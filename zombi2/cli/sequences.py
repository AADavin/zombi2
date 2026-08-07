"""``zombi2 sequences`` — evolve a sequence inside each gene, along its gene tree.

A sequence sees the species tree only through its gene tree, so this command takes a **prior genomes
run** (the run directory, or ``--from PATH``) and replays its gene genealogy: it reads that directory's
``genome_species_tree.nwk`` and ``genome_events.tsv``, rebuilds the ``{family: GeneTree}`` the run
produced, and evolves one sequence down each family's *complete* gene tree under a substitution
**model** (the menu — nucleotide ``jc69`` · ``k80`` · ``hky85`` · ``gtr``, or protein ``poisson`` ·
``jtt`` · ``dayhoff`` · ``wag`` · ``lg``) at a per-site substitution **rate**.

Long options are the API keyword names, and ``--substitution`` takes the written form of a rate
(SPEC §5): a bare number is the strict clock, and the uncorrelated ("relaxed") lineage clock is that
rate times a ``ByLineage`` modifier — ``--substitution "1.0 * ByLineage(spread=0.3)"``. The model's
physical parameters (``--kappa`` / ``--frequencies`` / ``--exchangeabilities``) are rejected for a model that
does not use them — including *every* protein model, which is empirical and takes none — so a
silently-ignored flag can't give a misleading run. See
`zombi2.sequences.simulate_sequences()`."""
from __future__ import annotations

import argparse
import os
import time

import numpy as np

from zombi2.genomes import FamilyGenomesResult
from zombi2.genomes.events import edges_from_tsv
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.rates.modifiers import DRAWN, INHERITED, DrivenBy, Modifier
from zombi2._runtime.report import write_run_report
from zombi2.sequences import (IMPLEMENTED_MODIFIERS, _calibrate, mean_pairwise_identity,
                              simulate_sequences)
from zombi2.sequences.substitution_models import (
    dayhoff, gtr, hky85, jc69, jtt, k80, lg, poisson, wag,
)
from zombi2.tree import read_newick
from zombi2.cli.framework import (_add_flat_arg, _add_force_arg, _add_quiet_arg, _add_parallel_arg, _add_from_arg,
                                  _add_params_arg, _add_run_arg, _rate, _rates_help, _write_params_log,
                                  conditioned_levels, default_outputs, signpost, level_dir, parallel_from_args,
                                  defaults_used, input_digests, record_conditioning, resolve_genomes,
                                  resolve_seed, warn)

#: the RATES block for ``zombi2 sequences -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    IMPLEMENTED_MODIFIERS, "--substitution",
    note="ByLineage draws one rate per species lineage, shared by every gene in it. spread is σ; "
         "dist is 'lognormal' (default) or 'gamma' (σ = the coefficient of variation). DrivenBy "
         "reads a trait grown first — the trait_events.tsv a 'zombi2 traits' run wrote, in this run "
         "or another: \"1.0 * DrivenBy('out/traits/trait_events.tsv', {'cave': 0.5, 'surface': "
         "1.0})\". A clock and a driver compose; a driver that switches mid-branch is integrated "
         "across the switch, not sampled once for the branch.")

# the write vocabulary, mirroring SequencesResult.write (there is no exported constant to import).
# The last two exist only for a nucleotide handoff, which is the only run with coordinates to lay a
# genome out in; asking for one otherwise writes nothing rather than failing.
_SEQUENCE_OUTPUTS = ("alignments", "phylograms", "ancestral", "founding", "species_phylogram",
                     "genomes", "initial_genome", "summary")

# the menu, by alphabet: the no-argument protein models are empirical (their exchangeabilities and
# frequencies come from the published matrices), so each is just its constructor.
_NUCLEOTIDE_MODELS = ("jc69", "k80", "hky85", "gtr")
#: the nucleotide constructors by name, for --intergene-model (which takes no parameters of its own:
#: the spacer's job is to be the unconstrained null, and a second set of knobs would only blur that)
_NUCLEOTIDE_CTORS = {"jc69": jc69, "k80": k80, "hky85": hky85, "gtr": gtr}
_PROTEIN_MODELS = {"poisson": poisson, "jtt": jtt, "dayhoff": dayhoff, "wag": wag, "lg": lg}

# which physical parameters each model reads; a knob given for a model that does not take it is
# rejected. Every model knob defaults to None, so "given" is simply "not None". A protein model
# reads none — its matrix is published, not tuned.
_MODEL_KNOBS = {
    "jc69": (),
    "k80": ("kappa",),
    "hky85": ("kappa", "frequencies"),
    "gtr": ("frequencies", "exchangeabilities"),
    **{name: () for name in _PROTEIN_MODELS},
}
_KNOB_FLAG = {"kappa": "--kappa", "frequencies": "--frequencies",
              "exchangeabilities": "--exchangeabilities"}


def _add_sequence_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "sequences evolve down the gene trees of the genomes run it holds")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    _add_from_arg(g, "the genomes run to replay")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")

    g = p.add_argument_group("substitution model")
    # validated in run() (not argparse-`required`) so a --params file can supply it — a required
    # argument is never satisfied by a default, which is what --params sets.
    g.add_argument("--model", default=None, metavar="MODEL",
                   choices=(*_NUCLEOTIDE_MODELS, *_PROTEIN_MODELS),
                   # which knob goes with which model is on the knobs themselves ([k80, hky85] …)
                   help="substitution model (default jc69). nucleotide: jc69, k80, hky85, gtr. "
                        "protein: poisson, jtt, dayhoff, wag, lg — empirical, no parameters to give")
    g.add_argument("--length", type=int, default=None, metavar="N",
                   help="alignment length in sites (default 1000). Rejected on a nucleotide run: "
                        "each block has its own length")
    g.add_argument("--intergene-model", default=None, metavar="MODEL", dest="intergene_model",
                   choices=_NUCLEOTIDE_MODELS,
                   help="[nucleotide runs] model for the spacer between genes (default jc69); "
                        "genes take --model")
    g.add_argument("--intergene-speed", type=float, default=3.0, metavar="X",
                   dest="intergene_speed",
                   help="[nucleotide runs] spacer rate as a multiple of the gene rate (default 3.0)")
    g.add_argument("--kappa", type=float, default=None, metavar="K",
                   help="[k80, hky85] transition/transversion ratio (default 2.0)")
    g.add_argument("--frequencies", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[hky85, gtr] equilibrium base frequencies A C G T (must be positive and "
                        "sum to 1; default equal)")
    g.add_argument("--exchangeabilities", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[gtr] the six exchangeabilities (default all 1)")
    # Across-site rate variation decorates ANY model on the menu, so these sit outside _MODEL_KNOBS
    # (which is the per-matrix physical parameters, and whose stray-knob check is per model).
    g.add_argument("--gamma-shape", type=float, default=None, metavar="A", dest="gamma_shape",
                   help="[any model] Gamma shape a for rate variation across sites (+G): smaller is "
                        "more unequal — 0.5 is strong, 5 nearly flat")
    g.add_argument("--invariant", type=float, default=None, metavar="P",
                   help="[any model] proportion of sites that never change (+I); 0 by default")
    g.add_argument("--rate-categories", type=int, default=None, metavar="N", dest="rate_categories",
                   help="[any model] how many discrete Gamma classes (default 4); needs --gamma-shape")

    g = p.add_argument_group("substitution rate & clock", "see RATES below")
    g.add_argument("--substitution", type=_rate, default=None, metavar="RATE",
                   help="substitutions per site per unit time (default 1.0, a strict clock); a "
                        "ByLineage modifier relaxes it, and a DrivenBy reads a trait grown first")
    g.add_argument("--divergence", type=float, default=None, metavar="D",
                   help="solve for the rate instead, so a site accrues D substitutions from root to "
                        "tip. Composes with --substitution: give the clock's shape alone "
                        "(\"ByLineage(spread=0.3)\") and this sets its scale")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_SEQUENCE_OUTPUTS, default=None, metavar="PART",
                   help="which outputs to write. default: alignments, phylograms, "
                        "species_phylogram, summary, and — on a nucleotide run — genomes (one "
                        "assembled FASTA per node, the big one) and initial_genome. also: "
                        "ancestral, founding")
    _add_flat_arg(g)
    _add_parallel_arg(g)
    g.add_argument("--stream", action="store_true",
                   help="[family/ordered] write each family's sequences straight to disk instead of "
                        "holding the whole run in memory. Not available on a nucleotide run")
    _add_quiet_arg(g)
    _add_force_arg(g)


def _resolve_model_knobs(args) -> dict:
    """The nucleotide substitution-model knobs with each default filled in — ``kappa`` (2.0),
    ``frequencies`` (uniform), ``exchangeabilities`` (all 1). Shared by `_build_model()`, which uses them,
    and `_effective_model_params()`, which logs them, so the two cannot drift."""
    # floats, not ints: --exchangeabilities parses as float, so a default logged as [1, 1, …] would not match
    # the [1.0, …] a reproduced run logs — and the run report is meant to reproduce byte-for-byte.
    return {"kappa": 2.0 if args.kappa is None else args.kappa,
            "frequencies": [0.25, 0.25, 0.25, 0.25] if args.frequencies is None else list(args.frequencies),
            "exchangeabilities": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] if args.exchangeabilities is None else list(args.exchangeabilities)}


def _effective_substitution(args, genome_run) -> dict:
    """The substitution rate the run actually used, for the log.

    With ``--divergence`` the base is solved for, so the rate written on the command line is not the
    rate that ran — and a log recording the flag rather than the resolved value would be a provenance
    record that disagrees with the run. Calls the same `_calibrate()` the engine
    does, so the logged number cannot drift from the one used."""
    if args.divergence is not None:
        return {"substitution": _calibrate(args.substitution, args.divergence,
                                           genome_run.complete_tree)}
    # No divergence: the flag may still have been left out, in which case the run used 1.0 and the
    # log should say so. A field logged None because a default applied is a field that stops meaning
    # anything the day the default changes — which is the whole reason to dump every parameter.
    return {} if args.substitution is not None else {"substitution": 1.0}


def _effective_model_params(args) -> dict:
    """The substitution-model params the run actually used, defaults filled — but only the knobs this
    model has (`_MODEL_KNOBS`), so the ``.log`` reproduces the exact model without the reader
    knowing each model's defaults. ``jc69`` and the protein models have no free knob, so this is empty.

    The across-site knobs are appended when either was given, resolved rather than as typed: someone
    reading the log should not have to know that the category count defaults to 4."""
    resolved = _resolve_model_knobs(args)
    params = {name: resolved[name] for name in _MODEL_KNOBS.get(args.model, ())}
    if args.gamma_shape is not None or args.invariant:
        params["invariant"] = args.invariant or 0.0
        if args.gamma_shape is not None:
            params["gamma_shape"] = args.gamma_shape
            params["rate_categories"] = args.rate_categories or 4
    return params


def _build_model(args: argparse.Namespace):
    """Build the substitution model from ``--model`` and its physical parameters (each knob falls back
    to the menu constructor's own default when not given; a protein model takes none), then decorate
    it with across-site rate variation if any was asked for.

    The decoration is applied last and to whatever matrix came out, because it is orthogonal to the
    chemistry: every model on the menu takes it, which is why the three flags are not per-model knobs."""
    if args.model in _PROTEIN_MODELS:
        model = _PROTEIN_MODELS[args.model]()
    else:
        knobs = _resolve_model_knobs(args)
        kappa = knobs["kappa"]
        frequencies = tuple(knobs["frequencies"])
        exchangeabilities = tuple(knobs["exchangeabilities"])
        if args.model == "jc69":
            model = jc69()
        elif args.model == "k80":
            model = k80(kappa=kappa)
        elif args.model == "hky85":
            model = hky85(kappa=kappa, frequencies=frequencies)
        else:
            model = gtr(exchangeabilities=exchangeabilities, frequencies=frequencies)
    if args.gamma_shape is not None or args.invariant:
        model = model.across_sites(gamma_shape=args.gamma_shape,
                                   invariant=args.invariant or 0.0,
                                   rate_categories=args.rate_categories or 4)
    return model


#: Below this much of the way from the random floor to identical, an alignment carries so little
#: signal that homology search and tree inference will fail on it. 0.15 sits well clear of both
#: sides in practice: a saturated DNA run lands near 0.04–0.09, a usable one above 0.7.
_SATURATED_BELOW = 0.15


def _identity_floor(model) -> float:
    """The identity two *unrelated* sequences would still show under this model.

    ``Σπ²`` — 25% for equal-frequency DNA, ~6% for a protein model — because two random draws from
    the stationary distribution agree that often. Under ``+I`` the floor is higher: the invariant
    sites match no matter how long the branch, so the floor is ``p + (1 - p)·Σπ²``. Without that
    term the saturation warning below would fire on a perfectly good ``+I`` run, whose identity
    cannot fall to ``Σπ²`` even in principle."""
    chance = float(np.sum(np.asarray(model.stationary) ** 2))
    invariant = model.site_shares[0] if model.site_rates[0] == 0.0 else 0.0
    return invariant + (1.0 - invariant) * chance


def _saturation_signal(identity: float, model) -> float:
    """How far the realised identity sits from random, as a fraction of the distance from the model's
    own floor (`_identity_floor`) to identical — so raw identity, which is not comparable across
    models, becomes a number that is."""
    floor = _identity_floor(model)
    return (identity - floor) / (1.0 - floor) if floor < 1.0 else 1.0


def run(args, parser):
    # validated here (not as argparse `required`) so a --params file can supply it
    if args.model is None:
        # jc69 has no free parameters, so a bare run needs nothing else to be well defined
        defaults_used(args, model="jc69")
    resolve_seed(args)                     # a run must be reproducible from its own log
    # reject a physical parameter given for a model that doesn't read it (e.g. --kappa with jc69),
    # so a silently-ignored flag can't give a misleading run — the genomes command's discipline
    allowed = set(_MODEL_KNOBS[args.model])
    stray = [_KNOB_FLAG[k] for k in ("kappa", "frequencies", "exchangeabilities")
             if getattr(args, k) is not None and k not in allowed]
    if stray:
        parser.error(f"these options don't apply to --model {args.model}: {', '.join(stray)}")
    # --rate-categories counts the classes of a Gamma, so on its own it asks for nothing; caught here
    # rather than in across_sites() so it reads as the flag error it is
    if args.rate_categories is not None and args.gamma_shape is None:
        parser.error("--rate-categories counts the discrete classes of the across-site Gamma, so it "
                     "needs --gamma-shape; --invariant on its own is a single never-changing class "
                     "and takes no count")

    handoff, tree_path = resolve_genomes(args.source or args.run)
    with open(tree_path, encoding="utf-8") as f:
        tree, _ = read_newick(f.read())

    # Which resolution wrote this handoff? blocks.tsv is the nucleotide resolution's and no other's,
    # so the run says what it is rather than needing a flag repeated from the genomes command.
    nucleotide = os.path.exists(os.path.join(handoff, "blocks.tsv"))
    if nucleotide:
        # Rebuild the whole run: at this resolution the sequences evolve down a tree per *block*, and
        # the blocks come from the genomes themselves, not from the event log alone.
        genome_run = read_nucleotide_genomes(handoff, tree)
        if args.length is not None:
            parser.error("--length does not apply to a nucleotide genome run: every block carries "
                         "its own length in bp, so one number here would contradict the coordinates "
                         "the genomes run wrote. Drop it — the genome sets the lengths.")
        extra = dict(intergene_speed=args.intergene_speed)
        if args.intergene_model is not None:
            extra["intergene_model"] = _NUCLEOTIDE_CTORS[args.intergene_model]()
    else:
        events_path = os.path.join(handoff, "genome_events.tsv")
        try:
            with open(events_path, encoding="utf-8") as f:
                events = edges_from_tsv(f.read())
        except FileNotFoundError:
            raise FileNotFoundError(
                f"{events_path} not found — re-run 'zombi2 genomes' with 'events' in --write so the "
                "gene genealogy can be replayed") from None
        # The genome run's spine from disk: its gene trees derive from (events, tree), and the species
        # tree drives the species phylogram. The sequence engine reads only .complete_tree and
        # .gene_trees, so an empty `genomes` map is the honest minimal shell (it never escapes here).
        genome_run = FamilyGenomesResult(complete_tree=tree, genomes={}, edges=events, seed=None)
        for flag, value in (("--intergene-model", args.intergene_model),):
            if value is not None:
                parser.error(f"{flag} applies to a nucleotide genome run, where blocks are genes or "
                             "spacer. This handoff has gene families only, so there is nothing for a "
                             "second model to evolve.")
        extra = dict(length=1000 if args.length is None else args.length)

    model = _build_model(args)

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "sequences", args.flat)
    streaming = args.stream and not nucleotide
    if args.stream and nucleotide:
        parser.error("--stream is not available on a nucleotide run: it reassembles every node's "
                     "genome, which needs every block's sequence at once. Drop --stream, or use "
                     "--resolution family at the genome level")

    t0 = time.perf_counter()
    result = simulate_sequences(genome_run, model=model, substitution=args.substitution,
                                divergence=args.divergence,
                                seed=args.seed, parallel=parallel_from_args(args, parser),
                                stream_to=out if streaming else None,
                                outputs=(tuple(args.write) if args.write else None)
                                if streaming else None,
                                flat=args.flat,
                                progress=not args.quiet, **extra)
    dt = time.perf_counter() - t0
    # the many-files-per-run outputs get a directory apiece (unless --flat): alignments and
    # phylograms are one file per family — per *block* on a nucleotide run, where a real genome has
    # thousands — and the assembled genome FASTAs are one per node. `initial_genome` is a single
    # file, but it is a whole-genome FASTA like the rest, so it lands in genomes/ with them.
    # The many-files-per-run outputs (alignments, phylograms, the assembled genome FASTAs) get a
    # subdirectory apiece, and `write` is where that is decided — so a run written from Python has
    # this layout too, and --flat is simply passed through. An output this run has none of writes
    # nothing and creates no directory, so a family run leaves no empty genomes/ behind.
    if streaming:                       # the engine wrote its own files and says which
        wanted = result.outputs
        n_families, n_seqs = result.n_families, result.n_sequences
    else:
        wanted = tuple(args.write) if args.write else default_outputs(result)
        result.write(out, outputs=wanted, flat=args.flat)
        n_families = sum(1 for aln in result.alignments.values() if aln)
        n_seqs = sum(len(aln) for aln in result.alignments.values())
    # the clock is now read off the rate itself: a ByLineage modifier is the relaxed clock
    # --substitution may now be a bare modifier (the clock's shape, with --divergence setting its
    # scale), which carries no `.modifiers` of its own — so look at the modifier itself as well, or a
    # relaxed run reports itself as strict.
    _sub = args.substitution
    _mods = (_sub,) if isinstance(_sub, Modifier) else getattr(_sub, "modifiers", ())
    # Both clock modifiers, not just the uncorrelated one: a FromParent run is autocorrelated, and
    # reporting it as "strict" was the same bug the ByLineage branch above was written to fix.
    clock = "strict clock"
    driven = []
    for m in _mods:
        if m.reads == (DRAWN, "lineage"):
            clock = f"{m.dist} lineage clock, spread {m.spread:g}"
        elif m.reads == (INHERITED, "lineage"):
            clock = (f"discrete-bin clock, {m.bins} bins, spread {m.spread:g}" if m.bins
                     else f"autocorrelated clock, spread {m.spread:g}")
        elif isinstance(m, DrivenBy):
            # a driver is a second factor, not a second clock — appended rather than replacing, or a
            # driven relaxed run would report itself as one or the other and never as both
            driven.append(os.path.basename(m.driver) if isinstance(m.driver, str)
                          else type(m.driver).__name__)
    if driven:
        clock += f", driven by {', '.join(driven)}"
    # What the run actually produced, not what was asked for: the rate is per unit time, so whether
    # it yields a usable alignment depends on the height of the tree it ran down, which the user has
    # no way to read off the flags. Reporting it turns a silent failure into a visible number.
    # a streamed run has no alignments to measure — it accumulated the same statistic as it wrote
    identity = result.identity if streaming else mean_pairwise_identity(result.alignments)
    realised = "" if identity is None else f", mean identity {identity:.1%}"
    if nucleotide:
        # the assembled genome of a node is exactly as long as its block layout (substitution keeps
        # length), so total bp comes from the genome run without assembling every node's sequence —
        # which, since `result.genomes` is now assembled lazily, would otherwise build them all just
        # to sum their lengths.
        bp = sum(g.length for g in genome_run.genomes.values())
        spacer = args.intergene_model or "jc69"
        summary = (f"{n_seqs} sequences across {n_families} blocks, {bp:,} bp assembled into "
                   f"{len(result.genomes)} genomes (every node), {model.name} genes / {spacer} spacer at "
                   f"{args.intergene_speed:g}x, {clock}{realised}")
    else:
        summary = (f"{n_seqs} sequences across {n_families} gene families, {model.name} "
                   f"{extra['length']} sites, {clock}{realised}")
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    if identity is not None and _saturation_signal(identity, model) < _SATURATED_BELOW:
        # the floor makes the identity readable: it is where unrelated sequences already sit under
        # this model (25% for equal-frequency DNA, ~6% for a protein one; higher under +I)
        floor = _identity_floor(model)
        warn(f"these sequences are close to saturated — mean identity {identity:.1%}, against a "
             f"{floor:.1%} floor. Set --divergence 0.2 instead.")
    if not args.flat:                             # record which same-run levels this run reads (if any),
        record_conditioning(out, conditioned_levels(  # so re-running one of them knows it orphans this
            args.run, (args.substitution,)))
    # the log is the run's parameters, not the parser's: the intergene knobs are for a nucleotide
    # handoff only (there is a spacer between genes to evolve); on a family/ordered run they never
    # applied, so recording them — and printing them in the reproduce command — is misleading.
    _write_params_log(os.path.join(out, "sequences.log"), args, summary,
                      omit=() if nucleotide else ("intergene_speed", "intergene_model"),
                      effective={"write": list(wanted), **_effective_model_params(args),
                                 **_effective_substitution(args, genome_run)},
                      inputs=input_digests(tree_path,
                                           os.path.join(handoff, "genome_events.tsv"),
                                           os.path.join(handoff, "blocks.tsv"),
                                           args.substitution))
    signpost(args, write_run_report(args.run), out)   # every file it wrote, then the run report
    return 0
