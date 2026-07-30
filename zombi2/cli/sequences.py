"""``zombi2 sequences`` — evolve a sequence inside each gene, along its gene tree.

A sequence sees the species tree only through its gene tree, so this command takes a **prior genomes
run** (``--genomes DIR``) and replays its gene genealogy: it reads that directory's
``genome_species_tree.nwk`` and ``genome_events.tsv``, rebuilds the ``{family: GeneTree}`` the run
produced, and evolves one sequence down each family's *complete* gene tree under a substitution
**model** (the menu — nucleotide ``jc69`` · ``k80`` · ``hky85`` · ``gtr``, or protein ``poisson`` ·
``jtt`` · ``dayhoff`` · ``wag`` · ``lg``) at a per-site substitution **rate**.

Long options are the API keyword names, and ``--substitution`` takes the written form of a rate
(SPEC §5): a bare number is the strict clock, and the uncorrelated ("relaxed") lineage clock is that
rate times a ``ByLineage`` modifier — ``--substitution "1.0 * ByLineage(spread=0.3)"``. The model's
physical parameters (``--kappa`` / ``--frequencies`` / ``--gtr-rates``) are rejected for a model that
does not use them — including *every* protein model, which is empirical and takes none — so a
silently-ignored flag can't give a misleading run. See
`zombi2.sequences.simulate_sequences()`."""
from __future__ import annotations

import argparse
import os
import time

import numpy as np

from zombi2.genomes import FamilyGenomesResult
from zombi2.genomes.events import events_from_tsv
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.rates.modifiers import ByLineage, Modifier
from zombi2.rates.parse import written_form
from zombi2.rates.rate import as_rate
from zombi2._runtime.report import write_run_report
from zombi2.rates.scope import PerSite
from zombi2.sequences import (WIRED_MODIFIERS, _calibrate, mean_pairwise_identity,
                              simulate_sequences)
from zombi2.sequences.substitution_models import (
    dayhoff, gtr, hky85, jc69, jtt, k80, lg, poisson, wag,
)
from zombi2.tree import read_newick
from zombi2.cli.framework import (_add_flat_arg, _add_force_arg, _add_quiet_arg, _add_parallel_arg, _add_from_arg,
                                  _add_params_arg, _add_run_arg, _rate, _rates_help, _write_params_log,
                                  default_outputs, guidance, level_dir, parallel_from_args,
                                  defaults_used, input_digests, resolve_genomes, resolve_seed, warn)

#: the RATES block for ``zombi2 sequences -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--substitution",
    note="ByLineage on --substitution IS the uncorrelated ('relaxed') clock: one i.i.d. multiplier "
         "per species lineage, shared across gene families. spread is σ; dist is 'lognormal' "
         "(default, σ = the log-scale) or 'gamma' (σ = the coefficient of variation).")

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
    "gtr": ("frequencies", "gtr_rates"),
    **{name: () for name in _PROTEIN_MODELS},
}
_KNOB_FLAG = {"kappa": "--kappa", "frequencies": "--frequencies", "gtr_rates": "--gtr-rates"}


def _add_sequence_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "sequences evolve down the gene trees of the genomes run it already holds")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    _add_from_arg(g, "the genomes run to replay — its genome_species_tree.nwk and "
                     "genome_events.tsv rebuild the gene trees")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")

    g = p.add_argument_group("substitution model", "the menu — one --model, its parameters below")
    # validated in run() (not argparse-`required`) so a --params file can supply it — a required
    # argument is never satisfied by a default, which is what --params sets.
    g.add_argument("--model", default=None, metavar="MODEL",
                   choices=(*_NUCLEOTIDE_MODELS, *_PROTEIN_MODELS),
                   help="substitution model. nucleotide (4 states, ACGT): jc69 (equal rates), "
                        "k80 (--kappa), hky85 (--kappa, --frequencies), gtr (--gtr-rates, "
                        "--frequencies). protein (20 states): poisson (equal rates), jtt, dayhoff, "
                        "wag, lg — empirical matrices, no parameters to give. Defaults to jc69, "
                        "announced when used")
    g.add_argument("--length", type=int, default=None, metavar="N",
                   help="alignment length in sites — residues under a protein model (default 1000). "
                        "Not for a nucleotide genome run: there every block carries its own length "
                        "in bp, so giving one here is an error rather than something ignored")
    g.add_argument("--intergene-model", default=None, metavar="MODEL", dest="intergene_model",
                   choices=_NUCLEOTIDE_MODELS,
                   help="[nucleotide runs] the model the spacer between genes evolves under "
                        "(default jc69 — flat, no free parameters). Genes take --model")
    g.add_argument("--intergene-speed", type=float, default=3.0, metavar="X",
                   dest="intergene_speed",
                   help="[nucleotide runs] how much faster the spacer evolves than the genes, as a "
                        "multiple of the substitution rate (default 3.0)")
    g.add_argument("--kappa", type=float, default=None, metavar="K",
                   help="[k80, hky85] transition/transversion ratio (default 2.0)")
    g.add_argument("--frequencies", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[hky85, gtr] equilibrium base frequencies A C G T (must be positive and "
                        "sum to 1; default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None, dest="gtr_rates",
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[gtr] the six exchangeabilities (default all 1)")

    g = p.add_argument_group("substitution rate & clock", "the per-site rate — see RATES below")
    g.add_argument("--substitution", type=_rate, default=None, metavar="RATE",
                   help="per-site substitution rate: a gene-tree branch of Δt time accrues "
                        "substitution·Δt substitutions/site (default 1.0 — the strict clock). A "
                        "ByLineage modifier makes it a relaxed clock: \"1.0 * ByLineage(spread=0.3)\"")
    g.add_argument("--divergence", type=float, default=None, metavar="D",
                   help="pick the rate for me, so that a site accrues D substitutions from the root "
                        "to a tip. The rate is per unit TIME, so what it produces depends on how "
                        "tall the tree is — 0.2 here is a readable alignment on any tree, where a "
                        "rate of 1.0 is fine on a short tree and pure noise on a tall one. Composes "
                        "with --substitution: give the clock's shape alone "
                        "(\"ByLineage(spread=0.3)\") and this sets its scale")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_SEQUENCE_OUTPUTS, default=None, metavar="PART",
                   help="which outputs to write (default: alignments, phylograms, "
                        "species_phylogram — the last written as clock_species_tree_*.nwk, the "
                        "species tree with its branches in substitutions/site — and, on a "
                        "nucleotide run, genomes and initial_genome: one assembled FASTA per node "
                        "of the complete tree, plus the genome the run started with. 'genomes' is "
                        "the big one, a whole genome times every node. also available: ancestral "
                        "(the sequence at every node that is not an extant tip) and founding (each "
                        "family's sequence at its origination, where the phylogram's root branch "
                        "starts)")
    _add_flat_arg(g)
    _add_parallel_arg(g)
    g.add_argument("--stream", action="store_true",
                   help="[family/ordered] write each family's sequences straight to disk instead of "
                        "building the whole run in memory. This is the level where a run's memory "
                        "goes — every alignment and every ancestral sequence live at once — so this "
                        "is the dial for a run that would not fit. Composes with --parallel and "
                        "--write, and writes the same files at the same seed, so a streamed run and "
                        "an in-memory one are the same dataset. Not available on a nucleotide run, "
                        "which reassembles whole genomes and so needs every block at once")
    _add_quiet_arg(g)
    _add_force_arg(g)


def _resolve_model_knobs(args) -> dict:
    """The nucleotide substitution-model knobs with each default filled in — ``kappa`` (2.0),
    ``frequencies`` (uniform), ``gtr_rates`` (all 1). Shared by `_build_model()`, which uses them,
    and `_effective_model_params()`, which logs them, so the two cannot drift."""
    # floats, not ints: --gtr-rates parses as float, so a default logged as [1, 1, …] would not match
    # the [1.0, …] a reproduced run logs — and the run report is meant to reproduce byte-for-byte.
    return {"kappa": 2.0 if args.kappa is None else args.kappa,
            "frequencies": [0.25, 0.25, 0.25, 0.25] if args.frequencies is None else list(args.frequencies),
            "gtr_rates": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] if args.gtr_rates is None else list(args.gtr_rates)}


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
    knowing each model's defaults. ``jc69`` and the protein models have no free knob, so this is empty."""
    resolved = _resolve_model_knobs(args)
    return {name: resolved[name] for name in _MODEL_KNOBS.get(args.model, ())}


def _build_model(args: argparse.Namespace):
    """Build the substitution model from ``--model`` and its physical parameters (each knob falls back
    to the menu constructor's own default when not given; a protein model takes none)."""
    if args.model in _PROTEIN_MODELS:
        return _PROTEIN_MODELS[args.model]()
    knobs = _resolve_model_knobs(args)
    kappa, freqs, rates = knobs["kappa"], tuple(knobs["frequencies"]), tuple(knobs["gtr_rates"])
    if args.model == "jc69":
        return jc69()
    if args.model == "k80":
        return k80(kappa=kappa)
    if args.model == "hky85":
        return hky85(kappa=kappa, freqs=freqs)
    return gtr(rates=rates, freqs=freqs)


#: Below this much of the way from the random floor to identical, an alignment carries so little
#: signal that homology search and tree inference will fail on it. 0.15 sits well clear of both
#: sides in practice: a saturated DNA run lands near 0.04–0.09, a usable one above 0.7.
_SATURATED_BELOW = 0.15


def _saturation_signal(identity: float, model) -> float:
    """How far the realised identity sits from random, as a fraction of the distance from the model's
    own random floor to identical. Two sequences related only by chance still match at ``Σπ²`` — 25%
    for equal-frequency DNA, ~6% for a protein model — so raw identity is not comparable across
    models and this is."""
    floor = float(np.sum(np.asarray(model.stationary) ** 2))
    return (identity - floor) / (1.0 - floor) if floor < 1.0 else 1.0


def run(args, parser):
    # validated here (not as argparse `required`) so a --params file can supply it
    if args.model is None:
        # jc69 has no free parameters, so a bare run needs nothing else to be well defined
        warn(defaults_used(args, model="jc69"))
    resolve_seed(args)                     # a run must be reproducible from its own log
    # reject a physical parameter given for a model that doesn't read it (e.g. --kappa with jc69),
    # so a silently-ignored flag can't give a misleading run — the genomes command's discipline
    allowed = set(_MODEL_KNOBS[args.model])
    stray = [_KNOB_FLAG[k] for k in ("kappa", "frequencies", "gtr_rates")
             if getattr(args, k) is not None and k not in allowed]
    if stray:
        parser.error(f"these options don't apply to --model {args.model}: {', '.join(stray)}")

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
                events = events_from_tsv(f.read())
        except FileNotFoundError:
            raise FileNotFoundError(
                f"{events_path} not found — re-run 'zombi2 genomes' with 'events' in --write so the "
                "gene genealogy can be replayed") from None
        # The genome run's spine from disk: its gene trees derive from (events, tree), and the species
        # tree drives the species phylogram. The sequence engine reads only .complete_tree and
        # .gene_trees, so an empty `genomes` map is the honest minimal shell (it never escapes here).
        genome_run = FamilyGenomesResult(complete_tree=tree, genomes={}, events=events, seed=None)
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
    clocks = [m for m in _mods if isinstance(m, ByLineage)]
    clock = (f"{clocks[0].dist} lineage clock, spread {clocks[0].spread:g}" if clocks
             else "strict clock")
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
        floor = float(np.sum(np.asarray(model.stationary) ** 2))
        # Name the rate that RAN, not the flag: --substitution is None when it was left out, and
        # "currently None" tells a reader nothing — and under --divergence the rate that ran is the
        # solved-for base, which is not on the command line at all. Then point at --divergence first:
        # it is the answer to "what should I put instead", which lowering a rate by guesswork is not.
        rate = _effective_substitution(args, genome_run).get("substitution", args.substitution)
        used = ("the default 1.0" if args.substitution is None and args.divergence is None
                else written_form(as_rate(rate, default_scope=PerSite)))
        warn(f"these sequences are close to saturated — mean pairwise identity is {identity:.1%}, "
             f"against {floor:.1%} for unrelated sequences under {model.name}. The substitution "
             f"rate is per unit time, so a tall tree accrues many substitutions per site and the "
             f"alignments keep little history: homology search and tree inference will both do "
             f"poorly on them. Say how diverged you want them instead — --divergence 0.2 is a "
             f"readable alignment on any tree — or lower the rate yourself (it ran at {used}).")
    guidance(args, f"alignments under {out}/")
    _write_params_log(os.path.join(out, "sequences.log"), args, summary,
                      effective={"write": list(wanted), **_effective_model_params(args),
                                 **_effective_substitution(args, genome_run)},
                      inputs=input_digests(tree_path,
                                           os.path.join(handoff, "genome_events.tsv"),
                                           os.path.join(handoff, "blocks.tsv"),
                                           args.substitution))
    if path := write_run_report(args.run):     # refresh the run's one-page report (grouped layout only)
        guidance(args, f"run report (one-page summary of the whole run): {path}")
    return 0
