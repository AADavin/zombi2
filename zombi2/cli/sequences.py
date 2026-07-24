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
:func:`zombi2.sequences.simulate_sequences`."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.genomes import GenomesResult
from zombi2.genomes.events import events_from_tsv
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.rates.modifiers import ByLineage
from zombi2.sequences import WIRED_MODIFIERS, simulate_sequences
from zombi2.sequences.substitution_models import (
    dayhoff, gtr, hky85, jc69, jtt, k80, lg, poisson, wag,
)
from zombi2.tree import read_newick
from zombi2.cli.framework import (_add_flat_arg, _add_quiet_arg, _add_parallel_arg, _add_from_arg,
                                  _add_params_arg, _add_run_arg, _rate, _rates_help, _write_params_log,
                                  default_outputs, level_dir, parallel_from_args, resolve_genomes)

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
                     "genomes", "initial_genome")

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
                        "wag, lg — empirical matrices, no parameters to give")
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
    g.add_argument("--substitution", type=_rate, default=1.0, metavar="RATE",
                   help="per-site substitution rate: a gene-tree branch of Δt time accrues "
                        "substitution·Δt substitutions/site (default 1.0 — the strict clock). A "
                        "ByLineage modifier makes it a relaxed clock: \"1.0 * ByLineage(spread=0.3)\"")

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
    _add_quiet_arg(g)


def _build_model(args: argparse.Namespace):
    """Build the substitution model from ``--model`` and its physical parameters (each knob falls back
    to the menu constructor's own default when not given; a protein model takes none)."""
    if args.model in _PROTEIN_MODELS:
        return _PROTEIN_MODELS[args.model]()
    kappa = 2.0 if args.kappa is None else args.kappa
    freqs = (0.25, 0.25, 0.25, 0.25) if args.frequencies is None else tuple(args.frequencies)
    rates = (1, 1, 1, 1, 1, 1) if args.gtr_rates is None else tuple(args.gtr_rates)
    if args.model == "jc69":
        return jc69()
    if args.model == "k80":
        return k80(kappa=kappa)
    if args.model == "hky85":
        return hky85(kappa=kappa, freqs=freqs)
    return gtr(rates=rates, freqs=freqs)


def run(args, parser):
    # validated here (not as argparse `required`) so a --params file can supply it
    if args.model is None:
        parser.error("--model is required (give it on the command line or in --params)")
    # reject a physical parameter given for a model that doesn't read it (e.g. --kappa with jc69),
    # so a silently-ignored flag can't give a misleading run — the genomes command's discipline
    allowed = set(_MODEL_KNOBS[args.model])
    stray = [_KNOB_FLAG[k] for k in ("kappa", "frequencies", "gtr_rates")
             if getattr(args, k) is not None and k not in allowed]
    if stray:
        parser.error(f"these options don't apply to --model {args.model}: {', '.join(stray)}")

    handoff, tree_path = resolve_genomes(args.source or args.run)
    with open(tree_path) as f:
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
            with open(events_path) as f:
                events = events_from_tsv(f.read())
        except FileNotFoundError:
            raise FileNotFoundError(
                f"{events_path} not found — re-run 'zombi2 genomes' with 'events' in --write so the "
                "gene genealogy can be replayed") from None
        # The genome run's spine from disk: its gene trees derive from (events, tree), and the species
        # tree drives the species phylogram. The sequence engine reads only .complete_tree and
        # .gene_trees, so an empty `genomes` map is the honest minimal shell (it never escapes here).
        genome_run = GenomesResult(complete_tree=tree, genomes={}, events=events, seed=None)
        for flag, value in (("--intergene-model", args.intergene_model),):
            if value is not None:
                parser.error(f"{flag} applies to a nucleotide genome run, where blocks are genes or "
                             "spacer. This handoff has gene families only, so there is nothing for a "
                             "second model to evolve.")
        extra = dict(length=1000 if args.length is None else args.length)

    model = _build_model(args)

    t0 = time.perf_counter()
    result = simulate_sequences(genome_run, model=model, substitution=args.substitution,
                                seed=args.seed, parallel=parallel_from_args(args, parser),
                                progress=not args.quiet, **extra)
    dt = time.perf_counter() - t0

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "sequences", args.flat)
    # the many-files-per-run outputs get a directory apiece (unless --flat): alignments and
    # phylograms are one file per family — per *block* on a nucleotide run, where a real genome has
    # thousands — and the assembled genome FASTAs are one per node. `initial_genome` is a single
    # file, but it is a whole-genome FASTA like the rest, so it lands in genomes/ with them.
    wanted = tuple(args.write) if args.write else default_outputs(result)
    own_dir = {"alignments": "alignments", "phylograms": "phylograms",
               "genomes": "genomes", "initial_genome": "genomes"}
    if rest := [o for o in wanted if o not in own_dir]:
        result.write(out, outputs=rest)
    for token, sub in own_dir.items():
        if token in wanted:
            result.write(level_dir(out, sub, args.flat), outputs=(token,))

    n_families = sum(1 for aln in result.alignments.values() if aln)
    n_seqs = sum(len(aln) for aln in result.alignments.values())
    # the clock is now read off the rate itself: a ByLineage modifier is the relaxed clock
    clocks = [m for m in getattr(args.substitution, "modifiers", ()) if isinstance(m, ByLineage)]
    clock = (f"{clocks[0].dist} lineage clock, spread {clocks[0].spread:g}" if clocks
             else "strict clock")
    if nucleotide:
        # the assembled genome of a node is exactly as long as its block layout (substitution keeps
        # length), so total bp comes from the genome run without assembling every node's sequence —
        # which, since `result.genomes` is now assembled lazily, would otherwise build them all just
        # to sum their lengths.
        bp = sum(g.length for g in genome_run.genomes.values())
        spacer = args.intergene_model or "jc69"
        summary = (f"{n_seqs} sequences across {n_families} blocks, {bp:,} bp assembled into "
                   f"{len(result.genomes)} genomes (every node), {model.name} genes / {spacer} spacer at "
                   f"{args.intergene_speed:g}x, {clock}")
    else:
        summary = (f"{n_seqs} sequences across {n_families} gene families, {model.name} "
                   f"{extra['length']} sites, {clock}")
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(out, "sequences.log"),
                      args, summary)
    return 0
