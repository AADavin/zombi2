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
from zombi2.rates.modifiers import ByLineage
from zombi2.sequences import WIRED_MODIFIERS, simulate_sequences
from zombi2.sequences.substitution_models import (
    dayhoff, gtr, hky85, jc69, jtt, k80, lg, poisson, wag,
)
from zombi2.species import read_newick
from zombi2.cli.framework import (_add_flat_arg, _add_quiet_arg, _add_from_arg, _add_params_arg, _add_run_arg,
                                  _rate, _rates_help, _write_params_log, level_dir,
                                  resolve_genomes)

#: the RATES block for ``zombi2 sequences -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--substitution",
    note="ByLineage on --substitution IS the uncorrelated ('relaxed') clock: one i.i.d. multiplier "
         "per species lineage, shared across gene families. spread is σ; dist is 'lognormal' "
         "(default, σ = the log-scale) or 'gamma' (σ = the coefficient of variation).")

# the write vocabulary, mirroring SequencesResult.write (there is no exported constant to import)
_SEQUENCE_OUTPUTS = ("alignments", "phylograms", "ancestral", "founding", "species_phylogram")

# the menu, by alphabet: the no-argument protein models are empirical (their exchangeabilities and
# frequencies come from the published matrices), so each is just its constructor.
_NUCLEOTIDE_MODELS = ("jc69", "k80", "hky85", "gtr")
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
    g.add_argument("--model", required=True, metavar="MODEL",
                   choices=(*_NUCLEOTIDE_MODELS, *_PROTEIN_MODELS),
                   help="substitution model. nucleotide (4 states, ACGT): jc69 (equal rates), "
                        "k80 (--kappa), hky85 (--kappa, --frequencies), gtr (--gtr-rates, "
                        "--frequencies). protein (20 states): poisson (equal rates), jtt, dayhoff, "
                        "wag, lg — empirical matrices, no parameters to give")
    g.add_argument("--length", type=int, default=1000, metavar="N",
                   help="alignment length in sites — residues under a protein model (default 1000)")
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
                   help="which outputs to write (default: alignments, phylograms). also available: "
                        "ancestral (internal-node sequences), founding (each family's sequence at "
                        "its origination, where the phylogram's root branch starts), "
                        "species_phylogram (the species tree "
                        "scaled by the clock)")
    _add_flat_arg(g)
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
    # reject a physical parameter given for a model that doesn't read it (e.g. --kappa with jc69),
    # so a silently-ignored flag can't give a misleading run — the genomes command's discipline
    allowed = set(_MODEL_KNOBS[args.model])
    stray = [_KNOB_FLAG[k] for k in ("kappa", "frequencies", "gtr_rates")
             if getattr(args, k) is not None and k not in allowed]
    if stray:
        parser.error(f"these options don't apply to --model {args.model}: {', '.join(stray)}")

    handoff, tree_path = resolve_genomes(args.source or args.run)
    events_path = os.path.join(handoff, "genome_events.tsv")
    with open(tree_path) as f:
        tree, _ = read_newick(f.read())
    try:
        with open(events_path) as f:
            events = events_from_tsv(f.read())
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{events_path} not found — re-run 'zombi2 genomes' with 'events' in --write so the "
            "gene genealogy can be replayed") from None

    # Rebuild the genome run's spine from disk: its gene trees derive from (events, tree), and the
    # species tree drives the species phylogram. The sequence engine reads only .complete_tree and
    # .gene_trees, so an empty `genomes` map is the honest minimal shell (never escapes this run).
    genome_run = GenomesResult(complete_tree=tree, genomes={}, events=events, seed=None)

    model = _build_model(args)

    t0 = time.perf_counter()
    result = simulate_sequences(genome_run, model=model, length=args.length,
                                substitution=args.substitution, seed=args.seed,
                                progress=not args.quiet)
    dt = time.perf_counter() - t0

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "sequences", args.flat)
    # alignments and phylograms are one file per family, so a hundred families is hundreds of files
    # each — they get a directory apiece unless --flat says otherwise
    wanted = args.write or ("alignments", "phylograms")   # SequencesResult.write's own default
    if rest := [o for o in wanted if o not in ("alignments", "phylograms")]:
        result.write(out, outputs=rest)
    for per_family in ("alignments", "phylograms"):
        if per_family in wanted:
            result.write(level_dir(out, per_family, args.flat), outputs=(per_family,))

    n_families = sum(1 for aln in result.alignments.values() if aln)
    n_seqs = sum(len(aln) for aln in result.alignments.values())
    # the clock is now read off the rate itself: a ByLineage modifier is the relaxed clock
    clocks = [m for m in getattr(args.substitution, "modifiers", ()) if isinstance(m, ByLineage)]
    clock = (f"{clocks[0].dist} lineage clock, spread {clocks[0].spread:g}" if clocks
             else "strict clock")
    summary = (f"{n_seqs} sequences across {n_families} gene families, {model.name} "
               f"{args.length} sites, {clock}")
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(out, "sequences.log"),
                      args, summary)
    return 0
