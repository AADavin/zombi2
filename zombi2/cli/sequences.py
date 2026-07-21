"""``zombi2 sequences`` — evolve a sequence inside each gene, along its gene tree.

A sequence sees the species tree only through its gene tree, so this command takes a **prior genomes
run** (``--genomes DIR``) and replays its gene genealogy: it reads that directory's
``genome_species_tree.nwk`` and ``genome_events.tsv``, rebuilds the ``{family: GeneTree}`` the run
produced, and evolves one sequence down each family's *complete* gene tree under a substitution
**model** (the menu — ``jc69`` · ``k80`` · ``hky85`` · ``gtr``) at a per-site substitution **rate**
(a bare number — the strict clock — optionally times a ``ByLineage`` lineage clock).

Long options are the API keyword names; the model's physical parameters (``--kappa`` / ``--frequencies``
/ ``--gtr-rates``) are rejected for a model that does not use them, so a silently-ignored flag can't
give a misleading run. See :func:`zombi2.sequences.simulate_sequences`."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.genomes import GenomesResult
from zombi2.genomes.events import events_from_tsv
from zombi2.rates.modifiers import ByLineage
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import gtr, hky85, jc69, k80
from zombi2.species import read_newick
from zombi2.cli.framework import _add_params_arg, _write_params_log

# the write vocabulary, mirroring SequencesResult.write (there is no exported constant to import)
_SEQUENCE_OUTPUTS = ("alignments", "phylograms", "ancestral", "species_phylogram")

# which physical parameters each model reads; a knob given for a model that does not take it is
# rejected. Every model knob defaults to None, so "given" is simply "not None".
_MODEL_KNOBS = {
    "jc69": (),
    "k80": ("kappa",),
    "hky85": ("kappa", "frequencies"),
    "gtr": ("frequencies", "gtr_rates"),
}
_KNOB_FLAG = {"kappa": "--kappa", "frequencies": "--frequencies", "gtr_rates": "--gtr-rates"}


def _add_sequence_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("--genomes", required=True, metavar="DIR",
                   help="a prior 'zombi2 genomes' output directory — reads its "
                        "genome_species_tree.nwk and genome_events.tsv and replays the gene "
                        "genealogy (sequences evolve down the complete gene trees, lost/extinct "
                        "lineages included)")
    g.add_argument("-o", "--output", required=True, metavar="DIR", dest="output",
                   help="output directory (created if needed)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")

    g = p.add_argument_group("substitution model", "the menu — one --model, its parameters below")
    g.add_argument("--model", required=True, choices=("jc69", "k80", "hky85", "gtr"), metavar="MODEL",
                   help="nucleotide substitution model: jc69 (equal rates), k80 (--kappa), "
                        "hky85 (--kappa, --frequencies), gtr (--gtr-rates, --frequencies)")
    g.add_argument("--length", type=int, default=1000, metavar="N",
                   help="alignment length in sites (default 1000)")
    g.add_argument("--kappa", type=float, default=None, metavar="K",
                   help="[k80, hky85] transition/transversion ratio (default 2.0)")
    g.add_argument("--frequencies", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[hky85, gtr] equilibrium base frequencies A C G T (must be positive and "
                        "sum to 1; default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None, dest="gtr_rates",
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[gtr] the six exchangeabilities (default all 1)")

    g = p.add_argument_group("substitution rate & clock", "bare-number per-site rate × a lineage clock")
    g.add_argument("--substitution", type=float, default=1.0, metavar="RATE",
                   help="per-site substitution rate: a gene-tree branch of Δt time accrues "
                        "substitution·Δt substitutions/site (default 1.0 — the strict clock)")
    g.add_argument("--clock-spread", type=float, default=0.0, metavar="SPREAD", dest="clock_spread",
                   help="lineage-clock spread (σ): one i.i.d. rate multiplier per species lineage, "
                        "shared across families — the uncorrelated ('relaxed') clock. 0 = strict "
                        "(default)")
    g.add_argument("--clock-dist", choices=("lognormal", "gamma"), default="lognormal",
                   metavar="DIST", dest="clock_dist",
                   help="lineage-clock distribution: lognormal (default, σ = log-scale) or gamma "
                        "(σ = coefficient of variation). Only used when --clock-spread > 0")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_SEQUENCE_OUTPUTS, default=None, metavar="PART",
                   help="which outputs to write (default: alignments, phylograms). also available: "
                        "ancestral (internal-node sequences), species_phylogram (the species tree "
                        "scaled by the clock)")


def _build_model(args: argparse.Namespace):
    """Build the substitution model from ``--model`` and its physical parameters (each knob falls back
    to the menu constructor's own default when not given)."""
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
    if args.clock_spread < 0:
        parser.error("--clock-spread must be >= 0")

    tree_path = os.path.join(args.genomes, "genome_species_tree.nwk")
    events_path = os.path.join(args.genomes, "genome_events.tsv")
    try:
        with open(tree_path) as f:
            tree, _ = read_newick(f.read())
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{tree_path} not found — is {args.genomes} a 'zombi2 genomes' output directory? "
            "(a genomes run writes genome_species_tree.nwk for this handoff)") from None
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
    substitution = args.substitution
    if args.clock_spread > 0:
        substitution = args.substitution * ByLineage(spread=args.clock_spread, dist=args.clock_dist)

    t0 = time.perf_counter()
    result = simulate_sequences(genome_run, model=model, length=args.length,
                                substitution=substitution, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.output, exist_ok=True)
    if args.write:
        result.write(args.output, outputs=args.write)
    else:
        result.write(args.output)               # SequencesResult.write's default: alignments + phylograms

    n_families = sum(1 for aln in result.alignments.values() if aln)
    n_seqs = sum(len(aln) for aln in result.alignments.values())
    clock = ("strict clock" if args.clock_spread == 0 else
             f"{args.clock_dist} lineage clock, spread {args.clock_spread:g}")
    summary = (f"{n_seqs} sequences across {n_families} gene families, {model.name} "
               f"{args.length} sites, {clock}")
    print(f"wrote {args.output}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(args.output, "sequences.log"), args, summary)
    return 0
