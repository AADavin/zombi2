"""zombi2 experimental command."""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from zombi2.tree import read_newick
from zombi2.cli.framework import ZombiHelpFormatter, _examples, _write_params_log

# --------------------------------------------------------------------------- #
# experimental: the zombi2.experimental layer (unstable, opt-in)
# --------------------------------------------------------------------------- #


def _add_experimental_args(p: argparse.ArgumentParser) -> None:
    """The ``experimental`` command groups unstable, not-yet-validated models (the
    ``zombi2.experimental`` layer). Currently one sub-subcommand: ``ils`` (multispecies
    coalescent)."""
    esub = p.add_subparsers(dest="experimental_command", metavar="<model>", required=True)
    sp_ils = esub.add_parser(
        "ils",
        help="incomplete lineage sorting: gene trees under the multispecies coalescent",
        description=(
            "Simulate gene trees under the multispecies coalescent, so gene lineages need not "
            "coalesce at the nodes they pass through; deep coalescence makes gene trees disagree with "
            "the containing tree -- incomplete lineage sorting (ILS). The amount of ILS is set by "
            "--population-size N, in the tree's own time units: it grows with branch_length / N.\n\n"
            "Two modes. Plain: the coalescent inside the species tree -t (single-copy orthologs). "
            "DTL + ILS: add --events-trace from a 'zombi2 genomes' run and the coalescent runs inside "
            "each gene family's locus tree (duplications/transfers/losses), one gene tree per family; "
            "a duplication's new copy, a transferred copy and the family origination are single-copy "
            "foundings (bounded coalescent), while speciations allow deep coalescence.\n\n"
            "EXPERIMENTAL: APIs and outputs may change. Pure numpy (no optional dependencies)."
        ),
        usage="zombi2 experimental ils -t FILE -N POP [--events-trace FILE] [-n R] [-k C] -o DIR",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # plain ILS: 1000 gene trees under the MSC on a species tree, one allele per species",
            "  zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 --seed 1 -o out/",
            "",
            "  # DTL + ILS: a coalescent gene tree per family from a 'genomes' run (write it with --write trace)",
            "  zombi2 experimental ils -t species_tree.nwk --events-trace run/events_trace.tsv -N 0.5 -o out/",
        ),
    )
    _add_experimental_ils_args(sp_ils)


def _add_experimental_ils_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="species-tree Newick (as written by 'zombi2 species') -- the coalescent "
                        "container, and (with --events-trace) the frame the locus trees live in")
    g.add_argument("--events-trace", default=None, metavar="FILE", dest="events_trace",
                   help="a 'genomes' run's events_trace.tsv (write it with 'zombi2 genomes ... "
                        "--write trace'). When given, run DTL + ILS: the coalescent within each gene "
                        "family's locus tree, one gene tree per family. Without it, plain species-tree ILS")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")
    g.add_argument("--seed", type=int, default=None, metavar="N", help="RNG seed for reproducibility")

    g = p.add_argument_group("coalescent")
    g.add_argument("-N", "--population-size", type=float, required=True, metavar="POP",
                   dest="population_size",
                   help="effective population size in the tree's time units: pairwise coalescence "
                        "rate 1/POP per unit time. Larger POP => more ILS (governed by branch / POP)")
    g.add_argument("-n", "--replicates", type=int, default=1, metavar="R",
                   help="independent gene trees to draw (default 1); with --events-trace, per family")
    g.add_argument("-k", "--samples", type=int, default=1, metavar="C",
                   help="alleles sampled per species tip / per extant gene copy (default 1)")


def _run_experimental_ils(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    print("zombi2: 'experimental ils' is unstable — APIs and outputs may change "
          "(zombi2.experimental).", file=sys.stderr)
    from zombi2.experimental.ils import MultispeciesCoalescent, is_concordant

    if args.replicates < 1:
        raise ValueError("--replicates must be >= 1")
    if args.samples < 1:
        raise ValueError("--samples must be >= 1")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    msc = MultispeciesCoalescent(population_size=args.population_size)
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick())

    if args.events_trace:                          # DTL + ILS: the coalescent within each locus tree
        from zombi2.genomes.reconciliation import extant_species_from_records
        from zombi2.genomes.simulation import read_events_trace
        with open(args.events_trace) as f:
            families = read_events_trace(f.read(), tree)
        if not families:
            raise ValueError(f"no gene-family events in {args.events_trace!r}")
        gid2species = extant_species_from_records(families, tree)
        fam_trees = msc._family_trees(families, gid2species, tree.total_age,
                                      args.samples, args.replicates, rng)
        gdir = os.path.join(args.out, "gene_trees")
        os.makedirs(gdir, exist_ok=True)
        for family, trees in sorted(fam_trees.items()):
            with open(os.path.join(gdir, f"{family}.nwk"), "w") as fh:
                for t in trees:
                    fh.write(t.to_newick(include_internal_names=False) + "\n")
        reps = "" if args.replicates == 1 else f", x{args.replicates} reps"
        summary = (f"experimental ils (DTL + ILS): coalescent gene trees for {len(fam_trees)} of "
                   f"{len(families)} families -> {args.out}/gene_trees/ (N={args.population_size:g}{reps})")
        print(summary)
        _write_params_log(os.path.join(args.out, "ils.log"), args, summary)
        return 0

    genes = msc.sample_gene_trees(tree, args.replicates, samples=args.samples, rng=rng)
    with open(os.path.join(args.out, "ils_gene_trees.nwk"), "w") as f:
        for g in genes:
            f.write(g.to_newick(include_internal_names=False) + "\n")
    copies = "1 copy/species" if args.samples == 1 else f"{args.samples} copies/species"
    summary = (f"experimental ils: {len(genes)} gene tree(s) under the multispecies coalescent "
               f"(N={args.population_size:g}, {copies}) -> {args.out}/gene_trees.nwk")
    if args.samples == 1:
        conc = sum(is_concordant(g, tree) for g in genes) / len(genes)
        summary += f"; {conc:.1%} match the species-tree topology"
    print(summary)
    _write_params_log(os.path.join(args.out, "ils.log"), args, summary)
    return 0


def run(args, parser):
    if args.experimental_command == "ils":
        return _run_experimental_ils(args, parser)
    parser.error(f"unknown experimental model {args.experimental_command!r}")   # unreachable
