"""``zombi2 tools`` — analyses that read a finished run and derive a new view of it.

Where the level commands *simulate*, the tools *read back* what a run wrote. Each tool is its own
sub-subcommand (``zombi2 tools <tool>``); the first is ``format``, which turns a genomes run into
analysis-ready tables. Its first table is the **homology** matrix: for each gene family, an n×n grid
(n the extant leaves) of ``O`` / ``P`` / ``X`` — ortholog, paralog, xenolog — read from the event at
each gene pair's most-recent common ancestor. Exact, because ZOMBI simulated the gene tree's embedding
rather than inferring it (see :mod:`zombi2.tools.homology`).
"""
from __future__ import annotations

import argparse
import os

from zombi2.genomes.events import events_from_tsv
from zombi2.genomes.gene_trees import gene_trees_from_events
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.species import read_newick
from zombi2.tools.homology import write_homology
from zombi2.cli.framework import (
    ZombiHelpFormatter, _add_flat_arg, _add_from_arg, _add_run_arg, _examples, level_dir,
    resolve_genomes,
)

#: the tables ``format`` can emit — ``name -> (subdirectory, writer, one-line gloss)``. One today; the
#: menu is declared so a second table is one entry plus its writer, and the ``--format`` help is built
#: from it, so it can never advertise a table that is not wired.
_FORMATS = {
    "homology": ("homology", write_homology,
                 "per-family n×n O/P/X table (ortholog / paralog / xenolog)"),
}

#: the tools description carries its own tool list (the house-style formatter hides argparse's auto
#: subcommand dump, exactly as the top-level help does), so ``zombi2 tools -h`` still names them.
_TOOLS_DESCRIPTION = (
    "Analyses that read a finished run and derive a new view of it. Run 'zombi2 tools <tool> -h' "
    "for a tool's options.\n\n"
    "Tools\n"
    "  format               turn a genomes run into analysis-ready tables (homology O/P/X, …)\n"
)


def _add_tools_args(p: argparse.ArgumentParser) -> None:
    """Wire the ``tools`` sub-subcommands. Adding a tool is a new ``add_parser`` here, its own
    ``_add_tools_<tool>_args`` builder, and one ``_TOOLS_RUN`` entry — the same three-touch shape the
    level commands have."""
    tsub = p.add_subparsers(dest="tools_command", metavar="<tool>", required=True)
    fp = tsub.add_parser(
        "format",
        help="turn a genomes run into analysis-ready tables (e.g. the homology matrix)",
        description=(
            "Read a finished 'zombi2 genomes' run and write analysis-ready tables derived from its "
            "gene trees. Today one --format is offered: 'homology' — for each family, an n×n table "
            "(n the extant leaves) of O/P/X (ortholog / paralog / xenolog), read from the event at "
            "each pair's most-recent common ancestor. Exact, not inferred: ZOMBI recorded the "
            "embedding as it simulated it. Works for every resolution; on a nucleotide run it is one "
            "table per declared gene (the intergenic spacer is not a gene and gets none). Tables land "
            "in the run's genomes/homology/."
        ),
        usage="zombi2 tools format DIR [--from PATH] [--format FORMAT ...] [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # O/P/X homology tables for a genomes run, written to its genomes/homology/",
            "  zombi2 tools format out/",
            "",
            "  # read a run that lives elsewhere",
            "  zombi2 tools format out/ --from other_run/",
        ),
    )
    _add_tools_format_args(fp)


def _add_tools_format_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "the genomes run whose gene trees the tables are derived from")
    g = p.add_argument_group("general")
    _add_from_arg(g, "the genomes run to read — its species tree and genome_events.tsv rebuild the "
                     "gene trees")
    g = p.add_argument_group("outputs")
    g.add_argument(
        "--format", nargs="+", choices=sorted(_FORMATS), default=["homology"], metavar="FORMAT",
        dest="formats",
        help="which tables to write (default: homology). " +
             "  ".join(f"{name}: {gloss}" for name, (_, _, gloss) in sorted(_FORMATS.items())))
    _add_flat_arg(g)


def _load_gene_trees(handoff, tree):
    """The run's ``{family: GeneTree}``, rebuilt from disk for either genome resolution.

    An unordered or ordered run derives its gene trees from the event log alone. A **nucleotide** run
    keys its events by ancestral interval, so its trees are recovered from the genome itself (the same
    ``read_nucleotide_genomes`` the sequence level replays) — and only its **declared genes** get a
    tree, never the intergenic spacer, which is what ``.gene_trees`` gives (the spacer's blocks live in
    ``.block_trees``). A nucleotide run that declared no genes is one long intergene, so there is
    nothing to relate."""
    if os.path.exists(os.path.join(handoff, "blocks.tsv")):     # the nucleotide resolution's mark
        genome_run = read_nucleotide_genomes(handoff, tree)
        if not genome_run.gene_spans:
            raise ValueError(
                "this nucleotide run declared no genes — it is one uninterrupted intergene, so there "
                "is nothing to relate. Re-run 'zombi2 genomes --resolution nucleotide' with --genes "
                "or --gff to lay down genes.")
        return genome_run.gene_trees                            # declared genes only; spacer excluded
    events_path = os.path.join(handoff, "genome_events.tsv")
    try:
        with open(events_path) as f:
            events = events_from_tsv(f.read())
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{events_path} not found — re-run 'zombi2 genomes' with 'events' in --write so the gene "
            "genealogy can be rebuilt") from None
    return gene_trees_from_events(events, tree)


def _run_format(args, parser) -> int:
    """``zombi2 tools format`` — rebuild the run's gene trees and write the requested tables."""
    handoff, tree_path = resolve_genomes(args.source or args.run)
    with open(tree_path) as f:
        tree, _ = read_newick(f.read())
    gene_trees = _load_gene_trees(handoff, tree)

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "genomes", args.flat)
    n_tables = sum(1 for gt in gene_trees.values() if gt.extant is not None)
    wrote = []
    for name in dict.fromkeys(args.formats):            # de-dupe, keep the order given
        subdir, writer, _ = _FORMATS[name]
        directory = level_dir(out, subdir, args.flat)
        writer(gene_trees, directory)
        wrote.append(f"{name} → {os.path.relpath(directory, args.run)}/")
    print(f"wrote {args.run}/ ({n_tables} table(s) per format: {', '.join(wrote)})")
    return 0


#: tool name -> handler; dispatch mirrors the level commands' ``_RUN``.
_TOOLS_RUN = {"format": _run_format}


def run(args, parser) -> int:
    return _TOOLS_RUN[args.tools_command](args, parser)
