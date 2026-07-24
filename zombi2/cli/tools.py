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
import sys

from zombi2 import tree as _tree
from zombi2.genomes.events import events_from_tsv
from zombi2.genomes.gene_trees import gene_trees_from_events
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.tree import read_newick
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
    "  tree                 transform one Newick tree (prune, round, stem, rescale, RED)\n"
    "  treedist             distance between two Newick trees (RF, branch-score)\n"
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

    trp = tsub.add_parser(
        "tree",
        help="transform one Newick tree (prune, round, stem, rescale, RED)",
        description=(
            "Apply one transform to a Newick tree and write the result (Newick to stdout, or to a file "
            "with -o). Exactly one action per call. Actions: --prune (drop dead/unsampled lineages), "
            "--round (snap a rounding-noisy dated tree to exactly ultrametric), --stem / --stem-add "
            "(set / extend the branch above the crown), --rescale-height / --rescale-factor (scale "
            "branch lengths), --red (the RED-rescaled tree; add --values for a per-node RED table). "
            "The RED-related actions and --stem/--rescale ignore tip fates, so any tree loads; --prune "
            "needs real fates (a ZOMBI tree, or an ultrametric one)."
        ),
        usage="zombi2 tools tree TREE (--prune | --round | --stem LEN | --rescale-height H | --red) [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # drop extinct lineages, to stdout",
            "  zombi2 tools tree out/species/species_complete.nwk --prune",
            "",
            "  # snap a rounding-noisy dated tree to ultrametric, to a file",
            "  zombi2 tools tree dated.nwk --round -o dated_ultrametric.nwk",
            "",
            "  # RED per node, as a table",
            "  zombi2 tools tree gtdb.nwk --red --values",
        ),
    )
    _add_tools_tree_args(trp)

    tdp = tsub.add_parser(
        "treedist",
        help="distance between two Newick trees (RF, branch-score)",
        description=(
            "Distance between two rooted Newick trees over their shared tips, printed as "
            "'<metric><TAB><value>' to stdout (or -o). --metric: rf (Robinson–Foulds), rf-normalized, "
            "branch-score (Kuhner–Felsenstein, uses branch lengths), or all. The two trees must carry "
            "the same tips, identically labelled; a mismatch is an error."
        ),
        usage="zombi2 tools treedist TREE_A TREE_B [--metric METRIC] [-o FILE]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # Robinson–Foulds between a true and an inferred tree",
            "  zombi2 tools treedist true.nwk inferred.nwk --metric rf",
            "",
            "  # every metric at once",
            "  zombi2 tools treedist true.nwk inferred.nwk --metric all",
        ),
    )
    _add_tools_treedist_args(tdp)


def _add_tools_tree_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("input", metavar="TREE", help="a Newick tree file (or - for stdin)")
    a = p.add_argument_group("action (exactly one)")
    m = a.add_mutually_exclusive_group(required=True)
    m.add_argument("--prune", action="store_true", help="drop dead/unsampled lineages → the extant tree")
    m.add_argument("--round", dest="round_", action="store_true",
                   help="snap a rounding-noisy dated tree to exactly ultrametric (tolerance --tol)")
    m.add_argument("--stem", type=float, metavar="LEN", help="set the stem (branch above the crown) to LEN")
    m.add_argument("--stem-add", type=float, metavar="LEN", dest="stem_add", help="extend the stem by LEN")
    m.add_argument("--rescale-height", type=float, metavar="H", dest="rescale_height",
                   help="scale branch lengths so root-to-tip = H")
    m.add_argument("--rescale-factor", type=float, metavar="F", dest="rescale_factor",
                   help="multiply every branch length by F")
    m.add_argument("--red", action="store_true",
                   help="the RED-rescaled tree (Relative Evolutionary Divergence on [0,1])")
    o = p.add_argument_group("options")
    o.add_argument("--tol", type=float, default=1e-3,
                   help="tolerance for --round, as a fraction of tree height (default 1e-3)")
    o.add_argument("--values", action="store_true",
                   help="with --red: write a node⇥RED table instead of a tree")
    o.add_argument("-o", "--output", metavar="FILE", help="write here instead of stdout")


def _add_tools_treedist_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("a", metavar="TREE_A", help="first Newick tree file")
    p.add_argument("b", metavar="TREE_B", help="second Newick tree file")
    p.add_argument("--metric", choices=["rf", "rf-normalized", "branch-score", "all"], default="rf",
                   help="which distance (default rf); 'all' prints every metric")
    p.add_argument("-o", "--output", metavar="FILE", help="write here instead of stdout")


def _emit(text: str, path: str | None) -> None:
    """stdout by default; a file with -o."""
    if path:
        with open(path, "w") as f:
            f.write(text.rstrip("\n") + "\n")
    else:
        print(text)


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


def _run_tree(args, parser) -> int:
    """``zombi2 tools tree`` — one transform, Newick in, Newick (or a RED table) out."""
    if args.values and not args.red:
        parser.error("--values only applies with --red")
    text = sys.stdin.read() if args.input == "-" else open(args.input).read()
    try:
        if args.prune:
            t, _ = _tree.read_newick(text)                      # prune needs real fates
            pruned = _tree.prune(t, keep="extant")
            if pruned is None:
                parser.error("no extant lineages to keep")
            out = pruned.to_newick()
        else:
            t, _ = _tree.read_newick(text, assume_extant=True)  # geometric: any tree, fates irrelevant
            if args.round_:
                out = _tree.make_ultrametric(t, tol=args.tol).to_newick()
            elif args.stem is not None:
                out = _tree.with_stem(t, args.stem).to_newick()
            elif args.stem_add is not None:
                out = _tree.with_stem(t, args.stem_add, mode="add").to_newick()
            elif args.rescale_height is not None:
                out = _tree.rescale(t, height=args.rescale_height).to_newick()
            elif args.rescale_factor is not None:
                out = _tree.rescale(t, factor=args.rescale_factor).to_newick()
            elif args.values:                                   # --red --values: the per-node table
                red = _tree.relative_evolutionary_divergence(t)
                out = "node\tRED\n" + "\n".join(f"n{i}\t{v:.6g}" for i, v in sorted(red.items()))
            else:                                               # --red: the RED-rescaled tree
                out = _tree.red_scaled(t).to_newick()
    except (ValueError, OSError) as e:
        parser.error(str(e))
    _emit(out, args.output)
    return 0


def _leaf_labels(tree, namemap: dict) -> dict:
    """``{leaf id: label}`` — the external name for an external tree, ``n<id>`` for a ZOMBI tree."""
    return {i: (namemap.get(i) or f"n{i}") for i, n in tree.nodes.items() if n.children is None}


def _relabel_leaves(tree, leaf_labels: dict, label_id: dict):
    """A copy whose LEAF ids are ``label_id[label]`` (so two trees share leaf ids **by label**);
    internal ids are shifted clear of the leaf range. Distance compares clades of leaf ids, so this
    makes treedist match tips by taxon rather than by the positionally-minted parse ids."""
    offset = len(label_id)
    new = {i: (label_id[leaf_labels[i]] if i in leaf_labels else i + offset) for i in tree.nodes}
    nodes = {new[i]: _tree.Node(new[i], None if n.parent is None else new[n.parent],
                                n.birth_time, n.end_time,
                                None if n.children is None else tuple(new[c] for c in n.children),
                                n.fate)
             for i, n in tree.nodes.items()}
    return _tree.Tree(nodes, new[tree.root])


def _run_treedist(args, parser) -> int:
    """``zombi2 tools treedist`` — a distance (or all) between two trees, to stdout. Tips are matched
    by **label** (the external name, or ``n<id>`` for a ZOMBI tree), not by parse order."""
    try:
        a, na = _tree.read_newick(open(args.a).read(), assume_extant=True)
        b, nb = _tree.read_newick(open(args.b).read(), assume_extant=True)
        la, lb = _leaf_labels(a, na), _leaf_labels(b, nb)
        sa, sb = set(la.values()), set(lb.values())
        if sa != sb:
            parser.error(f"the two trees have different leaf sets ({len(sa)} vs {len(sb)} tips, "
                         f"{len(sa ^ sb)} not shared) — treedist needs the same taxa on both")
        label_id = {lab: k for k, lab in enumerate(sorted(sa))}
        a, b = _relabel_leaves(a, la, label_id), _relabel_leaves(b, lb, label_id)
        metrics = ["rf", "rf-normalized", "branch-score"] if args.metric == "all" else [args.metric]
        lines = [f"{m}\t{_tree.distance(a, b, metric=m):g}" for m in metrics]
    except (ValueError, OSError) as e:
        parser.error(str(e))
    _emit("\n".join(lines), args.output)
    return 0


#: tool name -> handler; dispatch mirrors the level commands' ``_RUN``.
_TOOLS_RUN = {"format": _run_format, "tree": _run_tree, "treedist": _run_treedist}


def run(args, parser) -> int:
    return _TOOLS_RUN[args.tools_command](args, parser)
