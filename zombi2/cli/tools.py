"""``zombi2 tools`` — analyses that read a finished run and derive a new view of it.

Where the level commands *simulate*, the tools *read back* what a run wrote. Each tool is its own
sub-subcommand (``zombi2 tools <tool>``); the first is ``format``, which turns a genomes run into
analysis-ready files, all derived from the gene trees and all exact rather than inferred, because
ZOMBI simulated the embedding it is reporting: the **homology** matrix — for each family an n×n grid
(n the extant leaves) of ``S`` / ``D`` / ``T`` (+ ``x``), the event at each pair's common ancestor and
whether transfer is in their history since (`zombi2.tools.homology`); the **marker table**, one row
per family saying whether it can be trusted to recover the species tree, which is what to read if you
came looking for orthologs (`zombi2.tools.markers`); and **recPhyloXML**, each family's complete gene
tree written inside the complete species tree in the community format for that
(`zombi2.tools.recphylo`).
"""
from __future__ import annotations

from collections import Counter

import argparse
import os
import re
import sys

from zombi2 import tree as _tree
from zombi2.genomes.events import events_from_tsv
from zombi2.genomes.gene_trees import gene_trees_from_events
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.tree import read_newick
from zombi2.tools.homology import write_homology
from zombi2.tools.markers import write_markers
from zombi2.tools.recphylo import write_recphylo
from zombi2.cli.framework import (
    ZombiHelpFormatter, _add_flat_arg, _add_from_arg, _add_quiet_arg, _add_run_arg, _examples,
    level_dir, resolve_genomes, warn,
)

#: what ``format`` can emit — ``name -> (subdirectory, writer, one-line gloss)``. The menu is
#: declared so a new one is an entry plus its writer, and the ``--format`` help is built from it, so
#: it can never advertise something that is not wired. A writer takes ``(gene_trees, species tree,
#: directory)`` and gives back a short description of what it wrote, for the summary line.
_FORMATS = {
    "homology": ("homology", write_homology,
                 "per-family n×n table: S/D/T for the event at each pair's ancestor, x if a "
                 "transfer came after"),
    "recphylo": ("recphylo", write_recphylo,
                 "per-family recPhyloXML — the gene tree drawn inside the species tree"),
    "markers": ("markers", write_markers,
                "one row per family: single-copy, universal, and does its tree match the species tree"),
}

#: the tools description carries its own tool list (the house-style formatter hides argparse's auto
#: subcommand dump, exactly as the top-level help does), so ``zombi2 tools -h`` still names them.
_TOOLS_DESCRIPTION = (
    "Analyses that read a finished run and derive a new view of it. Run 'zombi2 tools <tool> -h' "
    "for a tool's options.\n\n"
    "Tools\n"
    "  format               turn a genomes run into analysis-ready files (homology, recPhyloXML)\n"
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
        help="turn a genomes run into analysis-ready files (homology, markers, recPhyloXML)",
        description=(
            "Read a finished 'zombi2 genomes' run and write analysis-ready files derived from its "
            "gene trees. Three --format choices: 'homology' — for each family, an n×n table (n the "
            "extant leaves) giving the event at each pair's common ancestor (S speciation, D "
            "duplication, T transfer) and whether transfer is in their history since (an x suffix); "
            "'markers' — one row per family: is it single-copy, is it universal, and does its true "
            "tree match the species tree, which is what to read if you are after genes to build a "
            "species tree from — and 'recphylo' — each family's complete gene tree "
            "written inside the complete species tree as recPhyloXML, the community format for that, "
            "ready for a viewer or for scoring a reconciliation method against. All three are exact, "
            "not inferred: ZOMBI recorded the embedding as it simulated it. All three work at every "
            "resolution; on a nucleotide run there is one per declared gene (the intergenic spacer is "
            "not a gene and gets none). Files land in the run's genomes/homology/, genomes/markers/ "
            "and genomes/recphylo/."
        ),
        usage="zombi2 tools format DIR [--from PATH] [--format FORMAT ...] [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # which families make trustworthy species-tree markers",
            "  zombi2 tools format out/ --format markers",
            "",
            "  # the per-pair table: the event at each pair's common ancestor",
            "  zombi2 tools format out/",
            "",
            "  # recPhyloXML instead — one file per family, for a viewer",
            "  zombi2 tools format out/ --format recphylo",
            "",
            "  # both, from a run that lives elsewhere",
            "  zombi2 tools format out/ --from other_run/ --format homology recphylo",
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
        with open(path, "w", encoding="utf-8") as f:
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
             "; ".join(f"{name}: {gloss}" for name, (_, _, gloss) in sorted(_FORMATS.items())))
    g.add_argument("--recphylo", choices=("complete", "extant", "both"), default="complete",
                   metavar="SCOPE",
                   help="which history --format recphylo writes (default complete). complete: the "
                        "whole simulated history inside the complete species tree. extant: it "
                        "projected onto what a dataset can hold — the extant gene tree inside the "
                        "extant species tree — written twice, 'true' (rooted where the family really "
                        "began, the answer key for ancestral gene content) and 'recoverable' (rooted "
                        "at the surviving copies' ancestor, the most any method could recover), plus "
                        "family_origins.tsv saying how each family entered. both: all of them")
    _add_flat_arg(g)
    _add_quiet_arg(g)


def _load_gene_trees(handoff, tree):
    """The run's ``{family: GeneTree}``, rebuilt from disk for either genome resolution.

    A family or ordered run derives its gene trees from the event log alone. A **nucleotide** run
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
        with open(events_path, encoding="utf-8") as f:
            events = events_from_tsv(f.read())
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{events_path} not found — re-run 'zombi2 genomes' with 'events' in --write so the gene "
            "genealogy can be rebuilt") from None
    return gene_trees_from_events(events, tree)


def _run_format(args, parser) -> int:
    """``zombi2 tools format`` — rebuild the run's gene trees and write the requested tables."""
    handoff, tree_path = resolve_genomes(args.source or args.run)
    with open(tree_path, encoding="utf-8") as f:
        tree, _ = read_newick(f.read())
    gene_trees = _load_gene_trees(handoff, tree)

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "genomes", args.flat)
    wrote = []
    for name in dict.fromkeys(args.formats):            # de-dupe, keep the order given
        subdir, writer, _ = _FORMATS[name]
        directory = level_dir(out, subdir, args.flat)
        # every writer takes (gene_trees, species tree, directory); recphylo alone has a choice to
        # make, so it is the one that takes an option too
        extra = {"scope": args.recphylo} if name == "recphylo" else {}
        what = writer(gene_trees, tree, directory, **extra)   # each writer says what it wrote
        wrote.append(f"{name}: {what} → {os.path.relpath(directory, args.run)}/")
    # unconditional, like every other command's completion line: --quiet takes away the progress bar,
    # not the one line saying what landed and where. This alone printed nothing under --quiet, so a
    # scripted run had to go and find the files to learn whether it had written any.
    print(f"wrote {args.run}/ ({'; '.join(wrote)})")
    return 0


def _run_tree(args, parser) -> int:
    """``zombi2 tools tree`` — one transform, Newick in, Newick (or a RED table) out."""
    if args.values and not args.red:
        parser.error("--values only applies with --red")
    text = sys.stdin.read() if args.input == "-" else open(args.input, encoding="utf-8").read()
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
                out = "node\tRED\n" + "\n".join(f"{_tree.node_label(i)}\t{v:.6g}"
                                                  for i, v in sorted(red.items()))
            else:                                               # --red: the RED-rescaled tree
                out = _tree.red_scaled(t).to_newick()
    except (ValueError, OSError) as e:
        parser.error(str(e))
    _emit(out, args.output)
    return 0


def _leaf_labels(tree, namemap: dict) -> dict:
    """``{leaf id: label}`` — the external name for an external tree, ``n<id>`` for a ZOMBI tree."""
    return {i: (namemap.get(i) or _tree.node_label(i))
            for i, n in tree.nodes.items() if n.children is None}


#: a ZOMBI gene-copy leaf, ``n<species>_g<copy>`` — the label every gene tree, alignment record and
#: homology header uses. What makes it safe to detect is that both halves are ``n``/``g`` + digits,
#: which a species label never is and a real taxon name (``E_coli``, ``Nostoc_sp_PCC7120``) never is.
_GENE_LEAF = re.compile(r"^(n\d+)_g\d+$")


def _species_behind(labels: dict) -> dict | None:
    """``{leaf id: its species label}`` when *every* leaf is a gene copy, else ``None``.

    A gene tree's tips are genes and a species tree's are species, so the two are not the same kind
    of thing and comparing them by label gives an empty intersection. They can still be compared —
    on the species each gene sits in — but only when the reader knows that is what happened, which is
    what this detection is for."""
    out = {}
    for leaf, label in labels.items():
        m = _GENE_LEAF.match(label)
        if m is None:
            return None
        out[leaf] = m.group(1)
    return out


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


def _match_gene_tree_to_species_tree(la: dict, lb: dict, parser) -> tuple[dict, dict]:
    """When exactly one of the two trees has gene-copy tips, relabel it by the species each gene sits
    in — and say so. Both trees the same kind ⇒ nothing to do.

    A gene tree and a species tree are not the same kind of object, and left alone they share no
    labels at all, so the comparison would fail as "different leaf sets" without saying why. They are
    comparable on the species, and often that is exactly the question (does this family's tree recover
    the species tree?) — but only when the mapping is one gene per species. A family with two copies
    somewhere has no well-defined answer, and a plausible number would be worse than a refusal."""
    ga, gb = _species_behind(la), _species_behind(lb)
    if (ga is None) == (gb is None):            # both gene trees, or both species trees
        return la, lb
    which, other = ("first", "second") if ga is not None else ("second", "first")
    species = ga if ga is not None else gb
    repeated = sorted({s for s, n in Counter(species.values()).items() if n > 1})
    if repeated:
        parser.error(
            f"the {which} tree's tips are gene copies (n<species>_g<copy>) and the {other} tree's are "
            f"species, so the two can only be compared on the species each gene sits in — but "
            f"{', '.join(repeated)} carry several copies each, so that is not one gene per species "
            f"and the distance is not defined. Compare two gene trees instead, or reduce the family "
            f"to one copy per species first.")
    warn(f"the {which} tree's tips are gene copies and the {other} tree's are species; comparing "
         f"them on the species each gene sits in. This family is single-copy, so that is one gene "
         f"per species — the distance is a gene tree scored against a species tree, not two trees "
         f"of the same kind.")
    return (species, lb) if ga is not None else (la, species)


def _run_treedist(args, parser) -> int:
    """``zombi2 tools treedist`` — a distance (or all) between two trees, to stdout. Tips are matched
    by **label** (the external name, or ``n<id>`` for a ZOMBI tree), not by parse order."""
    try:
        a, na = _tree.read_newick(open(args.a, encoding="utf-8").read(), assume_extant=True)
        b, nb = _tree.read_newick(open(args.b, encoding="utf-8").read(), assume_extant=True)
        la, lb = _leaf_labels(a, na), _leaf_labels(b, nb)
        la, lb = _match_gene_tree_to_species_tree(la, lb, parser)
        # Uniqueness first, and separately from the leaf-set check below — that one compares *sets*,
        # so a repeated label collapses into it and passes. The relabelling then maps both copies to
        # one id, and the distance comes back as a plausible number computed on a tree that is not a
        # tree. A wrong number from a scoring tool is worse than a refusal.
        for which, labels in (("first", la), ("second", lb)):
            repeated = sorted({lab for lab, n in Counter(labels.values()).items() if n > 1})
            if repeated:
                parser.error(
                    f"the {which} tree repeats the tip label(s) {', '.join(repeated)} — a distance "
                    f"between trees is only defined when each taxon appears once, so this cannot be "
                    f"scored. If these are gene-tree tips, several copies of a family in one genome "
                    f"share a species: compare gene trees to each other, or pick one copy per species "
                    f"first.")
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
