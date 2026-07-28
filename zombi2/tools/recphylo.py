"""recPhyloXML — a family's gene tree written *inside* the species tree it evolved in.

recPhyloXML (Duchemin et al. 2018) is the community format for a gene tree embedded in a species
tree: each gene-tree node carries the species branch it sat on and the event that ended it, so a
viewer can draw one inside the other and the whole history reads at a glance. It is normally the
output of an *inference* — an embedding recovered from an observed gene tree, which is what the word
reconciliation means. Nothing is recovered here: ZOMBI simulated the embedding, so what is written
is the true history, and it can be used as the answer key an inference is scored against.

The **complete** gene tree is what goes in, inside the **complete** species tree, because the events
this format exists to show are the losses — a gene that died leaves nothing in the extant tree to
put a `<loss>` on. Extinct and unsampled species are in the species tree for the same reason: a
transfer can arrive from a lineage that later died, and the edge has to land somewhere.

The mapping from ZOMBI2's event log is direct, one gene-tree node to one `<clade>`:

===============  ======================================================================
duplication      ``<duplication speciesLocation="n<species>">``
speciation       ``<speciation speciesLocation="n<species>">`` — the *parent* species, the
                 branch the gene was on when its species split
loss             ``<loss speciesLocation="n<species>">``
extant tip       ``<leaf speciesLocation="n<species>">``
extinct tip      ``<leaf …>`` as well — the gene reached the end of a species branch that
unsampled tip    happened to die; the species tree says which fate that branch had
transfer         two tags, the format's own two-step: ``<branchingOut speciesLocation=
                 "n<donor>">`` on the node the copy left from, and ``<transferBack
                 destinationSpecies="n<recipient>">`` opening the child that arrived
===============  ======================================================================

Origination has no tag, and needs none: a family founded mid-branch is simply a gene tree whose root
starts there. Branch lengths are left out — no reference file carries them and no viewer reads them;
the dated trees are next door in `gene_trees/` and `species_complete.nwk`.
"""

from __future__ import annotations

import pathlib

from ..genomes.events import copy_label, node_label
from ..genomes.gene_trees import GeneTree

#: gene-node kind → the recPhyloXML tag that ends its ``<eventsRec>``. Every kind a
#: `GeneNode` can carry is here, so an unmapped one is a real gap
#: rather than a silently dropped event.
_TAG = {"duplication": "duplication", "speciation": "speciation", "transfer": "branchingOut",
        "loss": "loss", "extant": "leaf", "extinct": "leaf", "unsampled": "leaf"}


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _species_lines(tree, indent: int) -> list[str]:
    """The complete species tree as nested ``<clade><name>n<id></name>`` — every node named, because
    that name is what every ``speciesLocation`` in the gene trees points at."""
    out: list[str] = []
    stack: list[tuple[int | None, int]] = [(tree.root, indent)]
    while stack:
        node_id, depth = stack.pop()
        pad = "  " * depth
        if node_id is None:                                  # a close marker pushed below
            out.append(f"{pad}</clade>")
            continue
        node = tree.nodes[node_id]
        out.append(f"{pad}<clade>")
        out.append(f"{pad}  <name>{node_label(node_id)}</name>")
        stack.append((None, depth))
        for child in reversed(node.children or ()):          # reversed: the stack restores the order
            stack.append((child, depth + 1))
    return out


def _gene_lines(root, indent: int) -> list[str]:
    """One gene tree as nested ``<clade>``, each with the ``<eventsRec>`` its node's event calls for.

    Iterative because gene-tree depth is unbounded — a duplication ladder runs past the interpreter's
    recursion limit, which is why every other walk over these trees is iterative too."""
    out: list[str] = []
    # (node, the species it was transferred INTO or None, depth); None in the node slot closes a clade
    stack: list[tuple[object, int | None, int]] = [(root, None, indent)]
    while stack:
        node, arrived_in, depth = stack.pop()
        pad = "  " * depth
        if node is None:
            out.append(f"{pad}</clade>")
            continue
        tag = _TAG[node.kind]
        out.append(f"{pad}<clade>")
        out.append(f"{pad}  <name>{_escape(copy_label(node.species, node.copy))}</name>")
        out.append(f"{pad}  <eventsRec>")
        if arrived_in is not None:   # the format's two-step transfer: where it landed, then what it did
            out.append(f"{pad}    <transferBack destinationSpecies=\"{node_label(arrived_in)}\">"
                       f"</transferBack>")
        out.append(f"{pad}    <{tag} speciesLocation=\"{node_label(node.species)}\"></{tag}>")
        out.append(f"{pad}  </eventsRec>")
        stack.append((None, None, depth))
        # A transfer's children are the copy that stayed and the copy that left, and the one that
        # left is the one now on another branch. A self-transfer puts both on the same branch, and
        # then it is the second child that arrived — the engine logs the donor's continuation first.
        moved = None
        if node.kind == "transfer" and node.children:
            moved = next((c for c in node.children if c.species != node.species), node.children[-1])
        for child in reversed(node.children):
            stack.append((child, child.species if child is moved else None, depth + 1))
    return out


def recphylo_xml(gene_trees: dict[int, GeneTree], tree) -> str:
    """One ``<recPhylo>`` document: the species tree once, then one ``<recGeneTree>`` per family.

    Handing it every family gives a single file a viewer can draw all of them in at once; handing it
    one gives that family's own file, which is what `write_recphylo()` writes."""
    lines = ["<recPhylo>", "  <spTree>", "    <phylogeny>"]
    lines += _species_lines(tree, 3)
    lines += ["    </phylogeny>", "  </spTree>"]
    for fam, gt in sorted(gene_trees.items()):
        # the family id as the phylogeny's own name, so a file holding several says which is which
        lines += ["  <recGeneTree>", "    <phylogeny rooted=\"true\">", f"      <name>fam{fam}</name>"]
        lines += _gene_lines(gt.complete, 3)
        lines += ["    </phylogeny>", "  </recGeneTree>"]
    lines.append("</recPhylo>")
    return "\n".join(lines) + "\n"


def write_recphylo(gene_trees: dict[int, GeneTree], tree, directory) -> str:
    """Write ``recphylo_fam<family>.xml`` — one family per file — into ``directory``.

    One file per family, like every other per-family output of a run, so a family can be handed to a
    viewer on its own. Unlike the homology tables these are written for **every** family, extinct
    ones included: the complete tree is the point, and a family that left no survivor still has a
    history worth looking at. For all of them in one document instead, call `recphylo_xml()` with
    the whole dictionary."""
    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    for fam, gt in sorted(gene_trees.items()):
        (d / f"recphylo_fam{fam}.xml").write_text(recphylo_xml({fam: gt}, tree), encoding="utf-8")
    return f"{len(gene_trees)} file(s)"


__all__ = ["recphylo_xml", "write_recphylo"]
