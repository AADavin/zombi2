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

from ..genomes.events import gene_label
from ..genomes.gene_trees import GeneTree
from ..tree import prune
from .reconciliation import extant_reconciliation, origins_tsv, visible_branches

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
    name = tree.labels()
    stack: list[tuple[int | None, int]] = [(tree.root, indent)]
    while stack:
        node_id, depth = stack.pop()
        pad = "  " * depth
        if node_id is None:                                  # a close marker pushed below
            out.append(f"{pad}</clade>")
            continue
        node = tree.nodes[node_id]
        out.append(f"{pad}<clade>")
        out.append(f"{pad}  <name>{name[node_id]}</name>")
        stack.append((None, depth))
        for child in reversed(node.children or ()):          # reversed: the stack restores the order
            stack.append((child, depth + 1))
    return out


def _gene_lines(root, indent: int, name) -> list[str]:
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
        out.append(f"{pad}  <name>{_escape(name[node.species])}_{gene_label(node.copy)}</name>")
        out.append(f"{pad}  <eventsRec>")
        if arrived_in is not None:   # the format's two-step transfer: where it landed, then what it did
            out.append(f"{pad}    <transferBack destinationSpecies=\"{name[arrived_in]}\">"
                       f"</transferBack>")
        out.append(f"{pad}    <{tag} speciesLocation=\"{name[node.species]}\"></{tag}>")
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
        lines += _gene_lines(gt.complete, 3, tree.labels())
        lines += ["    </phylogeny>", "  </recGeneTree>"]
    lines.append("</recPhylo>")
    return "\n".join(lines) + "\n"


def reconciliation_xml(root, extant_tree, fam: int) -> str:
    """One ``<recPhylo>`` for an **extant-only** reconciliation: the extant species tree, then the
    projected gene tree. Same serialiser as the complete case — the projection hands back an ordinary
    `GeneNode` tree whose species are extant-tree branches, so nothing here is special-cased."""
    lines = ["<recPhylo>", "  <spTree>", "    <phylogeny>"]
    lines += _species_lines(extant_tree, 3)
    lines += ["    </phylogeny>", "  </spTree>",
              "  <recGeneTree>", "    <phylogeny rooted=\"true\">", f"      <name>fam{fam}</name>"]
    lines += _gene_lines(root, 3, extant_tree.labels())
    lines += ["    </phylogeny>", "  </recGeneTree>", "</recPhylo>"]
    return "\n".join(lines) + "\n"


def write_recphylo(gene_trees: dict[int, GeneTree], tree, directory, *,
                   scope: str = "complete") -> str:
    """Write recPhyloXML into ``directory``, one family per file, as every other per-family output of
    a run does — so a family can be handed to a viewer on its own.

    ``scope`` picks which history:

    - ``"complete"`` → ``recphylo_fam<family>.xml``, the whole simulated history inside the complete
      species tree. Written for **every** family, extinct ones included: a family that left no
      survivor still has a history worth looking at.
    - ``"extant"`` → ``recphylo_fam<family>_{true,recoverable}.xml`` plus ``family_origins.tsv``, the
      history projected onto what a dataset can hold. Only families with a surviving copy, because
      the others are not in a dataset at all. See `zombi2.tools.reconciliation` for the two scopes and
      why both are written.
    - ``"both"`` → all of the above.
    """
    if scope not in ("complete", "extant", "both"):
        raise ValueError(f"scope is 'complete', 'extant' or 'both', got {scope!r}")
    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    wrote = []
    if scope in ("complete", "both"):
        for fam, gt in sorted(gene_trees.items()):
            (d / f"recphylo_fam{fam}.xml").write_text(recphylo_xml({fam: gt}, tree),
                                                      encoding="utf-8")
        wrote.append(f"{len(gene_trees)} file(s)")
    if scope in ("extant", "both"):
        extant_tree = prune(tree, keep="extant")
        if extant_tree is None:
            raise ValueError("no extant lineages, so there is no extant tree to reconcile against")
        visible = visible_branches(tree)
        origins: dict[int, object] = {}
        n = 0
        for fam, gt in sorted(gene_trees.items()):
            for which in ("true", "recoverable"):
                rec = extant_reconciliation(gt, tree, scope=which, visible=visible)
                if rec is None:                     # no surviving copy: not in a dataset at all
                    break
                (d / f"recphylo_fam{fam}_{which}.xml").write_text(
                    reconciliation_xml(rec.root, extant_tree, fam), encoding="utf-8")
                if which == "true":
                    origins[fam] = rec
                    n += 1
        (d / "family_origins.tsv").write_text(origins_tsv(origins), encoding="utf-8")
        wrote.append(f"{n} extant pair(s) + origins")
    return ", ".join(wrote)


__all__ = ["recphylo_xml", "reconciliation_xml", "write_recphylo"]
