"""A compact, likelihood-ready view of a (reconciled) gene tree.

ALE uses only the gene tree's **topology** and its tip → species assignment — gene-tree
branch lengths are irrelevant to the reconciliation likelihood (the DP integrates over where
on each species branch every event fell). So this structure keeps only the rooted binary
topology in post-order plus, for each tip, the species it was sampled from.

Tips are read from ZOMBI2's reconciled extant Newick, whose leaf labels are
``"<species>|<gid>"`` (see :func:`zombi2.reconciliation.reconcile`); the species is the part
before the first separator. An explicit ``tip_species`` map overrides the label parsing.
"""

from __future__ import annotations

from dataclasses import dataclass

from zombi2.tree import read_newick


@dataclass(slots=True)
class GNode:
    """One gene-tree node in post-order. Internal nodes carry child indices; leaves carry the
    species name they were sampled from."""

    is_leaf: bool
    species: str | None
    left: int | None
    right: int | None


class GeneTree:
    """A rooted binary gene tree: a post-order list of :class:`GNode`, root last."""

    def __init__(self, nodes: list[GNode]):
        self.nodes = nodes
        self.n = len(nodes)
        self.root = len(nodes) - 1  # post-order ⇒ the root is visited last

    @classmethod
    def from_newick(cls, newick: str, *, sep: str = "|",
                    tip_species: dict[str, str] | None = None) -> "GeneTree":
        """Parse a Newick string into a :class:`GeneTree`.

        ``sep`` splits ``"<species>|<gid>"`` tip labels (the ZOMBI2 reconciled-tree convention);
        pass ``sep="_"`` for ``build_gene_trees``' ``"<species>_<gid>"`` labels, or a
        ``tip_species={label: species}`` map to bypass label parsing entirely.
        """
        tree = read_newick(newick)
        nodes: list[GNode] = []
        idx: dict[int, int] = {}

        def visit(node) -> int:
            for c in node.children:
                visit(c)
            i = len(nodes)
            idx[id(node)] = i
            if not node.children:
                if node.name.startswith("LOSS"):
                    raise ValueError(
                        "gene tip labelled 'LOSS…' — ALElite consumes the EXTANT (survivors-only) "
                        "gene tree, where losses are marginalised via the extinction probability. "
                        "Pass reconciliation.extant, never .complete."
                    )
                if tip_species is not None and node.name in tip_species:
                    species = tip_species[node.name]
                else:
                    species = node.name.split(sep)[0]
                nodes.append(GNode(is_leaf=True, species=species, left=None, right=None))
            else:
                if len(node.children) != 2:
                    raise ValueError(
                        f"ALElite needs a strictly binary gene tree; a node has "
                        f"{len(node.children)} children."
                    )
                lc, rc = node.children
                nodes.append(GNode(is_leaf=False, species=None,
                                   left=idx[id(lc)], right=idx[id(rc)]))
            return i

        visit(tree.root)
        return cls(nodes)

    @classmethod
    def from_reconciliation(cls, recon, *, sep: str = "|",
                            tip_species: dict[str, str] | None = None) -> "GeneTree":
        """Build from a :class:`zombi2.reconciliation.Reconciliation`, using its **extant** tree.

        This is the intended entry point: it consumes the survivors-only tree (the observable
        gene tree ALE scores), so there is no way to accidentally feed the complete tree with
        its losses. Raises if the family left no surviving copies (``extant is None``).
        """
        extant = getattr(recon, "extant", None)
        if extant is None:
            raise ValueError(
                "reconciliation has no extant tree (the family is fully extinct); there is no "
                "observable gene tree to score."
            )
        return cls.from_newick(extant, sep=sep, tip_species=tip_species)

    def species_set(self) -> set[str]:
        return {g.species for g in self.nodes if g.is_leaf}

    def __repr__(self) -> str:
        return f"GeneTree(n={self.n}, tips={sum(g.is_leaf for g in self.nodes)})"
