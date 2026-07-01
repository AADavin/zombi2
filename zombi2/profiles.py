"""The phylogenetic profile matrix — families x extant species.

This is the key v1 output and the σ-sample dataset the future inverse-Potts / DCA
validation will consume. It is read directly off the final genome state at each extant
leaf; no gene trees are needed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from .genome import Genome
from .tree import TreeNode


def _natkey(name: str) -> tuple[int, str]:
    """Natural-ish sort key: order by the numeric run in the name, then the string."""
    digits = re.sub(r"\D", "", name)
    return (int(digits) if digits else 0, name)


@dataclass
class ProfileMatrix:
    """Copy-number matrix with row (family) and column (species) labels."""

    families: list[str]
    species: list[str]
    matrix: np.ndarray  # shape (n_families, n_species), integer copy numbers

    def presence(self) -> np.ndarray:
        """Binary presence/absence matrix (copy number > 0)."""
        return (self.matrix > 0).astype(np.int8)

    def to_tsv(self, presence: bool = False) -> str:
        data = self.presence() if presence else self.matrix
        lines = ["family\t" + "\t".join(self.species)]
        for i, family in enumerate(self.families):
            lines.append(family + "\t" + "\t".join(str(int(x)) for x in data[i]))
        return "\n".join(lines) + "\n"

    @classmethod
    def from_leaf_genomes(cls, leaf_genomes: dict[TreeNode, Genome]) -> "ProfileMatrix":
        species_nodes = sorted(leaf_genomes.keys(), key=lambda n: _natkey(n.name))
        species = [n.name for n in species_nodes]

        famset: set[str] = set()
        for genome in leaf_genomes.values():
            famset.update(genome.families())
        families = sorted(famset, key=_natkey)

        matrix = np.zeros((len(families), len(species)), dtype=int)
        fidx = {f: i for i, f in enumerate(families)}
        for j, node in enumerate(species_nodes):
            genome = leaf_genomes[node]
            for family in genome.families():
                matrix[fidx[family], j] = genome.copy_number(family)

        return cls(families=families, species=species, matrix=matrix)
