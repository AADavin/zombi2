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


def _read_text(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


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
    def from_tsv(cls, source: str) -> "ProfileMatrix":
        """Load a copy-number matrix written by :meth:`to_tsv`.

        ``source`` is either a path to a TSV file or the raw TSV text (anything
        containing a newline is treated as text). The header's first column label is
        ignored; the remaining labels are species, and each subsequent row is a family.
        """
        text = source if "\n" in source else _read_text(source)
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            raise ValueError("empty profile table")
        species = lines[0].split("\t")[1:]
        families: list[str] = []
        rows: list[list[int]] = []
        for line in lines[1:]:
            parts = line.split("\t")
            families.append(parts[0])
            rows.append([int(x) for x in parts[1:]])
        matrix = (np.array(rows, dtype=int) if rows
                  else np.zeros((0, len(species)), dtype=int))
        if matrix.shape[1:] and matrix.shape[1] != len(species):
            raise ValueError("row width does not match the number of species columns")
        return cls(families=families, species=species, matrix=matrix)

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
