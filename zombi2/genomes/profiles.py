"""Phyletic profiles — each gene family's copy count in each extant species.

The classic comparative-genomics view of a genome run: a families × extant-species matrix of copy
numbers (0 = absent). It is **derived** from the observed genomes (the extant tips of the complete
tree), so it costs nothing to keep the run lean and materialise it on access. Stored sparse — only
the nonzero ``(family, species)`` cells — with the dense matrix and presence/absence derived.
"""

from __future__ import annotations

import collections
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Profiles:
    """Gene-family copy counts across the extant species: rows are ``families``, columns are the
    extant ``species`` (node ids). Only nonzero cells are stored in ``counts``; ``matrix`` and
    ``presence`` densify on request."""

    families: tuple[int, ...]                  # row labels — families present at some extant tip, sorted
    species: tuple[int, ...]                   # column labels — extant species (node ids), sorted
    counts: dict[tuple[int, int], int]         # {(family, species): copy count}, nonzero only

    @property
    def shape(self) -> tuple[int, int]:
        """``(n_families, n_species)``."""
        return (len(self.families), len(self.species))

    @property
    def matrix(self) -> np.ndarray:
        """Dense ``families × species`` copy-count matrix (``int``; 0 where a family is absent)."""
        fi = {f: i for i, f in enumerate(self.families)}
        si = {s: j for j, s in enumerate(self.species)}
        m = np.zeros(self.shape, dtype=int)
        for (f, s), c in self.counts.items():
            m[fi[f], si[s]] = c
        return m

    @property
    def presence(self) -> np.ndarray:
        """Dense ``families × species`` presence/absence matrix (``0``/``1``)."""
        return (self.matrix > 0).astype(int)

    def to_tsv(self, *, presence: bool = False) -> str:
        """The matrix as TSV — a ``family`` column then one column per extant species (``n<id>``),
        one row per family. ``presence=True`` writes 0/1 instead of copy counts."""
        m = self.presence if presence else self.matrix
        header = "family\t" + "\t".join(f"n{s}" for s in self.species)
        rows = [f"{f}\t" + "\t".join(str(v) for v in m[i]) for i, f in enumerate(self.families)]
        return "\n".join([header, *rows]) + "\n"


def profiles_from_genomes(genomes: dict, extant_ids) -> Profiles:
    """Build the phyletic profiles from the per-node ``genomes``, over the extant species
    ``extant_ids`` — tally each family's copy count at each extant tip."""
    species = tuple(sorted(extant_ids))
    counts: dict[tuple[int, int], int] = {}
    families: set[int] = set()
    for s in species:
        for family, c in collections.Counter(g.family for g in genomes[s]).items():
            counts[(family, s)] = c
            families.add(family)
    return Profiles(families=tuple(sorted(families)), species=species, counts=counts)
