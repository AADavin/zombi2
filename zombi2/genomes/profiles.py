"""The phylogenetic profile matrix — families x extant species.

This is the key v1 output: a presence/copy-number panel over the extant leaves. It is read
directly off the final genome state at each extant leaf; no gene trees are needed.

**Sparse by construction.** A copy-number profile is overwhelmingly zero — each gene
family is present in only a handful of species — so the dense ``families x species``
array is O(N²) in the number of tips and becomes the memory wall long before the
simulation does. :class:`ProfileMatrix` therefore stores the data in **COO** form (three
parallel arrays: family index, species index, copy number, one entry per *present* cell),
which is O(non-zeros) ≈ O(N). The dense array is still available on demand via
:attr:`~ProfileMatrix.matrix` (it densifies when you ask), and every summary the library
computes — presence, genome sizes, frequency spectrum — is done straight off the sparse
arrays without ever materialising it.
"""

from __future__ import annotations

import re

import numpy as np

from zombi2.genomes.genome import Genome
from zombi2.tree import TreeNode


def _read_text(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _natkey(name: str) -> tuple[int, str]:
    """Natural-ish sort key: order by the numeric run in the name, then the string."""
    digits = re.sub(r"\D", "", name)
    return (int(digits) if digits else 0, name)


class ProfileMatrix:
    """Sparse copy-number profile with row (family) and column (species) labels.

    Stored as COO: ``rows`` (family indices), ``cols`` (species indices) and ``data``
    (positive copy numbers), one entry per present cell. Construct it from a dense array
    (``ProfileMatrix(families, species, dense)`` — the entry point stays backward
    compatible) or directly from COO arrays (``coo=(rows, cols, data)``), the path the
    engines use so a million-tip run never builds the dense matrix.
    """

    __slots__ = ("families", "species", "_rows", "_cols", "_data", "_shape")

    def __init__(self, families, species, matrix=None, *, coo=None):
        self.families = list(families)
        self.species = list(species)
        self._shape = (len(self.families), len(self.species))
        if coo is not None:
            rows = np.asarray(coo[0], dtype=np.int64)
            cols = np.asarray(coo[1], dtype=np.int64)
            data = np.asarray(coo[2], dtype=np.int64)
            if data.shape[0]:
                # coalesce duplicate (family, species) cells by summing — COO permits repeats, and
                # without this the dense view (last-write-wins) would disagree with the bincount
                # reductions (which sum). Only rewrites when duplicates exist, so a normal run
                # (no repeats) keeps its exact cell order and stays byte-identical.
                flat = rows * self._shape[1] + cols
                uniq, inv = np.unique(flat, return_inverse=True)
                if uniq.shape[0] != flat.shape[0]:
                    data = np.bincount(inv.ravel(), weights=data).astype(np.int64)
                    rows = (uniq // self._shape[1]).astype(np.int64)
                    cols = (uniq % self._shape[1]).astype(np.int64)
            self._rows, self._cols, self._data = rows, cols, data
        elif matrix is not None:
            m = np.asarray(matrix)
            if m.size:
                r, c = np.nonzero(m)
                self._rows, self._cols = r.astype(np.int64), c.astype(np.int64)
                self._data = np.asarray(m[r, c], dtype=np.int64)
            else:
                self._rows = self._cols = self._data = np.zeros(0, dtype=np.int64)
        else:
            self._rows = self._cols = self._data = np.zeros(0, dtype=np.int64)

    # --- shape / sparsity -------------------------------------------------
    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @property
    def nnz(self) -> int:
        """Number of present (non-zero) cells."""
        return int(self._data.shape[0])

    @property
    def coo(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The raw ``(family_index, species_index, copy_number)`` arrays."""
        return self._rows, self._cols, self._data

    # --- sparse-native reductions (never densify) -------------------------
    def presence_per_family(self) -> np.ndarray:
        """Number of species each family is present in (length = #families)."""
        return np.bincount(self._rows, minlength=self._shape[0]).astype(np.int64)

    def copies_per_species(self) -> np.ndarray:
        """Total copy number in each species, in ``self.species`` order (length = #species)."""
        return np.bincount(self._cols, weights=self._data,
                           minlength=self._shape[1]).astype(np.int64)

    def copies_per_family(self) -> np.ndarray:
        """Total copy number of each family across all species (length = #families)."""
        return np.bincount(self._rows, weights=self._data,
                           minlength=self._shape[0]).astype(np.int64)

    def copy_values(self) -> np.ndarray:
        """The copy number of every present cell (the non-zero values)."""
        return self._data

    # --- dense view (materialises on demand) ------------------------------
    @property
    def matrix(self) -> np.ndarray:
        """Dense ``(n_families, n_species)`` copy-number array — **densifies on access**.

        Memory is O(#families × #species); at large tip counts prefer the sparse
        reductions above (or :meth:`to_coo_tsv`), which stay O(non-zeros).
        """
        m = np.zeros(self._shape, dtype=np.int64)
        if self._data.shape[0]:
            m[self._rows, self._cols] = self._data
        return m

    def presence(self) -> np.ndarray:
        """Binary presence/absence matrix (copy number > 0). Densifies on access."""
        return (self.matrix > 0).astype(np.int8)

    # --- dense TSV (wide format; for small / full-output runs) -------------
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
        rows, cols, data = [], [], []
        for i, line in enumerate(lines[1:]):
            parts = line.split("\t")
            families.append(parts[0])
            values = parts[1:]
            if len(values) != len(species):
                raise ValueError("row width does not match the number of species columns")
            for j, v in enumerate(values):
                iv = int(v)
                if iv:
                    rows.append(i); cols.append(j); data.append(iv)
        return cls(families=families, species=species, coo=(rows, cols, data))

    # --- sparse TSV (long / COO format; the scalable on-disk output) ------
    def to_coo_tsv(self) -> str:
        """Serialise as a sparse long table: one ``family<TAB>species<TAB>copies`` row per
        present cell. Two header lines, ``#species`` and ``#families``, record the full label
        sets in order, so the round trip is lossless (even all-absent rows/columns survive).
        This is O(non-zeros + labels), the format to use when the dense table would be
        astronomically large."""
        out = ["#species\t" + "\t".join(self.species),
               "#families\t" + "\t".join(self.families),
               "family\tspecies\tcopies"]
        fam, sp = self.families, self.species
        for r, c, v in zip(self._rows.tolist(), self._cols.tolist(), self._data.tolist()):
            out.append(f"{fam[r]}\t{sp[c]}\t{v}")
        return "\n".join(out) + "\n"

    @classmethod
    def from_coo_tsv(cls, source: str) -> "ProfileMatrix":
        """Load a sparse profile written by :meth:`to_coo_tsv`."""
        text = source if "\n" in source else _read_text(source)
        lines = [ln for ln in text.splitlines() if ln.strip()]
        species: list[str] = []
        families: list[str] = []
        triples: list[tuple[str, str, int]] = []
        for ln in lines:
            if ln.startswith("#species\t"):
                species = ln.split("\t")[1:]
            elif ln.startswith("#families\t"):
                families = ln.split("\t")[1:]
            else:
                f, s, v = ln.split("\t")
                if f == "family" and s == "species" and v == "copies":   # the header row itself
                    continue
                triples.append((f, s, int(v)))
        if not families:  # header missing → derive rows from the data
            seen: dict[str, None] = {}
            for f, _, _ in triples:
                seen.setdefault(f, None)
            families = sorted(seen, key=_natkey)
        if not species:  # header missing → derive columns from the data
            seen = {}
            for _, s, _ in triples:
                seen.setdefault(s, None)
            species = sorted(seen, key=_natkey)
        fidx = {f: i for i, f in enumerate(families)}
        sidx = {s: j for j, s in enumerate(species)}
        rows = [fidx[f] for f, _, _ in triples]
        cols = [sidx[s] for _, s, _ in triples]
        data = [v for _, _, v in triples]
        return cls(families=families, species=species, coo=(rows, cols, data))

    # --- construction from a genome dict (sparse; no dense allocation) ----
    @classmethod
    def from_leaf_genomes(cls, leaf_genomes: dict[TreeNode, Genome]) -> "ProfileMatrix":
        species_nodes = sorted(leaf_genomes.keys(), key=lambda n: _natkey(n.name))
        species = [n.name for n in species_nodes]

        famset: set[str] = set()
        for genome in leaf_genomes.values():
            famset.update(genome.families())
        families = sorted(famset, key=_natkey)
        fidx = {f: i for i, f in enumerate(families)}

        rows, cols, data = [], [], []
        for j, node in enumerate(species_nodes):
            genome = leaf_genomes[node]
            for family in genome.families():
                cn = genome.copy_number(family)
                if cn:
                    rows.append(fidx[family]); cols.append(j); data.append(cn)
        return cls(families=families, species=species, coo=(rows, cols, data))
