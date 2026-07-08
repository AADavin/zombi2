"""GFF-driven genome selection (**P2.5**): evolve a *real annotated genome* down a tree, applying
codon-level selection (P2) to protein-coding genes and neutral nucleotide evolution to everything else.

You supply a root genome (one nucleotide contig) and its CDS annotations (from a GFF, with **strand**
and reading frame). :class:`GenomeSelection` partitions the genome into coding intervals and the
intergenic gaps between them; each CDS is translated and evolved under
:class:`~zombi2.experimental.codon_selection.CodonSelection` (mutation on DNA, selection on the encoded
protein -> emergent dN/dS), while intergenic DNA drifts neutrally under the same nucleotide model. The
segments are then reassembled, in coordinate order, into a whole evolved genome at every tree node.

This is the "specify a real genome at the root" scenario: vertical inheritance down one tree (no DTL
gene-family dynamics here -- that is a later slice). **Preconditions (validated -- unsupported inputs
fail loudly, they are never silently mishandled):**

* the genome is upper-cased on input (soft-masking / lowercase is not preserved);
* CDS are **single-exon**, in frame (**phase 0**), length a multiple of 3, over ``ACGT`` sense codons
  (a terminal stop codon is kept fixed and excluded from selection); phased/partial and multi-exon
  CDS are **rejected**;
* non-``ACGT`` bases (``N`` / assembly gaps / IUPAC codes) are allowed only in **non-coding** DNA,
  where they are preserved unchanged (they carry no signal); inside a CDS they are an error.

Like the rest of the feature, nothing here imports torch/esm at module load.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from zombi2.experimental import warn_experimental
from zombi2.experimental.codon_selection import CodonSelection, translate
from zombi2.experimental.selection import Critic
from zombi2.sequences.models import SubstitutionModel, evolve_on_tree, reverse_complement

__all__ = ["CDS", "GenomeSelection", "read_cds_gff"]

_ACGT = frozenset("ACGT")


@dataclass(frozen=True)
class CDS:
    """A single-exon coding interval: 0-based half-open ``[start, end)`` on the forward strand, with
    ``strand`` +1 (sense) or -1 (antisense), ``phase`` (GFF reading-frame offset; only 0 supported).
    ``end - start`` must be a multiple of 3."""

    start: int
    end: int
    strand: int
    name: str = ""
    phase: int = 0


def _attrs(field: str) -> dict:
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_cds_gff(path, *, seqid: str | None = None) -> list[CDS]:
    """Read single-exon **CDS** features (keeping strand + phase) from a GFF3 -> a sorted ``list[CDS]``.

    Coordinates convert from GFF 1-based inclusive to 0-based half-open. If ``seqid`` is given only that
    sequence is kept; if it is ``None`` and the GFF spans several sequences, that is an error (pass
    ``seqid=``). Multi-exon CDS (several CDS lines sharing an ``ID``/``Parent``) and malformed strands
    are rejected. (ZOMBI2's own :func:`~zombi2.genomes.gff.read_gff` drops strand/phase.)
    """
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    rows, group_counts = [], {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] != "CDS":
                continue
            if seqid is not None and f[0] != seqid:
                continue
            strand = {"+": 1, "-": -1}.get(f[6])
            if strand is None:
                raise ValueError(f"CDS at {f[0]}:{f[3]}-{f[4]} has unsupported strand {f[6]!r}")
            phase = int(f[7]) if f[7] in ("0", "1", "2") else 0
            attrs = _attrs(f[8]) if len(f) > 8 else {}
            gid = attrs.get("Parent") or attrs.get("ID")
            name = attrs.get("locus_tag") or attrs.get("Name") or gid or f"cds{len(rows)}"
            rows.append((f[0], CDS(int(f[3]) - 1, int(f[4]), strand, name, phase)))
            if gid:
                group_counts[gid] = group_counts.get(gid, 0) + 1
    multi = sorted(g for g, n in group_counts.items() if n > 1)
    if multi:
        raise ValueError(f"multi-exon CDS not supported: ID/Parent {multi} span multiple CDS lines")
    contigs = {c for c, _ in rows}
    if seqid is None and len(contigs) > 1:
        raise ValueError(f"GFF spans multiple sequences {sorted(contigs)}; pass seqid= to pick one")
    return sorted((cds for _, cds in rows), key=lambda c: c.start)


class GenomeSelection:
    """Overlay language-model-guided selection on a whole annotated genome (CDS = codon selection,
    non-coding = neutral). ``critic``/``beta``/``nuc_model`` are handed to an inner
    :class:`~zombi2.experimental.codon_selection.CodonSelection`; the same ``nuc_model`` drives the
    neutral intergenic evolution, so coding and non-coding share one mutation process."""

    def __init__(self, critic: Critic, *, beta: float = 1.0,
                 nuc_model: SubstitutionModel | None = None):
        warn_experimental("GenomeSelection")
        self.codon = CodonSelection(critic, beta=beta, nuc_model=nuc_model)
        self.nuc = self.codon.nuc

    def _segments(self, length: int, cds: list[CDS]):
        """Ordered, non-overlapping segments covering ``[0, length)``: coding CDS + intergenic gaps."""
        segs, pos = [], 0
        for c in sorted(cds, key=lambda c: c.start):
            if c.strand not in (1, -1):
                raise ValueError(f"CDS {c.name!r} strand must be +1 or -1, got {c.strand}")
            if c.phase != 0:
                raise ValueError(f"CDS {c.name!r} has phase {c.phase} != 0 "
                                 "(phased/partial CDS are not supported)")
            if not 0 <= c.start < c.end <= length:
                raise ValueError(f"CDS {c.name!r} interval [{c.start},{c.end}) out of [0,{length})")
            if c.start < pos:
                raise ValueError(f"CDS {c.name!r} overlaps a previous feature")
            if (c.end - c.start) % 3:
                raise ValueError(f"CDS {c.name!r} length {c.end - c.start} is not a multiple of 3")
            if c.start > pos:
                segs.append(("nc", pos, c.start, 1, None))
            segs.append(("cds", c.start, c.end, c.strand, c.name))
            pos = c.end
        if pos < length:
            segs.append(("nc", pos, length, 1, None))
        return segs

    def _evolve_cds(self, root, subst, coding: str, name: str, rng) -> dict:
        """Evolve one forward-strand ``ACGT`` coding sequence; a terminal stop codon is kept fixed."""
        bad = next((i for i, ch in enumerate(coding) if ch not in _ACGT), None)
        if bad is not None:
            raise ValueError(f"CDS {name!r} coding position {bad} is not ACGT ({coding[bad]!r})")
        stop = ""
        if len(coding) >= 3 and translate(coding[-3:]) == "*":
            stop, coding = coding[-3:], coding[:-3]
        if not coding:
            raise ValueError(f"CDS {name!r} has no sense codons to evolve")
        ev = self.codon.evolve_coding_family(root, subst, coding, rng=rng)
        return {gid: seq + stop for gid, seq in ev.items()}

    def _evolve_nc(self, root, subst, sub: str, rng) -> dict:
        """Evolve one intergenic segment neutrally, keeping any non-ACGT base (N/gap/IUPAC) unchanged."""
        ev = evolve_on_tree(root, subst, self.nuc, rng, root_seq=sub)
        fixed = {i for i, ch in enumerate(sub) if ch not in _ACGT}
        if fixed:                                     # ambiguity codes carry no signal -> freeze them
            ev = {g: "".join(sub[i] if i in fixed else s[i] for i in range(len(s)))
                  for g, s in ev.items()}
        return ev

    def evolve_genome(self, root, subst: dict, genome: str, cds: list[CDS], *,
                      rng: np.random.Generator) -> dict:
        """Evolve the ``genome`` (one nucleotide contig, 5'->3') down the node tree ``(root, subst)``:
        each CDS under codon selection (respecting strand), the rest neutrally. Returns
        ``{node.gid: genome}`` -- the reassembled genome at every node (length preserved; a zero-length
        root branch reproduces the upper-cased input, non-ACGT bases and all)."""
        if not genome:
            raise ValueError("genome is empty")
        genome = genome.upper()
        per_seg = []                                          # {gid: genomic_seq} per segment, in order
        for kind, a, b, strand, name in self._segments(len(genome), cds):
            sub = genome[a:b]
            if kind == "nc":
                per_seg.append(self._evolve_nc(root, subst, sub, rng))
            else:
                coding = sub if strand == 1 else reverse_complement(sub)
                ev = self._evolve_cds(root, subst, coding, name, rng)
                per_seg.append({g: (s if strand == 1 else reverse_complement(s))
                                for g, s in ev.items()})
        return {gid: "".join(seg[gid] for seg in per_seg) for gid in per_seg[0]}
