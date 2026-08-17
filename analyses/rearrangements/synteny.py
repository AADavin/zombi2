"""The gene-order summary statistics this study fits, and the divergence axis they are binned on.

A genome is represented as its signed gene order, one row per gene, in genomic order:

    {species: [(chromosome, gene, strand in {+1, -1}), ...]}

Three statistics are computed for every pair of genomes. Each is named here exactly as it is
named in ``REPORT.md``.

**Gene-order conservation.** Consecutive core genes on one chromosome form an *adjacency*, the
unordered pair ``{f, g}``. The *breakpoint distance* between two genomes is the number of
adjacencies present in one and absent in the other. Gene-order conservation is its normalised
complement: the fraction of adjacencies the two genomes share. Conservation is used rather than
the raw count only so that genomes with different numbers of core adjacencies stay comparable;
the two carry the same information.

**Conserved segment length.** The maximal runs of genes whose consecutive adjacencies all survive
in the other genome. This is the spacing between breakpoints, that is, the same breaks seen from
the other side. It is reported as the mean run length in genes.

**Cross-chromosome breaks.** Of the adjacencies that are broken, the fraction whose two genes lie
on different chromosomes in the other genome. An inversion acts inside one chromosome, so it
breaks adjacencies but leaves both genes on the same chromosome. A translocation moves a segment
between chromosomes, so it separates them. This statistic is therefore what tells the two event
types apart, and it is what makes the translocation rate readable in the mixed arm.

The divergence axis is the patristic distance between two tips, the time from one tip up to their
common ancestor and back down to the other. Rearrangements accumulate along both of those paths,
so patristic time is the quantity the decay of conservation is read against.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

# Divergence-time bins, in tree time units. The tree is scaled to a crown depth of
# CROWN_DEPTH (see experiment.py), so patristic distances run from 0 to 2 * CROWN_DEPTH.
DEFAULT_BINS = (0.0, 60.0, 120.0, 200.0)


# ------------------------------------------------------------------ the divergence axis
def pairwise_divergence(tree, namemap) -> dict[frozenset, float]:
    """``{frozenset({tip_a, tip_b}): patristic time}`` for every pair of extant tips.

    Computed from the ZOMBI2 tree directly. Every node carries ``birth_time`` and ``end_time``,
    so the common ancestor's ``end_time`` is the time the two lineages split.
    """
    tips = list(tree.extant_leaves())
    # ancestor chain of every tip, root first, so the last shared entry is the common ancestor
    chain: dict[int, list[int]] = {}
    for tip in tips:
        path, node = [], tip
        while node is not None:
            path.append(node)
            node = tree.nodes[node].parent
        chain[tip] = path[::-1]
    out: dict[frozenset, float] = {}
    for a, b in combinations(tips, 2):
        ca, cb = chain[a], chain[b]
        shared = 0
        for x, y in zip(ca, cb):
            if x != y:
                break
            shared += 1
        mrca = ca[shared - 1]
        split = tree.nodes[mrca].end_time
        depth_a = tree.nodes[a].end_time
        depth_b = tree.nodes[b].end_time
        out[frozenset((namemap[a], namemap[b]))] = (depth_a - split) + (depth_b - split)
    return out


# ------------------------------------------------------------------ the core gene set
def core_families(genomes: dict[str, list[tuple]]) -> set:
    """Genes present in exactly one copy in every genome.

    Under the nucleotide model a breakpoint never falls strictly inside a gene, and this study
    simulates no gain or loss, so in practice every gene qualifies and the core is the whole gene
    set. The filter is kept so that the statistics are computed over a one-to-one gene set by
    construction rather than by assumption.
    """
    counts: dict = {}
    for rows in genomes.values():
        seen: dict = {}
        for _chrom, gene, _strand in rows:
            seen[gene] = seen.get(gene, 0) + 1
        for gene, n in seen.items():
            counts.setdefault(gene, []).append(n)
    n_species = len(genomes)
    return {g for g, ns in counts.items() if len(ns) == n_species and max(ns) == 1}


# ------------------------------------------------------------------ per-genome precomputation
def adjacencies(rows, core):
    """``(adjacency set, signed adjacency set, per-chromosome gene order, gene -> chromosome)``.

    A signed adjacency records *which ends* of the two genes meet, so a pair that is still
    adjacent but in flipped relative orientation has the same unsigned adjacency and a different
    signed one.
    """
    plain: set = set()
    signed: set = set()
    chroms: dict = {}
    where: dict = {}
    previous = None
    for chrom, gene, strand in rows:
        if gene not in core:
            continue
        chroms.setdefault(chrom, []).append(gene)
        where[gene] = chrom
        if previous is not None and previous[0] == chrom and previous[1] != gene:
            _pc, pgene, pstrand = previous
            plain.add(frozenset((pgene, gene)))
            # a '+' gene is read tail then head, so its right end is its head
            right = (pgene, "h") if pstrand > 0 else (pgene, "t")
            left = (gene, "t") if strand > 0 else (gene, "h")
            signed.add(frozenset((right, left)))
        previous = (chrom, gene, strand)
    return plain, signed, list(chroms.values()), where


def precompute(genomes, core):
    return {sp: dict(zip(("plain", "signed", "chroms", "where"), adjacencies(rows, core)))
            for sp, rows in genomes.items()}


# ------------------------------------------------------------------ the three statistics
def gene_order_conservation(a, b) -> float:
    """Fraction of adjacencies the two genomes share, that is, the normalised complement of the
    breakpoint distance."""
    pa, pb = a["plain"], b["plain"]
    denominator = 0.5 * (len(pa) + len(pb))
    if denominator == 0:
        return np.nan
    return len(pa & pb) / denominator


def conserved_segment_lengths(a, b) -> list[int]:
    """Maximal runs of ``a``'s core genes whose consecutive adjacencies all survive in ``b``."""
    shared = b["plain"]
    sizes = []
    for chrom in a["chroms"]:
        if not chrom:
            continue
        run = 1
        for i in range(1, len(chrom)):
            if frozenset((chrom[i - 1], chrom[i])) in shared:
                run += 1
            else:
                sizes.append(run)
                run = 1
        sizes.append(run)
    return sizes


def cross_chromosome_breaks(a, b) -> float:
    """Of ``a``'s adjacencies that are broken in ``b``, the fraction whose two genes lie on
    different chromosomes in ``b``. Reads translocations, which inversions cannot produce."""
    shared, where = b["plain"], b["where"]
    broken = cross = 0
    for chrom in a["chroms"]:
        for i in range(1, len(chrom)):
            f, g = chrom[i - 1], chrom[i]
            if frozenset((f, g)) in shared:
                continue
            if f not in where or g not in where:
                continue
            broken += 1
            if where[f] != where[g]:
                cross += 1
    return cross / broken if broken else np.nan


# ------------------------------------------------------------------ the summary vector
KEYS = ("conservation", "segment", "cross")


def summary(genomes, divergence, bins=DEFAULT_BINS, core=None) -> dict:
    """Mean of each statistic in each divergence-time bin, plus the core size.

    Keys are ``conservation{i}``, ``segment{i}`` and ``cross{i}`` for bin ``i``.
    """
    if core is None:
        core = core_families(genomes)
    pc = precompute(genomes, core)
    n_bins = len(bins) - 1
    buckets = {k: [[] for _ in range(n_bins)] for k in KEYS}
    for a, b in combinations(list(genomes), 2):
        key = frozenset((a, b))
        if key not in divergence:
            continue
        time = divergence[key]
        index = None
        for i in range(n_bins):
            if bins[i] <= time < bins[i + 1]:
                index = i
                break
        if index is None:
            continue
        pa, pb = pc[a], pc[b]
        runs = conserved_segment_lengths(pa, pb)
        values = {"conservation": gene_order_conservation(pa, pb),
                  "segment": float(np.mean(runs)) if runs else np.nan,
                  "cross": cross_chromosome_breaks(pa, pb)}
        for k, v in values.items():
            if not np.isnan(v):
                buckets[k][index].append(v)
    out = {"n_core": float(len(core))}
    for k in KEYS:
        for i in range(n_bins):
            vals = buckets[k][i]
            out[f"{k}{i}"] = float(np.mean(vals)) if vals else np.nan
    return out


def fit_keys(bins=DEFAULT_BINS, use_cross=False) -> list[str]:
    """The keys the ABC distance is computed over.

    Gene-order conservation and conserved segment length always. Cross-chromosome breaks only
    in the mixed arm, where a translocation rate has to be read.
    """
    n = len(bins) - 1
    keys = [f"conservation{i}" for i in range(n)] + [f"segment{i}" for i in range(n)]
    if use_cross:
        keys += [f"cross{i}" for i in range(n)]
    return keys


def all_keys(bins=DEFAULT_BINS) -> list[str]:
    n = len(bins) - 1
    return [f"{k}{i}" for k in KEYS for i in range(n)]


def distance(candidate, target, keys, weights=None) -> float:
    """Scaled Euclidean distance between two summary vectors over ``keys``."""
    total = n = 0.0
    for k in keys:
        va, vb = candidate.get(k, np.nan), target.get(k, np.nan)
        if np.isnan(va) or np.isnan(vb):
            continue
        w = (weights or {}).get(k, 1.0)
        total += w * (va - vb) ** 2
        n += 1
    return float(np.sqrt(total / n)) if n else np.inf


def scales(summaries, keys) -> dict:
    """Per-key spread across the grid, so statistics on different scales weigh comparably.

    Conservation is a fraction and conserved segment length is a count of genes, so an unweighted
    distance would be dominated by the segment term.
    """
    out = {}
    for k in keys:
        col = np.array([s.get(k, np.nan) for s in summaries], float)
        col = col[~np.isnan(col)]
        sd = float(np.std(col)) if col.size else 0.0
        out[k] = sd if sd > 0 else 1.0
    return out
