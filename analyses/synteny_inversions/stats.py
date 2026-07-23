"""Signed, within-chromosome micro-synteny statistics — the observables that read the inversion rate.

Micro-synteny (local gene order + orientation) persists only among close relatives, and it is where
inversions leave a readable mark. Representation (real or simulated, multi-chromosome):

    genomes = {species: [(chromosome, family, orientation ∈ {+1,-1}), …]}   in genomic order

Between consecutive core genes on one chromosome there are two adjacency notions:

  * UNSIGNED  {f, g}                      — retained micro-synteny (any rearrangement breaks it)
  * SIGNED    {right_end(f), left_end(g)} — retained *and* same relative orientation

A pair adjacent in both genomes but joined by different ends is a **flipped adjacency** — an
inversion-boundary signature. Gene-order conservation (unsigned) and conserved block size are the
ABC fit targets; the flip fraction and reversed-block size are carried as checks.

Ported from the retired ``synteny/scripts/signed_chrom_stats.py`` + ``synteny_stats.py``; the
pairwise-divergence axis is computed from the dated tree with ete3.
"""
from __future__ import annotations

import json
from itertools import combinations

import numpy as np
from ete3 import Tree as ETree

DEFAULT_BINS = (0.0, 70.0, 130.0, 200.0)


def _extremities(fam, sign: int):
    """(left_end, right_end) of a gene read left→right: a '+' gene points tail→head."""
    h, t = (fam, "h"), (fam, "t")
    return (t, h) if sign > 0 else (h, t)


def load_signed(path: str) -> dict[str, list[tuple]]:
    raw = json.load(open(path))
    return {sp: [(str(c), str(f), int(o)) for c, f, o in rows] for sp, rows in raw.items()}


def pairwise_divergence(nwk_path: str) -> dict[frozenset, float]:
    """{frozenset({A, B}) -> path length in Myr} for every tip pair — the divergence axis.
    Rearrangements accumulate along both branches from the MRCA, so the patristic distance
    t_A + t_B − 2 t_MRCA is the natural axis."""
    t = ETree(open(nwk_path).read(), format=1)
    leaves = t.get_leaves()
    return {frozenset((a.name, b.name)): a.get_distance(b) for a, b in combinations(leaves, 2)}


def core_single_copy(genomes, min_frac=1.0) -> set:
    """Families present in exactly one copy in at least ``min_frac`` of species."""
    n = len(genomes)
    fam_species: dict = {}
    fam_multi: set = set()
    for sp, rows in genomes.items():
        seen: dict = {}
        for _c, f, _o in rows:
            seen[f] = seen.get(f, 0) + 1
        for f, c in seen.items():
            fam_species.setdefault(f, set()).add(sp)
            if c > 1:
                fam_multi.add(f)
    thr = min_frac * n
    return {f for f, sps in fam_species.items() if len(sps) >= thr and f not in fam_multi}


def within_chrom_adj(rows, core):
    """(unsigned adjacencies, signed adjacencies, per-chromosome core-gene order)."""
    uns, sgn = set(), set()
    chroms: dict = {}
    prev = None
    for c, f, o in rows:
        if f not in core:
            continue
        chroms.setdefault(c, []).append((f, o))
        if prev is not None and prev[0] == c and prev[1] != f:
            _pc, pf, po = prev
            uns.add(frozenset((pf, f)))
            right = _extremities(pf, po)[1]
            left = _extremities(f, o)[0]
            sgn.add(frozenset((right, left)))
        prev = (c, f, o)
    return uns, sgn, list(chroms.values())


def precompute(genomes, core):
    pc = {}
    for sp, rows in genomes.items():
        u, s, ch = within_chrom_adj(rows, core)
        pos = {}
        for ci, chrom in enumerate(ch):
            for j, (f, _o) in enumerate(chrom):
                pos[f] = (ci, j)
        pc[sp] = {"u": u, "s": s, "chroms": ch, "pos": pos}
    return pc


def pair_stats(a, b):
    """(unsigned conservation fraction, orientation-flip fraction) for a species pair."""
    ua, ub, sa, sb = a["u"], b["u"], a["s"], b["s"]
    denom = 0.5 * (len(ua) + len(ub))
    if denom == 0:
        return np.nan, np.nan
    u_int = len(ua & ub)
    unsigned_frac = u_int / denom
    flipped = (u_int - len(sa & sb)) / u_int if u_int > 0 else np.nan
    return unsigned_frac, flipped


def conserved_block_sizes(a, b):
    """Maximal runs of A's core genes whose consecutive adjacency also survives in B (in genes)."""
    bshared = b["u"]
    sizes = []
    for chrom in a["chroms"]:
        if not chrom:
            continue
        run = 1
        for i in range(1, len(chrom)):
            if frozenset((chrom[i - 1][0], chrom[i][0])) in bshared:
                run += 1
            else:
                sizes.append(run)
                run = 1
        sizes.append(run)
    return sizes


def pair_block_mean(a, b):
    bs = conserved_block_sizes(a, b)
    return float(np.mean(bs)) if bs else np.nan


def reversed_block_sizes(a, b):
    """Runs of A's genes that are consecutive in B but in reversed order — the footprint of a single
    inversion, so the run length reads the event size directly (the observable that sees *size*)."""
    posB = b["pos"]
    sizes = []
    for chrom in a["chroms"]:
        slopes = []
        for i in range(1, len(chrom)):
            pf = chrom[i - 1][0]; f = chrom[i][0]
            if pf in posB and f in posB:
                (cb1, jb1), (cb2, jb2) = posB[pf], posB[f]
                slopes.append(jb2 - jb1 if cb1 == cb2 and abs(jb2 - jb1) == 1 else 0)
            else:
                slopes.append(0)
        minority = -1 if slopes.count(1) >= slopes.count(-1) else 1
        run = 0
        for s in slopes:
            if s == minority:
                run += 1
            else:
                if run > 0:
                    sizes.append(run + 1)
                run = 0
        if run > 0:
            sizes.append(run + 1)
    return sizes


def pair_rev_mean(a, b):
    rs = reversed_block_sizes(a, b)
    return float(np.mean(rs)) if rs else np.nan


def summary_vector(genomes, divergence, bins=DEFAULT_BINS, core=None, min_frac=1.0):
    """Per divergence-time bin: mean gene-order conservation (uns), conserved block size (blk),
    orientation-flip fraction (flip), reversed-block size (rev). ``divergence`` is the
    {frozenset(pair): Myr} map from :func:`pairwise_divergence`."""
    if core is None:
        core = core_single_copy(genomes, min_frac)
    pc = precompute(genomes, core)
    names = list(genomes)
    ub = [[] for _ in range(len(bins) - 1)]
    bb = [[] for _ in range(len(bins) - 1)]
    fb = [[] for _ in range(len(bins) - 1)]
    rb = [[] for _ in range(len(bins) - 1)]
    for a, b in combinations(names, 2):
        key = frozenset((a, b))
        if key not in divergence:
            continue
        uf, ff = pair_stats(pc[a], pc[b])
        blk = pair_block_mean(pc[a], pc[b])
        rev = pair_rev_mean(pc[a], pc[b])
        d = divergence[key]
        for i in range(len(bins) - 1):
            if bins[i] <= d < bins[i + 1]:
                if not np.isnan(uf):
                    ub[i].append(uf)
                if not np.isnan(blk):
                    bb[i].append(blk)
                if not np.isnan(ff):
                    fb[i].append(ff)
                if not np.isnan(rev):
                    rb[i].append(rev)
                break
    out = {"n_core": float(len(core))}
    for i in range(len(bins) - 1):
        out[f"uns{i}"] = float(np.mean(ub[i])) if ub[i] else np.nan
        out[f"blk{i}"] = float(np.mean(bb[i])) if bb[i] else np.nan
        out[f"flip{i}"] = float(np.mean(fb[i])) if fb[i] else np.nan
        out[f"rev{i}"] = float(np.mean(rb[i])) if rb[i] else np.nan
    return out


def summary_keys(bins=DEFAULT_BINS):
    """Keys the ABC distance is computed over: conservation + conserved block size."""
    return ([f"uns{i}" for i in range(len(bins) - 1)]
            + [f"blk{i}" for i in range(len(bins) - 1)])


def flip_keys(bins=DEFAULT_BINS):
    return [f"flip{i}" for i in range(len(bins) - 1)]


def rev_keys(bins=DEFAULT_BINS):
    return [f"rev{i}" for i in range(len(bins) - 1)]


def distance(sa, sb, keys=None, weights=None):
    keys = keys or [k for k in sa if k.startswith("uns")]
    tot, n = 0.0, 0
    for k in keys:
        va, vb = sa.get(k, np.nan), sb.get(k, np.nan)
        if np.isnan(va) or np.isnan(vb):
            continue
        w = (weights or {}).get(k, 1.0)
        tot += w * (va - vb) ** 2
        n += 1
    return float(np.sqrt(tot / n)) if n else np.inf
