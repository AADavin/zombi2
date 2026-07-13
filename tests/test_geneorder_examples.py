"""Gene-aware rearrangement worked examples, ported in spirit from Krister Swenson's fork.

The fork (``thekswenson/Zombi``, ``root_genome`` branch, ``tests/test_events.py`` /
``test_genomes.py`` / ``test_divisions.py``) hand-verified rearrangements on tiny gene+intergene
genomes (``30_10.gff`` = 3 genes of length 5; ``30_6.gff`` = 5 genes of length 3), asserting the
*exact* gene order, orientation and segment boundaries after a specific scripted inversion or
transposition (cut-and-paste).

zombi2's :mod:`test_nucleotide_genome` already stress-tests the same geometry the way the fork
did in bulk — random ``(s, ell)`` inversions checked against an independent array oracle. What
those tests do *not* do is exercise the **gene annotation**: a seeded genome whose segments carry
``gene_id`` / gene-block structure. These golden cases fill that gap. The fork used a
two-breakpoint ``(bp1, bp2, direction)`` parameterization on its own intergenic coordinate space;
we re-express each scenario in zombi2's native ``_apply_inversion(start, length)`` (half-open arc
``[s, s+ell)`` on a circular nucleotide genome) and assert gene-level outcomes derived by hand —
so these are correctness checks, not regression pins of current output.

Fork concept -> zombi2 name: *division* -> maximal uncut segment/``Block``; *pieces* ->
ordered segment list; ``ch.genes`` order -> gene-carrying segments in genome order.
"""

import numpy as np

from zombi2.genomes.events import TargetParams
from zombi2.genomes.genome import IdManager
from zombi2.genomes.nucleotide_genome import NucleotideGenome, SegmentRegistry

# Gene layouts as 0-based half-open ancestral intervals, mirroring the fork's GFF fixtures.
G_30_10 = [("1", 0, 5), ("2", 10, 15), ("3", 20, 25)]                       # 30_10.gff
G_30_6 = [("1", 0, 3), ("2", 6, 9), ("3", 12, 15), ("4", 18, 21), ("5", 24, 27)]  # 30_6.gff


def _seed(genes, length):
    """A gene-annotated root nucleotide genome (segments tiled into gene / intergene)."""
    reg = SegmentRegistry(pending_genes=[(a, b, name) for (name, a, b) in genes])
    g = NucleotideGenome(IdManager(), root_length=length, extension=0.9, registry=reg)
    g.originate(np.random.default_rng(0), TargetParams())  # lays down the seed chromosome
    return g


def _gene_order(g):
    """(gene_id, strand) for each gene-carrying segment, in genome order (the fork's ch.genes)."""
    return [(s.gene_id, s.strand) for s in g._segments if s.gene_id is not None]


def _origins(g):
    """Sorted ancestral (source, position) multiset — invariant under any content-preserving event."""
    return sorted((src, p) for (src, p, _st) in g.to_cells())


def _divisions(g):
    """Sorted ancestral source-intervals of every segment — the fork's 'divisions'/'natural_cuts'."""
    return sorted((s.src_start, s.src_end) for s in g._segments)


# --------------------------------------------------------------------------- #
# Inversion enclosing a single gene: that gene flips, order is preserved.
# (fork: test_events.test_inversion_0 / _2 family)
# --------------------------------------------------------------------------- #
def test_inversion_flips_enclosed_gene():
    g = _seed(G_30_10, 30)
    before = _origins(g)
    g._apply_inversion(8, 12)  # arc [8, 20) fully contains gene "2" [10,15), cuts only intergenes
    assert _gene_order(g) == [("1", 1), ("2", -1), ("3", 1)]  # only the enclosed gene reversed
    assert _origins(g) == before                              # content conserved
    assert g.size() == 30


# --------------------------------------------------------------------------- #
# Inversion spanning several genes reverses their order and flips each.
# (fork: test_events.test_inversion_3 -> ['2','1','5','4','3'])
# --------------------------------------------------------------------------- #
def test_inversion_reverses_and_flips_spanned_genes():
    g = _seed(G_30_6, 30)
    before = _origins(g)
    g._apply_inversion(4, 18)  # arc [4, 22) encloses genes 2,3,4 -> become 4,3,2, all reversed
    assert _gene_order(g) == [("1", 1), ("4", -1), ("3", -1), ("2", -1), ("5", 1)]
    assert _origins(g) == before
    assert g.size() == 30


# --------------------------------------------------------------------------- #
# A wrapping inversion (arc crosses the origin) is the fork's trickiest case.
# (fork: test_events.test_inversion_1LEFT / _2 / _Adri_20_09_21_B, "this test wraps")
# --------------------------------------------------------------------------- #
def test_wrapping_inversion_conserves_content_and_flips_ends():
    g = _seed(G_30_6, 30)
    before = _origins(g)
    g._apply_inversion(24, 12)  # [24,30)+[0,6) wraps: gene "5" and gene "1" both flip
    order = _gene_order(g)
    assert order == [("1", -1), ("5", -1), ("2", 1), ("3", 1), ("4", 1)]
    assert {gid for gid, st in order if st == -1} == {"1", "5"}  # only the wrapped ends flipped
    assert _origins(g) == before                                 # nothing lost across the origin
    assert g.size() == 30


# --------------------------------------------------------------------------- #
# Whole-genome inversion reverses every gene and flips them all.
# --------------------------------------------------------------------------- #
def test_whole_genome_inversion_reverses_all_genes():
    g = _seed(G_30_6, 30)
    before = _origins(g)
    g._apply_inversion(0, 30)
    assert _gene_order(g) == [("5", -1), ("4", -1), ("3", -1), ("2", -1), ("1", -1)]
    assert _origins(g) == before
    assert g.size() == 30


# --------------------------------------------------------------------------- #
# Genes are indivisible: an intergenic inversion never splits a gene block.
# (fork invariant behind the whole gene-order model)
#
# ``_apply_inversion`` is the raw primitive and takes *physical* coordinates, so a fixed arc
# list would drift into a gene interior once earlier inversions permute the layout. We therefore
# apply one inversion per fresh genome, with arcs whose endpoints are legal on the known tiling
# (genes at [0,5) [10,15) [20,25); intergene interiors 6-9, 16-19, 26-29; boundaries 0/5/10/…).
# --------------------------------------------------------------------------- #
def test_genes_stay_indivisible_under_intergenic_inversion():
    expected = {name: (a, b) for name, a, b in G_30_10}  # gene id -> its one ancestral interval
    for s, ell in [(8, 12), (5, 20), (16, 9), (26, 9)]:  # each arc encloses whole genes only
        g = _seed(G_30_10, 30)
        before = _origins(g)
        g._apply_inversion(s, ell)
        by_id = {}
        for seg in g._segments:
            if seg.gene_id is not None:
                by_id.setdefault(seg.gene_id, []).append((seg.src_start, seg.src_end))
        # every gene survives as exactly one intact block spanning its full ancestral interval
        for gid, ivals in by_id.items():
            assert ivals == [expected[gid]], f"gene {gid} was split or moved: {ivals}"
        assert set(by_id) == set(expected)  # all three genes still present
        assert _origins(g) == before        # nothing lost


# --------------------------------------------------------------------------- #
# Divisions subdivide only at the cut, and tile the genome. (fork: test_divisions,
# init_divisions / natural_cuts).
# --------------------------------------------------------------------------- #
def test_inversion_subdivides_only_the_cut_intergene():
    g = _seed(G_30_10, 30)
    assert _divisions(g) == [(0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30)]
    g._apply_inversion(8, 12)  # cuts inside intergene [5,10) at 8; 20 is already a boundary
    # only [5,10) splits into [5,8) + [8,10); every other division is untouched
    assert _divisions(g) == [(0, 5), (5, 8), (8, 10), (10, 15), (15, 20), (20, 25), (25, 30)]
    # divisions still tile [0, 30) with no gap or overlap
    covered = sorted(_divisions(g))
    assert covered[0][0] == 0 and covered[-1][1] == 30
    assert all(a[1] == b[0] for a, b in zip(covered, covered[1:]))


# --------------------------------------------------------------------------- #
# Transposition (cut a block, paste it elsewhere) moves genes WITHOUT flipping
# them — the defining difference from an inversion. (fork: test_genomes.cut_and_paste)
# --------------------------------------------------------------------------- #
def test_transposition_moves_gene_block_forward_keeping_orientation():
    g = _seed(G_30_6, 30)
    before = _origins(g)
    g._apply_transposition(3, 6, 21)  # cut arc [3,9) (intergene + gene "2"), paste at physical 21
    order = _gene_order(g)
    assert order == [("1", 1), ("3", 1), ("4", 1), ("5", 1), ("2", 1)]  # "2" moved to the end
    assert all(st == 1 for _gid, st in order)                            # NO strand flip
    assert _origins(g) == before                                        # content conserved
    assert g.size() == 30


def test_transposition_moves_gene_block_backward():
    g = _seed(G_30_6, 30)
    before = _origins(g)
    g._apply_transposition(15, 6, 3)  # cut arc [15,21) (intergene + gene "4"), paste after gene "1"
    assert _gene_order(g) == [("1", 1), ("4", 1), ("2", 1), ("3", 1), ("5", 1)]
    assert _origins(g) == before
    assert g.size() == 30


def test_transposition_of_wrapping_arc_conserves_content():
    g = _seed(G_30_6, 30)
    before = _origins(g)
    g._apply_transposition(27, 6, 15)  # cut wrapping arc [27,30)+[0,3) (intergene + gene "1")
    order = _gene_order(g)
    assert order == [("2", 1), ("3", 1), ("1", 1), ("4", 1), ("5", 1)]  # "1" reinserted mid-genome
    assert all(st == 1 for _gid, st in order)
    assert _origins(g) == before
    assert g.size() == 30
