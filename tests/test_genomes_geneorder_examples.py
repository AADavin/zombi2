"""Gene-aware rearrangement worked examples, ported from Krister Swenson's fork.

The fork (``thekswenson/Zombi``, ``root_genome`` branch: ``tests/test_events.py``,
``tests/test_genomes.py``, ``tests/test_divisions.py``) hand-verified rearrangements on two tiny
genomes — ``30_6.gff`` (30 bp, 5 genes of 3 bp) and ``30_10.gff`` (30 bp, 3 genes of 5 bp) — and
asserted the *exact* gene order, orientation and segment boundaries after a scripted inversion.

zombi2 already stress-tests the same geometry in bulk: ``test_genomes_nucleotide`` fires random
inversions and checks each against an independent brute-force array oracle. What that leaves
untested is the **gene annotation** — a genome whose blocks are declared, indivisible genes
separated by spacer. These worked examples fill exactly that gap, and being hand-derived they pin
down *what the model should do*, which random testing never states.

The two genomes are rebuilt here as GFF fixtures, so they are the fork's, not a paraphrase. The
fork parameterised an inversion by two intergenic breakpoints and a direction
(``make_inversion_intergenic(ch, bp1, bp2, T_DIR)``); zombi2 takes a half-open arc,
``Chromosome.invert(start, length)``, on a circular nucleotide genome. Every expected value below is
re-derived in zombi2's coordinates rather than transcribed.

Concept map — fork → zombi2: *division* → :class:`~zombi2.genomes.nucleotide.Block` (a run of one
unbroken ancestry); *pieces* → the ordered block list; ``ch.genes`` → the gene-carrying blocks in
physical order; ``natural_cuts`` → block boundaries; ``Inversion.afterToBeforeT`` → ``trace_back()``.

Not ported: the fork's per-nucleotide ``afterToBeforeS``/``afterToBeforeT`` arithmetic cases, whose
content the oracle already proves for arbitrary inputs (one scenario is kept below, to fix the
*convention*); and its transfer worked examples, which need a donor and a recipient alive at the
same instant — a global timeline, not a scripted chromosome.

Every example runs the same way the fork's did: a run with all rates zero, so the root genome is
exactly the declaration, and then the events applied by hand. What the engine draws from the rng,
these call outright — the mutators take explicit coordinates (``invert``, ``duplicate``, ``delete``,
``originate``, ``excise``/``place``), and the engine calls those same ones after picking, so there
is no second code path to drift.
"""

import pytest

from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.genomes.nucleotide import _CutsGene
from zombi2.species import read_newick

#: the fork's fixtures: (total bp, gene length, gene count). 30_6 = 5 genes of 3 bp with 3 bp
#: spacers; 30_10 = 3 genes of 5 bp with 5 bp spacers.
FORK_30_6 = (30, 3, 5)
FORK_30_10 = (30, 5, 3)


def _gff(tmp_path, fixture, name="seed.gff"):
    """One of the fork's genomes as a GFF: genes evenly spaced, each followed by a spacer of its
    own length."""
    total, gene_len, n = fixture
    stride = total // n
    lines = ["##gff-version 3", f"##sequence-region c 1 {total}"]
    lines += [f"c\tfork\tgene\t{i * stride + 1}\t{i * stride + gene_len}\t.\t+\t.\tID=g{i + 1}"
              for i in range(n)]
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def _seed(tmp_path, fixture=FORK_30_6):
    """The seed chromosome, with every rate zero so the root genome is exactly the declaration."""
    tree, _ = read_newick("(A:1.0,B:1.0);")
    r = simulate_genomes_nucleotide(tree, gff=_gff(tmp_path, fixture), seed=1)
    return r.genomes[tree.root].chromosomes[0]


def _layout(chrom):
    """``[(position, length, gene, strand)]`` for every block, in physical order — the fork's view
    of a chromosome. ``gene`` is 0 for spacer."""
    out, at = [], 0
    for b in chrom.blocks:
        out.append((at, b.end - b.start, b.gene, b.strand))
        at += b.end - b.start
    return out


def _genes(chrom):
    """``[(gene, strand)]`` in physical order — the fork's ``ch.genes``, its signed gene order."""
    return [(b.gene, b.strand) for b in chrom.blocks if b.gene]


def _gene_cycle(chrom):
    """The signed gene order as a **cycle**, normalised to start at the lowest gene id.

    A circular chromosome has no privileged origin, and an inversion spanning it re-roots the
    physical coordinates (see ``test_a_wrapping_inversion_re_roots_the_circular_coordinates``), so
    order-up-to-rotation is the right comparison whenever the origin is not the point.
    """
    genes = _genes(chrom)
    if not genes:
        return []
    i = min(range(len(genes)), key=lambda k: genes[k][0])
    return genes[i:] + genes[:i]


# --- the fixtures are the fork's genomes ----------------------------------------------------------

def test_the_seed_reproduces_the_forks_30_6_genome(tmp_path):
    # 5 genes of 3 bp, each followed by a 3 bp spacer — the fork's 30_6.gff
    assert _layout(_seed(tmp_path, FORK_30_6)) == [
        (0, 3, 1, 1), (3, 3, 0, 1), (6, 3, 2, 1), (9, 3, 0, 1), (12, 3, 3, 1),
        (15, 3, 0, 1), (18, 3, 4, 1), (21, 3, 0, 1), (24, 3, 5, 1), (27, 3, 0, 1)]


def test_the_seed_reproduces_the_forks_30_10_genome(tmp_path):
    # 3 genes of 5 bp, each followed by a 5 bp spacer — the fork's 30_10.gff
    assert _layout(_seed(tmp_path, FORK_30_10)) == [
        (0, 5, 1, 1), (5, 5, 0, 1), (10, 5, 2, 1), (15, 5, 0, 1), (20, 5, 3, 1), (25, 5, 0, 1)]


# --- inversions: the fork's scenarios, re-derived --------------------------------------------------

def test_an_inversion_inside_one_intergene_moves_no_gene(tmp_path):
    # the arc [3, 6) is exactly the first spacer. No gene moves or turns — but the spacer itself is
    # now read backwards, which is the honest answer: its nucleotides really were reversed.
    chrom = _seed(tmp_path)
    before = _layout(chrom)
    chrom.invert(3, 3)
    assert _genes(chrom) == _genes(_seed(tmp_path))
    assert _layout(chrom) == [(p, n, g, -s if (p, n, g) == (3, 3, 0) else s)
                              for (p, n, g, s) in before]


def test_an_inversion_over_one_gene_flips_it_and_mirrors_its_place(tmp_path):
    # the arc [3, 9) holds spacer(3,6) then g2(6,9). Reversed, the arc reads g2 then spacer, so g2
    # lands at 3 pointing backwards and the spacer follows it — the rest of the genome is untouched.
    chrom = _seed(tmp_path)
    chrom.invert(3, 6)
    assert _layout(chrom) == [
        (0, 3, 1, 1), (3, 3, 2, -1), (6, 3, 0, -1), (9, 3, 0, 1), (12, 3, 3, 1),
        (15, 3, 0, 1), (18, 3, 4, 1), (21, 3, 0, 1), (24, 3, 5, 1), (27, 3, 0, 1)]


def test_a_multi_gene_inversion_reverses_the_order_and_flips_each(tmp_path):
    # the fork's headline case. The arc [3, 18) spans g2 and g3 with their spacers; reversed, the
    # two genes swap places *and* each points backwards — g1, g3-, g2-, g4, g5.
    chrom = _seed(tmp_path)
    chrom.invert(3, 15)
    assert _genes(chrom) == [(1, 1), (3, -1), (2, -1), (4, 1), (5, 1)]
    # ...and they sit where the reversal puts them, spacers included
    assert _layout(chrom) == [
        (0, 3, 1, 1), (3, 3, 0, -1), (6, 3, 3, -1), (9, 3, 0, -1), (12, 3, 2, -1),
        (15, 3, 0, -1), (18, 3, 4, 1), (21, 3, 0, 1), (24, 3, 5, 1), (27, 3, 0, 1)]


def test_inverting_the_whole_genome_reverses_every_gene(tmp_path):
    chrom = _seed(tmp_path)
    chrom.invert(0, 30)
    assert _genes(chrom) == [(5, -1), (4, -1), (3, -1), (2, -1), (1, -1)]
    assert chrom.length == 30


def test_a_wrapping_inversion_re_roots_the_circular_coordinates(tmp_path):
    # the fork flagged its wrapping cases specially, and they are where a convention shows itself.
    # The arc [27, 30) + [0, 3) is the last spacer followed by g1. Reversed it reads g1-, spacer -
    # and zombi2 re-roots the block list at the arc's start, so g1 stays at position 0 rather than
    # moving to 27. A circular chromosome has no privileged origin; only the cycle is meaningful.
    chrom = _seed(tmp_path)
    chrom.invert(27, 6)
    assert _layout(chrom)[:2] == [(0, 3, 1, -1), (3, 3, 0, -1)]
    assert _gene_cycle(chrom) == [(1, -1), (2, 1), (3, 1), (4, 1), (5, 1)]
    assert chrom.length == 30                       # and nothing is gained or lost by the wrap


def test_a_wrapping_inversion_can_carry_a_gene_across_the_origin(tmp_path):
    # the arc [24, 30) + [0, 3) covers g5, the last spacer and g1: both genes flip, and their order
    # around the circle swaps
    chrom = _seed(tmp_path)
    chrom.invert(24, 9)
    assert _gene_cycle(chrom) == [(1, -1), (5, -1), (2, 1), (3, 1), (4, 1)]


def test_two_inversions_in_a_row_compose(tmp_path):
    # the fork's multi-event cases (test_divisions' 2- and 3-inversion branches): each acts on the
    # genome the previous one left
    chrom = _seed(tmp_path)
    chrom.invert(3, 15)                             # g1, g3-, g2-, g4, g5
    chrom.invert(6, 3)                              # flip g3- back, alone in its arc
    assert _genes(chrom) == [(1, 1), (3, 1), (2, -1), (4, 1), (5, 1)]


def test_an_inversion_is_its_own_inverse(tmp_path):
    chrom = _seed(tmp_path)
    before = _layout(chrom)
    chrom.invert(6, 12)
    assert _layout(chrom) != before
    chrom.invert(6, 12)
    assert _layout(chrom) == before


# --- genes are indivisible ------------------------------------------------------------------------

def test_an_inversion_that_would_cut_a_gene_is_refused(tmp_path):
    # the fork only ever cut in intergenic space (make_inversion_*intergenic*); zombi2 states the
    # same rule as a guard — a breakpoint strictly inside a declared gene is not a legal cut
    chrom = _seed(tmp_path)
    with pytest.raises(_CutsGene):
        chrom.invert(4, 4)                          # ends at 8, inside g2 [6, 9)


@pytest.mark.parametrize("start, length", [(0, 3), (3, 3), (3, 6), (0, 30), (6, 12), (27, 6)])
def test_a_legal_inversion_never_splits_a_gene(tmp_path, start, length):
    chrom = _seed(tmp_path)
    chrom.invert(start, length)
    for b in chrom.blocks:
        if b.gene:
            assert b.end - b.start == 3, f"gene {b.gene} was split into a {b.end - b.start} bp piece"
    assert {b.gene for b in chrom.blocks if b.gene} == {1, 2, 3, 4, 5}


# --- divisions: the fork's natural_cuts ------------------------------------------------------------

def test_blocks_tile_the_chromosome_and_split_only_at_the_cuts(tmp_path):
    # the fork's test_divisions: an inversion subdivides the genome only where it cut, and the
    # pieces still tile it end to end
    chrom = _seed(tmp_path)
    assert len(chrom.blocks) == 10                  # 5 genes + 5 spacers, the natural boundaries
    chrom.invert(4 + 1, 4)                          # [5, 9): cuts the first spacer at 5, g2 at 9...

    at = 0
    for b in chrom.blocks:                          # ...and the pieces still tile the chromosome
        assert b.start < b.end
        at += b.end - b.start
    assert at == chrom.length == 30
    # one new boundary, mid-spacer at 5 — the gene boundaries were already there
    assert len(chrom.blocks) == 11
    assert sorted((b.end - b.start) for b in chrom.blocks if not b.gene) == [1, 2, 3, 3, 3, 3]


# --- the coordinate map: the fork's afterToBeforeT --------------------------------------------------

def test_trace_back_maps_every_position_home(tmp_path):
    # the fork asserted a coordinate map (Inversion.afterToBeforeT) after a scripted inversion:
    # which ancestral position does position x now hold? zombi2 answers with trace_back(). One case
    # is kept, to fix the convention — the oracle covers the arithmetic for arbitrary inputs.
    chrom = _seed(tmp_path)
    chrom.invert(3, 6)                              # spacer(3,6) + g2(6,9) reversed in place

    where = [pos for (_src, pos, _strand) in chrom.trace_back()]
    assert where[:3] == [0, 1, 2]                   # g1, untouched
    assert where[3:9] == [8, 7, 6, 5, 4, 3]         # the arc, read backwards
    assert where[9:] == list(range(9, 30))          # everything after the arc, untouched

    strands = [s for (_src, _pos, s) in chrom.trace_back()]
    assert strands[:3] == [1, 1, 1] and strands[3:9] == [-1] * 6 and strands[9:12] == [1, 1, 1]
    # the multiset of ancestral positions is conserved: an inversion moves material, never makes it
    assert sorted(where) == list(range(30))


# --- tandem duplication: the fork's test_tandemdup_* ----------------------------------------------

def _counter(start=100):
    """A copy-id minter, standing in for the engine's."""
    n = [start]

    def mint():
        n[0] += 1
        return n[0]
    return mint


def test_a_tandem_duplication_copies_a_gene_in_place(tmp_path):
    # the arc [3, 9) is spacer + g2. Duplicated in tandem the copy lands immediately after, so the
    # genome reads g1, spacer, g2, spacer, g2, then on as before — g2 twice, six bases longer.
    chrom = _seed(tmp_path)
    copied = chrom.duplicate(3, 6, _counter())
    assert chrom.length == 36
    assert _genes(chrom) == [(1, 1), (2, 1), (2, 1), (3, 1), (4, 1), (5, 1)]
    assert [(p, n, g) for (p, n, g, _s) in _layout(chrom)][:5] == [
        (0, 3, 1), (3, 3, 0), (6, 3, 2), (9, 3, 0), (12, 3, 2)]
    # the record names each block's parent copy and the fresh child it begot
    assert [(src, beg, end) for (_par, _child, src, beg, end) in copied] == [(0, 3, 6), (0, 6, 9)]
    # one fresh lineage per *copy lineage*, not per block: the seed genome is all one lineage, so
    # the two blocks of the arc beget a single child between them
    assert len({parent for (parent, *_rest) in copied}) == 1
    assert len({child for (_par, child, *_rest) in copied}) == 1


def test_a_duplicated_gene_keeps_its_orientation(tmp_path):
    # the contrast with inversion, which is the fork's whole point in separating the two: a
    # duplication moves nothing and turns nothing, it only adds
    chrom = _seed(tmp_path)
    chrom.duplicate(3, 6, _counter())
    assert all(s == 1 for (_p, _n, _g, s) in _layout(chrom))


def test_duplicating_a_multi_gene_arc_keeps_the_order(tmp_path):
    # [3, 18) spans g2 and g3; the copy repeats them in the same order, unlike an inversion
    chrom = _seed(tmp_path)
    chrom.duplicate(3, 15, _counter())
    assert _genes(chrom) == [(1, 1), (2, 1), (3, 1), (2, 1), (3, 1), (4, 1), (5, 1)]


def test_a_duplication_that_would_cut_a_gene_is_refused(tmp_path):
    chrom = _seed(tmp_path)
    with pytest.raises(_CutsGene):
        chrom.duplicate(4, 4, _counter())               # ends at 8, inside g2 [6, 9)


# --- loss: the fork's test_loss_* -----------------------------------------------------------------

def test_a_loss_removes_a_gene_and_closes_the_gap(tmp_path):
    # the arc [3, 9) is spacer + g2: deleting it leaves g1 abutting the spacer that followed g2
    chrom = _seed(tmp_path)
    lost = chrom.delete(3, 6)
    assert chrom.length == 24
    assert _genes(chrom) == [(1, 1), (3, 1), (4, 1), (5, 1)]
    assert [(src, beg, end) for (_cp, src, beg, end) in lost] == [(0, 3, 6), (0, 6, 9)]


def test_a_loss_inside_one_intergene_removes_no_gene(tmp_path):
    chrom = _seed(tmp_path)
    chrom.delete(3, 3)
    assert chrom.length == 27
    assert _genes(chrom) == [(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)]


def test_a_loss_that_would_take_the_last_gene_does_not_happen(tmp_path):
    # a chromosome never exists without a gene, so the whole-genome deletion is refused outright
    chrom = _seed(tmp_path)
    assert chrom.delete(0, 30) is None
    assert chrom.length == 30 and len(_genes(chrom)) == 5


def test_a_loss_may_take_every_gene_but_one(tmp_path):
    # ...but stripping it down to a single gene is allowed
    chrom = _seed(tmp_path)
    assert chrom.delete(3, 24) is not None               # [3, 27): g2..g5 and their spacers
    assert _genes(chrom) == [(1, 1)]
    assert chrom.length == 6


# --- transposition: the fork's test_transposition* and test_genomes.py's cut_and_paste ------------

def test_a_transposition_moves_a_gene_without_turning_it(tmp_path):
    # the fork's cut-and-paste. Lift [3, 9) (spacer + g2) out; the remainder is 24 bp reading
    # g1, spacer, g3, spacer, g4, spacer, g5, spacer. Drop the arc back at 12 — after g3's spacer —
    # and g2 sits between g3 and g4, still pointing forwards.
    chrom = _seed(tmp_path)
    assert chrom.transpose(3, 6, 12) is True
    assert chrom.length == 30
    assert _genes(chrom) == [(1, 1), (3, 1), (2, 1), (4, 1), (5, 1)]
    assert all(s == 1 for (_p, _n, g, s) in _layout(chrom) if g)


def test_a_flipped_transposition_turns_only_the_moved_block(tmp_path):
    chrom = _seed(tmp_path)
    chrom.transpose(3, 6, 12, flipped=True)
    assert _genes(chrom) == [(1, 1), (3, 1), (2, -1), (4, 1), (5, 1)]


def test_a_transposition_that_lands_inside_a_gene_leaves_the_genome_intact(tmp_path):
    # the rollback the engine relies on: a landing that would split a gene undoes the excision
    chrom = _seed(tmp_path)
    before = _layout(chrom)
    with pytest.raises(_CutsGene):
        chrom.transpose(3, 6, 13)                        # 13 is inside g3, which sits at [12, 15)
    assert _layout(chrom) == before                      # nothing lost, nothing moved


def test_transposing_a_whole_arc_back_to_where_it_was_is_a_no_op(tmp_path):
    chrom = _seed(tmp_path)
    before = _genes(chrom)
    chrom.transpose(3, 6, 3)
    assert _genes(chrom) == before


# --- origination: the fork's test_origination1 ----------------------------------------------------

def test_an_origination_lays_down_a_new_gene(tmp_path):
    # a de-novo gene of its own fresh source, inserted at a legal cut — indivisible from birth
    chrom = _seed(tmp_path)
    chrom.originate(9, 12, source=99, copy=500, family=42)
    assert chrom.length == 42
    assert _genes(chrom) == [(1, 1), (2, 1), (42, 1), (3, 1), (4, 1), (5, 1)]
    new = next(b for b in chrom.blocks if b.gene == 42)
    assert (new.source, new.start, new.end, new.strand) == (99, 0, 12, 1)


def test_an_origination_inside_a_gene_is_refused(tmp_path):
    chrom = _seed(tmp_path)
    with pytest.raises(_CutsGene):
        chrom.originate(7, 12, source=99, copy=500, family=42)   # 7 is inside g2 [6, 9)


def test_an_originated_gene_is_indivisible_like_a_declared_one(tmp_path):
    chrom = _seed(tmp_path)
    chrom.originate(9, 12, source=99, copy=500, family=42)
    with pytest.raises(_CutsGene):
        chrom.invert(10, 6)                              # would cut the new gene at [9, 21)
