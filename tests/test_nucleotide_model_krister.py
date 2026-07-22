"""The nucleotide genome model, tested with Krister Swenson's suite.

**These tests were designed by Krister Swenson for ZOMBI1** — for his fork of it,
`thekswenson/Zombi` (`root_genome` branch, `tests/`), built for gene-order work — and are ported
here at his suggestion. They are kept together, grouped by where they came from rather than by what
they test, so this file reads against his originals and can be handed to him whole. He wrote the
scenarios; the coordinates are re-derived, because ZOMBI2 is a rewrite rather than a descendant and
the two models express the same ideas differently (see the concept map below).

`docs/design/swenson-test-port.md` tabulates each of his test files, where it landed, and what was
deliberately not ported.

His argument, which is the one that shapes all of this: **ZOMBI2's output is a set of files, and
the files are what a reader consumes.** Verifying the in-memory structures says nothing about
whether the written coordinates mean what they claim — a convention slip yields a plausible file
that replays to the wrong genome, and every in-memory test still passes.

Four sections, in the order they were built:

1. **Worked examples** (his `test_events.py`, `test_genomes.py`, `test_divisions.py`) — hand-derived
   gene-order outcomes on his own tiny genomes, run his way: every rate zero, then events applied by
   hand. The oracle in `test_genomes_nucleotide` already stress-tests this geometry in bulk; what it
   leaves untested is the gene annotation, which is what these cover.
2. **File replay** (his `test_geneorder_events.py`) — read the written event files back, replay them
   forward with an independent implementation, and check they reproduce the written genomes. His
   crown jewel, and the reason three output features got built.
3. **Pipeline determinism** (his `test_randomization.py`) — same seed, same bytes, across the whole
   pipeline rather than one command.
4. **The golden pin** — not his, but built to make section 1 possible. Splitting the engine's
   mutators into a choosing half and an applying half is what let the worked examples script a
   duplication or a loss at all, and the whole risk of that refactor was disturbing the rng draw
   order. This holds four seeded runs against a stored fixture.

Concept map — fork → ZOMBI2: *division* → `Block` (a run of one unbroken ancestry); *pieces* → the
ordered block list; `ch.genes` → the gene-carrying blocks in physical order; `natural_cuts` → block
boundaries; `Inversion.afterToBeforeT` → `trace_back()`; `Geneorder_events_per_branch/` →
`genome_event_positions.tsv` + `rearrangements.tsv`; `All_genomes/` → `gene_order.tsv` / `blocks.tsv`.

Most of this exercises the **nucleotide** resolution, which is where ZOMBI2 keeps genes, spacer and
block ancestry — the model his gene-order tests are about. Section 3 is broader: it compares the
whole CLI pipeline, because that is what his `test_randomization.py` did.

Regenerate the golden fixture (deliberately — never to make a red test green):

    python tests/test_nucleotide_model_krister.py
"""

import pathlib

import pytest

from zombi2.genomes.events import node_from_label
from zombi2.cli.main import main
from zombi2.genomes import simulate_genomes_nucleotide, simulate_genomes_ordered
from zombi2.genomes.nucleotide import _CutsGene
from zombi2.species import read_newick, simulate_species_tree


# ==================================================================================================
# 1. WORKED EXAMPLES
#
# Krister Swenson's, from ZOMBI1: test_events.py, test_genomes.py, test_divisions.py.
# He hand-verified each outcome on two tiny genomes of his own (30_6.gff, 30_10.gff),
# which are rebuilt here as GFF fixtures. The scenarios are his; every expected value is
# re-derived in ZOMBI2's coordinates rather than transcribed, since the two models
# parameterise an event differently.
# ==================================================================================================

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

# ==================================================================================================
# 2. FILE REPLAY
#
# Krister Swenson's, from ZOMBI1: test_geneorder_events.py (checkEventsAgainstGenomes).
# The written event files must replay to the written genomes. The best test in his suite,
# and the one that required ZOMBI2 to grow three output features before it could run.
# ==================================================================================================

# --------------------------------------------------------------------------- #
# Reading the written files
# --------------------------------------------------------------------------- #
def _read_tsv(path):
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, row.split("\t"))) for row in lines[1:] if row]


def _read_gene_order(path):
    """``gene_order.tsv`` -> ``{node: [(family, strand), ...]}`` in genome order."""
    genomes = {}
    for r in _read_tsv(path / "gene_order.tsv"):
        genomes.setdefault(node_from_label(r["lineage"]), []).append((int(r["family"]), int(r["strand"])))
    return genomes


def _read_steps(path):
    """Every written event that moves genes, merged into one time-ordered stream.

    Rows sharing a timestamp keep the order they were written (a replacing transfer's displacements
    precede its arrival), so the sort key carries each row's index within its own file.
    """
    steps = []
    for i, r in enumerate(_read_tsv(path / "genome_event_positions.tsv")):
        steps.append((float(r["time"]), 0, i, r))
    for i, r in enumerate(_read_tsv(path / "rearrangements.tsv")):
        steps.append((float(r["time"]), 1, i, r))
    steps.sort(key=lambda s: s[:3])
    return [r for *_, r in steps]


# --------------------------------------------------------------------------- #
# The replay — an independent implementation of each operation
# --------------------------------------------------------------------------- #
def _flip(segment):
    return [(fam, -strand) for fam, strand in reversed(segment)]


def _replay(steps, tree):
    """Run the written history forward from an empty root genome, returning ``{node: genome}``.

    Speciations are interleaved by the tree's own timing: when a branch ends, both daughters start
    from a copy of what it had.
    """
    live = {tree.root: []}
    ended = {}
    in_flight = {}                    # a transfer's donor row, waiting for its recipient row
    pending = sorted((n.end_time, n.id) for n in tree.nodes.values() if n.end_time is not None)
    p = 0

    def settle(upto):
        """Retire every branch that ends at or before ``upto``, seeding its daughters."""
        nonlocal p
        while p < len(pending) and pending[p][0] <= upto:
            _, nid = pending[p]
            p += 1
            genome = live.pop(nid)
            ended[nid] = genome
            for c in (tree.nodes[nid].children or ()):
                live[c] = list(genome)

    for r in steps:
        t, kind, lineage = float(r["time"]), r["kind"], node_from_label(r["lineage"])
        settle(t)
        # a branch that has already ended takes no more events; the engine never emits such a row
        assert lineage in live, f"event at {t} on branch {lineage}, which is not alive then"
        g = live[lineage]
        start, length = int(r["start"]), int(r["length"])

        if kind == "origination":
            g.insert(start, (int(r["family"]), +1))
        elif kind == "duplication":
            at = int(r["dest_position"])
            g[at:at] = g[start:start + length]
        elif kind == "loss":
            del g[start:start + length]
        elif kind == "transfer_donor":
            # the donor branch is unchanged; hold what left until its arriving row turns up
            in_flight[(t, node_from_label(r["donor"]), node_from_label(r["recipient"]))] = g[start:start + length]
        elif kind == "transfer_recipient":
            block = in_flight.pop((t, node_from_label(r["donor"]), node_from_label(r["recipient"])))
            assert len(block) == length, "the two rows of a transfer disagree on its length"
            g[start:start] = block                      # the block arrives at start
        elif kind == "inversion":
            g[start:start + length] = _flip(g[start:start + length])
        elif kind in ("transposition", "translocation"):
            segment = g[start:start + length]
            del g[start:start + length]                       # excised first, then placed
            at = int(r["dest_position"])
            g[at:at] = _flip(segment) if int(r["flipped"]) else segment
        else:
            raise AssertionError(f"unhandled event kind {kind!r}")

    settle(float("inf"))
    assert not in_flight, f"{len(in_flight)} transfer(s) left without an arriving row"
    return ended | live


def _positions_run(tmp_path, *, seed, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.4, n_extant=12, seed=seed)
    params = dict(duplication=0.4, loss=0.3, origination=0.5, transfer=0.3, inversion=0.5,
                  transposition=0.4, chromosomes=1, initial_families=10, seed=seed)
    params.update(kw)
    r = simulate_genomes_ordered(sp, **params)
    r.write(tmp_path, outputs=("events", "gene_order", "rearrangements", "event_positions"))
    return r


# --------------------------------------------------------------------------- #
# The test
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_written_positions_replay_to_the_written_genomes(tmp_path, seed):
    r = _positions_run(tmp_path, seed=seed)
    replayed = _replay(_read_steps(tmp_path), r.complete_tree)
    written = _read_gene_order(tmp_path)

    assert any(written.values()), "the fixture should produce genes to compare"
    for node, genome in written.items():
        assert replayed[node] == genome, f"replay diverges at node {node}"


def test_replacing_transfers_displace_before_they_arrive(tmp_path):
    # the ordering rule the file format promises: rows sharing a timestamp apply as written, so a
    # replacing transfer's losses land before its block does
    r = _positions_run(tmp_path, seed=7, replacement=True, transfer=0.8)
    # the fixture must actually displace something, or the ordering rule goes untested
    arrivals = [i for i, p in enumerate(r.event_positions) if p.kind == "transfer_recipient"]
    assert arrivals, "fixture produced no transfers"
    displaced = [i for i, p in enumerate(r.event_positions)
                 if p.kind == "loss" and any(p.time == r.event_positions[j].time and i < j
                                             for j in arrivals)]
    assert displaced, "fixture produced no replacement displacements"

    replayed = _replay(_read_steps(tmp_path), r.complete_tree)
    for node, genome in _read_gene_order(tmp_path).items():
        assert replayed[node] == genome, f"replay diverges at node {node}"


def test_the_replay_is_sensitive_to_a_one_position_slip(tmp_path):
    # a negative control: if the test could not tell a correct coordinate from a wrong one, passing
    # would mean nothing. Nudge one written start by a single position and the replay must break.
    r = _positions_run(tmp_path, seed=1)
    steps = _read_steps(tmp_path)
    victim = next(s for s in steps if s["kind"] == "loss" and int(s["start"]) > 0)
    victim["start"] = str(int(victim["start"]) - 1)

    replayed = _replay(steps, r.complete_tree)
    assert replayed != _read_gene_order(tmp_path)


def test_the_replay_is_sensitive_to_a_misplaced_arrival(tmp_path):
    # the same control for the half a transfer's second row is responsible for. Perturbing *one*
    # arrival need not show: genomes are compared as (family, strand), and a block landing between
    # two genes of the same family reads identically either side of the boundary. So perturb each
    # arrival in turn and require that the replay notices at least one.
    r = _positions_run(tmp_path, seed=1)
    written = _read_gene_order(tmp_path)
    n_arrivals = sum(1 for s in _read_steps(tmp_path) if s["kind"] == "transfer_recipient")
    assert n_arrivals, "fixture produced no transfers"

    noticed = 0
    for i in range(n_arrivals):
        steps = _read_steps(tmp_path)                      # a fresh copy for each perturbation
        victim = [s for s in steps if s["kind"] == "transfer_recipient"][i]
        victim["start"] = str(int(victim["start"]) + 1)    # one position further along
        try:
            noticed += _replay(steps, r.complete_tree) != written
        except AssertionError:                             # a start past the end is also "noticed"
            noticed += 1
    assert noticed, f"none of the {n_arrivals} arrivals mattered — the replay is not reading them"


def test_every_gene_content_event_has_a_position(tmp_path):
    # the table is total over the events that change gene content: no silent gaps, and speciation
    # (which copies a genome wholesale) is correctly absent
    r = _positions_run(tmp_path, seed=3)
    # a transfer writes one row per branch, so every genealogy row — donor side and recipient side —
    # finds its position under its own lineage once the two transfer kinds are folded together
    positioned = {(p.time, p.lineage, p.kind.split("_")[0]) for p in r.event_positions}

    kinds = set()
    for e in r.events:
        if e.kind == "speciation":            # a genome is copied wholesale: no position to record
            continue
        kinds.add(e.kind)
        assert (e.time, e.lineage, e.kind) in positioned, \
            f"{e.kind} at {e.time} on {e.lineage} has no position"
    assert kinds == {"origination", "duplication", "loss", "transfer"}, \
        f"the fixture should exercise every gene-content kind, got {sorted(kinds)}"

# ==================================================================================================
# 3. PIPELINE DETERMINISM
#
# Krister Swenson's, from ZOMBI1: test_randomization.py. Same seed, same output —
# widened here from one mode to the whole pipeline, compared byte for byte.
# ==================================================================================================

PIPELINE_SEED = "20"


def _pipeline(root):
    """Run species → genomes (all three resolutions) → sequences → traits into ``root``."""
    tree = str(root / "species_complete.nwk")
    assert main(["species", "--birth", "1.0", "--death", "0.3", "--n-extant", "12",
                 "--seed", PIPELINE_SEED, "-o", str(root), "--flat"]) == 0

    assert main(["genomes", "-t", tree, "--duplication", "0.3", "--loss", "0.25",
                 "--origination", "0.6", "--seed", PIPELINE_SEED, "-o", str(root / "g_unordered"), "--flat"]) == 0
    assert main(["genomes", "-t", tree, "--resolution", "ordered", "--duplication", "0.3",
                 "--loss", "0.25", "--origination", "0.6", "--transfer", "0.2",
                 "--inversion", "0.4", "--transposition", "0.3", "--chromosomes", "2",
                 "--seed", PIPELINE_SEED, "-o", str(root / "g_ordered"),
                 "--write", "events", "profiles", "gene_order", "rearrangements",
                 "chromosome_events", "event_positions", "--flat"]) == 0
    assert main(["genomes", "-t", tree, "--resolution", "nucleotide", "--root-length", "600",
                 "--genes", "4", "--inversion", "0.8", "--duplication", "0.4", "--loss", "0.3",
                 "--seed", PIPELINE_SEED, "-o", str(root / "g_nucleotide"),
                 "--write", "events", "genes", "blocks", "rearrangements", "--flat"]) == 0

    assert main(["sequences", "--genomes", str(root / "g_unordered"), "--model", "hky85",
                 "--length", "150", "--seed", PIPELINE_SEED, "-o", str(root / "s"),
                 "--write", "alignments", "phylograms", "ancestral",
                 "species_phylogram", "--flat"]) == 0
    assert main(["traits", "-t", tree, "--rate", "1.0", "--seed", PIPELINE_SEED,
                 "-o", str(root / "t"), "--flat"]) == 0


def _artifacts(root):
    """``{relative path: bytes}`` for everything written, bar the run logs (they stamp the wall
    clock, so they differ between two runs however deterministic the simulation)."""
    return {str(p.relative_to(root)): p.read_bytes()
            for p in sorted(root.rglob("*")) if p.is_file() and p.suffix != ".log"}


@pytest.fixture(scope="module")
def two_runs(tmp_path_factory):
    base = tmp_path_factory.mktemp("determinism")
    first, second = base / "run1", base / "run2"
    for root in (first, second):
        root.mkdir()
        _pipeline(root)
    return _artifacts(first), _artifacts(second)


def test_the_pipeline_writes_the_same_files_twice(two_runs):
    first, second = two_runs
    assert sorted(first) == sorted(second)
    assert first, "the pipeline should have written something to compare"


def test_every_written_byte_is_identical(two_runs):
    first, second = two_runs
    differing = [name for name in first if first[name] != second[name]]
    assert not differing, f"same seed, different output: {differing}"


def test_the_comparison_covers_every_level(two_runs):
    # a guard on the fixture rather than the engine: if a command stopped writing, or a level were
    # dropped from the pipeline above, the two runs would still agree — vacuously
    names = " ".join(two_runs[0])
    for expected in ("species_complete.nwk", "genome_events.tsv", "gene_order.tsv",
                     "genome_event_positions.tsv", "blocks.tsv", "genes.tsv",
                     "sequences_alignment", "trait_values.tsv"):
        assert expected in names, f"{expected} is missing — the pipeline is not covering that level"


def test_excluding_the_run_logs_is_justified_and_narrow(tmp_path):
    # the exclusion has to be justified, not assumed. A log carries a wall-clock `timestamp`, which
    # is why it cannot join a byte comparison; everything else in it must still agree between two
    # same-seed runs. (Two runs a second apart would show the timestamp differing too, but asserting
    # that would make this test depend on how fast the machine is.)
    runs = []
    for tag in ("a", "b"):
        root = tmp_path / tag
        root.mkdir()
        main(["species", "--birth", "1.0", "--death", "0.3", "--n-extant", "10", "--seed", PIPELINE_SEED,
              "-o", str(root), "--flat"])
        runs.append((root / "species.log").read_text().splitlines())

    first, second = runs
    assert any(line.startswith("timestamp\t") for line in first), \
        "no timestamp in the run log — the exclusion would be unnecessary"
    # `output` differs by construction: the two runs write to different directories
    differing = [(a, b) for a, b in zip(first, second, strict=True) if a != b]
    assert all(a.startswith(("timestamp\t", "output\t")) for a, _ in differing), \
        f"a run log differs by more than its timestamp and path: {differing}"

# ==================================================================================================
# 4. THE GOLDEN PIN
#
# NOT Krister's — this one is ZOMBI2's own. It is the safety net for splitting the
# nucleotide mutators into a choosing half and an applying half, which is what let
# section 1 script a duplication or a loss at all. Kept here because it exists only
# because of that port, and it is what makes those tests safe to trust.
# ==================================================================================================

GOLDEN = pathlib.Path(__file__).parent / "data" / "nucleotide_golden.txt"

#: a fixed tree, so the fixture does not depend on the species engine's draws either
GOLDEN_TREE = "(((A:0.4,B:0.4):0.3,C:0.7):0.3,(D:0.5,E:0.5):0.5);"

#: runs chosen to fire every mutator: content events, rearrangements, and the chromosome tier
GOLDEN_RUNS = {
    "content": dict(duplication=0.6, loss=0.5, origination=0.7, transfer=0.4, root_length=400,
                    genes=3),
    "rearrangements": dict(inversion=0.8, transposition=0.6, translocation=0.5, chromosomes=3,
                           inversion_probability=0.5, root_length=300, genes=2),
    "tier": dict(fission=2.0, fusion=0.4, chromosome_origination=0.3, chromosome_loss=0.2,
                 chromosomes=2, root_length=300, genes=2),
    "everything": dict(duplication=0.4, loss=0.3, origination=0.5, transfer=0.3, inversion=0.5,
                       transposition=0.4, translocation=0.3, fission=1.0, fusion=0.4,
                       chromosomes=3, root_length=400, genes=3),
}


def _render_golden(name, params):
    """One run, flattened to text: every node's blocks, then every event and rearrangement."""
    tree, _ = read_newick(GOLDEN_TREE)
    r = simulate_genomes_nucleotide(tree, seed=11, **params)
    lines = [f"## {name}"]
    for node in sorted(r.genomes):
        for chrom in r.genomes[node].chromosomes:
            blocks = " ".join(f"{b.source}:{b.start}-{b.end}:{b.strand}:{b.copy}:{b.gene}"
                              for b in chrom.blocks)
            lines.append(f"node {node} chrom {chrom.id} {chrom.topology} {blocks}")
    lines += [f"event {e}" for e in r.events]
    lines += [f"rearrangement {e}" for e in r.rearrangements]
    lines += [f"chromosome_event {e}" for e in r.chromosome_events]
    return lines


def _render_all_golden():
    out = []
    for name, params in GOLDEN_RUNS.items():
        out += _render_golden(name, params)
    return "\n".join(out) + "\n"


def test_seeded_runs_match_the_recorded_output():
    assert GOLDEN.exists(), f"missing fixture — regenerate with: python {__file__}"
    assert _render_all_golden() == GOLDEN.read_text(), (
        "a seeded nucleotide run no longer matches the recorded output. If the engine's sampling "
        "genuinely changed, regenerate the fixture; if not, an rng draw was added, removed or "
        "reordered.")


def test_the_fixture_covers_every_mutator():
    # a guard on the fixture: if the runs above stopped firing an event kind, the pin above would
    # still pass while protecting nothing
    text = GOLDEN.read_text()
    for kind in ("Duplication", "Loss", "Origination", "Transfer", "Inversion", "Transposition",
                 "Translocation"):
        assert kind in text, f"the fixture fires no {kind} — it is not pinning that mutator"
    for tier in ("fission", "fusion"):
        assert tier in text, f"the fixture fires no {tier}"


if __name__ == "__main__":
    GOLDEN.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN.write_text(_render_all_golden())
    print(f"wrote {GOLDEN} ({len(GOLDEN.read_text().splitlines())} lines)")
