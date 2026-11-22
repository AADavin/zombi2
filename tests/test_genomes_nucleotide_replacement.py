"""Replacing transfers at the nucleotide resolution — ``replacement=`` on
`simulate_genomes_nucleotide`.

An additive transfer copies an arc of DNA into the recipient beside whatever it already has, so a
single-copy gene becomes two and the genome grows. A replacing one lands on the recipient's own copy
of that stretch and overwrites it, which is how bacteria exchange DNA with close relatives. Here it
is one event made of two the engine already had: a transfer placed where the homologue sits, and a
loss of that homologue. These check the rule that finds the homologue, that the genomes and the gene
trees still agree, and that the event is written and read as one.
"""

import collections

import pytest

from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.genomes.events import events_from_edges
from zombi2.genomes.nucleotide import (Block, Chromosome, NucleotideGenome, Transfer,
                                       _homologous_stretches, read_nucleotide_genomes)
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import jc69
from zombi2.species import simulate_species_tree

#: the run of the issue: ten single-copy genes of 1 kb in 20 kb, transfers to close relatives
ISSUE = dict(root_length=20000, genes=10, gene_length=1000, transfer=1.0, transfer_extent=1500,
             transfer_to="distance", seed=1)
#: rearranged genomes on a ring and a linear replicon, with duplicates and a donor that may be
#: its own recipient
BUSY = dict(chromosomes=[(40000, "circular"), (20000, "linear")], genes=15, gene_length=1000,
            transfer=3.0, transfer_extent=2500, transfer_to="distance", self_transfer=True,
            inversion=2.0, inversion_extent=6000, loss=0.3, loss_extent=1500, duplication=0.3,
            duplication_extent=1500, seed=2)


def _issue_tree():
    return simulate_species_tree(birth=1.0, death=0.0, n_extant=8, seed=1)


def _busy_tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=2)


def _walk(root):
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return out


def _shape(node):
    return (node.kind, node.species, round(node.time, 9),
            tuple(sorted(_shape(c) for c in node.children)))


def _copies_per_tip(g):
    out = []
    for leaf in sorted(g.complete_tree.extant_leaves()):
        out.append(collections.Counter(b.gene for c in g.node_genomes[leaf].chromosomes
                                       for b in c.blocks if b.is_gene))
    return out


def _trees_match_genomes(g) -> bool:
    for leaf in g.complete_tree.extant_leaves():
        blocks = [b for c in g.node_genomes[leaf].chromosomes for b in c.blocks]
        for trees, spans in ((g.gene_trees, g.gene_spans),
                             (g.block_trees, dict(enumerate(g.root_blocks)))):
            for key, tree in trees.items():
                src, lo, hi = spans[key]
                held = {b.copy for b in blocks if b.source == src and lo < b.end and b.start < hi}
                grown = sum(1 for n in _walk(tree.complete)
                            if n.kind == "extant" and n.species == leaf)
                if grown != len(held):
                    return False
    return True


# --- the run of the issue ---------------------------------------------------------------------------

def test_a_replacing_transfer_keeps_a_single_copy_gene_single_copy():
    sp = _issue_tree()
    added = simulate_genomes_nucleotide(sp, **ISSUE)
    replaced = simulate_genomes_nucleotide(sp, replacement=1.0, **ISSUE)
    # additive, as it always was: two tips hold two families twice and have grown by the arcs
    assert max(sum(v > 1 for v in c.values()) for c in _copies_per_tip(added)) == 2
    assert max(g.length for g in added.genomes.values()) > 20000
    # replacing: every tip still holds each of the ten genes once, in 20 kb
    assert all(c == {fam: 1 for fam in replaced.gene_spans} for c in _copies_per_tip(replaced))
    assert {g.length for g in replaced.genomes.values()} == {20000}
    transfers = [e for e in replaced.events if isinstance(e, Transfer)]
    assert transfers and all(e.replaced for e in transfers)
    for e in transfers:                       # what came in is as long as what it overwrote
        assert sum(b - a for (*_r, a, b) in e.transferred) == sum(b - a for (*_r, a, b) in e.replaced)
    assert replaced.replacements == {"found_resident": 2, "replaced": 2, "added_beside": 0}


def test_without_replacement_the_run_is_the_run_it_was():
    sp = _busy_tree()
    plain = simulate_genomes_nucleotide(sp, **BUSY)
    explicit = simulate_genomes_nucleotide(sp, replacement=0.0, **BUSY)
    assert plain.events == explicit.events and plain.rearrangements == explicit.rearrangements
    assert plain.replacements is None and "transfers" not in plain.summary()
    assert not any(e.replaced for e in plain.events if isinstance(e, Transfer))


# --- genomes, trees and the record still agree --------------------------------------------------------

@pytest.mark.parametrize("extra", [{}, dict(deletion=20.0, insertion=20.0, deletion_extent=10,
                                           insertion_extent=10)], ids=["no indels", "indels"])
@pytest.mark.parametrize("replacement", [1.0, 0.5])
def test_the_gene_trees_are_the_genomes_history(replacement, extra):
    g = simulate_genomes_nucleotide(_busy_tree(), replacement=replacement, **BUSY, **extra)
    assert sum(bool(e.replaced) for e in g.events if isinstance(e, Transfer)) > 5
    assert _trees_match_genomes(g)
    for node in g.node_genomes:               # and every node's genome can be put back together
        g.assembly(node)
    tally = g.replacements
    assert tally["found_resident"] == tally["replaced"] + tally["added_beside"]
    assert g.summary()["transfers"] == tally
    if replacement == 1.0:
        assert tally["added_beside"] == 0     # every homologue that was found was overwritten
    else:
        assert tally["replaced"] > 0 and tally["added_beside"] > 0


def test_the_overwritten_gene_dies_at_the_transfer_and_the_transfer_names_it():
    g = simulate_genomes_nucleotide(_busy_tree(), replacement=1.0, **BUSY)
    edges = g.genealogy
    died = {(e.family, e.copy): e for e in edges if e.kind == "loss"}
    named = [e for e in edges if e.kind == "transfer" and e.replaced is not None]
    assert len(named) > 10 and len(named) % 2 == 0          # two rows a transfer, both name it
    for e in named:
        loss = died[(e.family, e.replaced)]
        assert loss.time == e.time                            # one event, not two that coincide
        assert loss.lineage == next(x.recipient for x in edges
                                    if x.kind == "transfer" and x.parent == e.parent
                                    and x.recipient is not None)
    # the events table writes the pair as one row, in the family resolution's own word
    rows = collections.Counter(e.kind for e in events_from_edges(edges))
    assert rows["transfer_replacing"] == len(named) // 2 and rows["transfer_additive"] > 0
    # and the overwritten gene has no loss row of its own there: it is that row's second parent
    assert rows["loss"] == len(died) - len(named) // 2


def test_a_replacing_transfer_is_written_and_read_back_as_one_event(tmp_path):
    g = simulate_genomes_nucleotide(_busy_tree(), replacement=1.0, **BUSY)
    g.write(tmp_path)
    kinds = collections.Counter(line.split("\t")[1] for line in
                                (tmp_path / "block_events.tsv").read_text().splitlines()[1:])
    assert kinds["transfer_replacing"] > 0 and kinds["transfer"] > 0
    table = (tmp_path / "genome_events.tsv").read_text()
    assert "transfer_replacing" in table and "transfer_additive" in table
    back = read_nucleotide_genomes(tmp_path, g.complete_tree)
    assert [e for e in back.events if isinstance(e, Transfer)] == \
           [e for e in g.events if isinstance(e, Transfer)]
    assert set(back.gene_trees) == set(g.gene_trees)
    for fam, tree in g.gene_trees.items():
        assert _shape(back.gene_trees[fam].complete) == _shape(tree.complete)
    assert _trees_match_genomes(back)


def test_a_sequence_run_follows_the_replaced_stretch():
    # at rate zero every genome is the starting DNA, moved about by its own history. A replaced
    # stretch is the donor's copy, so a wrong tree there would put the wrong letters in the genome.
    g = simulate_genomes_nucleotide(_issue_tree(), replacement=1.0, **ISSUE)
    r = simulate_sequences(g, model=jc69(), substitution=0.0, seed=1)
    start = r.initial_genome[g.initial_genome.chromosomes[0].id]
    names = g.complete_tree.labels()
    complement = str.maketrans("ACGT", "TGCA")
    for node, genome in g.node_genomes.items():
        for chrom in genome.chromosomes:
            expected = "".join(start[b.start:b.end] if b.strand == 1
                               else start[b.start:b.end].translate(complement)[::-1]
                               for b in chrom.blocks)
            assert r.node_genomes[names[node]][chrom.id] == expected


# --- the rule that finds the homologue ----------------------------------------------------------------
#
# Source 0 is 1000 bases: spacer [0, 200), gene 1 [200, 500), spacer [500, 1000). The donor's arc is
# [100, 600): spacer, the gene whole, spacer.

def _arc(strand=1, copy=9):
    blocks = [Block(0, 100, 200, 1, copy), Block(0, 200, 500, 1, copy, 1), Block(0, 500, 600, 1, copy)]
    if strand == -1:
        blocks = [Block(b.source, b.start, b.end, -1, b.copy, b.gene) for b in reversed(blocks)]
    return blocks


def _genome(blocks, topology="linear"):
    return NucleotideGenome([Chromosome(1, topology, blocks)])


def test_the_homologue_lies_between_the_arcs_two_end_bases():
    recipient = _genome([Block(0, 0, 200, 1, 3), Block(0, 200, 500, 1, 3, 1), Block(0, 500, 1000, 1, 3)])
    (chrom, start, length, flipped), = _homologous_stretches(recipient, _arc())
    assert (start, length, flipped) == (100, 500, False)
    # the same stretch whichever way the donor happened to hold it
    (_c, start, length, flipped), = _homologous_stretches(recipient, _arc(strand=-1))
    assert (start, length, flipped) == (100, 500, True)


def test_a_recipient_that_holds_the_stretch_reversed_gives_it_flipped():
    recipient = _genome([Block(0, 500, 1000, -1, 3), Block(0, 200, 500, -1, 3, 1), Block(0, 0, 200, -1, 3)])
    (_c, start, length, flipped), = _homologous_stretches(recipient, _arc())
    assert (start, length, flipped) == (400, 500, True)


@pytest.mark.parametrize("blocks, why", [
    ([Block(0, 0, 150, 1, 3), Block(0, 700, 1000, 1, 3)], "the region was lost"),
    ([Block(0, 0, 200, 1, 3), Block(0, 200, 500, 1, 3, 1), Block(0, 500, 550, 1, 3),
      Block(0, 700, 1000, 1, 3)], "the last base was lost"),
    ([Block(0, 500, 1000, 1, 3), Block(0, 0, 200, 1, 3), Block(0, 200, 500, 1, 3, 1)],
     "the two ends are out of order"),
    ([Block(0, 0, 200, 1, 3), Block(5, 0, 300, 1, 4, 7), Block(0, 200, 500, 1, 3, 1),
      Block(0, 500, 1000, 1, 3)], "another gene lies between them"),
])
def test_a_transfer_with_no_homologue_to_land_on_is_additive(blocks, why):
    assert _homologous_stretches(_genome(blocks), _arc()) == [], why


def test_spacer_an_indel_put_there_does_not_stand_in_the_way():
    blocks = [Block(0, 0, 200, 1, 3), Block(8, 0, 12, 1, 6), Block(0, 200, 500, 1, 3, 1),
              Block(0, 500, 1000, 1, 3)]
    assert _homologous_stretches(_genome(blocks), _arc()) == []              # foreign, for all it knows
    (_c, start, length, _f), = _homologous_stretches(_genome(blocks), _arc(), inserted={8})
    assert (start, length) == (100, 512)                                     # and the insert goes too


def test_a_region_held_twice_gives_two_stretches_not_one_over_both():
    twice = [Block(0, 0, 200, 1, 3), Block(0, 200, 500, 1, 3, 1), Block(0, 500, 1000, 1, 3),
             Block(0, 0, 200, 1, 4), Block(0, 200, 500, 1, 4, 1), Block(0, 500, 1000, 1, 4)]
    found = _homologous_stretches(_genome(twice), _arc())
    assert sorted((start, length) for _c, start, length, _f in found) == [(100, 500), (1100, 500)]


def test_a_self_transfer_never_lands_on_the_stretch_it_copied():
    arc = _arc(copy=3)
    alone = _genome([Block(0, 0, 100, 1, 3), *arc, Block(0, 600, 1000, 1, 3)])
    assert _homologous_stretches(alone, arc) == []
    other = [Block(0, 0, 200, 1, 4), Block(0, 200, 500, 1, 4, 1), Block(0, 500, 1000, 1, 4)]
    beside = _genome([Block(0, 0, 100, 1, 3), *arc, Block(0, 600, 1000, 1, 3), *other])
    (_c, start, length, _f), = _homologous_stretches(beside, arc)
    assert (start, length) == (1100, 500)                 # the other copy of the region


def test_a_stretch_may_wrap_the_origin_of_a_ring():
    # the ring's origin sits inside the gene's neighbourhood: [500, 1000) then [0, 500)
    ring = _genome([Block(0, 500, 1000, 1, 3), Block(0, 0, 200, 1, 3), Block(0, 200, 500, 1, 3, 1)],
                   topology="circular")
    (_c, start, length, flipped), = _homologous_stretches(ring, _arc())
    assert (start, length, flipped) == (600, 500, False)          # 600..1000, then 0..100


# --- the argument and the command line -----------------------------------------------------------------

def test_replacement_is_a_probability():
    sp = _issue_tree()
    for bad in (1.5, -0.1, "often"):
        with pytest.raises(ValueError, match=r"probability in \[0, 1\]"):
            simulate_genomes_nucleotide(sp, replacement=bad, **ISSUE)
    assert simulate_genomes_nucleotide(sp, replacement=True, **ISSUE).replacements["replaced"] == 2


def test_the_command_line_takes_replacement_at_the_nucleotide_resolution(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--death", "0", "--n-extant", "8",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run), "--resolution", "nucleotide", "--root-length", "20000",
                 "--genes", "10", "--gene-length", "1000", "--transfer", "1.0",
                 "--transfer-extent", "1500", "--transfer-to", "distance", "--replacement",
                 "--seed", "1", "--quiet"]) == 0
    table = (run / "genomes" / "genome_events.tsv").read_text(encoding="utf-8")
    assert "transfer_replacing" in table and "transfer_additive" not in table
    summary = (run / "genomes" / "genome_summary.json").read_text(encoding="utf-8")
    assert '"replaced": 2' in summary
