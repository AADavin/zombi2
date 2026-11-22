"""Overlapping genes at the nucleotide resolution — kept, on one indivisible block.

A real annotation overlaps: a base or four where genes abut in an operon, a gene nested inside
another on the opposite strand. The engine has one rule, *a breakpoint never falls strictly inside a
gene*, and inside a run of overlapping genes every position is inside one. So the run is one block,
never cut, and the genes on it are never separated. These check that, and the bookkeeping around it:
each gene keeps its own coordinates, its own row in a GFF and its own tree.
"""

import collections
import random

import pytest

from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.genomes.gff import read_gff
from zombi2.genomes.nucleotide import read_nucleotide_genomes
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import jc69
from zombi2.species import simulate_species_tree

#: name: (start, end, strand), 1-based inclusive as a GFF writes it
GENES = {
    "a": (101, 400, "+"),      # a, b and c chain: b starts 4 bases inside a,
    "b": (397, 700, "+"),      # and c 1 base inside b, on the other strand
    "c": (700, 900, "-"),
    "d": (1201, 1500, "+"),    # on its own
    "e": (1800, 2400, "-"),    # f is nested inside e, on the other strand
    "f": (1900, 2000, "+"),
    "g": (2401, 2600, "+"),    # touches e, which is not overlapping it
}
LENGTH = 3000
#: the blocks those genes make, 0-based half-open, with the genes on each
BLOCKS = {(100, 900): ["a", "b", "c"], (1200, 1500): ["d"], (1799, 2400): ["e", "f"],
          (2400, 2600): ["g"]}
DNA = "".join(random.Random(7).choice("ACGT") for _ in range(LENGTH))


@pytest.fixture
def files(tmp_path):
    gff = tmp_path / "genes.gff"
    gff.write_text(f"##gff-version 3\n##sequence-region c 1 {LENGTH}\n" + "".join(
        f"c\tt\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={name}\n"
        for name, (s, e, strand) in GENES.items()), encoding="utf-8")
    fasta = tmp_path / "genes.fasta"
    fasta.write_text(f">c\n{DNA}\n", encoding="utf-8")
    return gff, fasta


def _run(files, *, seed=3, n_extant=8, **params):
    gff, fasta = files
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=seed)
    return simulate_genomes_nucleotide(sp, gff=gff, fasta=fasta, seed=seed, **params)


#: events of every kind that moves, copies or removes an arc, long enough to take whole blocks
BUSY = dict(inversion=6.0, inversion_extent=900, transposition=3.0, transposition_extent=900,
            duplication=2.0, duplication_extent=1200, loss=2.0, loss_extent=1200,
            transfer=2.0, transfer_extent=1200)


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


def _gff_rows(text):
    """``[(start, end, strand, attributes)]`` of a GFF's gene rows."""
    rows = []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) == 9:
            rows.append((int(cols[3]), int(cols[4]), cols[6],
                         dict(x.split("=", 1) for x in cols[8].split(";"))))
    return rows


def _listed_and_grown(g, leaf):
    """For one tip: how many copies of each family its GFF lists, and how many extant leaves each
    family's gene tree has there. The two must agree."""
    copies = collections.defaultdict(set)
    for _s, _e, _strand, attrs in _gff_rows(g._gff("x", g.node_genomes[leaf])):
        copies[int(attrs["family"])].add(attrs["copy"])
    listed = {fam: len(c) for fam, c in copies.items()}
    grown = {}
    for fam, tree in g.gene_trees.items():
        n = sum(1 for node in _walk(tree.complete)
                if node.kind == "extant" and node.species == leaf)
        if n:
            grown[fam] = n
    return listed, grown


def test_read_gff_keeps_overlapping_genes(files):
    _lengths, genes = read_gff(files[0])
    assert {g.name: (g.start + 1, g.end, "+" if g.strand == 1 else "-") for g in genes} == GENES
    _lengths, trimmed = read_gff(files[0], trim_overlaps=True)     # the old way, still on request
    assert [g.name for g in trimmed] == ["a", "b", "c", "d", "e", "g"]     # f was swallowed whole
    assert all(x.end <= y.start for x, y in zip(trimmed, trimmed[1:]))


def test_overlapping_genes_are_laid_down_as_one_block(files):
    g = _run(files)
    family = g.gene_names
    blocks = [(b.start, b.end) for b in g.initial_genome.chromosomes[0].blocks if b.is_gene]
    assert blocks == list(BLOCKS)
    # every gene keeps its own coordinates and strand, whatever block it sits on
    assert {name: g.gene_spans[fam][1:] for name, fam in family.items()} == \
           {name: (s - 1, e) for name, (s, e, _strand) in GENES.items()}
    assert {name: g.gene_strands[fam] for name, fam in family.items()} == \
           {name: 1 if strand == "+" else -1 for name, (_s, _e, strand) in GENES.items()}
    # and the spacer is everything else: gene blocks and spacer tile the replicon
    assert sum(b.length for b in g.initial_genome.chromosomes[0].blocks) == LENGTH


def test_the_initial_gff_is_the_annotation_that_was_given(files):
    g = _run(files)
    rows = _gff_rows(g._gff("initial", g.initial_genome))
    assert {attrs["Name"]: (s, e, strand) for s, e, strand, attrs in rows} == GENES
    assert [attrs["Name"] for *_r, attrs in rows] == sorted(GENES, key=lambda n: GENES[n][0])


def test_describe_names_every_gene_on_a_block(files):
    g = _run(files)
    text = g.describe(g.complete_tree.root)
    for names in BLOCKS.values():
        assert any(line.endswith("  " + ", ".join(names)) for line in text.splitlines())


def test_no_event_cuts_a_block_of_overlapping_genes(files):
    # without indels a gene block is whole wherever it turns up: at every node, every block that
    # holds genes is one of the four that were laid down, never a part of one
    g = _run(files, **BUSY)
    assert len(g.rearrangements) > 20 and len(g.events) > 20
    seen = collections.Counter()
    for genome in g.node_genomes.values():
        for chrom in genome.chromosomes:
            for b in chrom.blocks:
                if b.is_gene:
                    assert (b.start, b.end) in BLOCKS
                    seen[(b.start, b.end)] += 1
    assert set(seen) == set(BLOCKS)


def test_genes_sharing_a_block_share_a_history(files):
    # an event takes a block whole, so without indels the genes on it have the same tree — one tree
    # per gene, each under its own family
    g = _run(files, **BUSY)
    trees = g.gene_trees
    assert set(trees) == set(g.gene_names.values())
    for names in BLOCKS.values():
        shapes = {_shape(trees[g.gene_names[n]].complete) for n in names}
        assert len(shapes) == 1
    assert len({_shape(t.complete) for t in trees.values()}) > 1      # and the blocks do differ
    for leaf in g.complete_tree.extant_leaves():
        listed, grown = _listed_and_grown(g, leaf)
        assert listed == grown


def test_a_moved_block_writes_each_gene_where_it_now_sits(files):
    # every row of every tip's GFF, read back through the block it sits on: the stretch it names is
    # that gene's own interval, and it reads on the gene's strand turned by the block's
    g = _run(files, **BUSY)
    flipped = 0
    for leaf in g.complete_tree.extant_leaves():
        genome = g.node_genomes[leaf]
        trace = {c.id: g.trace_back(leaf)[c.id] for c in genome.chromosomes}
        chrom_of = {f"x_chr{c.id}": c.id for c in genome.chromosomes}
        for line in g._gff("x", genome).splitlines():
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            attrs = dict(x.split("=", 1) for x in cols[8].split(";"))
            _src, lo, hi = g.gene_spans[int(attrs["family"])]
            here = trace[chrom_of[cols[0]]][int(cols[3]) - 1:int(cols[4])]
            positions = [p for (_s, p, _strand) in here]
            turned = here[0][2]
            assert sorted(positions) == list(range(lo, hi))
            assert positions == (list(range(lo, hi)) if turned == 1 else list(range(hi - 1, lo - 1, -1)))
            coding = g.gene_strands[int(attrs["family"])]
            assert cols[6] == ("+" if coding * turned == 1 else "-")
            flipped += turned == -1
    assert flipped > 0


def test_a_deletion_can_end_one_gene_and_spare_its_neighbour(files):
    # with indels a block is no longer taken whole: a deletion can remove one gene from it and leave
    # the gene beside it. So the two trees part ways, and each tip's GFF lists exactly the genes
    # whose trees reach it.
    parted = 0
    for seed in range(1, 9):
        g = _run(files, seed=seed, loss=1.0, loss_extent=900, inversion=3.0, inversion_extent=900,
                 deletion=60.0, deletion_extent=150)
        for leaf in g.complete_tree.extant_leaves():
            listed, grown = _listed_and_grown(g, leaf)
            assert listed == grown
        for names in BLOCKS.values():
            fams = [g.gene_names[n] for n in names]
            parted += len({_shape(g.gene_trees[f].complete) if f in g.gene_trees else None
                           for f in fams}) > 1
    assert parted > 0


def test_a_run_read_back_from_its_files_groups_the_same_way(files, tmp_path):
    # genes.tsv holds each gene's own interval and nothing about blocks: overlapping intervals share
    # one, so the grouping is read off the file
    g = _run(files, **BUSY)
    g.write(tmp_path / "run")
    back = read_nucleotide_genomes(tmp_path / "run", g.complete_tree)
    assert back.gene_spans == g.gene_spans and back.gene_strands == g.gene_strands
    assert set(back.gene_trees) == set(g.gene_trees)
    for fam, tree in g.gene_trees.items():
        assert _shape(back.gene_trees[fam].complete) == _shape(tree.complete)
    for node in g.node_genomes:
        assert back.describe(node) == g.describe(node)
        assert back._gff("x", back.node_genomes[node]) == g._gff("x", g.node_genomes[node])


def test_a_gene_is_a_column_range_of_its_block(files):
    g = _run(files, **BUSY)
    r = simulate_sequences(g, model=jc69(), substitution=0.0, seed=1)
    for name, (start, end, _strand) in GENES.items():
        fam = g.gene_names[name]
        block, lo, hi = g.gene_span_in_block(fam)
        assert block == g.block_of(fam)
        assert r.founding[block][lo:hi] == DNA[start - 1:end]
    # genes on one block are slices of one alignment
    assert g.block_of(g.gene_names["a"]) == g.block_of(g.gene_names["c"])
    assert g.block_of(g.gene_names["e"]) != g.block_of(g.gene_names["g"])


def test_a_block_of_overlapping_genes_evolves_under_the_gene_model(files):
    # spacer all but frozen, genes not: a block of several genes must change like a gene
    g = _run(files, n_extant=12)
    r = simulate_sequences(g, model=jc69(), intergene_model=jc69(), intergene_speed=1e-9,
                           divergence=0.3, seed=1)
    genic = {g.block_of(fam) for fam in g.gene_spans}
    assert len(genic) == len(BLOCKS)
    for i in range(len(g.root_blocks)):
        varies = len(set(r.alignments[i].values())) > 1
        assert varies == (i in genic)


def test_the_command_line_keeps_overlapping_genes(files, tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--death", "0", "--n-extant", "4",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run), "--resolution", "nucleotide", "--gff", str(files[0]),
                 "--inversion", "2", "--inversion-extent", "900", "--seed", "1", "--quiet"]) == 0
    rows = _gff_rows((run / "genomes" / "gff" / "genome_initial.gff").read_text(encoding="utf-8"))
    assert {attrs["Name"]: (s, e, strand) for s, e, strand, attrs in rows} == GENES
