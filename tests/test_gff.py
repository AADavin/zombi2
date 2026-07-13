"""Reading genome architecture (length + gene coordinates) from a GFF3 annotation.

A compact synthetic GFF exercises the parser (coordinate conversion, overlap trimming, seqid
selection, length fallbacks, gzip); an integration test feeds a GFF straight into the genic
nucleotide simulation.
"""

import gzip

import pytest

from zombi2 import BirthDeath, read_gff, simulate_nucleotide_genomes, simulate_species_tree
from zombi2.genomes.gff import GffGenome, read_gff_all

# chr1 (1000 bp, circular, 4 genes incl. an overlap + a nested one) and a smaller plasmid.
GFF = """\
##gff-version 3
##sequence-region chr1 1 1000
chr1\tx\tregion\t1\t1000\t.\t+\t.\tID=chr1;Is_circular=true
chr1\tx\tgene\t10\t50\t.\t+\t.\tID=gene-g1;locus_tag=g1;Name=alpha
chr1\tx\tCDS\t10\t50\t.\t+\t0\tID=cds1;Parent=gene-g1
chr1\tx\tgene\t45\t80\t.\t-\t.\tID=gene-g2;locus_tag=g2;Name=beta
chr1\tx\tgene\t100\t200\t.\t+\t.\tID=gene-g3;locus_tag=g3
chr1\tx\tgene\t120\t150\t.\t+\t.\tID=gene-g4;locus_tag=g4
##sequence-region plasmid1 1 300
plasmid1\tx\tregion\t1\t300\t.\t+\t.\tID=plasmid1;Is_circular=true
plasmid1\tx\tgene\t20\t60\t.\t+\t.\tID=gene-p1;locus_tag=p1
"""


def _write(tmp_path, text=GFF, name="genome.gff"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_reads_length_genes_and_circular(tmp_path):
    g = read_gff(_write(tmp_path))
    assert isinstance(g, GffGenome)
    assert g.seqid == "chr1" and g.length == 1000 and g.circular is True
    assert g.n_features == 4                         # 4 gene features on chr1 (CDS ignored)


def test_coordinates_are_zero_based_half_open(tmp_path):
    g = read_gff(_write(tmp_path))
    names = {n: (a, b) for a, b, n in g.genes}
    assert names["g1"] == (9, 50)                    # GFF [10,50] -> [9,50)
    assert names["g3"] == (99, 200)


def test_overlaps_trimmed_and_nested_dropped(tmp_path):
    g = read_gff(_write(tmp_path))
    got = {n: (a, b) for a, b, n in g.genes}
    assert got == {"g1": (9, 50), "g2": (50, 80), "g3": (99, 200)}  # g2 start clipped, g4 gone
    assert g.n_trimmed == 1 and g.n_dropped == 1
    for (a, b, _), (c, _d, _) in zip(g.genes, g.genes[1:]):         # disjoint, sorted
        assert c >= b


def test_seqid_selection(tmp_path):
    path = _write(tmp_path)
    assert read_gff(path).seqid == "chr1"            # default: most-annotated sequence
    p = read_gff(path, seqid="plasmid1")
    assert p.length == 300 and [n for _a, _b, n in p.genes] == ["p1"]
    with pytest.raises(ValueError, match="not found"):
        read_gff(path, seqid="chrZ")


def test_gzip_is_supported(tmp_path):
    raw = _write(tmp_path)
    gz = tmp_path / "genome.gff.gz"
    with gzip.open(gz, "wt") as f:
        f.write(GFF)
    a, b = read_gff(raw), read_gff(str(gz))
    assert (a.length, a.genes) == (b.length, b.genes)


def test_length_falls_back_to_region_then_max_end(tmp_path):
    # no ##sequence-region pragma -> length from the `region` feature
    no_pragma = "\n".join(l for l in GFF.splitlines() if not l.startswith("##sequence-region"))
    assert read_gff(_write(tmp_path, no_pragma, "a.gff")).length == 1000
    # neither pragma nor region -> largest gene end (chr1 max end = 200)
    bare = "\n".join(l for l in no_pragma.splitlines() if "\tregion\t" not in l)
    assert read_gff(_write(tmp_path, bare, "b.gff")).length == 200


def test_no_gene_features_raises(tmp_path):
    only_region = "##gff-version 3\nchr1\tx\tregion\t1\t100\t.\t+\t.\tID=chr1\n"
    with pytest.raises(ValueError, match="no gene"):
        read_gff(_write(tmp_path, only_region, "c.gff"))


def test_names_unique_even_when_annotation_repeats(tmp_path):
    dup = ("##gff-version 3\n##sequence-region c 1 100\n"
           "c\tx\tgene\t1\t10\t.\t+\t.\tName=dup\n"
           "c\tx\tgene\t20\t30\t.\t+\t.\tName=dup\n")
    g = read_gff(_write(tmp_path, dup, "d.gff"))
    assert len({n for _a, _b, n in g.genes}) == len(g.genes)   # disambiguated


def test_gff_feeds_the_genic_simulation(tmp_path):
    g = read_gff(_write(tmp_path))
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=5, age=1.0, seed=1)
    res = simulate_nucleotide_genomes(
        tree, inversion=0.003, loss=0.003, duplication=0.002, transfer=0.002,
        root_length=g.length, gene_intervals=g.genes, extension=0.95,
        pseudogenization=0.3, replacement=0.3, seed=2)
    # the GFF genes are recovered as gene blocks named by their locus tags
    gene_ids = {a.gene_id for a in res.gene_blocks()}
    assert {"g1", "g2", "g3"} <= gene_ids
    assert res.gene_trees() and res.intergene_trees()


# --- multi-sequence GFF -> multi-chromosome genome ------------------------------------

def test_read_gff_all_returns_every_sequence(tmp_path):
    gs = read_gff_all(_write(tmp_path))
    assert [g.seqid for g in gs] == ["chr1", "plasmid1"]   # most-annotated (the chromosome) first
    assert gs[0].length == 1000 and gs[1].length == 300
    assert [n for _a, _b, n in gs[1].genes] == ["p1"]      # the plasmid's gene survives (not dropped)


def test_multisequence_gff_feeds_a_multichromosome_sim(tmp_path):
    gs = read_gff_all(_write(tmp_path))
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=5, age=1.0, seed=1)
    res = simulate_nucleotide_genomes(
        tree, inversion=0.003, loss=0.002, extension=0.95,
        root_chromosomes=[(g.length, g.genes) for g in gs], seed=2)
    # every leaf descends from both replicons (no fission), each its own source
    assert all(len(leaf.chromosomes) == 2 for leaf in res.leaf_genomes.values())
    # genes from BOTH the chromosome and the plasmid are recovered (p1 would be dropped by read_gff)
    gene_ids = {a.gene_id for a in res.gene_blocks()}
    assert {"g1", "g3", "p1"} <= gene_ids
    assert len(res.block_gene_trees()) == len(res.blocks)
