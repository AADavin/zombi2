"""Tests for the GFF reader — declaring the genes of a seed nucleotide genome.

GFF is 1-based inclusive; blocks are 0-based half-open, so the reader converts on the way in. Only
``##sequence-region`` (a replicon's extent) and ``gene`` features are read; everything else is ignored.
"""

import pytest

from zombi2.genomes.gff import GffGene, read_gff


def _write(tmp_path, text):
    path = tmp_path / "genes.gff"
    path.write_text(text)
    return path


HEADER = "##gff-version 3\n##sequence-region chrom1 1 3000\n"


def test_reads_replicons_genes_and_converts_coordinates(tmp_path):
    gff = _write(tmp_path, HEADER +
                 "chrom1\tZOMBI2\tgene\t201\t500\t.\t+\t.\tID=dnaA\n"
                 "chrom1\tZOMBI2\tgene\t900\t1400\t.\t-\t.\tID=recA;Name=recombinase\n")
    lengths, genes = read_gff(gff)
    assert lengths == {"chrom1": 3000}
    # 1-based inclusive [201, 500] -> 0-based half-open [200, 500)
    assert genes == [GffGene("chrom1", 200, 500, 1, "dnaA"),
                     GffGene("chrom1", 899, 1400, -1, "recA")]


def test_only_gene_features_are_declared(tmp_path):
    gff = _write(tmp_path, HEADER +
                 "chrom1\tZOMBI2\tgene\t201\t500\t.\t+\t.\tID=dnaA\n"
                 "chrom1\tZOMBI2\tCDS\t201\t500\t.\t+\t0\tParent=dnaA\n"
                 "chrom1\tZOMBI2\texon\t201\t500\t.\t+\t.\tParent=dnaA\n")
    _lengths, genes = read_gff(gff)
    assert [g.name for g in genes] == ["dnaA"]               # CDS / exon are ignored


def test_several_replicons_and_deterministic_order(tmp_path):
    gff = _write(tmp_path, "##sequence-region b 1 500\n##sequence-region a 1 900\n"
                           "b\t.\tgene\t10\t20\t.\t+\t.\tID=g2\n"
                           "a\t.\tgene\t30\t40\t.\t+\t.\tID=g1\n")
    lengths, genes = read_gff(gff)
    assert lengths == {"a": 900, "b": 500}
    assert [(g.seqid, g.name) for g in genes] == [("a", "g1"), ("b", "g2")]   # sorted by replicon, start


def test_name_falls_back_from_ID_to_Name_to_coordinates(tmp_path):
    gff = _write(tmp_path, HEADER +
                 "chrom1\t.\tgene\t10\t20\t.\t+\t.\tName=only_name\n"
                 "chrom1\t.\tgene\t30\t40\t.\t+\t.\t\n")
    _lengths, genes = read_gff(gff)
    assert [g.name for g in genes] == ["only_name", "chrom1:30-40"]


def test_length_is_inferred_when_sequence_region_is_missing(tmp_path):
    gff = _write(tmp_path, "chrom1\t.\tgene\t10\t900\t.\t+\t.\tID=g\n")
    lengths, _genes = read_gff(gff)
    assert lengths == {"chrom1": 900}                        # the replicon ends at its last gene


def test_accepts_an_iterable_of_lines():
    lengths, genes = read_gff(["##sequence-region c 1 100", "c\t.\tgene\t1\t10\t.\t+\t.\tID=x"])
    assert lengths == {"c": 100} and genes[0].name == "x"


def test_rejects_overlapping_genes(tmp_path):
    gff = _write(tmp_path, HEADER +
                 "chrom1\t.\tgene\t100\t200\t.\t+\t.\tID=a\n"
                 "chrom1\t.\tgene\t150\t300\t.\t+\t.\tID=b\n")
    with pytest.raises(ValueError, match="overlap"):
        read_gff(gff)


def test_touching_genes_are_allowed(tmp_path):
    # laid end to end: [100,200) then [200,300) share a boundary but do not overlap
    gff = _write(tmp_path, HEADER +
                 "chrom1\t.\tgene\t101\t200\t.\t+\t.\tID=a\n"
                 "chrom1\t.\tgene\t201\t300\t.\t+\t.\tID=b\n")
    _lengths, genes = read_gff(gff)
    assert (genes[0].end, genes[1].start) == (200, 200)


def test_rejects_a_gene_beyond_its_replicon(tmp_path):
    gff = _write(tmp_path, "##sequence-region chrom1 1 100\n"
                           "chrom1\t.\tgene\t50\t200\t.\t+\t.\tID=a\n")
    with pytest.raises(ValueError, match="beyond replicon"):
        read_gff(gff)


@pytest.mark.parametrize("line, match", [
    ("chrom1\t.\tgene\tx\t500\t.\t+\t.\tID=a", "must be integers"),
    ("chrom1\t.\tgene\t500\t200\t.\t+\t.\tID=a", "start <= end"),
    ("chrom1\t.\tgene\t0\t200\t.\t+\t.\tID=a", "start <= end"),
    ("chrom1\t.\tgene\t1", "at least 8"),
])
def test_rejects_malformed_lines(tmp_path, line, match):
    with pytest.raises(ValueError, match=match):
        read_gff(_write(tmp_path, HEADER + line + "\n"))


def test_rejects_an_empty_declaration(tmp_path):
    with pytest.raises(ValueError, match="no replicon and no gene"):
        read_gff(_write(tmp_path, "##gff-version 3\n"))
