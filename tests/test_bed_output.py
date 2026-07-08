"""BED gene annotations for the nucleotide genic model (``--write bed``).

``genes.bed`` is the root (seed) genome's gene annotation; ``BED/<node>.bed`` is every node's
annotation after rearrangements (genes at their coordinates on that node's chromosome). All are
standard BED6 (chrom, chromStart, chromEnd, name, score, strand), 0-based half-open. BED is
genic-only and drives the Python engine, so these tests do not need the Rust extension.
"""

from zombi2.cli import main
from zombi2.tree import read_newick


def _sp(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--tips", "6", "--seed", "1", "-o", str(sp)])
    return sp / "species_tree.nwk"


def _genes(tmp_path):
    p = tmp_path / "genes.tsv"
    p.write_text("100\t400\tgeneA\n600\t900\tgeneB\n1200\t1500\tgeneC\n")
    return p


def _bed_rows(path):
    return [ln.split("\t") for ln in path.read_text().strip().splitlines() if ln]


def _run_bed(tmp_path, out, extra=()):
    return main(["genomes", "-t", str(_sp(tmp_path)), "--genome-model", "nucleotide",
                 "--genes", str(_genes(tmp_path)), "--root-length", "2000",
                 "--transposition", "0.001", "--inversion", "0.001",
                 "--seed", "3", "--write", "bed", *extra, "-o", str(out)])


def test_genes_bed_is_valid_bed6_matching_seed(tmp_path):
    out = tmp_path / "g"
    assert _run_bed(tmp_path, out) == 0
    rows = _bed_rows(out / "genes.bed")
    assert [(r[3], int(r[1]), int(r[2])) for r in rows] == [
        ("geneA", 100, 400), ("geneB", 600, 900), ("geneC", 1200, 1500)]  # the seed intervals
    for chrom, start, end, name, score, strand in rows:
        assert chrom == "root_chromosome"
        assert int(start) < int(end)
        assert score == "0"
        assert strand in ("+", "-")


def test_per_node_bed_written_for_every_node(tmp_path):
    out = tmp_path / "g"
    assert _run_bed(tmp_path, out) == 0
    bdir = out / "BED"
    tree = read_newick((out / "species_tree.nwk").read_text())
    node_names = [n.name for n in tree.nodes_preorder()]
    assert node_names  # sanity
    for name in node_names:
        assert (bdir / f"{name}.bed").exists()
    # root.bed carries the same intervals as genes.bed (only the chromosome name differs)
    root_iv = [r[1:] for r in _bed_rows(bdir / "root.bed")]
    genes_iv = [r[1:] for r in _bed_rows(out / "genes.bed")]
    assert root_iv == genes_iv


def test_bed_coordinates_monotonic_and_named(tmp_path):
    out = tmp_path / "g"
    assert _run_bed(tmp_path, out) == 0
    seed_genes = {"geneA", "geneB", "geneC"}
    for bed in (out / "BED").glob("*.bed"):
        rows = _bed_rows(bed)
        starts = [int(r[1]) for r in rows]
        assert starts == sorted(starts)                    # features in genome order
        for a, b in zip(rows, rows[1:]):
            assert int(a[2]) <= int(b[1])                  # non-overlapping (intergenes between)
        assert {r[3] for r in rows} <= seed_genes          # every feature is a known gene


def test_bed_requires_gene_coordinates(tmp_path, capsys):
    rc = main(["genomes", "-t", str(_sp(tmp_path)), "--genome-model", "nucleotide",
               "--root-length", "500", "--write", "bed", "-o", str(tmp_path / "e")])
    assert rc == 1
    assert "needs gene coordinates" in capsys.readouterr().err


def test_bed_and_ancestral_together(tmp_path):
    out = tmp_path / "g"
    assert _run_bed(tmp_path, out, extra=("ancestral",)) == 0
    assert (out / "genes.bed").exists()
    assert (out / "BED" / "root.bed").exists()
    assert list((out / "Genomes").glob("*.fasta.gz"))     # ancestral genomes still produced
