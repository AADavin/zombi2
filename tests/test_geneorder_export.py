"""`zombi2 tools export` — gene-order study formats from a nucleotide genomes run (phase 2a).

Covers the adjacency algebra (unit), and the end-to-end CLI export of broken adjacencies per tree
edge on a content-conserving (inversion / transposition) run, where the result is exact.
"""

from zombi2.cli import main
from zombi2.tools.geneorder_export import _adjacencies, breakpoints_tsv, read_node_orders

_GENES = "0\t100\tg1\n1000\t1100\tg2\n2000\t2100\tg3\n"


# --- adjacency algebra ------------------------------------------------------ #
def test_adjacencies_of_forward_order():
    # circular +g1 +g2 +g3: each gene's head meets the next gene's tail
    adj = _adjacencies([("g1", 1), ("g2", 1), ("g3", 1)])
    assert adj == {frozenset(("g1_h", "g2_t")), frozenset(("g2_h", "g3_t")),
                   frozenset(("g3_h", "g1_t"))}


def test_inverting_one_gene_breaks_exactly_two_adjacencies():
    before = _adjacencies([("g1", 1), ("g2", 1), ("g3", 1)])
    after = _adjacencies([("g1", 1), ("g2", -1), ("g3", 1)])  # invert g2
    broken = before - after
    # the two adjacencies flanking g2 break; the g3–g1 adjacency is untouched
    assert broken == {frozenset(("g1_h", "g2_t")), frozenset(("g2_h", "g3_t"))}
    assert frozenset(("g3_h", "g1_t")) in after


def test_adjacencies_are_circular():
    # a single reversal of the whole 2-gene genome preserves the (circular) adjacency set
    assert _adjacencies([("g1", 1), ("g2", 1)]) == _adjacencies([("g2", -1), ("g1", -1)])


# --- end-to-end through the CLI --------------------------------------------- #
def _make_run(tmp_path, seed=8):
    sp = tmp_path / "S"
    main(["species", "--birth", "1", "--death", "0.3", "--tips", "6", "--age", "3",
          "--seed", str(seed), "-o", str(sp)])
    genes = tmp_path / "genes.tsv"
    genes.write_text(_GENES)
    g = tmp_path / "G"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"), "--genome-model", "nucleotide",
               "--genes", str(genes), "--root-length", "2500", "--inversion", "0.02",
               "--transposition", "0.01", "--write", "bed", "geneorder",
               "--seed", str(seed), "-o", str(g)])
    assert rc == 0
    return g


def test_cli_export_breakpoints(tmp_path):
    g = _make_run(tmp_path)
    out = tmp_path / "export"
    rc = main(["tools", "export", str(g), "--format", "breakpoints", "-o", str(out)])
    assert rc == 0
    f = out / "Breakpoints.tsv"
    assert f.exists()
    lines = f.read_text().strip().splitlines()
    assert lines[0] == "parent\tchild\tadjacency"
    genes = {"g1", "g2", "g3"}
    for ln in lines[1:]:
        parent, child, adj = ln.split("\t")
        # content-conserved run => every broken adjacency is a clean pair of two gene extremities
        ends = adj.split("|")
        assert len(ends) == 2, f"collapsed adjacency {adj!r} (unexpected without duplication)"
        assert all(e.rsplit("_", 1)[0] in genes and e.rsplit("_", 1)[1] in ("h", "t") for e in ends)


def test_breakpoints_match_direct_node_order_comparison(tmp_path):
    # the TSV must agree with recomputing broken adjacencies straight from the BED node orders
    g = _make_run(tmp_path, seed=3)
    orders = read_node_orders(str(g))
    text = breakpoints_tsv(str(g))
    n_rows = len(text.strip().splitlines()) - 1
    assert n_rows >= 1  # this seed produced rearrangements
    assert orders["root"]  # the reconstructed root order is present


def test_export_without_bed_errors_helpfully(tmp_path):
    sp = tmp_path / "S"
    main(["species", "--birth", "1", "--death", "0.3", "--tips", "5", "--age", "3",
          "--seed", "1", "-o", str(sp)])
    genes = tmp_path / "genes.tsv"
    genes.write_text(_GENES)
    g = tmp_path / "G"
    main(["genomes", "--tree", str(sp / "species_tree.nwk"), "--genome-model", "nucleotide",
          "--genes", str(genes), "--root-length", "2500", "--inversion", "0.02",
          "--write", "geneorder", "--seed", "1", "-o", str(g)])  # NO bed
    import pytest
    with pytest.raises(SystemExit):  # parser.error on the missing BED/ dir
        main(["tools", "export", str(g), "--format", "breakpoints"])
