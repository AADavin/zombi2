"""geneorder_events.tsv — the structural-event log with physical breakpoints (design phase 1).

Covers the two additive pieces: the event ``region`` now survives into the ``EventLog`` (for both
the nucleotide and ordered models), and the nucleotide ``--write geneorder`` serialises it as
``geneorder_events.tsv`` (native half-open coordinates). See docs/design/geneorder-export.md.
"""

from zombi2 import BirthDeath, OrderedGenome, Rates, simulate_genomes, simulate_species_tree
from zombi2.cli import main
from zombi2.genomes.events import EventType, Region
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes
from zombi2.genomes.simulation import GENEORDER_EVENTS_HEADER, geneorder_events_from_log

_GENES = "0\t50\tg1\n100\t150\tg2\n200\t250\tg3\n"


def _nuc(seed=7):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    return simulate_nucleotide_genomes(tree, output="genomes", root_length=300,
                                       inversion=0.03, transposition=0.02, loss=0.02, seed=seed)


# --- region persistence (the one-line genome_sim change) -------------------- #
def test_inversion_records_carry_region():
    res = _nuc()
    inv = [r for r in res.event_log if r.event is EventType.INVERSION]
    assert inv, "the seed should have produced inversions"
    for r in inv:
        assert isinstance(r.region, Region)
        assert r.region.length >= 1 and r.region.start >= 0


def test_ordered_model_also_carries_region():
    # the user asked for both models; the ordered genome fills Region the same way
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=2)
    rates = Rates(inversion=0.5, transposition=0.3, origination=0.0)
    g = simulate_genomes(tree, rates, initial_families=10, seed=3,
                         genome_factory=lambda ids: OrderedGenome(ids, extension=0.5))
    inv = [r for r in g.event_log if r.event is EventType.INVERSION]
    assert inv, "the seed should have produced inversions"
    assert all(r.region is not None for r in inv)


# --- phase 1.5: the paste / insert dest is captured too --------------------- #
def test_transposition_records_carry_paste_dest():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    res = simulate_nucleotide_genomes(tree, output="genomes", root_length=300,
                                      transposition=0.05, seed=4)
    tps = [r for r in res.event_log if r.event is EventType.TRANSPOSITION]
    assert tps, "the seed should have produced transpositions"
    assert all(isinstance(r.region.dest, int) for r in tps)  # paste position logged


def test_transfer_records_carry_recipient_dest():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=41)
    res = simulate_nucleotide_genomes(tree, output="genomes", inversion=0.004, transfer=0.006,
                                      loss=0.002, root_length=300, extension=0.95, seed=41)
    trs = [r for r in res.event_log if r.event is EventType.TRANSFER]
    assert trs, "the seed should have produced transfers"
    for r in trs:  # donor arc + integer recipient insert position (not a (chrom,pos) tuple)
        assert r.region is not None and isinstance(r.region.dest, int)
        assert r.donor and r.recipient


def test_novel_origination_carries_region_but_seed_does_not():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    res = simulate_nucleotide_genomes(tree, output="genomes", root_length=300,
                                      origination=0.5, inversion=0.02, seed=2)
    os_ = [r for r in res.event_log if r.event is EventType.ORIGINATION]
    novel = [r for r in os_ if r.region is not None]
    seed = [r for r in os_ if r.region is None]
    assert novel and seed  # both kinds occur
    assert all(r.region.length >= 1 for r in novel)  # insert at start, gene of length


# --- serialisation ---------------------------------------------------------- #
def test_writer_header_and_skips_speciation():
    text = geneorder_events_from_log(_nuc().event_log)
    lines = text.strip().splitlines()
    assert lines[0] == GENEORDER_EVENTS_HEADER
    for ln in lines[1:]:
        assert ln.split("\t")[1] not in ("S", "F")  # markers filtered out


# --- CLI ------------------------------------------------------------------- #
def _run_genomes(tmp_path, tag, seed):
    sp = tmp_path / f"S_{tag}"
    main(["species", "--birth", "1", "--death", "0.3", "--tips", "6", "--age", "3",
          "--seed", str(seed), "-o", str(sp)])
    genes = tmp_path / "genes.tsv"
    genes.write_text(_GENES)
    g = tmp_path / f"G_{tag}"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"), "--genome-model", "nucleotide",
               "--genes", str(genes), "--root-length", "300", "--inversion", "0.05",
               "--transposition", "0.03", "--write", "geneorder", "--seed", str(seed), "-o", str(g)])
    assert rc == 0
    return g / "geneorder_events.tsv"


def test_cli_writes_geneorder_with_populated_breakpoints(tmp_path):
    f = _run_genomes(tmp_path, "a", seed=5)
    assert f.exists()
    lines = f.read_text().strip().splitlines()
    assert lines[0] == GENEORDER_EVENTS_HEADER
    rows = [ln.split("\t") for ln in lines[1:]]
    invs = [r for r in rows if r[1] == "I"]
    assert invs, "expected inversion rows"
    for row in invs:  # chrom, start, length, strand populated for a structural inversion
        assert row[4] != "" and row[5] != "" and row[6] != ""
        assert int(row[6]) >= 1  # arc length
    tps = [r for r in rows if r[1] == "P"]
    assert tps, "expected transposition rows"
    for row in tps:  # dest (col 8) = paste position, populated for a transposition
        assert row[8] != "" and int(row[8]) >= 0


def test_geneorder_file_is_deterministic(tmp_path):
    assert _run_genomes(tmp_path, "x", seed=9).read_text() == \
           _run_genomes(tmp_path, "y", seed=9).read_text()
