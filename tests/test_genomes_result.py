"""Tests for the genome Result spine: phyletic profiles + write (zombi2.genomes)."""

from zombi2.species import simulate_species_tree
from zombi2.genomes import Profiles, simulate_genomes_unordered


def _run(seed=1, n_extant=12, death=0.3):
    sp = simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)
    g = simulate_genomes_unordered(sp, duplication=0.3, transfer=0.1, loss=0.2, origination=0.6,
                                   initial_families=8, seed=seed)
    return sp, g


# --- profiles: the families × extant-species copy-count matrix -------------

def test_profiles_columns_are_the_extant_tips():
    sp, g = _run(seed=2)
    assert set(g.profiles.species) == {n.id for n in sp.complete_tree.extant()}
    assert len(g.profiles.species) == sp.n_extant


def test_profiles_matrix_sum_is_the_extant_copy_total():
    # every copy at an extant tip is counted exactly once
    sp, g = _run(seed=3, death=0.5)
    total = sum(len(g.genomes[n.id]) for n in sp.complete_tree.extant())
    assert g.profiles.matrix.sum() == total


def test_profiles_column_matches_family_counts_at_that_tip():
    sp, g = _run(seed=4)
    p = g.profiles
    fi = {f: i for i, f in enumerate(p.families)}
    for j, s in enumerate(p.species):
        col = p.matrix[:, j]
        for family, count in g.family_counts(s).items():
            assert col[fi[family]] == count            # the cell is that family's copy count
        assert col.sum() == len(g.genomes[s])          # nothing else in the column


def test_profiles_presence_is_binary_and_tracks_the_matrix():
    _, g = _run(seed=5)
    p = g.profiles
    assert set(p.presence.flatten().tolist()) <= {0, 1}
    assert ((p.presence == 1) == (p.matrix > 0)).all()


def test_profiles_excludes_extinct_nodes():
    # extinct/internal nodes never appear as columns, only observed extant tips
    sp, g = _run(seed=3, death=0.7)
    extinct = {n.id for n in sp.complete_tree.extinct()}
    assert extinct and not (extinct & set(g.profiles.species))


def test_profiles_is_cached():
    _, g = _run(seed=6)
    assert g.profiles is g.profiles                    # derived once, reused


def test_empty_run_has_empty_profiles():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1)
    g = simulate_genomes_unordered(sp, seed=1)          # no families, no rates
    assert g.profiles.families == ()
    assert g.profiles.matrix.shape == (0, 10)           # no families, still the 10 extant columns


# --- write: materialise the chosen outputs ---------------------------------

def test_write_produces_events_and_profiles(tmp_path):
    _, g = _run(seed=7)
    g.write(tmp_path)
    # the default writes what the run computed: the log, the profiles, the genomes at every node,
    # and a gene tree per family
    written = {p.name for p in tmp_path.iterdir()}
    assert {"genome_events.tsv", "profiles.tsv", "genomes.tsv"} <= written
    assert any(n.startswith("gene_tree_fam") for n in written)
    ev = (tmp_path / "genome_events.tsv").read_text().splitlines()
    assert ev[0].split("\t") == ["time", "kind", "lineage", "family", "copy", "parent", "recipient",
                                 "donor"]
    assert len(ev) - 1 == len(g.events)                 # one row per event
    pr = (tmp_path / "profiles.tsv").read_text().splitlines()
    assert len(pr) - 1 == len(g.profiles.families)      # one row per family


def test_write_is_selective(tmp_path):
    _, g = _run(seed=8)
    g.write(tmp_path, outputs=("profiles",))
    assert [p.name for p in tmp_path.iterdir()] == ["profiles.tsv"]


def test_write_presence_tsv_is_binary():
    _, g = _run(seed=9)
    cells = [c for row in g.profiles.to_tsv(presence=True).splitlines()[1:] for c in row.split("\t")[1:]]
    assert set(cells) <= {"0", "1"}


def test_profiles_is_the_public_type():
    _, g = _run(seed=1)
    assert isinstance(g.profiles, Profiles)


# --- the initial genome ----------------------------------------------------------------------------

def test_the_initial_genome_is_the_one_the_run_started_with(tmp_path):
    """The genome at the START of the root branch. It is not genomes[root]: a node sits at the END of
    its branch, and the root branch is real simulated time, so events happen along it."""
    from zombi2.species import simulate_species_tree
    from zombi2.genomes import simulate_genomes_unordered

    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=5, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.5, loss=0.5, initial_families=6, seed=1)
    assert len(g.initial_genome) == 6                     # one copy per seeded family
    assert sorted(c.family for c in g.initial_genome) == list(range(6))
    root = sp.complete_tree.root
    assert g.initial_genome != g.genomes[root], "the stem was quiet — pick another seed"

    g.write(tmp_path)
    rows = (tmp_path / "initial_genome.tsv").read_text().splitlines()
    assert rows[0] == "family\tcopy" and len(rows) == 7    # its own file, no lineage column
    assert "lineage" not in rows[0]


def test_the_initial_genome_survives_a_run_that_loses_everything():
    # the run starts with what it starts with, whatever later becomes of it
    from zombi2.species import simulate_species_tree
    from zombi2.genomes import simulate_genomes_unordered

    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=3)
    g = simulate_genomes_unordered(sp, loss=6.0, initial_families=5, seed=3)
    assert len(g.initial_genome) == 5
    assert sum(len(g.genomes[n.id]) for n in sp.complete_tree.extant()) < 5 * 4


def test_a_transfer_names_both_ends_of_its_edge_on_every_row():
    """A transfer is an edge between two branches, and each of its two rows names both. Without
    `donor` the arriving row said only where the copy landed — twice over, since `lineage` and
    `recipient` are the same branch there — so reading who donated to whom meant pairing the rows on
    (time, parent). It is also what makes a self-transfer visible: donor == recipient."""
    from zombi2.species import simulate_species_tree
    from zombi2.genomes import simulate_genomes_unordered

    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=4)
    r = simulate_genomes_unordered(sp, initial_families=8, transfer=1.0, seed=4)
    transfers = [e for e in r.events if e.kind == "transfer"]
    assert transfers
    for e in transfers:
        assert e.donor is not None                       # on both rows, arriving and departing
        assert e.lineage in (e.donor, e.recipient)
        if e.recipient is None:                          # the donor's own continuation
            assert e.lineage == e.donor
        else:                                            # the copy that arrived
            assert e.lineage == e.recipient
    assert all(e.donor is None for e in r.events if e.kind != "transfer")


def test_self_transfers_are_readable_from_one_row():
    from zombi2.genomes._transfer import Distance
    from zombi2.species import simulate_species_tree
    from zombi2.genomes import simulate_genomes_unordered

    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=4)
    r = simulate_genomes_unordered(sp, initial_families=8, transfer=1.0, self_transfer=True,
                                   transfer_to=Distance(decay=10.0), seed=4)
    arrived = [e for e in r.events if e.kind == "transfer" and e.recipient is not None]
    selfies = [e for e in arrived if e.donor == e.recipient]
    assert selfies, "a steep distance decay with self_transfer should give plenty of self-transfers"
    assert len(selfies) < len(arrived), "and not all of them"
