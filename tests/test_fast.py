"""Optional Rust fast-path engine (profiles-only). Skipped if zombi2_core isn't built."""

import numpy as np
import pytest

import zombi2 as z

pytestmark = pytest.mark.skipif(not z.rust_available(),
                                reason="zombi2_core (Rust extension) not built")

RATES = dict(duplication=0.15, transfer=0.1, loss=0.2, origination=0.5,
             initial_size=20, max_family_size=0.5)


def _tree(n=30, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n, age=5.0, seed=seed)


def test_reproducible_same_seed():
    tree = _tree()
    a = z.simulate_profiles_fast(tree, seed=7, **RATES)
    b = z.simulate_profiles_fast(tree, seed=7, **RATES)
    assert np.array_equal(a.matrix, b.matrix)


def test_shape_and_species_are_extant_leaves():
    tree = _tree(n=30)
    pm = z.simulate_profiles_fast(tree, seed=3, **RATES)
    assert pm.matrix.shape[1] == len(tree.extant_leaves())
    assert set(pm.species) == {n.name for n in tree.extant_leaves()}
    assert (pm.matrix >= 0).all()


def test_hard_cap_respected():
    tree = _tree(n=40)
    cap = 5
    pm = z.simulate_profiles_fast(tree, duplication=0.6, transfer=0.2, loss=0.05,
                                  origination=0.4, initial_size=20, max_family_size=cap,
                                  seed=11)
    assert pm.matrix.max() <= cap


def test_accepts_uniform_rates_object():
    tree = _tree()
    obj = z.simulate_profiles_fast(tree, z.UniformRates(0.15, 0.1, 0.2, 0.5),
                                   initial_size=20, max_family_size=0.5, seed=7)
    kw = z.simulate_profiles_fast(tree, seed=7, **RATES)
    assert np.array_equal(obj.matrix, kw.matrix)


def test_rejects_unsupported_models():
    tree = _tree()
    with pytest.raises(ValueError):
        z.simulate_profiles_fast(tree, z.UniformRates(0.2, 0, 0.1, 0.3, carrying_capacity=10),
                                 seed=1)
    with pytest.raises(TypeError):
        z.simulate_profiles_fast(tree, z.GenomeWiseRates(0.2, 0.1, 0.2, 0.5), seed=1)


FULL = dict(duplication=0.15, transfer=0.1, loss=0.2, origination=0.5,
            initial_size=20, max_family_size=0.5)


def _extant_leaves(newick):
    # gene trees are strictly binary, so #leaves = #commas + 1
    return 0 if newick is None else newick.count(",") + 1


def test_full_log_reproducible():
    tree = _tree()
    a = z.simulate_genomes_fast(tree, seed=7, **FULL)
    b = z.simulate_genomes_fast(tree, seed=7, **FULL)
    assert len(a.event_log) == len(b.event_log)
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)


def test_full_log_gene_tree_invariant():
    # reconciliation invariant: extant gene-tree leaves == family's extant copy count
    tree = _tree(n=40, seed=3)
    g = z.simulate_genomes_fast(tree, seed=5, **FULL)
    fam_row = {f: i for i, f in enumerate(g.profiles.families)}
    checked = 0
    for fam, (_complete, extant) in g.gene_trees().items():
        copies = int(g.profiles.matrix[fam_row[fam]].sum()) if fam in fam_row else 0
        assert _extant_leaves(extant) == copies, (fam, _extant_leaves(extant), copies)
        checked += 1
    assert checked > 0


def test_full_log_writes_all_outputs(tmp_path):
    tree = _tree()
    g = z.simulate_genomes_fast(tree, seed=1, **FULL)
    g.write(tmp_path)
    for f in ["species_tree.nwk", "species_nodes.tsv", "Transfers.tsv",
              "Gene_family_summary.tsv", "Profiles.tsv", "Presence.tsv"]:
        assert (tmp_path / f).exists()
    assert (tmp_path / "gene_trees").is_dir()
    assert (tmp_path / "gene_family_events").is_dir()


def test_full_log_event_types_present():
    tree = _tree(n=40, seed=2)
    g = z.simulate_genomes_fast(tree, seed=9, **FULL)
    kinds = {r.event.value for r in g.event_log}
    assert {"O", "S", "D", "L"}.issubset(kinds)  # transfers may be absent by chance


def test_full_log_rejects_unsupported():
    tree = _tree()
    with pytest.raises(TypeError):
        z.simulate_genomes_fast(tree, z.GenomeWiseRates(0.2, 0.1, 0.2, 0.5), seed=1)


def test_full_log_family_counts_match_profiles_fast():
    # both fast paths run the same model; mean family count should agree
    tree = _tree(n=60, seed=2)
    a = np.mean([len(z.simulate_genomes_fast(tree, duplication=0.2, loss=0.1, origination=0.4,
                                             initial_size=30, max_family_size=0.5,
                                             seed=100 + s).profiles.families) for s in range(8)])
    b = np.mean([z.simulate_profiles_fast(tree, duplication=0.2, loss=0.1, origination=0.4,
                                          initial_size=30, max_family_size=0.5,
                                          seed=200 + s).matrix.shape[0] for s in range(8)])
    assert abs(a - b) / b < 0.2


def _read_profiles(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    return rows[0][1:], {r[0]: list(map(int, r[1:])) for r in rows[1:]}


def test_write_fast_produces_all_files(tmp_path):
    tree = _tree()
    s = z.simulate_and_write_fast(tree, tmp_path, seed=7, **FULL)
    assert s["n_species"] == len(tree.extant_leaves())
    for f in ["species_tree.nwk", "species_nodes.tsv", "Transfers.tsv",
              "Gene_family_summary.tsv", "Profiles.tsv", "Presence.tsv"]:
        assert (tmp_path / f).exists()
    assert (tmp_path / "gene_trees").is_dir()
    assert (tmp_path / "gene_family_events").is_dir()


def test_write_fast_gene_tree_invariant(tmp_path):
    tree = _tree(n=40, seed=3)
    z.simulate_and_write_fast(tree, tmp_path, seed=5, **FULL)
    _, prof = _read_profiles(tmp_path / "Profiles.tsv")
    checked = 0
    for fam, row in prof.items():
        rowsum = sum(row)
        ep = tmp_path / "gene_trees" / f"{fam}_extant.nwk"
        n_leaves = (ep.read_text().count(",") + 1) if ep.exists() else 0
        assert n_leaves == rowsum
        checked += 1
    assert checked > 0


def test_write_fast_presence_matches_profiles(tmp_path):
    z.simulate_and_write_fast(_tree(), tmp_path, seed=1, **FULL)
    _, prof = _read_profiles(tmp_path / "Profiles.tsv")
    _, pres = _read_profiles(tmp_path / "Presence.tsv")
    for fam in prof:
        assert pres[fam] == [1 if x > 0 else 0 for x in prof[fam]]


def test_write_fast_reproducible(tmp_path):
    tree = _tree()
    z.simulate_and_write_fast(tree, tmp_path / "a", seed=9, **FULL)
    z.simulate_and_write_fast(tree, tmp_path / "b", seed=9, **FULL)
    assert (tmp_path / "a" / "Profiles.tsv").read_text() == (tmp_path / "b" / "Profiles.tsv").read_text()
    assert (tmp_path / "a" / "Transfers.tsv").read_text() == (tmp_path / "b" / "Transfers.tsv").read_text()


def test_write_fast_rejects_unsupported(tmp_path):
    with pytest.raises(TypeError):
        z.simulate_and_write_fast(_tree(), tmp_path, z.GenomeWiseRates(0.2, 0.1, 0.2, 0.5), seed=1)


def test_statistically_matches_python_engine():
    # mean copy-number over the matrix should agree within Monte-Carlo error
    tree = _tree(n=60, seed=2)
    r = np.mean([z.simulate_profiles_fast(tree, duplication=0.2, loss=0.1, origination=0.4,
                                          initial_size=30, max_family_size=0.5,
                                          seed=1000 + s).matrix.mean() for s in range(15)])
    p = np.mean([z.simulate_genomes(tree, duplication=0.2, loss=0.1, origination=0.4,
                                    initial_size=30, max_family_size=0.5,
                                    seed=2000 + s).profiles.matrix.mean() for s in range(15)])
    assert abs(r - p) / p < 0.15
