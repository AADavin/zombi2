"""Replicate-level parallel runner."""

from zombi2 import BirthDeath, SharedRates, run_replicates


def _summary_keys(rows):
    return [(r["replicate"], r["seed"], r["n_species"], r["n_families"], r["n_events"])
            for r in rows]


def test_run_replicates_serial_writes_outputs(tmp_path):
    rows = run_replicates(
        4, str(tmp_path), BirthDeath(1.0, 0.3), n_tips=8, age=4.0,
        duplication=0.15, transfer=0.1, loss=0.2, origination=0.5,
        initial_families=8, seed=1, processes=1,
    )
    assert len(rows) == 4
    for i, r in enumerate(rows):
        assert r["replicate"] == i
        assert r["n_species"] == 8
        d = tmp_path / f"replicate_{i:04d}"
        assert (d / "species_tree.nwk").exists()
        assert (d / "Profiles.tsv").exists()
        assert (d / "gene_trees").is_dir()


def test_replicates_are_independent():
    # different seeds -> the replicates are not all identical
    rows = run_replicates(
        5, "/tmp/zombi2_indep", BirthDeath(1.0, 0.3), n_tips=10, age=4.0,
        duplication=0.2, transfer=0.1, loss=0.2, origination=0.5,
        initial_families=10, seed=7, processes=1,
    )
    assert len({r["n_events"] for r in rows}) > 1


def test_reproducible_by_base_seed(tmp_path):
    kw = dict(n_tips=8, age=4.0, duplication=0.15, transfer=0.1, loss=0.2,
              origination=0.5, initial_families=8, seed=42, processes=1)
    a = run_replicates(3, str(tmp_path / "a"), BirthDeath(1.0, 0.3), **kw)
    b = run_replicates(3, str(tmp_path / "b"), BirthDeath(1.0, 0.3), **kw)
    assert _summary_keys(a) == _summary_keys(b)


def test_parallel_matches_serial(tmp_path):
    # results depend only on the (base seed, replicate index), not the process count
    kw = dict(species_model=BirthDeath(1.0, 0.3), n_tips=10, age=4.0,
              rates=SharedRates(0.2, 0.1, 0.2, 0.5), initial_families=10, seed=5)
    serial = run_replicates(4, str(tmp_path / "s"), **kw, processes=1)
    parallel = run_replicates(4, str(tmp_path / "p"), **kw, processes=2)
    assert _summary_keys(serial) == _summary_keys(parallel)
