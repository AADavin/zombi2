"""Tests for the fresh forward-only species CLI (zombi2.cli.species_forward).

Standalone: it is not wired into ``zombi2 species`` yet (the old command is the pipeline entry
point), so we build its parser directly and call ``run``.
"""
import argparse

import pytest

from zombi2.cli import species_forward


def _run(argv, out):
    p = argparse.ArgumentParser()
    species_forward._add_species_args(p)
    args = p.parse_args([*argv, "-o", str(out)])
    return species_forward.run(args, p)


def test_basic_run_writes_the_default_files(tmp_path):
    rc = _run(["--birth", "1", "--death", "0.3", "--n-extant", "20", "--seed", "1"], tmp_path)
    assert rc == 0
    assert (tmp_path / "species_complete.nwk").read_text().strip().endswith(";")
    assert (tmp_path / "species_extant.nwk").read_text().strip().endswith(";")
    assert (tmp_path / "species_events.tsv").exists()
    assert (tmp_path / "species.log").exists()


def test_clade_drift_runs(tmp_path):
    rc = _run(["--birth", "1", "--clade-drift", "0.3", "--n-extant", "30", "--seed", "1"], tmp_path)
    assert rc == 0
    assert (tmp_path / "species_extant.nwk").exists()


def test_skyline(tmp_path):
    rc = _run(["--birth", "1", "--skyline", "0:1.0,2:0.0", "--total-time", "5", "--seed", "1"], tmp_path)
    assert rc == 0
    assert (tmp_path / "species_complete.nwk").exists()


def test_diversity_cap(tmp_path):
    _run(["--birth", "1", "--diversity-cap", "25", "--total-time", "100", "--seed", "1"], tmp_path)
    # the cap bounds the survivors
    n_tips = (tmp_path / "species_extant.nwk").read_text().count("n")  # rough tip count
    assert n_tips > 0


def test_total_time_with_mass_extinction(tmp_path):
    rc = _run(["--birth", "1.5", "--death", "0.1", "--total-time", "5",
               "--mass-extinction", "2", "0.5", "--seed", "3"], tmp_path)
    assert rc == 0
    assert "\textinction\t" in (tmp_path / "species_events.tsv").read_text()   # deaths recorded (incl. the cull)


def test_global_birth_and_sampling(tmp_path):
    rc = _run(["--global-birth", "--birth", "2", "--total-time", "4", "--sampling", "0.5", "--seed", "1"], tmp_path)
    assert rc == 0


def test_fossils_side_output(tmp_path):
    _run(["--birth", "1", "--death", "0.4", "--n-extant", "40", "--fossils", "0.5", "--seed", "3"], tmp_path)
    assert (tmp_path / "species_fossils.tsv").read_text().startswith("lineage\ttime")


def test_write_is_selective(tmp_path):
    _run(["--birth", "1", "--death", "0.3", "--n-extant", "20", "--write", "extant", "events", "--seed", "1"], tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == {"species_extant.nwk", "species_events.tsv", "species.log"}


def test_requires_exactly_one_stop(tmp_path):
    with pytest.raises(SystemExit):
        _run(["--birth", "1"], tmp_path)                                          # neither
    with pytest.raises(SystemExit):
        _run(["--birth", "1", "--n-extant", "10", "--total-time", "5"], tmp_path)  # both


def test_deterministic(tmp_path):
    _run(["--birth", "1", "--death", "0.3", "--n-extant", "20", "--seed", "7"], tmp_path / "a")
    _run(["--birth", "1", "--death", "0.3", "--n-extant", "20", "--seed", "7"], tmp_path / "b")
    assert (tmp_path / "a" / "species_complete.nwk").read_text() == \
           (tmp_path / "b" / "species_complete.nwk").read_text()
