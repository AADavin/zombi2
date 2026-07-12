"""Regression tests for the CLI/docs contract fixes (2026-07 audit).

* ``zombi2 genomes --write reconciliations`` writes ``Reconciled_complete/extant.nwk`` (tips
  ``<species>|<gid>``) — the truth input ``tools recon-accuracy`` and ``tools reconcile`` read.
  Before the fix ``reconciliations`` was not a valid ``--write`` value, so the documented
  benchmark workflow could not be started from the CLI.
* ``--threads > 1`` is rejected up front with a clear, flag-level message for every
  non-parallelisable configuration, instead of crashing deep in the engine (per-genome/family
  rates) or silently running serial (``--conversion`` / ``--genome-model ordered``).
"""

import pytest

from zombi2 import rust_available
from zombi2.cli import main

needs_rust = pytest.mark.skipif(not rust_available(),
                                reason="zombi2_core (Rust extension) not built")

SPECIES_TREE = "((a:1,b:1)i:1,(c:1,d:1)j:1)R:0;"


def _species_tree(tmp_path):
    p = tmp_path / "sp.nwk"
    p.write_text(SPECIES_TREE + "\n")
    return str(p)


@needs_rust
def test_write_reconciliations_produces_scorable_truth(tmp_path, capsys):
    sp = _species_tree(tmp_path)
    out = tmp_path / "rc"
    rc = main(["genomes", "-t", sp, "--dup", "0.3", "--trans", "0.1", "--loss", "0.2",
               "--initial-families", "6", "--seed", "7", "--write", "reconciliations",
               "-o", str(out)])
    assert rc == 0
    for name in ("Reconciled_complete.nwk", "Reconciled_extant.nwk", "Reconciliation_events.tsv"):
        assert (out / name).exists(), name
    extant_file = out / "Reconciled_extant.nwk"
    lines = extant_file.read_text().strip().splitlines()
    assert lines, "no reconciled gene trees written"
    assert lines[0].endswith(";") and "|" in lines[0]      # tips are <species>|<gid>

    # the file feeds `tools reconcile` (the documented workflow that previously failed)
    capsys.readouterr()
    rc = main(["tools", "reconcile", "-g", str(extant_file), "-t", sp,
               "--dup", "0.3", "--trans", "0.1", "--loss", "0.2"])
    assert rc == 0

    # scoring the truth against itself is a perfect reconciliation
    capsys.readouterr()
    rc = main(["tools", "recon-accuracy", "-t", str(extant_file), "-i", str(extant_file)])
    assert rc == 0
    assert "joint_acc=1.0000" in capsys.readouterr().out


@needs_rust
def test_write_reconciliations_nucleotide(tmp_path):
    sp = _species_tree(tmp_path)
    out = tmp_path / "rcn"
    rc = main(["genomes", "-t", sp, "--genome-model", "nucleotide",
               "--dup", "0.2", "--loss", "0.2", "--seed", "3",
               "--write", "reconciliations", "-o", str(out)])
    assert rc == 0
    assert (out / "Reconciled_extant.nwk").exists()


@pytest.mark.parametrize("extra, needle", [
    (["--rate-model", "per-genome"], "rate-model shared"),
    (["--conversion", "0.1"], "--conversion"),
    (["--genome-model", "ordered"], "unordered"),
])
def test_threads_gt1_rejected_with_flag_level_message(tmp_path, capsys, extra, needle):
    sp = _species_tree(tmp_path)
    with pytest.raises(SystemExit):
        main(["genomes", "-t", sp, "--dup", "0.2", "--loss", "0.2",
              "--write", "profiles", "--threads", "4", "-o", str(tmp_path / "x")] + extra)
    err = capsys.readouterr().err
    assert "counts-only path" in err and needle in err


def test_threads_gt1_rejected_for_non_profiles_write(tmp_path, capsys):
    sp = _species_tree(tmp_path)
    with pytest.raises(SystemExit):
        main(["genomes", "-t", sp, "--dup", "0.2", "--loss", "0.2",
              "--write", "trees", "--threads", "4", "-o", str(tmp_path / "x")])
    assert "exactly --write profiles" in capsys.readouterr().err


@needs_rust
def test_threads_gt1_valid_combo_still_runs(tmp_path):
    sp = _species_tree(tmp_path)
    rc = main(["genomes", "-t", sp, "--dup", "0.2", "--loss", "0.2",
               "--write", "profiles", "--threads", "4", "-o", str(tmp_path / "ok")])
    assert rc == 0
    assert (tmp_path / "ok" / "Profiles.tsv").exists()
