"""``run.zombi2`` — the one-page run report (:mod:`zombi2._runtime.report`).

The report is a derived projection of the per-level records a run already writes, so the tests run real
pipelines and check what the report says against them: that every level appears, that every output file
is documented, that the reproduce block actually reproduces, and that a changed upstream is flagged.
"""
from __future__ import annotations

import os
import shlex

from zombi2._runtime.report import (RUN_REPORT_NAME, _GLOSS, _RECORD_SUFFIXES, build_run_report,
                                    write_run_report)
from zombi2.cli.main import main

_LEVELS = ("species", "genomes", "sequences", "traits")


def _pipeline(run, *, kind: str = "continuous") -> None:
    """A small four-level run into ``run`` (grouped layout), fast enough for a test."""
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1"])
    main(["genomes", str(run), "--initial-families", "5", "--duplication", "0.3", "--transfer", "0.1",
          "--loss", "0.2", "--seed", "2"])
    main(["sequences", str(run), "--model", "hky85", "--length", "30", "--seed", "1"])
    if kind == "discrete":
        main(["traits", str(run), "--kind", "discrete", "--states", "a,b", "--switch", "0.2", "--seed", "1"])
    else:
        main(["traits", str(run), "--kind", "continuous", "--rate", "1.0", "--seed", "1"])


def test_a_full_pipeline_writes_one_report_at_the_run_root(tmp_path):
    _pipeline(tmp_path)
    report = tmp_path / RUN_REPORT_NAME
    assert report.exists()
    text = report.read_text(encoding="utf-8")
    for section in ("SPECIES", "GENOMES", "SEQUENCES", "TRAITS", "TO REPRODUCE"):
        assert section in text, f"the report is missing its {section} section"
    # the software environment, for a reproduction or a bug report
    assert "built with zombi2" in text and "numpy" in text and "python" in text
    # the report belongs to the run, not to any one level's directory
    assert [p.name for p in tmp_path.iterdir() if p.is_file()] == [RUN_REPORT_NAME]


def test_every_output_file_is_documented(tmp_path):
    """The report must account for every file a run writes — the guarantee this feature was asked for.
    A new output with no gloss in ``_GLOSS`` (and that is not a record) fails here, so the report can
    never silently fall behind what the levels produce."""
    _pipeline(tmp_path, kind="discrete")            # discrete so trait_events.tsv is written too
    undocumented = []
    for level in _LEVELS:
        for entry in os.listdir(tmp_path / level):
            if entry == RUN_REPORT_NAME or entry.endswith(_RECORD_SUFFIXES):
                continue                            # the report itself, or a .log / _summary.json record
            if entry not in _GLOSS:
                undocumented.append(f"{level}/{entry}")
    assert not undocumented, ("outputs with no description — add a gloss in report._GLOSS: "
                              f"{undocumented}")


def test_per_family_directories_are_summarised_not_enumerated(tmp_path):
    _pipeline(tmp_path)
    text = (tmp_path / RUN_REPORT_NAME).read_text(encoding="utf-8")
    assert "gene_trees/" in text and "files)" in text          # the directory, with a count
    assert "gene_tree_fam0" not in text                        # not one line per family


def test_the_reproduce_block_actually_reproduces(tmp_path):
    """The strongest check: run each command the report prints, into a fresh directory, and confirm it
    succeeds and lands the same tree — the reproduce block is a promise, so it is executed here."""
    a, b = tmp_path / "a", tmp_path / "b"
    _pipeline(a)
    ran = 0
    for raw in (a / RUN_REPORT_NAME).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line.startswith("zombi2 "):
            continue
        argv = [str(b) if tok == str(a) else tok for tok in shlex.split(line)[1:]]
        assert main(argv) == 0, f"a reproduce command failed: {line}"
        ran += 1
    assert ran == 4, "expected one reproduce command per level"
    for level in _LEVELS:                                       # deterministic: same seed, same output
        for name in os.listdir(a / level):
            if name.endswith(".nwk"):
                assert (a / level / name).read_bytes() == (b / level / name).read_bytes()


def test_a_changed_upstream_is_flagged_as_stale(tmp_path):
    _pipeline(tmp_path)
    tree = tmp_path / "species" / "species_complete.nwk"
    tree.write_text(tree.read_text(encoding="utf-8") + "\n", encoding="utf-8")   # change its bytes
    write_run_report(str(tmp_path))
    text = (tmp_path / RUN_REPORT_NAME).read_text(encoding="utf-8")
    assert "changed since" in text and "GENOMES was computed on" in text


def test_no_report_for_a_flat_or_empty_run(tmp_path):
    assert build_run_report(str(tmp_path)) is None      # nothing there yet
    flat = tmp_path / "flat"
    main(["species", str(flat), "--birth", "1", "--death", "0.3", "--n-extant", "6", "--seed", "1", "--flat"])
    assert build_run_report(str(flat)) is None          # --flat keeps no per-level records to report on
    assert not (flat / RUN_REPORT_NAME).exists()


def test_ordered_and_nucleotide_outputs_are_documented(tmp_path):
    """The other two genome resolutions write different files (gene_order, chromosome_events; blocks,
    genes, bed/, gff/) and nucleotide writes no summary. Every file must still be documented, and the
    section must still render — from the log's result line when there is no summary."""
    for res, extra in (("ordered", ["--initial-families", "8", "--inversion", "0.2", "--chromosomes", "2"]),
                       ("nucleotide", ["--root-length", "400"])):
        run = tmp_path / res
        main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "6", "--seed", "1"])
        main(["genomes", str(run), "--resolution", res, "--duplication", "0.1", "--loss", "0.1",
              "--origination", "0.3", "--seed", "3", *extra])
        text = (run / RUN_REPORT_NAME).read_text(encoding="utf-8")
        assert "GENOMES" in text, f"{res}: no genomes section in the report"
        undocumented = [f"genomes/{e}" for e in os.listdir(run / "genomes")
                        if e != RUN_REPORT_NAME and not e.endswith(_RECORD_SUFFIXES) and e not in _GLOSS]
        assert not undocumented, f"{res}: undocumented outputs {undocumented}"


def test_staleness_survives_a_move_and_a_different_cwd(tmp_path, monkeypatch):
    """The staleness warning must hold when the report is regenerated from another working directory, or
    after the run directory is moved — the inputs are located relative to the run dir, not the CWD."""
    monkeypatch.chdir(tmp_path)                              # record with a RELATIVE run-arg
    main(["species", "run", "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1"])
    main(["genomes", "run", "--initial-families", "4", "--duplication", "0.2", "--seed", "1"])
    tree = tmp_path / "run" / "species" / "species_complete.nwk"
    tree.write_text(tree.read_text(encoding="utf-8") + "\n", encoding="utf-8")   # change the upstream
    moved = tmp_path / "moved"
    (tmp_path / "run").rename(moved)                         # move the whole run
    monkeypatch.chdir(tmp_path.parent)                       # and rebuild from a different cwd
    assert "GENOMES was computed on" in build_run_report(str(moved)), \
        "a changed upstream must be flagged even after a move / from another cwd"


def test_a_joint_run_reports_species_and_its_driver(tmp_path):
    main(["joint", str(tmp_path), "--birth", "1.0 * DrivenBy('trait', {'a': 1.0, 'b': 2.0})",
          "--death", "0.2", "--states", "a,b", "--switch", "0.3", "--n-extant", "20", "--seed", "1"])
    text = (tmp_path / RUN_REPORT_NAME).read_text(encoding="utf-8")
    assert "JOINT" in text and "driven by a trait" in text
    assert "traits/trait_values.tsv" in text            # the driver's files listed with the species'
    assert "zombi2 joint" in text                       # a reproduce line for the joint command
