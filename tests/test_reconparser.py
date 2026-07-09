"""Tests for the vendored ``reconparser`` interop tool (zombi2.tools.reconparser).

The parsers need the optional ``reconparser`` extra (ete3, pandas); every test is skipped if
those are absent. We exercise the ALE parser against a tiny hand-built ALE output fixture (a
full AleRax run directory is too large to fixture here), plus the ``zombi2 tools parse`` CLI.
"""

import textwrap

import pytest

pytest.importorskip("ete3")
pytest.importorskip("pandas")

from zombi2.cli import main  # noqa: E402
from zombi2.tools.reconparser import ALEParser, AleRaxFamily, AleRaxRun  # noqa: E402


# A minimal but structurally-complete ALE output triple (v1.0-style .uml_rec).
UCONS_TREE = "#constructor string\n((A:1,B:1):1,C:1);\n"

UTS = "#from\tto\tfreq\nA\tB\t0.5\nB\tC\t0.3\n"

UML_REC = textwrap.dedent("""\
    #ALErec v1.0
    S:\t((A:1,B:1):1,C:1):0;

    Input ale from:\ttest.ale
    >logl: -42.500000
    rate of\tDuplications\tTransfers\tLosses
    ML \t0.1\t0.05\t0.15
    2\treconciled G-s:
    (A:1,B:1);
    ((A:1,B:1):1,C:1);

    # of\tDuplications\tTransfers\tLosses\tSpeciations
    Total \t3\t2\t5\t4
    # of\tDuplications\tTransfers\tLosses\tOriginations\tcopies
    S_terminal_branch\tA(0)\t1\t0\t1\t1\t1
    S_terminal_branch\tB(1)\t0\t1\t2\t1\t1
    S_internal_branch\t3\t2\t1\t2\t0\t1
    """)


@pytest.fixture
def ale_base(tmp_path):
    """Write a tiny ALE output triple and return its base path (``.../test.ale``)."""
    base = tmp_path / "test.ale"
    (tmp_path / "test.ale.ucons_tree").write_text(UCONS_TREE)
    (tmp_path / "test.ale.uTs").write_text(UTS)
    (tmp_path / "test.ale.uml_rec").write_text(UML_REC)
    return base


# --------------------------------------------------------------------------- #
# import surface
# --------------------------------------------------------------------------- #
def test_public_surface():
    # the three documented classes are what the package exposes, and nothing here leaks
    # into the top-level zombi2 namespace
    import zombi2
    import zombi2.tools.reconparser as rp

    assert set(rp.__all__) == {"ALEParser", "AleRaxRun", "AleRaxFamily"}
    for cls in (ALEParser, AleRaxRun, AleRaxFamily):
        assert isinstance(cls, type)
    assert not hasattr(zombi2, "ALEParser")


# --------------------------------------------------------------------------- #
# ALE parser
# --------------------------------------------------------------------------- #
def test_ale_accepts_base_or_file_path(tmp_path, ale_base):
    # a specific-file path is normalised back to the base path
    from_file = ALEParser(str(ale_base) + ".uml_rec")
    assert from_file.base_path == ale_base
    assert ALEParser(ale_base).base_path == ale_base


def test_ale_files_exist(ale_base):
    assert ALEParser(ale_base).files_exist() == {
        "consensus_tree": True,
        "transfers": True,
        "reconciliation": True,
    }


def test_ale_ml_rates_and_loglik(ale_base):
    p = ALEParser(ale_base)
    assert p.get_log_likelihood() == pytest.approx(-42.5)
    r = p.get_ml_rates()
    assert r == pytest.approx({"duplications": 0.1, "transfers": 0.05, "losses": 0.15})


def test_ale_transfers(ale_base):
    tr = ALEParser(ale_base).get_transfers()
    assert list(tr.columns) == ["from", "to", "freq"]
    assert len(tr) == 2
    top = tr.nlargest(1, "freq").iloc[0]
    assert (top["from"], top["to"]) == ("A", "B")
    assert top["freq"] == pytest.approx(0.5)


def test_ale_summary_and_branch_stats(ale_base):
    p = ALEParser(ale_base)
    s = p.get_summary_statistics()
    assert s["total_transfers"] == pytest.approx(2)
    bs = p.get_branch_statistics()
    # terminal branch id "A(0)" is cleaned to "A"
    assert set(bs["branch_id"]) == {"A", "B", "3"}
    assert bs.loc[bs["branch_id"] == "B", "transfers"].iloc[0] == pytest.approx(1)


def test_ale_consensus_tree(ale_base):
    tree = ALEParser(ale_base).get_consensus_tree()
    assert sorted(leaf.name for leaf in tree.get_leaves()) == ["A", "B", "C"]


def test_ale_missing_files_raise(tmp_path):
    p = ALEParser(tmp_path / "nope.ale")
    assert p.files_exist() == {
        "consensus_tree": False, "transfers": False, "reconciliation": False}
    with pytest.raises(FileNotFoundError):
        p.get_log_likelihood()


# --------------------------------------------------------------------------- #
# AleRax parser (surface only — a full run directory is out of scope to fixture)
# --------------------------------------------------------------------------- #
def test_alerax_run_requires_directory(tmp_path):
    with pytest.raises(NotADirectoryError):
        AleRaxRun(tmp_path / "does_not_exist")


# --------------------------------------------------------------------------- #
# CLI: zombi2 tools parse
# --------------------------------------------------------------------------- #
def test_cli_parse_ale(ale_base, capsys):
    rc = main(["tools", "parse", str(ale_base)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "ALE reconciliation" in out
    assert "-42.500000" in out
    assert "D=0.1" in out and "T=0.05" in out and "L=0.15" in out
    assert "A -> B" in out


def test_cli_parse_ale_explicit_file_and_top(ale_base, capsys):
    # a specific file path works too, and --top bounds the printed transfers
    rc = main(["tools", "parse", str(ale_base) + ".uml_rec", "--top", "1"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "A -> B" in out
    assert "B -> C" not in out


def test_cli_parse_ale_writes_tsv(ale_base, tmp_path, capsys):
    outdir = tmp_path / "out"
    rc = main(["tools", "parse", str(ale_base), "-o", str(outdir)])
    assert rc == 0
    assert (outdir / "ale_transfers.tsv").exists()
    assert (outdir / "ale_branch_statistics.tsv").exists()


def test_cli_parse_no_ale_files_errors(tmp_path, capsys):
    rc = main(["tools", "parse", str(tmp_path / "missing.ale"), "--tool", "ale"])
    assert rc == 1
    assert "no ALE output files" in capsys.readouterr().err


def test_cli_parse_alerax_bad_dir_errors(tmp_path, capsys):
    rc = main(["tools", "parse", str(tmp_path / "nodir"), "--tool", "alerax"])
    assert rc == 1
    assert "AleRax output directory not found" in capsys.readouterr().err
