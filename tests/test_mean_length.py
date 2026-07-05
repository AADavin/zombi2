"""The structural-event length knob is a mean length (--mean-length), converted to the engine's
geometric continuation probability."""

import math

import pytest

from zombi2.cli import _extension_from_mean_length, main
from zombi2.species import BirthDeath, simulate_species_tree


def _tree_file(tmp_path):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=1.0, seed=1)
    p = tmp_path / "sp.nwk"
    p.write_text(tree.to_newick() + "\n")
    return str(p)


def test_mean_length_to_extension_conversion():
    assert _extension_from_mean_length(None) is None            # per-level default preserved
    assert _extension_from_mean_length(1.0) == 0.0              # single-element events
    assert _extension_from_mean_length(2.0) == 0.5
    assert math.isclose(_extension_from_mean_length(100.0), 0.99)
    with pytest.raises(ValueError):
        _extension_from_mean_length(0.5)                        # a segment spans at least one unit


def test_cli_ordered_accepts_mean_length(tmp_path):
    tree = _tree_file(tmp_path)
    out = tmp_path / "ord"
    rc = main(["genomes", "-t", tree, "--genome-model", "ordered", "--rate-model", "shared",
               "--inversion", "0.3", "--mean-length", "3", "--write", "profiles",
               "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "Profiles.tsv").exists()


def test_cli_nucleotide_accepts_mean_length(tmp_path):
    tree = _tree_file(tmp_path)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--genome-model", "nucleotide", "--inversion", "0.004",
               "--root-length", "800", "--mean-length", "50", "--write", "profiles",
               "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "Profiles.tsv").exists()


def test_cli_mean_length_below_one_is_a_clean_error(tmp_path, capsys):
    tree = _tree_file(tmp_path)
    rc = main(["genomes", "-t", tree, "--genome-model", "ordered", "--rate-model", "shared",
               "--inversion", "0.3", "--mean-length", "0.5", "-o", str(tmp_path / "x")])
    assert rc == 1
    assert "mean-length" in capsys.readouterr().err
