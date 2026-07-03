"""End-to-end tests for the command-line wrapper (``python -m zombi2 ...``)."""

from __future__ import annotations

import os

import pytest

from zombi2 import rust_available
from zombi2.cli import main

# The `genomes` / `all` subcommands run the built-in model, which uses the Rust engine.
needs_rust = pytest.mark.skipif(not rust_available(),
                                reason="zombi2_core (Rust extension) not built")


def test_species_writes_newick(tmp_path):
    out = tmp_path / "sp"
    rc = main(["species", "--birth", "1", "--death", "0.3",
               "--tips", "15", "--age", "5", "--seed", "1", "-o", str(out)])
    assert rc == 0
    newick = (out / "species_tree.nwk").read_text().strip()
    assert newick.endswith(";")


@needs_rust
def test_genomes_on_supplied_tree(tmp_path):
    """`species` output feeds straight into `genomes` — the round-trip the CLI enables."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "20", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--max-family-size", "0.5", "--seed", "42", "-o", str(gen)])
    assert rc == 0
    for f in ("species_tree.nwk", "Profiles.tsv", "Transfers.tsv"):
        assert (gen / f).exists()
    assert os.listdir(gen / "gene_trees")


@needs_rust
def test_genomes_writes_full_output(tmp_path):
    """The default `genomes` run writes the full ZOMBI-1 output (Rust engine)."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "-o", str(gen)])
    assert rc == 0
    for f in ("species_tree.nwk", "Profiles.tsv", "Presence.tsv", "Transfers.tsv",
              "Gene_family_summary.tsv"):
        assert (gen / f).exists()
    assert os.listdir(gen / "gene_trees")
    assert os.listdir(gen / "gene_family_events")


@needs_rust
def test_genomes_profiles_only(tmp_path):
    """`--profiles-only` writes just the profile matrices — no gene trees / event log."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "--profiles-only", "-o", str(gen)])
    assert rc == 0
    assert (gen / "Profiles.tsv").exists()
    assert (gen / "Presence.tsv").exists()
    assert not (gen / "gene_trees").exists()
    assert not (gen / "Transfers.tsv").exists()


@needs_rust
def test_all_runs_end_to_end(tmp_path):
    out = tmp_path / "all"
    rc = main(["all", "--birth", "1", "--death", "0.2", "--tips", "15", "--age", "5",
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "-o", str(out)])
    assert rc == 0
    assert (out / "Profiles.tsv").exists()


def test_max_family_size_parses_int_and_fraction(tmp_path):
    """Integer -> absolute cap; a decimal -> fraction of the number of species."""
    from zombi2.cli import _int_or_float

    assert _int_or_float("40") == 40 and isinstance(_int_or_float("40"), int)
    assert _int_or_float("0.5") == 0.5 and isinstance(_int_or_float("0.5"), float)


@needs_rust
def test_seed_makes_genomes_reproducible(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "20", "--age", "5", "--seed", "1", "-o", str(sp)])
    tree = str(sp / "species_tree.nwk")

    def run(name):
        dest = tmp_path / name
        main(["genomes", "--tree", tree, "--dup", "0.2", "--loss", "0.2",
              "--orig", "0.5", "--seed", "7", "-o", str(dest)])
        return (dest / "Profiles.tsv").read_text()

    assert run("a") == run("b")
