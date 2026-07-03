"""End-to-end tests for the command-line wrapper (``python -m zombi2 ...``)."""

from __future__ import annotations

import os

import pytest

from zombi2 import rust_available
from zombi2.cli import main

# The `genomes` subcommand runs the built-in model, which uses the Rust engine.
needs_rust = pytest.mark.skipif(not rust_available(),
                                reason="zombi2_core (Rust extension) not built")


def test_species_writes_newick(tmp_path):
    out = tmp_path / "sp"
    rc = main(["species", "--birth", "1", "--death", "0.3",
               "--tips", "15", "--age", "5", "--seed", "1", "-o", str(out)])
    assert rc == 0
    newick = (out / "species_tree.nwk").read_text().strip()
    assert newick.endswith(";")


def test_species_bare_defaults(tmp_path):
    """`zombi2 species -o out` runs with defaults — no required rate/size args."""
    out = tmp_path / "sp"
    rc = main(["species", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")


def test_species_forward_keeps_extinct(tmp_path):
    """Forward mode grows the complete tree, keeping extinct lineages (not in backward)."""
    out = tmp_path / "fwd"
    rc = main(["species", "--model", "forward", "--age", "6", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()


def test_species_forward_requires_exactly_one_of_tips_age(tmp_path):
    """Forward needs exactly one of --tips/--age; neither or both is an error."""
    with pytest.raises(SystemExit):
        main(["species", "--model", "forward", "-o", str(tmp_path / "a")])
    with pytest.raises(SystemExit):
        main(["species", "--model", "forward", "--tips", "20", "--age", "5",
              "-o", str(tmp_path / "b")])


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


# --- trait command: overlay a trait on a GIVEN tree (pure-Python, no Rust) --
def _tree_file(tmp_path, tips=10, seed=1):
    sp = tmp_path / "sp"
    main(["species", "--tips", str(tips), "--seed", str(seed), "-o", str(sp)])
    return str(sp / "species_tree.nwk")


def test_trait_writes_tips_and_ancestral_values(tmp_path):
    tree = _tree_file(tmp_path, tips=10)
    dest = tmp_path / "tr"
    rc = main(["trait", "-t", tree, "--model", "ou", "--alpha", "3", "--theta", "5",
               "--seed", "2", "-o", str(dest)])
    assert rc == 0
    rows = (dest / "traits.tsv").read_text().splitlines()
    assert rows[0] == "node\ttrait"
    names = {r.split("\t")[0] for r in rows[1:]}
    assert any(n.startswith("n") for n in names)             # tips present
    assert "root" in names and any(n.startswith("i") for n in names)   # ancestral values present
    newick = (dest / "trait_tree.nwk").read_text()
    assert newick.strip().endswith(";") and "[&trait=" in newick


def test_trait_requires_tree_and_out(tmp_path):
    with pytest.raises(SystemExit):                          # -t and -o are required
        main(["trait", "--model", "bm"])


def test_trait_mk_writes_discrete_states(tmp_path):
    tree = _tree_file(tmp_path, tips=8)
    dest = tmp_path / "tr"
    main(["trait", "-t", tree, "--model", "mk", "--states", "3", "--rate", "0.6",
          "--seed", "1", "-o", str(dest)])
    values = {r.split("\t")[1] for r in (dest / "traits.tsv").read_text().splitlines()[1:]}
    assert values <= {"0", "1", "2"}


def test_trait_reproducible(tmp_path):
    tree = _tree_file(tmp_path, tips=8)

    def run(name):
        dest = tmp_path / name
        main(["trait", "-t", tree, "--model", "bm", "--sigma2", "0.5", "--seed", "3",
              "-o", str(dest)])
        return (dest / "traits.tsv").read_text()

    assert run("a") == run("b")


def test_trait_replicates_write_wide_table(tmp_path):
    tree = _tree_file(tmp_path, tips=8)
    dest = tmp_path / "rep"
    rc = main(["trait", "-t", tree, "--model", "bm", "--sigma2", "0.5",
               "--replicates", "5", "--seed", "1", "-o", str(dest)])
    assert rc == 0
    rows = (dest / "traits.tsv").read_text().splitlines()
    assert rows[0] == "node\trep_1\trep_2\trep_3\trep_4\trep_5"
    assert all(len(r.split("\t")) == 6 for r in rows[1:])   # one column per replicate
    assert not (dest / "trait_tree.nwk").exists()           # wide-table mode: no annotated tree


def test_trait_mk_ordered(tmp_path):
    tree = _tree_file(tmp_path, tips=10)
    dest = tmp_path / "ord"
    rc = main(["trait", "-t", tree, "--model", "mk", "--states", "4", "--rate", "0.6",
               "--ordered", "--seed", "1", "-o", str(dest)])
    assert rc == 0
    vals = {r.split("\t")[1] for r in (dest / "traits.tsv").read_text().splitlines()[1:]}
    assert vals <= {"0", "1", "2", "3"}


def test_trait_mk_q_matrix(tmp_path):
    q = tmp_path / "q.tsv"
    q.write_text("0 2 0\n1 0 1\n0 3 0\n")               # a 3-state arbitrary Markov chain
    tree = _tree_file(tmp_path, tips=10)
    dest = tmp_path / "qm"
    rc = main(["trait", "-t", tree, "--model", "mk", "--q-matrix", str(q),
               "--seed", "1", "-o", str(dest)])
    assert rc == 0
    vals = {r.split("\t")[1] for r in (dest / "traits.tsv").read_text().splitlines()[1:]}
    assert vals <= {"0", "1", "2"}


def test_trait_dec_writes_ranges(tmp_path):
    tree = _tree_file(tmp_path, tips=10)
    dest = tmp_path / "dec"
    rc = main(["trait", "-t", tree, "--model", "dec", "--areas", "A,B,C",
               "--dispersal", "0.3", "--extinction", "0.1", "--root-range", "A",
               "--seed", "1", "-o", str(dest)])
    assert rc == 0
    vals = {r.split("\t")[0]: r.split("\t")[1]
            for r in (dest / "traits.tsv").read_text().splitlines()[1:]}
    assert vals["root"] == "{A}"                             # root range respected
    assert all(v.startswith("{") and v.endswith("}") for v in vals.values())   # ranges
    assert "[&trait={" in (dest / "trait_tree.nwk").read_text()
