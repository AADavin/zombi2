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
    for f in ("species_tree.nwk", "Profiles.tsv"):     # default output = profiles + trees
        assert (gen / f).exists()
    assert os.listdir(gen / "gene_trees")


@needs_rust
def test_genomes_output_all_writes_full(tmp_path):
    """`--output all` writes the full ZOMBI-1 output (Rust engine)."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "--output", "all", "-o", str(gen)])
    assert rc == 0
    for f in ("species_tree.nwk", "Profiles.tsv", "Presence.tsv", "Transfers.tsv",
              "Gene_family_summary.tsv"):
        assert (gen / f).exists()
    assert os.listdir(gen / "gene_trees")
    assert os.listdir(gen / "gene_family_events")


@needs_rust
def test_genomes_output_profiles(tmp_path):
    """`--output profiles` writes just the profile matrices — no gene trees / event log."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "--output", "profiles", "-o", str(gen)])
    assert rc == 0
    assert (gen / "Profiles.tsv").exists()
    assert (gen / "Presence.tsv").exists()
    assert not (gen / "gene_trees").exists()
    assert not (gen / "Transfers.tsv").exists()


@needs_rust
def test_genomes_output_selection(tmp_path):
    """`--output` writes exactly the requested components and nothing else."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "30", "--seed", "1", "-o", str(sp)])
    gen = tmp_path / "gen"
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.2", "--trans", "0.2",
               "--loss", "0.2", "--orig", "0.5", "--seed", "1",
               "--output", "trees", "transfers", "-o", str(gen)])
    assert rc == 0
    assert (gen / "gene_trees").exists() and (gen / "Transfers.tsv").exists()
    assert not (gen / "Profiles.tsv").exists()          # profiles not requested
    assert not (gen / "gene_family_events").exists()    # events not requested


def test_sparse_requires_profiles(tmp_path):
    """--sparse without 'profiles' in --output is a clean error (exit 1)."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.2", "--loss", "0.2",
               "--orig", "0.5", "--output", "trees", "--sparse", "-o", str(tmp_path / "g")])
    assert rc == 1


def test_forward_extinction_returns_clean_error(tmp_path, capsys):
    """A forward run that keeps going extinct exits 1 with a clean message (no traceback)."""
    rc = main(["species", "--model", "forward", "--age", "8", "--birth", "1", "--death", "20",
               "--max-attempts", "5", "-o", str(tmp_path / "x")])
    assert rc == 1
    err = capsys.readouterr().err
    assert err.startswith("zombi2: error:") and "extinct" in err


def test_genomes_missing_tree_file_returns_clean_error(tmp_path, capsys):
    """A missing --tree file exits 1 with a clean message, not a traceback."""
    rc = main(["genomes", "--tree", str(tmp_path / "nope.nwk"),
               "--dup", "0.2", "-o", str(tmp_path / "g")])
    assert rc == 1
    assert "zombi2: error:" in capsys.readouterr().err


def test_species_episodic(tmp_path):
    """Multiple --birth/--death with --shifts builds an episodic (skyline) model."""
    out = tmp_path / "ep"
    rc = main(["species", "--birth", "1", "2", "--death", "0.3", "0.1", "--shifts", "2",
               "--age", "5", "--tips", "30", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")


def test_species_ghosts_adds_extinct_tips(tmp_path):
    """--ghosts un-prunes the backward tree, adding extinct (e*) leaves."""
    import zombi2 as z
    out = tmp_path / "gh"
    rc = main(["species", "--tips", "30", "--death", "0.6", "--ghosts", "--seed", "1",
               "-o", str(out)])
    assert rc == 0
    leaves = z.read_newick((out / "species_tree.nwk").read_text()).leaves()
    assert any(n.name.startswith("e") for n in leaves)   # extinct/ghost tips are named e*


def test_species_forward_fossils(tmp_path):
    """Forward + fossilization runs (fossilized birth–death)."""
    rc = main(["species", "--model", "forward", "--age", "6", "--fossilization", "0.3",
               "--seed", "1", "-o", str(tmp_path / "fbd")])
    assert rc == 0


def test_species_backward_fossils_is_error(tmp_path):
    """Fossil / removal / sampling flags require forward mode."""
    with pytest.raises(SystemExit):
        main(["species", "--fossilization", "0.2", "-o", str(tmp_path / "x")])


def test_log_written_by_default(tmp_path):
    """Every run always writes species_tree.log with the full set of parameters."""
    out = tmp_path / "sp"
    main(["species", "--tips", "15", "--seed", "3", "-o", str(out)])
    log = (out / "species_tree.log").read_text()
    assert "zombi2_version" in log and "seed\t3" in log and "model\tbackward" in log


@needs_rust
def test_genomes_annotate_species(tmp_path):
    """--annotate-species labels internal gene nodes <gid>|<species-branch>."""
    import re
    sp = tmp_path / "sp"
    main(["species", "--tips", "15", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "g"
    main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.3", "--loss", "0.3",
          "--orig", "0.5", "--max-family-size", "0.5", "--annotate-species", "--seed", "1",
          "-o", str(out)])
    trees = "".join(p.read_text() for p in (out / "gene_trees").glob("*.nwk"))
    assert re.search(r"g\d+\|", trees)                     # gid|branch labels present
    # without the flag, internal nodes are bare gids
    out2 = tmp_path / "g2"
    main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.3", "--loss", "0.3",
          "--orig", "0.5", "--max-family-size", "0.5", "--seed", "1", "-o", str(out2)])
    trees2 = "".join(p.read_text() for p in (out2 / "gene_trees").glob("*.nwk"))
    assert "|" not in trees2


def test_genomes_genome_wise_rate_model(tmp_path):
    """--rate-model genome-wise runs (Python engine) and is recorded in the log."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "gw"
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.5", "--loss", "0.4",
               "--orig", "0.5", "--rate-model", "genome-wise", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert "rate_model\tgenome-wise" in (out / "genomes.log").read_text()


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


def test_trait_writes_log(tmp_path):
    """`trait` always writes trait.log with the command line and the full parameter set."""
    tree = _tree_file(tmp_path, tips=10)
    dest = tmp_path / "tr"
    rc = main(["trait", "-t", tree, "--model", "bm", "--seed", "5", "-o", str(dest)])
    assert rc == 0
    log = (dest / "trait.log").read_text()
    assert "zombi2_version" in log and "command_line\t" in log
    assert "model\tbm" in log and "seed\t5" in log


# --- nucleotide genome model (--rate-model nucleotide): structural events, atoms as genes --

def test_genomes_nucleotide_model(tmp_path):
    """`--rate-model nucleotide` evolves nucleotide genomes; profiles+trees writes the atom
    table, the emergent profile, per-atom gene trees, and reconciliations (Python engine)."""
    tree = _tree_file(tmp_path, tips=12)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide",
               "--dup", "0.0006", "--loss", "0.0006", "--seed", "1", "-o", str(out)])
    assert rc == 0
    for f in ("species_tree.nwk", "atoms.tsv", "Profiles.tsv", "Presence.tsv",
              "Mosaics.tsv", "Reconciled_complete.nwk", "Reconciliation_events.tsv"):
        assert (out / f).exists()
    assert os.listdir(out / "gene_trees")                    # one gene tree per atom
    assert "rate_model\tnucleotide" in (out / "genomes.log").read_text()
    assert "initial_size\t1" in (out / "genomes.log").read_text()   # nucleotide default = 1


def test_genomes_nucleotide_profiles_only(tmp_path):
    """Nucleotide `--output profiles` writes only the emergent profile — no gene trees."""
    tree = _tree_file(tmp_path, tips=12)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide",
               "--dup", "0.0006", "--loss", "0.0006", "--output", "profiles",
               "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "Profiles.tsv").exists() and (out / "atoms.tsv").exists()
    assert not (out / "gene_trees").exists()
    assert not (out / "Reconciled_complete.nwk").exists()


def test_genomes_nucleotide_sparse(tmp_path):
    """Nucleotide profiles honour --sparse (single Profiles_sparse.tsv long table)."""
    tree = _tree_file(tmp_path, tips=12)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide", "--loss", "0.0008",
               "--output", "profiles", "--sparse", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "Profiles_sparse.tsv").exists()
    assert not (out / "Profiles.tsv").exists()


# --- abc: fit gene-family rates to an empirical profile by ABC inference --

def _profile_file(tmp_path, tree, **rates):
    """Simulate an empirical profile with known rates (uniform model) -> its Profiles.tsv."""
    dest = tmp_path / "truth"
    args = ["genomes", "-t", tree, "--output", "profiles", "--seed", "7", "-o", str(dest)]
    for k, v in rates.items():
        args += [f"--{k}", str(v)]
    main(args)
    return str(dest / "Profiles.tsv")


@needs_rust
def test_abc_fits_and_writes_posterior(tmp_path):
    """`abc` fits the rates given as priors and writes posterior + summary + spectra + log."""
    tree = _tree_file(tmp_path, tips=25)
    prof = _profile_file(tmp_path, tree, dup=0.3, loss=0.6, orig=2.0)
    out = tmp_path / "fit"
    rc = main(["abc", "-t", tree, "--profiles", prof,
               "--dup", "0", "1", "--loss", "0", "1.5", "--orig", "0", "4",
               "--n-sims", "200", "--seed", "1", "-o", str(out)])
    assert rc == 0
    for f in ("posterior.tsv", "summary.tsv", "spectra.tsv", "abc.log"):
        assert (out / f).exists()
    summ = (out / "summary.tsv").read_text()
    assert "duplication" in summ and "loss" in summ and "origination" in summ
    assert "transfer" not in summ                          # not given a prior -> not fitted
    header = (out / "posterior.tsv").read_text().splitlines()[0].split("\t")
    assert header == ["duplication", "loss", "origination"]
    assert len((out / "posterior.tsv").read_text().splitlines()) > 1     # >= 1 accepted draw


def test_abc_requires_a_range_to_fit(tmp_path, capsys):
    """All-fixed priors (no range) is a clean error, exit 1 — nothing to fit."""
    tree = _tree_file(tmp_path, tips=15)
    prof = tmp_path / "p.tsv"
    prof.write_text("family\tn1\tn2\nf1\t1\t0\n")
    rc = main(["abc", "-t", tree, "--profiles", str(prof), "--dup", "0.3",
               "-o", str(tmp_path / "o")])
    assert rc == 1
    assert "range" in capsys.readouterr().err


@needs_rust
def test_abc_regression_adjust(tmp_path):
    """--regression-adjust adds adjusted rows to summary.tsv."""
    tree = _tree_file(tmp_path, tips=25)
    prof = _profile_file(tmp_path, tree, dup=0.3, loss=0.6, orig=2.0)
    out = tmp_path / "ra"
    main(["abc", "-t", tree, "--profiles", prof, "--dup", "0", "1", "--loss", "0", "1.5",
          "--orig", "0", "4", "--n-sims", "200", "--regression-adjust", "--seed", "1",
          "-o", str(out)])
    assert "duplication_adj" in (out / "summary.tsv").read_text()


@needs_rust
def test_abc_smc(tmp_path):
    """--smc runs the sequential sampler and writes a posterior."""
    tree = _tree_file(tmp_path, tips=25)
    prof = _profile_file(tmp_path, tree, dup=0.3, loss=0.6)
    out = tmp_path / "smc"
    rc = main(["abc", "-t", tree, "--profiles", prof, "--dup", "0", "1", "--loss", "0", "1.5",
               "--smc", "--rounds", "2", "--particles", "60", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "posterior.tsv").exists()
