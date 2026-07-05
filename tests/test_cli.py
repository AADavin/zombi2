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
    rc = main(["species", "--mode", "forward", "--age", "6", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()


def test_species_forward_requires_exactly_one_of_tips_age(tmp_path):
    """Forward needs exactly one of --tips/--age; neither or both is an error."""
    with pytest.raises(SystemExit):
        main(["species", "--mode", "forward", "-o", str(tmp_path / "a")])
    with pytest.raises(SystemExit):
        main(["species", "--mode", "forward", "--tips", "20", "--age", "5",
              "-o", str(tmp_path / "b")])


def test_species_forward_mass_extinction(tmp_path):
    """Forward mode accepts one or more --mass-extinction pulses and logs them."""
    out = tmp_path / "me"
    rc = main(["species", "--mode", "forward", "--birth", "1.2", "--death", "0.2",
               "--age", "6", "--mass-extinction", "2", "0.8",
               "--mass-extinction", "4", "0.5", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")
    log = (out / "species_tree.log").read_text()
    assert "mass_extinction" in log and "[2.0, 0.8]" in log


def test_species_mass_extinction_rejected_backward(tmp_path):
    """--mass-extinction is forward-only; on the default backward model it errors."""
    with pytest.raises(SystemExit):
        main(["species", "--mass-extinction", "2", "0.5", "-o", str(tmp_path / "a")])


def test_species_mass_extinction_needs_age_not_tips(tmp_path):
    """A pulse age needs a fixed present, so --mass-extinction requires --age (not --tips)."""
    with pytest.raises(SystemExit):
        main(["species", "--mode", "forward", "--tips", "20",
              "--mass-extinction", "2", "0.5", "-o", str(tmp_path / "b")])


def test_species_clads(tmp_path):
    """ClaDS (per-lineage rates) is a forward diversification process."""
    out = tmp_path / "clads"
    rc = main(["species", "--mode", "forward", "--diversification", "clads",
               "--birth", "1.0", "--clads-alpha", "0.9", "--clads-sigma", "0.2",
               "--turnover", "0.1", "--age", "5", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")


def test_species_diversity_dependent(tmp_path):
    """Diversity-dependent BD needs a carrying capacity and forward mode."""
    out = tmp_path / "dd"
    rc = main(["species", "--mode", "forward", "--diversification", "diversity-dependent",
               "--birth", "3", "--death", "0.2", "-K", "40", "--age", "20",
               "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")


def test_species_diversity_dependent_needs_K(tmp_path):
    """diversity-dependent without --carrying-capacity errors."""
    with pytest.raises(SystemExit):
        main(["species", "--mode", "forward", "--diversification", "diversity-dependent",
              "--age", "10", "-o", str(tmp_path / "a")])


def test_species_clads_requires_forward(tmp_path):
    """clads/diversity-dependent are forward-only; backward errors."""
    with pytest.raises(SystemExit):
        main(["species", "--diversification", "clads", "--tips", "20",
              "-o", str(tmp_path / "b")])


def test_species_clade_shift(tmp_path):
    """Forward mode accepts one or more --clade-shift regimes and logs them."""
    out = tmp_path / "cs"
    rc = main(["species", "--mode", "forward", "--birth", "0.9", "--death", "0.4",
               "--age", "5", "--clade-shift", "3.0", "1.6", "0.2",
               "--clade-shift", "1.5", "0.4", "0.6", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")
    assert "clade_shift" in (out / "species_tree.log").read_text()


def test_species_clade_shift_needs_age_not_tips(tmp_path):
    """Clade-shift times need a fixed present, so --clade-shift requires --age (not --tips)."""
    with pytest.raises(SystemExit):
        main(["species", "--mode", "forward", "--tips", "20",
              "--clade-shift", "3.0", "2.0", "0.1", "-o", str(tmp_path / "a")])


def test_species_clade_shift_rejected_backward(tmp_path):
    """--clade-shift is forward-only; on the default backward model it errors."""
    with pytest.raises(SystemExit):
        main(["species", "--clade-shift", "3.0", "2.0", "0.1", "-o", str(tmp_path / "b")])


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
def test_genomes_score_likelihoods_writes_table(tmp_path):
    """`genomes --score-likelihoods` writes the per-family ALE likelihood table."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.2", "--tips", "8", "--age", "1",
          "--seed", "5", "-o", str(sp)])
    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.3", "--trans", "0.2", "--loss", "0.4", "--initial-families", "12",
               "--seed", "5", "--write", "trees",
               "--score-likelihoods", "--score-model", "dated", "undated", "-o", str(gen)])
    assert rc == 0
    table = gen / "Reconciliation_likelihoods.tsv"
    assert table.exists()
    lines = table.read_text().strip().split("\n")
    assert lines[0] == "family\textant_copies\tdated_loglik\tundated_loglik"
    assert len(lines) >= 2  # header + at least one scored family


@needs_rust
def test_genomes_output_all_writes_full(tmp_path):
    """`--output all` writes the full ZOMBI-1 output (Rust engine)."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "--write", "all", "-o", str(gen)])
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
               "--seed", "42", "--write", "profiles", "-o", str(gen)])
    assert rc == 0
    assert (gen / "Profiles.tsv").exists()
    assert (gen / "Presence.tsv").exists()
    assert not (gen / "gene_trees").exists()
    assert not (gen / "Transfers.tsv").exists()


@needs_rust
def test_genomes_output_trace(tmp_path):
    """`--output trace` (optionally + profiles) writes the compact single-file event log via the
    fast path — Events_trace.tsv, no per-family event files, no gene trees."""
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3",
          "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])

    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.1", "--loss", "0.25", "--orig", "0.5",
               "--seed", "42", "--write", "trace", "profiles", "--sparse", "-o", str(gen)])
    assert rc == 0
    assert (gen / "Events_trace.tsv").exists()
    assert (gen / "Profiles_sparse.tsv").exists()
    assert not (gen / "gene_family_events").exists()   # the per-family dir is not written
    assert not (gen / "gene_trees").exists()            # trees are only reconstructed on request
    # header + at least the root originations
    lines = (gen / "Events_trace.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["time", "event", "branch", "donor", "recipient",
                                     "family", "parent", "child1", "child2"]
    assert len(lines) > 1


@needs_rust
def test_genomes_output_trace_with_trees(tmp_path):
    """`--output trace trees` writes the trace file *and* reconstructs the gene trees."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "30", "--seed", "1", "-o", str(sp)])
    gen = tmp_path / "gen"
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.2", "--trans", "0.1",
               "--loss", "0.25", "--orig", "0.5", "--seed", "1",
               "--write", "trace", "trees", "-o", str(gen)])
    assert rc == 0
    assert (gen / "Events_trace.tsv").exists()
    assert os.listdir(gen / "gene_trees")


@needs_rust
def test_genomes_output_selection(tmp_path):
    """`--output` writes exactly the requested components and nothing else."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "30", "--seed", "1", "-o", str(sp)])
    gen = tmp_path / "gen"
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.2", "--trans", "0.2",
               "--loss", "0.2", "--orig", "0.5", "--seed", "1",
               "--write", "trees", "transfers", "-o", str(gen)])
    assert rc == 0
    assert (gen / "gene_trees").exists() and (gen / "Transfers.tsv").exists()
    assert not (gen / "Profiles.tsv").exists()          # profiles not requested
    assert not (gen / "gene_family_events").exists()    # events not requested


def test_sparse_requires_profiles(tmp_path):
    """--sparse without 'profiles' in --output is a clean error (exit 1)."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.2", "--loss", "0.2",
               "--orig", "0.5", "--write", "trees", "--sparse", "-o", str(tmp_path / "g")])
    assert rc == 1


def test_forward_extinction_returns_clean_error(tmp_path, capsys):
    """A forward run that keeps going extinct exits 1 with a clean message (no traceback)."""
    rc = main(["species", "--mode", "forward", "--age", "8", "--birth", "1", "--death", "20",
               "--max-attempts", "5", "-o", str(tmp_path / "x")])
    assert rc == 1
    err = capsys.readouterr().err
    # the run banner is printed to stderr first; the clean one-line error follows it
    assert "zombi2: error:" in err and "extinct" in err


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
    rc = main(["species", "--mode", "forward", "--age", "6", "--fossilization", "0.3",
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


def _genomes_run_with_trace(tmp_path):
    """A helper: species -> genomes (genome-wise, with a written Events_trace.tsv). No Rust."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    run = tmp_path / "run"
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"), "--dup", "0.4", "--trans", "0.1",
               "--loss", "0.3", "--orig", "0.6", "--rate-model", "genome-wise",
               "--write", "trace", "profiles", "--seed", "2", "-o", str(run)])
    assert rc == 0
    assert (run / "Events_trace.tsv").exists()
    return run


def test_sequence_command_writes_phylograms(tmp_path):
    """`zombi2 sequence` replays a genomes run's trace and writes substitution-unit gene trees."""
    run = _genomes_run_with_trace(tmp_path)
    out = tmp_path / "seq"
    rc = main(["sequence", "--genomes", str(run), "--family-speed", "0.5",
               "--branch-speed", "0.4", "--seed", "7", "-o", str(out)])
    assert rc == 0
    subst = list((out / "gene_trees").glob("*_extant_subst.nwk"))
    assert subst
    assert all(p.read_text().strip().endswith(";") for p in subst)
    assert (out / "gene_family_speeds.tsv").read_text().startswith("family\tspeed")
    assert (out / "branch_rates.tsv").read_text().startswith("species_branch\trate")
    assert (out / "sequence.log").exists()


def test_sequence_command_discrete_bins(tmp_path):
    """`--branch-bins` selects the discrete-bin (GTDB) lineage clock."""
    run = _genomes_run_with_trace(tmp_path)
    out = tmp_path / "seq"
    rc = main(["sequence", "--genomes", str(run), "--family-speed", "0.5",
               "--branch-bins", "0.25,0.5,1,2,4", "--branch-switch-rate", "1.5",
               "--seed", "7", "-o", str(out)])
    assert rc == 0
    assert list((out / "gene_trees").glob("*_extant_subst.nwk"))
    assert (out / "branch_rates.tsv").exists()


def test_sequence_command_reproducible(tmp_path):
    """Same genomes run + same seed -> identical phylograms."""
    run = _genomes_run_with_trace(tmp_path)
    a, b = tmp_path / "a", tmp_path / "b"
    for out in (a, b):
        main(["sequence", "--genomes", str(run), "--family-speed", "0.5",
              "--branch-speed", "0.4", "--seed", "7", "-o", str(out)])
    fa = sorted((a / "gene_trees").glob("*_extant_subst.nwk"))
    fb = sorted((b / "gene_trees").glob("*_extant_subst.nwk"))
    assert fa and [p.name for p in fa] == [p.name for p in fb]
    assert all(x.read_text() == y.read_text() for x, y in zip(fa, fb))


def test_sequence_two_clocks_rejected(tmp_path):
    """--branch-speed and --branch-bins are two lineage clocks — giving both is an error."""
    run = _genomes_run_with_trace(tmp_path)
    rc = main(["sequence", "--genomes", str(run), "--branch-speed", "0.4",
               "--branch-bins", "0.5,1,2", "-o", str(tmp_path / "seq")])
    assert rc == 1


def test_sequence_missing_trace_is_clean_error(tmp_path):
    """A genomes run without a written trace gives a clear error, not a traceback."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "15", "--seed", "1", "-o", str(sp)])
    # a genomes dir with only a species tree (no Events_trace.tsv)
    rc = main(["sequence", "--genomes", str(sp), "--branch-speed", "0.4",
               "-o", str(tmp_path / "seq")])
    assert rc == 1


def _alignment_letters(aln_dir):
    """Read every FASTA in a `sequence` run's alignments/ dir; return the set of all letters seen."""
    from zombi2.sequence_sim import read_fasta
    files = sorted(aln_dir.glob("*.fasta"))
    letters = set()
    for f in files:
        for seq in read_fasta(str(f)).values():
            letters |= set(seq)
    return files, letters


def test_sequence_command_dna_alignments(tmp_path):
    """`zombi2 sequence --subst-model hky85` writes per-family FASTA over the ACGT alphabet."""
    run = _genomes_run_with_trace(tmp_path)
    out = tmp_path / "seq"
    rc = main(["sequence", "--genomes", str(run), "--subst-model", "hky85",
               "--seq-length", "60", "--branch-speed", "0.3", "--seed", "7", "-o", str(out)])
    assert rc == 0
    files, letters = _alignment_letters(out / "alignments")
    assert files                                            # at least one family alignment
    assert letters <= set("ACGT") and letters                # only nucleotides
    # the substitution-unit gene trees are still written alongside
    assert list((out / "gene_trees").glob("*_extant_subst.nwk"))


def test_sequence_command_protein_alignments(tmp_path):
    """`zombi2 sequence --subst-model lg` writes per-family FASTA over the 20-AA alphabet."""
    run = _genomes_run_with_trace(tmp_path)
    out = tmp_path / "seq"
    rc = main(["sequence", "--genomes", str(run), "--subst-model", "lg",
               "--seq-length", "50", "--gamma-shape", "0.5", "--seed", "7", "-o", str(out)])
    assert rc == 0
    files, letters = _alignment_letters(out / "alignments")
    assert files
    assert letters <= set("ARNDCQEGHILKMFPSTWYV")           # amino acids
    assert not letters <= set("ACGT")                       # genuinely protein
    # every sequence in a family alignment is --seq-length long
    from zombi2.sequence_sim import read_fasta
    recs = read_fasta(str(files[0]))
    assert recs and all(len(s) == 50 for s in recs.values())


def test_sequence_without_model_writes_no_alignments(tmp_path):
    """Omitting --subst-model keeps the old rescale-only behaviour (no alignments/ dir)."""
    run = _genomes_run_with_trace(tmp_path)
    out = tmp_path / "seq"
    rc = main(["sequence", "--genomes", str(run), "--branch-speed", "0.3",
               "--seed", "7", "-o", str(out)])
    assert rc == 0
    assert not (out / "alignments").exists()
    assert list((out / "gene_trees").glob("*_extant_subst.nwk"))


def test_sequence_gamma_without_model_rejected(tmp_path):
    """--gamma-shape without --subst-model is a clean error (nothing to apply it to)."""
    run = _genomes_run_with_trace(tmp_path)
    rc = main(["sequence", "--genomes", str(run), "--gamma-shape", "0.5",
               "-o", str(tmp_path / "seq")])
    assert rc == 1


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
    assert "initial_chromosomes\t1" in (out / "genomes.log").read_text()   # nucleotide default = 1


def test_genomes_nucleotide_profiles_only(tmp_path):
    """Nucleotide `--output profiles` writes only the emergent profile — no gene trees."""
    tree = _tree_file(tmp_path, tips=12)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide",
               "--dup", "0.0006", "--loss", "0.0006", "--write", "profiles",
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
               "--write", "profiles", "--sparse", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "Profiles_sparse.tsv").exists()
    assert not (out / "Profiles.tsv").exists()


def test_genomes_nucleotide_genic(tmp_path):
    """`--genes` enables genic mode: genes.tsv, split Gene/Intergene trees, Pseudogenizations,
    and atoms.tsv carrying the gene/intergene classification."""
    tree = _tree_file(tmp_path, tips=10)
    genes = tmp_path / "genes.tsv"
    genes.write_text("100\t180\tgeneA\n300\t360\tgeneB\n500\t620\tgeneC\n750\t800\tgeneD\n")
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide",
               "--inversion", "0.004", "--loss", "0.003", "--dup", "0.002", "--trans", "0.002",
               "--root-length", "1000", "--extension", "0.96", "--genes", str(genes),
               "--pseudogenization", "0.4", "--replacement", "0.6",
               "--write", "profiles", "trees", "--seed", "9", "-o", str(out)])
    assert rc == 0
    for f in ("genes.tsv", "atoms.tsv", "Profiles.tsv", "Pseudogenizations.tsv"):
        assert (out / f).exists()
    assert os.listdir(out / "Gene_trees") and os.listdir(out / "Intergene_trees")
    # atoms.tsv carries the classification; the four seed genes are present in genes.tsv
    assert "kind\tgene_id" in (out / "atoms.tsv").read_text()
    genes_txt = (out / "genes.tsv").read_text()
    for name in ("geneA", "geneB", "geneC", "geneD"):
        assert name in genes_txt


def test_genomes_nucleotide_genes_reject_profiles_only(tmp_path):
    """Genic mode needs the Python engine — a profiles-only run still writes the genic outputs
    (it never silently takes the Rust path)."""
    tree = _tree_file(tmp_path, tips=8)
    genes = tmp_path / "genes.tsv"
    genes.write_text("100 200 g1\n400 500 g2\n")
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide", "--loss", "0.002",
               "--root-length", "800", "--genes", str(genes), "--write", "profiles",
               "--seed", "2", "-o", str(out)])
    assert rc == 0
    assert (out / "genes.tsv").exists() and (out / "atoms.tsv").exists()


def test_genomes_nucleotide_from_gff(tmp_path):
    """`--gff` starts genic mode from a real annotation: its length + (trimmed) gene coordinates."""
    gff = tmp_path / "genome.gff"
    gff.write_text("##gff-version 3\n##sequence-region c 1 600\n"
                   "c\tx\tregion\t1\t600\t.\t+\t.\tID=c;Is_circular=true\n"
                   "c\tx\tgene\t50\t150\t.\t+\t.\tlocus_tag=a\n"
                   "c\tx\tgene\t140\t260\t.\t-\t.\tlocus_tag=b\n"       # overlaps a -> trimmed
                   "c\tx\tgene\t400\t500\t.\t+\t.\tlocus_tag=d\n")
    tree = _tree_file(tmp_path, tips=8)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide", "--gff", str(gff),
               "--loss", "0.003", "--dup", "0.002", "--pseudogenization", "0.3",
               "--write", "profiles", "trees", "--seed", "5", "-o", str(out)])
    assert rc == 0
    genes_txt = (out / "genes.tsv").read_text()
    for locus in ("a", "b", "d"):
        assert locus in genes_txt
    assert (out / "atoms.tsv").exists() and os.listdir(out / "Gene_trees")
    assert "root_length\t600" in (out / "genomes.log").read_text()   # length came from the GFF


def test_genomes_gff_and_genes_conflict(tmp_path, capsys):
    """--gff and --genes both set the gene coordinates — giving both is a clean error."""
    gff = tmp_path / "g.gff"
    gff.write_text("##sequence-region c 1 100\nc\tx\tgene\t1\t10\t.\t+\t.\tlocus_tag=a\n")
    genes = tmp_path / "genes.tsv"
    genes.write_text("1 10 a\n")
    tree = _tree_file(tmp_path, tips=6)
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide", "--gff", str(gff),
               "--genes", str(genes), "--loss", "0.002", "-o", str(tmp_path / "nt")])
    assert rc == 1
    assert "either --gff or --genes" in capsys.readouterr().err


def test_genomes_nucleotide_ancestral(tmp_path):
    """`--output ancestral` simulates DNA and writes the genome (architecture + gzipped FASTA) at
    every node; with --genome-fasta the reconstructed root == the input genome."""
    import gzip
    import numpy as np
    genes = tmp_path / "genes.tsv"
    genes.write_text("20 60 gA\n90 130 gB\n160 200 gC\n230 280 gD\n")
    genome = "".join(np.random.default_rng(1).choice(list("ACGT"), size=300))
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">seq\n" + genome + "\n")
    tree = _tree_file(tmp_path, tips=6)
    out = tmp_path / "nt"
    rc = main(["genomes", "-t", tree, "--rate-model", "nucleotide", "--genes", str(genes),
               "--root-length", "300", "--inversion", "0.01", "--loss", "0.006", "--dup", "0.005",
               "--subst-model", "hky85", "--subst-rate", "0.4", "--gamma-shape", "0.5",
               "--genome-fasta", str(fasta), "--write", "ancestral", "--seed", "5", "-o", str(out)])
    assert rc == 0
    assert os.listdir(out / "Architecture") and os.listdir(out / "Genomes")
    assert os.listdir(out / "Gene_alignments")
    # the root genome FASTA reconstructs the input exactly
    with gzip.open(out / "Genomes" / "root.fasta.gz", "rt") as fh:
        root_seq = "".join(l.strip() for l in fh if not l.startswith(">"))
    assert root_seq == genome
    # root architecture keeps the four genes intact and in order
    arch = (out / "Architecture" / "root.tsv").read_text()
    assert arch.count("\tgene\t") == 4


# --- abc: fit gene-family rates to an empirical profile by ABC inference --

def _profile_file(tmp_path, tree, **rates):
    """Simulate an empirical profile with known rates (uniform model) -> its Profiles.tsv."""
    dest = tmp_path / "truth"
    args = ["genomes", "-t", tree, "--write", "profiles", "--seed", "7", "-o", str(dest)]
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


# --------------------------------------------------------------------------- coevolve (traits:species = SSE)
def test_coevolve_traits_species_bisse(tmp_path):
    """`coevolve --couple traits:species` grows a tree jointly with a binary trait (BiSSE)."""
    out = tmp_path / "cv"
    rc = main(["coevolve", "--couple", "traits:species", "--sse-model", "bisse",
               "--lambda0", "1", "--lambda1", "2", "--mu0", "0.2", "--mu1", "0.2",
               "--q01", "0.1", "--q10", "0.1", "--tips", "20", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")
    assert (out / "trait_tree.nwk").exists()
    assert (out / "coevolve.log").exists()
    # traits.tsv covers every node (tips + ancestors), values in {0, 1}
    rows = [ln.split("\t") for ln in (out / "traits.tsv").read_text().splitlines()[1:]]
    assert rows and all(v in ("0", "1") for _, v in rows)


def test_coevolve_default_edge_is_traits_species(tmp_path):
    """With no --couple, the default edge is traits:species."""
    out = tmp_path / "cv"
    rc = main(["coevolve", "--age", "3", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()


def test_coevolve_quasse(tmp_path):
    """--sse-model quasse grows a tree with a continuous trait; traits are floats."""
    out = tmp_path / "cv"
    rc = main(["coevolve", "--couple", "traits:species", "--sse-model", "quasse",
               "--spec-low", "0.5", "--spec-high", "2", "--qmu", "0.2", "--diffusion", "0.5",
               "--age", "3", "--seed", "2", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()


def test_coevolve_musse_needs_q_matrix(tmp_path):
    """musse requires --birth/--death/--q-matrix."""
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "--sse-model", "musse",
              "--birth", "1", "1", "--death", "0.1", "0.1", "--tips", "10", "-o", str(tmp_path / "a")])


def test_coevolve_reproducible(tmp_path):
    """Same seed -> identical tree."""
    a, b = tmp_path / "a", tmp_path / "b"
    args = ["coevolve", "--couple", "traits:species", "--tips", "15", "--seed", "7", "-o"]
    main(args + [str(a)])
    main(args + [str(b)])
    assert (a / "species_tree.nwk").read_text() == (b / "species_tree.nwk").read_text()


def test_coevolve_unbuilt_edge_errors(tmp_path):
    """A planned-but-unbuilt edge errors clearly (does not silently run)."""
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "species:genes", "--age", "3", "-o", str(tmp_path / "a")])


def test_coevolve_traits_genes_on_given_tree(tmp_path):
    """traits:genes (folded in from the old coevolve-genetrait command) evolves a trait, then a
    gene-family panel conditioned on it, along a GIVEN tree."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "24", "--age", "4", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "tg"
    rc = main(["coevolve", "--couple", "traits:genes", "-t", str(sp / "species_tree.nwk"),
               "--trait-model", "mk", "--states", "2", "--trait-center", "--panel", "24",
               "--responsive", "0.5", "--effect-loss", "3", "--write", "profiles", "trees",
               "--seed", "7", "-o", str(out)])
    assert rc == 0
    for name in ("Profiles.tsv", "traits.tsv", "trait_tree.nwk", "coupling.tsv",
                 "coevolve.log", "species_tree.nwk"):
        assert (out / name).exists(), name
    assert (out / "gene_trees").is_dir()


def test_coevolve_traits_genes_rejects_grow_flags(tmp_path):
    """traits:genes runs on a given tree; --age/--tips (which grow a tree) are an error."""
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:genes", "--age", "3", "-o", str(tmp_path / "a")])


def test_coevolve_needs_exactly_one_stop_condition(tmp_path):
    """traits:species grows the tree: neither or both of --age/--tips is an error."""
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "-o", str(tmp_path / "a")])
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "--age", "3", "--tips", "10",
              "-o", str(tmp_path / "b")])


def test_coevolve_bad_edge_name_errors(tmp_path):
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "foo:bar", "--age", "3", "-o", str(tmp_path / "a")])


# --------------------------------------------------------------------------- coevolve Phase 2: species:traits / ClaSSE
def test_coevolve_species_traits_on_given_tree(tmp_path):
    """species:traits alone evolves a cladogenetic trait along a GIVEN tree (no diversification)."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "30", "--age", "4", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "cv"
    rc = main(["coevolve", "--couple", "species:traits", "-t", str(sp / "species_tree.nwk"),
               "--sse-model", "bisse", "--q01", "0", "--q10", "0", "--clado-shift", "0.4",
               "--seed", "2", "-o", str(out)])
    assert rc == 0
    assert (out / "traits.tsv").exists() and (out / "trait_tree.nwk").exists()
    assert (out / "species_tree.nwk").exists()       # the given tree, copied for provenance


def test_coevolve_classe_both_arrows(tmp_path):
    """traits:species + species:traits = ClaSSE (grows a tree with cladogenetic jumps)."""
    out = tmp_path / "cv"
    rc = main(["coevolve", "--couple", "traits:species", "--couple", "species:traits",
               "--lambda0", "1", "--lambda1", "3", "--q01", "0.05", "--q10", "0.05",
               "--clado-shift", "0.3", "--tips", "80", "--seed", "3", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")


def test_coevolve_couple_space_separated_list(tmp_path):
    """--couple also accepts a space-separated list in a single flag."""
    out = tmp_path / "cv"
    rc = main(["coevolve", "--couple", "traits:species", "species:traits",
               "--tips", "60", "--seed", "4", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()


def test_coevolve_species_traits_alone_needs_tree(tmp_path):
    """species:traits alone runs on a given tree; without -t it errors."""
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "species:traits", "--clado-shift", "0.4",
              "-o", str(tmp_path / "a")])


def test_coevolve_into_species_rejects_input_tree(tmp_path):
    """An into-species edge grows the tree, so passing -t is an error."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "-t", str(sp / "species_tree.nwk"),
              "--tips", "20", "-o", str(tmp_path / "a")])


def test_coevolve_quasse_classe(tmp_path):
    """Continuous ClaSSE: quasse + a cladogenetic jump."""
    out = tmp_path / "cv"
    rc = main(["coevolve", "--couple", "traits:species", "species:traits", "--sse-model", "quasse",
               "--spec-low", "0.4", "--spec-high", "3", "--qmu", "0.2", "--diffusion", "0",
               "--clado-jump", "1.0", "--tips", "60", "--seed", "5", "-o", str(out)])
    assert rc == 0
    assert (out / "trait_tree.nwk").exists()


# --------------------------------------------------------------------------- coevolve Phase 3: genes:species
def test_coevolve_genes_species(tmp_path):
    """genes:species grows a tree driven by key-innovation gene families; writes drivers + manifest."""
    out = tmp_path / "gs"
    rc = main(["coevolve", "--couple", "genes:species", "--drivers", "2",
               "--lambda0", "1", "--mu0", "0.2", "--driver-speciation", "1.2",
               "--driver-transfer", "0.8", "--driver-loss", "0.3", "--root-drivers", "1",
               "--tips", "120", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").read_text().strip().endswith(";")
    assert (out / "drivers_manifest.tsv").exists()
    header = (out / "drivers.tsv").read_text().splitlines()[0]
    assert header == "node\tD0\tD1"


def test_coevolve_genes_species_then_genomes_overlay(tmp_path):
    """The neutral genome overlays on the grown tree with the ordinary genomes command."""
    out = tmp_path / "gs"
    main(["coevolve", "--couple", "genes:species", "--drivers", "1", "--root-drivers", "1",
          "--tips", "60", "--seed", "2", "-o", str(out)])
    ov = tmp_path / "ov"
    rc = main(["genomes", "-t", str(out / "species_tree.nwk"), "--trans", "1", "--loss", "0.5",
               "--write", "profiles", "-o", str(ov)])
    assert rc == 0
    assert (ov / "Profiles.tsv").exists()


def test_coevolve_genes_species_rejects_tree(tmp_path):
    """genes:species grows the tree, so -t is an error."""
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "genes:species", "-t", str(sp / "species_tree.nwk"),
              "--tips", "20", "-o", str(tmp_path / "a")])


def test_coevolve_genes_species_not_combinable_yet(tmp_path):
    """Combining genes:species with another edge (the joint model) is not in this phase."""
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "genes:species", "--couple", "traits:species",
              "--tips", "20", "-o", str(tmp_path / "a")])


# ── the redesigned CLI surface: --version, grouped help, and the renamed flags ──────
def test_version_flag_prints_version(capsys):
    """`zombi2 --version` prints the version and exits 0 (argparse version action)."""
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert "ZOMBI2" in capsys.readouterr().out


def test_top_level_help_lists_commands_grouped(capsys):
    """The top-level help carries the banner and the curated, theme-grouped command list
    (not a duplicate argparse positional dump)."""
    with pytest.raises(SystemExit):
        main(["--help"])
    out = capsys.readouterr().out
    assert "a simulator of species trees, genomes, traits and sequences" in out   # banner
    assert "Species trees" in out and "Inference" in out
    for cmd in ("species", "genomes", "trait", "abc", "coevolve", "sequence"):
        assert cmd in out
    assert "==SUPPRESS==" not in out                   # the auto command dump is hidden


def test_genomes_help_is_sectioned_by_model(capsys):
    """`genomes --help` groups its options into UPPERCASE model sections (IQ-TREE style),
    and advertises --write (not the old --output)."""
    with pytest.raises(SystemExit):
        main(["genomes", "--help"])
    out = capsys.readouterr().out
    for section in ("GENERAL", "GENE-FAMILY RATES", "OUTPUT",
                    "NUCLEOTIDE MODEL", "GENES & INTERGENES"):
        assert section in out
    assert "--write" in out and "--output" not in out


def test_renamed_flags_work_and_old_names_rejected(tmp_path):
    """The renamed flags parse; the old names are gone (a clean break, as designed)."""
    out = tmp_path / "sp"
    assert main(["species", "--mode", "backward", "--tips", "8", "--seed", "1",
                 "-o", str(out)]) == 0
    with pytest.raises(SystemExit):                    # old species --model is gone
        main(["species", "--model", "backward", "-o", str(tmp_path / "old")])
