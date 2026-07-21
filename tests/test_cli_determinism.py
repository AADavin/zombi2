"""Same seed, same files — the whole pipeline, end to end.

Ported from Krister Swenson's fork (``thekswenson/Zombi``, ``tests/test_randomization.py``), which
runs each mode twice with a fixed seed and compares the two output directories with
``filecmp.dircmp``.

zombi2 already checks reproducibility per command (``test_cli`` has a
``test_<level>_is_deterministic_given_the_seed`` for species, sequences and traits). What it does
not have — and what the fork's test is really about — is the **whole pipeline** compared **as
files**: run species → genomes → sequences → traits twice, and diff every byte written. That is
additive in two ways. It covers the handoffs, where one command reads another's output and a
non-determinism would only show downstream; and it compares the artifacts a user actually keeps,
so a stray timestamp or a set iterated in hash order is caught wherever it lands, not only in the
outputs a per-command test happens to look at.

``.log`` files are excluded: they record a wall-clock ``timestamp`` by design, so they differ
between two runs however deterministic the simulation is. Everything else must match exactly.
"""

import pytest

from zombi2.cli.main import main

SEED = "20"


def _pipeline(root):
    """Run species → genomes (all three resolutions) → sequences → traits into ``root``."""
    tree = str(root / "species_complete.nwk")
    assert main(["species", "--birth", "1.0", "--death", "0.3", "--n-extant", "12",
                 "--seed", SEED, "-o", str(root)]) == 0

    assert main(["genomes", "-t", tree, "--duplication", "0.3", "--loss", "0.25",
                 "--origination", "0.6", "--seed", SEED, "-o", str(root / "g_unordered")]) == 0
    assert main(["genomes", "-t", tree, "--resolution", "ordered", "--duplication", "0.3",
                 "--loss", "0.25", "--origination", "0.6", "--transfer", "0.2",
                 "--inversion", "0.4", "--transposition", "0.3", "--chromosomes", "2",
                 "--seed", SEED, "-o", str(root / "g_ordered"),
                 "--write", "events", "profiles", "gene_order", "rearrangements",
                 "chromosome_events", "event_positions"]) == 0
    assert main(["genomes", "-t", tree, "--resolution", "nucleotide", "--root-length", "600",
                 "--genes", "4", "--inversion", "0.8", "--duplication", "0.4", "--loss", "0.3",
                 "--seed", SEED, "-o", str(root / "g_nucleotide"),
                 "--write", "events", "genes", "blocks", "rearrangements"]) == 0

    assert main(["sequences", "--genomes", str(root / "g_unordered"), "--model", "hky85",
                 "--length", "150", "--seed", SEED, "-o", str(root / "s"),
                 "--write", "alignments", "phylograms", "ancestral",
                 "species_phylogram"]) == 0
    assert main(["traits", "-t", tree, "--rate", "1.0", "--seed", SEED,
                 "-o", str(root / "t")]) == 0


def _artifacts(root):
    """``{relative path: bytes}`` for everything written, bar the run logs (they stamp the wall
    clock, so they differ between two runs however deterministic the simulation)."""
    return {str(p.relative_to(root)): p.read_bytes()
            for p in sorted(root.rglob("*")) if p.is_file() and p.suffix != ".log"}


@pytest.fixture(scope="module")
def two_runs(tmp_path_factory):
    base = tmp_path_factory.mktemp("determinism")
    first, second = base / "run1", base / "run2"
    for root in (first, second):
        root.mkdir()
        _pipeline(root)
    return _artifacts(first), _artifacts(second)


def test_the_pipeline_writes_the_same_files_twice(two_runs):
    first, second = two_runs
    assert sorted(first) == sorted(second)
    assert first, "the pipeline should have written something to compare"


def test_every_written_byte_is_identical(two_runs):
    first, second = two_runs
    differing = [name for name in first if first[name] != second[name]]
    assert not differing, f"same seed, different output: {differing}"


def test_the_comparison_covers_every_level(two_runs):
    # a guard on the fixture rather than the engine: if a command stopped writing, or a level were
    # dropped from the pipeline above, the two runs would still agree — vacuously
    names = " ".join(two_runs[0])
    for expected in ("species_complete.nwk", "genome_events.tsv", "gene_order.tsv",
                     "genome_event_positions.tsv", "blocks.tsv", "genes.tsv",
                     "sequences_alignment", "trait_values.tsv"):
        assert expected in names, f"{expected} is missing — the pipeline is not covering that level"


def test_excluding_the_run_logs_is_justified_and_narrow(tmp_path):
    # the exclusion has to be justified, not assumed. A log carries a wall-clock `timestamp`, which
    # is why it cannot join a byte comparison; everything else in it must still agree between two
    # same-seed runs. (Two runs a second apart would show the timestamp differing too, but asserting
    # that would make this test depend on how fast the machine is.)
    runs = []
    for tag in ("a", "b"):
        root = tmp_path / tag
        root.mkdir()
        main(["species", "--birth", "1.0", "--death", "0.3", "--n-extant", "10", "--seed", SEED,
              "-o", str(root)])
        runs.append((root / "species.log").read_text().splitlines())

    first, second = runs
    assert any(line.startswith("timestamp\t") for line in first), \
        "no timestamp in the run log — the exclusion would be unnecessary"
    # `output` differs by construction: the two runs write to different directories
    differing = [(a, b) for a, b in zip(first, second, strict=True) if a != b]
    assert all(a.startswith(("timestamp\t", "output\t")) for a, _ in differing), \
        f"a run log differs by more than its timestamp and path: {differing}"
