"""Tests for the run summary — what came out, beside the log that says what went in.

Three things nobody could get from a run's output, each of which a reader reconstructed by hand and
one of them wrongly:

- the event log holds one row per gene-tree **edge**, so duplications, transfers and speciations each
  write two. Counting by row inflates them exactly 2×;
- ``gene_trees/`` holds a pair of files per family that *ever existed* while the summary line counts
  the families that **survived**, and neither number explained the other;
- the family-size cap was invisible: when it binds it discards events, so realised rates fall below
  the declared ones, and nothing said so.

The first is the one that matters most, because it produces a plausible wrong number rather than an
error, so it is checked against the log's own ``event`` column rather than against itself.
"""

import collections
import csv
import json

import pytest

from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85
from zombi2.species import simulate_species_tree


def _run(**kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1)
    return sp, simulate_genomes_family(sp, **{"duplication": 0.4, "transfer": 0.2, "loss": 0.2,
                                              "origination": 0.4, "initial_families": 15,
                                              "seed": 1, **kw})


def test_the_event_counts_are_deduplicated_against_the_logs_own_event_column(tmp_path):
    # the whole point. A duplication, a transfer and a speciation each write TWO rows, so a reader
    # counting by row doubles them — the mistake that produces a wrong number instead of a crash.
    # Checked against the `event` column the log itself writes, not against the summary's own logic.
    sp, g = _run()
    g.write(tmp_path, outputs=("events", "summary"))
    rows = list(csv.DictReader((tmp_path / "genome_events.tsv").read_text(
        encoding="utf-8").splitlines(), delimiter="\t"))
    from_log = collections.defaultdict(set)
    for r in rows:
        from_log[r["kind"]].add(r["event"])

    s = json.loads((tmp_path / "genome_summary.json").read_text(encoding="utf-8"))
    for kind in ("duplication", "transfer", "loss", "speciation"):
        assert s["events"][kind] == len(from_log[kind]), kind
    # origination is split: the initial genome is logged as origination at the root's own start time,
    # so a bare count of that kind is de-novo arrivals PLUS initial_families
    assert s["events"]["initial"] + s["events"]["origination"] == len(from_log["origination"])
    assert s["events"]["initial"] == 15                      # what --initial-families asked for
    # and the row count is reported too, so the 2x is visible rather than a trap
    assert s["event_rows"] == len(rows) > sum(s["events"].values())


def test_the_family_counts_explain_the_gene_tree_file_count(tmp_path):
    # "96 gene families" beside 213 files, with nothing to reconcile them: born vs surviving is the
    # missing sentence, so the summary says both.
    sp, g = _run()
    g.write(tmp_path, outputs=("gene_trees", "summary"))
    s = json.loads((tmp_path / "genome_summary.json").read_text(encoding="utf-8"))
    trees = tmp_path / "gene_trees"
    assert len(list(trees.glob("*_complete.nwk"))) == s["families"]["born"]
    assert len(list(trees.glob("*_extant.nwk"))) == s["families"]["surviving"]
    assert s["families"]["died_out"] == s["families"]["born"] - s["families"]["surviving"]


def test_the_cap_reports_the_families_sitting_at_it():
    # the cap was invisible. A family at the ceiling had duplications and arriving transfers
    # discarded, so its realised rates are below the declared ones — this is how a reader finds out.
    _, tight = _run(duplication=1.2, loss=0.05, max_family_size=4)
    s = tight.summary()["family_size_cap"]
    assert s["cap"] == 4
    assert s["families_at_cap"] > 0 and s["cells_at_cap"] >= s["families_at_cap"]
    assert len(s["family_ids_at_cap"]) == s["families_at_cap"]
    # and every family it names really is at the cap in some genome
    for fam in s["family_ids_at_cap"]:
        assert max(tight.family_counts(n.id)[fam]
                   for n in tight.complete_tree.extant()) == 4

    _, free = _run(duplication=1.2, loss=0.05, max_family_size=None)
    lifted = free.summary()["family_size_cap"]
    assert lifted == {"cap": None, "families_at_cap": 0, "cells_at_cap": 0,
                      "family_ids_at_cap": []}


def test_a_cap_that_never_bites_says_so():
    # the common case, and the one that must not cry wolf: realistic rates nowhere near the ceiling
    _, g = _run(duplication=0.02, loss=0.1, max_family_size=10)
    assert g.summary()["family_size_cap"]["families_at_cap"] == 0


def test_the_species_summary_counts_every_tip_once():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, sampling=0.7, fossils=0.2, seed=1)
    s = r.summary()
    tips = s["tips"]
    assert tips["extant"] + tips["extinct"] + tips["unsampled"] == tips["total"]
    assert tips["extant"] == r.n_extant
    assert s["events"]["speciation"] == sum(1 for e in r.events if e.kind == "speciation")
    assert s["fossils"] == len(r.fossils)
    # the realised rates are events over the exposure that produced them, which is what a declared
    # per-lineage rate means — the cheapest check anyone can make on a tree
    assert s["realised_rates"]["birth"] == pytest.approx(
        s["events"]["speciation"] / s["tree"]["total_branch_length"], abs=1e-6)


def test_the_sequences_summary_reports_what_the_run_printed(tmp_path):
    sp, g = _run()
    r = simulate_sequences(g, model=hky85(2.0), length=200, divergence=0.2, seed=1)
    r.write(tmp_path, outputs=("alignments", "summary"))
    s = json.loads((tmp_path / "sequences_summary.json").read_text(encoding="utf-8"))
    assert s["families_with_sequences"] == len(list((tmp_path / "alignments").glob("*.fasta")))
    assert s["sequences"] == sum(len(a) for a in r.alignments.values())
    assert s["sites"] == {"min": 200, "max": 200}
    assert 0.0 < s["mean_pairwise_identity"] < 1.0


def test_a_summary_is_reproducible_byte_for_byte(tmp_path):
    # it is an output like any other, so the same seed must give the same file
    for name in ("a", "b"):
        sp, g = _run()
        g.write(tmp_path / name, outputs=("summary",))
    assert (tmp_path / "a" / "genome_summary.json").read_bytes() == \
           (tmp_path / "b" / "genome_summary.json").read_bytes()


def test_the_cli_writes_a_summary_at_every_level(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "10",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run), "--duplication", "0.3", "--transfer", "0.1", "--loss", "0.2",
                 "--origination", "0.4", "--initial-families", "10", "--seed", "1", "--quiet"]) == 0
    assert main(["sequences", str(run), "--model", "jc69", "--length", "100", "--divergence", "0.2",
                 "--seed", "1", "--quiet"]) == 0
    for level, name in (("species", "species_summary.json"), ("genomes", "genome_summary.json"),
                        ("sequences", "sequences_summary.json")):
        path = run / level / name
        assert path.exists(), f"no summary for {level}"
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["level"] == level and payload["seed"] == 1


def test_the_summary_can_be_left_out(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--n-extant", "8", "--seed", "1",
                 "--write", "complete", "extant", "--quiet"]) == 0
    assert not (run / "species" / "species_summary.json").exists()
