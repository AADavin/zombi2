"""Tests for the run summary — what came out, beside the log that says what went in.

Three things nobody could get from a run's output, each of which a reader reconstructed by hand and
one of them wrongly:

- how many events there were. The log now writes one row per event, so counting rows answers it —
  and the summary must agree with those rows, which is what the first test checks;
- ``gene_trees/`` holds a pair of files per family that *ever existed* while the summary line counts
  the families that **survived**, and neither number explained the other;
- the family-size cap was invisible: when it binds it discards events, so realised rates fall below
  the declared ones, and nothing said so.

The first is the one that matters most, because it produces a plausible wrong number rather than an
error, so it is checked against the log itself rather than against the summary's own logic.
"""

import collections
import csv
import json

import pytest

from zombi2.params import PerLineage
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85
from zombi2.species import simulate_species_tree


def _run(**kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1)
    return sp, simulate_genomes_family(sp, **{"duplication": 0.4, "transfer": 0.2, "loss": 0.2,
                                              "origination": 0.4, "initial_families": 15,
                                              "seed": 1, **kw})


@pytest.mark.parametrize("replacement", [False, True])
def test_the_event_counts_agree_with_the_logs_rows(tmp_path, replacement):
    # the whole point. The log is one row per event, so the summary's counts and a `cut -f2 | sort |
    # uniq -c` on the file are the same numbers — checked against the file, not against the summary's
    # own logic, because a wrong count here is a plausible number rather than a crash.
    sp, g = _run(replacement=replacement)
    g.write(tmp_path, outputs=("events", "summary"))
    rows = list(csv.DictReader((tmp_path / "genome_events.tsv").read_text(
        encoding="utf-8").splitlines(), delimiter="\t"))
    from_log = collections.Counter(r["kind"] for r in rows)

    s = json.loads((tmp_path / "genome_summary.json").read_text(encoding="utf-8"))
    assert s["events"]["duplication"] == from_log["duplication"]
    assert s["events"]["speciation"] == from_log["speciation"]
    # a transfer is written as one kind or the other by whether the arriving copy replaced a resident
    assert s["events"]["transfer"] == from_log["transfer_additive"] + from_log["transfer_replacing"]
    # ...and a replacing transfer kills the copy it overwrote without spending a `loss` row on it:
    # that death is the row's second parent, and the summary counts it, because the gene tree has it
    assert s["events"]["loss"] == from_log["loss"] + from_log["transfer_replacing"]
    assert bool(from_log["transfer_replacing"]) is replacement
    # origination is split: the initial genome is logged as origination at the root's own start time,
    # so a bare count of that kind is de-novo arrivals PLUS initial_families
    assert s["events"]["initial"] + s["events"]["origination"] == from_log["origination"]
    assert s["events"]["initial"] == 15                      # what --initial-families asked for
    assert len(rows) == sum(s["events"].values()) - from_log["transfer_replacing"]


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
                   for n in tight.complete_tree.extant_leaves()) == 4

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


def test_a_discrete_trait_summary_describes_its_states():
    # what you look at first on a discrete run: a trait whose tips all share one state has told you
    # nothing, and `most_common_share` is the number that says so
    from zombi2 import traits

    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=1)
    r = traits.simulate_discrete(sp, states=("cave", "surface"), switch=0.5, seed=1)
    s = r.summary()
    assert s["kind"] == "discrete" and s["level"] == "traits"
    assert sum(s["states"].values()) == s["tips"] == len(r.values)
    assert s["states_at_tips"] == len(s["states"])
    assert s["most_common_share"] == pytest.approx(max(s["states"].values()) / s["tips"])
    assert s["events"]["on_branch"] == sum(1 for e in r.events if e.kind == "on_branch")
    assert "values" not in s                      # a state is not a number to take a mean of


def test_a_continuous_trait_summary_describes_where_the_values_got_to():
    from zombi2 import traits

    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=1)
    r = traits.simulate_continuous(sp, rate=2.0, seed=1)
    s = r.summary()
    assert s["kind"] == "continuous"
    vals = list(r.values.values())
    assert s["values"]["min"] == min(vals) and s["values"]["max"] == max(vals)
    # the ROOT NODE, i.e. after the trait diffused along the stem — not the value the run started
    # from, which belongs to no node at all
    assert s["value_at_root_node"] == r.node_values[sp.complete_tree.root]
    assert s["value_at_root_node"] != 0.0         # it diffused; `start` was 0
    assert "states" not in s


def test_a_joint_summary_holds_both_levels_it_grew():
    from zombi2.joint import simulate_joint
    from zombi2.traits import DiscreteTrait

    r = simulate_joint(birth=PerLineage(1.0).scaled_by("trait", {"small": 1.0, "large": 3.0}), death=0.2,
                       trait=DiscreteTrait(states=("small", "large"), switch=0.3),
                       n_extant=20, seed=1)
    s = r.summary()
    assert s["level"] == "joint" and s["driver"] == "trait"
    # the same payloads the two levels would write alone — one vocabulary, not a third
    assert s["species"] == r.species.summary()
    assert s["trait"] == r.trait.summary()
    assert "genome" not in s                      # exactly one driver per run
    # the tree is an OUTPUT here, so its realised birth rate is what the driver did: a trait whose
    # tips are mostly in the 3x state must have pushed it above the declared base of 1.0
    assert s["species"]["realised_rates"]["birth"] > 1.0


def test_the_cli_writes_a_joint_summary_at_the_run_root(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "j"
    assert main(["joint", str(run), "--birth",
                 "PerLineage(1.0).scaled_by('trait', {'small': 1.0, 'large': 3.0})",
                 "--death", "0.2",
                 "--states", "small,large", "--switch", "0.3", "--n-extant", "20", "--seed", "1",
                 "--quiet"]) == 0
    # at the root, because it describes both levels and belongs to neither
    payload = json.loads((run / "joint_summary.json").read_text(encoding="utf-8"))
    assert payload["level"] == "joint"
    assert payload["species"]["tips"]["extant"] == 20
    # and each level still has its own, so a joint run is not a special case downstream
    assert (run / "species" / "species_summary.json").exists()
    assert (run / "traits" / "trait_summary.json").exists()


def test_a_threshold_trait_is_summarised_as_states_not_numbers():
    # a threshold trait reads a discrete state off a continuous liability, so its values are states.
    # The first version of this took a mean of them and raised, which took the whole CLI run down.
    from zombi2 import traits

    sp = simulate_species_tree(birth=1.0, n_extant=10, seed=1)
    r = traits.simulate_discrete(sp, states=("low", "high"), liability=1.0, threshold=0.0, seed=1)
    assert r.kind == "threshold"
    s = r.summary()
    assert s["kind"] == "threshold"
    assert sum(s["states"].values()) == s["tips"]
    assert "values" not in s and "value_at_root_node" not in s
