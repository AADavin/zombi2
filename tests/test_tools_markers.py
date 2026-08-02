"""Tests for :mod:`zombi2.tools.markers` — can a family be trusted as a phylogenetic marker?

The question is about a family, not a pair: single-copy, universal, and does its true tree match the
species tree over the genomes it occupies. The last is the one that matters, because it is invisible
in real data — a family can pass every filter and still give the wrong topology.
"""

import collections
import csv

import pytest

from zombi2.genomes import simulate_genomes_family
from zombi2.species import simulate_species_tree
from zombi2.tools.markers import marker_row, markers_tsv


def _rows(text):
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def test_with_no_duplication_or_transfer_every_family_is_a_perfect_marker():
    # the control: nothing can make a gene tree disagree with the species tree, so the congruence
    # column must say so for every family — otherwise it is measuring something else
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=12, seed=3)
    g = simulate_genomes_family(sp, initial_families=40, loss=0.05, seed=5)
    rows = _rows(markers_tsv(g.gene_trees, sp.complete_tree))
    assert rows and all(r["single_copy"] == "yes" for r in rows)
    assert all(r["congruent"] == "yes" and r["rf"] == "0" for r in rows)


def test_a_family_can_pass_every_filter_and_still_give_the_wrong_tree():
    # the whole point. Replacing transfers keep a family single-copy and universal — it survives every
    # filter a phylogenomicist applies — while moving its topology away from the species tree.
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=15, seed=3)
    g = simulate_genomes_family(sp, duplication=0.02, transfer=0.25, loss=0.05,
                                replacement=True, initial_families=300, seed=5)
    rows = _rows(markers_tsv(g.gene_trees, sp.complete_tree))
    clean = [r for r in rows if r["single_copy"] == "yes" and r["universal"] == "yes"]
    assert len(clean) > 20, "the run should leave plenty of families that look usable"
    misleading = [r for r in clean if r["congruent"] == "no"]
    assert misleading, "replacing transfer should have moved some of them off the species tree"
    assert all(int(r["transfers"]) > 0 for r in misleading)     # and it is the transfers that did it


def test_rf_is_left_empty_where_it_would_mean_nothing():
    # several copies in one genome: no one-to-one gene -> genome map, so no distance to a species tree
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=10, seed=2)
    g = simulate_genomes_family(sp, duplication=0.8, initial_families=25, seed=4)
    rows = _rows(markers_tsv(g.gene_trees, sp.complete_tree))
    multi = [r for r in rows if r["single_copy"] == "no"]
    assert multi and all(r["rf"] == "" and r["congruent"] == "" for r in multi)


def test_universal_means_present_in_every_extant_genome():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1)
    g = simulate_genomes_family(sp, loss=0.4, origination=0.5, initial_families=30, seed=2)
    n_extant = len(sp.complete_tree.extant_leaves())
    for family, gt in g.gene_trees.items():
        if gt.extant is None:
            continue
        row = marker_row(family, gt, sp.complete_tree)
        assert row["universal"] == (row["genomes"] == n_extant)
        assert row["genomes"] <= n_extant


def test_the_counts_are_the_family_s_own_history():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=7)
    g = simulate_genomes_family(sp, duplication=0.3, transfer=0.2, loss=0.3, initial_families=20,
                                seed=3)
    # every duplication / transfer / loss in the log belongs to exactly one family, so the table's
    # columns must add up to the log
    from_log = collections.Counter()
    for e in g.events:
        if e.kind in ("duplication", "transfer_additive", "transfer_replacing", "loss"):
            from_log[(e.family, "transfer" if e.kind.startswith("transfer") else e.kind)] += 1
    for family, gt in g.gene_trees.items():
        row = marker_row(family, gt, sp.complete_tree)
        assert row["duplications"] == from_log[(family, "duplication")]
        assert row["transfers"] == from_log[(family, "transfer")]
        assert row["losses"] == from_log[(family, "loss")]


def test_a_family_in_part_of_the_tree_is_judged_against_the_part_it_occupies():
    # a family present in five genomes must not be scored against the whole species tree; the
    # comparison is the species tree restricted to those five
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=14, seed=5)
    g = simulate_genomes_family(sp, initial_families=60, loss=0.25, seed=6)
    rows = _rows(markers_tsv(g.gene_trees, sp.complete_tree))
    partial = [r for r in rows if r["universal"] == "no" and r["rf"]]
    assert partial, "the loss rate should leave families in only part of the tree"
    # no duplication and no transfer, so every one of them still recovers the species tree exactly
    assert all(r["congruent"] == "yes" for r in partial)


def test_the_cli_writes_one_table_for_the_whole_run(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--n-extant", "8", "--seed", "1",
                 "--quiet"]) == 0
    assert main(["genomes", str(run), "--duplication", "0.2", "--transfer", "0.1",
                 "--initial-families", "10", "--seed", "1", "--quiet"]) == 0
    assert main(["tools", "format", str(run), "--format", "markers", "--quiet"]) == 0
    out = run / "genomes" / "markers.tsv"
    assert out.exists()                                    # one file, not one per family
    assert not (run / "genomes" / "markers").exists()      # and so no directory of its own
    rows = _rows(out.read_text(encoding="utf-8"))
    assert rows and set(rows[0]) == {"family", "genomes", "copies", "single_copy", "universal",
                                     "duplications", "transfers", "losses", "rf", "congruent"}


@pytest.mark.parametrize("seed", range(1, 4))
def test_congruence_agrees_with_treedist(seed, tmp_path, capsys):
    # the number must mean what `zombi2 tools treedist` reports, or two parts of one tool disagree.
    # Run the command itself: it already knows how to compare a single-copy gene tree to a species
    # tree, on the species each gene sits in, which is exactly this measurement.
    from zombi2.cli.main import main

    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=9, seed=seed)
    g = simulate_genomes_family(sp, transfer=0.25, replacement=True, initial_families=60,
                                seed=seed * 3)
    species_file = tmp_path / "sp.nwk"
    species_file.write_text(sp.extant_tree.to_newick() + "\n", encoding="utf-8")
    n_extant = len(sp.complete_tree.extant_leaves())
    checked, disagreeing = 0, 0
    for family, gt in g.gene_trees.items():
        if gt.extant is None:
            continue
        row = marker_row(family, gt, sp.complete_tree)
        if not row["single_copy"] or row["genomes"] != n_extant:
            continue                              # treedist needs the same genomes on both trees
        gene_file = tmp_path / f"g{family}.nwk"
        gene_file.write_text(gt.to_newick("extant") + "\n", encoding="utf-8")
        capsys.readouterr()
        assert main(["tools", "treedist", str(gene_file), str(species_file), "--metric", "rf"]) == 0
        reported = float(capsys.readouterr().out.strip().split("\t")[1])
        assert row["rf"] == reported, family
        checked += 1
        disagreeing += row["congruent"] is False
    assert checked > 3 and disagreeing > 0        # both agreements and disagreements were seen
