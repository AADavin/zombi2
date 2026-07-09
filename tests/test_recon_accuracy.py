"""Tests for reconciliation accuracy (``zombi2.tools.recon_accuracy``)."""
from __future__ import annotations

import pytest

from zombi2 import Yule, simulate_genomes, simulate_species_tree
from zombi2.genomes.reconciliation import Reconciliation
from zombi2.tools.recon_accuracy import reconciliation_accuracy

# A pure-speciation truth over four tips: three internal nodes (root, sX, sY), all S.
TRUTH = "((A|1,B|2)sX|S,(C|3,D|4)sY|S)root|S;"


# --------------------------------------------------------------------------- perfect / identity

def test_identical_reconciliation_is_perfect():
    a = reconciliation_accuracy(TRUTH, TRUTH)
    assert a.n_nodes == 3
    assert a.event_accuracy == 1.0 and a.mapping_accuracy == 1.0 and a.joint_accuracy == 1.0
    assert a.per_event["S"].precision == 1.0 and a.per_event["S"].recall == 1.0
    assert a.transfer.n_true == 0


def test_child_order_does_not_matter():
    # Same tree with the two root subtrees swapped and cherries flipped — identical topology.
    reordered = "((D|4,C|3)sY|S,(B|2,A|1)sX|S)root|S;"
    a = reconciliation_accuracy(TRUTH, reordered)
    assert a.event_accuracy == 1.0 and a.mapping_accuracy == 1.0


def test_accepts_reconciliation_object():
    recon = Reconciliation(complete=None, extant=TRUTH, events=[])
    a = reconciliation_accuracy(recon, recon)
    assert a.n_nodes == 3 and a.joint_accuracy == 1.0


# --------------------------------------------------------------------------- event errors

def test_one_wrong_event_type():
    # sX inferred as a duplication instead of a speciation.
    inferred = "((A|1,B|2)sX|D,(C|3,D|4)sY|S)root|S;"
    a = reconciliation_accuracy(TRUTH, inferred)
    assert a.event_accuracy == pytest.approx(2 / 3)
    assert a.mapping_accuracy == 1.0                 # species branches still all correct
    assert a.joint_accuracy == pytest.approx(2 / 3)
    s = a.per_event["S"]
    assert (s.tp, s.support_true, s.support_pred) == (2, 3, 2)
    assert s.recall == pytest.approx(2 / 3) and s.precision == 1.0
    d = a.per_event["D"]
    assert (d.tp, d.support_true, d.support_pred) == (0, 0, 1)
    assert d.precision == 0.0 and d.recall == 0.0


def test_wrong_species_mapping_only():
    # right event (S) but sX mapped to the wrong species branch.
    inferred = "((A|1,B|2)sZ|S,(C|3,D|4)sY|S)root|S;"
    a = reconciliation_accuracy(TRUTH, inferred)
    assert a.event_accuracy == 1.0
    assert a.mapping_accuracy == pytest.approx(2 / 3)
    assert a.joint_accuracy == pytest.approx(2 / 3)


# --------------------------------------------------------------------------- transfers

TRUTH_T = "((A|1,B|2)sX|S,(C|3,D|4)donor|T>recip)root|S;"


def test_transfer_fully_recovered():
    a = reconciliation_accuracy(TRUTH_T, TRUTH_T)
    assert a.transfer.n_true == 1
    assert a.transfer == (1, 1, 1, 1, 1)             # detected, donor, recipient, both all correct


def test_transfer_missed():
    inferred = "((A|1,B|2)sX|S,(C|3,D|4)donor|S)root|S;"   # called a speciation, not a transfer
    a = reconciliation_accuracy(TRUTH_T, inferred)
    assert a.transfer.n_true == 1 and a.transfer.detected == 0
    assert a.transfer.both_correct == 0
    assert a.per_event["T"].recall == 0.0


def test_transfer_wrong_recipient():
    inferred = "((A|1,B|2)sX|S,(C|3,D|4)donor|T>WRONG)root|S;"
    a = reconciliation_accuracy(TRUTH_T, inferred)
    assert a.transfer == (1, 1, 1, 0, 0)             # detected + donor right, recipient wrong
    assert a.event_accuracy == 1.0                   # event type (T) is still correct


def test_transfer_wrong_donor():
    inferred = "((A|1,B|2)sX|S,(C|3,D|4)OTHER|T>recip)root|S;"
    a = reconciliation_accuracy(TRUTH_T, inferred)
    assert a.transfer == (1, 1, 0, 1, 0)             # detected + recipient right, donor wrong
    assert a.mapping_accuracy == pytest.approx(2 / 3)  # the transfer node's species branch is wrong


# --------------------------------------------------------------------------- errors & edge cases

def test_topology_mismatch_raises():
    other = "((A|1,C|3)x|S,(B|2,D|4)y|S)root|S;"     # regroups the tips
    with pytest.raises(ValueError, match="topology"):
        reconciliation_accuracy(TRUTH, other)


def test_tip_label_mismatch_raises():
    other = "((A|1,B|2)sX|S,(C|3,E|9)sY|S)root|S;"   # D|4 -> E|9
    with pytest.raises(ValueError, match="same gene tree"):
        reconciliation_accuracy(TRUTH, other)


def test_single_leaf_tree_has_no_nodes():
    a = reconciliation_accuracy("A|1;", "A|1;")
    assert a.n_nodes == 0
    assert a.event_accuracy == 0.0 and a.transfer.n_true == 0


# --------------------------------------------------------------------------- real simulated data

def _sim_reconciliations():
    tree = simulate_species_tree(Yule(1.0), n_tips=8, age=2.0, seed=11)
    g = simulate_genomes(tree, duplication=0.2, transfer=0.15, loss=0.2,
                         origination=0.0, initial_families=25, seed=11)
    return {f: r for f, r in g.reconciliations().items() if r.extant is not None}


def test_real_reconciliations_perfect_against_themselves():
    recons = _sim_reconciliations()
    assert len(recons) > 0
    for r in recons.values():
        a = reconciliation_accuracy(r, r)
        if a.n_nodes == 0:
            continue
        assert a.event_accuracy == 1.0 and a.mapping_accuracy == 1.0 and a.joint_accuracy == 1.0
        assert a.transfer.n_true == a.transfer.both_correct


def test_corrupting_an_event_lowers_accuracy():
    # flip one speciation label to a duplication in a real reconciled tree.
    recons = _sim_reconciliations()
    fam = next(r for r in recons.values() if "|S:" in r.extant and r.extant.count("|S") >= 1)
    corrupted = fam.extant.replace("|S:", "|D:", 1)
    assert corrupted != fam.extant
    a = reconciliation_accuracy(fam, corrupted)
    assert a.event_accuracy < 1.0
