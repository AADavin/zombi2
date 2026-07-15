"""The genome opportunity knob: ``Rates(per="copy"|"lineage")`` and ``FamilySampledRates(per=…)``.

``per`` selects the unit the D/T/L clock rides on. ``per="copy"`` (default) scales the total rate by
copy number → exponential families (the built-in Rust model); ``per="lineage"`` is a constant rate
per genome → linear families. It is orthogonal to per-family heterogeneity. Also covers the deprecated
``PerCopyRates`` / ``PerLineageRates`` presets and Rust eligibility.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

import zombi2 as z
from zombi2 import _rust
from zombi2.genomes.genome import UnorderedGenome


def _tree(seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.2), age=4.0, direction="forward", seed=seed)


def _max_copies(rates, seed=3):
    g = z.simulate_genomes(_tree(), rates, initial_families=1, seed=seed)
    return max(gm.size() for gm in g.leaf_genomes.values())


def test_per_copy_is_exponential_per_lineage_is_linear():
    # same duplication rate: per-copy explodes; per-lineage stays small (linear)
    assert _max_copies(z.Rates(duplication=0.6, per="copy")) > 3 * _max_copies(z.Rates(duplication=0.6, per="lineage"))
    assert z.Rates(duplication=0.3).per == "copy"   # default


def test_per_copy_uses_rust_per_lineage_does_not():
    assert _rust.eligible(z.Rates(duplication=0.2, transfer=0.1, loss=0.2), UnorderedGenome, None)
    assert not _rust.eligible(z.Rates(duplication=0.2, per="lineage"), UnorderedGenome, None)


def test_deprecated_presets_are_byte_identical_to_the_knob():
    tree = _tree()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # per="copy" == PerCopyRates (Rust path)
        a = z.simulate_genomes(tree, z.Rates(duplication=0.5, transfer=0.1, loss=0.2), seed=7).profiles.to_tsv()
        b = z.simulate_genomes(tree, z.PerCopyRates(duplication=0.5, transfer=0.1, loss=0.2), seed=7).profiles.to_tsv()
        # per="lineage" == PerLineageRates (Python path)
        c = z.simulate_genomes(tree, z.Rates(duplication=0.5, transfer=0.1, loss=0.2, per="lineage"), seed=7).profiles.to_tsv()
        d = z.simulate_genomes(tree, z.PerLineageRates(duplication=0.5, transfer=0.1, loss=0.2), seed=7).profiles.to_tsv()
    assert a == b
    assert c == d


def test_presets_warn_and_left_public_api():
    for preset in (z.PerCopyRates, z.PerLineageRates):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            preset(duplication=0.1)
        assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert "Rates" in z.__all__
    assert "PerCopyRates" not in z.__all__ and "PerLineageRates" not in z.__all__
    assert issubclass(z.PerCopyRates, z.Rates) and issubclass(z.PerLineageRates, z.Rates)


def test_validation():
    with pytest.raises(ValueError, match="per must be"):
        z.Rates(duplication=0.3, per="bogus")
    with pytest.raises(ValueError, match="supports only"):
        z.Rates(duplication=0.3, per="lineage", inversion=0.1)   # rearrangements are per-copy
    with pytest.raises(ValueError, match="transfer"):
        z.Rates(duplication=0.3, transfer=0.1, per="shared")     # transfer not yet supported for shared


def test_per_shared_is_one_tree_wide_clock():
    # per="shared": a single tree-wide duplication clock → #dup events ≈ base·age, INDEPENDENT of how
    # many lineages the family spans (contrast per-lineage, which scales with the lineage count).
    from zombi2.genomes.events import EventType
    tree = z.simulate_species_tree(z.BirthDeath(1.5, 0.2), age=5.0, direction="forward", seed=1)
    assert len(tree.extant_leaves()) > 20   # a bushy tree, so "independent of lineage count" bites

    def n_dups(rates, seed):
        g = z.simulate_genomes(tree, rates, initial_families=1, seed=seed)
        return sum(1 for e in g.event_log.records if e.event == EventType.DUPLICATION)

    shared = [n_dups(z.Rates(duplication=0.6, per="shared"), s) for s in range(120)]
    lineage = [n_dups(z.Rates(duplication=0.6, per="lineage"), s) for s in range(40)]
    assert 2.0 < np.mean(shared) < 4.5          # ≈ 0.6 × 5 = 3, not scaled by the many lineages
    assert np.mean(lineage) > 3 * np.mean(shared)  # per-lineage scales with the lineage count → far more


def test_modifiers_reject_a_shared_base():
    # a modifier scales per-branch weights, but per="shared" routes dup/loss through a tree-wide pool
    # a ModifiedRates never sees — wrapping one would silently switch the shared clock off, so refuse it
    with pytest.raises(ValueError, match="shared"):
        z.LineageRates(z.Rates(duplication=0.5, per="shared"), factors={"n1": 2.0})
    with pytest.raises(ValueError, match="shared"):
        z.ModifiedRates(z.Rates(duplication=0.5, per="shared"), [z.FamilyModifier(factors={"1": 2.0})])


def test_per_shared_is_deterministic_and_python_engine():
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), age=4.0, direction="forward", seed=1)
    a = z.simulate_genomes(tree, z.Rates(duplication=0.5, loss=0.1, per="shared"), initial_families=2, seed=9)
    b = z.simulate_genomes(tree, z.Rates(duplication=0.5, loss=0.1, per="shared"), initial_families=2, seed=9)
    assert a.profiles.to_tsv() == b.profiles.to_tsv()          # reproducible
    assert not _rust.eligible(z.Rates(duplication=0.5, per="shared"), UnorderedGenome, None)  # Python


def test_family_sampled_rates_per_lineage_closes_the_wrinkle():
    # per-family rates AND a per-lineage opportunity — the combination that needed the modifier route
    tree = _tree()
    fam_copy = z.FamilySampledRates(rates={"1": (0.8, 0.0, 0.1)}, per="copy")
    fam_lin = z.FamilySampledRates(rates={"1": (0.8, 0.0, 0.1)}, per="lineage")
    # copy scales with copies → bigger families than the fixed-per-genome lineage clock
    big = max(gm.size() for gm in z.simulate_genomes(tree, fam_copy, initial_families=1, seed=5).leaf_genomes.values())
    small = max(gm.size() for gm in z.simulate_genomes(tree, fam_lin, initial_families=1, seed=5).leaf_genomes.values())
    assert big >= small
    with pytest.raises(ValueError):
        z.FamilySampledRates(duplication=0.3, per="lineage", carrying_capacity=10)
