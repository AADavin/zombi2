"""The deprecation window for renamed public names (see docs/design/naming-consolidation.md).

Five rate names were renamed; the old spellings keep working for one minor version but are
*marked*: they resolve through a PEP-562 ``__getattr__`` with a ``DeprecationWarning`` and are
absent from ``__all__``/``dir()`` (so they leave the API reference and tab-completion). The silent
deep-module anchor (``zombi2.genomes.rates``) is unchanged so existing deep imports do not warn.
"""

from __future__ import annotations

import warnings

import pytest

import zombi2 as z
import zombi2.genomes as zg

# old -> new canonical name
DEPRECATED = {
    "SharedRates": "PerCopyRates",
    "PerGenomeRates": "PerLineageRates",
    "BranchRates": "LineageRates",
    "BranchModifier": "LineageModifier",
    "read_branch_rates": "read_lineage_rates",
}


@pytest.mark.parametrize("old,new", list(DEPRECATED.items()))
def test_top_level_alias_warns_and_resolves(old, new):
    """``zombi2.<old>`` warns once and returns the *same object* as ``zombi2.<new>``."""
    with pytest.warns(DeprecationWarning, match=old):
        obj = getattr(z, old)
    assert obj is getattr(z, new)


@pytest.mark.parametrize("old,new", list(DEPRECATED.items()))
def test_genomes_namespace_alias_warns_and_resolves(old, new):
    """``zombi2.genomes.<old>`` warns and resolves to the canonical namespace object."""
    with pytest.warns(DeprecationWarning, match=old):
        obj = getattr(zg, old)
    assert obj is getattr(zg, new)


@pytest.mark.parametrize("old", list(DEPRECATED))
def test_aliases_are_absent_from_the_public_catalog(old):
    """The old names must not appear in ``__all__`` or ``dir()`` (gone from docs + completion)."""
    assert old not in z.__all__
    assert old not in zg.__all__
    assert old not in dir(z)
    assert old not in dir(zg)


def test_unknown_attribute_still_raises_attributeerror():
    """``__getattr__`` only intercepts the known aliases; everything else is a normal error."""
    with pytest.raises(AttributeError):
        z.NoSuchThing
    with pytest.raises(AttributeError):
        zg.NoSuchThing


def test_deep_module_anchor_stays_silent():
    """Deep imports from the implementation modules keep working with NO warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning here fails the test
        from zombi2.genomes.rates import (
            SharedRates, PerGenomeRates, BranchRates, BranchModifier,
        )
        from zombi2.genomes.read_rates import read_branch_rates
    assert SharedRates is z.PerCopyRates
    assert PerGenomeRates is z.PerLineageRates
    assert BranchRates is z.LineageRates
    assert BranchModifier is z.LineageModifier
    assert read_branch_rates is z.read_lineage_rates
