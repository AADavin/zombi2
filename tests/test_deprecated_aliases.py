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


def test_singular_command_aliases_still_work_and_warn(tmp_path, capsys):
    """The commands are plural (`traits`, `sequences`); the singular spellings are accepted but
    deprecated — they warn and produce the canonical (plural) run-manifest name."""
    from zombi2.cli import main
    tree = tmp_path / "sp.nwk"
    tree.write_text("((A:1,B:1):1,C:2);")
    rc = main(["trait", "-t", str(tree), "--model", "bm", "--seed", "1", "-o", str(tmp_path / "o")])
    assert rc == 0
    err = capsys.readouterr().err
    assert "deprecated" in err and "traits" in err
    assert (tmp_path / "o" / "traits.log").exists()   # canonical output, whichever spelling was used


# C8: the traits:genomes edge's TraitLinked* names were unified onto the TraitGene* stem.
COEVOLVE_DEPRECATED = {
    "TraitLinkedRates": "TraitGeneRates",
    "TraitLinkedResult": "TraitGeneResult",
    "simulate_trait_linked_genomes": "simulate_trait_conditioned_genomes",
}


@pytest.mark.parametrize("old,new", list(COEVOLVE_DEPRECATED.items()))
def test_coevolve_traitgene_alias_warns_and_resolves(old, new):
    """The old TraitLinked* names warn and resolve to the TraitGene* stem, at both the top level
    and the coevolve namespace; the deep implementation module stays a silent anchor."""
    import zombi2.coevolve as cv
    with pytest.warns(DeprecationWarning, match=old):
        assert getattr(z, old) is getattr(z, new)
    with pytest.warns(DeprecationWarning, match=old):
        assert getattr(cv, old) is getattr(cv, new)
    assert old not in z.__all__ and old not in cv.__all__ and old not in dir(cv)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # the deep module must not warn
        from zombi2.coevolve import trait_coupling as tc
        assert getattr(tc, old) is getattr(tc, new)
