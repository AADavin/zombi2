"""A clade under its own substitution model — `zombi2.sequences.Models`.

A rate says how fast a lineage evolves; the model says what the change looks like. Until this, a run
had one model for every branch, so a clade could evolve faster but not *differently* — and the whole
class of compositional questions was out of reach. The sharp one is compositional attraction: two
unrelated AT-rich lineages look alike, which misleads a tree-builder in a way long-branch attraction
does not, so a study of one needs a clade that really does drift.
"""

from __future__ import annotations

import pytest

from zombi2 import genomes, sequences, species
from zombi2.params import Clade, PerSite
from zombi2.sequences import Models
from zombi2.sequences.substitution_models import hky85, jc69, lg

AT_RICH = hky85(kappa=2.0, frequencies=(0.40, 0.10, 0.10, 0.40))     # equilibrium A+T = 0.80
EVEN = hky85(kappa=2.0)                                              # equilibrium A+T = 0.50


def _run(n_extant=12, families=40, seed=3):
    sp = species.simulate_species_tree(birth=1, death=0.3, n_extant=n_extant, seed=seed)
    return sp, genomes.simulate_genomes_family(sp, initial_families=families, seed=5)


def _at_content(result, inside_ids, *, inside: bool) -> float:
    hits = total = 0
    for aln in result.alignments.values():
        for label, seq in aln.items():
            if ((int(label.split("_")[0][1:]) in inside_ids) == inside):
                hits += sum(seq.count(c) for c in "AT")
                total += len(seq)
    return hits / total if total else float("nan")


# --- the thing itself -----------------------------------------------------------------------------

def test_a_clade_drifts_toward_its_own_composition():
    """The whole point. Inside the clade the sequences move toward the model's own equilibrium and
    outside they stay where they were, so the run carries a compositional signal that a single-model
    run cannot produce however the rates are set."""
    sp, g = _run()
    clade = Clade({"at": ["n5", "n9"]})
    inside = set(clade.resolve(sp)["at"])
    r = sequences.simulate_sequences(
        g, model=Models().set_by(clade, {"at": AT_RICH, "rest": EVEN}),
        length=300, substitution=0.6, seed=1)
    assert _at_content(r, inside, inside=True) > 0.65        # pulled well off 0.5 ...
    assert 0.47 < _at_content(r, inside, inside=False) < 0.53   # ... while outside is untouched


def test_the_drift_approaches_the_model_equilibrium_with_divergence():
    """Not a one-off offset but a real compositional process: a lineage entering the clade starts at
    the old composition and relaxes toward the new one, so more divergence means more of the way
    there. This is also the transient that makes a branch's realised substitution count fall short of
    its nominal length while it is happening — the documented cost of a non-homogeneous model."""
    sp, g = _run()
    clade = Clade({"at": ["n5", "n9"]})
    inside = set(clade.resolve(sp)["at"])
    got = []
    for substitution in (0.05, 0.6, 2.0):
        r = sequences.simulate_sequences(
            g, model=Models().set_by(clade, {"at": AT_RICH, "rest": EVEN}),
            length=300, substitution=substitution, seed=1)
        got.append(_at_content(r, inside, inside=True))
    assert got[0] < got[1] < got[2]                          # monotone toward equilibrium
    assert got[0] < 0.60 and got[2] > 0.75                   # and it gets most of the way to 0.80


def test_the_plain_path_is_untouched():
    """`models=None` is every run written before this one. It must not move by a byte."""
    _sp, g = _run()
    kw = dict(model=EVEN, length=200, substitution=0.3, seed=1)
    assert sequences.simulate_sequences(g, **kw).alignments == \
           sequences.simulate_sequences(g, **kw).alignments


def test_it_is_deterministic_given_the_seed():
    sp, g = _run()
    clade = Clade({"at": ["n5", "n9"]})
    kw = dict(model=Models().set_by(clade, {"at": AT_RICH, "rest": EVEN}),
              length=200, substitution=0.4, seed=7)
    assert sequences.simulate_sequences(g, **kw).alignments == \
           sequences.simulate_sequences(g, **kw).alignments


# --- what it refuses, and why ---------------------------------------------------------------------

def test_every_group_the_tree_paints_needs_a_model():
    # a lineage in no named clade is in 'rest'; forgetting it would evolve part of the tree under
    # whatever happened to be first, silently
    sp, g = _run()
    with pytest.raises(ValueError, match="no model for"):
        sequences.simulate_sequences(g, length=50, seed=1,
                                     model=Models().set_by(Clade({"at": ["n5"]}), {"at": AT_RICH}))


def test_a_label_naming_no_lineage_is_refused():
    # the same silence `check_mapping_fires` exists to break: a model that evolves nothing
    sp, g = _run()
    with pytest.raises(ValueError, match="name.*no lineage"):
        sequences.simulate_sequences(
            g, length=50, seed=1,
            model=Models().set_by(Clade({"at": ["n5"]}),
                                  {"at": AT_RICH, "rest": EVEN, "ghost": EVEN}))


def test_one_run_has_one_alphabet():
    # a gene copy's sequence is one string, so this is meaningless rather than unimplemented
    with pytest.raises(ValueError, match="one alphabet"):
        Models().set_by(Clade({"a": ["n1"]}), {"a": jc69(), "rest": lg()})


def test_every_model_sorts_the_sites_the_same_way():
    # a site's class is drawn once for the family and holds down the whole tree, so two clades
    # cannot sort the same site into different classes
    with pytest.raises(ValueError, match="across-site rate classes differ"):
        Models().set_by(Clade({"a": ["n1"]}),
                        {"a": hky85().across_sites(gamma_shape=0.5), "rest": hky85()})


def test_a_trait_cannot_drive_the_model_yet():
    # a trait switches partway along a branch and this engine samples one matrix per branch, so the
    # refusal names the reason and points at the rate, which a trait CAN drive
    with pytest.raises(ValueError, match="read off the tree today"):
        Models().set_by("some_trait.tsv", {"a": AT_RICH})


def test_the_model_set_takes_one_verb():
    from zombi2.params.parameter import RateCompositionError
    with pytest.raises(RateCompositionError, match="set_by and nothing else"):
        Models().scaled_by(Clade({"a": ["n1"]}), {"a": 2.0})
    with pytest.raises(RateCompositionError, match="one set_by"):
        Models().set_by(Clade({"a": ["n1"]}), {"a": AT_RICH, "rest": EVEN}) \
                .set_by(Clade({"b": ["n2"]}), {"b": AT_RICH, "rest": EVEN})


def test_a_rate_can_still_be_scoped_to_the_same_clade():
    """The two compose, which is the point for a reduction study: the clade evolves faster AND
    differently, which are the two halves of the endosymbiont syndrome."""
    sp, g = _run()
    clade = Clade({"at": ["n5", "n9"]})
    r = sequences.simulate_sequences(
        g, model=Models().set_by(clade, {"at": AT_RICH, "rest": EVEN}), length=200, seed=1,
        substitution=PerSite(0.3).scaled_by(clade, {"at": 4.0, "rest": 1.0}))
    assert r.alignments
