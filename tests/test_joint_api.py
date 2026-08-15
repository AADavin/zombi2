"""`joint.simulate` — the one front door for a run that simulates two levels at once.

The design note's §3 shape. Every participant is a process spec, you **give what you are not
simulating**, and a rate reads the other participant by name. These cover the two built cells through
the new door, the naming, and the refusals that keep the shape honest.
"""

import pytest

from zombi2 import genomes, joint, species, traits
from zombi2.genomes import family
from zombi2.joint import JointResult
from zombi2.params import PerLineage
from zombi2.species import BirthDeath


def test_birth_death_is_an_unexecuted_spec():
    spec = species.birth_death(birth=1.0, death=0.2, n_extant=50)
    assert isinstance(spec, BirthDeath)
    assert (spec.birth, spec.death, spec.n_extant, spec.total_time) == (1.0, 0.2, 50, None)


def test_a_trait_drives_speciation_through_the_new_door():
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(1.0).scaled_by("traits:size", {"small": 1.0, "large": 4.0}),
            death=0.2, n_extant=120),
        traits.discrete(name="size", states=["small", "large"], switch=0.15), seed=1)
    assert isinstance(r, JointResult) and r.n_extant == 120
    fast = sum(1 for v in r.trait.values.values() if v == "large") / r.n_extant
    assert fast > 0.6, "the fast state should be over-represented at the tips"


def test_gene_content_drives_speciation_through_the_new_door():
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(0.8).scaled_by("genomes:toxin", {"present": 3.0, "absent": 1.0}),
            death=0.15, n_extant=50),
        genomes.genome(origination=0.2, loss=0.3, families=[family('toxin')], initial_families=5),
        seed=3)
    assert r.genome is not None and r.trait is None
    assert r.n_extant == 50


def test_the_retired_spelling_names_its_replacement():
    """`simulate_joint(birth=…, trait=…)` was the whole API. Removing it without a sentence would
    leave a port with "unexpected keyword argument", which names the mistake and not the fix."""
    with pytest.raises(TypeError, match="joint.simulate"):
        joint.simulate_joint(birth=1.0, n_extant=10, seed=1)


def test_an_unnamed_trait_answers_to_the_bare_name():
    """A run holds one trait, so ``"trait"`` names it. A name is what lets a run hold two, and is not
    required for the common case of one."""
    r = joint.simulate(
        species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"a": 1.0, "b": 2.0}),
                            n_extant=20),
        traits.discrete(states=["a", "b"], switch=0.2), seed=1)
    assert r.n_extant == 20


def test_a_named_trait_answers_to_its_name():
    r = joint.simulate(
        species.birth_death(birth=PerLineage(1.0).scaled_by("traits:habitat", {"a": 1.0, "b": 2.0}),
                            n_extant=20),
        traits.discrete(name="habitat", states=["a", "b"], switch=0.2), seed=1)
    assert r.n_extant == 20


def test_a_driver_naming_the_wrong_trait_is_refused():
    with pytest.raises(ValueError, match="traits:habitat"):
        joint.simulate(
            species.birth_death(birth=PerLineage(1.0).scaled_by("traits:size", {"a": 2.0}),
                                n_extant=10),
            traits.discrete(name="habitat", states=["a", "b"], switch=0.2), seed=1)


# --- the shape ------------------------------------------------------------------------------------

def test_a_finished_result_is_not_a_participant():
    """A finished driver is conditioning, and belongs in the driven level's own run."""
    tree = species.simulate_species_tree(birth=1.0, n_extant=10, seed=1)
    with pytest.raises(TypeError, match="conditioning"):
        joint.simulate(species.birth_death(birth=1.0, n_extant=10), tree, seed=1)


def test_the_tree_is_not_given_when_it_is_being_simulated():
    tree = species.simulate_species_tree(birth=1.0, n_extant=10, seed=1).complete_tree
    with pytest.raises(ValueError, match="drop tree="):
        joint.simulate(
            species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"a": 2.0}), n_extant=10),
            traits.discrete(states=["a", "b"], switch=0.2), tree=tree, seed=1)


def test_a_run_on_a_supplied_tree_needs_the_tree():
    with pytest.raises(ValueError, match="pass tree="):
        joint.simulate(genomes.genome(loss=0.1),
                       traits.discrete(name="h", states=["a", "b"], switch=0.2), seed=1)


def test_a_pair_that_has_no_engine_says_so():
    tree = species.simulate_species_tree(birth=1.0, n_extant=10, seed=1).complete_tree
    with pytest.raises(NotImplementedError, match="one genome and one trait"):
        joint.simulate(traits.discrete(name="a", states=["x", "y"], switch=0.2),
                       traits.discrete(name="b", states=["x", "y"], switch=0.2), tree=tree, seed=1)


def test_one_level_must_ride_with_the_tree():
    with pytest.raises(ValueError, match="exactly one level"):
        joint.simulate(species.birth_death(birth=1.0, n_extant=10), seed=1)


def test_one_tree_per_run():
    with pytest.raises(ValueError, match="one species.birth_death"):
        joint.simulate(species.birth_death(birth=1.0, n_extant=10),
                       species.birth_death(birth=1.0, n_extant=10),
                       traits.discrete(states=["a", "b"], switch=0.2), seed=1)
