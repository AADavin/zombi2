"""Joint slice 3 — gene content drives speciation, grown jointly (P(Species, Genomes)).

The genome half of joint: `.scaled_by("genomes:count", …)` / `.scaled_by("genomes:<family>", …)`
with a live genome simulated by `joint.simulate(species.birth_death(...), genomes.genome(...))`. Covers named
families (the referenceable handle), the genome process spec, the result shape, determinism, and the
two gene-content-dependent-diversification signals (count and named-presence).
"""

import statistics

import pytest

from zombi2 import genomes, joint, species
from zombi2.genomes import FamilyGenome, family
from zombi2.joint import JointResult
from zombi2.params import PerLineage
from zombi2.species import simulate_species_tree


# --- named families in the genome level (the referenceable handle) --------------------------------

def test_named_families_seed_and_track():
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=3).complete_tree
    res = genomes.simulate_genomes_family(tree, duplication=0.1, loss=0.3,
                                          initial_families=3, families=[family("toxin"), family("operon")], seed=1)
    assert set(res.family_names) == {"toxin", "operon"}
    assert res.has_family(tree.root, "toxin")           # seeded at the crown
    with pytest.raises(KeyError):
        res.has_family(tree.root, "absent_name")


def test_named_families_default_path_byte_identical():
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    a = genomes.simulate_genomes_family(tree, loss=0.2, initial_families=4, seed=7)
    b = genomes.simulate_genomes_family(tree, loss=0.2, initial_families=4, families=[], seed=7)
    assert [(e.kind, e.copy, e.family) for e in a.edges] == [(e.kind, e.copy, e.family) for e in b.edges]


def test_named_families_validate():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError):
        genomes.simulate_genomes_family(tree, families=[family("a"), family("a")], seed=1)     # duplicate
    with pytest.raises(ValueError):
        genomes.simulate_genomes_family(tree, families=[family("")], seed=1)           # empty name


# --- the genome process spec ----------------------------------------------------------------------

def test_family_spec_is_unexecuted_bundle():
    spec = genomes.genome(duplication=0.1, loss=0.2, origination=0.3, families=[family("toxin")])
    assert isinstance(spec, FamilyGenome)
    assert spec.family_names == ("toxin",)


def test_family_spec_validates():
    with pytest.raises(ValueError):
        genomes.genome(families=[family("a"), family("a")])
    with pytest.raises(ValueError):
        genomes.genome(initial_families=-1)


# --- the result shape -----------------------------------------------------------------------------

def _count_joint(curve=lambda n: 1.0 + 0.2 * n, n_extant=150, seed=1):
    return joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:count", curve), death=0.1, n_extant=n_extant), genomes.genome(duplication=0.3, loss=0.3, origination=0.3, initial_families=3), seed=seed)


def test_joint_genome_result_shape():
    res = _count_joint(n_extant=120, seed=2)
    assert isinstance(res, JointResult)
    assert res.n_extant == 120
    assert res.trait is None and res.genome is not None
    assert res.genome.__class__.__name__ == "FamilyGenomesResult"
    # a genome is recorded at every node (the profiles derive from the extant tips)
    assert set(res.genome.node_genomes) == set(res.complete_tree.nodes)
    assert res.genome.profiles is not None


def test_joint_genome_writes_both_levels(tmp_path):
    res = joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:toxin", {"present": 3.0, "absent": 1.0}), death=0.1, n_extant=80), genomes.genome(duplication=0.3, loss=0.3, families=[family("toxin")]), seed=3)
    res.write(tmp_path)
    for f in ("species_complete.nwk", "species_extant.nwk", "species_events.tsv",
              "genome_events.tsv", "profiles.tsv"):
        assert (tmp_path / f).exists(), f"missing {f}"


def test_joint_genome_is_deterministic():
    a, b = _count_joint(seed=9), _count_joint(seed=9)
    assert [(e.time, e.kind, e.node) for e in a.events] == \
           [(e.time, e.kind, e.node) for e in b.events]
    assert [(e.time, e.kind, e.lineage, e.copy) for e in a.genome.edges] == \
           [(e.time, e.kind, e.lineage, e.copy) for e in b.genome.edges]


# --- the gene-content-dependent-diversification signals -------------------------------------------

def test_gene_count_drives_diversification():
    # bigger genomes speciate faster → driven tips carry larger genomes than a flat (neutral) curve
    def mean_size(curve):
        sizes = []
        for s in (1, 2, 3):
            r = _count_joint(curve, seed=s)
            sizes.append(statistics.mean(len(r.genome.node_genomes[n.id]) for n in (r.complete_tree.nodes[_i] for _i in r.complete_tree.extant_leaves())))
        return statistics.mean(sizes)
    driven = mean_size(lambda n: 1.0 + 0.3 * n)
    neutral = mean_size(lambda n: 2.0)              # flat curve: driven but no effect = neutral null
    assert driven > neutral, f"count driving gave no size signal: {driven:.2f} vs {neutral:.2f}"


def test_named_family_presence_drives_diversification():
    def frac_toxin(birth, seed):
        r = joint.simulate(species.birth_death(birth=birth, death=0.1, n_extant=200), genomes.genome(duplication=0.4, loss=0.4, origination=0.1, families=[family("toxin")]), seed=seed)
        tips = list(r.complete_tree.extant_leaves())
        return sum(r.genome.has_family(i, "toxin") for i in tips) / len(tips)
    driven = statistics.mean(
        frac_toxin(PerLineage(1.0).scaled_by("genomes:toxin", {"present": 5.0, "absent": 1.0}), s) for s in (1, 2, 3))
    neutral = statistics.mean(
        frac_toxin(PerLineage(1.0).scaled_by("genomes:toxin", {"present": 1.0, "absent": 1.0}), s) for s in (1, 2, 3))
    assert driven > neutral + 0.15, f"toxin driving gave no signal: {driven:.2f} vs {neutral:.2f}"


def test_total_time_mode():
    res = joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:count", lambda n: 1.0 + 0.1 * n), death=0.1, total_time=3.0), genomes.genome(origination=0.4, loss=0.2, initial_families=2), seed=4)
    assert all(n.end_time == pytest.approx(3.0) for n in (res.complete_tree.nodes[_i] for _i in res.complete_tree.extant_leaves()))


# --- validation -----------------------------------------------------------------------------------

def test_exactly_one_level_rides_with_the_tree():
    """A joint run simulates the tree with one other level. Neither none nor two."""
    from zombi2 import traits

    bare = species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:count", lambda n: n),
                               n_extant=10)
    with pytest.raises(ValueError, match="exactly one level"):
        joint.simulate(bare, seed=1)                                          # neither
    with pytest.raises(ValueError, match="exactly one level"):
        joint.simulate(bare, traits.discrete(states=["a", "b"], switch=0.1),
                       genomes.genome(origination=0.1), seed=1)               # both


def test_undeclared_named_family_rejected():
    with pytest.raises(ValueError, match="not.*declared"):
        joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("genomes:toxin", {"present": 2.0, "absent": 1.0}), n_extant=10), genomes.genome(origination=0.1), seed=1)


def test_trait_source_on_genome_joint_rejected():
    with pytest.raises(ValueError, match="genomes:"):
        joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"a": 2.0}), n_extant=10), genomes.genome(origination=0.1, families=[family("a")]), seed=1)
