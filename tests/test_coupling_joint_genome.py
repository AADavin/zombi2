"""Coupling slice 3 — gene content drives speciation, grown jointly (P(Species, Genomes)).

The genome half of joint: `mod.DrivenBy("genomes:count", …)` / `mod.DrivenBy("genomes:<family>", …)`
with a live genome grown by `joint.simulate_joint(genome=genomes.unordered(...))`. Covers named
families (the referenceable handle), the genome process spec, the result shape, determinism, and the
two gene-content-dependent-diversification signals (count and named-presence).
"""

import statistics

import pytest

from zombi2 import genomes, joint
from zombi2.genomes import UnorderedGenome
from zombi2.joint import JointResult
from zombi2.rates import modifiers as mod
from zombi2.species import simulate_species_tree


# --- named families in the genome level (the referenceable handle) --------------------------------

def test_named_families_seed_and_track():
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=3).complete_tree
    res = genomes.simulate_genomes_unordered(tree, duplication=0.1, loss=0.3,
                                             initial_families=3, families=["toxin", "operon"], seed=1)
    assert set(res.family_names) == {"toxin", "operon"}
    assert res.has_family(tree.root, "toxin")           # seeded at the crown
    with pytest.raises(KeyError):
        res.has_family(tree.root, "absent_name")


def test_named_families_default_path_byte_identical():
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    a = genomes.simulate_genomes_unordered(tree, loss=0.2, initial_families=4, seed=7)
    b = genomes.simulate_genomes_unordered(tree, loss=0.2, initial_families=4, families=None, seed=7)
    assert [(e.kind, e.copy, e.family) for e in a.events] == [(e.kind, e.copy, e.family) for e in b.events]


def test_named_families_validate():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError):
        genomes.simulate_genomes_unordered(tree, families=["a", "a"], seed=1)     # duplicate
    with pytest.raises(ValueError):
        genomes.simulate_genomes_unordered(tree, families=[""], seed=1)           # empty name


# --- the genome process spec ----------------------------------------------------------------------

def test_unordered_spec_is_unexecuted_bundle():
    spec = genomes.unordered(duplication=0.1, loss=0.2, origination=0.3, families=["toxin"])
    assert isinstance(spec, UnorderedGenome)
    assert spec.families == ("toxin",)


def test_unordered_spec_validates():
    with pytest.raises(ValueError):
        genomes.unordered(families=["a", "a"])
    with pytest.raises(ValueError):
        genomes.unordered(initial_families=-1)


# --- the result shape -----------------------------------------------------------------------------

def _count_joint(curve=lambda n: 1.0 + 0.2 * n, n_extant=150, seed=1):
    return joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("genomes:count", curve), death=0.1,
        genome=genomes.unordered(duplication=0.3, loss=0.3, origination=0.3, initial_families=3),
        n_extant=n_extant, seed=seed)


def test_joint_genome_result_shape():
    res = _count_joint(n_extant=120, seed=2)
    assert isinstance(res, JointResult)
    assert res.n_extant == 120
    assert res.trait is None and res.genome is not None
    assert res.genome.__class__.__name__ == "GenomesResult"
    # a genome is recorded at every node (the profiles derive from the extant tips)
    assert set(res.genome.genomes) == set(res.complete_tree.nodes)
    assert res.genome.profiles is not None


def test_joint_genome_writes_both_levels(tmp_path):
    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("genomes:toxin", {"present": 3.0, "absent": 1.0}), death=0.1,
        genome=genomes.unordered(duplication=0.3, loss=0.3, families=["toxin"]),
        n_extant=80, seed=3)
    res.write(tmp_path)
    for f in ("species_complete.nwk", "species_extant.nwk", "species_events.tsv",
              "genome_events.tsv", "profiles.tsv"):
        assert (tmp_path / f).exists(), f"missing {f}"


def test_joint_genome_is_deterministic():
    a, b = _count_joint(seed=9), _count_joint(seed=9)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]
    assert [(e.time, e.kind, e.lineage, e.copy) for e in a.genome.events] == \
           [(e.time, e.kind, e.lineage, e.copy) for e in b.genome.events]


# --- the gene-content-dependent-diversification signals -------------------------------------------

def test_gene_count_drives_diversification():
    # bigger genomes speciate faster → coupled tips carry larger genomes than a flat (neutral) curve
    def mean_size(curve):
        sizes = []
        for s in (1, 2, 3):
            r = _count_joint(curve, seed=s)
            sizes.append(statistics.mean(len(r.genome.genomes[n.id]) for n in r.complete_tree.extant()))
        return statistics.mean(sizes)
    coupled = mean_size(lambda n: 1.0 + 0.3 * n)
    neutral = mean_size(lambda n: 2.0)              # flat curve: driven but no effect = neutral null
    assert coupled > neutral, f"count coupling gave no size signal: {coupled:.2f} vs {neutral:.2f}"


def test_named_family_presence_drives_diversification():
    def frac_toxin(birth, seed):
        r = joint.simulate_joint(birth=birth, death=0.1,
            genome=genomes.unordered(duplication=0.4, loss=0.4, origination=0.1, families=["toxin"]),
            n_extant=200, seed=seed)
        tips = [n.id for n in r.complete_tree.extant()]
        return sum(r.genome.has_family(i, "toxin") for i in tips) / len(tips)
    coupled = statistics.mean(
        frac_toxin(1.0 * mod.DrivenBy("genomes:toxin", {"present": 5.0, "absent": 1.0}), s) for s in (1, 2, 3))
    neutral = statistics.mean(
        frac_toxin(1.0 * mod.DrivenBy("genomes:toxin", {"present": 1.0, "absent": 1.0}), s) for s in (1, 2, 3))
    assert coupled > neutral + 0.15, f"toxin coupling gave no signal: {coupled:.2f} vs {neutral:.2f}"


def test_total_time_mode():
    res = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("genomes:count", lambda n: 1.0 + 0.1 * n), death=0.1,
        genome=genomes.unordered(origination=0.4, loss=0.2, initial_families=2),
        total_time=3.0, seed=4)
    assert all(n.end_time == pytest.approx(3.0) for n in res.complete_tree.extant())


# --- validation -----------------------------------------------------------------------------------

def test_exactly_one_driver():
    with pytest.raises(TypeError, match="exactly one driver"):
        joint.simulate_joint(birth=1.0 * mod.DrivenBy("trait", {"a": 1.0}), n_extant=10, seed=1)  # neither
    with pytest.raises(TypeError, match="exactly one driver"):
        joint.simulate_joint(
            birth=1.0 * mod.DrivenBy("genomes:count", lambda n: n),
            trait=__import__("zombi2.traits", fromlist=["discrete"]).discrete(states=["a", "b"], switch=0.1),
            genome=genomes.unordered(origination=0.1), n_extant=10, seed=1)  # both


def test_unseeded_named_family_rejected():
    with pytest.raises(ValueError, match="not.*seeded"):
        joint.simulate_joint(
            birth=1.0 * mod.DrivenBy("genomes:toxin", {"present": 2.0, "absent": 1.0}),
            genome=genomes.unordered(origination=0.1),   # no families=["toxin"]
            n_extant=10, seed=1)


def test_trait_source_on_genome_joint_rejected():
    with pytest.raises(ValueError, match="genomes:"):
        joint.simulate_joint(
            birth=1.0 * mod.DrivenBy("trait", {"a": 2.0}),
            genome=genomes.unordered(origination=0.1, families=["a"]),
            n_extant=10, seed=1)
