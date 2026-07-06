"""Executable proof that the v1 seams hold.

A new *rate model* and a new *genome representation* both run through the **unchanged**
simulator, rate interface, sampler and profile matrix.
"""

import numpy as np

from zombi2 import (
    BirthDeath,
    EventWeight,
    Exponential,
    FamilySampledRates,
    Gamma,
    GenomeSimulator,
    RateModel,
    SharedRates,
    UnorderedGenome,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes.events import EventType


# --- axis 1: a new RateModel (genome-wise, size-independent totals) ----------
class PerGenomeRates(RateModel):
    def __init__(self, duplication, transfer, loss, origination):
        self.d, self.t, self.l, self.o = duplication, transfer, loss, origination

    def event_weights(self, genome, branch, time):
        out = []
        if genome.size() > 0:
            out += [
                EventWeight(EventType.DUPLICATION, None, self.d),
                EventWeight(EventType.TRANSFER, None, self.t),
                EventWeight(EventType.LOSS, None, self.l),
            ]
        out.append(EventWeight(EventType.ORIGINATION, None, self.o))
        return out


def test_genome_wise_rate_model_swap():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=4)
    genomes = simulate_genomes(tree, PerGenomeRates(0.5, 0.2, 0.5, 0.4),
                               initial_families=10, seed=4)
    assert genomes.profiles.matrix.shape[1] == 10
    assert len(genomes.event_log) > 0


# --- per-family sampled rates (a stateful rate model) -----------------------
def test_family_sampled_rates_differ_and_cache():
    fs = FamilySampledRates(duplication=Exponential(0.5), transfer=0.0,
                            loss=Exponential(0.3), origination=0.0)
    fs.bind(np.random.default_rng(0))
    r1, r2 = fs.rates_for("1"), fs.rates_for("2")
    assert r1 != r2                      # different families -> different rates
    assert fs.rates_for("1") == r1       # cached and stable for the life of a family


def test_family_sampled_accepts_float_callable_and_dist():
    fs = FamilySampledRates(duplication=lambda rng: 0.1, transfer=0.05,
                            loss=Gamma(2, 0.1), origination=0.3)
    fs.bind(np.random.default_rng(0))
    d, t, l = fs.rates_for("x")
    assert d == 0.1 and t == 0.05 and l > 0


def test_family_sampled_full_run_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    make = lambda: FamilySampledRates(duplication=Exponential(0.2), transfer=Exponential(0.1),
                                      loss=Exponential(0.25), origination=0.5)
    a = simulate_genomes(tree, make(), initial_families=15, seed=2)
    b = simulate_genomes(tree, make(), initial_families=15, seed=2)
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)
    assert a.profiles.matrix.shape[1] == 10


# --- axis 2: a new Genome representation with extra state --------------------
class OrderedListGenome(UnorderedGenome):
    """A different representation that also tracks an explicit gene order.

    It only overrides the low-level add/remove hooks; every operation (including the
    inherited ``clone_reminting`` used at speciation) then maintains the extra state, so
    the simulator never needs to know about it.
    """

    def __init__(self, ids):
        super().__init__(ids)
        self.order: list[str] = []

    def _add(self, gene):
        super()._add(gene)
        self.order.append(gene.gid)

    def _remove(self, gene):
        super()._remove(gene)
        self.order.remove(gene.gid)


def test_alternative_genome_representation_swap():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=0)
    gr = GenomeSimulator().simulate(
        tree, SharedRates(0.15, 0.1, 0.2, 0.4),
        np.random.default_rng(1), initial_size=8, genome_factory=OrderedListGenome,
    )
    assert gr.leaf_genomes
    for genome in gr.leaf_genomes.values():
        assert isinstance(genome, OrderedListGenome)
        assert len(genome.order) == genome.size()
        assert sorted(genome.order) == sorted(g.gid for g in genome.genes())
