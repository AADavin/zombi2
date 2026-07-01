"""Executable proof that the v1 seams hold.

A new *rate model* and a new *genome representation* both run through the **unchanged**
simulator, rate interface, sampler and profile matrix. This is the compile-time
guarantee behind the spec's extensibility claim.
"""

import numpy as np

from zombi2 import (
    EventRates,
    GenomeSimulator,
    RateModel,
    Simulation,
    SpeciesTreeModel,
    SpeciesTreeSimulator,
    UnorderedGenome,
)
from zombi2.events import EventType
from zombi2.genome import Gene


# --- axis 1: a new RateModel (genome-wise, size-independent) ----------------
class GenomeWiseRates(RateModel):
    def propensities(self, genome, branch, time):
        n = genome.size()
        r = self.rates
        return {
            EventType.DUPLICATION: r.duplication if n > 0 else 0.0,
            EventType.TRANSFER: r.transfer if n > 0 else 0.0,
            EventType.LOSS: r.loss if n > 0 else 0.0,
            EventType.ORIGINATION: r.origination,  # size-independent
        }


def test_genome_wise_rate_model_swap():
    sp = SpeciesTreeModel(1.0, 0.2, 10, age=3.0)
    rates = GenomeWiseRates(EventRates(0.5, 0.2, 0.5, 0.4))
    res = Simulation(sp, rates, seed=4, initial_size=10).run()
    assert res.profiles.matrix.shape[1] == 10
    assert len(res.event_log) > 0


# --- axis 2: a new Genome representation with extra state --------------------
class OrderedListGenome(UnorderedGenome):
    """A different representation that also tracks an explicit gene order.

    A naive design would forget the extra state on clone(); doing it correctly proves
    the simulator (which only calls the interface) never needs to know.
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

    def clone(self):
        new = OrderedListGenome(self.ids)
        for lst in self._genes.values():
            for g in lst:
                new._add(Gene(g.gid, g.family))
        return new


def test_alternative_genome_representation_swap():
    sp = SpeciesTreeModel(1.0, 0.2, 8, age=3.0)
    rates = RateModel(EventRates(0.15, 0.1, 0.2, 0.4))
    tree = SpeciesTreeSimulator().simulate(sp, np.random.default_rng(0))
    gr = GenomeSimulator().simulate(
        tree, rates, np.random.default_rng(1), initial_size=8, genome_factory=OrderedListGenome
    )
    assert gr.leaf_genomes  # produced extant genomes
    for genome in gr.leaf_genomes.values():
        assert isinstance(genome, OrderedListGenome)
        # the extra state stayed consistent with gene content through the whole run
        assert len(genome.order) == genome.size()
        assert sorted(genome.order) == sorted(g.gid for g in genome.genes())
