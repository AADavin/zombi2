"""Transfer-model tests: replacement, growth cap (incl. transfers), distance, self."""

from zombi2 import (
    BirthDeath,
    TransferModel,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.events import EventType


def _mrca_time(a, b):
    ancestors = set()
    node = a
    while node is not None:
        ancestors.add(node)
        node = node.parent
    node = b
    while node not in ancestors:
        node = node.parent
    return node.time


def test_replacement_transfers_curb_growth():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=10, age=3.0, seed=1)
    # one seed family, present in every lineage (no loss); transfers only.
    add = simulate_genomes(tree, transfer=1.0, initial_families=1, seed=2,
                           transfers=TransferModel(replacement=0.0))
    rep = simulate_genomes(tree, transfer=1.0, initial_families=1, seed=2,
                           transfers=TransferModel(replacement=1.0))
    # additive transfers accumulate copies; replacement keeps every recipient at 1
    assert add.profiles.matrix.max() >= 2
    assert rep.profiles.matrix.max() == 1
    # replacement transfers into an occupied recipient must log a compensating loss
    assert any(r.event is EventType.LOSS for r in rep.event_log)


def test_max_family_size_fraction_of_species_bounds_transfers():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=20, age=4.0, seed=1)
    g = simulate_genomes(tree, transfer=1.0, duplication=0.5, initial_families=2,
                         max_family_size=0.25, seed=3)  # cap = round(0.25 * 20) = 5
    assert g.profiles.matrix.max() <= 5


def test_distance_decay_favors_close_recipients():
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=30, age=5.0, seed=1)
    name2node = {n.name: n for n in tree.nodes_preorder()}

    def mean_transfer_distance(decay):
        g = simulate_genomes(tree, transfer=0.8, initial_families=20, seed=5,
                             transfers=TransferModel(distance_decay=decay))
        d = [2.0 * (r.time - _mrca_time(name2node[r.donor], name2node[r.recipient]))
             for r in g.event_log if r.event is EventType.TRANSFER]
        return sum(d) / len(d), len(d)

    m_uniform, n_uniform = mean_transfer_distance(None)
    m_local, n_local = mean_transfer_distance(5.0)
    assert n_uniform > 10 and n_local > 10
    assert m_local < m_uniform  # locality pulls transfers to nearer lineages


def test_self_transfer_acts_as_duplication():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=10, age=3.0, seed=1)
    # no duplication at all; growth within a lineage must come from self-transfers
    g = simulate_genomes(tree, transfer=1.0, duplication=0.0, initial_families=3, seed=4,
                         transfers=TransferModel(allow_self=True))
    self_transfers = [r for r in g.event_log
                      if r.event is EventType.TRANSFER and r.donor == r.recipient]
    assert len(self_transfers) > 0
    assert g.profiles.matrix.max() >= 2  # families gained copies without duplications
