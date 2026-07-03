"""Ghost (dead) lineage grafting — un-pruning the reconstructed species tree."""

import numpy as np
import pytest

import zombi2 as z
from zombi2.tree import Tree, TreeNode


def _recon(birth, death, n=25, age=5.0, seed=1):
    return z.simulate_species_tree(z.BirthDeath(birth, death), n_tips=n, age=age, seed=seed)


def _prune_to_extant(node):
    """Rebuild the subtree keeping only paths to extant leaves, suppressing degree-two nodes
    and preserving names/times — should reproduce the reconstructed tree."""
    if node.is_leaf():
        if node.is_extant:
            return TreeNode(name=node.name, time=node.time, is_extant=True)
        return None
    kept = [k for k in (_prune_to_extant(c) for c in node.children) if k is not None]
    if not kept:
        return None
    if len(kept) == 1:  # suppress this node
        return kept[0]
    m = TreeNode(name=node.name, time=node.time)
    for k in kept:
        m.add_child(k)
    return m


def _extant_newick(tree):
    root = _prune_to_extant(tree.root)
    return Tree(root, tree.total_age).to_newick()


def _dead_leaves(tree):
    return [n for n in tree.leaves() if not n.is_extant]


# --- sanity check 1: Yule / no extinction -> no ghosts -----------------------

def test_yule_adds_no_ghosts():
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=20, age=5.0, seed=3)
    before = tree.to_newick()
    n_before = len(tree.nodes())
    z.add_ghost_lineages(tree, z.Yule(1.0), seed=7)
    assert tree.to_newick() == before
    assert len(tree.nodes()) == n_before
    assert _dead_leaves(tree) == []


# --- sanity check 2: pruning invariant --------------------------------------

def test_pruning_recovers_reconstructed_tree():
    tree = _recon(1.0, 0.6, n=30, seed=5)
    original = tree.to_newick()
    z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.6), seed=11)
    assert len(_dead_leaves(tree)) > 0            # ghosts were actually added
    assert _extant_newick(tree) == original        # ...but pruning them off is exact


# --- sanity check 3: extant tips untouched ----------------------------------

def test_extant_tips_untouched():
    tree = _recon(1.0, 0.5, n=25, seed=2)
    before = {n.name: n.time for n in tree.extant_leaves()}
    z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.5), seed=9)
    after = {n.name: n.time for n in tree.extant_leaves()}
    assert after == before
    assert len(tree.extant_leaves()) == 25
    assert all(not n.is_extant for n in _dead_leaves(tree))


# --- sanity check 4: density scales with extinction -------------------------

def test_ghost_density_increases_with_extinction():
    def mean_ghosts(death):
        counts = []
        for s in range(8):
            tree = _recon(1.0, death, n=30, seed=100 + s)
            z.add_ghost_lineages(tree, z.BirthDeath(1.0, death), seed=200 + s)
            counts.append(len(_dead_leaves(tree)))
        return np.mean(counts)

    assert mean_ghosts(0.0) == 0
    assert mean_ghosts(0.3) < mean_ghosts(0.7)


# --- sanity check 5: reproducibility ----------------------------------------

def test_reproducible_augmentation():
    def run():
        tree = _recon(1.0, 0.5, n=25, seed=4)
        z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.5), seed=42)
        return tree.to_newick()
    assert run() == run()


# --- structural checks ------------------------------------------------------

def test_ghost_nodes_are_non_extant_and_bounded_in_time():
    tree = _recon(1.0, 0.6, n=30, seed=6)
    z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.6), seed=13)
    for n in tree.nodes():
        assert n.time <= tree.total_age + 1e-9
        if n.name.startswith("e"):          # grafted extinct leaves are named e*
            assert not n.is_extant
    assert any(n.name.startswith("e") for n in tree.leaves())  # some ghosts were added
    # every internal node stays binary (the forward sim relies on this)
    assert all(len(n.children) == 2 for n in tree.internal_nodes())
    # extinct ghost leaves die strictly before the present
    assert all(leaf.time < tree.total_age for leaf in _dead_leaves(tree))


def _episodic():
    return z.EpisodicBirthDeath(birth=[1.0, 1.6], death=[0.4, 0.7], shifts=[2.0],
                                sampling_fraction=0.7)


def test_episodic_pruning_invariant():
    tree = z.simulate_species_tree(_episodic(), n_tips=30, age=6.0, seed=5)
    original = tree.to_newick()
    z.add_ghost_lineages(tree, _episodic(), seed=11)
    assert len(_dead_leaves(tree)) > 0
    assert _extant_newick(tree) == original


def test_episodic_incomplete_sampling_adds_unsampled_extant_ghosts():
    # with ρ<1 some ghosts survive to the present but are unsampled (dead leaf at total_age)
    tree = z.simulate_species_tree(_episodic(), n_tips=40, age=6.0, seed=3)
    z.add_ghost_lineages(tree, _episodic(), seed=8)
    at_present = [g for g in _dead_leaves(tree) if abs(g.time - tree.total_age) < 1e-9]
    before_present = [g for g in _dead_leaves(tree) if g.time < tree.total_age - 1e-9]
    assert len(at_present) > 0        # unsampled-extant ghosts (from ρ<1)
    assert len(before_present) > 0    # extinct-before-present ghosts (from μ>0)


def test_episodic_reproducible():
    def run():
        tree = z.simulate_species_tree(_episodic(), n_tips=30, age=6.0, seed=4)
        z.add_ghost_lineages(tree, _episodic(), seed=42)
        return tree.to_newick()
    assert run() == run()


def _mean_dead(method, make_model, seeds):
    counts = []
    for s in seeds:
        tree = z.simulate_species_tree(make_model(), n_tips=30, age=6.0, seed=s)
        z.add_ghost_lineages(tree, make_model(), method=method, seed=1000 + s)
        counts.append(len(_dead_leaves(tree)))
    return np.mean(counts)


def test_htransform_matches_rejection_constant():
    make = lambda: z.BirthDeath(1.0, 0.7)
    rej = _mean_dead("rejection", make, range(30))
    ht = _mean_dead("htransform", make, range(30))
    assert abs(rej - ht) / rej < 0.25   # statistically equivalent samplers


def test_htransform_matches_rejection_episodic():
    make = lambda: z.EpisodicBirthDeath([1.0, 1.6], [0.3, 0.8], [3.0], sampling_fraction=0.6)
    rej = _mean_dead("rejection", make, range(30))
    ht = _mean_dead("htransform", make, range(30))
    assert abs(rej - ht) / rej < 0.25


def test_htransform_pruning_invariant():
    tree = _recon(1.0, 0.7, n=30, seed=5)
    original = tree.to_newick()
    z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.7), method="htransform", seed=11)
    assert len(_dead_leaves(tree)) > 0
    assert _extant_newick(tree) == original


def test_htransform_reproducible():
    def run():
        tree = _recon(1.0, 0.6, n=25, seed=4)
        z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.6), method="htransform", seed=42)
        return tree.to_newick()
    assert run() == run()


def test_htransform_episodic_incomplete_sampling():
    m = z.EpisodicBirthDeath([1.0, 1.6], [0.3, 0.8], [3.0], sampling_fraction=0.6)
    tree = z.simulate_species_tree(m, n_tips=40, age=6.0, seed=3)
    z.add_ghost_lineages(tree, m, method="htransform", seed=8)
    dead = _dead_leaves(tree)
    assert any(abs(g.time - tree.total_age) < 1e-9 for g in dead)   # unsampled-extant
    assert any(g.time < tree.total_age - 1e-9 for g in dead)        # extinct


def test_invalid_method_raises():
    tree = _recon(1.0, 0.5, seed=1)
    with pytest.raises(ValueError):
        z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.5), method="nope", seed=1)


def test_unsupported_model_raises():
    tree = _recon(1.0, 0.5, seed=1)

    class NotAModel:
        pass

    with pytest.raises(NotImplementedError):
        z.add_ghost_lineages(tree, NotAModel(), seed=1)


# --- integration: the forward gene sim uses ghosts as transfer partners -----

def test_forward_sim_runs_on_augmented_tree():
    tree = _recon(1.0, 0.6, n=30, seed=8)
    z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.6), seed=8)
    ghost_branches = {n.name for n in tree.nodes() if not n.is_extant}
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.4, loss=0.15,
                           origination=0.5, initial_size=30, max_family_size=0.5, seed=42)
    # profiles only over sampled (extant) species
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}
    # with transfers on a ghost-laden tree, at least one transfer should involve a ghost branch
    involved = any(
        (r.donor in ghost_branches or r.recipient in ghost_branches or r.branch in ghost_branches)
        for r in g.event_log if r.event is z.EventType.TRANSFER
    )
    assert involved
