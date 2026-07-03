"""Forward-in-time species-tree simulation (complete trees with extinct lineages)."""

import numpy as np
import pytest

import zombi2 as z


def _fwd(model, **kw):
    return z.simulate_species_tree(model, direction="forward", **kw)


def _dead_names(tree):
    """Names of branches with no extant descendant (the dead part of the tree)."""
    dead = set()

    def mark(node):
        alive = node.is_extant if node.is_leaf() else any([mark(c) for c in node.children])
        if not alive:
            dead.add(node.name)
        return alive

    mark(tree.root)
    return dead


def test_age_mode_complete_tree():
    tree = _fwd(z.BirthDeath(1.0, 0.4), age=5.0, seed=1)
    assert abs(tree.total_age - 5.0) < 1e-9
    assert tree.root.time == 0.0
    assert len(tree.extant_leaves()) >= 2
    assert all(len(n.children) == 2 for n in tree.internal_nodes())  # binary
    # extant leaves at the present; extinct leaves strictly before it
    for leaf in tree.leaves():
        if leaf.is_extant:
            assert abs(leaf.time - tree.total_age) < 1e-9
        else:
            assert leaf.time < tree.total_age - 1e-12


def test_n_tips_mode_hits_target():
    tree = _fwd(z.BirthDeath(1.0, 0.5), n_tips=20, seed=2)
    assert len(tree.extant_leaves()) == 20
    assert tree.total_age > 0.0


def test_extinction_produces_extinct_leaves():
    n_dead = []
    for s in range(15):
        t = _fwd(z.BirthDeath(1.0, 0.6), age=6.0, seed=s)
        n_dead.append(sum(1 for leaf in t.leaves() if not leaf.is_extant))
    assert np.mean(n_dead) > 0  # extinction leaves dead lineages


def test_yule_has_no_extinction_and_matches_theory():
    lam, age = 1.0, 2.0
    counts = [len(_fwd(z.Yule(lam), age=age, seed=s).extant_leaves())
              for s in range(400)]
    # Yule crown (2 lineages) grown for `age`: E[extant] = 2 e^{λ·age}
    assert abs(np.mean(counts) - 2 * np.exp(lam * age)) / (2 * np.exp(lam * age)) < 0.15
    # no extinction under Yule
    t = _fwd(z.Yule(lam), age=age, seed=1)
    assert all(leaf.is_extant for leaf in t.leaves())


def test_reproducible():
    a = _fwd(z.BirthDeath(1.0, 0.4), age=5.0, seed=7).to_newick()
    b = _fwd(z.BirthDeath(1.0, 0.4), age=5.0, seed=7).to_newick()
    assert a == b


def test_prune_recovers_reconstructed():
    tree = _fwd(z.BirthDeath(1.0, 0.6), n_tips=25, seed=3)
    recon = z.prune(tree)
    assert len(recon.extant_leaves()) == len(tree.extant_leaves()) == 25
    assert all(leaf.is_extant for leaf in recon.leaves())
    assert all(len(n.children) == 2 for n in recon.internal_nodes())


def test_argument_validation():
    with pytest.raises(ValueError):  # neither
        _fwd(z.BirthDeath(1.0, 0.3))
    with pytest.raises(ValueError):  # both
        _fwd(z.BirthDeath(1.0, 0.3), age=5.0, n_tips=10)
    with pytest.raises(ValueError):  # n_tips too small
        _fwd(z.BirthDeath(1.0, 0.3), n_tips=1)
    with pytest.raises(ValueError):  # bad age
        _fwd(z.BirthDeath(1.0, 0.3), age=0.0)
    with pytest.raises(NotImplementedError):  # episodic needs a fixed present -> age mode only
        _fwd(z.EpisodicBirthDeath([1.0], [0.3], []), n_tips=10)


def test_episodic_single_epoch_matches_constant():
    # a one-epoch episodic model == constant BirthDeath forward (same mean extant count)
    def mean_extant(model):
        return np.mean([len(_fwd(model, age=4.0, seed=s).extant_leaves())
                        for s in range(300)])
    epi = z.EpisodicBirthDeath([1.0], [0.4], [])
    const = z.BirthDeath(1.0, 0.4)
    a, b = mean_extant(epi), mean_extant(const)
    assert abs(a - b) / b < 0.12


def test_episodic_recent_mass_extinction_reduces_tips():
    # a high-extinction recent epoch (last 1.0 before present) should leave fewer extant tips
    calm = z.EpisodicBirthDeath([1.0, 1.0], [0.2, 0.2], [1.0])
    crash = z.EpisodicBirthDeath([1.0, 1.0], [2.5, 0.2], [1.0])  # μ=2.5 in the recent epoch
    m_calm = np.mean([len(_fwd(calm, age=5.0, seed=s).extant_leaves())
                      for s in range(200)])
    m_crash = np.mean([len(_fwd(crash, age=5.0, seed=s).extant_leaves())
                       for s in range(200)])
    assert m_crash < m_calm


def test_episodic_incomplete_sampling_marks_unsampled_extant():
    m = z.EpisodicBirthDeath([1.0], [0.3], [], sampling_fraction=0.5)
    # collect trees; with ρ=0.5 some present-day lineages should be unsampled (is_extant=False)
    saw_unsampled = False
    for s in range(30):
        t = _fwd(m, age=5.0, seed=s)
        present = [n for n in t.leaves() if abs(n.time - t.total_age) < 1e-9]
        if any(not n.is_extant for n in present):
            saw_unsampled = True
            break
    assert saw_unsampled


def test_episodic_reproducible():
    m = z.EpisodicBirthDeath([1.0, 1.6], [0.3, 0.6], [2.0], sampling_fraction=0.8)
    a = _fwd(m, age=5.0, seed=7).to_newick()
    b = _fwd(m, age=5.0, seed=7).to_newick()
    assert a == b


# --- FBD / serial (through-time) sampling ------------------------------------

def _fossils(tree):
    return [n for n in tree.leaves() if n.sampled and not n.is_extant]


def test_fbd_no_fossils_when_psi_zero():
    t = _fwd(z.BirthDeath(1.0, 0.4, fossilization=0.0),
                                        age=5.0, seed=1)
    assert _fossils(t) == []
    assert all(len(n.children) == 2 for n in t.internal_nodes())  # still binary


def test_fbd_fossils_scale_with_psi():
    def mean_fossils(psi):
        return np.mean([len(_fossils(_fwd(
            z.BirthDeath(1.0, 0.4, fossilization=psi), age=5.0, seed=s)))
            for s in range(40)])
    assert mean_fossils(0.0) == 0
    assert mean_fossils(0.2) < mean_fossils(0.6)


def test_fbd_fossils_are_dated_and_sampled():
    t = _fwd(
        z.BirthDeath(1.0, 0.5, fossilization=0.6), age=6.0, seed=2)
    fossils = _fossils(t)
    assert fossils
    for f in fossils:
        assert f.sampled and not f.is_extant
        assert f.time < t.total_age - 1e-9  # a past sample, before the present


def test_fbd_sampled_tree_extraction():
    t = _fwd(
        z.BirthDeath(1.0, 0.5, fossilization=0.5, sampling_fraction=0.9), age=6.0, seed=1)
    n_sampled = sum(1 for n in t.leaves() if n.sampled)
    samp = z.prune(t, keep="sampled")
    assert len(samp.leaves()) == n_sampled           # fossils + extant samples
    assert all(len(nd.children) == 2 for nd in samp.internal_nodes())
    # prune keeps only extant sampled tips (no fossils)
    recon = z.prune(t)
    assert len(recon.leaves()) == len(t.extant_leaves())
    assert len(recon.leaves()) <= len(samp.leaves())


def test_fbd_n_tips_mode_allowed():
    t = _fwd(
        z.BirthDeath(1.0, 0.3, fossilization=0.3), n_tips=15, seed=3)
    assert len(t.extant_leaves()) == 15  # constant-rate FBD supports n_tips


def test_fbd_reproducible():
    m = z.BirthDeath(1.0, 0.4, fossilization=0.4, sampling_fraction=0.9)
    a = _fwd(m, age=5.0, seed=7).to_newick()
    b = _fwd(m, age=5.0, seed=7).to_newick()
    assert a == b


# --- sampled ancestors (removal r < 1) ---------------------------------------

def _sampled_ancestors(tree):
    return [n for n in tree.nodes() if n.sampled and len(n.children) == 1]


def test_sampled_ancestors_only_when_removal_below_one():
    common = dict(fossilization=0.6, sampling_fraction=1.0)
    t_removed = _fwd(
        z.BirthDeath(1.0, 0.4, removal=1.0, **common), age=6.0, seed=1)
    t_kept = _fwd(
        z.BirthDeath(1.0, 0.4, removal=0.0, **common), age=6.0, seed=1)
    assert _sampled_ancestors(t_removed) == []
    assert len(_sampled_ancestors(t_kept)) > 0


def test_gene_sim_passes_through_sampled_ancestors():
    tree = _fwd(
        z.BirthDeath(1.0, 0.4, fossilization=0.6, removal=0.0), age=6.0, seed=1)
    assert _sampled_ancestors(tree)  # the tree really has degree-two nodes
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.2, loss=0.15,
                           origination=0.5, initial_size=20, max_family_size=0.5, seed=42)
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}


def test_prune_sampled_keeps_sampled_ancestors():
    tree = _fwd(
        z.BirthDeath(1.0, 0.4, fossilization=0.6, removal=0.0, sampling_fraction=0.9),
        age=6.0, seed=1)
    samp = z.prune(tree, keep="sampled")
    n_leaf_samples = sum(1 for n in tree.leaves() if n.sampled)
    assert len(samp.leaves()) >= n_leaf_samples          # fossil/extant leaf samples kept
    assert any(len(n.children) == 1 and n.sampled for n in samp.nodes())  # SAs preserved
    # the extant-only reconstructed tree suppresses sampled ancestors
    recon = z.prune(tree)
    assert all(len(n.children) == 2 for n in recon.internal_nodes())


def test_removal_validation():
    with pytest.raises(ValueError):
        _fwd(
            z.BirthDeath(1.0, 0.4, fossilization=0.5, removal=1.5), age=5.0)


# --- episodic FBD ------------------------------------------------------------

def test_episodic_fbd_produces_fossils():
    m = z.EpisodicBirthDeath(
        birth=[1.0, 1.4], death=[0.3, 0.5], fossilization=[0.4, 0.4], shifts=[2.0],
        sampling_fraction=0.9, removal=0.5)
    t = _fwd(m, age=6.0, seed=3)
    assert len(_fossils(t)) > 0
    assert len(t.extant_leaves()) >= 2


def test_episodic_fbd_requires_age_mode():
    m = z.EpisodicBirthDeath([1.0], [0.3], [], fossilization=[0.3])
    with pytest.raises(NotImplementedError):
        _fwd(m, n_tips=10)


def test_episodic_fbd_reproducible():
    m = z.EpisodicBirthDeath(
        birth=[1.0, 1.4], death=[0.3, 0.5], fossilization=[0.3, 0.5], shifts=[2.0],
        sampling_fraction=0.8, removal=0.5)
    a = _fwd(m, age=5.0, seed=7).to_newick()
    b = _fwd(m, age=5.0, seed=7).to_newick()
    assert a == b


def test_forward_tree_feeds_gene_sim_with_ghost_transfers():
    tree = _fwd(z.BirthDeath(1.0, 0.6), n_tips=40, seed=8)
    dead = _dead_names(tree)
    assert dead  # forward tree has a dead part
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.4, loss=0.15,
                           origination=0.5, initial_size=30, max_family_size=0.5, seed=42)
    # profiles only over extant species
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}
    # at least one transfer should involve a dead (extinct) branch — transfer from the dead,
    # for free, with no ghost-grafting step
    involved = any(
        (r.donor in dead or r.recipient in dead or r.branch in dead)
        for r in g.event_log if r.event is z.EventType.TRANSFER
    )
    assert involved
