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
                           origination=0.5, initial_families=20, max_family_size=0.5, seed=42)
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


# --- mass extinctions (instantaneous tree-wide survival pulses) --------------

def _pulse_deaths(tree, pulse_time):
    """Extinct leaves that died exactly at ``pulse_time`` — the mass-extinction victims
    (background extinction lands at continuous random times, never exactly on the pulse)."""
    return [n for n in tree.leaves()
            if not n.is_extant and abs(n.time - pulse_time) < 1e-9]


def test_mass_extinction_reduces_extant_tips():
    # a severe pulse partway through should leave fewer extant tips than no pulse
    calm = z.BirthDeath(1.0, 0.3)
    crash = z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.5, 0.9)])  # 90% die at age 2.5
    m_calm = np.mean([len(_fwd(calm, age=5.0, seed=s).extant_leaves()) for s in range(200)])
    m_crash = np.mean([len(_fwd(crash, age=5.0, seed=s).extant_leaves()) for s in range(200)])
    assert m_crash < m_calm


def test_mass_extinction_kills_at_the_pulse_time():
    # victims become extinct leaves exactly at the pulse instant (present - age = 5 - 2 = 3)
    m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.0, 0.8)])
    saw_victims = False
    for s in range(30):
        t = _fwd(m, age=5.0, seed=s)
        if _pulse_deaths(t, 3.0):
            saw_victims = True
            # every victim is a dead leaf sitting strictly before the present
            assert all((not v.is_extant) and v.time < t.total_age for v in _pulse_deaths(t, 3.0))
            break
    assert saw_victims


def test_mass_extinction_multiple_pulses():
    # two pulses (ages 3.0 and 1.0 before a present at 5.0 -> tree-times 2.0 and 4.0)
    m = z.BirthDeath(1.2, 0.2, mass_extinctions=[(3.0, 0.7), (1.0, 0.7)])
    saw_first = saw_second = False
    for s in range(40):
        t = _fwd(m, age=5.0, seed=s)
        saw_first = saw_first or bool(_pulse_deaths(t, 2.0))
        saw_second = saw_second or bool(_pulse_deaths(t, 4.0))
        if saw_first and saw_second:
            break
    assert saw_first and saw_second


def test_mass_extinction_full_wipe_is_rejected():
    # fraction == 1.0 kills every lineage -> the run always dies out -> conditioning fails
    m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.5, 1.0)])
    with pytest.raises(RuntimeError):
        _fwd(m, age=5.0, seed=1, max_attempts=30)


def test_mass_extinction_requires_age_mode():
    m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.0, 0.5)])
    with pytest.raises(NotImplementedError):  # a pulse age needs a fixed present
        _fwd(m, n_tips=10)


def test_mass_extinction_age_must_precede_crown():
    m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(5.0, 0.5)])
    with pytest.raises(ValueError):  # pulse at/after the crown age is undefined
        _fwd(m, age=5.0, seed=1)


@pytest.mark.parametrize("frac", [0.0, 1.5, -0.1])
def test_mass_extinction_fraction_validation(frac):
    with pytest.raises(ValueError):
        _fwd(z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.0, frac)]), age=5.0)


def test_mass_extinction_backward_is_rejected():
    m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.0, 0.5)])
    with pytest.raises(ValueError):  # not represented in the reconstructed tree
        z.simulate_species_tree(m, direction="backward", n_tips=10, age=5.0)


def test_mass_extinction_reproducible():
    m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(2.5, 0.6), (1.0, 0.5)])
    a = _fwd(m, age=5.0, seed=7).to_newick()
    b = _fwd(m, age=5.0, seed=7).to_newick()
    assert a == b


def test_mass_extinction_composes_with_episodic():
    # a pulse layered on an episodic (skyline) background still fires at its instant
    m = z.EpisodicBirthDeath([1.0, 1.0], [0.2, 0.2], [2.0], mass_extinctions=[(1.0, 0.8)])
    saw = False
    for s in range(30):
        t = _fwd(m, age=5.0, seed=s)   # pulse tree-time = 5 - 1 = 4
        if _pulse_deaths(t, 4.0):
            saw = True
            break
    assert saw


def test_mass_extinction_feeds_gene_sim():
    # the pulse's dead lineages are ordinary extinct leaves the gene simulator can use
    tree = _fwd(z.BirthDeath(1.2, 0.3, mass_extinctions=[(2.5, 0.8)]), age=5.0, seed=4)
    assert _dead_names(tree)  # the tree has a dead part
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.3, loss=0.15,
                           origination=0.5, initial_families=30, max_family_size=0.5, seed=42)
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}


# --- ClaDS (per-lineage rates that shift at each speciation) -----------------

def test_clads_age_mode_complete_tree():
    t = _fwd(z.ClaDS(1.0, alpha=0.9, sigma=0.2, turnover=0.1), age=5.0, seed=1)
    assert abs(t.total_age - 5.0) < 1e-9
    assert t.root.time == 0.0
    assert len(t.extant_leaves()) >= 2
    assert all(len(n.children) == 2 for n in t.internal_nodes())  # binary
    for leaf in t.extant_leaves():
        assert abs(leaf.time - t.total_age) < 1e-9


def test_clads_n_tips_mode_hits_target():
    t = _fwd(z.ClaDS(1.0, turnover=0.1), n_tips=25, seed=2)
    assert len(t.extant_leaves()) == 25
    assert t.total_age > 0.0


def test_clads_turnover_zero_has_no_extinction():
    t = _fwd(z.ClaDS(1.0, turnover=0.0), age=4.0, seed=3)
    assert all(leaf.is_extant for leaf in t.leaves())  # pure birth with shifts


def test_clads_turnover_produces_extinct_leaves():
    n_dead = [sum(1 for leaf in _fwd(z.ClaDS(1.2, turnover=0.4), age=6.0, seed=s).leaves()
                  if not leaf.is_extant) for s in range(15)]
    assert np.mean(n_dead) > 0


def test_clads_reproducible():
    m = z.ClaDS(1.0, alpha=0.95, sigma=0.2, turnover=0.1)
    a = _fwd(m, age=5.0, seed=7).to_newick()
    b = _fwd(m, age=5.0, seed=7).to_newick()
    assert a == b


def test_clads_backward_is_rejected():
    with pytest.raises(ValueError):  # no closed-form reconstructed CDF -> forward-only
        z.simulate_species_tree(z.ClaDS(1.0), direction="backward", n_tips=10, age=5.0)


@pytest.mark.parametrize("kw", [dict(lambda_0=0.0), dict(lambda_0=1.0, alpha=0.0),
                                dict(lambda_0=1.0, sigma=-0.1), dict(lambda_0=1.0, turnover=1.0),
                                dict(lambda_0=1.0, sampling_fraction=1.5)])
def test_clads_validation(kw):
    with pytest.raises(ValueError):
        _fwd(z.ClaDS(**kw), age=5.0)


def test_clads_composes_with_mass_extinction():
    m = z.ClaDS(1.2, sigma=0.2, turnover=0.1, mass_extinctions=[(2.0, 0.8)])
    saw = False
    for s in range(30):
        t = _fwd(m, age=5.0, seed=s)   # pulse tree-time = 5 - 2 = 3
        if _pulse_deaths(t, 3.0):
            saw = True
            break
    assert saw


def test_clads_feeds_gene_sim():
    tree = _fwd(z.ClaDS(1.2, sigma=0.2, turnover=0.2), age=5.0, seed=4)
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.3, loss=0.15,
                           origination=0.5, initial_families=30, max_family_size=0.5, seed=42)
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}


# --- diversity-dependent (density-dependent) birth–death ---------------------

def test_diversity_dependent_saturates_at_K():
    # with μ=0 and a long age, the tree fills to exactly its carrying capacity K
    K = 20
    for s in range(6):
        t = _fwd(z.DiversityDependent(5.0, 0.0, carrying_capacity=K), age=30.0, seed=s)
        assert len(t.extant_leaves()) == K


def test_diversity_dependent_larger_K_more_tips():
    def mean_extant(K):
        return np.mean([len(_fwd(z.DiversityDependent(3.0, 0.2, carrying_capacity=K),
                                 age=25.0, seed=s).extant_leaves()) for s in range(15)])
    assert mean_extant(15) < mean_extant(60)


def test_diversity_dependent_n_tips_mode():
    t = _fwd(z.DiversityDependent(3.0, 0.3, carrying_capacity=50), n_tips=25, seed=3)
    assert len(t.extant_leaves()) == 25


def test_diversity_dependent_n_tips_above_K_rejected():
    with pytest.raises(ValueError):  # cannot grow past carrying capacity
        _fwd(z.DiversityDependent(3.0, 0.0, carrying_capacity=20), n_tips=40)


def test_diversity_dependent_reproducible():
    m = z.DiversityDependent(3.0, 0.3, carrying_capacity=40)
    a = _fwd(m, age=8.0, seed=7).to_newick()
    b = _fwd(m, age=8.0, seed=7).to_newick()
    assert a == b


def test_diversity_dependent_backward_is_rejected():
    with pytest.raises(ValueError):
        z.simulate_species_tree(z.DiversityDependent(2.0, carrying_capacity=30),
                                direction="backward", n_tips=10, age=5.0)


@pytest.mark.parametrize("kw", [dict(lambda_0=0.0, carrying_capacity=10),
                                dict(lambda_0=1.0, death=-0.1, carrying_capacity=10),
                                dict(lambda_0=1.0, carrying_capacity=0.0),
                                dict(lambda_0=1.0, carrying_capacity=10, sampling_fraction=0.0)])
def test_diversity_dependent_validation(kw):
    with pytest.raises(ValueError):
        _fwd(z.DiversityDependent(**kw), age=5.0)


def test_diversity_dependent_feeds_gene_sim():
    tree = _fwd(z.DiversityDependent(3.0, 0.3, carrying_capacity=40), age=12.0, seed=4)
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.3, loss=0.15,
                           origination=0.5, initial_families=30, max_family_size=0.5, seed=42)
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}


# --- clade-specific rate shifts ----------------------------------------------

def test_clade_shift_runs_and_is_binary():
    m = z.CladeShiftBirthDeath(0.8, 0.5, clade_shifts=[(3.5, 1.4, 0.2)])
    t = _fwd(m, age=4.0, seed=1)
    assert abs(t.total_age - 4.0) < 1e-9
    assert len(t.extant_leaves()) >= 2
    assert all(len(n.children) == 2 for n in t.internal_nodes())


def test_clade_shift_reproducible():
    m = z.CladeShiftBirthDeath(0.9, 0.4, clade_shifts=[(3.0, 1.5, 0.2), (2.0, 0.4, 0.6)])
    a = _fwd(m, age=4.0, seed=7).to_newick()
    b = _fwd(m, age=4.0, seed=7).to_newick()
    assert a == b


def test_clade_shift_to_fast_regime_raises_tip_count():
    # an early shift of one clade to a faster regime lifts the average number of extant tips
    base = z.BirthDeath(0.9, 0.5)
    fast = z.CladeShiftBirthDeath(0.9, 0.5, clade_shifts=[(3.5, 1.5, 0.2)])  # shift near the crown
    m_base = np.mean([len(_fwd(base, age=4.0, seed=s).extant_leaves()) for s in range(40)])
    m_fast = np.mean([len(_fwd(fast, age=4.0, seed=s).extant_leaves()) for s in range(40)])
    assert m_fast > m_base


def test_clade_shift_requires_age_mode():
    m = z.CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(2.0, 2.0, 0.1)])
    with pytest.raises(NotImplementedError):  # shift ages need a fixed present
        _fwd(m, n_tips=10)


def test_clade_shift_age_must_precede_crown():
    m = z.CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(5.0, 2.0, 0.1)])
    with pytest.raises(ValueError):
        _fwd(m, age=5.0, seed=1)


def test_clade_shift_backward_is_rejected():
    m = z.CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(2.0, 2.0, 0.1)])
    with pytest.raises(ValueError):
        z.simulate_species_tree(m, direction="backward", n_tips=10, age=5.0)


@pytest.mark.parametrize("kw", [dict(birth=0.0, clade_shifts=[(2.0, 1.0, 0.1)]),
                                dict(birth=1.0, clade_shifts=[]),
                                dict(birth=1.0, clade_shifts=[(2.0, 0.0, 0.1)]),   # shift birth <= 0
                                dict(birth=1.0, clade_shifts=[(0.0, 1.0, 0.1)])])  # shift age <= 0
def test_clade_shift_validation(kw):
    with pytest.raises(ValueError):
        _fwd(z.CladeShiftBirthDeath(**kw), age=5.0)


def test_clade_shift_composes_with_mass_extinction():
    m = z.CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(3.0, 1.8, 0.1)],
                               mass_extinctions=[(2.0, 0.7)])
    saw = False
    for s in range(30):
        t = _fwd(m, age=5.0, seed=s)   # pulse tree-time = 5 - 2 = 3
        if _pulse_deaths(t, 3.0):
            saw = True
            break
    assert saw


def test_clade_shift_feeds_gene_sim():
    tree = _fwd(z.CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(3.5, 1.6, 0.2)]), age=4.0, seed=4)
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.3, loss=0.15,
                           origination=0.5, initial_families=30, max_family_size=0.5, seed=42)
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}


def test_forward_tree_feeds_gene_sim_with_ghost_transfers():
    tree = _fwd(z.BirthDeath(1.0, 0.6), n_tips=40, seed=8)
    dead = _dead_names(tree)
    assert dead  # forward tree has a dead part
    g = z.simulate_genomes(tree, duplication=0.1, transfer=0.4, loss=0.15,
                           origination=0.5, initial_families=30, max_family_size=0.5, seed=42)
    # profiles only over extant species
    assert set(g.profiles.species) == {n.name for n in tree.extant_leaves()}
    # at least one transfer should involve a dead (extinct) branch — transfer from the dead,
    # for free, with no ghost-grafting step
    involved = any(
        (r.donor in dead or r.recipient in dead or r.branch in dead)
        for r in g.event_log if r.event is z.EventType.TRANSFER
    )
    assert involved
