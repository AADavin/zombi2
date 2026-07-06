"""Tests for the gene-family coupling model (:mod:`zombi2.coupling`).

Three layers:

1. **Unit** — the coupling spec constructors and the exact loss-rate formula
   ``base_loss·exp(-β·f_i)`` with ``f_i = h_i + Σ_j J_ij σ_j`` (partners only).
2. **Ground-truth recovery** — inject a known ``J`` and confirm the generated profiles show
   the prescribed co-occurrence structure: positive ``J`` → co-occurrence, negative ``J`` →
   avoidance, zero ``J`` → no structure, and uncoupled families stay uncorrelated.
3. **Driver** — :func:`simulate_coupled` shape/panel/reproducibility contract.

The recovery tests run on a *near-star* tree (all lineages split just below the root, then
evolve independently for the whole age) which isolates the injected coupling from the
phylogenetic confounding that inflates co-occurrence on a normal birth–death tree — the very
"shared ancestry" trap the design note (``docs/non_independence.tex``) warns about, and the
reason real inference (Fukunaga & Iwasaki 2022) corrects for the tree.
"""

import itertools
import math

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.coevolve.coupling import (
    CouplingSpec,
    PottsRates,
    pathway_blocks,
    simulate_coupled,
)
from zombi2.genomes.events import EventType
from zombi2.genomes.genome import Gene, IdManager, UnorderedGenome
from zombi2.tree import Tree, TreeNode


# ── fixtures ─────────────────────────────────────────────────────────────────
def near_star_tree(k: int, age: float, delta: float = 1e-3) -> Tree:
    """A balanced binary tree of ``2**k`` tips whose internal nodes all sit within
    ``[0, delta·k]`` of the root — so every terminal branch is ~``age`` long and the tips
    are near-independent draws (minimal shared ancestry)."""
    counter = itertools.count()
    root = TreeNode(name="root", time=0.0)

    def split(node, level):
        if level == k:
            return
        for _ in range(2):
            child = TreeNode(name=f"i{next(counter)}", time=min(age, delta * (level + 1)))
            node.add_child(child)
            split(child, level + 1)

    split(root, 0)
    stack, leaves = [root], []
    while stack:
        n = stack.pop()
        leaves.append(n) if not n.children else stack.extend(n.children)
    for idx, leaf in enumerate(leaves):
        leaf.time, leaf.name, leaf.is_extant = age, f"n{idx}", True
    return Tree(root, age)


def _corr(a, b) -> float:
    a, b = a.astype(float), b.astype(float)
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _make_genome(families):
    ids = IdManager()
    g = UnorderedGenome(ids)
    for fam in families:
        g._add(Gene(ids.new_gene(), fam))
    return g


# canonical recovery-test regime (see the parameter sweep in the design work): moderate
# occupancy (~0.78) so both co-occurrence and avoidance register.
_REGIME = dict(base_loss=1.0, transfer=0.2, beta=1.0, h=2.0)


# ── 1. unit: spec constructors ───────────────────────────────────────────────
def test_from_dense_and_from_edges_agree():
    J = np.array([[0.0, 1.2, -0.4],
                  [1.2, 0.0, 0.0],
                  [-0.4, 0.0, 0.0]])
    a = CouplingSpec.from_dense(J)
    b = CouplingSpec.from_edges(3, {(0, 1): 1.2, (0, 2): -0.4})
    assert np.allclose(a.dense_J(), J)
    assert np.allclose(b.dense_J(), J)          # edges are symmetrised
    assert np.allclose(a.dense_J(), b.dense_J())


def test_from_edges_rejects_self_and_out_of_range():
    with pytest.raises(ValueError):
        CouplingSpec.from_edges(3, {(1, 1): 1.0})          # self-coupling
    with pytest.raises(ValueError):
        CouplingSpec.from_edges(3, {(0, 5): 1.0})          # out of range


def test_pathway_blocks_structure():
    spec = pathway_blocks([2, 2], within=3.0, between=-1.0)
    J = spec.dense_J()
    assert np.all(np.diag(J) == 0.0)                       # zero diagonal
    assert J[0, 1] == 3.0 and J[2, 3] == 3.0               # within-block
    assert J[0, 2] == -1.0 and J[1, 3] == -1.0             # between-block
    assert np.allclose(J, J.T)                             # symmetric
    assert spec.panel_ids == ["F0", "F1", "F2", "F3"]


def test_bad_spec_rejected():
    with pytest.raises(ValueError):
        CouplingSpec.from_dense(np.zeros((3, 3)), base_loss=-1.0)
    with pytest.raises(ValueError):
        CouplingSpec.from_dense(np.zeros((3, 3)), h=np.zeros(2))   # wrong h length


# ── 1. unit: the loss-rate formula ───────────────────────────────────────────
def test_loss_rate_matches_field_formula():
    spec = CouplingSpec.from_edges(
        3, {(0, 1): 1.2, (0, 2): -0.4},
        h=[0.5, -0.3, 0.0], base_loss=1.0, beta=0.5, transfer=0.7,
    )
    rates = PottsRates(spec)
    g = _make_genome(["F0", "F1", "F2"])            # all three present
    ws = rates.event_weights(g, "b", 0.0)
    loss = {e.family: e.rate for e in ws if e.event is EventType.LOSS}

    # f0 = 0.5 + J01·1 + J02·1 = 0.5 + 1.2 - 0.4 = 1.3
    assert loss["F0"] == pytest.approx(math.exp(-0.5 * 1.3))
    # f1 = -0.3 + J01·1 = 0.9
    assert loss["F1"] == pytest.approx(math.exp(-0.5 * 0.9))
    # f2 =  0.0 + J02·1 = -0.4   → field < 0 raises loss above base
    assert loss["F2"] == pytest.approx(math.exp(-0.5 * -0.4))
    assert loss["F2"] > spec.base_loss

    # transfer is the field-blind gain channel: per-copy rate × genome size
    trans = [e for e in ws if e.event is EventType.TRANSFER]
    assert len(trans) == 1 and trans[0].family is None
    assert trans[0].rate == pytest.approx(0.7 * g.size())


def test_present_partner_changes_the_field():
    spec = CouplingSpec.from_edges(3, {(0, 1): 1.2, (0, 2): -0.4},
                                   h=[0.5, -0.3, 0.0], base_loss=1.0, beta=0.5)
    rates = PottsRates(spec)
    # drop F2: now f0 = 0.5 + 1.2 = 1.7 (only the F1 partner remains)
    g = _make_genome(["F0", "F1"])
    loss = {e.family: e.rate for e in rates.event_weights(g, "b", 0.0)
            if e.event is EventType.LOSS}
    assert loss["F0"] == pytest.approx(math.exp(-0.5 * 1.7))


def test_positive_partner_protects_negative_partner_exposes():
    spec = CouplingSpec.from_edges(3, {(0, 1): 2.0, (0, 2): -2.0}, base_loss=1.0, beta=1.0)
    rates = PottsRates(spec)
    solo = {e.family: e.rate for e in rates.event_weights(_make_genome(["F0"]), "b", 0)
            if e.event is EventType.LOSS}["F0"]
    with_pos = {e.family: e.rate for e in rates.event_weights(_make_genome(["F0", "F1"]), "b", 0)
                if e.event is EventType.LOSS}["F0"]
    with_neg = {e.family: e.rate for e in rates.event_weights(_make_genome(["F0", "F2"]), "b", 0)
                if e.event is EventType.LOSS}["F0"]
    assert with_pos < solo < with_neg                     # protection vs exposure
    assert solo == pytest.approx(1.0)                     # field 0 → base_loss


def test_non_panel_family_is_uncoupled():
    spec = CouplingSpec.from_edges(2, {(0, 1): 3.0}, base_loss=0.7, beta=1.0)
    rates = PottsRates(spec)
    g = _make_genome(["F0", "F1", "X"])                   # X is not in the panel
    loss = {e.family: e.rate for e in rates.event_weights(g, "b", 0.0)
            if e.event is EventType.LOSS}
    assert loss["X"] == pytest.approx(0.7)               # base_loss, no field


def test_origination_channel_optional():
    spec = CouplingSpec.from_edges(2, {(0, 1): 1.0}, origination=0.3)
    ws = PottsRates(spec).event_weights(_make_genome(["F0"]), "b", 0.0)
    orig = [e for e in ws if e.event is EventType.ORIGINATION]
    assert len(orig) == 1 and orig[0].rate == pytest.approx(0.3)
    # default: no origination channel (closed panel)
    ws0 = PottsRates(CouplingSpec.from_edges(2, {(0, 1): 1.0})).event_weights(
        _make_genome(["F0"]), "b", 0.0)
    assert not [e for e in ws0 if e.event is EventType.ORIGINATION]


# ── 2. ground-truth recovery ─────────────────────────────────────────────────
def test_positive_coupling_creates_cooccurrence():
    """Injected +J → coupled pair co-occurs, while an uncoupled pair does not."""
    tree = near_star_tree(8, age=6.0)
    spec = pathway_blocks([2, 1, 1], within=3.0, between=0.0, **_REGIME)  # (0,1) coupled
    P = simulate_coupled(tree, spec, seed=3).profiles.presence()
    coupled = _corr(P[0], P[1])
    uncoupled = _corr(P[2], P[3])
    assert coupled > 0.25
    assert coupled > uncoupled + 0.2


def test_zero_coupling_is_null():
    """J = 0 → the same panel positions show no injected co-occurrence."""
    tree = near_star_tree(8, age=6.0)
    spec = pathway_blocks([2, 1, 1], within=0.0, between=0.0, **_REGIME)
    P = simulate_coupled(tree, spec, seed=3).profiles.presence()
    assert abs(_corr(P[0], P[1])) < 0.25
    assert abs(_corr(P[2], P[3])) < 0.25


def test_negative_coupling_creates_avoidance():
    """Injected -J → coupled pair avoids each other (anti-correlated presence)."""
    tree = near_star_tree(8, age=6.0)
    spec = pathway_blocks([2, 1, 1], within=-3.0, between=0.0, **_REGIME)
    P = simulate_coupled(tree, spec, seed=3).profiles.presence()
    assert _corr(P[0], P[1]) < -0.25


def test_coupling_beats_no_coupling_on_the_same_pair():
    """Phylogeny-controlled: on one tree, the coupled pair's correlation exceeds the same
    pair's correlation under J = 0 (differencing out the shared-ancestry baseline)."""
    tree = near_star_tree(8, age=6.0)
    coupled = simulate_coupled(
        tree, pathway_blocks([2, 1, 1], within=3.0, **_REGIME), seed=3).profiles.presence()
    null = simulate_coupled(
        tree, pathway_blocks([2, 1, 1], within=0.0, **_REGIME), seed=3).profiles.presence()
    assert _corr(coupled[0], coupled[1]) > _corr(null[0], null[1]) + 0.2


# ── 3. driver contract ───────────────────────────────────────────────────────
def test_profile_shape_and_panel_rows():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=25, age=4.0, seed=5)
    spec = pathway_blocks([3, 2], within=2.0, **_REGIME)      # 5-family panel
    res = simulate_coupled(tree, spec, seed=1)
    # every panel family is a row (even any that went globally extinct), species = tips
    assert res.profiles.families == spec.panel_ids
    assert res.profiles.shape == (5, 25)
    assert res.profiles.shape[1] == len(tree.extant_leaves())


def test_reproducible_with_seed():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=20, age=4.0, seed=2)
    spec = pathway_blocks([2, 2], within=2.0, **_REGIME)
    a = simulate_coupled(tree, spec, seed=7).profiles
    b = simulate_coupled(tree, spec, seed=7).profiles
    assert np.array_equal(a.matrix, b.matrix)


def test_initial_presence_mask():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=15, age=1.0, seed=9)
    spec = pathway_blocks([2, 2], within=2.0, base_loss=0.0, transfer=0.0, h=0.0)
    # no loss, no transfer, tiny age → the initial presence pattern is preserved at the tips
    mask = np.array([1, 0, 1, 0])
    res = simulate_coupled(tree, spec, seed=1, initial_presence=mask)
    present_rows = set(np.unique(res.profiles.coo[0]))
    assert present_rows == {0, 2}                            # only the seeded families


def test_runs_on_realistic_tree():
    """Smoke test on an ordinary birth–death tree (with phylogenetic structure)."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=40, age=5.0, seed=4)
    spec = pathway_blocks([4, 4], within=2.5, between=-1.0, **_REGIME)
    res = simulate_coupled(tree, spec, seed=2)
    assert res.profiles.shape == (8, 40)
    assert len(res.event_log) > 0


# ── 3b. per-family internal-node origin seeding (``origins=``) ────────────────
def two_clade_tree(age: float = 1.0) -> Tree:
    """``root → {cladeA → (a0,a1)}, {cladeB → (b0,b1)}`` — cladeA's subtree is exactly
    ``{a0, a1}``, so a family born at cladeA must appear only there when nothing erases it."""
    root = TreeNode(name="root", time=0.0)
    a = TreeNode(name="cladeA", time=0.5 * age)
    b = TreeNode(name="cladeB", time=0.5 * age)
    root.add_child(a)
    root.add_child(b)
    for parent, names in ((a, ("a0", "a1")), (b, ("b0", "b1"))):
        for nm in names:
            leaf = TreeNode(name=nm, time=age, is_extant=True)
            parent.add_child(leaf)
    return Tree(root, age)


def _presence_by_leaf(res, family):
    """Leaf names where ``family`` is present (queries the genomes, not the natkey column order)."""
    return {node.name for node, g in res.leaf_genomes.items() if g.copy_number(family) > 0}


_FROZEN = dict(base_loss=0.0, transfer=0.0, h=0.0)  # no loss, no gain → seeds are frozen at the tips


def test_origin_at_internal_node_confines_family_to_that_clade():
    tree = two_clade_tree()
    spec = pathway_blocks([2, 2], within=2.0, **_FROZEN)  # panel F0..F3
    node_a = next(n for n in tree.nodes_preorder() if n.name == "cladeA")
    res = simulate_coupled(tree, spec, seed=1, origins={"F0": node_a})  # F0 born at cladeA
    assert _presence_by_leaf(res, "F0") == {"a0", "a1"}                 # confined to the clade
    for fam in ("F1", "F2", "F3"):                                      # the rest default to root
        assert _presence_by_leaf(res, fam) == {"a0", "a1", "b0", "b1"}


def test_origin_logs_origination_at_the_birth_node():
    tree = two_clade_tree()
    spec = pathway_blocks([2, 2], within=2.0, **_FROZEN)
    res = simulate_coupled(tree, spec, seed=1, origins={"F0": "cladeA"})  # by node name
    origins = [(r.branch, r.family) for r in res.event_log.records
               if r.event is EventType.ORIGINATION]
    assert ("cladeA", "F0") in origins           # F0's birth is recorded at cladeA
    assert ("root", "F0") not in origins         # and not at the root


def test_origins_all_at_root_matches_default():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=15, age=2.0, seed=9)
    spec = pathway_blocks([2, 2], within=2.0, **_REGIME)
    default = simulate_coupled(tree, spec, seed=4).profiles
    at_root = simulate_coupled(tree, spec, seed=4,
                               origins={f: "root" for f in spec.panel_ids}).profiles
    assert np.array_equal(default.matrix, at_root.matrix)


def test_origin_none_keeps_family_absent_everywhere():
    tree = two_clade_tree()
    spec = pathway_blocks([2, 2], within=2.0, **_FROZEN)
    res = simulate_coupled(tree, spec, seed=1, origins={"F0": None})
    assert _presence_by_leaf(res, "F0") == set()          # never seeded
    assert _presence_by_leaf(res, "F1") == {"a0", "a1", "b0", "b1"}


def test_origins_unknown_node_rejected():
    tree = two_clade_tree()
    spec = pathway_blocks([2, 2], within=2.0, **_FROZEN)
    with pytest.raises(ValueError, match="unknown node name"):
        simulate_coupled(tree, spec, seed=1, origins={"F0": "no_such_node"})


def test_origins_non_panel_family_rejected():
    tree = two_clade_tree()
    spec = pathway_blocks([2, 2], within=2.0, **_FROZEN)
    with pytest.raises(ValueError, match="non-panel families"):
        simulate_coupled(tree, spec, seed=1, origins={"F9": "cladeA"})


def test_origins_and_initial_presence_mutually_exclusive():
    tree = two_clade_tree()
    spec = pathway_blocks([2, 2], within=2.0, **_FROZEN)
    with pytest.raises(ValueError, match="not both"):
        simulate_coupled(tree, spec, seed=1,
                         initial_presence=np.ones(4), origins={"F0": "cladeA"})


# ── 4. gain-side coupling: field-biased HGT establishment (option b) ──────────
def _sel(family):
    """A minimal stand-in for a transfer selection: the hook reads ``selection.genes[0].family``."""
    from types import SimpleNamespace
    return SimpleNamespace(genes=[SimpleNamespace(family=family)])


def test_establishment_probability_formula():
    """p = exp(-β·g·(f_maxᵢ - f_i)); deficit measured from the best-possible field."""
    spec = CouplingSpec.from_edges(3, {(0, 1): 2.0, (0, 2): -1.0},
                                   h=[0.5, 0.0, 0.0], beta=0.5, gain_coupling=2.0)
    rates = PottsRates(spec)
    # f_max[0] = h0 + max(J01,0) + max(J02,0) = 0.5 + 2 + 0 = 2.5
    assert spec.f_max[0] == pytest.approx(2.5)
    # only the positive partner present → field == f_max → establishes freely (p = 1)
    assert rates.establishment_probability(_sel("F0"), _make_genome(["F1"]), 0.0) == pytest.approx(1.0)
    # both partners present → f0 = 0.5+2-1 = 1.5, deficit 1.0 → exp(-0.5·2·1.0)
    assert rates.establishment_probability(_sel("F0"), _make_genome(["F1", "F2"]), 0.0) \
        == pytest.approx(math.exp(-1.0))
    # only the negative partner → f0 = -0.5, deficit 3.0 → exp(-0.5·2·3.0)
    assert rates.establishment_probability(_sel("F0"), _make_genome(["F2"]), 0.0) \
        == pytest.approx(math.exp(-3.0))
    # non-panel family is uncoupled → p = 1
    assert rates.establishment_probability(_sel("X"), _make_genome(["F1"]), 0.0) == 1.0


def test_establishment_probability_off_by_default():
    """gain_coupling = 0 (the default) → every transfer establishes (p = 1), field ignored."""
    spec = CouplingSpec.from_edges(3, {(0, 1): 2.0, (0, 2): -1.0}, h=[0.5, 0.0, 0.0], beta=0.5)
    assert spec.gain_coupling == 0.0
    rates = PottsRates(spec)
    for g in ([], ["F1"], ["F2"], ["F1", "F2"]):
        assert rates.establishment_probability(_sel("F0"), _make_genome(g), 0.0) == 1.0


def test_gain_coupling_inert_without_transfers():
    """With no transfers there is nothing to gate: gain_coupling leaves output byte-identical
    (and consumes no RNG draw), confirming gain stays purely donor-limited."""
    tree = near_star_tree(6, age=3.0)
    common = dict(within=3.0, between=0.0, base_loss=1.0, transfer=0.0, beta=1.0, h=2.0)
    a = simulate_coupled(tree, pathway_blocks([2, 2], **common, gain_coupling=0.0), seed=5).profiles
    b = simulate_coupled(tree, pathway_blocks([2, 2], **common, gain_coupling=5.0), seed=5).profiles
    assert np.array_equal(a.matrix, b.matrix)


def test_gain_coupling_reproducible():
    """The extra establishment draw is deterministic under a fixed seed."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=20, age=4.0, seed=2)
    spec = pathway_blocks([2, 2], within=2.0, gain_coupling=3.0, **_REGIME)
    a = simulate_coupled(tree, spec, seed=7).profiles
    b = simulate_coupled(tree, spec, seed=7).profiles
    assert np.array_equal(a.matrix, b.matrix)


def test_gain_coupling_strengthens_cooccurrence():
    """Field-biased establishment adds signal on the *gain* side: on the same tree/seed, the
    coupled pair co-occurs more strongly with gain_coupling > 0 than with the loss coupling
    alone (differential retention of transferred copies reinforces co-occurrence)."""
    tree = near_star_tree(8, age=6.0)
    base = dict(within=2.0, between=0.0, **_REGIME)
    off = simulate_coupled(
        tree, pathway_blocks([2, 1, 1], **base, gain_coupling=0.0), seed=3).profiles.presence()
    on = simulate_coupled(
        tree, pathway_blocks([2, 1, 1], **base, gain_coupling=4.0), seed=3).profiles.presence()
    assert _corr(on[0], on[1]) > _corr(off[0], off[1])
    assert _corr(on[2], on[3]) < 0.25          # the uncoupled pair stays unstructured


def test_coupled_result_reconstructs_gene_trees():
    """The coupled event log carries a root origination for each panel family, so the standard
    Genomes machinery reconstructs a gene tree per family — the fix that lets --rate-model
    coupled write gene_trees/ like every other model."""
    from zombi2.genomes.simulation import Genomes
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=15, age=3.0, seed=1)
    spec = pathway_blocks([3, 3], within=2.0, **_REGIME)
    res = simulate_coupled(tree, spec, seed=1)

    fam_events = res.event_log.by_family()
    assert set(spec.panel_ids) <= set(fam_events)
    for fam in spec.panel_ids:                  # each panel family has its root birth record
        assert any(r.event is EventType.ORIGINATION for r in fam_events[fam])

    genomes = Genomes(species_tree=tree, leaf_genomes=res.leaf_genomes,
                      event_log=res.event_log, profiles=res.profiles)
    trees = genomes.gene_trees()
    assert set(trees) == set(spec.panel_ids)
    assert any(complete for complete, _extant in trees.values())  # non-empty genealogies


# --- ABC for the coupled model: co-occurrence summary + match_coupled -------------

from zombi2 import ProfileMatrix                                             # noqa: E402
from zombi2.matching import (                                               # noqa: E402
    cooccurrence_features, cooccurrence_summary, match_coupled,
)

_ABC_COMMON = dict(between=0.0, h=-0.5, base_loss=1.0, transfer=0.4)


def test_cooccurrence_features_detect_modules():
    """At matched marginal prevalence, coupled families (modules) close far more triangles
    than independent ones. (Prevalence must be matched — the raw triangle count is otherwise
    confounded, since sparse independent families also produce spurious co-occurrence.)"""
    tree = simulate_species_tree(BirthDeath(1.0, 0.5), n_tips=45, age=1.0, seed=1)
    # h calibrated so both sit at prevalence ~0.5 (coupled needs a lower field to offset the
    # co-retention its coupling induces): independent h=+0.2, coupled(J=1.2) h=-2.4.
    coupled = simulate_coupled(
        tree, pathway_blocks([4]*15, within=1.2, between=0.0, h=-2.4, base_loss=1.0, transfer=0.4), seed=1)
    indep = simulate_coupled(
        tree, pathway_blocks([4]*15, within=0.0, between=0.0, h=+0.2, base_loss=1.0, transfer=0.4), seed=1)
    assert abs(coupled.profiles.matrix.mean() - indep.profiles.matrix.mean()) < 0.12  # prevalence matched
    fc = cooccurrence_features(coupled.profiles)
    fi = cooccurrence_features(indep.profiles)
    assert fc.shape == (3,)                                   # [edges, triangles, transitivity]
    assert fc[1] > 2 * fi[1]                                  # coupling closes many more triangles


def test_cooccurrence_features_constant_is_zero():
    """Families with no presence variance -> zero features, no crash / no nan."""
    pm = ProfileMatrix(["a", "b", "c"], ["s1", "s2"], np.ones((3, 2), dtype=int))
    assert np.array_equal(cooccurrence_features(pm), np.zeros(3))


def test_cooccurrence_summary_structure():
    """Summary = default marginal block + 3 co-occurrence features, with balancing weights."""
    species = [f"s{i}" for i in range(8)]
    summ = cooccurrence_summary(species, max_copies=4)
    pm = ProfileMatrix(["f1", "f2"], species, np.random.default_rng(0).integers(0, 2, (2, 8)))
    v = summ(pm)
    assert len(v) == 2 * 8 + 4 + 3                            # freq + sizes + copyspec + 3 co-occ
    assert len(summ.feature_weights) == len(v)
    assert summ.feature_weights[-1] > 1.0                     # small co-occ block is up-weighted


def test_match_coupled_recovers_coupling_strength():
    """match_coupled fits the within-pathway coupling J: a coupled target yields a clearly
    higher fitted J than an uncoupled one, and returns a well-formed posterior."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.5), n_tips=40, age=1.0, seed=1)

    def builder(p):
        return pathway_blocks([4]*12, within=p["within"], **_ABC_COMMON)

    def fit_median(true_J):
        emp = simulate_coupled(tree, builder({"within": true_J}), seed=777).profiles
        fit = match_coupled(tree, emp, builder, {"within": (0.0, 1.5)},
                            n_sims=120, accept=0.1, seed=1)
        assert fit.param_names == ["within"]
        assert len(fit.posterior["within"]) >= 1
        return fit.summary()["within"]["median"]

    assert fit_median(1.0) > fit_median(0.0) + 0.2           # coupling strength is identified


def test_match_coupled_threads_origins_through():
    """origins= reaches the simulator: fitting the same target with families born at an internal
    clade vs all-at-root explores different simulations (different accepted summaries)."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.5), n_tips=24, age=1.0, seed=1)
    spec = pathway_blocks([3, 3], within=1.0, **_ABC_COMMON)
    node = next(n for n in tree.internal_nodes() if n is not tree.root)
    origins = {f: node for f in spec.panel_ids}
    emp = simulate_coupled(tree, spec, seed=5, origins=origins).profiles

    def builder(p):
        return pathway_blocks([3, 3], within=p["within"], **_ABC_COMMON)

    common = dict(n_sims=25, accept=0.2, seed=1)
    at_root = match_coupled(tree, emp, builder, {"within": (0.0, 1.5)}, **common)
    at_node = match_coupled(tree, emp, builder, {"within": (0.0, 1.5)}, origins=origins, **common)
    # same empirical target, but the origins run simulates a different (clade-restricted) process
    assert np.array_equal(at_root.empirical_summary, at_node.empirical_summary)
    assert not np.array_equal(at_root.accepted_summaries, at_node.accepted_summaries)
