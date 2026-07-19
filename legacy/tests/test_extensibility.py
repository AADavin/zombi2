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
    Rates,
    UnorderedGenome,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes import PerGenomeRates as RealPerGenomeRates
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
        tree, Rates(0.15, 0.1, 0.2, 0.4),
        np.random.default_rng(1), initial_size=8, genome_factory=OrderedListGenome,
    )
    assert gr.leaf_genomes
    for genome in gr.leaf_genomes.values():
        assert isinstance(genome, OrderedListGenome)
        assert len(genome.order) == genome.size()
        assert sorted(genome.order) == sorted(g.gid for g in genome.genes())


# --- Poisson oracle for the real (library) PerGenomeRates model --------------
def _total_branch_length(tree):
    return sum(n.branch_length() for n in tree.nodes_preorder() if n.parent is not None)


def _mean_event_count(tree, event, *, duplication=0.0, loss=0.0,
                      initial_size, reps, seed):
    """Mean number of ``event`` events over ``reps`` genome simulations on ``tree``."""
    gs = GenomeSimulator()
    master = np.random.default_rng(seed)
    counts = []
    for _ in range(reps):
        rng = np.random.default_rng(master.integers(1 << 62))
        gr = gs.simulate(
            tree,
            RealPerGenomeRates(duplication=duplication, transfer=0.0, loss=loss, origination=0.0),
            rng, initial_size=initial_size,
        )
        counts.append(sum(1 for r in gr.event_log if r.event is event))
    return np.array(counts, dtype=float)


def test_per_genome_rates_event_counts_match_poisson_oracle():
    """Realized D and L counts equal rate * total_branch_length (Poisson mean),
    and doubling/tripling the per-genome rate scales the mean count by the same factor.

    Under PerGenomeRates each event type fires at a *constant per-genome rate*,
    independent of genome size, so along the whole species tree (as long as the
    genome never empties) the number of duplication/loss events is
    Poisson(rate * total_branch_length). Duplication-only and loss-only are run
    separately, each seeded from a large, non-empty genome so the constant-rate
    assumption holds for the whole run.
    """
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    tbl = _total_branch_length(tree)

    reps = 3000

    # --- duplication-only oracle: mean == d * total_branch_length -------------
    d = 0.3
    dup = _mean_event_count(tree, EventType.DUPLICATION, duplication=d,
                            initial_size=50, reps=reps, seed=11)
    pred_d = d * tbl
    se_d = dup.std(ddof=1) / np.sqrt(reps)
    assert abs(dup.mean() - pred_d) < 5 * se_d, (dup.mean(), pred_d, se_d)

    # --- loss-only oracle: mean == l * total_branch_length --------------------
    l = 0.4
    loss = _mean_event_count(tree, EventType.LOSS, loss=l,
                             initial_size=200, reps=reps, seed=7)
    pred_l = l * tbl
    se_l = loss.std(ddof=1) / np.sqrt(reps)
    assert abs(loss.mean() - pred_l) < 5 * se_l, (loss.mean(), pred_l, se_l)

    # --- factor-r relationship: 3x the rate -> ~3x the mean count -------------
    r = 3.0
    dup_r = _mean_event_count(tree, EventType.DUPLICATION, duplication=d * r,
                              initial_size=50, reps=reps, seed=13)
    ratio = dup_r.mean() / dup.mean()
    assert abs(ratio - r) < 0.08 * r, (ratio, r, dup.mean(), dup_r.mean())
    assert abs(dup_r.mean() - r * pred_d) < 5 * (dup_r.std(ddof=1) / np.sqrt(reps))


# --- per-family sampled loss calibrated against a closed-form survival oracle -
def test_family_sampled_loss_calibrates_to_per_family_rate():
    # A per-family loss-only model is calibrated against a closed-form oracle: on an
    # ULTRAMETRIC extant tree every root->leaf path has the same total time T, so a
    # family that starts as one copy at the root and can only be *lost* survives to a
    # given extant leaf iff no loss fired on that length-T path. Loss is memoryless and
    # speciation reminting keeps copy number 1 (and the family's cached rate), so the
    # marginal presence probability at any extant leaf is exactly exp(-loss * T),
    # independent of topology. Families differ only through their sampled loss rate, so
    # binning by rate turns "sim tracks each family's rate" into a per-bin moment check
    # P_present = mean(exp(-l * T)). Fully-lost families are dropped from the profile, so
    # we score presence over ALL initial families (missing => present in zero leaves) to
    # avoid survivorship bias.
    T = 3.0
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=60, age=T, seed=7)
    extant = [n for n in tree.nodes_preorder() if n.is_leaf() and n.is_extant]
    assert {round(n.time, 9) for n in extant} == {T}   # oracle needs ultrametricity
    n_leaves = len(extant)

    fs = FamilySampledRates(duplication=0.0, transfer=0.0,
                            loss=Exponential(3.0), origination=0.0)
    n_families = 4000
    gr = simulate_genomes(tree, fs, initial_families=n_families, seed=11)

    # presence count per family across extant leaves, over EVERY initial family
    present_in = {fam: int(c) for fam, c in
                  zip(gr.profiles.families, (gr.profiles.matrix > 0).sum(axis=1))}
    all_families = [str(i) for i in range(1, n_families + 1)]
    present_count = np.array([present_in.get(fam, 0) for fam in all_families])
    loss_rate = np.array([fs.rates_for(fam)[2] for fam in all_families])  # cached draws

    edges = [0.0, 0.1, 0.2, 0.35, 0.6, 1.0, np.inf]
    checked = 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        idx = np.where((loss_rate >= lo) & (loss_rate < hi))[0]
        if len(idx) < 100:
            continue
        # oracle: mean per-leaf survival probability over the families in this bin
        expected = float(np.mean(np.exp(-loss_rate[idx] * T)))
        observed = present_count[idx].sum() / (len(idx) * n_leaves)
        # families are independent draws; bin mean has se ~ sqrt(p(1-p)/n_fam). A generous
        # multiple keeps it several-sigma non-flaky while still catching a real miscalibration.
        se = np.sqrt(max(expected * (1 - expected), 1e-6) / len(idx))
        assert abs(observed - expected) < 4 * se + 0.01, (
            f"bin [{lo},{hi}): observed {observed:.4f} vs oracle {expected:.4f} "
            f"(tol={4 * se + 0.01:.4f}, n_fam={len(idx)})")
        checked += 1
    assert checked >= 5  # several well-separated rate bins were actually validated

    # monotonicity: low-rate families are present in more leaves than high-rate ones
    order = np.argsort(loss_rate)
    low_third = present_count[order[: len(order) // 3]].mean()
    high_third = present_count[order[-len(order) // 3:]].mean()
    assert low_third > 3 * high_third
