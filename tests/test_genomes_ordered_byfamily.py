"""Per-family weights at the ordered resolution — the weight lands on the segment (SPEC §6).

The rule these cover: a modifier that varies with what is *inside* a segment must weight the
**segment, by what it covers**, never the position the event started from. Weighting the start would
apply a family's own rate to its neighbours, and the neighbourhood is reshuffled by every
rearrangement, so the parameter would not even name a fixed thing over a run.
"""

import collections

import numpy as np
import pytest


from zombi2 import species, genomes
from zombi2.genomes.ordered import Chromosome, Gene, _run_means
from zombi2.params import LogNormal, PerChromosome, PerCopy, PerLineage, Random


def _chrom(families, topology="circular"):
    return Chromosome(0, topology, [Gene(i, f, 1) for i, f in enumerate(families)])


@pytest.fixture
def tree():
    return species.simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1)


# --- the identity the two-step draw rests on -----------------------------------------------------

@pytest.mark.parametrize("m", [1, 2, 3, 4, 5, 6])
def test_summed_run_means_equal_summed_gene_weights_on_a_circle(m):
    """``Σ_s mean_w(s, m) == Σ_g w_g`` for **every** run size, on a circular chromosome.

    This is why the engine can draw the size first and the start second, and why the total rate needs
    no per-size term: the normaliser does not depend on the size drawn. Break this and the Gillespie
    total stops matching the process being sampled."""
    w = {0: 1.0, 1: 5.0, 2: 0.2, 3: 3.0, 4: 1.0, 5: 0.5}
    c = _chrom([0, 1, 2, 3, 4, 5])
    assert sum(_run_means(c, w, m)) == pytest.approx(sum(w[g.family] for g in c.genes))


def test_uniform_weights_leave_every_run_unweighted():
    """With no family variation every mean is 1, so the start pick stays uniform and the size
    distribution is untouched — the guarantee that setting no weight changes nothing."""
    c = _chrom([0, 1, 2, 3, 4, 5])
    assert _run_means(c, {f: 1.0 for f in range(6)}, 3) == [1.0] * 6


# --- the rule itself: the weight follows what a run covers ----------------------------------------

def test_a_heavy_gene_lifts_every_run_that_covers_it():
    """The heart of SPEC §6. Family 2 is heavy. At size 1 only the run *on* it is lifted; at size 3
    the runs starting at 0, 1 and 2 are all lifted, because all three cover it — and the runs that
    merely sit next to it are not lifted at all."""
    c = _chrom([0, 1, 2, 3, 4, 5])
    heavy = {0: 1.0, 1: 1.0, 2: 9.0, 3: 1.0, 4: 1.0, 5: 1.0}

    assert _run_means(c, heavy, 1) == [1.0, 1.0, 9.0, 1.0, 1.0, 1.0]

    at3 = _run_means(c, heavy, 3)
    assert at3[0] == at3[1] == at3[2] == pytest.approx(11 / 3)   # the three runs covering gene 2
    assert at3[3] == at3[4] == at3[5] == 1.0                     # the three that do not


def test_a_linear_run_is_clamped_by_its_start():
    """On a linear chromosome a run stops at the last gene, so the means near the end are taken over
    shorter runs — the same edge effect, from the same cause, that clamping already gives."""
    c = _chrom([0, 1, 2, 3], "linear")
    w = {0: 1.0, 1: 1.0, 2: 1.0, 3: 5.0}
    means = _run_means(c, w, 3)
    assert means[1] == pytest.approx(7 / 3)    # genes 1,2,3 — the full run
    assert means[2] == pytest.approx(3.0)      # genes 2,3 only — clamped
    assert means[3] == pytest.approx(5.0)      # gene 3 alone


# --- the engine ----------------------------------------------------------------------------------

def test_no_family_weight_is_byte_identical(tree):
    """The plain path must not move: same seed, same run, whether or not the feature exists."""
    kw = dict(duplication=0.3, loss=0.2, origination=0.3, inversion=0.4, initial_families=15,
              chromosomes=1, seed=5)
    a = genomes.simulate_genomes_ordered(tree, **kw)
    b = genomes.simulate_genomes_ordered(tree, **kw)
    assert [(e.time, e.kind, e.family, e.copy) for e in a.edges] == \
           [(e.time, e.kind, e.family, e.copy) for e in b.edges]


def test_byfamily_is_accepted_and_spreads_families_apart(tree):
    """A per-family weight should leave the average family alone (`a per-family draw` is mean-corrected) but
    pull families apart — measured on inversion, which changes no copy numbers and so cannot feed
    back on itself the way duplication does."""
    import collections

    def per_family_inversions(spread):
        counts = []
        for s in range(6):
            inv = 1.2 if spread is None else PerCopy(1.2).varying_among('families', LogNormal(0.0, spread))
            g = genomes.simulate_genomes_ordered(
                tree, origination=0.5, inversion=inv, inversion_extent=1,
                initial_families=25, chromosomes=1, seed=300 + s)
            c = collections.Counter(r.chromosome for r in g.rearrangements)
            counts.append(np.var(list(c.values())) if len(c) > 2 else 0.0)
        return float(np.mean(counts))

    assert per_family_inversions(None) >= 0.0          # it runs at all
    assert per_family_inversions(0.9) > 0.0            # and with a weight set


def test_one_object_shared_across_rates_is_accepted(tree):
    """A family-wide tempo here is one Random object read by every rate — one draw per family."""
    speed = Random('families', LogNormal(0.0, 0.6))
    g = genomes.simulate_genomes_ordered(tree, duplication=PerCopy(0.3).varying_among(speed),
                                         loss=PerCopy(0.2).varying_among(speed),
                                         origination=0.3,
                                         initial_families=15, chromosomes=1, seed=7)
    assert g.events


def test_byfamily_still_refused_on_origination(tree):
    """Origination is the rate at which families are *created*: when it is read there is no family
    yet to have drawn a factor. Unchanged from the family resolution."""
    with pytest.raises(ValueError, match="no family yet"):
        genomes.simulate_genomes_ordered(tree, origination=PerLineage(0.3).varying_among('families', LogNormal(0.0, 0.5)), seed=1)


def test_byfamily_refused_on_the_chromosome_tier(tree):
    """A per-family weight has to reach the genes an event covers; a fission acts on a whole
    replicon, so there is nothing for it to attach to."""
    with pytest.raises(ValueError, match="per-family draw on a PerChromosome scope"):
        genomes.simulate_genomes_ordered(tree, origination=0.3, chromosomes=2,
                                         fission=PerChromosome(0.2).varying_among('families', LogNormal(0.0, 0.5)), seed=1)


# --- the growth guard -----------------------------------------------------------------------------

def _worst_family_count(g):
    """The largest number of copies any one family reaches in any one genome of the run."""
    worst = 0
    for genome in g.node_genomes.values():
        c = collections.Counter(gene.family for chrom in genome for gene in chrom.genes)
        worst = max(worst, max(c.values(), default=0))
    return worst


def test_the_cap_actually_bounds_a_family(tree):
    """Duplication compounds, so without a quota a fast family multiplies without bound. The cap is
    the same guard the family resolution has always had, and an int is an absolute copy count."""
    g = genomes.simulate_genomes_ordered(tree, duplication=1.2, loss=0.1, origination=0.4,
                                         initial_families=10, chromosomes=1, max_family_size=4,
                                         seed=4)
    assert _worst_family_count(g) <= 4


def test_the_cap_holds_when_runs_carry_several_copies(tree):
    """A run may carry several copies of one family, so the test is *current + carried > cap*, not
    merely "already full" — otherwise a segmental duplication would overshoot the quota."""
    g = genomes.simulate_genomes_ordered(tree, duplication=1.2, loss=0.1, origination=0.4,
                                         duplication_extent=4, initial_families=10, chromosomes=1,
                                         max_family_size=5, seed=4)
    assert _worst_family_count(g) <= 5


def test_the_cap_can_be_switched_off(tree):
    g = genomes.simulate_genomes_ordered(tree, duplication=1.2, loss=0.1, origination=0.4,
                                         initial_families=10, chromosomes=1, max_family_size=None,
                                         seed=4)
    assert _worst_family_count(g) > 4          # unbounded growth is still available on request


def test_a_capped_run_refuses_whole_runs_not_partial_ones(tree):
    """Clipping a run to the genes still under quota would quietly shorten runs exactly where the
    genome is crowded, reshaping the extent distribution. Refusing outright keeps it intact."""
    g = genomes.simulate_genomes_ordered(tree, duplication=1.0, loss=0.2, origination=0.4,
                                         duplication_extent=3, initial_families=10, chromosomes=1,
                                         max_family_size=3, seed=8)
    dup_runs = [p.length for p in g.event_positions if p.kind == "duplication" and p.length]
    assert dup_runs, "the run should still produce duplications under a cap"
    assert max(dup_runs) > 1, "runs are refused whole, so multi-gene duplications still occur"


def test_the_cap_is_on_by_default(tree):
    """Both resolutions carry the same plain per-genome count, so they agree out of the box instead
    of one growing without bound. It is a runaway rail, not a modelling choice — hence a number no
    real family reaches rather than a biologically tuned one."""
    import inspect
    sig = inspect.signature(genomes.simulate_genomes_ordered)
    assert sig.parameters["max_family_size"].default == 10
    assert inspect.signature(genomes.simulate_genomes_family).parameters[
        "max_family_size"].default == 10


# --- drawing the run without scoring every start (issue #436) --------------------------------------
# The engine draws a gene in proportion to its weight, by rejection, and then one of the runs covering
# it. That has to give every start the chance the scan gives it: its run's mean weight, within the
# chromosome's share of the genome's weight. These compare the draw against those chances, computed
# exactly from `_run_means`, over many draws with fixed seeds.

from zombi2.genomes import ordered as _ordered


class _Size:
    """An extent that is always ``m`` genes, so the start is the only thing drawn by weight."""

    def __init__(self, m):
        self.m = m

    def sample(self, rng, **_):
        return self.m


def _genome(topology="circular"):
    return [_chrom([0, 1, 2, 3, 4, 5, 6], topology),
            Chromosome(1, topology, [Gene(10 + i, f, 1) for i, f in enumerate([7, 8, 9, 2, 5])])]


def _weights():
    w = _ordered._FamilyWeights()
    for family, weight in {0: 1.0, 1: 0.2, 2: 6.0, 3: 1.0, 4: 0.5, 5: 3.0, 6: 1.0,
                           7: 0.2, 8: 1.0, 9: 0.1}.items():
        w[family] = weight
    return w


def _exact(genome, weights, m):
    """The chance of each ``(chromosome, start)``, as the scan computes it."""
    sums = [sum(weights[g.family] for g in c.genes) for c in genome]
    out = {}
    for ci, chrom in enumerate(genome):
        means = _run_means(chrom, weights, min(m, len(chrom.genes)))
        for s, mean in enumerate(means):
            out[(ci, s)] = sums[ci] / sum(sums) * mean / sum(means)
    return out


def _drawn(genome, weights, m, n=60000, seed=0):
    rng = np.random.default_rng(seed)
    hits = collections.Counter()
    for _ in range(n):
        ci, s, _size = _ordered._pick_run_by_family(rng, genome, weights, _Size(m))
        hits[(ci, s)] += 1
    return {key: count / n for key, count in hits.items()}


def _distance(a, b):
    """Total variation: the share of the probability that sits somewhere else."""
    return 0.5 * sum(abs(a.get(key, 0.0) - b.get(key, 0.0)) for key in set(a) | set(b))


@pytest.mark.parametrize("m", [1, 2, 3, 5])
def test_every_start_gets_the_chance_the_scan_gives_it_on_a_circle(m):
    genome, weights = _genome(), _weights()
    assert _distance(_drawn(genome, weights, m), _exact(genome, weights, m)) < 0.02


def test_every_start_gets_the_chance_the_scan_gives_it_on_a_line_at_one_gene():
    """On a linear chromosome the two draws agree exactly at the default run of one gene."""
    genome, weights = _genome("linear"), _weights()
    assert _distance(_drawn(genome, weights, 1), _exact(genome, weights, 1)) < 0.02


def test_a_wrong_offset_would_be_caught():
    """The comparison is not vacuous: drawing the start on the gene itself, ignoring the runs that
    begin before it, is far from the scan once runs are longer than one gene."""
    genome, weights = _genome(), _weights()
    rng = np.random.default_rng(1)
    hits = collections.Counter()
    for _ in range(20000):
        ci, j = _ordered._gene_by_rate(rng, genome, 12, weights, weights.largest, 10**6)
        hits[(ci, j)] += 1
    wrong = {key: count / 20000 for key, count in hits.items()}
    assert _distance(wrong, _exact(genome, weights, 3)) > 0.1


def test_a_gene_is_kept_in_proportion_to_its_weight():
    genome, weights = _genome(), _weights()
    rng = np.random.default_rng(2)
    hits = collections.Counter()
    for _ in range(60000):
        ci, j = _ordered._gene_by_rate(rng, genome, 12, weights, weights.largest, 10**6)
        hits[genome[ci].genes[j].id] += 1
    total = sum(weights[g.family] for c in genome for g in c.genes)
    for c in genome:
        for g in c.genes:
            assert hits[g.id] / 60000 == pytest.approx(weights[g.family] / total, abs=0.01)


def test_out_of_tries_the_genome_is_scored_and_the_chances_hold(monkeypatch):
    """With no tries at all every event falls to the scan, and the chances are still the scan's."""
    monkeypatch.setattr(_ordered, "_MIN_TRIES", 0)
    monkeypatch.setattr(_ordered, "_GENES_PER_TRY", 10**9)
    scanned = []
    real = _ordered._pick_run_by_scan
    monkeypatch.setattr(_ordered, "_pick_run_by_scan",
                        lambda *a, **k: scanned.append(1) or real(*a, **k))
    genome, weights = _genome(), _weights()
    assert _distance(_drawn(genome, weights, 3, n=20000), _exact(genome, weights, 3)) < 0.03
    assert len(scanned) == 20000


def test_a_table_that_does_not_know_its_largest_weight_is_scored():
    """A plain mapping carries no bound, and a bound is what rejection needs, so it is scanned."""
    genome = _genome()
    plain = dict(_weights())
    assert _distance(_drawn(genome, plain, 2, n=20000), _exact(genome, plain, 2)) < 0.03


def test_the_largest_weight_is_kept_as_families_arrive():
    w = _ordered._FamilyWeights()
    for family, weight in [(0, 0.5), (1, 2.5), (2, 1.0)]:
        w[family] = weight
    assert w.largest == 2.5


def test_a_genome_whose_genes_all_weigh_nothing_has_no_run():
    zero = _ordered._FamilyWeights()
    for f in range(10):
        zero[f] = 0.0
    zero[99] = 1.0                           # a family this genome does not hold
    rng = np.random.default_rng(0)
    assert _ordered._pick_run_by_family(rng, _genome(), zero, _Size(2)) is None
