"""Validation — the simulator against theory, not against itself.

Every other test in this suite asks whether the code does what the code intends. These ask something
different and harder: whether the process ZOMBI2 samples from is the process it claims. A simulator
can be internally consistent, deterministic, well covered and still be sampling the wrong
distribution — and nothing that compares the code to itself will ever notice.

So each test here computes a quantity with a **closed form** and checks the run against it:

- a Yule tree's size is geometric with mean ``e^(λT)``;
- a birth–death run's mean size is ``e^((λ-μ)T)`` and its extinction probability is known exactly;
- ``n_extant`` conditioning is **length-weighted**, which is what makes the terminal branches
  non-degenerate — the naive "stop at the Nth birth" scheme biases every downstream estimate, and is
  ruled out here rather than assumed against;
- origination is counted **per lineage**, so its expected count is the rate times the tree's total
  branch length;
- duplication, transfer and loss are counted **per copy**, so by the compensator identity their
  expected counts are the rate times the total gene-copy-time — one test that pins the units *and*
  the scope of all three at once;
- inversion breakpoints fall **uniformly** around a chromosome, which on a ring is the statement that
  the arcs they cut are ``Dirichlet(1, …, 1)``, and their extent has the geometric *shape* its mean
  implies rather than merely the right average.

**Every test is deterministic.** The seeds are fixed, so a "statistical" test here cannot flake: it
computes one number and compares it to one expectation. The tolerances are wide (``|z| < 4``) because
their job is to catch a model that has changed, not to detect a rate that is off by a percent — and
the hypotheses these rule out are rejected by tens of standard deviations, not by twos.

What that buys, measured rather than assumed: scaling a rate by 1.15 behind the engine's back fails
these tests, and scaling it by 1.02 does not. Resolution is around **ten percent**. So this file is
evidence that the model is the one advertised, not a calibration of it — a rate quietly off by a few
percent would still get past, and catching that needs more replicates than a test suite should spend.

If one of these fails, the arithmetic below is not what to doubt first. Something about what the
engine samples has changed.
"""
from __future__ import annotations

import math

import numpy as np

from zombi2.genomes import simulate_genomes_family, simulate_genomes_ordered
from zombi2.genomes.ordered import Inversion
from zombi2.species import simulate_species_tree

#: How far an observation may sit from its closed form before the test fails, in standard errors.
#: Generous on purpose: with fixed seeds these numbers are deterministic, so this is not guarding
#: against bad luck, it is the width at which a real change in the model becomes undeniable.
Z_MAX = 4.0


def _z(sample: np.ndarray, expected: float) -> float:
    """How many standard errors the sample mean sits from its expectation."""
    se = sample.std(ddof=1) / math.sqrt(len(sample))
    return (sample.mean() - expected) / se


def _n_extant(result) -> int:
    return len([n for n in result.complete_tree.nodes.values()
                if n.children is None and n.fate == "extant"])


# --- the species level: a birth-death process ---------------------------------------------------

def test_yule_tree_size_matches_the_geometric_distribution():
    """A pure-birth run from one lineage for time T leaves ``N(T)`` lineages, and ``N(T)`` is
    geometric: ``E[N] = e^(λT)`` and ``P(N = 1) = e^(-λT)``. Checking the mean alone would pass for
    any process with the right average, so the single-lineage probability is checked too — it pins
    the shape, not just the centre."""
    lam, T, reps = 1.0, 2.0, 2000
    counts = np.array([_n_extant(simulate_species_tree(birth=lam, death=0.0, total_time=T, seed=s))
                       for s in range(reps)])

    assert abs(_z(counts, math.exp(lam * T))) < Z_MAX, (
        f"mean tree size {counts.mean():.4f}, expected e^(λT) = {math.exp(lam * T):.4f}")

    p1 = math.exp(-lam * T)                       # P(no split ever happened)
    se = math.sqrt(p1 * (1 - p1) / reps)
    assert abs(((counts == 1).mean() - p1) / se) < Z_MAX, (
        f"P(N=1) {(counts == 1).mean():.4f}, expected {p1:.4f}")


def test_birth_death_mean_and_extinction_probability_match_theory():
    """With birth λ and death μ over time T, the *unconditioned* mean size is ``e^((λ-μ)T)`` and the
    probability the whole run dies is ``μ(e^{rT}-1)/(λe^{rT}-μ)``.

    The extinction probability is the load-bearing half. It only comes out right if a time-limited
    run is **not** quietly conditioned on survival — if extinct runs were retried instead of raised,
    the observed rate would collapse toward zero and the mean would be inflated."""
    lam, mu, T, reps = 1.0, 0.4, 2.0, 3000
    counts, died = [], 0
    for s in range(reps):
        try:
            counts.append(_n_extant(simulate_species_tree(birth=lam, death=mu, total_time=T, seed=s)))
        except RuntimeError:                       # the run went extinct: that is an outcome, not a failure
            counts.append(0)
            died += 1
    counts = np.array(counts)

    r = lam - mu
    assert abs(_z(counts, math.exp(r * T))) < Z_MAX, (
        f"mean size {counts.mean():.4f}, expected e^(rT) = {math.exp(r * T):.4f}")

    p_ext = mu * (math.exp(r * T) - 1) / (lam * math.exp(r * T) - mu)
    se = math.sqrt(p_ext * (1 - p_ext) / reps)
    assert abs((died / reps - p_ext) / se) < Z_MAX, (
        f"P(extinct) {died / reps:.4f}, expected {p_ext:.4f}")


def test_n_extant_conditioning_is_length_weighted_not_stop_at_the_nth_birth():
    """The one most simulators get wrong, so it is ruled out rather than assumed.

    Stopping a Yule run the instant the Nth lineage is born (the naive scheme) leaves the last split
    *at* the present, so every tree ends with two zero-length terminal branches and any estimate
    reading terminal branch lengths is biased. Sampling the stopping time correctly — weighting by
    how long the tree spends at each size — instead makes the time from the last speciation to the
    present ``Exp(Nλ)``, with mean ``1/(Nλ)``.

    Three schemes are in use in the wild and they predict three different answers, so this one
    quantity separates them: 0 for the naive scheme, ``1/(2Nλ)`` for unweighted sampling, and
    ``1/(Nλ)`` for the length-weighted one."""
    lam, n, reps = 1.0, 20, 1500
    gaps = []
    for s in range(reps):
        nodes = simulate_species_tree(birth=lam, death=0.0, n_extant=n, seed=s).complete_tree.nodes
        present = max(x.end_time for x in nodes.values())
        gaps.append(present - max(x.end_time for x in nodes.values() if x.children))
    gaps = np.array(gaps)

    # the naive scheme would put every last split exactly at the present
    assert (gaps > 0).all(), "some tree ended at its last speciation — terminal branches are degenerate"

    assert abs(_z(gaps, 1 / (n * lam))) < Z_MAX, (
        f"mean time since the last split {gaps.mean():.5f}, expected 1/(Nλ) = {1 / (n * lam):.5f}")

    # and the unweighted alternative is not merely less good, it is far away
    assert _z(gaps, 1 / (2 * n * lam)) > 10, "unweighted sampling is not ruled out"

    # Exp has sd == mean; a different stopping rule would not
    assert abs(gaps.std(ddof=1) / gaps.mean() - 1.0) < 0.15, "the gap is not exponential"


# --- the genome level: what each rate is counted per ---------------------------------------------

def _fixed_tree():
    """One tree, reused so the expectations below are exact rather than averaged over trees."""
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=1).complete_tree


def test_origination_is_counted_per_lineage():
    """Origination is per lineage, so over a fixed tree of total branch length L the count is
    ``Poisson(o·L)``. Both moments are checked: a per-*copy* origination would scale with genome
    size instead, and would not stay Poisson."""
    tree = _fixed_tree()
    length = sum(n.end_time - n.birth_time for n in tree.nodes.values())
    for rate, reps in ((0.5, 400), (2.0, 400)):
        counts = np.array([
            sum(1 for e in simulate_genomes_family(tree, origination=rate, initial_families=0,
                                                   seed=s).events if e.kind == "origination")
            for s in range(reps)])
        expected = rate * length
        assert abs(_z(counts, expected)) < Z_MAX, (
            f"origination={rate}: {counts.mean():.3f} events, expected o·L = {expected:.3f}")
        # Poisson: variance equals the mean
        assert 0.7 < counts.var(ddof=1) / expected < 1.4, (
            f"origination={rate}: variance {counts.var(ddof=1):.2f} against mean {expected:.2f}")


def test_inversion_breakpoints_are_uniform_around_the_chromosome():
    """Where an inversion lands, and how much it takes.

    A rearrangement model is only neutral if its breakpoints fall **uniformly**: a sampler that
    favoured the middle of the array, or its ends, would still produce plausible-looking genomes
    while quietly making some gene orders unreachable, and no invariant on the resulting genome would
    notice. Uniform breakpoints on a ring are the statement that the arcs they cut are
    ``Dirichlet(1, …, 1)``; the chromosome here is discrete, so the same statement is that the counts
    per starting position are multinomial with equal probabilities — a chi-square.

    The **extent** is checked for its shape rather than its mean, which
    ``test_realised_extent_on_a_circle_matches_the_nominal_one`` already pins: a geometric of mean M
    puts ``1/M`` of its weight on 1 and ``(1-1/M)/M`` on 2, and a distribution with the right average
    but the wrong spread would pass on the mean alone and fail here.

    Circular, and with no duplication or loss, so the gene count stays fixed and every position is
    equally available for the whole run."""
    genes, mean_extent = 60, 6.0
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=10, seed=1).complete_tree
    starts, lengths = [], []
    for s in range(6):
        run = simulate_genomes_ordered(tree, initial_families=genes, chromosomes=1,
                                       topology="circular", inversion=3.0,
                                       inversion_extent=mean_extent, seed=s)
        for r in run.rearrangements:
            if isinstance(r, Inversion):
                starts.append(r.start)
                lengths.append(r.length)
    starts, lengths = np.array(starts), np.array(lengths)
    assert len(starts) > 5000, "too few inversions to say anything"

    observed = np.bincount(starts, minlength=genes)
    assert len(observed) == genes, "an inversion started outside the chromosome"
    expected = len(starts) / genes
    chi2 = float(((observed - expected) ** 2 / expected).sum())
    df = genes - 1
    assert abs((chi2 - df) / math.sqrt(2 * df)) < Z_MAX, (
        f"breakpoints are not uniform around the ring: chi2 {chi2:.1f} on {df} df")

    for k, p in ((1, 1 / mean_extent), (2, (1 - 1 / mean_extent) / mean_extent)):
        seen = (lengths == k).mean()
        se = math.sqrt(p * (1 - p) / len(lengths))
        assert abs((seen - p) / se) < Z_MAX, (
            f"P(extent = {k}) is {seen:.4f}, expected {p:.4f} for a geometric of mean {mean_extent}")


def test_duplication_transfer_and_loss_are_counted_per_copy():
    """The compensator identity, which pins the units *and* the scope of all three rates at once.

    For a counting process whose rate is ``r`` per live copy, the expected number of events is
    ``r × ∫(live copies)dt`` — the total gene-copy-time, which the event log is enough to
    reconstruct. A rate counted per *lineage* rather than per copy, or one off by a constant, shows
    up here immediately; a factor of two shows up as a ratio of two.

    Note the events are counted by their identity, not by their rows: a duplication or transfer
    writes one row per descendant, so counting rows would report exactly twice the truth."""
    tree = _fixed_tree()
    rates = {"duplication": 0.2, "transfer": 0.35, "loss": 0.15}
    reps = 300
    copy_time, counts = [], {kind: [] for kind in rates}

    for s in range(reps):
        g = simulate_genomes_family(tree, initial_families=20, seed=s, **rates)
        born, ended = {}, {}
        for e in g.events:
            born.setdefault(e.copy, (e.time, e.lineage))
            if e.parent is not None:
                ended.setdefault(e.parent, e.time)
            if e.kind == "loss":
                ended[e.copy] = e.time
        copy_time.append(sum(ended.get(c, tree.nodes[lin].end_time) - t
                             for c, (t, lin) in born.items()))
        for kind in rates:
            counts[kind].append(len({e.event for e in g.events if e.kind == kind}))

    total = float(np.mean(copy_time))
    for kind, rate in rates.items():
        observed = np.array(counts[kind])
        assert abs(_z(observed, rate * total)) < Z_MAX, (
            f"{kind}={rate}: {observed.mean():.2f} events against rate × gene-copy-time "
            f"= {rate * total:.2f} (ratio {observed.mean() / (rate * total):.4f})")
