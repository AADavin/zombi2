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
  implies rather than merely the right average;
- a mass extinction culls the fraction it was given, and ``sampling`` observes the fraction it was
  given;
- a transfer picks its recipient uniformly among the lineages alive at that instant;
- a **substitution model** puts the distance between an ancestor and its descendant where its own
  transition probability does: Jukes–Cantor's ``3/4·(1 − e^(−4d/3))`` for JC69 and Kimura's separate
  transition and transversion probabilities for K80 — written out from the 1969 and 1980 papers, not
  taken from ``p_matrix()``, which is the code under test. HKY85 is checked on the two things that
  follow from building ``Q`` out of ``π`` correctly and are invisible otherwise: the composition
  stays at ``π`` at the tips of a deep tree, and changes ``i→j`` and ``j→i`` balance, which is
  reversibility itself. A **protein** matrix gets the same two checks over twenty states, where a
  20x20 eigendecomposition gives an indexing error more room to hide — but note what they cannot
  reach: ``Q`` is assembled as ``S_ij * pi_j``, so any pi handed in is stationary by construction and
  neither check can detect a mis-transcribed published table;
- a **site profile** puts row *i* on site *i*. Over an equal-exchangeability model a profiled site is
  Felsenstein's F81, whose ``S·(1 − e^(−d/S))`` with ``S = 1 − Σπ²`` gives each block of a profile its
  own closed form — and the total over sites cannot catch a shuffled profile, because summing over
  positions is permutation-invariant. The profile's other promise, that it changes *where* residues
  belong and not *which pairs interchange*, is checked separately: changes still follow
  ``π_i·π_j·S_ij`` with the base model's exchangeabilities, which rules out a rebuild that keeps the
  frequencies and drops the chemistry;
- **across-site rate classes** reach the sites: ``+I+Γ`` gives the Jukes–Cantor curve *averaged over*
  the classes rather than evaluated at their mean, which by Jensen's inequality is strictly fewer
  differences — so a set of classes that is computed, normalised and then never applied is ruled out;
- the **lineage clocks** are mean 1, so a relaxed clock does not inflate every branch in the tree.
  This is where the historical lognormal bug lived, and the uncorrected draw is the alternative
  ruled out;
- a **Brownian trait** moves by ``Normal(0, σ²·Δt)`` on each branch, and its tips separate by ``σ²×``
  the path between them — Felsenstein's covariance, the statement that shared ancestry makes
  relatives similar; **Ornstein–Uhlenbeck** hits its exact transition moments, mean *and* variance,
  which is what distinguishes a pull from a diffusion that happens to be near an optimum;
- an **Mk trait** fires each transition at the rate that *direction* was given — the compensator
  identity again, and the check that catches a transposed rate dict — and its branches end in the
  other state as often as the chain's ``(1 − e^(−2qΔt))/2`` says, which pins the saturation a rate
  alone does not;
- **correlated traits** come out at the ρ they were declared with. The null a comparative method is
  graded against is a run with the correlation switched off, so if the overlay were not reaching the
  draw, the signal and the null would be the same run;
- **a driven rate realises the multiplier it was declared with** — on the conditioned path, on the
  joint one, and on the sequence level, where it is not statistical at all: a species-phylogram
  branch is the driven rate *integrated* over the branch, so when the driver switches part-way along,
  the length the engine wrote can be checked against the trait's own segments to floating point, and
  the two "sample the driver once per branch" wirings are ruled out by how far they would land from
  it. These matter most and are checked last, because a mis-wired driver is the one error
  that leaves every output well formed: a tree, a genome and a log that are all internally consistent,
  with only the strength of the association wrong, which is precisely what the run was made to
  measure.

Where the sequence and trait checks pool across a tree, note **which standard error applies**. Under
JC69, K80 and a symmetric Mk chain, whether a site or a lineage ends up somewhere else does not
depend on where it started, so those indicators are independent across branches even though the
sequences and the states are not, and the binomial standard error is the right one. Base composition,
reversibility and the rate classes are *not* state-independent in that way — sites within one run are
correlated through their shared ancestry — so those are measured across independent replicates, each
its own tree, with the spread taken from the replicates themselves. Using a binomial standard error
there would understate the spread and read ordinary noise as a broken model.

**Every test is deterministic.** The seeds are fixed, so a "statistical" test here cannot flake: it
computes one number and compares it to one expectation. The tolerances are wide (``|z| < 4``) because
their job is to catch a model that has changed, not to detect a rate that is off by a percent — and
the hypotheses these rule out are rejected by tens of standard deviations, not by twos.

What that buys, measured by mutating the model rather than assumed. Scaling a rate by 1.15 behind the
engine's back fails these; by 1.02 it does not. Biasing 10% of inversions toward one half of the ring
fails; 3% does not. Mis-wiring the driver strength by 5% is caught on the conditioned path and by
10% on the joint one. At the sequence and trait levels, scaling the substitution rate by 1.02 fails,
κ by 1.03, the Brownian variance-rate by 1.05, the OU pull by 1.02, an Mk transition rate by 1.03 and
a declared correlation by 1.05; shifting one stationary frequency by 0.005 fails by 15 standard
errors, moving a profile row's dominant frequency from 0.80 to 0.81 by 3.4, and inflating the lineage
clock — the historical bug's signature — by 1.10. So resolution sits between a fraction of a percent
and about ten percent, depending on the check.

The **exchangeabilities under a profile are the loose one**, and knowing which check is weakest is
part of reading this file: one exchangeability has to be off by about 20% before it shows, because
pinning six of them needs branches short enough that a difference is still a substitution, which caps
how many changes there are to count. What it does catch decisively is the failure that path invites —
rebuilding each site from its frequencies and losing the chemistry altogether, which it rejects by
270 standard errors.

That figure is a property of *these* trees and replicate counts, not of the method: more material, or
more replicates, tightens it. Which is the honest way to read this file — it is evidence that the
model is the one advertised, not a calibration of it. A rate quietly off by a percent still gets
past, and catching that costs more replicates than a test suite should spend on every commit.

If one of these fails, the arithmetic below is not what to doubt first. Something about what the
engine samples has changed.
"""
from __future__ import annotations

import math
import re
from collections import Counter

import numpy as np

from zombi2.genomes import simulate_genomes_family, simulate_genomes_ordered
from zombi2.genomes.ordered import Inversion
from zombi2.joint import simulate_joint
from zombi2.rates.modifiers import ByLineage, DrivenBy, FromParent
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import BASES, gtr, hky85, jc69, k80, lg, poisson
from zombi2.species import simulate_species_tree
from zombi2.traits import discrete, simulate_continuous, simulate_discrete

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


def test_a_mass_extinction_culls_the_fraction_it_was_given():
    """At the pulse each standing lineage is lost independently with probability ``f``.

    The tree has to be **large** at the pulse, and that is not a convenience. A pulse can wipe out a
    small tree entirely, and a run that goes extinct raises rather than returning — so measuring on
    the runs that came back conditions on survival, and runs where the pulse bit hardest are exactly
    the ones missing. Measured on a tree of a few standing lineages that bias shows up as roughly
    four standard errors of apparent *under*-culling, which reads convincingly like a real defect.
    With hundreds of lineages standing, total wipeout has no measurable probability and the estimate
    is honest."""
    pulse, fraction = 4.0, 0.75
    standing = culled = 0
    for s in range(15):
        # tall enough that no seed here is wiped out, so nothing is excluded and nothing is biased
        run = simulate_species_tree(birth=2.0, death=0.0, total_time=pulse + 0.1,
                                    mass_extinctions=[(pulse, fraction)], seed=s)
        alive = [n for n in run.complete_tree.nodes.values()
                 if n.birth_time < pulse <= n.end_time]
        standing += len(alive)
        culled += sum(1 for n in alive
                      if not n.children and abs(n.end_time - pulse) < 1e-9)

    assert standing > 20_000, f"only {standing} lineages stood at the pulse — too few to be unbiased"
    observed = culled / standing
    se = math.sqrt(fraction * (1 - fraction) / standing)
    assert abs((observed - fraction) / se) < Z_MAX, (
        f"the pulse culled {observed:.4f} of standing lineages, not {fraction}")


def test_incomplete_sampling_observes_the_fraction_it_was_given():
    """``sampling=ρ`` observes each survivor independently with probability ρ; the rest are
    ``unsampled`` — present in the complete tree, absent from the extant one. A ρ applied to the
    process rather than to the observation would change the tree's shape instead of its labelling,
    so what is checked is the labelling of survivors."""
    rho = 0.4
    observed = unsampled = 0
    for s in range(120):
        run = simulate_species_tree(birth=1.0, death=0.2, n_extant=30, sampling=rho, seed=s)
        for n in run.complete_tree.nodes.values():
            if n.children is None:
                observed += n.fate == "extant"
                unsampled += n.fate == "unsampled"

    total = observed + unsampled
    seen = observed / total
    se = math.sqrt(rho * (1 - rho) / total)
    assert abs((seen - rho) / se) < Z_MAX, (
        f"{seen:.4f} of survivors were observed, not ρ = {rho}")


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
        for e in g.edges:
            born.setdefault(e.copy, (e.time, e.lineage))
            if e.parent is not None:
                ended.setdefault(e.parent, e.time)
            if e.kind == "loss":
                ended[e.copy] = e.time
        copy_time.append(sum(ended.get(c, tree.nodes[lin].end_time) - t
                             for c, (t, lin) in born.items()))
        for kind in rates:
            counts[kind].append(len({e.event for e in g.edges if e.kind == kind}))

    total = float(np.mean(copy_time))
    for kind, rate in rates.items():
        observed = np.array(counts[kind])
        assert abs(_z(observed, rate * total)) < Z_MAX, (
            f"{kind}={rate}: {observed.mean():.2f} events against rate × gene-copy-time "
            f"= {rate * total:.2f} (ratio {observed.mean() / (rate * total):.4f})")


def test_transfer_picks_its_recipient_uniformly():
    """Who receives a transfer, under ``transfer_to="uniform"``.

    The same concern as the inversion breakpoints, one level up: a recipient sampler that quietly
    favoured some lineages — the deepest, the most recently born, the lowest-numbered — would still
    produce a plausible genome for everyone and would still satisfy every invariant on the event log.
    Only the distribution shows it.

    Which lineages are even eligible changes through the run, so a flat count per lineage is not the
    expectation. Each transfer offers ``1/k`` to each of the ``k`` lineages alive at that instant, so
    a lineage's expected receipts are the sum of ``1/k`` over the transfers it was exposed to."""
    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=1).complete_tree
    spans = {n.id: (n.birth_time, n.end_time) for n in tree.nodes.values()}
    received, exposure, total = Counter(), Counter(), 0

    for s in range(12):
        for e in simulate_genomes_family(tree, initial_families=30, transfer=0.6, seed=s).edges:
            if e.kind != "transfer" or e.recipient is None:
                continue
            eligible = [i for i, (b, d) in spans.items() if b < e.time <= d and i != e.lineage]
            if len(eligible) < 2:
                continue
            received[e.recipient] += 1
            total += 1
            for i in eligible:
                exposure[i] += 1 / len(eligible)

    assert total > 10_000, f"only {total} transfers to judge by"
    keys = [k for k in exposure if exposure[k] >= 20]        # enough exposure for chi-square to apply
    chi2 = sum((received[k] - exposure[k]) ** 2 / exposure[k] for k in keys)
    df = len(keys) - 1
    assert abs((chi2 - df) / math.sqrt(2 * df)) < Z_MAX, (
        f"recipients are not uniform: chi2 {chi2:.1f} on {df} df")


# --- driven rates: a rate that reads another level ------------------------------------------------

def _state_at(history, node, tree, when):
    """Which state a lineage held at ``when``, walking its own segment history."""
    at = tree.nodes[node].birth_time
    segments = history.get(node) or []
    for state, duration in segments:
        if when <= at + duration + 1e-12:
            return state
        at += duration
    return segments[-1][0] if segments else None


def test_a_conditioned_rate_realises_the_multiplier_it_was_given():
    """The driving mechanism itself, on the conditioned path (the driver is a finished level).

    This is the one a downstream invariant cannot catch. If ``DrivenBy`` attached its factors to the
    wrong branches, or applied them a step late, every output would still be well formed — a tree, a
    genome, an event log, all internally consistent — and only the *strength* of the association
    would be wrong, which is exactly the quantity a user is measuring when they reach for a
    conditioned run.

    Origination is per lineage, so the compensator is simply lineage-time: a lineage sitting in state
    ``s`` originates at ``o × m(s)`` per unit time, and each state's realised rate is checked against
    its own declared value rather than against the other's — a ratio would pass even if both were
    scaled by the same wrong constant."""
    factors, base = {"hot": 5.0, "cold": 1.0}, 1.4
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=25, seed=1).complete_tree
    seen, lineage_time = Counter(), Counter()

    for s in range(120):
        habitat = simulate_discrete(tree, states=list(factors), switch=0.5, start="hot", seed=s)
        for segments in habitat.history.values():
            for state, duration in segments:
                lineage_time[state] += duration
        run = simulate_genomes_family(tree, initial_families=0, seed=s,
                                      origination=base * DrivenBy(habitat, factors))
        for e in run.edges:
            if e.kind == "origination":
                seen[_state_at(habitat.history, e.lineage, tree, e.time)] += 1

    for state, factor in factors.items():
        realised = seen[state] / lineage_time[state]
        expected = base * factor
        se = math.sqrt(seen[state]) / lineage_time[state]      # Poisson count, known exposure
        assert abs((realised - expected) / se) < Z_MAX, (
            f"{state}: origination ran at {realised:.4f} per unit time, but was told to run at "
            f"base × {factor} = {expected:.4f}")


def test_the_parallel_engine_realises_a_driven_multiplier_too():
    """The same check on the **per-family** engine, which computes the driven rate a different way.

    The serial loop sums a driven rate over every living lineage at once; the parallel engine evolves
    one family at a time and sums it over that family's footprint — the lineages that family happens
    to occupy. Those are two different arithmetics for the same process, and only one of them was
    ever checked. If the footprint version attached the driver to the wrong branch, or missed a
    switch because the horizon only looks at occupied lineages, the run would still be a well-formed
    genome and only the driver *strength* would be wrong.

    Origination again, for the same reason: the compensator is plain lineage-time."""
    factors, base = {"hot": 5.0, "cold": 1.0}, 1.4
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=25, seed=1).complete_tree
    seen, lineage_time = Counter(), Counter()

    for s in range(60):
        habitat = simulate_discrete(tree, states=list(factors), switch=0.5, start="hot", seed=s)
        for segments in habitat.history.values():
            for state, duration in segments:
                lineage_time[state] += duration
        run = simulate_genomes_family(tree, initial_families=0, seed=s, parallel=1,
                                      origination=base * DrivenBy(habitat, factors))
        for e in run.edges:
            if e.kind == "origination":
                seen[_state_at(habitat.history, e.lineage, tree, e.time)] += 1

    for state, factor in factors.items():
        realised = seen[state] / lineage_time[state]
        expected = base * factor
        se = math.sqrt(seen[state]) / lineage_time[state]
        assert abs((realised - expected) / se) < Z_MAX, (
            f"{state}: the parallel engine originated at {realised:.4f} per unit time, but was told "
            f"to run at base × {factor} = {expected:.4f}")


def test_a_joint_rate_realises_the_multiplier_it_was_given():
    """The same check on the joint path, where the driver is growing at the same time.

    Harder than the conditioned case, and the more likely to be subtly wrong: the trait is evolving
    along the very tree its own speciation rate is building, so the rate has to be re-read as the
    state changes, mid-race. Getting that wrong — reading a stale state, or missing the re-read at a
    split — gives a tree that looks entirely ordinary and a state-dependence of the wrong strength.

    Speciation is per lineage, so again the exposure is lineage-time in each state, and each state's
    realised rate is checked against its own ``λ × m(s)``."""
    factors, base = {"fast": 3.0, "slow": 1.0}, 1.0
    splits, lineage_time = Counter(), Counter()

    for s in range(60):
        run = simulate_joint(birth=base * DrivenBy("trait", factors),
                             trait=discrete(states=list(factors), switch=0.4),
                             n_extant=120, seed=s)
        for segments in run.trait.history.values():
            for state, duration in segments:
                lineage_time[state] += duration
        for node in run.complete_tree.nodes.values():
            if node.children:                                   # it split: credit the state it held
                segments = run.trait.history.get(node.id)
                if segments:
                    splits[segments[-1][0]] += 1

    for state, factor in factors.items():
        realised = splits[state] / lineage_time[state]
        expected = base * factor
        se = math.sqrt(splits[state]) / lineage_time[state]
        assert abs((realised - expected) / se) < Z_MAX, (
            f"{state}: lineages split at {realised:.4f} per unit time, but the trait was told to "
            f"make them split at λ × {factor} = {expected:.4f}")


# --- the sequence level: a substitution model down a gene tree ------------------------------------
#
# A branch length here is in **substitutions/site**, and every model on the menu is normalised to one
# expected substitution per site per unit branch length. So a branch of length d has a transition
# probability the model itself fixes, and each test below writes that probability out by hand — from
# Jukes and Cantor (1969) and Kimura (1980), not from `p_matrix()`, which is the code under test.

#: no duplication, transfer, loss or origination, so every family's gene tree *is* the species tree
#: and each lineage carries exactly one copy — which makes ``n<species>_g<copy>`` an unambiguous
#: handle on "the sequence this lineage ended with".
_NO_GENOME_EVENTS = {"duplication": 0.0, "transfer": 0.0, "loss": 0.0, "origination": 0.0}

def _coder(alphabet: str) -> np.ndarray:
    """A lookup from byte to state index, for counting with ``bincount`` rather than a Python loop
    over sites. Out-of-alphabet bytes map to 255, which indexes out of bounds rather than quietly
    counting as a state."""
    table = np.full(128, 255, dtype=np.uint8)
    for index, letter in enumerate(alphabet):
        table[ord(letter)] = index
    return table


_ACGT = _coder(BASES)


def _codes(sequence: str, table: np.ndarray = _ACGT) -> np.ndarray:
    return table[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]


def _one_copy_per_lineage(*, n_extant: int, families: int, seed: int):
    """A species tree and a genome run with no D/T/L/O — one gene copy per lineage per family."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=seed).complete_tree
    return tree, simulate_genomes_family(tree, initial_families=families, seed=seed,
                                         **_NO_GENOME_EVENTS)


def _sequence_branches(tree, run, family: int, substitution: float, table: np.ndarray = _ACGT):
    """``(ancestor, descendant, branch length)`` for every branch of one family.

    The root branch counts too, and it is the one that quietly goes missing: its ancestor is the
    family's ``founding`` sequence — drawn from the model's stationary frequencies at the origination
    point, ``t = 0`` for a family the run started with — and ``founding`` is deliberately *not* in
    ``ancestral``, whose keys pair one-to-one with phylogram nodes. So a walk that only pairs a node
    with its parent silently drops one branch per family."""
    seqs = {int(label.split("_", 1)[0][1:]): _codes(s, table)
            for label, s in {**run.ancestral[family], **run.alignments[family]}.items()}
    founding = _codes(run.founding[family], table)
    for node in tree.nodes.values():
        ancestor = founding if node.parent is None else seqs[node.parent]
        yield ancestor, seqs[node.id], substitution * (node.end_time - node.birth_time)


def test_jc69_site_differences_match_the_jukes_cantor_distance():
    """Under JC69 a site differs across a branch of length ``d`` with probability
    ``3/4·(1 − e^(−4d/3))`` — the Jukes–Cantor formula, written out here rather than taken from
    `p_matrix()`.

    One number checks the whole chain at once: that ``Q`` is normalised to one substitution per site
    per unit branch length, that the exponential is applied over the right ``d``, and that a
    gene-tree branch is ``substitution × Δt``. A ``Q`` scaled by any constant fails, because the
    saturating curve pins the *scale*, not just the ordering.

    The binomial standard error is the right one despite the shared ancestry: under JC69 whether a
    site changes over a branch does not depend on which base it started from, so the indicators are
    independent across sites **and** across branches, even though the sequences are not.

    The alternative ruled out is the one every textbook warns about — reading ``d`` itself as the
    expected fraction of differing sites, i.e. ignoring multiple hits at the same site."""
    substitution = 0.6
    tree, genomes = _one_copy_per_lineage(n_extant=12, families=8, seed=1)
    run = simulate_sequences(genomes, model=jc69(), length=400, substitution=substitution, seed=1)

    observed = expected = naive = variance = 0.0
    for family in genomes.gene_trees:
        for ancestor, descendant, d in _sequence_branches(tree, run, family, substitution):
            sites = len(ancestor)
            observed += float((ancestor != descendant).sum())
            p = 0.75 * (1 - math.exp(-4 * d / 3))
            expected += p * sites
            variance += p * (1 - p) * sites
            naive += min(d, 1.0) * sites          # the uncorrected p-distance, capped to stay a probability

    assert observed > 10_000, "too few substitutions to say anything"
    assert abs((observed - expected) / math.sqrt(variance)) < Z_MAX, (
        f"{observed:.0f} sites differ, but Jukes–Cantor expects {expected:.0f} "
        f"(ratio {observed / expected:.4f})")
    assert abs((observed - naive) / math.sqrt(variance)) > 10, (
        "the uncorrected p-distance is not ruled out — multiple hits may not be accumulating")


def test_k80_transitions_and_transversions_match_kimuras_two_probabilities():
    """Under K80 the two kinds of change have separate closed forms (Kimura 1980), and κ is what
    sets them apart. With ``Q`` normalised, the transversion rate is ``β = 1/(κ+2)`` and the
    transition rate ``α = κ/(κ+2)``, giving over a branch of length ``d``

        ``P(transition)   = 1/4 + 1/4·e^(−4βd) − 1/2·e^(−2(α+β)d)``
        ``P(transversion) = 1/2 − 1/2·e^(−4βd)``     (both transversions together)

    Checking the two separately, each against its own expression, is what pins κ. The ts/tv *ratio*
    alone would not: a model that got both probabilities wrong by the same factor would keep the
    ratio and fail here. And at κ = 1 the pair collapses to Jukes–Cantor, so a κ that was parsed but
    never reached the matrix is the alternative ruled out at the end."""
    kappa, substitution = 4.0, 0.6
    alpha, beta = kappa / (kappa + 2), 1 / (kappa + 2)
    # more material than the JC69 check above: κ is the hardest of these numbers to pin, because it
    # only shows in *which* of two changes happened rather than in whether one did at all
    tree, genomes = _one_copy_per_lineage(n_extant=18, families=24, seed=1)
    run = simulate_sequences(genomes, model=k80(kappa), length=1200, substitution=substitution,
                             seed=1)

    #: purine (A, G) = 0, pyrimidine (C, T) = 1 — a change within a class is a transition
    ring = np.array([0, 1, 0, 1], dtype=np.uint8)
    seen = {"transition": 0.0, "transversion": 0.0}
    expected = {"transition": 0.0, "transversion": 0.0}
    variance = {"transition": 0.0, "transversion": 0.0}
    flat = 0.0                                     # what κ = 1 (Jukes–Cantor) would predict for transitions

    for family in genomes.gene_trees:
        for ancestor, descendant, d in _sequence_branches(tree, run, family, substitution):
            changed = ancestor != descendant
            same_ring = ring[ancestor] == ring[descendant]
            sites = len(ancestor)
            counts = {"transition": float((changed & same_ring).sum()),
                      "transversion": float((changed & ~same_ring).sum())}
            p = {"transition": 0.25 + 0.25 * math.exp(-4 * beta * d)
                               - 0.5 * math.exp(-2 * (alpha + beta) * d),
                 "transversion": 0.5 - 0.5 * math.exp(-4 * beta * d)}
            for kind in seen:
                seen[kind] += counts[kind]
                expected[kind] += p[kind] * sites
                variance[kind] += p[kind] * (1 - p[kind]) * sites
            flat += (0.25 + 0.25 * math.exp(-4 * d / 3) - 0.5 * math.exp(-4 * d / 3)) * sites

    for kind in seen:
        assert abs((seen[kind] - expected[kind]) / math.sqrt(variance[kind])) < Z_MAX, (
            f"κ={kappa}: {seen[kind]:.0f} {kind}s against Kimura's {expected[kind]:.0f} "
            f"(ratio {seen[kind] / expected[kind]:.4f})")

    assert abs((seen["transition"] - flat) / math.sqrt(variance["transition"])) > 10, (
        f"κ={kappa} is not distinguishable from κ=1 — the ratio may never reach the matrix")


def test_hky85_holds_the_composition_it_was_given_and_stays_reversible():
    """Two consequences of building ``Q`` from ``π`` correctly, neither visible in a well-formed run.

    **The composition stays at π.** Every sequence starts drawn from π, so if π really is stationary
    for ``Q`` the composition is still π at the tips of a deep tree (which is where it is measured
    here — the extant alignments). Build ``Q`` with π on the wrong axis —
    a transpose, the classic — and π stops being its stationary distribution, so the sequences drift
    away from the frequencies the run was asked for while every output stays perfectly well formed.

    **The change matrix is symmetric.** Reversibility says ``π_i·P_ij = π_j·P_ji``, so ancestor→
    descendant changes ``i→j`` and ``j→i`` are equally common *whatever* the branch lengths are.
    That is a statement about ``Q`` alone, so it holds pooled over a whole tree with no expectation
    to compute.

    Unlike the JC69 and K80 checks, neither quantity is state-independent, so sites within a run are
    correlated through their shared ancestry and a binomial standard error would be too small. Both
    are therefore measured across **independent replicates**, each its own tree and its own
    sequences, with the spread taken from the replicates themselves."""
    frequencies, substitution, reps = (0.4, 0.3, 0.2, 0.1), 1.0, 60
    composition, asymmetry = [], {pair: [] for pair in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))}

    for seed in range(reps):
        tree, genomes = _one_copy_per_lineage(n_extant=8, families=4, seed=seed)
        run = simulate_sequences(genomes, model=hky85(2.5, frequencies), length=400,
                                 substitution=substitution, seed=seed)
        bases = np.zeros(4)
        changes = np.zeros((4, 4))
        for family in genomes.gene_trees:
            for sequence in run.alignments[family].values():
                bases += np.bincount(_codes(sequence), minlength=4)
            for ancestor, descendant, _ in _sequence_branches(tree, run, family, substitution):
                changed = ancestor != descendant
                np.add.at(changes, (ancestor[changed], descendant[changed]), 1)
        composition.append(bases / bases.sum())
        for i, j in asymmetry:
            both = changes[i, j] + changes[j, i]
            asymmetry[(i, j)].append((changes[i, j] - changes[j, i]) / both)

    composition = np.array(composition)
    for i, pi in enumerate(frequencies):
        assert abs(_z(composition[:, i], pi)) < Z_MAX, (
            f"{'ACGT'[i]} settles at {composition[:, i].mean():.4f} of the sites, but the model was "
            f"given π = {pi}")

    for (i, j), sample in asymmetry.items():
        assert abs(_z(np.array(sample), 0.0)) < Z_MAX, (
            f"{'ACGT'[i]}→{'ACGT'[j]} and {'ACGT'[j]}→{'ACGT'[i]} do not balance "
            f"(relative asymmetry {np.mean(sample):+.4f}) — the model is not reversible")


def test_an_empirical_protein_matrix_is_simulated_at_its_own_stationary_distribution():
    """A twenty-state model, checked end to end: the run stays at the model's own π and its changes
    balance.

    **Be precise about what this can and cannot catch**, because the obvious reading is wrong.
    `reversible()` assembles ``Q_ij = S_ij · π_j`` from a symmetric ``S``, so whatever π it is handed
    is stationary for the ``Q`` it builds, and detailed balance holds by construction. Neither half
    of this test can therefore detect a **transcription** error in the published table: swap two
    entries of ``_LG_PI`` and the model simply believes the swapped numbers, which are stationary for
    the matrix built from them, and everything here still passes (verified). Checking 190 published
    numbers needs an independent copy of them, which this suite does not have.

    What it does check is the **engine**, over twenty states rather than four: that a founding
    sequence is drawn from the model's frequencies rather than uniformly; that ``exp(Qt)`` — computed
    by eigendecomposing ``diag(√π)·Q·diag(1/√π)`` — preserves that distribution instead of drifting
    away from it down a deep tree; and that the simulation applies the transition matrix in the right
    orientation, since applying ``P`` transposed would break the ``i→j`` / ``j→i`` balance while
    leaving every sequence well formed. On the four-state models the same properties ride on
    `test_hky85_holds_the_composition_it_was_given_and_stays_reversible`; here they are exercised on
    a 20×20 eigendecomposition, where an indexing error has far more room to hide.

    With 190 pairs the asymmetry is pooled into one χ² per degree of freedom rather than tested pair
    by pair, and pairs too rare to have a χ² distribution are dropped rather than trusted. The
    composition is a maximum over 20 residues, so it is a wider net than the ``|z| < 4`` on any one
    of them suggests: under a correct matrix the largest of 20 deviations sits near 2.5 standard
    errors by construction."""
    model, reps = lg(), 40
    table = _coder(model.alphabet)
    composition, asymmetry = [], []

    for seed in range(reps):
        tree, genomes = _one_copy_per_lineage(n_extant=8, families=4, seed=seed)
        run = simulate_sequences(genomes, model=model, length=300, substitution=1.0, seed=seed)
        residues = np.zeros(model.k)
        changes = np.zeros((model.k, model.k))
        for family in genomes.gene_trees:
            for sequence in run.alignments[family].values():
                residues += np.bincount(_codes(sequence, table), minlength=model.k)
            for ancestor, descendant, _ in _sequence_branches(tree, run, family, 1.0, table):
                changed = ancestor != descendant
                np.add.at(changes, (ancestor[changed], descendant[changed]), 1)
        composition.append(residues / residues.sum())

        upper = np.triu_indices(model.k, 1)
        forward, backward = changes[upper], changes.T[upper]
        common = (forward + backward) >= 10           # rarer pairs have no χ² distribution to speak of
        asymmetry.append(float((((forward - backward) ** 2)[common]
                                / (forward + backward)[common]).sum()) / int(common.sum()))

    composition = np.array(composition)
    for index, pi in enumerate(model.stationary):
        assert abs(_z(composition[:, index], pi)) < Z_MAX, (
            f"{model.name}: {model.alphabet[index]} settles at {composition[:, index].mean():.4f} of "
            f"the sites, but the published frequency is {pi:.4f}")

    assert abs(_z(np.array(asymmetry), 1.0)) < Z_MAX, (
        f"{model.name}: forward and reverse changes give χ²/df = {np.mean(asymmetry):.4f} rather "
        f"than 1 — the published matrix is not being used reversibly")


def test_a_site_profile_governs_the_site_it_belongs_to():
    """``profiles=``: row *i* of the array is the equilibrium frequencies at position *i*, and the
    thing worth proving is that row *i* really lands on site *i*.

    A profile that was transposed, shifted by one, or collapsed so that one row governed the whole
    gene would still produce a well-formed alignment of exactly the right length, over the right
    alphabet, with a sensible-looking overall composition. Only a per-site expectation catches it, and
    the *total* over sites cannot: summing over positions is permutation-invariant, so a shuffled
    profile has the same total divergence as the right one. So the profile here is **blocked** — a
    flat half and a sharply peaked half — and each half is checked against its own closed form.

    That closed form exists because the base model is `poisson()`, whose exchangeabilities are all
    equal. A site of an equal-exchangeability model over frequencies ``π`` is Felsenstein's F81, and
    F81 has a closed form for the chance a site differs across a branch of length ``d``:

        ``P(differ) = S·(1 − e^(−d/S))``   with ``S = 1 − Σπ²``

    — which at uniform ``π`` is ``3/4·(1 − e^(−4d/3))``, Jukes–Cantor, as it must be. ``S`` is what a
    profile changes: a peaked site has a low ``S`` and saturates early, because most of the time the
    residue it mutates to is the one it already was. Here the two blocks sit at ``S = 0.95`` and
    ``S = 0.36``, so they saturate in visibly different places and each block's expectation is
    strongly wrong for the other — which is the assertion at the end.

    The formula also pins the **normalisation**. Every per-site model is rebuilt from the base
    model's exchangeabilities and its own row, and rebuilt models are renormalised to one expected
    substitution per site per unit branch length; if that step were dropped, a peaked profile would
    quietly run at a different rate and the phylogram would stop meaning substitutions per site."""
    model, block, substitution, reps = poisson(), 120, 1.0, 12
    table = _coder(model.alphabet)
    k = model.k

    flat = np.full(k, 1 / k)
    peaked = np.full(k, 0.2 / (k - 1))
    peaked[0] = 0.8
    rows = np.vstack([np.tile(flat, (block, 1)), np.tile(peaked, (block, 1))])
    halves = {"flat": (np.arange(block), flat),
              "peaked": (np.arange(block, 2 * block), peaked)}
    spread = {name: 1 - float((pi ** 2).sum()) for name, (_, pi) in halves.items()}

    seen = dict.fromkeys(halves, 0.0)
    expected = dict.fromkeys(halves, 0.0)
    variance = dict.fromkeys(halves, 0.0)
    swapped = dict.fromkeys(halves, 0.0)        # what the *other* half's row would predict
    composition = {name: [] for name in halves}

    for seed in range(reps):
        tree, genomes = _one_copy_per_lineage(n_extant=10, families=2, seed=seed)
        run = simulate_sequences(genomes, model=model, length=2 * block,
                                 substitution=substitution,
                                 profiles={f: rows for f in genomes.gene_trees}, seed=seed)
        residues = {name: np.zeros(k) for name in halves}
        for family in genomes.gene_trees:
            for sequence in run.alignments[family].values():
                coded = _codes(sequence, table)
                for name, (sites, _) in halves.items():
                    residues[name] += np.bincount(coded[sites], minlength=k)
            for ancestor, descendant, d in _sequence_branches(tree, run, family, substitution, table):
                for name, (sites, _) in halves.items():
                    other = spread["peaked" if name == "flat" else "flat"]
                    p = spread[name] * (1 - math.exp(-d / spread[name]))
                    seen[name] += float((ancestor[sites] != descendant[sites]).sum())
                    expected[name] += p * len(sites)
                    variance[name] += p * (1 - p) * len(sites)
                    swapped[name] += other * (1 - math.exp(-d / other)) * len(sites)
        for name in halves:
            composition[name].append(residues[name] / residues[name].sum())

    for name, (_, pi) in halves.items():
        assert abs((seen[name] - expected[name]) / math.sqrt(variance[name])) < Z_MAX, (
            f"the {name} half of the profile differs at {seen[name]:.0f} sites against the "
            f"F81 curve its own row implies, {expected[name]:.0f} (S = {spread[name]:.3f})")
        assert abs((seen[name] - swapped[name]) / math.sqrt(variance[name])) > 10, (
            f"the {name} half is not distinguishable from the other half's row — the rows may not be "
            f"reaching the sites they belong to")

        # Where the mass sits, on the residue the peaked row concentrates it on: 0.8 in that half and
        # 1/20 in the flat one. One well-conditioned number per half rather than all 20 frequencies —
        # the rare residues of a peaked row carry too few counts for their spread to be estimated from
        # a dozen replicates, and `test_an_empirical_protein_matrix…` already pins a full frequency
        # vector where the counts support it. Between this and the F81 curve above the whole row is
        # constrained anyway: S fixes how concentrated the row is, this fixes where.
        dominant = np.array([replicate[0] for replicate in composition[name]])
        assert abs(_z(dominant, pi[0])) < Z_MAX, (
            f"in the {name} half {model.alphabet[0]} settles at {dominant.mean():.4f} of the sites "
            f"but its row says {pi[0]:.4f}")
        other = halves["peaked" if name == "flat" else "flat"][1]
        assert abs(_z(dominant, other[0])) > 20, (
            f"the {name} half's composition is not distinguishable from the other half's row")


def test_a_profile_keeps_the_chemistry_of_the_model_it_decorates():
    """The other half of what a profile promises: it says **where** residues belong, not **which
    pairs interchange**.

    A profile replaces a model's frequencies site by site while keeping its exchangeabilities, so a
    run under a profile should still show the base model's chemistry. The failure this rules out is
    the plausible one — rebuilding each site from its frequencies alone and losing the exchangeability
    matrix on the way, which turns every model into F81. Nothing downstream would look wrong: the
    composition would still match the profile (the test above would still pass), the sequences would
    still be well formed, and only *which* substitutions happen would be homogenised.

    Under a reversible model at stationarity, changes between ``i`` and ``j`` happen at a rate
    proportional to ``π_i·π_j·S_ij``, so with one profile shared by every site the pooled change
    counts should follow that product — with ``S`` the exchangeabilities **this test supplied** to
    ``gtr()``, not any matrix read back out of the model.

    ``π_i·π_j·S_ij`` is a statement about *rates*, and observed differences are only substitutions
    while branches are short enough that a site rarely changes twice. So the substitution rate here
    is deliberately small: at ``0.04`` the second-order correction sits below the counting noise,
    while at the ``0.6`` the tests above use it does not, and this same check would report a large
    and entirely spurious misfit. Four states rather than twenty for the same reason — six
    exchangeabilities need far less material to pin than 190, and the run has to stay short."""
    exchangeabilities = (1.0, 3.0, 0.5, 2.0, 4.0, 1.5)      # AC AG AT CG CT GT, deliberately uneven
    profile = np.array([0.4, 0.35, 0.15, 0.10])
    substitution, length, families, reps = 0.04, 300, 3, 6

    upper = np.triu_indices(4, 1)
    S = np.zeros((4, 4))
    S[upper] = exchangeabilities
    S = S + S.T

    rows = np.tile(profile, (length, 1))
    changes = np.zeros((4, 4))
    for seed in range(reps):
        tree, genomes = _one_copy_per_lineage(n_extant=10, families=families, seed=seed)
        run = simulate_sequences(genomes, model=gtr(exchangeabilities, (0.25, 0.25, 0.25, 0.25)),
                                 length=length, substitution=substitution,
                                 profiles={f: rows for f in genomes.gene_trees}, seed=seed)
        for family in genomes.gene_trees:
            for ancestor, descendant, _ in _sequence_branches(tree, run, family, substitution):
                changed = ancestor != descendant
                np.add.at(changes, (ancestor[changed], descendant[changed]), 1)

    observed = (changes + changes.T)[upper]                 # changes between i and j, either way
    exposure = (profile[:, None] * profile[None, :])[upper]
    assert observed.sum() > 1500 and observed.min() > 50, (
        f"only {observed.sum():.0f} changes, rarest pair {observed.min():.0f} — too few to judge")

    def misfit(weights) -> float:
        expected = observed.sum() * weights / weights.sum()
        chi2 = float((((observed - expected) ** 2) / expected).sum())
        return (chi2 - (len(observed) - 1)) / math.sqrt(2 * (len(observed) - 1))

    assert abs(misfit(exposure * S[upper])) < Z_MAX, (
        f"changes under a profile do not follow π_i·π_j·S_ij with the exchangeabilities the model "
        f"was built from ({misfit(exposure * S[upper]):+.1f} standard errors) — the profile is not "
        f"keeping the base model's chemistry")
    assert misfit(exposure) > 20, (
        "equal exchangeabilities are not ruled out — a profile may be collapsing its model to F81")


def test_a_trait_driving_the_substitution_rate_is_integrated_across_the_branch():
    """The Traits→Sequences driver, which is exact rather than statistical.

    A branch of the species phylogram is the substitution rate **integrated** over the branch, so
    when a trait drives that rate and the trait switches part-way along, the branch length is
    ``base × Σ(factor of each state × how long it was held)`` — not the base rate times one factor
    sampled for the branch. Both give plausible phylograms; they differ by how much of the branch was
    spent in which state, which is precisely the quantity a conditioned run was set up to express.

    Because this is an integral of a step function and not a random draw, it can be checked exactly:
    the trait's own ``history`` gives the segments, and the phylogram gives what the engine used.
    They should agree to floating point. The two sampled-once alternatives — take the state the
    branch started in, take the one it ended in — are ruled out at the end, and the run is set up so
    that a good share of branches carry a switch, without which all three agree and the test would
    pass on a wiring that is wrong."""
    factors, base, switch, reps = {"hot": 3.0, "cold": 1.0}, 0.5, 1.2, 20
    label = re.compile(r"[ne](\d+):([0-9.eE+-]+)")
    checked = 0
    worst = 0.0
    #: for the branches the trait switched on, how far the two sampled-once wirings would land from
    #: what the engine actually wrote, relative to the branch
    off_if_sampled = {"the state it started in": [], "the state it ended in": []}

    for seed in range(reps):
        tree, genomes = _one_copy_per_lineage(n_extant=20, families=1, seed=seed)
        habitat = simulate_discrete(tree, states=list(factors), switch=switch, start="hot", seed=seed)
        run = simulate_sequences(genomes, model=jc69(), length=1,
                                 substitution=base * DrivenBy(habitat, factors), seed=seed)
        lengths = {int(i): float(v)
                   for i, v in label.findall(run.species_phylogram["complete"])}

        for node, length in lengths.items():
            segments = habitat.history[node]
            integrated = base * sum(factors[state] * held for state, held in segments)
            worst = max(worst, abs(length - integrated) / integrated)
            checked += 1
            if len(segments) > 1:                     # the trait switched part-way along this branch
                elapsed = tree.nodes[node].end_time - tree.nodes[node].birth_time
                for wiring, state in (("the state it started in", segments[0][0]),
                                      ("the state it ended in", segments[-1][0])):
                    sampled_once = base * factors[state] * elapsed
                    off_if_sampled[wiring].append(abs(sampled_once - length) / length)

    assert checked > 500, f"only {checked} branches to check"
    switched = len(off_if_sampled["the state it started in"])
    assert switched > 100, (
        f"only {switched} branches carried a switch — without them every wiring agrees and this test "
        f"proves nothing")
    assert worst < 1e-9, (
        f"a driven branch is off its integrated rate by {worst:.2e} relative — the driver is not "
        f"being integrated across the branch")

    for wiring, gaps in off_if_sampled.items():
        assert np.mean(gaps) > 0.1, (
            f"sampling the driver once per branch — taking {wiring} — would land within "
            f"{np.mean(gaps):.1%} of the integral on average, so this test would not tell the two "
            f"apart")


def test_across_site_rate_classes_are_applied_per_site_and_average_to_one():
    """``+I+Γ``: the sites of one gene run at different speeds, and the classes average to 1.

    With classes ``r_c`` in shares ``s_c``, a site differs across a branch of length ``d`` with
    probability ``Σ_c s_c · 3/4·(1 − e^(−4·d·r_c/3))`` — the Jukes–Cantor curve averaged over the
    classes, *not* evaluated at their mean. The difference is Jensen's inequality and it is the whole
    point of the model: because the curve saturates, spreading rates across sites gives **fewer**
    differences than one rate at the same average, which is why ignoring rate variation
    underestimates distance.

    So the single-rate expectation is the alternative ruled out here. It is what the run would
    produce if the classes were computed, normalised, reported on the model — and then never reached
    the sites. That is a failure mode with no other symptom: the invariant ``Σ r_c·s_c = 1`` still
    holds, the phylograms still read as mean substitutions per site, and the alignments still look
    like alignments.

    Classes are drawn per site, once per family, so the realised shares are multinomial rather than
    exact and the sites of a family are correlated across branches through their class. Replicates
    again, with the spread taken from them."""
    shape, invariant, substitution, reps = 0.5, 0.2, 0.8, 40
    model = jc69().across_sites(gamma_shape=shape, invariant=invariant)
    assert abs(sum(r * s for r, s in zip(model.site_rates, model.site_shares)) - 1.0) < 1e-12

    against_classes, against_one_rate = [], []
    for seed in range(reps):
        tree, genomes = _one_copy_per_lineage(n_extant=8, families=4, seed=seed)
        run = simulate_sequences(genomes, model=model, length=500, substitution=substitution,
                                 seed=seed)
        observed = mixture = single = 0.0
        for family in genomes.gene_trees:
            for ancestor, descendant, d in _sequence_branches(tree, run, family, substitution):
                sites = len(ancestor)
                observed += float((ancestor != descendant).sum())
                mixture += sites * sum(share * 0.75 * (1 - math.exp(-4 * d * rate / 3))
                                       for rate, share in zip(model.site_rates, model.site_shares))
                single += sites * 0.75 * (1 - math.exp(-4 * d / 3))
        against_classes.append(observed / mixture)
        against_one_rate.append(observed / single)

    assert abs(_z(np.array(against_classes), 1.0)) < Z_MAX, (
        f"observed differences are {np.mean(against_classes):.4f}× what the rate classes predict")
    assert abs(_z(np.array(against_one_rate), 1.0)) > 10, (
        "a single rate across sites is not ruled out — the classes may not be reaching the sites")


def test_the_lineage_clocks_are_mean_one_so_the_tree_is_not_inflated():
    """The relaxed clocks, at the one place they have been wrong before.

    A lineage clock multiplies the substitution rate by a random factor per species branch, and the
    factor has to have **mean 1** — otherwise every branch in the tree is systematically longer than
    the rate the run was given, and the whole phylogram inflates. Drawing ``exp(Normal(0, σ))``
    instead of ``exp(Normal(−σ²/2, σ))`` does exactly that, silently: the tree still looks like a
    tree, the alignments still look like alignments, and only the *scale* is wrong — by ``e^(σ²/2)``,
    which at σ = 0.5 is 13%. That is the historical lognormal-clock bug, and it is what the last
    assertion in each half rules out.

    The realised factor is recoverable from the run itself: a species-phylogram branch is
    ``base × factor × Δt``, so dividing out the base rate and the branch's time gives the factor the
    engine actually used. `ByLineage` draws one per branch independently, so the factors themselves
    are the sample; `FromParent` drifts parent→child, so the *ratios* down each branch are."""
    label = re.compile(r"[ne](\d+):([0-9.eE+-]+)")

    def realised_factors(modifier, base=1.0, reps=25):
        """``{seed: (tree, {node id: the clock factor the engine used on its branch})}``."""
        out = {}
        for seed in range(reps):
            tree, genomes = _one_copy_per_lineage(n_extant=25, families=1, seed=seed)
            run = simulate_sequences(genomes, model=jc69(), length=1, substitution=base * modifier,
                                     seed=seed)
            lengths = {int(i): float(v)
                       for i, v in label.findall(run.species_phylogram["complete"])}
            out[seed] = (tree, {i: length / (base * (tree.nodes[i].end_time - tree.nodes[i].birth_time))
                                for i, length in lengths.items()})
        return out

    spread = 0.5
    drawn = np.array([f for _, factors in realised_factors(ByLineage(spread=spread)).values()
                      for f in factors.values()])
    assert len(drawn) > 1000, "too few branches to judge the clock by"
    assert abs(_z(drawn, 1.0)) < Z_MAX, (
        f"ByLineage(spread={spread}) factors average {drawn.mean():.4f}, not 1 — every branch in the "
        f"phylogram is scaled by that")
    assert abs(np.log(drawn).std(ddof=1) - spread) < 0.05, (
        f"the log-scale spread is {np.log(drawn).std(ddof=1):.4f}, not the {spread} it was given")
    assert _z(drawn, math.exp(spread ** 2 / 2)) < -5, (
        "an uncorrected lognormal is not ruled out — the mean correction may be missing")

    spread = 0.4
    ratios = []
    for tree, factors in realised_factors(FromParent(spread=spread)).values():
        for i, factor in factors.items():
            parent = tree.nodes[i].parent
            if parent is not None:
                ratios.append(factor / factors[parent])
    ratios = np.array(ratios)
    assert abs(_z(ratios, 1.0)) < Z_MAX, (
        f"FromParent(spread={spread}) drifts by {ratios.mean():.4f} per branch on average, not 1 — "
        f"the rate ratchets down the tree")
    assert abs(np.log(ratios).std(ddof=1) - spread) < 0.05, (
        f"the log-scale spread of the drift is {np.log(ratios).std(ddof=1):.4f}, not {spread}")
    assert _z(ratios, math.exp(spread ** 2 / 2)) < -5, (
        "an uncorrected lognormal drift is not ruled out")


# --- the trait level: a value riding the tree -----------------------------------------------------

def _trait_tree():
    """One tree, reused so every expectation below is exact for *this* tree rather than averaged."""
    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=16, seed=1).complete_tree
    # a node's depth below is read straight off its end_time, which is only the elapsed time since
    # the origin if the origin is 0 — asserted rather than assumed, because a tree that started
    # elsewhere would shift every expectation here by a constant and still look entirely ordinary
    assert tree.nodes[tree.root].birth_time == 0.0
    return tree


def _mrca(tree, a: int, b: int) -> int:
    """The most recent node both ``a`` and ``b`` descend from."""
    above, node = set(), a
    while node is not None:
        above.add(node)
        node = tree.nodes[node].parent
    node = b
    while node not in above:
        node = tree.nodes[node].parent
    return node


def _trait_branches(tree, values, start):
    """``(ancestor value, node value, Δt)`` for every branch, the root's included.

    ``node_values`` is the value at each node's ``end_time``, and the root diffuses over its own
    branch from ``start`` at ``t = 0`` (convention B), so the root branch is a branch like any other
    and its ancestor is ``start``."""
    for node in tree.nodes.values():
        ancestor = start if node.parent is None else values[node.parent]
        yield ancestor, values[node.id], node.end_time - node.birth_time


def test_brownian_motion_moves_by_the_variance_it_was_given():
    """The diffusion, checked twice: once per branch, once across the tree.

    **Per branch** the increment is ``Normal(0, σ²·Δt)``, so dividing every increment by
    ``√(σ²·Δt)`` should give one standard normal sample pooled over branches of every length. Mean
    and variance are both checked — a σ² read as a standard deviation rather than a variance, or a
    ``Δt`` dropped from the draw, moves the variance while leaving the mean at 0.

    **Across the tree** the tips are what a user actually gets, and their law is Felsenstein's: two
    tips have variance ``σ²×`` their own depth and covariance ``σ²×`` their shared path, so the
    expected squared difference between them is ``σ²×`` the path length *between* them — down from
    one, up to the other, with everything above their MRCA cancelling. That is the statement that
    shared ancestry makes relatives similar, and it is the one a wrong tree traversal breaks while
    leaving every individual branch correct. It is measured as a squared separation rather than as a
    covariance because a covariance of near-independent tips is a badly conditioned thing to
    estimate, and the alternative it rules out — tips that ignore their shared ancestry — is the same
    either way."""
    sigma2, start, reps = 2.5, 0.0, 400
    tree = _trait_tree()
    runs = [simulate_continuous(tree, start=start, rate=sigma2, seed=s) for s in range(reps)]

    standardised = np.array([(value - ancestor) / math.sqrt(sigma2 * dt)
                             for run in runs
                             for ancestor, value, dt in _trait_branches(tree, run.node_values,
                                                                           start)])
    assert abs(_z(standardised, 0.0)) < Z_MAX, (
        f"branch increments average {standardised.mean():+.4f} rather than 0 — the diffusion drifts")
    variance_z = (standardised.var(ddof=1) - 1) / math.sqrt(2 / len(standardised))
    assert abs(variance_z) < Z_MAX, (
        f"branch increments have variance {standardised.var(ddof=1):.4f} against the σ²·Δt they were "
        f"drawn with, {variance_z:+.2f} standard errors away")

    tips = [n.id for n in tree.extant_leaves()]
    depth = {i: tree.nodes[i].end_time for i in tips}      # the root is born at t = 0
    pairs = [(a, b, depth[a] + depth[b] - 2 * tree.nodes[_mrca(tree, a, b)].end_time)
             for i, a in enumerate(tips) for b in tips[i + 1:]]

    separation = np.array([np.mean([(run.node_values[a] - run.node_values[b]) ** 2 / (sigma2 * path)
                                    for a, b, path in pairs]) for run in runs])
    assert abs(_z(separation, 1.0)) < Z_MAX, (
        f"tips sit {separation.mean():.4f}× as far apart as the path between them implies")

    ignoring_ancestry = np.array([np.mean([(run.node_values[a] - run.node_values[b]) ** 2
                                           / (sigma2 * (depth[a] + depth[b])) for a, b, _ in pairs])
                                  for run in runs])
    assert abs(_z(ignoring_ancestry, 1.0)) > 5, (
        "tips that ignored their shared ancestry are not ruled out — relatives may not be inheriting "
        "their common history")


def test_ornstein_uhlenbeck_realises_its_exact_transition_moments():
    """The OU pull, against the transition density it is supposed to have.

    Over a branch of length ``Δt`` from a value ``x``, an OU process ends at

        ``Normal(θ + (x − θ)·e^(−α·Δt),  σ²/(2α)·(1 − e^(−2α·Δt)))``

    — both moments depending on α, and neither reducing to the Brownian one except in the limit
    ``α → 0``. Standardising every branch by its own mean and variance turns the whole run into one
    standard normal sample, so a pull applied to the wrong quantity, or applied once per branch
    instead of continuously, moves the mean; a variance that forgets the ``(1 − e^(−2αΔt))`` factor
    — the part that keeps an OU from wandering off the way a Brownian motion does — moves the spread.

    The run starts far from the optimum (``start`` well below θ) on purpose: near θ the pull term is
    small and a broken pull would be nearly invisible. The alternative ruled out is exactly that —
    reading the run as if there were no pull at all."""
    sigma2, theta, pull, start, reps = 1.5, 4.0, 2.0, -3.0, 2000
    tree = _trait_tree()
    runs = [simulate_continuous(tree, start=start, rate=sigma2, reverts_to=theta, pull=pull, seed=s)
            for s in range(reps)]

    standardised, as_brownian = [], []
    for run in runs:
        for ancestor, value, dt in _trait_branches(tree, run.node_values, start):
            mean = theta + (ancestor - theta) * math.exp(-pull * dt)
            variance = sigma2 / (2 * pull) * (1 - math.exp(-2 * pull * dt))
            standardised.append((value - mean) / math.sqrt(variance))
            as_brownian.append((value - ancestor) / math.sqrt(sigma2 * dt))
    standardised = np.array(standardised)

    assert abs(_z(standardised, 0.0)) < Z_MAX, (
        f"standardised increments average {standardised.mean():+.4f} rather than 0 — the pull is not "
        f"landing where the OU transition density puts it")
    variance_z = (standardised.var(ddof=1) - 1) / math.sqrt(2 / len(standardised))
    assert abs(variance_z) < Z_MAX, (
        f"standardised increments have variance {standardised.var(ddof=1):.4f}, {variance_z:+.2f} "
        f"standard errors from the σ²/(2α)·(1−e^(−2αΔt)) they were drawn with")

    assert abs(_z(np.array(as_brownian), 0.0)) > 10, (
        "a run with no pull at all is not ruled out — reverts_to/pull may not be reaching the draw")


def test_an_mk_trait_switches_at_the_rate_each_direction_was_given():
    """The compensator identity at the trait level, one rate at a time.

    A lineage in state ``a`` leaves for ``b`` at ``q(a→b)`` per unit time, so over the whole run the
    expected number of ``a→b`` transitions is ``q(a→b) × (time spent in a)`` — and both quantities
    are in the result already: the events log has the transitions, and ``history`` has the
    ``(state, duration)`` segments that say how long each state was occupied. Nothing has to reach
    equilibrium for this to hold, which matters, because a tree this size does not.

    The rates are deliberately **asymmetric**, and each direction is checked against its own declared
    value. That is what catches a transposed rate dict — ``{"a->b": …}`` read as ``b→a`` — which is
    a mistake with no other symptom: the trait still switches, the log is still consistent, the map
    still derives, and only the direction of the asymmetry is reversed."""
    forward, backward, reps = 0.8, 0.3, 1200
    tree = _trait_tree()
    fired, occupancy = Counter(), Counter()

    for seed in range(reps):
        run = simulate_discrete(tree, states=["a", "b"], start="a",
                                switch={"a->b": forward, "b->a": backward}, seed=seed)
        for segments in run.history.values():
            for state, duration in segments:
                occupancy[state] += duration
        for change in run.events:
            if change.kind == "on_branch":
                fired[(change.from_state, change.to_state)] += 1

    for (frm, to), rate in ((("a", "b"), forward), (("b", "a"), backward)):
        count = fired[(frm, to)]
        assert count > 2000, f"only {count} {frm}→{to} transitions to judge by"
        realised = count / occupancy[frm]
        se = math.sqrt(count) / occupancy[frm]           # a Poisson count over a known exposure
        assert abs((realised - rate) / se) < Z_MAX, (
            f"{frm}→{to} fired at {realised:.4f} per unit time spent in {frm}, but was given {rate}")

    transposed = abs((fired[("a", "b")] / occupancy["a"] - backward)
                     / (math.sqrt(fired[("a", "b")]) / occupancy["a"]))
    assert transposed > 10, "a transposed rate dict is not ruled out"


def test_an_mk_branch_ends_in_the_other_state_as_often_as_the_chain_says():
    """The two-state chain's transition probability, which pins the *shape* the rate alone does not.

    For a symmetric two-state chain at rate ``q``, a branch of length ``Δt`` ends in the other state
    with probability ``(1 − e^(−2qΔt))/2`` — saturating at ½ however long the branch, because the
    chain forgets where it started. Counting transitions (the test above) pins the rate; this pins
    what the chain *does* with it, and rules out the naive alternative that treats every switch as
    permanent, ``q·Δt``, which has the same slope at zero and no saturation at all.

    As with JC69, the end state's disagreement with the start does not depend on which state that
    was, so branches are independent and the binomial standard error applies despite the tree."""
    switch, reps = 0.9, 150
    tree = _trait_tree()
    observed = expected = naive = variance = 0.0

    for seed in range(reps):
        run = simulate_discrete(tree, states=["a", "b"], switch=switch, start="a", seed=seed)
        for ancestor, value, dt in _trait_branches(tree, run.node_values, "a"):
            observed += ancestor != value
            p = (1 - math.exp(-2 * switch * dt)) / 2
            expected += p
            variance += p * (1 - p)
            naive += min(switch * dt, 1.0)

    assert observed > 500, "too few branches ended elsewhere to say anything"
    assert abs((observed - expected) / math.sqrt(variance)) < Z_MAX, (
        f"{observed:.0f} branches ended in the other state, against the chain's {expected:.1f}")
    assert abs((observed - naive) / math.sqrt(variance)) > 5, (
        "a chain whose switches never reverse is not ruled out")


def test_correlated_traits_realise_the_correlation_they_were_given():
    """The ``correlation=`` overlay, which is the whole reason to evolve two traits in one call.

    Two traits diffusing jointly take their branch increment from ``MVN(0, Σ·Δt)`` with
    ``Σ = D·R·D``, so after dividing each trait's increment by its own ``√(σ²·Δt)`` the pair has
    correlation exactly ρ — whatever the branch lengths, and whatever the two σ² are. The σ² here
    differ by a factor of four on purpose: a Σ assembled without ``D`` on both sides would still give
    a correlated-looking pair, just not at ρ.

    ρ = 0 is the alternative ruled out, and it is the one that matters. A comparative method is
    usually being asked whether an apparent association between two traits is real, and the null it
    is graded against is a run with the correlation switched off. If the overlay were not reaching
    the draw, both the signal and the null would be the same run and the grading would be
    meaningless."""
    rho, sigma_x, sigma_y, reps = 0.6, 1.0, 4.0, 60
    tree = _trait_tree()
    increments = []

    for seed in range(reps):
        run = simulate_continuous(tree, start={"x": 0.0, "y": 0.0},
                                  rate={"x": sigma_x, "y": sigma_y},
                                  correlation={("x", "y"): rho}, seed=seed)
        origin = {"x": 0.0, "y": 0.0}
        for ancestor, value, dt in _trait_branches(tree, run.node_values, origin):
            increments.append(((value["x"] - ancestor["x"]) / math.sqrt(sigma_x * dt),
                               (value["y"] - ancestor["y"]) / math.sqrt(sigma_y * dt)))

    increments = np.array(increments)
    realised = float(np.corrcoef(increments.T)[0, 1])
    se = (1 - rho ** 2) / math.sqrt(len(increments))
    assert abs((realised - rho) / se) < Z_MAX, (
        f"the two traits came out correlated at {realised:.4f}, not the ρ = {rho} they were given")
    assert abs(realised / (1 / math.sqrt(len(increments)))) > 10, (
        "uncorrelated traits are not ruled out — the overlay may not be reaching the draw")


def test_the_manual_modifier_table_matches_what_the_engines_wire():
    """Appendix A's "which level accepts which" table, checked against the engines it describes.

    It had drifted: `ByFamily` was missing from the genome row and `FromParent` from the sequence
    row, so a reader who trusted the appendix would think two working things were unsupported. A
    tester found the contradiction against `--help`, trusted the CLI, and was right — but said the
    appendix had lost her confidence, which is the real cost of a doc that disagrees with the code."""
    import pathlib
    import re

    from zombi2.genomes import IMPLEMENTED_MODIFIERS as GENOMES
    from zombi2.genomes.nucleotide import IMPLEMENTED_MODIFIERS as NUCLEOTIDE
    from zombi2.genomes.ordered import IMPLEMENTED_MODIFIERS as ORDERED
    from zombi2.joint import IMPLEMENTED_MODIFIERS as JOINT
    from zombi2.sequences import IMPLEMENTED_MODIFIERS as SEQUENCES
    from zombi2.species import IMPLEMENTED_MODIFIERS as SPECIES
    from zombi2.traits import IMPLEMENTED_MODIFIERS as TRAITS
    from zombi2.traits.discrete import IMPLEMENTED_MODIFIERS as TRAITS_DISCRETE

    appendix = pathlib.Path(__file__).resolve().parent.parent / "manual" / "book" / "appendix-a.md"
    if not appendix.exists():                       # the manual is not shipped in every checkout
        import pytest
        pytest.skip("manual/book/appendix-a.md not present")
    # scope the search to the modifier table — an earlier table in the same appendix lists scopes
    # and also has a "| Species |" row
    text = appendix.read_text(encoding="utf-8").split("### Which level accepts which", 1)[-1]

    # One row covers two engines — "Genomes, family and ordered" — and only family's tuple was ever
    # checked against it, so the row was free to be false for ordered, and for a while it was: it
    # listed DrivenBy while the ordered engine refused it. Assert the two agree, and the shared row
    # cannot lie again about either.
    assert set(ORDERED) == set(GENOMES), (
        f"appendix A lists family and ordered under one row, but they wire "
        f"{sorted(m.__name__ for m in GENOMES)} and {sorted(m.__name__ for m in ORDERED)} — either "
        f"bring them back into line or split the row")

    for row, wired in (("Species", SPECIES),
                       ("Genomes, family and ordered", GENOMES),
                       ("Genomes, nucleotide", NUCLEOTIDE),
                       ("Sequences", SEQUENCES),
                       ("Traits, continuous `rate`", TRAITS),
                       ("Traits, discrete `switch`", TRAITS_DISCRETE),
                       ("Joint, `birth` / `death`", JOINT)):
        line = next((ln for ln in text.splitlines() if ln.startswith(f"| {row} |")), None)
        assert line, f"appendix A has no row for {row!r}"
        # the SECOND column only: a row label may itself carry a backticked word (the rate a row is
        # about — `rate`, `switch`, `birth`), and scanning the whole line would read those as
        # modifiers.
        listed = set(re.findall(r"`(\w+)`", line.split("|")[2]))
        assert listed == {m.__name__ for m in wired}, (
            f"appendix A's {row!r} row lists {sorted(listed)}, but the engine wires "
            f"{sorted(m.__name__ for m in wired)}")
