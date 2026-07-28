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
- **a driven rate realises the multiplier it was declared with** — on the conditioned path and on the
  joint one. These matter most and are checked last, because a mis-wired coupling is the one error
  that leaves every output well formed: a tree, a genome and a log that are all internally consistent,
  with only the strength of the association wrong, which is precisely what the run was made to
  measure.

**Every test is deterministic.** The seeds are fixed, so a "statistical" test here cannot flake: it
computes one number and compares it to one expectation. The tolerances are wide (``|z| < 4``) because
their job is to catch a model that has changed, not to detect a rate that is off by a percent — and
the hypotheses these rule out are rejected by tens of standard deviations, not by twos.

What that buys, measured by mutating the model rather than assumed. Scaling a rate by 1.15 behind the
engine's back fails these; by 1.02 it does not. Biasing 10% of inversions toward one half of the ring
fails; 3% does not. Mis-wiring the coupling strength by 5% is caught on the conditioned path and by
10% on the joint one. So resolution sits between a few and ten percent, depending on the check.

That figure is a property of *these* trees and replicate counts, not of the method: more material, or
more replicates, tightens it. Which is the honest way to read this file — it is evidence that the
model is the one advertised, not a calibration of it. A rate quietly off by a percent still gets
past, and catching that costs more replicates than a test suite should spend on every commit.

If one of these fails, the arithmetic below is not what to doubt first. Something about what the
engine samples has changed.
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np

from zombi2.genomes import simulate_genomes_family, simulate_genomes_ordered
from zombi2.genomes.ordered import Inversion
from zombi2.joint import simulate_joint
from zombi2.rates.modifiers import DrivenBy
from zombi2.species import simulate_species_tree
from zombi2.traits import discrete, simulate_discrete

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
        for e in simulate_genomes_family(tree, initial_families=30, transfer=0.6, seed=s).events:
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


# --- the coupled models: a rate that reads another level ------------------------------------------

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
    """The coupling mechanism itself, on the conditioned path (the driver is a finished level).

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
        for e in run.events:
            if e.kind == "origination":
                seen[_state_at(habitat.history, e.lineage, tree, e.time)] += 1

    for state, factor in factors.items():
        realised = seen[state] / lineage_time[state]
        expected = base * factor
        se = math.sqrt(seen[state]) / lineage_time[state]      # Poisson count, known exposure
        assert abs((realised - expected) / se) < Z_MAX, (
            f"{state}: origination ran at {realised:.4f} per unit time, but was told to run at "
            f"base × {factor} = {expected:.4f}")


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


def test_the_manual_modifier_table_matches_what_the_engines_wire():
    """Appendix A's "which level accepts which" table, checked against the engines it describes.

    It had drifted: `ByFamily` was missing from the genome row and `FromParent` from the sequence
    row, so a reader who trusted the appendix would think two working things were unsupported. A
    tester found the contradiction against `--help`, trusted the CLI, and was right — but said the
    appendix had lost her confidence, which is the real cost of a doc that disagrees with the code."""
    import pathlib
    import re

    from zombi2.genomes import WIRED_MODIFIERS as GENOMES
    from zombi2.genomes.nucleotide import WIRED_MODIFIERS as NUCLEOTIDE
    from zombi2.sequences import WIRED_MODIFIERS as SEQUENCES
    from zombi2.species import WIRED_MODIFIERS as SPECIES
    from zombi2.traits import WIRED_MODIFIERS as TRAITS

    appendix = pathlib.Path(__file__).resolve().parent.parent / "manual" / "book" / "appendix-a.md"
    if not appendix.exists():                       # the manual is not shipped in every checkout
        import pytest
        pytest.skip("manual/book/appendix-a.md not present")
    # scope the search to the modifier table — an earlier table in the same appendix lists scopes
    # and also has a "| Species |" row
    text = appendix.read_text().split("### Which level accepts which", 1)[-1]

    for row, wired in (("Species", SPECIES),
                       ("Genomes — family, ordered", GENOMES),
                       ("Genomes — nucleotide", NUCLEOTIDE),
                       ("Sequences", SEQUENCES),
                       ("Traits", TRAITS)):
        line = next((ln for ln in text.splitlines() if ln.startswith(f"| {row} |")), None)
        assert line, f"appendix A has no row for {row!r}"
        listed = set(re.findall(r"`(\w+)`", line))
        assert listed == {m.__name__ for m in wired}, (
            f"appendix A's {row!r} row lists {sorted(listed)}, but the engine wires "
            f"{sorted(m.__name__ for m in wired)}")
