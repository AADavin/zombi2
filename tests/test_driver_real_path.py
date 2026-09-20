"""A continuous driver is read along its real path, not the straight line (issue #454).

ZOMBI2 records a continuous trait only at the nodes, and a reader used to take the straight line
between a branch's two node values. That line is the mean of the bridge between them with the
excursions dropped, so under a non-linear mapping it is biased, and a smaller reading step does not
remove the bias. The path is now drawn instead — after the run, off its own random stream — and
written to ``trait_path.tsv`` so a file-driven run reads the same path an in-memory one does.

These tests check the law against closed forms rather than against itself: the Brownian bridge's
pointwise variance and its integral, the Ornstein–Uhlenbeck bridge's ``sinh`` form, the exact
integrated rate under a log link, and a brute-force fine simulation under a ``changing_at``
schedule.
"""

import math
import warnings

import numpy as np
import pytest

from zombi2 import traits
from zombi2.params import PerLineage
from zombi2.params.conditioned import (branch_starts, driver_from_continuous_result,
                                       load_driver, path_points)
from zombi2.species import simulate_species_tree
from zombi2.tree import Node, Tree


def one_branch(length: float = 1.0) -> Tree:
    """A tree that is a single branch from 0 to ``length`` — the shape every closed form here is
    written for."""
    return Tree({0: Node(0, None, 0.0, length, (), "extant")}, 0)


def pinned(tree, *, value: float = 0.0, rate=1.0, pull=None, reverts_to=None, seed: int = 1):
    """A continuous run on ``tree`` whose node values are then forced to ``value``.

    The bridge law does not read the node values — it is the rate, the pull and the optimum — so
    replacing them afterwards gives a branch with known endpoints under the real law, which is what
    the closed forms below need. ``start`` is pinned too, through the ``initial`` event."""
    result = traits.simulate_continuous(tree, start=value, rate=rate, pull=pull,
                                        reverts_to=reverts_to, seed=seed)
    result.node_values = {i: float(value) for i in result.node_values}
    return result


def midpoint_values(result, *, step, node_id=0, seeds=range(400)):
    """The value each of many independent draws puts at the branch's midpoint.

    Independent draws come from independent runs: the path is keyed by the run's seed, so the seed
    is what is varied. ``step`` must divide the branch into an odd number of stretches for the
    branch midpoint to be one of them."""
    out = []
    for s in seeds:
        result.seed = int(s)
        traj = driver_from_continuous_result(result, step=step)
        node = result.complete_tree.nodes[node_id]
        mid = 0.5 * (node.birth_time + node.end_time)
        out.append(float(traj.value(node_id, mid)))
    return np.array(out)


# --- 1. the feature changes no node value --------------------------------------------------------

def test_node_values_are_untouched_by_the_path():
    """The path is drawn after the run, off its own stream, so the engine's draw order is the same
    and every node value is the number it was."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=7).complete_tree
    a = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=11)
    before = dict(a.node_values)

    traj = driver_from_continuous_result(a, step=0.05)       # draws every branch's path
    assert traj.value(ct.root, ct.nodes[ct.root].birth_time) is not None
    assert a.node_values == before

    b = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=11)
    assert b.node_values == before                            # the same seed, element for element

    # and an OU run with a modified σ², which walks a different code path
    rate = PerLineage(1.0).changing_at({0: 1.0, 1.0: 4.0})
    c = traits.simulate_continuous(ct, start=0.0, rate=rate, reverts_to=2.0, pull=0.5, seed=13)
    kept = dict(c.node_values)
    driver_from_continuous_result(c, step=0.05)
    assert c.node_values == kept


# --- 2. Brownian motion: the exact pointwise variance and the exact ∫X² ---------------------------

def test_brownian_bridge_matches_its_closed_form():
    """One branch, T=1, σ²=1, both ends 0. The bridge's variance at the midpoint is σ²T/4 = 0.25 and
    E[∫X² dt] over the branch is σ²T²/6 = 1/6."""
    result = pinned(one_branch(1.0), rate=1.0)
    mids = midpoint_values(result, step=1.0)                  # one stretch: its midpoint is t=0.5
    assert len(mids) == 400
    var = float(mids.var())
    assert var == pytest.approx(0.25, rel=0.15), var          # n=400 ⇒ ~7% s.e. on a variance

    # E[∫ X² dt] = ∫ Var(X_t) dt = σ²T²/6. Read on a fine grid, so the stretch values are the path.
    step = 1 / 64
    total = []
    for s in range(400):
        result.seed = s
        traj = driver_from_continuous_result(result, step=step)
        xs = [traj.value(0, (k + 0.5) * step) for k in range(64)]
        total.append(sum(x * x for x in xs) * step)
    got = float(np.mean(total))
    assert got == pytest.approx(1 / 6, rel=0.12), got


# --- 3. Ornstein–Uhlenbeck: the sinh closed form --------------------------------------------------

def test_ou_bridge_matches_the_sinh_closed_form():
    """The OU bridge's variance at an interior time is σ²·sinh(αs)·sinh(αu)/(α·sinh(αT)).

    Checked twice: the general V()-based formula agrees with that closed form to 1e-12, and sampled
    paths agree with it to sampling error."""
    alpha, sigma2, T = 1.3, 2.0, 1.0
    tree = one_branch(T)
    result = pinned(tree, rate=sigma2, pull=alpha, reverts_to=0.0)
    law = result.path_law
    assert law is not None and law.pull == alpha

    def closed(s, u):
        return sigma2 * math.sinh(alpha * s) * math.sinh(alpha * u) / (alpha * math.sinh(alpha * (s + u)))

    # the general formula, from V() alone — var = V(0,t) − C²/V(0,T)
    for t in (0.1, 0.25, 0.5, 0.75, 0.9):
        v_t = law.variance(0, 0.0, t)
        v_full = law.variance(0, 0.0, T)
        c = math.exp(-alpha * (T - t)) * v_t
        assert v_t - c * c / v_full == pytest.approx(closed(t, T - t), abs=1e-12)

    mids = midpoint_values(result, step=T)
    assert float(mids.var()) == pytest.approx(closed(0.5, 0.5), rel=0.15), float(mids.var())
    assert abs(float(mids.mean())) < 0.15                     # both ends at θ=0, so the bridge is centred


# --- 4. the bias the issue is about ---------------------------------------------------------------

def test_a_log_link_is_biased_on_the_line_and_right_on_the_path():
    """One branch, T=1, σ²=1, both node values 0, rate = exp(X).

    Along the straight line X is 0 throughout, so the integrated rate is exactly 1.0. Along the real
    path E[exp(X_t)] = exp(σ²t(1−t)/2T), whose integral over [0,1] is 1.087653. The line is wrong by
    8.8%, halving the step does not move it, and the drawn path lands on the closed form."""
    exact = 1.0876530389043015
    result = pinned(one_branch(1.0), rate=1.0)

    def integrated(step, *, path, seed=0):
        result.seed = seed
        traj = driver_from_continuous_result(result, step=step, path=path)
        n = round(1.0 / step)
        return sum(math.exp(traj.value(0, (k + 0.5) * step)) for k in range(n)) * step

    line = integrated(1 / 64, path=False)
    assert line == pytest.approx(1.0, abs=1e-12)
    assert abs(line - exact) / exact > 0.05                   # the line is wrong by more than 5%
    finer = integrated(1 / 128, path=False)
    assert finer == pytest.approx(1.0, abs=1e-12)             # …and halving the step does not help
    assert abs(finer - exact) / exact > 0.05

    drawn = float(np.mean([integrated(1 / 64, path=True, seed=s) for s in range(600)]))
    assert drawn == pytest.approx(exact, rel=0.02), (drawn, exact)


# --- 5. the same path every time ------------------------------------------------------------------

def test_a_path_is_the_same_however_it_is_asked_for():
    """The path is keyed by the run's seed and the branch, so resolving twice gives the same
    trajectory and resolving another driver first changes nothing."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=12, seed=3).complete_tree
    a = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=21)
    b = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=22)

    first = driver_from_continuous_result(a, step=0.05)
    again = driver_from_continuous_result(a, step=0.05)
    driver_from_continuous_result(b, step=0.05)               # another driver in between
    third = driver_from_continuous_result(a, step=0.05)

    for i in ct.nodes:
        assert first._states[i] == again._states[i]
        assert first._states[i] == third._states[i]
    assert any(first._states[i] != driver_from_continuous_result(b, step=0.05)._states[i]
               for i in ct.nodes)                             # two runs are two paths


# --- 6. the at_speciation endpoint --------------------------------------------------------------

def test_a_branch_starts_from_the_post_jump_value():
    """With ``at_speciation`` the engine jumps at the split and then diffuses from the POST-jump
    value. The branch's left endpoint is that value, not the parent's node value — reading the
    parent's pulled every branch of such a run toward the parent."""
    ct = simulate_species_tree(birth=1.0, n_extant=8, seed=5).complete_tree
    result = traits.simulate_continuous(ct, start=0.0, rate=0.5, at_speciation=4.0, seed=9)
    jumps = {e.lineage: (e.from_state, e.to_state)
             for e in result.events if e.kind == "on_speciation"}
    assert jumps, "this run is meant to have on-speciation jumps"

    points = path_points(ct, result.node_values, 1e9,          # one stretch: endpoints only
                         law=None, seed=None, starts=branch_starts(result))
    for i, (pre, post) in jumps.items():
        if ct.nodes[i].end_time <= ct.nodes[i].birth_time:
            continue
        assert points[i][0][1] == pytest.approx(float(post))   # the branch starts after the jump
        parent_value = float(result.node_values[ct.nodes[i].parent])
        assert parent_value == pytest.approx(float(pre))        # the jump's `from` IS the parent
        if abs(float(post) - parent_value) > 1e-9:
            assert points[i][0][1] != pytest.approx(parent_value)   # …and NOT from the parent's value


def test_the_root_branch_starts_from_the_initial_value():
    """The root's left endpoint is the ``initial`` row, the trait at t=0. ``node_values[root]`` is
    the value at the first split, so the root branch has no earlier node to read."""
    ct = simulate_species_tree(birth=1.0, n_extant=6, seed=5).complete_tree
    result = traits.simulate_continuous(ct, start=3.5, rate=1.0, seed=9)
    assert branch_starts(result)[ct.root] == pytest.approx(3.5)


# --- 7. a modified σ² ----------------------------------------------------------------------------

def test_a_changing_sigma2_still_gives_the_exact_bridge():
    """σ² steps from 1 to 9 halfway along the branch. The bridge is built on the accrued variance,
    so it stays exact; compared against a brute-force fine-grained forward simulation conditioned on
    the same two endpoints, to a tolerance of 8% on the variance."""
    T = 1.0
    rate = PerLineage(1.0).changing_at({0: 1.0, 0.5: 9.0})
    result = pinned(one_branch(T), rate=rate)

    # V(0, t) the bridge uses: ∫σ², i.e. t for t ≤ 0.5 and 0.5 + 9(t − 0.5) after
    law = result.path_law
    assert law.variance(0, 0.0, 0.25) == pytest.approx(0.25)
    assert law.variance(0, 0.0, 1.0) == pytest.approx(0.5 + 4.5)

    # brute force: simulate the unconditioned process on a fine grid, keep the paths that end near 0
    rng = np.random.default_rng(0)
    n, keep = 4000, []
    grid = np.linspace(0.0, T, n + 1)
    sig2 = np.where(grid[:-1] < 0.5, 1.0, 9.0)
    for _ in range(120_000):
        steps = rng.normal(0.0, np.sqrt(sig2 * (T / n)))
        walk = np.concatenate([[0.0], np.cumsum(steps)])
        if abs(walk[-1]) < 0.08:                               # condition on ending at 0
            keep.append(walk[n // 2])
    assert len(keep) > 2000, len(keep)
    brute = float(np.var(keep))

    mids = midpoint_values(result, step=T, seeds=range(3000))
    assert float(mids.var()) == pytest.approx(brute, rel=0.08), (float(mids.var()), brute)


# --- the written path ----------------------------------------------------------------------------

def written(tmp_path, result, *, step=None):
    """Write ``result`` and return the directory — values plus path, the default pair."""
    result.write(tmp_path, step=step)
    return tmp_path


def test_the_written_path_round_trips(tmp_path):
    """``trait_path.tsv`` is written at full precision, so what is read back equals what was
    written, bit for bit."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=4).complete_tree
    result = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=8)
    step = 0.05
    written(tmp_path, result, step=step)

    from zombi2.params.conditioned import _load_path_file
    law, starts, _ = _load_path_file(tmp_path / "trait_values.tsv")
    wanted = path_points(ct, result.node_values, step, law=result.path_law, seed=result.seed,
                         starts=branch_starts(result))
    assert set(law.points) == set(wanted)
    for i in wanted:
        assert law.points[i] == wanted[i]                      # exact, not approx


def test_a_file_driver_reads_the_same_path_as_the_in_memory_one(tmp_path):
    """The point of the file. Written and read back at the step it was written at, a file-backed
    driver gives the same trajectory, value for value, as the driver used in memory."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=6).complete_tree
    result = traits.simulate_continuous(ct, start=0.0, rate=1.0, at_speciation=0.5, seed=12)
    step = 0.05
    written(tmp_path, result, step=step)

    memory = driver_from_continuous_result(result, step=step)
    from_file = load_driver(tmp_path / "trait_values.tsv", ct, step=step)
    for i in ct.nodes:
        assert from_file._starts[i] == pytest.approx(memory._starts[i])
        for a, b in zip(from_file._states[i], memory._states[i]):
            assert a == pytest.approx(b, abs=1e-12), i


def test_a_finer_reader_agrees_with_the_written_path(tmp_path):
    """A reader at half the written step refines between written points rather than redrawing, so
    at the times the two share — the written points themselves — it gives the written values."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=6).complete_tree
    result = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=15)
    step = 0.05
    written(tmp_path, result, step=step)

    from zombi2.params.conditioned import _load_path_file
    law, _, seed = _load_path_file(tmp_path / "trait_values.tsv")
    rng = np.random.default_rng(0)
    for i, pts in law.points.items():
        times = [t for t, _, _ in pts]
        got = law.sample(i, times[0], times[-1], pts[0][1], pts[-1][1], times, rng)
        for (t, v, _), (x, _) in zip(pts, got):
            assert x == pytest.approx(v, abs=1e-12), (i, t)


def test_a_values_file_with_no_path_file_falls_back_and_says_so(tmp_path):
    """An older run's output has no ``trait_path.tsv``. The reader then takes the straight line, and
    warns, naming the file that is missing."""
    ct = simulate_species_tree(birth=1.0, n_extant=6, seed=6).complete_tree
    result = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=16)
    result.write(tmp_path, outputs=("values",))
    assert not (tmp_path / "trait_path.tsv").exists()

    with pytest.warns(RuntimeWarning, match="trait_path.tsv"):
        traj = load_driver(tmp_path / "trait_values.tsv", ct, step=0.05)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        line = driver_from_continuous_result(result, step=0.05, path=False)
    # every branch but the root: the values file carries the value at each node, and the line
    # between two of them is what both readers take. The ROOT is the exception — its left endpoint
    # is the trait at t=0, which is in the run and in trait_path.tsv but in no values file — so a
    # directory with neither can only hold the root branch at its node value.
    for i in ct.nodes:
        if i == ct.root:
            continue
        for a, b in zip(traj._states[i], line._states[i]):
            assert a == pytest.approx(b, rel=1e-5)


# --- the written form round-trips ------------------------------------------------------------------

def test_path_false_round_trips_through_the_written_form():
    """``path=False`` is a different model from the default, so the written form records it and
    reparses to the same rate."""
    from zombi2.params.parse import parse_rate

    r = PerLineage(1.0).scaled_by("t.tsv", {"a": 2.0}, path=False)
    assert ", path=False" in repr(r)
    assert parse_rate(repr(r)) == r
    d = PerLineage(1.0).scaled_by("t.tsv", {"a": 2.0})
    assert "path=" not in repr(d)
    assert parse_rate(repr(d)) == d
    assert d != r
