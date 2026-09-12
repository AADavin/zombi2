"""`WeightedIndex` — one weight per living lineage, their total and a draw from them (issue #436).

A step adds the living lineages' weights up to race the events and an event draws one of them by
weight, and both used to walk the lineages, which is what made a driven run grow with the tree. The
weights are the leaves of a binary tree now. What has to hold is that the tree says what the weights
say: the total is their sum, the draw is the draw a scan of them would make, and both survive
lineages entering and being retired.
"""

import numpy as np
import pytest

from zombi2._runtime.draw import WeightedIndex, weighted_index


def _of(values):
    w = WeightedIndex()
    for v in values:
        w.append(float(v))
    return w


def test_the_total_is_the_weights_summed():
    w = _of([3.0, 1.0, 4.0, 1.0, 5.0])
    assert w.total == pytest.approx(14.0)
    w.set(2, 0.5)
    assert w.total == pytest.approx(10.5)


def test_the_weights_read_back_as_they_were_set():
    w = _of(range(1, 40))
    w.set(7, 100.0)
    assert list(w) == [100.0 if i == 7 else float(i + 1) for i in range(39)]
    assert len(w) == 39 and w[7] == 100.0


def test_retiring_moves_the_last_lineage_into_the_gap():
    """The swap-remove `zombi2.genomes._live.retire` makes, so the two stay in step."""
    w = _of([1.0, 2.0, 3.0, 4.0])
    w.remove(1)
    assert list(w) == [1.0, 4.0, 3.0] and w.total == pytest.approx(8.0)
    w.remove(2)                                  # the last one: nothing moves
    assert list(w) == [1.0, 4.0] and w.total == pytest.approx(5.0)


def test_the_total_does_not_drift_over_many_changes():
    """It is recomputed from the weights rather than added to and subtracted from, so it stays the
    sum of what is there however long the run."""
    rng = np.random.default_rng(4)
    w = _of(rng.random(64))
    for _ in range(20000):
        w.set(int(rng.integers(len(w))), float(rng.random()))
    assert w.total == pytest.approx(sum(w), rel=0, abs=1e-12)


@pytest.mark.parametrize("n", [1, 2, 3, 17, 64, 100])
def test_the_draw_is_the_draw_a_scan_would_make(n):
    """The same random value picks the same lineage either way, so replacing the scan changes which
    lineage acts in no run — only the total, which a tree adds in pairs."""
    rng = np.random.default_rng(n)
    values = [float(x) for x in np.round(rng.random(n) * 10, 6)]
    w = _of(values)
    for seed in range(200):
        assert w.pick(np.random.default_rng(seed)) == weighted_index(
            np.random.default_rng(seed), values, sum(values))


def test_a_draw_lands_on_each_lineage_in_proportion_to_its_weight():
    rng = np.random.default_rng(11)
    values = [1.0, 3.0, 0.0, 6.0]
    w = _of(values)
    hits = [0, 0, 0, 0]
    for _ in range(40000):
        hits[w.pick(rng)] += 1
    assert hits[2] == 0                          # a lineage weighing nothing is never drawn
    for hit, value in zip(hits, values):
        assert hit / 40000 == pytest.approx(value / sum(values), abs=0.01)


def test_it_grows_past_its_first_capacity():
    """The tree doubles as lineages enter, and the weights already in it come across."""
    w = _of(range(1, 130))
    assert len(w) == 129 and w.total == pytest.approx(129 * 130 / 2)
    assert list(w) == [float(i) for i in range(1, 130)]


# --- the same tree read as counts: which lineage holds the j-th element ---------------------------
# A uniform pick over the whole gene or chromosome pool used to walk the lineages subtracting each
# one's count. `find` answers the same question from the tree, and counts are whole numbers, so it
# reaches the element the walk reached.

def _walk(counts, j):
    """The walk `find` replaces: the lineage the j-th element falls in, and its offset inside it."""
    for k, c in enumerate(counts):
        if j < c:
            return k, j
        j -= c
    raise AssertionError("j is past the end")


@pytest.mark.parametrize("counts", [
    [5], [0, 3, 0, 7], [1, 1, 1, 1, 1], [12, 0, 0, 40, 3, 9, 1], list(range(1, 30)),
], ids=["one", "with empty ones", "all equal", "mixed", "many"])
def test_find_reaches_the_element_the_walk_reached(counts):
    w = _of(counts)
    for j in range(sum(counts)):
        assert w.find(float(j)) == _walk(counts, j)


def test_find_on_a_pool_of_one_lineage_is_the_offset_itself():
    w = _of([9])
    assert [w.find(float(j)) for j in range(9)] == [(0, float(j)) for j in range(9)]
