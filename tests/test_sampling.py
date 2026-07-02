"""The Fenwick tree used for O(log n) weighted event selection."""

import numpy as np

from zombi2._sampling import Fenwick


def test_fenwick_total_and_find_match_brute_force():
    rng = np.random.default_rng(0)
    for _ in range(200):
        n = int(rng.integers(1, 40))
        vals = np.zeros(n)
        f = Fenwick(n)
        for _ in range(2 * n):  # random point assignments
            i = int(rng.integers(n))
            v = float(rng.random() * 5)
            vals[i] = v
            f.set(i, v)
        assert abs(f.total - vals.sum()) < 1e-9
        if vals.sum() <= 0:
            continue
        cum = np.cumsum(vals)
        for _ in range(50):
            r = float((1.0 - rng.random()) * vals.sum())   # in (0, total]
            got = f.find(r)
            assert got == int(np.searchsorted(cum, r, side="left"))
            assert vals[got] > 0                            # always a live slot


def test_fenwick_set_can_zero_a_slot():
    f = Fenwick(5)
    for i, v in enumerate([1.0, 2.0, 3.0, 4.0, 5.0]):
        f.set(i, v)
    assert f.total == 15.0
    f.set(2, 0.0)                 # "remove" a branch
    assert f.total == 12.0
    # a draw can never land on the zeroed slot
    assert all(f.find(r) != 2 for r in np.linspace(0.01, 12.0, 200))
