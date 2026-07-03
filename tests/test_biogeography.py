"""Tests for the DEC biogeographic model (``zombi2.biogeography``)."""

import numpy as np
import pytest

import zombi2 as z


def _tree(n_tips=20, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n_tips, age=5.0, seed=seed)


# --------------------------------------------------------------------------- ranges & rates
def test_range_enumeration():
    d = z.DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1)
    assert len(d.states) == 7                                    # all non-empty subsets of 3
    assert d.states[:3] == [("A",), ("B",), ("C",)]             # ordered by size
    capped = z.DEC(areas=3, dispersal=0.1, extinction=0.1, max_range_size=2)
    assert len(capped.states) == 6                              # full 3-area range excluded
    assert all(len(s) <= 2 for s in capped.states)


def test_dispersal_and_extinction_rates():
    d = z.DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.2)
    Q = d.Q
    a, ab, abc, b = (d.encode({"A"}), d.encode({"A", "B"}),
                     d.encode({"A", "B", "C"}), d.encode({"B"}))
    assert Q[a, ab] == pytest.approx(0.1)                       # gain B: D[A,B]
    assert Q[ab, a] == pytest.approx(0.2)                       # lose B: E[B]
    assert Q[ab, b] == pytest.approx(0.2)                       # lose A: E[A]
    assert Q[ab, abc] == pytest.approx(0.2)                     # gain C: D[A,C]+D[B,C]
    # a single-area range cannot go locally extinct (no shrinking transitions)
    assert Q[a, a] == pytest.approx(-0.2)                       # total out-rate = 2 dispersals


def test_max_range_size_blocks_expansion():
    d = z.DEC(areas=3, dispersal=0.5, extinction=0.1, max_range_size=2)
    assert max(len(s) for s in d.states) == 2
    ab = d.encode({0, 1})
    # a size-2 range is at the cap: no dispersal (would exceed it), only extinction transitions
    for j, s in enumerate(d.states):
        if d.Q[ab, j] > 0:
            assert len(s) == 1                                 # only shrinking is allowed


# --------------------------------------------------------------------------- cladogenesis
def test_cladogenesis_single_area_is_narrow_sympatry():
    d = z.DEC(areas=3, dispersal=0.1, extinction=0.1)
    a = d.encode({0})
    rng = np.random.default_rng(0)
    assert d.cladogenesis(a, rng) == (a, a)


def test_cladogenesis_widespread_gives_a_single_area_daughter():
    d = z.DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1)
    R = d.encode({"A", "B", "C"})
    rng = np.random.default_rng(0)
    for _ in range(50):
        i1, i2 = d.cladogenesis(R, rng)
        r1, r2 = set(d.states[i1]), set(d.states[i2])
        assert r1 and r2                                        # both non-empty
        assert r1 <= {"A", "B", "C"} and r2 <= {"A", "B", "C"}  # both within the ancestor
        assert len(r1) == 1 or len(r2) == 1                    # one daughter is a single area


# --------------------------------------------------------------------------- simulation
def test_simulate_biogeography_basic():
    tree = _tree()
    d = z.DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1)
    res = z.simulate_biogeography(tree, d, root_state={"A"}, seed=1)
    assert res.kind == "discrete"
    assert res.label(res.node_values[tree.root]) == ("A",)      # root range respected
    for leaf in tree.extant_leaves():
        rng_labels = res.labeled_values()[leaf]
        assert 1 <= len(rng_labels) <= 3 and set(rng_labels) <= {"A", "B", "C"}
    assert set(res.ancestral_states()) == set(tree.internal_nodes())


def test_root_state_as_index_or_range():
    tree = _tree()
    d = z.DEC(areas=3, dispersal=0.1, extinction=0.1)
    by_range = z.simulate_biogeography(tree, d, root_state={1}, seed=1)
    by_index = z.simulate_biogeography(tree, d, root_state=d.encode({1}), seed=1)
    assert by_range.label(by_range.node_values[tree.root]) == (1,)
    assert {k.name: v for k, v in by_index.labeled_values().items()} == \
           {k.name: v for k, v in by_range.labeled_values().items()}


def test_dispersal_grows_ranges():
    tree = _tree()

    def mean_tip_range(disp, ext, reps=40):
        m = z.DEC(areas=4, dispersal=disp, extinction=ext)
        rng = np.random.default_rng(0)
        sizes = []
        for _ in range(reps):
            r = z.simulate_biogeography(tree, m, root_state={0}, rng=rng)
            sizes += [len(r.label(r.node_values[leaf])) for leaf in tree.extant_leaves()]
        return np.mean(sizes)

    assert mean_tip_range(0.4, 0.1) > mean_tip_range(0.02, 0.4) + 0.3


def test_biogeography_reproducible():
    tree = _tree()
    d = z.DEC(areas=3, dispersal=0.2, extinction=0.15)
    a = z.simulate_biogeography(tree, d, root_state={0}, seed=7).labeled_values()
    b = z.simulate_biogeography(tree, d, root_state={0}, seed=7).labeled_values()
    assert {k.name: v for k, v in a.items()} == {k.name: v for k, v in b.items()}


def test_biogeography_validation():
    with pytest.raises(ValueError):
        z.DEC(areas=2, dispersal=[[0, 1], [1, 0], [0, 0]], extinction=0.1)  # dispersal shape
    with pytest.raises(ValueError):
        z.DEC(areas=2, dispersal=-0.1, extinction=0.1)                      # negative rate
    with pytest.raises(ValueError):
        z.DEC(areas=3, dispersal=0.1, extinction=0.1, max_range_size=5)     # out of range
    with pytest.raises(ValueError):
        z.DEC(areas=3, dispersal=0.1, extinction=0.1).encode([])           # empty range


def test_dec_requires_binary_tree():
    tree = z.read_newick("((A:1,B:1,C:1):1,D:2);")                 # a hard polytomy (3-child node)
    with pytest.raises(ValueError):
        z.simulate_biogeography(tree, z.DEC(areas=3, dispersal=0.1, extinction=0.1), seed=1)


def test_dec_needs_two_areas():
    with pytest.raises(ValueError):
        z.DEC(areas=1, dispersal=0.1, extinction=0.1)


def test_dec_invalid_root_range_is_clean_error():
    with pytest.raises(ValueError):                                # not a cryptic KeyError
        z.DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1,
              max_range_size=2, root={"A", "B", "C"})
