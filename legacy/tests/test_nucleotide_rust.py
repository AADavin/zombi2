"""Rust engine for the nucleotide model, reached via ``simulate_nucleotide_genomes(...,
output="profiles")`` (leaf segments -> blocks / profile / mosaic / trace-back).

Skipped if zombi2_core isn't built. Validates internal consistency (the same invariants the
pure-Python tests assert) on the Rust result; cross-engine equivalence is statistical (the
RNG differs) and checked outside the suite.
"""

import pytest

import zombi2 as z

pytestmark = pytest.mark.skipif(not z.rust_available(),
                                reason="zombi2_core (Rust extension) not built")

FULL = dict(inversion=0.02, duplication=0.01, loss=0.015, transfer=0.008,
            transposition=0.01, origination=0.3, root_length=500, extension=0.9)


def _sim(tree, **kw):
    return z.simulate_nucleotide_genomes(tree, output="profiles", **kw)


def _tree(n=10, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=n, age=3.0, seed=seed)


def test_reproducible_same_seed():
    tree = _tree()
    a = _sim(tree, seed=7, **FULL)
    b = _sim(tree, seed=7, **FULL)
    assert [a2.block_id for a2 in a.blocks] == [b2.block_id for b2 in b.blocks]
    for leaf, ga in a.leaf_genomes.items():
        assert ga.to_cells() == b.leaf_genomes[leaf].to_cells()


def test_species_are_extant_leaves():
    tree = _tree(n=12)
    res = _sim(tree, seed=3, **FULL)
    _ids, species, M = res.profile_matrix()
    assert M.shape[1] == len(tree.extant_leaves())
    assert set(species) == {n.name for n in tree.extant_leaves()}


def test_mosaic_reassembles_each_leaf():
    tree = _tree(n=8)
    res = _sim(tree, seed=2, **FULL)
    amap = {a.block_id: a for a in res.blocks}
    for leaf, g in res.leaf_genomes.items():
        cells = []
        for aid, strand in res.leaf_mosaic(leaf):
            a = amap[aid]
            if strand == 1:
                cells.extend((a.source, p, 1) for p in range(a.start, a.end))
            else:
                cells.extend((a.source, p, -1) for p in range(a.end - 1, a.start - 1, -1))
        assert cells == g.to_cells()


def test_blocks_cover_exactly_surviving_positions():
    tree = _tree(n=8)
    res = _sim(tree, seed=4, **FULL)
    block_positions = {(a.source, p) for a in res.blocks for p in range(a.start, a.end)}
    surviving = set()
    for g in res.leaf_genomes.values():
        surviving.update((src, p) for src, p, _st in g.to_cells())
    assert block_positions == surviving


def test_profile_matches_leaf_coverage_and_grows():
    tree = _tree(n=10)
    res = _sim(tree, seed=5, **FULL)
    _ids, _sp, M = res.profile_matrix()
    assert M.max() >= 2                         # duplication/transfer -> paralogs
    assert M.min() == 0                         # loss / patchy novel sources
    assert len({a.source for a in res.blocks}) > 1  # origination made novel sources


def test_gene_trees_unavailable_on_profiles_path():
    tree = _tree(n=6)
    res = _sim(tree, seed=1, **FULL)
    with pytest.raises(NotImplementedError):
        res.block_gene_trees()
    with pytest.raises(NotImplementedError):
        res.block_histories()


def test_profiles_requires_numeric_extension():
    tree = _tree(n=6)
    with pytest.raises(ValueError):
        z.simulate_nucleotide_genomes(tree, output="profiles", extension=None, **{
            k: v for k, v in FULL.items() if k != "extension"})


def test_inversion_only_conserves_content():
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.0), n_tips=6, age=2.0, seed=9)
    res = _sim(tree, inversion=0.05, root_length=200, extension=0.9, seed=9)
    root = {("1", i) for i in range(200)}
    for g in res.leaf_genomes.values():
        origins = {(src, p) for src, p, _st in g.to_cells()}
        assert origins == root                  # inversion loses nothing
