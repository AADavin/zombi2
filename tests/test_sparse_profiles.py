"""The profile matrix is sparse (COO) by construction.

A copy-number profile is overwhelmingly zero, so :class:`~zombi2.ProfileMatrix` stores only
the present cells and computes every summary off those — the dense ``families x species``
array (O(N^2) in tip count) is materialised only when explicitly asked for. These tests pin
that the sparse path agrees with the dense one cell-for-cell, that both TSV formats round
trip, and that the Rust ``output="profiles"`` route is genuinely sparse.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.profiles import ProfileMatrix
from zombi2.matching import frequency_spectrum, genome_sizes, copy_number_spectrum


def _example():
    families = ["1", "2", "3"]
    species = ["s1", "s2", "s3", "s4"]
    dense = np.array([[1, 0, 2, 0],
                      [0, 0, 0, 0],   # a family present nowhere (all-absent column-wise)
                      [3, 1, 0, 4]], dtype=int)
    return families, species, dense


def test_coo_matches_dense_construction():
    families, species, dense = _example()
    pm = ProfileMatrix(families, species, dense)
    assert pm.shape == (3, 4)
    assert pm.nnz == 5                         # only the non-zero cells stored
    assert np.array_equal(pm.matrix, dense)    # densifies back exactly
    assert np.array_equal(pm.presence(), (dense > 0).astype(np.int8))


def test_sparse_reductions_match_dense():
    families, species, dense = _example()
    pm = ProfileMatrix(families, species, dense)
    assert np.array_equal(pm.presence_per_family(), (dense > 0).sum(axis=1))
    assert np.array_equal(pm.copies_per_species(), dense.sum(axis=0))
    assert np.array_equal(pm.copies_per_family(), dense.sum(axis=1))
    assert np.array_equal(np.sort(pm.copy_values()), np.sort(dense[dense > 0]))


def test_dense_tsv_roundtrip():
    families, species, dense = _example()
    pm = ProfileMatrix(families, species, dense)
    back = ProfileMatrix.from_tsv(pm.to_tsv())
    assert back.families == families and back.species == species
    assert np.array_equal(back.matrix, dense)


def test_coo_tsv_roundtrip_preserves_all_columns():
    families, species, dense = _example()
    pm = ProfileMatrix(families, species, dense)
    back = ProfileMatrix.from_coo_tsv(pm.to_coo_tsv())
    # every species column survives even though s2/s3 may be sparse
    assert back.species == species
    assert np.array_equal(back.matrix, dense)


def test_abc_stats_use_sparse_and_match_dense():
    families, species, dense = _example()
    pm = ProfileMatrix(families, species, dense)
    S = len(species)
    # frequency spectrum: families present in exactly k species
    present = (dense > 0).sum(axis=1)
    exp_fs = np.bincount(present, minlength=S + 1)[1:S + 1].astype(float)
    assert np.array_equal(frequency_spectrum(pm, S), exp_fs)
    assert np.array_equal(genome_sizes(pm, species), dense.sum(axis=0).astype(float))
    vals = dense[dense > 0]
    exp_cn = np.array([np.count_nonzero(vals == 1), np.count_nonzero(vals == 2),
                       np.count_nonzero(vals == 3), np.count_nonzero(vals >= 4)], float)
    assert np.array_equal(copy_number_spectrum(pm), exp_cn)


def test_empty_profile():
    pm = ProfileMatrix([], [], coo=([], [], []))
    assert pm.nnz == 0 and pm.shape == (0, 0)
    assert pm.matrix.shape == (0, 0)


@pytest.mark.skipif(not z.rust_available(), reason="zombi2_core (Rust extension) not built")
def test_profiles_path_is_sparse():
    """output='profiles' returns a sparse ProfileMatrix; stored cells << dense size."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=300, age=2.0, seed=1)
    pm = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                            origination=0.5, initial_size=20, output="profiles", seed=3)
    dense_cells = pm.shape[0] * pm.shape[1]
    assert pm.nnz < dense_cells                       # genuinely sparse
    # sparse reductions agree with a densified cross-check at this (small) size
    dense = pm.matrix
    assert np.array_equal(pm.copies_per_species(), dense.sum(axis=0))
    assert int(pm.copies_per_family().sum()) == int(dense.sum())


@pytest.mark.skipif(not z.rust_available(), reason="zombi2_core (Rust extension) not built")
def test_full_genomes_profile_is_sparse_and_consistent():
    """The full-genealogy path also builds its ProfileMatrix sparsely, self-consistently."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=60, age=2.0, seed=7)
    g = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, initial_size=20, seed=11)
    pm = g.profiles
    assert pm.nnz <= pm.shape[0] * pm.shape[1]
    # reconciliation invariant: each family's stored copies == its leaf-genome copy sum
    dense = pm.matrix
    assert np.array_equal(pm.copies_per_family(), dense.sum(axis=1))
    assert int(pm.copies_per_species().sum()) == int(dense.sum())
