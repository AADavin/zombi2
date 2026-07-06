"""Gene-tree reconstruction and growth-regulation tests."""


from zombi2 import (
    BirthDeath,
    SharedRates,
    Yule,
    simulate_genomes,
    simulate_species_tree,
)


def _n_leaves(newick: str) -> int:
    # In any Newick tree the number of leaves equals the number of commas + 1.
    return newick.count(",") + 1


def test_extant_tree_leaf_count_matches_profile():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=2)
    g = simulate_genomes(tree, duplication=0.15, transfer=0.1, loss=0.2,
                         origination=0.5, initial_families=10, seed=2)
    trees = g.gene_trees()
    fam_row = {f: i for i, f in enumerate(g.profiles.families)}
    checked = 0
    for family, (complete, extant) in trees.items():
        expected = int(g.profiles.matrix[fam_row[family]].sum()) if family in fam_row else 0
        if extant is None:
            assert expected == 0
        else:
            assert _n_leaves(extant) == expected, (family, _n_leaves(extant), expected)
            checked += 1
    assert checked > 0  # we actually reconstructed some extant trees


def test_complete_tree_has_more_or_equal_leaves_than_extant():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=12, age=4.0, seed=3)
    g = simulate_genomes(tree, duplication=0.2, transfer=0.05, loss=0.3,
                         origination=0.5, initial_families=12, seed=3)
    for complete, extant in g.gene_trees().values():
        assert complete  # always produced
        if extant is not None:
            assert _n_leaves(complete) >= _n_leaves(extant)  # complete keeps losses too


def test_carrying_capacity_bounds_family_size():
    tree = simulate_species_tree(Yule(1.0), n_tips=8, age=3.0, seed=1)  # no extinction
    g = simulate_genomes(
        tree,
        SharedRates(duplication=2.0, transfer=0.0, loss=0.1, origination=0.0, carrying_capacity=5),
        initial_families=3, seed=1,
    )
    # duplication stops at n=K and there are no transfers, so no family exceeds K
    assert g.profiles.matrix.max() <= 5


def test_max_family_size_caps_family_size():
    tree = simulate_species_tree(Yule(1.0), n_tips=8, age=3.0, seed=1)
    g = simulate_genomes(
        tree, SharedRates(duplication=3.0, transfer=0.0, loss=0.05, origination=0.0),
        initial_families=3, max_family_size=4, seed=1,
    )
    assert g.profiles.matrix.max() <= 4


def test_write_outputs(tmp_path):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=5)
    g = simulate_genomes(tree, duplication=0.15, transfer=0.1, loss=0.2,
                         origination=0.5, initial_families=8, seed=5)
    g.write(tmp_path)
    assert (tmp_path / "species_tree.nwk").exists()
    assert (tmp_path / "Profiles.tsv").exists()
    assert (tmp_path / "Transfers.tsv").exists()
    assert (tmp_path / "Gene_family_summary.tsv").exists()
    assert any((tmp_path / "gene_trees").glob("*_complete.nwk"))
