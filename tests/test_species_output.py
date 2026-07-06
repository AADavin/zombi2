"""Species-output polish: distinct extinct (e*) vs unsampled (u*) leaf names, the reconstructed
`species_tree_extant.nwk` companion file, and the banner showing only on --help."""

import pytest

from zombi2.cli import main
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.tree import read_newick

BANNER = "a simulator of species trees"


def _leaf_kinds(tree):
    n = [x for x in tree.leaves() if x.name.startswith("n")]
    u = [x for x in tree.leaves() if x.name.startswith("u")]
    e = [x for x in tree.leaves() if x.name.startswith("e")]
    return n, u, e


def test_forward_splits_extinct_and_unsampled_names():
    """A forward run with extinction (death>0) and incomplete sampling (ρ<1) names the two dead
    fates differently: e* for lineages gone before the present, u* for unsampled-extant ghosts."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.5, sampling_fraction=0.5),
                                 age=6.0, direction="forward", seed=1)
    n, u, e = _leaf_kinds(tree)
    assert n and u and e                                   # all three fates present in this scenario
    assert len(n) + len(u) + len(e) == len(tree.leaves())  # every leaf named n/u/e, no overlap
    present = tree.total_age
    for leaf in u:                                         # u* = alive at the present, not sampled
        assert not leaf.is_extant and not leaf.sampled
        assert abs(leaf.time - present) <= 1e-9 * max(1.0, present)
    for leaf in e:                                         # e* = gone before the present
        assert not leaf.is_extant
        assert leaf.time < present - 1e-9
    for leaf in n:                                         # n* = sampled extant
        assert leaf.is_extant


def test_cli_writes_extant_tree_when_dead_tips(tmp_path, capsys):
    """Forward mode with dead tips also writes the pruned reconstructed tree, holding exactly the
    sampled-extant (n*) tips — no e*/u*."""
    out = tmp_path / "fwd"
    rc = main(["species", "--mode", "forward", "--age", "6", "--birth", "1", "--death", "0.5",
               "--sampling-fraction", "0.5", "--seed", "1", "-o", str(out)])
    assert rc == 0
    complete = read_newick((out / "species_tree.nwk").read_text())
    assert (out / "species_tree_extant.nwk").exists()
    extant = read_newick((out / "species_tree_extant.nwk").read_text())
    n_complete, _, _ = _leaf_kinds(complete)
    assert len(extant.leaves()) == len(n_complete)         # pruned == the n* tips
    assert all(leaf.is_extant for leaf in extant.leaves())
    assert not any(leaf.name.startswith(("e", "u")) for leaf in extant.leaves())
    assert "species_tree_extant.nwk" in capsys.readouterr().out


def test_cli_no_extant_file_when_all_tips_extant(tmp_path):
    """A backward run (complete sampling, no ghosts) has no dead tips, so no companion file."""
    out = tmp_path / "bwd"
    rc = main(["species", "--mode", "backward", "--tips", "12", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()
    assert not (out / "species_tree_extant.nwk").exists()


def test_banner_only_on_help(tmp_path, capsys):
    """The banner prints on --help, never on a normal run."""
    out = tmp_path / "s"
    rc = main(["species", "--mode", "backward", "--tips", "6", "--seed", "1", "-o", str(out)])
    assert rc == 0
    captured = capsys.readouterr()
    assert BANNER not in captured.out and BANNER not in captured.err   # no banner on a real run

    with pytest.raises(SystemExit):
        main(["--help"])
    assert BANNER in capsys.readouterr().out                            # but it is on --help


def test_cli_writes_species_nodes_table(tmp_path):
    """`species` writes species_nodes.tsv — one row per node (name, time, is_leaf, is_extant) —
    matching the simulated tree, on both forward and backward runs."""
    out = tmp_path / "fwd"
    rc = main(["species", "--mode", "forward", "--age", "6", "--birth", "1", "--death", "0.5",
               "--sampling-fraction", "0.5", "--seed", "1", "-o", str(out)])
    assert rc == 0
    # the same seed reproduces the exact tree the CLI built the table from
    tree = simulate_species_tree(BirthDeath(1.0, 0.5, sampling_fraction=0.5),
                                 age=6.0, direction="forward", seed=1)

    rows = (out / "species_nodes.tsv").read_text().splitlines()
    assert rows[0] == "name\ttime\tis_leaf\tis_extant"
    table = {}
    for line in rows[1:]:
        name, t, is_leaf, is_extant = line.split("\t")
        table[name] = (float(t), is_leaf == "True", is_extant == "True")
    assert len(table) == len(list(tree.nodes()))               # one row per node, names unique
    for node in tree.nodes():
        t, is_leaf, is_extant = table[node.name]
        assert is_leaf == (not node.children)                  # leaf iff no children
        assert is_extant == bool(node.is_extant)               # extant flag faithful to the sim
        assert abs(t - node.time) <= 1e-6 * max(1.0, abs(node.time))

    # written on a plain backward run too (no dead tips, no companion tree)
    out2 = tmp_path / "bwd"
    assert main(["species", "--mode", "backward", "--tips", "10", "--seed", "1", "-o", str(out2)]) == 0
    assert (out2 / "species_nodes.tsv").exists()
