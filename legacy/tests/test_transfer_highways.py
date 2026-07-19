"""Transfer highways — donor→recipient bias at the recipient seam (TransferModel.pair).

A PairModifier multiplies each candidate recipient's selection weight by a factor for the
(donor, recipient) pair, specified as explicit branch pairs and/or clade→clade blocks. It runs on
the Python engine.
"""

import numpy as np
import pytest

from zombi2 import (
    BirthDeath,
    PairModifier,
    Rates,
    TransferModel,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes.events import EventType
from zombi2.tree import Tree, TreeNode


def _mini_tree() -> Tree:
    """((A,B)ab,(C,D)cd)root with all internal/leaf nodes named."""
    root = TreeNode(name="root", time=0.0)
    ab, cd = TreeNode(name="ab", time=1.0), TreeNode(name="cd", time=1.0)
    root.add_child(ab); root.add_child(cd)
    for nm in ("A", "B"):
        ab.add_child(TreeNode(name=nm, time=2.0))
    for nm in ("C", "D"):
        cd.add_child(TreeNode(name=nm, time=2.0))
    return Tree(root, total_age=2.0)


def _subtree_names(node) -> set:
    names, stack = set(), [node]
    while stack:
        n = stack.pop()
        names.add(n.name)
        stack.extend(n.children)
    return names


# --- unit: PairModifier.factor (explicit pairs, clade blocks, MRCA, direction) -------

def test_pair_factor_explicit_and_block():
    tree = _mini_tree()
    by = {n.name: n for n in tree.nodes_preorder()}

    # explicit A->C (x3) AND a clade block ab->cd (x5), specified by tip sets (MRCA = ab, cd)
    pm = PairModifier(pairs={("A", "C"): 3.0}, blocks=[({"A", "B"}, {"C", "D"}, 5.0)])
    pm.bind(tree)
    assert pm.factor(by["A"], by["C"]) == 15.0   # explicit x block = 3 * 5
    assert pm.factor(by["A"], by["D"]) == 5.0    # block only
    assert pm.factor(by["C"], by["A"]) == 1.0    # reverse direction: nothing applies

    # a block by node name resolves to that node's whole subtree
    pm2 = PairModifier(blocks=[("ab", "cd", 4.0)])
    pm2.bind(tree)
    assert pm2.factor(by["B"], by["C"]) == 4.0
    assert pm2.factor(by["A"], by["A"]) == 1.0   # A (recipient) is not under cd


def test_pair_modifier_validation():
    tree = _mini_tree()
    with pytest.raises(ValueError):
        PairModifier()                                   # no pairs and no blocks
    with pytest.raises(ValueError):
        PairModifier(pairs={("A", "B"): -1.0})           # negative pair factor
    with pytest.raises(ValueError):
        PairModifier(blocks=[("A", "B", -2.0)])          # negative block factor
    with pytest.raises(ValueError):
        PairModifier(blocks=[("nope", "cd", 2.0)]).bind(tree)      # unknown node name
    with pytest.raises(ValueError):
        PairModifier(blocks=[({"A", "Z"}, "cd", 2.0)]).bind(tree)  # unknown tip in a clade spec


# --- simulation: a highway biases recipient choice ----------------------------------

def _sim_tree():
    return simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)


def _transfers_into(tree, target_names, pair, seeds=range(6)):
    """(#transfers whose recipient is in target_names, #transfers total) over several seeds."""
    tm = TransferModel(pair=pair) if pair is not None else TransferModel()
    into = total = 0
    for s in seeds:
        g = simulate_genomes(tree, Rates(transfer=0.4, loss=0.4, origination=0.2),
                             transfers=tm, initial_families=6, seed=s, max_family_size=25)
        for r in g.event_log:
            if r.event is EventType.TRANSFER:
                total += 1
                into += r.recipient in target_names
    return into, total


def _a_clade(tree):
    root = tree.root
    node = next(n for n in tree.nodes_preorder() if n is not root and n.children)
    return node, _subtree_names(node)


def test_zero_block_forbids_transfers_into_a_clade():
    tree = _sim_tree()
    clade, names = _a_clade(tree)
    # a highway from everything (root) into the clade with factor 0 -> the clade never receives
    forbid = PairModifier(blocks=[(tree.root.name, clade.name, 0.0)])
    into, total = _transfers_into(tree, names, forbid)
    assert total > 0            # transfers did happen (elsewhere)
    assert into == 0            # but none landed inside the forbidden clade


def test_block_highway_concentrates_transfers_into_a_clade():
    tree = _sim_tree()
    clade, names = _a_clade(tree)
    boost = PairModifier(blocks=[(tree.root.name, clade.name, 50.0)])
    into_hw, tot_hw = _transfers_into(tree, names, boost)
    into_no, tot_no = _transfers_into(tree, names, None)
    frac_hw = into_hw / tot_hw
    frac_no = into_no / tot_no
    assert frac_hw > frac_no + 0.1   # the highway clearly concentrates transfers into the clade


def test_highway_reproducible_given_seed():
    tree = _sim_tree()
    clade, _ = _a_clade(tree)
    make = lambda: TransferModel(pair=PairModifier(blocks=[(tree.root.name, clade.name, 8.0)]))
    a = simulate_genomes(tree, Rates(transfer=0.4, loss=0.3, origination=0.2),
                         transfers=make(), initial_families=6, seed=4, max_family_size=25)
    b = simulate_genomes(tree, Rates(transfer=0.4, loss=0.3, origination=0.2),
                         transfers=make(), initial_families=6, seed=4, max_family_size=25)
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)
