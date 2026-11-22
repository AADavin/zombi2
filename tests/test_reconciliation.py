"""Tests for :mod:`zombi2.tools.reconciliation` — the history projected onto what a dataset holds.

The complete reconciliation is the truth and is not what a method is scored against, because a method
sees only survivors. These check the projection: that it keeps exactly the observable events, that
what it does with the unobservable ones is defensible, and — the part most likely to go wrong — that
the result is still a *valid* reconciliation once nodes have been moved.
"""

import collections
import xml.etree.ElementTree as ET

import pytest

from zombi2.genomes import simulate_genomes_family
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.species import simulate_species_tree
from zombi2.tools.reconciliation import (ORIGINATION, TRANSFER, extant_reconciliation,
                                         origins_tsv, visible_branches)
from zombi2.tree import prune

#: regimes that differ in what is *invisible*: none, dead lineages, and unobserved survivors. The
#: three with no extinction differ in the transfers: none, transfers between branches, and transfers
#: that may land on the branch they left and may replace the copy they meet
_REGIMES = {
    "no extinction": (dict(birth=1.0, death=0.0, n_extant=25, seed=3),
                      dict(duplication=0.2, loss=0.3, origination=0.4, initial_families=30, seed=5)),
    "transfers, no extinction": (dict(birth=1.0, death=0.0, n_extant=20, seed=1),
                                 dict(duplication=0.1, transfer=0.5, loss=0.5, initial_families=40,
                                      seed=2)),
    "self-transfers, no extinction": (dict(birth=1.0, death=0.0, n_extant=20, seed=4),
                                      dict(duplication=0.1, transfer=0.6, loss=0.4,
                                           initial_families=40, self_transfer=True, replacement=0.5,
                                           seed=7)),
    "heavy extinction": (dict(birth=1.2, death=0.6, n_extant=40, seed=11),
                         dict(duplication=0.15, transfer=0.3, loss=0.2, origination=0.3,
                              initial_families=30, max_family_size=None, seed=11)),
    "unsampled survivors": (dict(birth=1.2, death=0.6, n_extant=30, sampling=0.6, seed=5),
                            dict(duplication=0.2, transfer=0.4, loss=0.3, initial_families=25,
                                 max_family_size=None, seed=5)),
}
#: the regimes where every branch is one the extant tree names
_ALL_VISIBLE = [name for name in _REGIMES if name.endswith("no extinction")]


def _run(name):
    sp_kw, g_kw = _REGIMES[name]
    sp = simulate_species_tree(**sp_kw)
    return sp, simulate_genomes_family(sp, **g_kw)


def _walk(root):
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return out


def _extant_leaves(root):
    return sorted((n.species, n.copy) for n in _walk(root) if n.kind == "extant")


def _survives(node):
    return any(n.kind == "extant" for n in _walk(node))


def _two_tips():
    """The smallest species tree, for histories written by hand: ``(tree, root, a, b, at)``, where
    ``at(x)`` is a time ``x`` of the way down the two tip branches."""
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=2, seed=1).complete_tree
    a, b = tree.nodes[tree.root].children
    split, now = tree.nodes[tree.root].end_time, tree.nodes[a].end_time
    return tree, tree.root, a, b, lambda x: split + x * (now - split)


@pytest.mark.parametrize("regime", list(_REGIMES))
@pytest.mark.parametrize("scope", ["true", "recoverable"])
def test_the_projection_is_a_valid_reconciliation(regime, scope):
    # the part most likely to be wrong: nodes MOVE when their own branch is not one the extant tree
    # names, and a moved node can easily end up somewhere its parent cannot reach. Checked here
    # against the trees themselves rather than against the builder's own bookkeeping.
    sp, g = _run(regime)
    tree = sp.complete_tree
    et = prune(tree, keep="extant")
    vis = visible_branches(tree)
    ancestors, tips = {}, {i for i, n in et.nodes.items() if not n.children}
    for i in sorted(et.nodes):
        p = et.nodes[i].parent
        ancestors[i] = (ancestors[p] | {p}) if p is not None else set()

    checked = 0
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, tree, scope=scope, visible=vis)
        if rec is None:
            continue
        for n in _walk(rec.root):
            checked += 1
            assert n.species in et.nodes, "every branch named must be one the extant tree has"
            if n.is_leaf:
                assert n.kind in ("extant", "loss")
                if n.kind == "extant":
                    assert n.species in tips
            elif n.kind in ("speciation", "duplication"):
                # vertical descent only: a child sits on its parent's branch or below it. (A transfer
                # is the exception by definition — that is what makes it a transfer.)
                for c in n.children:
                    assert c.species == n.species or n.species in ancestors[c.species]
    assert checked > 100


@pytest.mark.parametrize("regime", list(_REGIMES))
def test_the_reconciliation_is_of_the_gene_tree_a_dataset_holds(regime):
    # whatever the projection does to internal nodes, the observable genes must come through
    # untouched — otherwise it is a reconciliation of some other tree
    sp, g = _run(regime)
    vis = visible_branches(sp.complete_tree)
    seen = 0
    for fam, gt in g.gene_trees.items():
        if gt.extant is None:
            continue
        seen += 1
        for scope in ("true", "recoverable"):
            rec = extant_reconciliation(gt, sp.complete_tree, scope=scope, visible=vis)
            assert _extant_leaves(rec.root) == _extant_leaves(gt.extant)
    assert seen > 5


def test_true_contains_recoverable():
    # the whole reason both are written: 'true' carries the family's ancestral presence and the
    # losses that narrowed it, 'recoverable' is that trimmed back to the surviving copies' ancestor.
    # Trimming is possible, inventing is not, so the richer one is the one to ship.
    sp, g = _run("heavy extinction")
    vis = visible_branches(sp.complete_tree)
    extra, pairs = 0, 0
    for fam, gt in g.gene_trees.items():
        t = extant_reconciliation(gt, sp.complete_tree, scope="true", visible=vis)
        r = extant_reconciliation(gt, sp.complete_tree, scope="recoverable", visible=vis)
        if t is None:
            continue
        pairs += 1
        assert t.losses >= r.losses                       # only ever adds
        assert len(_walk(t.root)) >= len(_walk(r.root))
        extra += t.losses - r.losses
    assert pairs > 5 and extra > 0, "some family must be ancestrally present above its survivors"


def test_a_family_that_arrived_from_a_dead_lineage_says_so():
    # the marker exists because the arrival is INVISIBLE in the reconciliation itself: with the donor
    # gone the copy simply appears, which a reader would call an origination and be wrong about
    sp, g = _run("heavy extinction")
    vis = visible_branches(sp.complete_tree)
    kinds = collections.Counter()
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, sp.complete_tree, visible=vis)
        if rec is not None:
            kinds[rec.entered_by] += 1
    assert kinds[ORIGINATION] > 0 and kinds[TRANSFER] > 0
    assert set(kinds) == {ORIGINATION, TRANSFER}


@pytest.mark.parametrize("regime", _ALL_VISIBLE)
def test_with_no_extinction_there_is_nothing_to_move(regime):
    # the control. Every lineage survives, so every branch is one the extant tree names and each is
    # its own image: no node can be relocated and no transfer synthesised, which is what isolates
    # those behaviours to the regimes that actually have invisible branches. That holds with
    # transfers too — a donor's copy dying is not a donor's lineage dying (issue #490).
    #
    # The two scopes still differ here, and that is the point worth pinning: 'true' vs 'recoverable'
    # is about ancestral presence above the surviving copies, NOT about extinction. A family present
    # at the root and lost outside one clade has history no method can reach even when every lineage
    # is sampled.
    sp, g = _run(regime)
    tree = sp.complete_tree
    kept, image = visible_branches(tree)
    assert kept == set(tree.nodes)
    assert image == {i: i for i in tree.nodes}
    events = 0
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, tree, visible=(kept, image))
        if rec is None:
            continue
        assert rec.entered_by == ORIGINATION          # no dead donor to arrive from
        real = {n.copy: n for n in _walk(gt.complete) if not n.is_leaf}
        for n in _walk(rec.root):
            if not n.is_leaf:                         # every event kept is the event that happened
                events += 1
                assert (n.kind, n.species, n.time) == (real[n.copy].kind, real[n.copy].species,
                                                       real[n.copy].time)
    assert events > 50


@pytest.mark.parametrize("regime", _ALL_VISIBLE[1:])
def test_a_donor_copy_that_died_on_a_visible_branch_leaves_a_transfer_and_a_loss(regime):
    # issue #490. The donor's branch is in the extant tree and only the donor's own copy died, so
    # the transfer stays where it happened, with a loss beside the copy that left. The exception is
    # a copy that landed on the branch it left: one copy on one branch, as after a dead duplicate.
    sp, g = _run(regime)
    tree = sp.complete_tree
    vis = visible_branches(tree)
    kept, came_back = 0, 0
    for fam, gt in g.gene_trees.items():
        rec = extant_reconciliation(gt, tree, visible=vis)
        if rec is None:
            continue
        projected = {n.copy: n for n in _walk(rec.root) if not n.is_leaf}
        for n in _walk(gt.complete):
            if n.kind != "transfer":
                continue
            left = next((c for c in n.children if c.species != n.species), n.children[-1])
            stayed = next(c for c in n.children if c is not left)
            if _survives(stayed) or not _survives(left):
                continue
            if left.species == n.species:
                came_back += 1
                assert n.copy not in projected
                continue
            kept += 1
            node = projected[n.copy]
            assert (node.kind, node.species, node.time) == ("transfer", n.species, n.time)
            loss, arrived = node.children             # the donor's side first
            assert (loss.kind, loss.species, loss.copy) == ("loss", n.species, stayed.copy)
            assert arrived.species == left.species
    assert kept >= 10
    assert came_back >= 1 or not regime.startswith("self")


def test_a_dead_donor_copy_by_hand():
    # the gene follows both daughters. On b it sends a copy to a, then b's own copy is lost. Every
    # branch is visible, so the projection is the history itself: a transfer from b, and a loss on b
    tree, root, a, b, at = _two_tips()
    transfer = GeneNode("transfer", b, at(0.3), 2, [GeneNode("loss", b, at(0.6), 3),
                                                     GeneNode("extant", a, at(1.0), 4)])
    whole = GeneNode("speciation", root, at(0.0), 0, [GeneNode("extant", a, at(1.0), 1), transfer])
    for scope in ("true", "recoverable"):
        rec = extant_reconciliation(GeneTree(0, whole, 0.0), tree, scope=scope)
        assert rec.root == whole
        assert (rec.entered_by, rec.branch, rec.losses) == (ORIGINATION, root, 1)


def test_a_family_whose_only_survivor_is_the_copy_that_left():
    # 'true' roots the family on b, where it began, with the transfer and the loss; 'recoverable'
    # trims both, because the copy on a is all that is left to recover
    tree, root, a, b, at = _two_tips()
    whole = GeneNode("transfer", b, at(0.3), 0, [GeneNode("loss", b, at(0.6), 1),
                                                 GeneNode("extant", a, at(1.0), 2)])
    true = extant_reconciliation(GeneTree(0, whole, at(0.1)), tree, scope="true")
    assert true.root == whole
    assert (true.entered_by, true.branch, true.losses) == (ORIGINATION, b, 1)
    recoverable = extant_reconciliation(GeneTree(0, whole, at(0.1)), tree, scope="recoverable")
    assert recoverable.root == GeneNode("extant", a, at(1.0), 2)
    assert (recoverable.branch, recoverable.losses) == (a, 0)


def test_a_copy_that_lands_on_the_branch_it_left_leaves_nothing_to_see():
    tree, root, a, b, at = _two_tips()
    whole = GeneNode("transfer", b, at(0.3), 0, [GeneNode("loss", b, at(0.6), 1),
                                                 GeneNode("extant", b, at(1.0), 2)])
    rec = extant_reconciliation(GeneTree(0, whole, at(0.1)), tree)
    assert rec.root == GeneNode("extant", b, at(1.0), 2)
    assert (rec.entered_by, rec.branch, rec.losses) == (ORIGINATION, b, 0)


def test_scope_is_validated():
    sp, g = _run("no extinction")
    gt = next(iter(g.gene_trees.values()))
    with pytest.raises(ValueError, match="'true' or 'recoverable'"):
        extant_reconciliation(gt, sp.complete_tree, scope="best")


def test_origins_tsv_names_every_family_once():
    sp, g = _run("heavy extinction")
    vis = visible_branches(sp.complete_tree)
    recs = {f: r for f, gt in g.gene_trees.items()
            if (r := extant_reconciliation(gt, sp.complete_tree, visible=vis)) is not None}
    rows = origins_tsv(recs).splitlines()
    assert rows[0] == "family\tentered_by\tbranch\tlosses"
    assert len(rows) == len(recs) + 1
    assert all(r.split("\t")[1] in (ORIGINATION, TRANSFER) for r in rows[1:])
    assert all(r.split("\t")[2].startswith("n") for r in rows[1:])   # a branch of the EXTANT tree


def test_the_cli_writes_both_scopes_and_the_origins_table(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1.2", "--death", "0.6", "--n-extant", "20",
                 "--seed", "11", "--quiet"]) == 0
    assert main(["genomes", str(run), "--duplication", "0.15", "--transfer", "0.3", "--loss", "0.2",
                 "--origination", "0.3", "--initial-families", "10", "--seed", "11", "--quiet"]) == 0
    assert main(["tools", "format", str(run), "--format", "recphylo", "--recphylo", "extant",
                 "--quiet"]) == 0
    out = run / "genomes" / "recphylo"
    assert (out / "family_origins.tsv").exists()
    assert not list(out.glob("recphylo_fam*[0-9].xml")), "extant scope writes no complete files"
    pairs = sorted(p.name for p in out.glob("*_true.xml"))
    assert pairs and len(pairs) == len(list(out.glob("*_recoverable.xml")))

    for p in out.glob("*_true.xml"):
        root = ET.fromstring(p.read_text(encoding="utf-8"))       # it is well-formed XML...
        names = {c.text for c in root.find("spTree").iter("name")}
        assert names and not any(n.startswith("e") for n in names)   # ...against the EXTANT tree
