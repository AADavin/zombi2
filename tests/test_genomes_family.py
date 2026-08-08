"""Tests for the D/L/O gene-family core (zombi2.genomes.family)."""

import collections
import math
import pathlib
import tempfile

import pytest

from zombi2.rates.scope import Global, PerCopy, PerLineage

from zombi2.rates import ScaledBy, modifiers as mod
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.tree import Node, Tree


def _tree(seed=1, n_extant=12, death=0.3):
    return simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)


def _transfers(events):
    """Reconstruct each transfer as ``(donor_lineage, recipient_lineage, time)`` from its two rows
    (same time + parent: the continuation on the donor, and the copy on the recipient — the one that
    names a recipient)."""
    rows = collections.defaultdict(list)
    for e in events:
        if e.kind == "transfer":
            rows[(e.time, e.parent)].append(e)
    out = []
    for pair in rows.values():
        xfer = next(r for r in pair if r.recipient is not None)     # the transferred copy, on the recipient
        cont = next(r for r in pair if r.recipient is None)         # the donor's continuation
        out.append((cont.lineage, xfer.lineage, xfer.time))
    return out


# --- the walk covers the whole complete tree -------------------------------

def test_genomes_on_every_node_including_extinct():
    sp = _tree(seed=3, death=0.5)
    g = simulate_genomes_family(sp, duplication=0.2, loss=0.2, origination=0.5, initial_families=4, seed=1)
    assert set(g.genomes) == set(sp.complete_tree.nodes)          # every node has a genome
    extinct = {n.id for n in sp.complete_tree.extinct_leaves()}
    assert extinct and extinct <= set(g.genomes)                 # extinct lineages included


def test_accepts_a_result_or_a_bare_tree():
    sp = _tree(seed=7)
    a = simulate_genomes_family(sp, origination=0.5, initial_families=3, seed=1)
    b = simulate_genomes_family(sp.complete_tree, origination=0.5, initial_families=3, seed=1)
    assert [(e.time, e.kind, e.copy) for e in a.edges] == [(e.time, e.kind, e.copy) for e in b.edges]


# --- determinism -----------------------------------------------------------

def test_deterministic_given_seed():
    sp = _tree(seed=2)
    kw = dict(duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=9)
    a, b = simulate_genomes_family(sp, **kw), simulate_genomes_family(sp, **kw)
    assert [(e.time, e.kind, e.lineage, e.copy, e.parent) for e in a.edges] == \
           [(e.time, e.kind, e.lineage, e.copy, e.parent) for e in b.edges]
    assert a.genomes == b.genomes


def test_different_seeds_differ():
    sp = _tree(seed=2)
    a = simulate_genomes_family(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=1)
    b = simulate_genomes_family(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=2)
    assert len(a.edges) != len(b.edges) or \
        [e.copy for e in a.edges] != [e.copy for e in b.edges]


# --- the three events behave --------------------------------------------------

def test_initial_families_seed_originations_at_the_crown():
    sp = _tree(seed=1)
    root = sp.complete_tree.root
    t0 = sp.complete_tree.nodes[root].birth_time
    g = simulate_genomes_family(sp, initial_families=6, seed=1)          # no D/L/O rates
    crown = [e for e in g.events if e.time == t0]
    assert len(crown) == 6 and all(e.kind == "origination" for e in crown)
    assert len(g.genomes[root]) == 6                               # the root carries all 6


def test_no_rates_means_pure_inheritance():
    # with every rate 0, nothing happens beyond the crown seeding + the speciation re-ids: every node
    # carries the root's families (per-node ids differ, but the family multiset is unchanged)
    sp = _tree(seed=4)
    g = simulate_genomes_family(sp, initial_families=3, seed=1)
    root_counts = g.family_counts(sp.complete_tree.root)
    assert all(g.family_counts(i) == root_counts for i in g.genomes)   # families inherited unchanged
    assert {e.kind for e in g.events} <= {"origination", "speciation"}  # only crown births + splits
    assert all(e.time == 0.0 for e in g.events if e.kind == "origination")


def test_origination_only_families_never_exceed_one_copy():
    sp = _tree(seed=6)
    g = simulate_genomes_family(sp, origination=0.8, initial_families=2, seed=1)  # no duplication
    for node_id in g.genomes:
        assert all(count == 1 for count in g.family_counts(node_id).values())


def test_duplication_grows_a_family():
    sp = _tree(seed=6)
    g = simulate_genomes_family(sp, duplication=0.8, initial_families=3, seed=1)  # no loss, no origination
    biggest = max((max(g.family_counts(i).values(), default=0)) for i in g.genomes)
    assert biggest > 1                                            # some family reached >1 copy
    # duplication never introduces a new family (only origination does); speciation re-ids at splits
    assert {e.kind for e in g.events} <= {"origination", "duplication", "speciation"}


def test_loss_can_shrink_and_empty_a_genome():
    # high loss, no origination/duplication except the seeded families -> some lineage loses all
    sp = _tree(seed=6)
    g = simulate_genomes_family(sp, loss=2.0, initial_families=4, seed=1)
    sizes = [len(g.genomes[i]) for i in g.genomes]
    assert min(sizes) < 4                                          # at least one node shrank
    assert any(e.kind == "loss" for e in g.events)


def test_summary_counts_the_genomes_that_emptied():
    # An emptied genome writes no row in profiles.tsv and leaves no gene tree, so nothing in the
    # outputs says it happened; `empty_genomes` is where a reader finds out.
    g = simulate_genomes_family(_tree(seed=6), loss=5.0, initial_families=10, seed=3)
    s = g.summary()
    assert s["empty_genomes"] == s["extant_genomes"] == 12          # every one of them
    assert s["families"]["surviving"] == 0 and s["genes_per_genome"]["max"] == 0
    # ... and the process itself is untouched: all ten founding families really were lost, no event
    # was refused. This is what tells "we made it visible" apart from "we added a floor".
    assert sum(1 for e in g.events if e.kind == "loss") == 10


def test_no_floor_means_a_last_copy_is_an_ordinary_copy():
    # A statistical guard against a floor being added here by accident. Loss is per copy and the
    # last copy is a copy like any other, so on a lone branch of length T a single gene survives
    # with probability exp(-λT) — exactly the survival of a pure death process, with no reflecting
    # boundary at one copy. A floor would push the empty fraction to zero; a partial one would pull
    # it below this.
    T, lam, n = 1.0, 0.5, 400
    tree = Tree({0: Node(0, None, 0.0, T, None, "extant")}, 0)
    empty = sum(1 for s in range(n)
                if not simulate_genomes_family(tree, loss=lam, initial_families=1, seed=s).genomes[0])
    p = 1 - math.exp(-lam * T)
    sd = math.sqrt(p * (1 - p) / n)                                # binomial, n independent runs
    assert abs(empty / n - p) < 4 * sd


def test_duplication_bifurcates_into_two_same_family_children():
    # ZOMBI1 model: a duplication ends the gene and starts two fresh ids descending from it
    sp = _tree(seed=8)
    g = simulate_genomes_family(sp, duplication=0.6, initial_families=3, seed=1)
    fam_of = {e.copy: e.family for e in g.edges if e.kind != "loss"}   # every gene's birth family
    kids = collections.defaultdict(list)
    for e in g.edges:
        if e.kind == "duplication":
            assert e.parent in fam_of and e.family == fam_of[e.parent]  # parent is a real, same-family gene
            kids[e.parent].append(e.copy)
    assert kids and all(len(cs) == 2 for cs in kids.values())          # each duplication has two descendants


def test_every_born_copy_id_is_unique():
    sp = _tree(seed=8, death=0.5)
    g = simulate_genomes_family(sp, duplication=0.4, loss=0.3, origination=0.6, initial_families=5, seed=2)
    born = [e.copy for e in g.edges if e.kind in ("origination", "duplication")]
    assert len(born) == len(set(born))


def test_family_counts_matches_the_genome():
    sp = _tree(seed=9)
    g = simulate_genomes_family(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=1)
    for node_id in g.genomes:
        assert sum(g.family_counts(node_id).values()) == len(g.genomes[node_id])


# --- empty / validation ----------------------------------------------------

def test_empty_run_has_no_events_or_content():
    sp = _tree(seed=1)
    g = simulate_genomes_family(sp, initial_families=0, seed=1)                            # no families, no rates
    assert g.events == []
    assert all(genome == () for genome in g.genomes.values())


def test_validation():
    sp = _tree(seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_family(sp, initial_families=-1, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_family(sp, initial_families=True, seed=1)     # bool is not a valid count
    with pytest.raises(ValueError):
        simulate_genomes_family(sp, origination=-1.0, initial_families=1, seed=1)   # negative rate (via scope)


# --- modifiers: OnTime (skyline) is wired; the rest are rejected, not silently dropped ---

def test_time_skyline_modifier_is_supported():
    # OnTime reads only `time`, which the walk supplies, so a skyline origination works: the rate
    # drops to 0 at t=1.5, so no family originates after it
    sp = simulate_species_tree(birth=1.0, death=0.2, total_time=4.0, seed=3)
    r = simulate_genomes_family(sp, origination=1.0 * mod.OnTime({0: 1.0, 1.5: 0.0}), seed=1)
    orig_times = [e.time for e in r.events if e.kind == "origination"]
    assert orig_times and max(orig_times) < 1.5


def test_unsupported_modifiers_are_rejected_not_silently_dropped():
    sp = _tree(seed=1)
    # clade drift would need per-lineage threading the walk doesn't do → reject, don't no-op
    with pytest.raises(ValueError, match="does not support"):
        simulate_genomes_family(sp, duplication=0.5 * mod.Inherited(per='lineage', spread=0.8), initial_families=3, seed=1)
    # OnTotalDiversity reads a `diversity` context the genome walk doesn't supply → reject, don't crash raw
    with pytest.raises(ValueError, match="does not support"):
        simulate_genomes_family(sp, loss=0.25 * mod.OnTotalDiversity(cap=100), initial_families=3, seed=1)


def test_non_default_scope_is_rejected_not_silently_mismatched():
    # a non-default scope sets the total rate one way while the engine still picks the affected
    # copy/lineage the default way — reject it (a PerCopy origination would be base×0 copies, a no-op)
    from zombi2.rates import scope
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="scope overrides are a later slice"):
        simulate_genomes_family(sp, origination=scope.PerCopy(2.0), seed=1)
    with pytest.raises(ValueError, match="scope overrides are a later slice"):
        simulate_genomes_family(sp, duplication=scope.PerLineage(0.5), initial_families=3, seed=1)
    with pytest.raises(ValueError, match="scope overrides are a later slice"):
        simulate_genomes_family(sp, loss=scope.Global(0.3), initial_families=3, seed=1)
    # the defaults — bare number and the explicit default scope — are accepted
    simulate_genomes_family(sp, origination=scope.PerLineage(0.5), duplication=scope.PerCopy(0.5),
                            initial_families=1, seed=1)


# --- transfer: horizontal moves between contemporaneous lineages -----------

def _alive_at(tree, node_id, t):
    n = tree.nodes[node_id]
    return n.birth_time <= t <= n.end_time


def test_transfer_events_are_contemporaneous_donor_to_recipient():
    sp = _tree(seed=3, death=0.4)
    g = simulate_genomes_family(sp, transfer=0.4, origination=0.5, initial_families=4, seed=1)
    xfers = _transfers(g.edges)
    assert xfers
    for donor, recipient, t in xfers:
        assert donor != recipient                               # a different lineage (default)
        assert _alive_at(sp.complete_tree, donor, t)            # donor alive at t
        assert _alive_at(sp.complete_tree, recipient, t)        # recipient alive at t


def test_transfer_copy_descends_from_a_real_donor_copy():
    sp = _tree(seed=8)
    g = simulate_genomes_family(sp, transfer=0.5, initial_families=5, seed=2)
    born = {e.copy for e in g.edges if e.kind != "loss"}       # every gene id that was born
    xfer_rows = [e for e in g.edges if e.kind == "transfer"]
    assert xfer_rows
    for e in xfer_rows:
        assert e.parent in born and e.copy in born             # the donor gene and the new copy are real


def test_only_transfer_events_carry_a_recipient():
    sp = _tree(seed=2)
    g = simulate_genomes_family(sp, duplication=0.2, transfer=0.3, loss=0.2, origination=0.4,
                                initial_families=4, seed=1)
    for e in g.edges:
        if e.recipient is not None:
            assert e.kind == "transfer"                        # only a transfer names a recipient


def test_no_transfer_at_zero_rate():
    sp = _tree(seed=1)
    g = simulate_genomes_family(sp, duplication=0.3, loss=0.2, origination=0.4, initial_families=5, seed=1)
    assert all(e.kind != "transfer" for e in g.edges)


def test_transfer_is_deterministic():
    sp = _tree(seed=2)
    kw = dict(duplication=0.2, transfer=0.4, loss=0.2, origination=0.4, initial_families=5, seed=9)
    a, b = simulate_genomes_family(sp, **kw), simulate_genomes_family(sp, **kw)
    assert [str(e) for e in a.events] == [str(e) for e in b.events]
    assert a.genomes == b.genomes


def test_replacement_can_displace_a_resident():
    # with loss=0 the only way a copy is lost is a replacement transfer overwriting a homologous copy;
    # each such loss sits at the same instant as a transfer
    sp = _tree(seed=6)
    g = simulate_genomes_family(sp, transfer=1.0, loss=0.0, replacement=True, initial_families=8, seed=1)
    losses = [e for e in g.edges if e.kind == "loss"]
    xfer_times = {e.time for e in g.edges if e.kind == "transfer"}
    assert losses                                               # replacement did displace some copies
    assert all(e.time in xfer_times for e in losses)           # every loss co-occurs with a transfer


def test_additive_transfer_never_loses():
    sp = _tree(seed=6)
    g = simulate_genomes_family(sp, transfer=1.0, loss=0.0, replacement=False, initial_families=8, seed=1)
    assert all(e.kind != "loss" for e in g.events)             # additive transfer only ever adds


def test_default_transfer_is_never_self_but_self_transfer_runs():
    sp = _tree(seed=4)
    g = simulate_genomes_family(sp, transfer=0.6, initial_families=5, seed=1)
    assert _transfers(g.edges)
    assert all(donor != recipient for donor, recipient, _ in _transfers(g.edges))
    s = simulate_genomes_family(sp, transfer=0.6, self_transfer=True, initial_families=5, seed=1)
    assert any(e.kind == "transfer" for e in s.edges)         # runs (self donor==recipient now allowed)


def test_distance_mode_runs_and_is_deterministic():
    sp = _tree(seed=7, death=0.4)
    from zombi2.genomes import Distance
    a = simulate_genomes_family(sp, transfer=0.5, transfer_to="distance", initial_families=5, seed=3)
    b = simulate_genomes_family(sp, transfer=0.5, transfer_to=Distance(decay=1.0), initial_families=5, seed=3)
    assert [str(e) for e in a.edges] == [str(e) for e in b.edges]  # "distance" == Distance(decay=1.0)
    assert any(e.kind == "transfer" for e in a.edges)


def test_transfer_can_come_from_the_dead():
    # high death + transfer: some transfers are donated by a lineage that later goes extinct
    sp = _tree(seed=3, death=0.7)
    g = simulate_genomes_family(sp, transfer=0.8, origination=0.6, initial_families=4, seed=2)
    donor_fates = {sp.complete_tree.nodes[e.lineage].fate for e in g.edges if e.kind == "transfer"}
    assert "extinct" in donor_fates


def test_transfer_to_validation():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="transfer_to"):
        simulate_genomes_family(sp, transfer=0.3, transfer_to="bogus", initial_families=3, seed=1)


# --- transfer_to = Clades: weight recipients by the donor's and recipient's named clade -----------

def _two_clades(sp):
    """Two disjoint internal clades of the complete tree, plus their tip lists and a labeler that maps
    a node id to ``"A"`` / ``"B"`` / ``"rest"``. Chosen to leave a real 'rest' remainder."""
    tree = sp.complete_tree

    def desc(r):
        out, st = set(), [r]
        while st:
            i = st.pop(); out.add(i)
            if tree.nodes[i].children:
                st.extend(tree.nodes[i].children)
        return out

    internals = [i for i, n in tree.nodes.items() if n.children is not None and i != tree.root]
    sub = {i: desc(i) for i in internals}
    for a in internals:
        if not (4 <= len(sub[a]) <= len(tree.nodes) // 3):
            continue
        for b in internals:
            if a != b and not (sub[a] & sub[b]) and 4 <= len(sub[b]) <= len(tree.nodes) // 3:
                lab = lambda i, A=sub[a], B=sub[b]: "A" if i in A else ("B" if i in B else "rest")
                tips_a = [i for i in sub[a] if tree.nodes[i].children is None]
                tips_b = [i for i in sub[b] if tree.nodes[i].children is None]
                return a, b, tips_a, tips_b, lab
    raise AssertionError("no two disjoint clades of a usable size")


def _clade_pairs(events, lab):
    """Every fired transfer as a ``(donor clade, recipient clade)`` label pair."""
    return [(lab(e.donor), lab(e.lineage)) for e in events
            if e.kind == "transfer" and e.recipient is not None]


# bounded params: a small tree and a low per-genome family cap keep additive transfers from growing
# copies to the cap (which would make the per-copy transfer rate — and the run — explode), while still
# firing enough transfers for an exact who-sends-to-whom assertion.
_CLADE_TREE = dict(seed=7, n_extant=20)
_CLADE_KW = dict(transfer=1.0, initial_families=15, max_family_size=8, seed=11)


def test_clades_between_only_excludes_within_and_rest():
    """The case a per-recipient weight cannot express: transfer runs strictly BETWEEN A and B —
    never A→A, B→B, or anything touching the rest of the tree."""
    from zombi2.genomes import Between, Clades
    sp = _tree(**_CLADE_TREE)
    a, b, _, _, lab = _two_clades(sp)
    g = simulate_genomes_family(
        sp, transfer_to=Clades({"A": a, "B": b}, Between({("A", "B"): 1.0, ("B", "A"): 1.0},
                                                         default=0.0)), **_CLADE_KW)
    pairs = _clade_pairs(g.edges, lab)
    assert pairs
    assert all(p in {("A", "B"), ("B", "A")} for p in pairs)


def test_clades_directional_a_to_b_only():
    """A donates, B receives, nothing else: every transfer is A→B."""
    from zombi2.genomes import Between, Clades
    sp = _tree(**_CLADE_TREE)
    a, b, _, _, lab = _two_clades(sp)
    g = simulate_genomes_family(
        sp, transfer_to=Clades({"A": a, "B": b}, Between({("A", "B"): 1.0}, default=0.0)), **_CLADE_KW)
    pairs = _clade_pairs(g.edges, lab)
    assert pairs and all(p == ("A", "B") for p in pairs)


def test_clades_named_by_tips_equals_named_by_node_id():
    """A clade named by all its tips (its MRCA's subtree) is the same clade as its node id — the two
    spellings give a byte-identical run."""
    from zombi2.genomes import Between, Clades
    sp = _tree(**_CLADE_TREE)
    a, b, tips_a, tips_b, _ = _two_clades(sp)
    kernel = Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
    by_id = simulate_genomes_family(sp, transfer_to=Clades({"A": a, "B": b}, kernel), **_CLADE_KW)
    by_tips = simulate_genomes_family(sp, transfer_to=Clades({"A": tips_a, "B": tips_b}, kernel),
                                      **_CLADE_KW)
    assert [str(e) for e in by_id.events] == [str(e) for e in by_tips.events]


def test_clades_default_baseline_allows_other_pairs():
    """default=1.0 is a baseline: an up-weighted pair is enriched, but unlisted pairs still happen —
    unlike default=0.0, which forbids them."""
    from zombi2.genomes import Between, Clades
    sp = _tree(**_CLADE_TREE)
    a, b, _, _, lab = _two_clades(sp)
    g = simulate_genomes_family(
        sp, transfer_to=Clades({"A": a, "B": b}, Between({("A", "B"): 8.0})), **_CLADE_KW)  # default 1.0
    kinds = set(_clade_pairs(g.edges, lab))
    assert ("A", "B") in kinds
    assert any(p not in {("A", "B"), ("B", "A")} for p in kinds)   # within/rest still occur at baseline


def test_clades_is_deterministic():
    from zombi2.genomes import Between, Clades
    sp = _tree(**_CLADE_TREE)
    a, b, _, _, _ = _two_clades(sp)
    kw = dict(_CLADE_KW, transfer_to=Clades({"A": a, "B": b},
                                            Between({("A", "B"): 1.0, ("B", "A"): 1.0})))
    x, y = simulate_genomes_family(sp, **kw), simulate_genomes_family(sp, **kw)
    assert [str(e) for e in x.events] == [str(e) for e in y.events]


def test_clades_construction_validation():
    from zombi2.genomes import Between, Clades
    sp = _tree(seed=1)
    root_kid = sp.complete_tree.nodes[sp.complete_tree.root].children[0]
    with pytest.raises(ValueError, match="reserved"):
        Clades({"rest": root_kid}, Between({("rest", "rest"): 1.0}))
    with pytest.raises(ValueError, match="not defined clades"):
        Clades({"A": root_kid}, Between({("A", "Z"): 1.0}))
    with pytest.raises(ValueError, match="Between kernel"):
        Clades({"A": root_kid}, "nope")
    with pytest.raises(ValueError, match="non-empty"):
        Clades({}, Between({("A", "B"): 1.0}))


def test_clades_runtime_validation():
    from zombi2.genomes import Between, Clades
    sp = _tree(seed=1)
    tree = sp.complete_tree
    root_kid = tree.nodes[tree.root].children[0]
    grandkid = tree.nodes[root_kid].children[0]           # nested inside root_kid
    with pytest.raises(ValueError, match="disjoint"):     # nested clades overlap
        simulate_genomes_family(sp, transfer=1.0, initial_families=3,
                                transfer_to=Clades({"A": root_kid, "B": grandkid},
                                                      Between({("A", "B"): 1.0})))
    with pytest.raises(ValueError, match="not a lineage"):
        simulate_genomes_family(sp, transfer=1.0, initial_families=3,
                                transfer_to=Clades({"A": 999999}, Between({("A", "A"): 1.0})))


def test_distance_decay_validation():
    from zombi2.genomes import Distance
    Distance(decay=0.0)                       # zero is fine — the uniform limit
    for bad in (-1.0, float("inf"), float("nan"), True):
        with pytest.raises(ValueError, match="Distance decay"):
            Distance(decay=bad)


# --- the written outputs ---------------------------------------------------

def _rows(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    cols = lines[0].split("\t")
    return cols, [dict(zip(cols, ln.split("\t"))) for ln in lines[1:] if ln]


def test_written_node_columns_carry_a_lineage_label():
    # species_events.tsv and trait_values.tsv have always written n<id>; the genome tables used to
    # write bare ints, so the same node read two ways in one output directory. A lineage that went
    # extinct is e<id>, so the test is that a node carries a LETTER, not that it carries "n".
    # The event log has no lineage column any more: a participant is n<species>_g<copy>, carrying
    # the branch it lived on inside the token, so the letter has to be there.
    sp = _tree(seed=2)
    g = simulate_genomes_family(sp, duplication=0.3, transfer=0.3, loss=0.2, origination=0.5,
                                initial_families=4, seed=7)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d)
        g.write(out, outputs=("events", "genomes"))
        _, events = _rows(out / "genome_events.tsv")
        copies = [tok for r in events for col in ("parents", "children")
                  for tok in r[col].split(";") if tok]
        assert copies
        for tok in copies:
            species, _, gene = tok.partition("_")
            assert species[:1] in "ne" and species[1:].isdigit(), tok
            assert gene[:1] == "g" and gene[1:].isdigit(), tok
        _, genomes = _rows(out / "genomes.tsv")
        assert genomes and all(r["lineage"][:1] in "ne" for r in genomes)


def test_written_log_round_trips_through_the_reader():
    # the written row is one EVENT and an Event is one gene-tree EDGE, so the reader has to expand
    # them again — exactly, or a downstream level replaying the log builds a different history
    from zombi2.genomes.events import edges_from_tsv, events_tsv
    sp = _tree(seed=2)
    for kw in ({}, {"replacement": True}, {"replacement": True, "self_transfer": True}):
        g = simulate_genomes_family(sp, duplication=0.3, transfer=0.3, loss=0.2, origination=0.5,
                                    initial_families=4, seed=7, **kw)
        kinds = {e.kind for e in g.edges}
        assert kinds == {"origination", "duplication", "loss", "transfer", "speciation"}, kw
        names = sp.complete_tree.labels()
        assert edges_from_tsv(events_tsv(g.edges, names)) == g.edges, kw


def test_write_genomes_covers_every_node_where_profiles_covers_only_tips():
    sp = _tree(seed=5, death=0.6)
    g = simulate_genomes_family(sp, duplication=0.3, loss=0.3, origination=0.6,
                                initial_families=5, seed=3)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d)
        g.write(out, outputs=("genomes", "profiles"))
        _, rows = _rows(out / "genomes.tsv")
        written = {r["lineage"] for r in rows}
        names = sp.complete_tree.labels()
        assert written == {names[s] for s in g.genomes if g.genomes[s]}
        internal = {names[n.id] for n in sp.complete_tree.nodes.values() if n.children is not None}
        assert written & internal, "ancestral genomes must be in there, not just the tips"
        # profiles is the extant-only view
        tips = {names[n.id] for n in sp.complete_tree.extant_leaves()}
        assert set((out / "profiles.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")[1:]) == tips


def test_write_gene_trees_emits_one_newick_per_family():
    sp = _tree(seed=5, death=0.6)
    g = simulate_genomes_family(sp, duplication=0.3, loss=0.3, origination=0.6,
                                initial_families=5, seed=3)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d)
        g.write(out, outputs=("gene_trees",))
        for fam, gt in g.gene_trees.items():
            complete = out / "gene_trees" / f"gene_tree_fam{fam}_complete.nwk"
            # the written tree names a dead lineage e<id>, which needs the run's tree to know
            assert complete.read_text(encoding="utf-8").strip() == \
                gt.to_newick("complete", names=sp.complete_tree.labels())
            extant = out / "gene_trees" / f"gene_tree_fam{fam}_extant.nwk"
            # a family with no survivor has no extant tree, and writes no file for it
            assert extant.exists() == (gt.to_newick("extant") is not None)


# --- ByFamily: per-family rate heterogeneity -------------------------------

def _dup_per_family(g, n_families):
    counts = collections.Counter(e.family for e in g.events if e.kind == "duplication")
    return [counts.get(f, 0) for f in range(n_families)]


def test_by_family_spreads_the_rates_without_moving_their_mean():
    # the point of the modifier: families stop being interchangeable. The draw is mean-corrected,
    # so widening the spread must widen the spread of outcomes without inflating the average.
    sp = _tree(seed=7, n_extant=20, death=0.0)
    flat = simulate_genomes_family(sp, duplication=0.25, loss=0.25, initial_families=150, seed=3)
    varied = simulate_genomes_family(sp, duplication=0.25 * mod.Drawn(per='family', spread=0.5),
                                     loss=0.25, initial_families=150, seed=3)
    f, v = _dup_per_family(flat, 150), _dup_per_family(varied, 150)
    import statistics
    assert statistics.pstdev(v) > 1.5 * statistics.pstdev(f)      # families genuinely differ
    assert statistics.mean(v) == pytest.approx(statistics.mean(f), rel=0.35)   # mean is held


def test_a_run_with_no_by_family_is_untouched():
    # the weighted path costs something, so it must only be taken when it is asked for
    sp = _tree(seed=2, n_extant=12)
    a = simulate_genomes_family(sp, duplication=0.2, transfer=0.1, loss=0.2,
                                initial_families=10, seed=5)
    b = simulate_genomes_family(sp, duplication=0.2, transfer=0.1, loss=0.2,
                                initial_families=10, seed=5)
    assert [(e.time, e.kind, e.copy) for e in a.edges] == [(e.time, e.kind, e.copy) for e in b.edges]


def test_by_family_is_deterministic_given_the_seed():
    sp = _tree(seed=2, n_extant=12)
    kw = dict(duplication=0.2 * mod.Drawn(per='family', spread=0.6), loss=0.2, initial_families=20, seed=5)
    a = simulate_genomes_family(sp, **kw)
    b = simulate_genomes_family(sp, **kw)
    assert [(e.time, e.kind, e.copy) for e in a.edges] == [(e.time, e.kind, e.copy) for e in b.edges]


def test_one_shared_draw_moves_every_rate_of_a_family_together():
    # the other placement: ONE ByFamily object read by both rates, so one draw per family scales
    # them together. A family that duplicates a lot should also be losing a lot — which is exactly
    # what two separately built ByFamily draws do NOT give.
    speed = mod.Drawn(per='family', spread=0.6)
    sp = _tree(seed=1, n_extant=20, death=0.0)
    g = simulate_genomes_family(sp, duplication=0.25 * speed, loss=0.25 * speed,
                                initial_families=150, seed=3)
    dup = collections.Counter(e.family for e in g.events if e.kind == "duplication")
    los = collections.Counter(e.family for e in g.events if e.kind == "loss")
    fams = [f for f in range(150) if dup.get(f, 0) + los.get(f, 0) > 4]
    assert len(fams) > 20
    xs = [dup.get(f, 0) for f in fams]
    ys = [los.get(f, 0) for f in fams]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sx = sum((a - mx) ** 2 for a in xs) ** 0.5
    sy = sum((b - my) ** 2 for b in ys) ** 0.5
    assert cov / (sx * sy) > 0.3          # a fast family is fast at everything


def test_the_carried_family_weights_match_a_full_recompute(monkeypatch):
    # The per-lineage multiplier sums are carried across events and only the touched lineage is
    # rebuilt, which is what keeps the weighted path off a quadratic. That is only sound while the
    # carried arrays say what a full recompute would say — so check exactly that, on every event of
    # a real run, with every mutation in play (duplication, loss, origination, transfer with
    # replacement, and a cap that makes some events no-ops).
    from zombi2.genomes.family import _FamilyWeights, _sum_mult   # (the module, not genomes.family)

    real, checked = _FamilyWeights.current, []

    def current(self, gen):
        out = real(self, gen)
        for m, _keys, arr in self._groups:
            assert arr == [_sum_mult(m, g) for g in gen]
        checked.append(len(gen))
        return out

    monkeypatch.setattr(_FamilyWeights, "current", current)
    speed = mod.Drawn(per='family', spread=0.7)
    sp = _tree(seed=4, n_extant=16, death=0.3)
    g = simulate_genomes_family(sp, duplication=0.3 * speed, transfer=0.2 * speed,
                                loss=0.3 * speed, origination=0.4, replacement=True,
                                max_family_size=4, initial_families=40, seed=6)
    assert len(checked) > 500                      # the run really did exercise the loop
    assert {e.kind for e in g.edges} == {"origination", "duplication", "loss", "transfer",
                                          "speciation"}


def test_by_family_is_refused_on_origination():
    sp = _tree(seed=1, n_extant=8)
    with pytest.raises(ValueError, match="families are CREATED"):
        simulate_genomes_family(sp, origination=0.5 * mod.Drawn(per='family', spread=0.3), seed=1)


def test_by_family_with_driven_by_is_refused_for_now():
    sp = _tree(seed=1, n_extant=8)
    with pytest.raises(ValueError, match="later slice"):
        simulate_genomes_family(sp, loss=0.2 * mod.Drawn(per='family', spread=0.3),
                                duplication=0.2 * ScaledBy("x.tsv", {"a": 2.0}), seed=1)


# --- max_family_size: a per-genome ceiling on a family's copies ------------

def _biggest_family(g):
    return max((collections.Counter(c.family for c in gen).most_common(1) or [(None, 0)])[0][1]
               for gen in g.genomes.values() if gen)


def test_max_family_size_binds_exactly():
    # duplication far above loss compounds without bound; the cap is what stops it, and it stops it
    # AT the number given rather than somewhere near it
    sp = _tree(seed=1, n_extant=12, death=0.0)
    for cap in (2, 5, 9):
        g = simulate_genomes_family(sp, duplication=0.9, loss=0.05, initial_families=6,
                                    max_family_size=cap, seed=2)
        assert _biggest_family(g) == cap


def test_max_family_size_none_lifts_the_ceiling():
    sp = _tree(seed=1, n_extant=12, death=0.0)
    capped = simulate_genomes_family(sp, duplication=0.9, loss=0.05, initial_families=6,
                                     max_family_size=5, seed=2)
    free = simulate_genomes_family(sp, duplication=0.9, loss=0.05, initial_families=6,
                                   max_family_size=None, seed=2)
    assert _biggest_family(free) > _biggest_family(capped)


def test_the_cap_is_a_plain_count_and_does_not_move_with_the_tree():
    # the bug this replaced: PerLineage(n) multiplied the cap by the node count of the species tree,
    # so a per-GENOME bound grew every time you added species — the shipped PerLineage(10) was 1470
    # copies of one family in one genome on a 147-node tree
    from zombi2.genomes import resolve_max_family_size
    assert resolve_max_family_size(7) == 7
    assert resolve_max_family_size(None) is None
    small = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1)
    big = simulate_species_tree(birth=1.0, death=0.0, n_extant=40, seed=1)
    kw = dict(duplication=0.9, loss=0.05, initial_families=6, max_family_size=4, seed=2)
    assert _biggest_family(simulate_genomes_family(small, **kw)) == 4      # the same number...
    assert _biggest_family(simulate_genomes_family(big, **kw)) == 4        # ...on either tree


def test_the_old_scope_spelling_is_refused_with_the_arithmetic():
    # a stale script must fail loudly rather than run at a cap different from the one it reads
    from zombi2.genomes import resolve_max_family_size
    with pytest.raises(ValueError, match="max_family_size=10.*size of the species tree"):
        resolve_max_family_size(PerLineage(10))
    with pytest.raises(ValueError, match="max_family_size=10.*exactly n"):
        resolve_max_family_size(Global(10))


def test_max_family_size_is_validated():
    from zombi2.genomes import resolve_max_family_size
    for bad in (0, -1):
        with pytest.raises(ValueError, match="at least 1"):
            resolve_max_family_size(bad)
    # a float is refused rather than rounded: 10 against 10.0 used to be two different caps, and
    # True is an int to Python but not a copy count to anybody else
    for bad in (10.0, 2.5, True, "big"):
        with pytest.raises(ValueError, match="whole number"):
            resolve_max_family_size(bad)
    with pytest.raises(ValueError, match="plain count"):      # a scope of any kind
        resolve_max_family_size(PerCopy(5))


def test_the_cap_also_holds_when_a_transfer_arrives():
    # a transfer adds a copy to the recipient, so the ceiling has to hold there too or a family
    # could be pushed past it sideways
    sp = _tree(seed=3, n_extant=15, death=0.0)
    g = simulate_genomes_family(sp, duplication=0.4, transfer=0.6, loss=0.05,
                                initial_families=8, max_family_size=4, seed=1)
    assert _biggest_family(g) <= 4


def test_the_carried_family_counts_match_a_full_scan(monkeypatch):
    # The cap's question — does this family already fill its quota here? — is answered from a counter
    # rather than by scanning the lineage's genome, which is what keeps the level linear in genome
    # size. Sound only while the counter says what the scan would, so check it against the scan on
    # every call of a real run, with every mutation in play.
    from zombi2.genomes.family import _FamilyCounts, _at_cap

    real_init, real_at_cap, checked, lent = _FamilyCounts.__init__, _FamilyCounts.at_cap, [], {}

    def init(self, gen):                    # the class does not hold `gen`; lend it for the check
        real_init(self, gen)
        lent[id(self)] = gen                # the same list objects the engine mutates in place

    def at_cap(self, k, family, cap):
        out = real_at_cap(self, k, family, cap)
        assert out == _at_cap(lent[id(self)][k], family, cap), (k, family, cap)
        checked.append(k)
        return out

    monkeypatch.setattr(_FamilyCounts, "__init__", init)
    monkeypatch.setattr(_FamilyCounts, "at_cap", at_cap)

    sp = _tree(seed=4, n_extant=14, death=0.3)
    g = simulate_genomes_family(sp, duplication=0.5, transfer=0.4, loss=0.3, origination=0.4,
                                initial_families=25, max_family_size=3, replacement=True,
                                self_transfer=True, seed=9)
    assert len(checked) > 200                       # the run really did exercise the cap
    assert {e.kind for e in g.edges} == {"origination", "duplication", "loss", "transfer",
                                          "speciation"}


def test_the_cap_binds_at_the_number_given_however_a_copy_arrives():
    # duplication is not the only way a family grows: a transfer arrives into a genome too, and the
    # cap has to hold for both
    sp = _tree(seed=2, n_extant=16, death=0.2)
    for cap in (1, 2, 5):
        g = simulate_genomes_family(sp, duplication=0.6, transfer=0.5, loss=0.1, origination=0.3,
                                    initial_families=20, max_family_size=cap, seed=3)
        biggest = max((max(collections.Counter(c.family for c in genome).values(), default=0)
                       for genome in g.genomes.values()), default=0)
        assert biggest == cap, f"cap {cap} gave a family of {biggest}"
