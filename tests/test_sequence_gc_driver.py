"""Sequence composition as a conditioning driver — a finished sequence making a trait's rate faster.

``gc()`` is the named door and ``composition(letters)`` the general one, so an amino-acid frequency is
the same driver with different letters; the second section below is where that is checked.

Three things are worth reading. The **cross-check**: the driver pools sequences by the lineage named
in each gene's label, so it is graded against a recount off the alignments. The **direction**: a
genome reading a sequence must be refused by name. And the **gap**: a lineage no family reached
carries no DNA, and has to answer anyway.
"""

from __future__ import annotations

import pytest

from zombi2 import traits
from zombi2.genomes import (simulate_genomes_family, simulate_genomes_nucleotide,
                            simulate_genomes_ordered)
from zombi2.params import PerCopy, PerLineage, PerSite, Recipients
from zombi2.params.mapping import Curve
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, jc69, lg
from zombi2.species import simulate_species_tree

#: skewed towards A/T, so GC is somewhere to move *to* and the tips are not all one number
_SKEWED = hky85(2.0, frequencies=(0.35, 0.15, 0.15, 0.35))
_RESPONSE = Curve(lambda x: 30.0 ** (x - 0.4))       # steep, so a driven run is visibly not the plain one


def _run(*, loss=0.05, families=6, length=200, model=None, n_extant=20, seed=6, genome_seed=5):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=4)
    g = simulate_genomes_family(sp.complete_tree, initial_families=families, duplication=0.1,
                                loss=loss, seed=genome_seed)
    return sp.complete_tree, simulate_sequences(g, model=model or _SKEWED, length=length, seed=seed)


def _recount(result, tree, letters: str = "GC") -> dict[int, tuple[int, int]]:
    """``{lineage: (hits, sites)}``, counted here rather than by the driver — the independently
    derived fact its per-node values are graded against."""
    out: dict[int, list[int]] = {}
    for by_unit in (result.alignments, result.ancestral):
        for sequences in by_unit.values():
            for label, seq in sequences.items():
                lineage = int(label.rsplit("_g", 1)[0][1:])   # n7_g12 / e7_g12 -> 7
                assert lineage in tree.nodes
                counts = out.setdefault(lineage, [0, 0])
                counts[0] += sum(seq.count(x) for x in letters)
                counts[1] += len(seq)
    return {k: (a, b) for k, (a, b) in out.items()}


def test_a_node_s_gc_is_its_whole_complement_pooled():
    """Checked against a recount: every family's sequence at that node, extant tips and ancestors
    alike, over the total sites — not one family, and not the tips alone."""
    tree, seqs = _run()
    values = seqs.gc()._node_values(tree)
    counted = _recount(seqs, tree)
    assert len(counted) > 40, "this run is meant to put sequences on most of the tree"
    for lineage, (gc, sites) in counted.items():
        assert values[lineage] == pytest.approx(gc / sites)
    assert min(values.values()) < max(values.values()), "GC is meant to vary across the tree"


def test_it_is_pooled_over_families_not_taken_from_one():
    """Were it reading one family, adding five more would not move the number."""
    tree, one = _run(families=1)
    _t, many = _run(families=8)
    a, b = one.gc()._node_values(tree), many.gc()._node_values(tree)
    root = tree.root
    assert a[root] != b[root]


def test_every_lineage_gets_an_answer():
    tree, seqs = _run()
    traj = seqs.gc().as_driver_trajectory(tree)
    assert set(traj._starts) == set(tree.nodes)
    for node in tree.nodes.values():
        assert 0.0 <= float(traj.value(node.id, node.birth_time)) <= 1.0


def test_a_lineage_with_no_dna_holds_its_parent_s_value():
    """One family, lost early, so whole clades have no sequence at all. A driver still has to answer
    for those branches, and 0.0 would drive them as though the genome had turned pure AT."""
    tree, seqs = _run(families=1, loss=0.2, genome_seed=8)
    values = seqs.gc()._node_values(tree)
    counted = _recount(seqs, tree)
    empty = [i for i in tree.nodes if i not in counted]
    assert len(empty) > 5, "this run is meant to leave some lineages with no sequence"
    assert len(counted) > 10, "...and to leave most of them with one"
    for i in empty:
        parent = tree.nodes[i].parent
        assert parent is not None
        assert values[i] == values[parent]


def test_the_path_between_two_nodes_is_the_line_between_them():
    """The stretches lie between the branch's endpoint values, and a finer `step` cuts more."""
    tree, seqs = _run()
    values = seqs.gc()._node_values(tree)
    node = max((n for n in tree.nodes.values() if n.parent is not None),
               key=lambda n: n.end_time - n.birth_time)
    coarse = seqs.gc().as_driver_trajectory(tree, step=10.0)
    fine = seqs.gc().as_driver_trajectory(tree, step=0.01)
    lo, hi = sorted((values[node.parent], values[node.id]))
    for traj in (coarse, fine):
        for value in traj._states[node.id]:
            assert lo <= float(value) <= hi
    assert len(fine._starts[node.id]) > len(coarse._starts[node.id])


def test_gc_drives_a_trait():
    tree, seqs = _run()
    kw = dict(states=["mesophile", "thermophile"], start="mesophile", seed=2)
    plain = traits.simulate_discrete(tree, switch=0.3, **kw)
    driven = traits.simulate_discrete(tree, switch=PerLineage(0.3).scaled_by(seqs.gc(), _RESPONSE), **kw)
    switches = lambda r: sum(1 for e in r.events if e.kind == "on_branch")
    assert switches(driven) != switches(plain)
    # ...and the continuous engine reads it the same way
    a = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=3)
    b = traits.simulate_continuous(tree, start=0.0, seed=3,
                                   rate=PerLineage(1.0).scaled_by(seqs.gc(), _RESPONSE))
    assert a.values != b.values


def test_one_sequence_run_drives_another():
    """Conditioning within a level (SPEC §3): allowed when the driver can be finished first."""
    tree, seqs = _run()
    g = simulate_genomes_family(tree, initial_families=2, duplication=0.0, loss=0.0, seed=4)
    plain = simulate_sequences(g, model=jc69(), length=60, seed=4)
    driven = simulate_sequences(g, model=jc69(), length=60, seed=4,
                                substitution=PerSite(1.0).scaled_by(seqs.gc(), _RESPONSE))
    assert plain.species_phylogram["complete"] != driven.species_phylogram["complete"]


# --- what it refuses, and why -------------------------------------------------------------------

def test_no_genome_can_be_driven_by_a_sequence():
    """All three resolutions produce the gene trees a sequence is grown along, so reading one back
    conditions a run on its own output. The refusal has to say that, not fail downstream."""
    tree, seqs = _run()
    gc = seqs.gc()
    driven = PerCopy(0.2).scaled_by(gc, _RESPONSE)
    for run in (lambda: simulate_genomes_family(tree, initial_families=3, loss=driven, seed=7),
                lambda: simulate_genomes_ordered(tree, initial_families=3, loss=driven, seed=7),
                lambda: simulate_genomes_nucleotide(
                    tree, root_length=2000, genes=2,
                    loss=PerLineage(0.2).scaled_by(gc, _RESPONSE), seed=7)):
        with pytest.raises(ValueError, match="genome level cannot be driven by a sequence"):
            run()
    # who RECEIVES a transfer is a read of the driver too, so the same refusal covers it
    with pytest.raises(ValueError, match="genome level cannot be driven by a sequence"):
        simulate_genomes_family(tree, initial_families=3, transfer=0.2, seed=7,
                                transfer_to=Recipients().weighted_by(gc, _RESPONSE))


def test_a_joint_run_and_the_species_level_refuse_it_too():
    """Neither reaches the driver, and both say why on their own terms."""
    from zombi2.joint import simulate_joint

    tree, seqs = _run()
    driven = PerLineage(1.0).scaled_by(seqs.gc(), _RESPONSE)
    with pytest.raises(TypeError, match="live level"):
        simulate_joint(birth=driven, trait=traits.discrete(states=["a", "b"], switch=0.3),
                       n_extant=10, seed=4)
    with pytest.raises(ValueError, match="species engine does not support"):
        simulate_species_tree(birth=driven, death=0.2, n_extant=8, seed=4)


def test_gc_is_a_number_so_a_state_table_is_refused():
    tree, seqs = _run()
    with pytest.raises(ValueError, match="driver is CONTINUOUS"):
        traits.simulate_discrete(tree, states=["a", "b"], seed=2,
                                 switch=PerLineage(0.3).scaled_by(seqs.gc(), {"0.5": 2.0}))


def test_one_family_s_gc_is_refused_rather_than_offered():
    tree, seqs = _run()
    with pytest.raises(ValueError, match="pooled over every family"):
        seqs.gc("fam3")


# --- the general door: any letters of the run's own alphabet ------------------------------------

def test_gc_is_composition_over_g_and_c():
    """One statistic, two front doors — so the named one cannot drift from the general one."""
    tree, seqs = _run()
    assert seqs.gc() == seqs.composition("GC")
    assert seqs.gc()._node_values(tree) == seqs.composition("GC")._node_values(tree)


def test_an_amino_acid_frequency_is_the_same_driver():
    """The plug-in case: a residue set on a protein run, counted and driven exactly as GC is."""
    tree, seqs = _run(model=lg(), length=300)
    basic = seqs.composition("KR")
    counted = _recount(seqs, tree, letters="KR")
    for lineage, (hits, sites) in counted.items():
        assert basic._node_values(tree)[lineage] == pytest.approx(hits / sites)
    kw = dict(states=["mesophile", "thermophile"], start="mesophile", seed=2)
    plain = traits.simulate_discrete(tree, switch=0.3, **kw)
    driven = traits.simulate_discrete(
        tree, switch=PerLineage(0.3).scaled_by(basic, Curve(lambda x: 1.0 + 40.0 * x)), **kw)
    switches = lambda r: sum(1 for e in r.events if e.kind == "on_branch")
    assert switches(driven) != switches(plain)


def test_a_composition_refuses_letters_the_run_cannot_contain():
    """A set that occurs nowhere would read 0.0 on every lineage — an undriven run wearing a driven
    rate, which is the silence SPEC §5 exists to break."""
    tree, dna = _run()
    _t, protein = _run(model=lg(), length=100)
    with pytest.raises(ValueError, match="not in this run's alphabet"):
        dna.composition("KR")
    with pytest.raises(ValueError, match="at least one letter"):
        dna.composition("")
    with pytest.raises(TypeError, match="as a string"):
        dna.composition(["K"])
    assert protein.composition("kr").letters == "KR"        # case is not the user's problem


def test_gc_is_refused_on_a_protein_run():
    """G and C are glycine and cysteine there, so the call is ambiguous rather than wrong."""
    tree, seqs = _run(model=lg(), length=100)
    with pytest.raises(ValueError, match="needs DNA"):
        seqs.gc()
    seqs.composition("GC")                                  # ...and saying so by name is allowed


def test_it_refuses_a_tree_it_was_not_grown_on():
    tree, seqs = _run()
    other = simulate_species_tree(birth=1.0, death=0.2, n_extant=9, seed=77).complete_tree
    with pytest.raises(ValueError, match="different species trees"):
        seqs.gc().as_driver_trajectory(other)
