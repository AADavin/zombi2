"""A genome and a gene's sequence, each driving the other — the last cell of the map.

See `docs/design/genomes-sequences.md`. The genome decides which sequences exist; the sequences
decide how fast the genome changes. Species time is sliced, because a composition moves with every
substitution and so no genome rate reading it is ever constant.
"""

import statistics

import pytest

from zombi2 import joint, species
from zombi2.genomes import family, genome as genome_spec
from zombi2.params import Curve, PerCopy
from zombi2.sequences import composition, gene, hky85


_AT_RICH = hky85(2.0, frequencies=(0.40, 0.10, 0.10, 0.40))
_EVEN = hky85(2.0)


def _tree(n=25, seed=1):
    return species.simulate_species_tree(birth=1.0, n_extant=n, seed=seed).complete_tree


def _run(ct, *, steep=True, seed=1, step=0.05, transfer=0.0, **kw):
    """hisA arrives AT-rich and ameliorates. While it is AT-rich its lineage loses genes fast, so how
    far it has got sets how much of the genome survives."""
    f = Curve((lambda gc: 30.0 ** ((0.35 - gc) / 0.2)) if steep else (lambda gc: 1.0))
    return joint.simulate(
        genome_spec(duplication=0.15, origination=0.0, initial_families=4, transfer=transfer,
                    loss=PerCopy(0.15).scaled_by("sequences:hisA", f, step=step),
                    families=[family("hisA")]),
        gene(name="hisA", model=_EVEN, length=250, start=_AT_RICH, substitution=0.8,
             offers=composition("GC", absent=0.35)),
        tree=ct, seed=seed, **kw)


def _copies(r):
    fam = r.genome.family_names["hisA"]
    return sum(1 for lin in r.genome.genomes.values() for c in lin if c.family == fam)


# --- what comes back ------------------------------------------------------------------------------

def test_it_returns_both_levels():
    ct = _tree()
    r = _run(ct)
    assert r.genome is not None and r.sequences is not None and r.trait is None
    assert r.complete_tree is ct
    assert r.species.events == []          # the tree came in, so its log is not this run's to write
    assert r.sequences.families == ("hisA",)


def test_every_node_of_the_gene_tree_has_its_sequence():
    """The engine keeps a copy's sequence the moment its branch ends, whatever ended it — a
    duplication, a transfer, a loss, a speciation, or the present."""
    ct = _tree()
    r = _run(ct, transfer=0.1)
    fam = r.genome.family_names["hisA"]
    gt = r.genome.gene_trees[fam]
    nodes, stack = 0, [gt.complete]
    while stack:
        n = stack.pop()
        nodes += 1
        stack.extend(n.children)
    assert nodes == len(r.sequences.alignments[fam]) + len(r.sequences.ancestral[fam])


def test_the_founding_sequence_is_kept():
    ct = _tree()
    r = _run(ct)
    fam = r.genome.family_names["hisA"]
    assert len(r.sequences.founding[fam]) == 250
    assert set(r.sequences.founding[fam]) <= set("ACGT")


def test_only_the_declared_gene_gets_sequences():
    """Every other family races exactly as it does today. That is what keeps the cost proportional
    to what the model reads."""
    ct = _tree()
    r = _run(ct)
    assert len(r.genome.gene_trees) > 1
    assert sorted(r.sequences.alignments) == [r.genome.family_names["hisA"]]


def test_it_is_deterministic():
    ct = _tree()
    a, b = _run(ct, seed=6), _run(ct, seed=6)
    assert a.sequences.alignments == b.sequences.alignments
    assert [(e.time, e.kind) for e in a.genome.events] == \
           [(e.time, e.kind) for e in b.genome.events]


def test_both_levels_write():
    import pathlib
    import tempfile
    ct = _tree()
    with tempfile.TemporaryDirectory() as d:
        _run(ct, seed=3).write(d)
        names = {p.name for p in pathlib.Path(d).rglob("*") if p.is_file()}
    assert "genome_events.tsv" in names and "joint_summary.json" in names


def test_it_records_when_asked():
    ct = _tree()
    kept = _run(ct, seed=4, record=True)
    assert kept.sequences.events
    assert {e.kind for e in kept.sequences.events} == {"substitution"}
    assert [e.time for e in kept.sequences.events] == sorted(e.time for e in kept.sequences.events)
    assert _run(ct, seed=4).sequences.events == []


# --- the loop does something -----------------------------------------------------------------------

def test_the_gene_s_composition_changes_what_the_genome_keeps():
    """hisA starts AT-rich, where the loss rate is high, and ameliorates toward the model's own even
    composition, where it is low. So how fast the gene turns over decides how much of the genome
    survives. Against the same run with the response flattened to 1.0, more of it does.
    """
    ct = _tree()
    seeds = [s for s in range(1, 9) if _copies(_run(ct, steep=False, seed=s)) > 0]
    assert len(seeds) >= 5, "too few runs kept the family to measure anything"
    loop = statistics.fmean(_copies(_run(ct, seed=s)) for s in seeds)
    flat = statistics.fmean(_copies(_run(ct, steep=False, seed=s)) for s in seeds)
    assert loop > flat * 1.1, f"loop {loop:.1f} against flat {flat:.1f}"


def test_halving_the_step_does_not_move_the_answer():
    ct = _tree()
    coarse = statistics.fmean(_copies(_run(ct, seed=s, step=0.1)) for s in range(1, 9))
    fine = statistics.fmean(_copies(_run(ct, seed=s, step=0.025)) for s in range(1, 9))
    assert abs(coarse - fine) < max(3.0, 0.25 * fine), f"{coarse:.1f} against {fine:.1f}"


def test_transfer_works_here():
    """Unlike the tree-growing joint models, and for the reason they refuse it: on a tree handed to
    the run, the set of lineages alive at an instant is already known."""
    ct = _tree()
    r = _run(ct, transfer=0.15, seed=2)
    assert [e for e in r.genome.events if e.kind.startswith("transfer")]


# --- what is refused ------------------------------------------------------------------------------

def _gene(**kw):
    kw.setdefault("name", "hisA")
    kw.setdefault("model", _EVEN)
    kw.setdefault("length", 100)
    return gene(**kw)


def _spec(**kw):
    kw.setdefault("duplication", 0.1)
    kw.setdefault("initial_families", 3)
    kw.setdefault("families", [family("hisA")])
    return genome_spec(**kw)


def test_it_needs_a_tree():
    with pytest.raises(ValueError, match="pass tree="):
        joint.simulate(
            _spec(loss=PerCopy(0.2).scaled_by("sequences:hisA", Curve(lambda x: 1.0), step=0.05)),
            _gene(offers=composition("GC", absent=0.5)), seed=1)


def test_a_finished_genome_run_is_not_the_handover_here():
    ct = _tree(8)
    from zombi2.genomes import simulate_genomes_family
    g = simulate_genomes_family(ct, initial_families=2, duplication=0.05, loss=0.1, seed=2)
    with pytest.raises(ValueError, match="one of the things being simulated"):
        joint.simulate(
            _spec(loss=PerCopy(0.2).scaled_by("sequences:hisA", Curve(lambda x: 1.0), step=0.05)),
            _gene(offers=composition("GC", absent=0.5)), genomes=g, seed=1)


def test_an_undeclared_gene_is_refused():
    ct = _tree(8)
    with pytest.raises(ValueError, match="names no family"):
        joint.simulate(
            _spec(families=[family("other")],
                  loss=PerCopy(0.2).scaled_by("sequences:hisA", Curve(lambda x: 1.0), step=0.05)),
            _gene(offers=composition("GC", absent=0.5)), tree=ct, seed=1)


def test_the_gene_has_to_say_what_it_offers():
    ct = _tree(8)
    with pytest.raises(ValueError, match="say what it offers"):
        joint.simulate(
            _spec(loss=PerCopy(0.2).scaled_by("sequences:hisA", Curve(lambda x: 1.0), step=0.05)),
            _gene(), tree=ct, seed=1)


def test_neither_reading_the_other_is_two_runs():
    ct = _tree(8)
    with pytest.raises(ValueError, match="two independent runs"):
        joint.simulate(_spec(loss=0.2), _gene(offers=composition("GC", absent=0.5)),
                       tree=ct, seed=1)


def test_a_reading_rate_needs_a_step():
    ct = _tree(8)
    with pytest.raises(ValueError, match="needs a step"):
        joint.simulate(
            _spec(loss=PerCopy(0.2).scaled_by("sequences:hisA", Curve(lambda x: 1.0))),
            _gene(offers=composition("GC", absent=0.5)), tree=ct, seed=1)


def test_a_genome_rate_reading_something_else_is_refused():
    ct = _tree(8)
    with pytest.raises(ValueError, match="reads the gene in this same run"):
        joint.simulate(
            _spec(loss=PerCopy(0.2).scaled_by("trait", {"a": 1.0}, step=0.05)),
            _gene(offers=composition("GC", absent=0.5)), tree=ct, seed=1)


def test_letters_outside_the_alphabet_are_refused():
    ct = _tree(8)
    with pytest.raises(ValueError, match="not in this model's alphabet"):
        joint.simulate(
            _spec(loss=PerCopy(0.2).scaled_by("sequences:hisA", Curve(lambda x: 1.0), step=0.05)),
            _gene(offers=composition("KR", absent=0.5)), tree=ct, seed=1)
