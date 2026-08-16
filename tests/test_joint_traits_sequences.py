"""A trait and a gene's sequence, each driving the other, on a tree the run is handed.

Step 10 of the joining design note — the cross-level join whose two ends are furthest apart. Both
directions already ran as conditioning; written at once they are a cycle, and the run walks species
time in slices to carry both.
"""

import statistics

import numpy as np
import pytest

from zombi2 import joint, species, traits
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, PerLineage, PerSite
from zombi2.sequences import composition, gene, lg
from zombi2.sequences.substitution_models import AMINO_ACIDS, reversible


def _poor(scale=0.12):
    """LG's chemistry over KR-depleted frequencies — where the gene starts, so its composition has
    somewhere to go and the trait rate reading it reads something that moves."""
    m = lg()
    S = m.Q / m.stationary[None, :]
    S = (S + S.T) / 2.0
    np.fill_diagonal(S, 0.0)
    pi = m.stationary.copy()
    pi[[AMINO_ACIDS.index(c) for c in "KR"]] *= scale
    return reversible(S, pi / pi.sum(), name="KR-poor LG", alphabet=AMINO_ACIDS)


def _genomes(n=30, loss=0.0, seed=2):
    ct = species.simulate_species_tree(birth=1.0, n_extant=n, seed=1).complete_tree
    return ct, simulate_genomes_family(ct, initial_families=3, duplication=0.0, loss=loss,
                                       origination=0.0, seed=seed, families=[family("rpoB")])


def _run(g, *, curve=lambda x: 0.05 + 30.0 * x, seed=1, step=0.05, absent=0.02, **trait_kw):
    """A KR-rich rpoB turns a lineage hot, and a hot lineage substitutes four times faster — which is
    what makes rpoB KR-rich. Neither can be finished first."""
    return joint.simulate(
        traits.discrete(name="habitat", states=["cold", "hot"], start="cold",
                        switch={"cold->hot": PerLineage(0.5).scaled_by(
                                    "sequences:rpoB", Curve(curve), step=step),
                                "hot->cold": 0.3},
                        **trait_kw),
        gene(name="rpoB", model=lg(), length=250, start=_poor(),
             offers=composition("KR", absent=absent),
             substitution=PerSite(0.6).scaled_by("trait", {"hot": 4.0, "cold": 1.0})),
        genomes=g, seed=seed)


def _hot_share(r):
    tips = sorted(r.complete_tree.extant_leaves())
    return sum(1 for t in tips if r.trait.node_values[t] == "hot") / len(tips)


# --- what comes back ------------------------------------------------------------------------------

def test_it_returns_both_levels():
    ct, g = _genomes()
    r = _run(g)
    assert r.trait is not None and r.sequences is not None and r.genome is None
    assert set(r.trait.node_values) == set(ct.nodes)
    fam = g.family_names["rpoB"]
    assert sorted(r.sequences.alignments) == [fam]
    assert r.sequences.families == ("rpoB",)
    assert r.complete_tree is ct


def test_the_tree_came_in_so_its_log_is_not_this_run_s():
    ct, g = _genomes()
    r = _run(g)
    assert r.species.events == []
    assert r.n_extant == sum(1 for n in ct.nodes.values() if n.fate == "extant")


def test_the_trait_log_opens_with_the_initial_row():
    ct, g = _genomes()
    r = _run(g)
    first = r.trait.events[0]
    assert (first.kind, first.lineage, first.from_state) == ("initial", ct.root, None)
    assert [c.time for c in r.trait.events] == sorted(c.time for c in r.trait.events)


def test_it_is_deterministic():
    _ct, g = _genomes()
    a, b = _run(g, seed=5), _run(g, seed=5)
    assert a.sequences.alignments == b.sequences.alignments
    assert a.trait.node_values == b.trait.node_values


def test_both_levels_write():
    import pathlib
    import tempfile
    _ct, g = _genomes()
    with tempfile.TemporaryDirectory() as d:
        _run(g, seed=3).write(d)
        names = {p.name for p in pathlib.Path(d).rglob("*") if p.is_file()}
    assert "trait_events.tsv" in names and "joint_summary.json" in names
    assert not any(n.startswith("clock_species_tree") for n in names)


def test_a_jump_at_the_split_still_works():
    _ct, g = _genomes()
    r = _run(g, at_speciation=0.4, seed=2)
    assert [c for c in r.trait.events if c.kind == "on_speciation"]


# --- the cycle does something ---------------------------------------------------------------------

def test_the_cycle_pushes_the_trait():
    """Both arrows point the same way round the loop: a KR-rich gene turns a lineage hot, and a hot
    lineage substitutes faster, which carries the gene further toward LG's own KR share. Against the
    same run with the trait's response flattened to 1.0, that leaves far more of the tree hot.
    """
    _ct, g = _genomes()
    loop = statistics.fmean(_hot_share(_run(g, seed=s)) for s in range(1, 9))
    flat = statistics.fmean(_hot_share(_run(g, seed=s, curve=lambda x: 1.0)) for s in range(1, 9))
    assert loop > flat + 0.10, f"loop {loop:.3f} against flat {flat:.3f}"


def test_halving_the_step_does_not_move_the_answer():
    _ct, g = _genomes()
    coarse = statistics.fmean(_hot_share(_run(g, seed=s, step=0.10)) for s in range(1, 9))
    fine = statistics.fmean(_hot_share(_run(g, seed=s, step=0.025)) for s in range(1, 9))
    assert abs(coarse - fine) < 0.12, f"step=0.10 gave {coarse:.3f}, step=0.025 gave {fine:.3f}"


def test_a_lineage_without_the_gene_reads_the_declared_absent():
    _ct, g = _genomes(loss=0.4)
    r = _run(g, absent=0.30)
    assert r.trait.events and r.sequences.alignments


# --- what is refused ------------------------------------------------------------------------------

def _gene(**kw):
    kw.setdefault("name", "rpoB")
    kw.setdefault("model", lg())
    kw.setdefault("length", 100)
    return gene(**kw)


def test_it_needs_the_genome_run_rather_than_a_bare_tree():
    ct, g = _genomes()
    with pytest.raises(ValueError, match="the genome run that"):
        joint.simulate(traits.discrete(name="h", states=["a", "b"], switch=0.1),
                       _gene(offers=composition("KR", absent=0.02)), tree=ct, seed=1)


def test_passing_both_a_tree_and_a_genome_run_is_refused():
    ct, g = _genomes()
    with pytest.raises(ValueError, match="two answers to one question"):
        joint.simulate(
            traits.discrete(name="h", states=["a", "b"],
                            switch=PerLineage(0.1).scaled_by("sequences:rpoB",
                                                             Curve(lambda x: 1.0), step=0.05)),
            _gene(offers=composition("KR", absent=0.02)), tree=ct, genomes=g, seed=1)


def test_a_species_spec_beside_them_is_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="inputs rather than participants"):
        joint.simulate(species.birth_death(birth=1.0, n_extant=10),
                       traits.discrete(name="h", states=["a", "b"], switch=0.1),
                       _gene(offers=composition("KR", absent=0.02)), genomes=g, seed=1)


def test_a_gene_the_genome_run_never_declared_is_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="names no family"):
        joint.simulate(
            traits.discrete(name="h", states=["a", "b"],
                            switch=PerLineage(0.1).scaled_by("sequences:nope",
                                                             Curve(lambda x: 1.0), step=0.05)),
            _gene(name="nope", offers=composition("KR", absent=0.02)), genomes=g, seed=1)


def test_the_gene_has_to_say_what_it_publishes():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="say what it publishes"):
        joint.simulate(
            traits.discrete(name="h", states=["a", "b"],
                            switch=PerLineage(0.1).scaled_by("sequences:rpoB",
                                                             Curve(lambda x: 1.0), step=0.05)),
            _gene(), genomes=g, seed=1)


def test_neither_reading_the_other_is_two_runs():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="two independent runs"):
        joint.simulate(traits.discrete(name="h", states=["a", "b"], switch=0.1),
                       _gene(offers=composition("KR", absent=0.02)), genomes=g, seed=1)


def test_the_trait_reading_a_composition_needs_a_step():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="needs a step"):
        joint.simulate(
            traits.discrete(name="h", states=["a", "b"],
                            switch=PerLineage(0.1).scaled_by("sequences:rpoB",
                                                             Curve(lambda x: 1.0))),
            _gene(offers=composition("KR", absent=0.02)), genomes=g, seed=1)


def test_a_continuous_trait_here_says_it_is_not_built():
    _ct, g = _genomes()
    with pytest.raises(NotImplementedError, match="CONTINUOUS trait"):
        joint.simulate(traits.continuous(name="size", rate=1.0),
                       _gene(offers=composition("KR", absent=0.02)), genomes=g, seed=1)


def test_a_genome_run_with_no_sequence_participant_is_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="only a sequence participant"):
        joint.simulate(species.birth_death(birth=1.0, n_extant=10),
                       traits.discrete(name="h", states=["a", "b"],
                                       switch=0.1), genomes=g, seed=1)
