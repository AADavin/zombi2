"""One family's sequences, and one family's composition as a driver (design note §9).

Two pieces that only work together. `families=` restricts a sequence run to named families, so its
pooled `composition` is that family's; `absent=` says what a branch reads where the family is not
there, which is what a pooled statistic never had to answer. This is **conditioning** — the driver
run finishes before the driven one starts — and it is here because the sequence loop needs it.
"""

import re

import pytest

from zombi2 import species
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, PerSite
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, lg


def _run(n=25, loss=0.25, seed=2, names=("chaperone", "client")):
    ct = species.simulate_species_tree(birth=1.0, n_extant=n, seed=1).complete_tree
    g = simulate_genomes_family(ct, initial_families=8, duplication=0.05, loss=loss,
                                origination=0.1, seed=seed,
                                families=[family(x) for x in names])
    return ct, g


def _phylogram_length(result, fam):
    return sum(float(x) for x in re.findall(r":([0-9.eE+-]+)", result.phylograms[fam]["complete"]))


# --- restricting a run ----------------------------------------------------------------------------

def test_families_evolves_only_the_named_ones():
    _ct, g = _run()
    s = simulate_sequences(g, families=["chaperone"], model=lg(), length=120, seed=3)
    assert sorted(s.alignments) == [g.family_names["chaperone"]]
    assert s.families == ("chaperone",)


def test_leaving_it_off_still_evolves_everything():
    _ct, g = _run()
    s = simulate_sequences(g, model=lg(), length=120, seed=3)
    assert len(s.alignments) == len(g.gene_trees) > 1
    assert s.families == ()


def test_a_restricted_run_is_the_same_family_the_whole_run_would_have_given():
    """Restricting picks families out of the run rather than changing how one evolves: the seed is
    the run's, so the family's own draws are whatever position it holds in the loop — this checks
    the tree it was evolved along, which is the genome's and cannot move."""
    _ct, g = _run()
    whole = simulate_sequences(g, model=lg(), length=120, seed=3)
    one = simulate_sequences(g, families=["client"], model=lg(), length=120, seed=3)
    fam = g.family_names["client"]
    assert one.phylograms[fam]["complete"].count(",") == whole.phylograms[fam]["complete"].count(",")


def test_it_is_deterministic():
    _ct, g = _run()
    a = simulate_sequences(g, families=["client"], model=lg(), length=120, seed=9)
    b = simulate_sequences(g, families=["client"], model=lg(), length=120, seed=9)
    assert a.alignments == b.alignments


# --- one family's composition ---------------------------------------------------------------------

def test_absent_is_what_a_branch_without_the_family_reads():
    ct, g = _run()
    s = simulate_sequences(g, families=["chaperone"], model=lg(), length=300, seed=3)
    values = s.composition("KR", absent=0.08)._node_values(ct)
    assert set(values) == set(ct.nodes)
    assert 0.08 in set(values.values()), "no branch was missing the family, so nothing is tested"
    assert all(0.0 <= v <= 1.0 for v in values.values())


def test_a_restricted_run_with_a_gap_and_no_absent_is_refused():
    """SPEC §5, on a driver: carrying the parent's value forward would drive the branch as though the
    family were still there, which is a different model from the one asked for."""
    ct, g = _run()
    s = simulate_sequences(g, families=["chaperone"], model=lg(), length=300, seed=3)
    with pytest.raises(ValueError, match="carries none of it"):
        s.composition("KR")._node_values(ct)


def test_a_family_present_everywhere_needs_no_absent():
    """The gate is about a gap, not about `families=`: a family that never goes missing has a value
    on every branch already."""
    ct, g = _run(n=12, loss=0.0, names=("core",))
    s = simulate_sequences(g, families=["core"], model=lg(), length=200, seed=3)
    values = s.composition("KR")._node_values(ct)
    assert set(values) == set(ct.nodes)


def test_a_pooled_run_still_carries_the_parent_value_forward():
    """Unchanged, and deliberately so: a gap in a whole-genome run is a lineage that momentarily
    held nothing, not a statistic that stopped existing."""
    ct, g = _run()
    s = simulate_sequences(g, model=lg(), length=200, seed=3)
    values = s.composition("KR")._node_values(ct)
    assert set(values) == set(ct.nodes)


def test_the_composition_says_what_it_is():
    _ct, g = _run()
    s = simulate_sequences(g, families=["chaperone"], model=lg(), length=120, seed=3)
    assert repr(s.composition("KR", absent=0.08)) == "composition('KR', absent=0.08)"
    assert repr(s.composition("KR")) == "composition('KR')"


# --- one gene's sequence driving another's rate ---------------------------------------------------

def test_one_family_drives_another_families_substitution_rate():
    """The pair `families=` exists for. Both runs read the same genome run, so the driver's history
    is the one on disk — without the restriction the second run would re-evolve the chaperone beside
    its target and read a third history that nothing records.
    """
    ct, g = _run()
    chaperone = simulate_sequences(g, families=["chaperone"], model=lg(), length=300, seed=3)
    basic = chaperone.composition("KR", absent=0.08)
    spread = chaperone.composition("KR", absent=0.08)._node_values(ct)
    assert max(spread.values()) - min(spread.values()) > 0.02, "the driver barely varies"

    def client(curve):
        return _phylogram_length(
            simulate_sequences(g, families=["client"], model=lg(), length=300, seed=4,
                               substitution=PerSite(1.0).scaled_by(basic, Curve(curve))),
            g.family_names["client"])

    flat = client(lambda x: 1.0)
    steep = client(lambda x: 20.0 ** (x - 0.1))
    assert steep > flat * 1.03, f"the driver left no mark: flat {flat:.2f}, steep {steep:.2f}"


# --- what is refused ------------------------------------------------------------------------------

def test_an_undeclared_name_is_refused():
    _ct, g = _run()
    with pytest.raises(ValueError, match="names no family"):
        simulate_sequences(g, families=["nope"], model=lg(), length=120, seed=3)


def test_one_name_is_a_list_not_a_string():
    _ct, g = _run()
    with pytest.raises(TypeError, match="list of names"):
        simulate_sequences(g, families="chaperone", model=lg(), length=120, seed=3)


def test_an_empty_list_is_refused():
    _ct, g = _run()
    with pytest.raises(ValueError, match="would evolve nothing"):
        simulate_sequences(g, families=[], model=lg(), length=120, seed=3)


def test_absent_must_be_a_fraction():
    _ct, g = _run()
    s = simulate_sequences(g, families=["chaperone"], model=lg(), length=120, seed=3)
    with pytest.raises(ValueError, match=r"fraction in \[0, 1\]"):
        s.composition("KR", absent=1.4)


def test_gc_still_takes_no_family_and_now_says_where_to_go():
    ct = species.simulate_species_tree(birth=1.0, n_extant=8, seed=1).complete_tree
    g = simulate_genomes_family(ct, initial_families=4, duplication=0.05, loss=0.1, seed=2,
                                families=[family("rrn")])
    s = simulate_sequences(g, model=hky85(2.0), length=120, seed=3)
    with pytest.raises(ValueError, match="restrict the run to it"):
        s.gc("rrn")
    assert repr(s.gc(absent=0.5)) == "composition('GC', absent=0.5)"


def test_a_nucleotide_run_has_blocks_rather_than_families():
    from zombi2.genomes import simulate_genomes_nucleotide
    ct = species.simulate_species_tree(birth=1.0, n_extant=6, seed=1).complete_tree
    g = simulate_genomes_nucleotide(ct, genes=3, gene_length=90, root_length=600,
                                    duplication=0.5, duplication_extent=40, loss=0.5,
                                    loss_extent=40, seed=2)
    with pytest.raises(ValueError, match="evolves blocks"):
        simulate_sequences(g, families=["anything"], model=hky85(2.0), seed=3)
