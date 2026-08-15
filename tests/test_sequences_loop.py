"""The sequence level joined to itself — two genes, each one's rate reading the other's composition.

Step 9 of the joining design note. The walk is by **time** rather than by family: species time is cut
into slices, and inside a slice every living gene copy of every family advances at a rate read off
the other family's composition as it stood at the top of the slice.
"""

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, PerSite
from zombi2.sequences import composition, gene, lg, simulate_sequences
from zombi2.sequences.substitution_models import AMINO_ACIDS, reversible

_KR = [AMINO_ACIDS.index(c) for c in "KR"]


def _poor(scale=0.12):
    """LG's chemistry over KR-depleted frequencies — a founding composition far from equilibrium, so
    the genes have somewhere to go and their compositions mean something while they go there."""
    m = lg()
    S = m.Q / m.stationary[None, :]
    S = (S + S.T) / 2.0
    np.fill_diagonal(S, 0.0)
    pi = m.stationary.copy()
    pi[_KR] *= scale
    return reversible(S, pi / pi.sum(), name="KR-poor LG", alphabet=AMINO_ACIDS)


def _genomes(n=25, loss=0.0, names=("A", "B"), seed=2):
    ct = species.simulate_species_tree(birth=1.0, n_extant=n, seed=1).complete_tree
    return ct, simulate_genomes_family(ct, initial_families=4, duplication=0.0, loss=loss,
                                       origination=0.0, seed=seed,
                                       families=[family(x) for x in names])


def _pair(g, *, curve=lambda x: 0.10 + 25.0 * x, seed=1, step=0.05, length=250, absent=0.02):
    start = _poor()
    f = Curve(curve)
    return simulate_sequences(g, joint=True, seed=seed, genes=[
        gene(name=n, model=lg(), length=length, start=start,
             offers=composition("KR", absent=absent),
             substitution=PerSite(0.5).scaled_by(f"sequences:{other}", f, step=step))
        for n, other in (("A", "B"), ("B", "A"))])


def _tip_share(result, fam):
    hits = total = 0
    for seq in result.alignments[fam].values():
        hits += sum(seq.count(c) for c in "KR")
        total += len(seq)
    return hits / total


def _both(result):
    a, b = sorted(result.alignments)
    return (_tip_share(result, a) + _tip_share(result, b)) / 2


# --- what comes back ------------------------------------------------------------------------------

def test_one_result_holding_both_genes():
    _ct, g = _genomes()
    r = _pair(g)
    assert sorted(r.alignments) == sorted(g.family_names[n] for n in ("A", "B"))
    assert r.families == ("A", "B")
    assert set(r.phylograms) == set(r.alignments)
    for fam in r.alignments:
        assert r.alignments[fam] and r.ancestral[fam]
        assert len(r.founding[fam]) == 250


def test_the_phylograms_are_the_trees_the_alignments_were_drawn_along():
    """A driven rate is not constant down a branch here — it changes at every slice boundary — so
    each branch length is the walk's own accumulation rather than one sample of the rate."""
    _ct, g = _genomes()
    r = _pair(g)
    for fam in r.phylograms:
        nwk = r.phylograms[fam]["complete"]
        lengths = [float(x) for x in __import__("re").findall(r":([0-9.eE+-]+)", nwk)]
        assert lengths and all(v >= 0.0 for v in lengths)
        assert sum(lengths) > 0.0


def test_it_is_deterministic():
    _ct, g = _genomes()
    assert _pair(g, seed=7).alignments == _pair(g, seed=7).alignments


def test_a_joint_run_has_no_species_phylogram_and_says_so():
    """It is the clock made visible, and here every gene runs at a rate the other gene sets — so no
    one set of branch lengths belongs to the run rather than to a gene."""
    import pathlib
    import tempfile
    _ct, g = _genomes()
    r = _pair(g)
    assert r.species_phylogram == {"complete": None, "extant": None}
    with tempfile.TemporaryDirectory() as d:
        r.write(d)
        names = {p.name for p in pathlib.Path(d).rglob("*") if p.is_file()}
    assert not any(n.startswith("clock_species_tree") for n in names)
    assert any(n.startswith("alignment") or "fam" in n for n in names)


# --- the loop does something ----------------------------------------------------------------------

def test_the_loop_changes_how_far_the_pair_gets():
    """Both genes are founded far from LG's equilibrium and ameliorate toward it. Each one's speed is
    read off the other's composition, and the response rises with it — so the pair is slow while both
    are still KR-poor and accelerates as they are not. Against the same run with the response
    flattened to 1.0, that carries them measurably further.
    """
    _ct, g = _genomes()
    loop = [_both(_pair(g, seed=s)) for s in range(1, 5)]
    flat = [_both(_pair(g, seed=s, curve=lambda x: 1.0)) for s in range(1, 5)]
    assert min(loop) > max(flat), f"loop {['%.3f' % v for v in loop]} against flat " \
                                  f"{['%.3f' % v for v in flat]}"


def test_halving_the_step_does_not_move_the_answer():
    """The check the manual asks for. Slicing is the approximation here, so a run whose answer moves
    when `step` is halved was read at too coarse a resolution."""
    _ct, g = _genomes()
    coarse = [_both(_pair(g, seed=s, step=0.10)) for s in range(1, 5)]
    fine = [_both(_pair(g, seed=s, step=0.025)) for s in range(1, 5)]
    assert abs(sum(coarse) / 4 - sum(fine) / 4) < 0.012


def test_a_gene_a_lineage_does_not_carry_reads_the_declared_absent():
    """The gate from step 8, on the live side: a lineage with no copy of the family has no sequence
    to count, so the driver reads what the gene declared."""
    _ct, g = _genomes(n=25, loss=0.35)
    r = _pair(g, absent=0.30)
    assert sum(len(a) for a in r.alignments.values()) > 0


# --- what is refused ------------------------------------------------------------------------------

def test_joint_needs_the_genes():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="needs the genes"):
        simulate_sequences(g, joint=True, seed=1)


def test_one_gene_is_an_ordinary_run():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="One gene is"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100,
                 substitution=PerSite(1.0).scaled_by("sequences:A", Curve(lambda x: 1.0),
                                                     step=0.05))])


def test_reading_another_gene_needs_joint_true():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="joint=True"):
        simulate_sequences(g, seed=1, genes=[
            gene(name="A", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100, offers=composition("KR", absent=0.02))])


def test_joint_true_needs_something_to_read():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="none reads another"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100),
            gene(name="B", model=lg(), length=100)])


def test_a_gene_cannot_read_itself():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="its own composition"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:A", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100)])


def test_reading_a_gene_that_offers_nothing_is_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="offers nothing"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100,
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100)])


def test_a_reading_rate_needs_a_step():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="needs a step"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0))),
            gene(name="B", model=lg(), length=100, offers=composition("KR", absent=0.02))])


def test_two_resolutions_are_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="have to agree"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:A", Curve(lambda x: 1.0),
                                                     step=0.01))])


def test_a_step_as_long_as_the_tree_is_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="not shorter than the tree"):
        _pair(g, step=50.0)


def test_an_undeclared_gene_name_is_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="name no family"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="nope", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100, offers=composition("KR", absent=0.02))])


def test_the_run_wide_arguments_have_nothing_to_apply_to():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="nothing to apply to"):
        simulate_sequences(g, joint=True, seed=1, model=lg(), length=100, genes=[
            gene(name="A", model=lg(), length=100, offers=composition("KR", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100, offers=composition("KR", absent=0.02))])


def test_an_offered_statistic_needs_an_absent():
    with pytest.raises(TypeError):
        composition("KR")                       # absent is keyword-only and required


def test_letters_outside_the_alphabet_are_refused():
    _ct, g = _genomes()
    with pytest.raises(ValueError, match="not in this run's alphabet"):
        simulate_sequences(g, joint=True, seed=1, genes=[
            gene(name="A", model=lg(), length=100, offers=composition("BJZ", absent=0.02),
                 substitution=PerSite(1.0).scaled_by("sequences:B", Curve(lambda x: 1.0),
                                                     step=0.05)),
            gene(name="B", model=lg(), length=100, offers=composition("KR", absent=0.02))])
