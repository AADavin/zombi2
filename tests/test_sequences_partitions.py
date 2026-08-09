"""Sequences: **partitions** — a family's sites split into blocks, each under its own model.

``partitions=[(hky85(kappa=2.0), 600), (jc69(), 400)]`` replaces ``model=`` and ``length=``: the
blocks are evolved down the *same* gene tree at the *same* rate and concatenated, in order, into one
sequence per gene copy. Three things are worth testing here and they are not the obvious ones.

The first is that the split is real — that the second block is not quietly evolving under the first
block's matrix, which nothing about the output would reveal. So a long two-partition run under models
with very different base composition is checked half by half against each model's own stationary
frequencies.

The second is that **one phylogram is still exact**. Every model on the menu is normalised to one
expected substitution per site per unit branch length and every set of across-site rate classes to a
mean of 1, so a branch of ``Δt`` accrues ``rate · Δt`` substitutions per site in each partition alike.
That is what lets the family keep a single tree, and it is checked by measuring the realised
divergence in each half of a run under two genuinely different matrices.

The third is that nothing moved. A single-partition run must be the same bytes as the plain
``model=``/``length=`` run it is spelled differently from — that is the regression test on the default
path, which everything else in the suite pins from the other side.
"""

from __future__ import annotations

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import FamilyGenomesResult, simulate_genomes_family, simulate_genomes_nucleotide
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.rates import LogNormal, modifiers as mod
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, jc69, lg
from zombi2.species import simulate_species_tree


def _genome_run(gene_trees, *, t_split: float = 1.0, t_now: float = 2.0) -> FamilyGenomesResult:
    """The genome run the sequence level requires, around hand-built gene trees: a three-lineage
    species tree (root 0 splits at ``t_split`` into extant tips 1 and 2 at ``t_now``) carrying the
    given families. The gene trees are attached directly — these families are written by hand, so
    there is no event log for the run to derive them from."""
    tree = species.Tree({0: species.Node(0, None, 0.0, t_split, (1, 2), "speciation"),
                         1: species.Node(1, 0, t_split, t_now, None, "extant"),
                         2: species.Node(2, 0, t_split, t_now, None, "extant")}, 0)
    run = FamilyGenomesResult(complete_tree=tree, genomes={}, edges=[], seed=None)
    run.gene_trees = dict(gene_trees)      # a cached_property: the instance dict wins
    return run


def _pair_run(t_spec: float, t_tip: float) -> FamilyGenomesResult:
    """One family: the founding gene (id 0) on species 0 speciates at ``t_spec``, its two daughters
    reach extant tips at ``t_tip``. Two tips, one branch each — so a p-distance between them is a
    measurement of ``2 · rate · (t_tip - t_spec)`` and nothing else."""
    root = GeneNode("speciation", 0, t_spec, 0)
    root.children = [GeneNode("extant", 1, t_tip, 1), GeneNode("extant", 2, t_tip, 2)]
    return _genome_run({0: GeneTree(0, root, 0.0)}, t_split=t_spec, t_now=t_tip)


def _tips(result, fam: int = 0) -> list[str]:
    return list(result.alignments[fam].values())


def _composition(seq: str) -> np.ndarray:
    return np.array([seq.count(c) for c in "ACGT"], dtype=float) / len(seq)


def _p_distance(a: str, b: str) -> float:
    return sum(1 for x, y in zip(a, b) if x != y) / len(a)


def _family_run(seed: int = 2):
    """A real (small) genome run, for the tests that need many families rather than one hand-built
    tree — streaming, and the coverage of every output."""
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1)
    return simulate_genomes_family(sp, duplication=0.3, transfer=0.1, loss=0.2, origination=0.1,
                                   initial_families=8, seed=seed)


# --- the default path did not move -----------------------------------------------------------------

def test_one_partition_is_bit_identical_to_the_plain_model():
    """The regression test. A single partition is the un-partitioned run spelled differently, so it
    must consume the same randomness and produce the same bytes — alignments, ancestral sequences,
    founding sequences and phylograms alike."""
    run = _family_run()
    split = simulate_sequences(run, partitions=[(jc69(), 200)], seed=7)
    plain = simulate_sequences(run, model=jc69(), length=200, seed=7)
    assert split.alignments == plain.alignments
    assert split.ancestral == plain.ancestral
    assert split.founding == plain.founding
    assert split.phylograms == plain.phylograms
    assert split.species_phylogram == plain.species_phylogram


# --- the split is real -----------------------------------------------------------------------------

def test_each_partition_evolves_under_its_own_model():
    """20 000 sites of JC69 followed by 20 000 of an HKY85 with strongly skewed frequencies, run far
    enough to saturate. Each half must sit at *its own* model's stationary composition: if the second
    block were evolving under the first block's matrix, both halves would come out uniform."""
    frequencies = (0.1, 0.1, 0.4, 0.4)
    n = 20000
    r = simulate_sequences(_pair_run(0.0, 5.0), substitution=2.0, seed=4,
                           partitions=[(jc69(), n), (hky85(kappa=1.0, frequencies=frequencies), n)])
    seq = _tips(r)[0]
    assert np.allclose(_composition(seq[:n]), 0.25, atol=0.02)
    assert np.allclose(_composition(seq[n:]), frequencies, atol=0.02)


def test_the_realised_divergence_is_the_same_in_every_partition():
    """The phylogram guarantee, measured. Two genuinely different matrices (JC69 and a strongly
    transition-biased HKY85) at the same rate must accrue the same *number* of substitutions per
    site, because both are normalised to one per unit branch length. So the two halves of the
    alignment have the same p-distance, and it is the one JC69 theory predicts for the branch
    lengths in the phylogram: ``p = 3/4 · (1 - exp(-4d/3))`` over ``d`` = 0.1 substitutions/site
    between the tips."""
    n = 20000
    r = simulate_sequences(_pair_run(0.0, 1.0), divergence=0.05, seed=9,
                           partitions=[(jc69(), n), (hky85(kappa=5.0), n)])
    a, b = _tips(r)
    first, second = _p_distance(a[:n], b[:n]), _p_distance(a[n:], b[n:])
    expected = 0.75 * (1 - np.exp(-4 * 0.1 / 3))               # two branches of 0.05 subs/site each
    assert abs(first - second) < 0.01
    assert abs(first - expected) < 0.01


def test_the_phylogram_does_not_depend_on_the_partitioning():
    """Branch lengths are ``base × clock × Δt`` and know nothing about which matrix ran where, so a
    partitioned run and a single-model run at the same seed carry the *same* trees."""
    run = _family_run()
    split = simulate_sequences(run, partitions=[(hky85(kappa=2.0), 60), (jc69(), 40)], seed=5)
    plain = simulate_sequences(run, model=jc69(), length=100, seed=5)
    assert split.phylograms == plain.phylograms
    assert split.species_phylogram == plain.species_phylogram


def test_a_sequence_is_as_long_as_the_partitions_sum_to():
    parts = [(jc69(), 30), (hky85(kappa=2.0), 45), (jc69(), 25)]
    r = simulate_sequences(_family_run(), partitions=parts, seed=1)
    total = sum(n for _, n in parts)
    for fam, aln in r.alignments.items():
        assert all(len(s) == total for s in aln.values())
        assert all(len(s) == total for s in r.ancestral[fam].values())
        assert len(r.founding[fam]) == total


def test_a_partition_carries_its_own_across_site_variation():
    """The two axes compose: rate classes belong to a model (SPEC §5 — they are not a rate modifier),
    so each partition brings its own. A partition under ``+I`` at 90% invariant changes far less than
    a flat one beside it, at the same rate, in the same run."""
    n = 4000
    flat, frozen = jc69(), jc69().across_sites(invariant=0.9)
    r = simulate_sequences(_pair_run(0.0, 1.0), substitution=0.5, seed=6,
                           partitions=[(flat, n), (frozen, n)])
    a, b = _tips(r)
    # the +I half concentrates all its change in a tenth of the columns, so most columns match
    assert _p_distance(a[n:], b[n:]) < 0.5 * _p_distance(a[:n], b[:n])


def test_partitions_compose_with_a_lineage_clock():
    """A clock is a per-lineage factor on the one rate every partition shares, so it reshapes the
    whole sequence rather than any part of it — and the run still works, which is the point."""
    run = _family_run()
    parts = [(hky85(kappa=2.0), 50), (jc69(), 50)]
    strict = simulate_sequences(run, partitions=parts, seed=3)
    relaxed = simulate_sequences(run, partitions=parts, seed=3,
                                 substitution=1.0 * mod.Drawn(per='lineage', dist=LogNormal(0.0, 0.6)))
    assert set(strict.alignments) == set(relaxed.alignments)
    assert strict.alignments != relaxed.alignments
    assert strict.phylograms != relaxed.phylograms


# --- what is refused -------------------------------------------------------------------------------

def test_partitions_and_model_or_length_cannot_both_be_given():
    run = _pair_run(1.0, 2.0)
    with pytest.raises(ValueError, match="second answer to the same question"):
        simulate_sequences(run, model=jc69(), partitions=[(jc69(), 10)], seed=1)
    with pytest.raises(ValueError, match="the partitions set the length"):
        simulate_sequences(run, partitions=[(jc69(), 10)], length=10, seed=1)
    with pytest.raises(ValueError, match="no model"):
        simulate_sequences(run, seed=1)                        # neither given


def test_partitions_must_share_one_alphabet():
    """Meaningless rather than unimplemented, and the message has to say so: the partitions are
    concatenated into one sequence, and there is no string that is half DNA and half protein."""
    with pytest.raises(ValueError, match="one alphabet") as e:
        simulate_sequences(_pair_run(1.0, 2.0), partitions=[(jc69(), 10), (lg(), 10)], seed=1)
    assert "JC69" in str(e.value) and "LG" in str(e.value)     # both partitions are named
    assert "meaningless" in str(e.value)


@pytest.mark.parametrize("parts, match", [
    ([(jc69(), 0)], "positive whole number"),
    ([(jc69(), 1.5)], "positive whole number"),
    ([(jc69(), 10), (jc69(), -3)], "partition 1"),             # the offending index, not the first
    ([("jc69", 10)], "not a SubstitutionModel"),
    ([jc69()], "on its own"),                                  # a bare model, no site count
    ([(jc69(), 10, 2)], "is a \\(model, sites\\) pair"),
    ([], "would evolve no sites"),
])
def test_a_partition_needs_a_model_and_a_positive_site_count(parts, match):
    with pytest.raises(ValueError, match=match):
        simulate_sequences(_pair_run(1.0, 2.0), partitions=parts, seed=1)


def test_a_nucleotide_run_refuses_partitions():
    """A nucleotide genome has already said which stretch takes which model — genes under ``model``,
    spacer under ``intergene_model`` — and each block carries its own length in bp. Partitions would
    be a second, contradicting answer."""
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=3)
    gen = simulate_genomes_nucleotide(sp, loss=0.5, loss_extent=40, duplication=0.5,
                                      duplication_extent=40, root_length=400, genes=2,
                                      gene_length=90, seed=3)
    with pytest.raises(ValueError, match="the genome sets both the lengths and the split"):
        simulate_sequences(gen, partitions=[(jc69(), 100)], seed=1)


# --- streaming and the parallel engine -------------------------------------------------------------

def test_a_streamed_partitioned_run_is_the_same_dataset(tmp_path):
    """Partitions are concatenated before anything downstream sees them, so streaming composes for
    free — but "for free" is the kind of claim that is worth checking once."""
    run = _family_run()
    parts = [(hky85(kappa=2.0), 60), (jc69(), 40)]
    streamed = simulate_sequences(run, partitions=parts, seed=11, stream_to=tmp_path / "streamed")
    memory = simulate_sequences(run, partitions=parts, seed=11)
    memory.write(tmp_path / "memory", outputs=streamed.outputs)
    for path in sorted((tmp_path / "streamed" / "alignments").glob("*.fasta")):
        twin = tmp_path / "memory" / "alignments" / path.name
        assert path.read_text(encoding="utf-8") == twin.read_text(encoding="utf-8")
    assert streamed.n_sequences == sum(len(a) for a in memory.alignments.values() if a)
