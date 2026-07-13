"""Regression tests: a structural LOSS must honour the min-genome floor (2026-07 audit).

Before the fix, ``_draw_length`` returned the whole chromosome at ``extension >= 1.0`` and the LOSS
path never clamped to ``MIN_GENOME_LENGTH`` (only intergenic deletions did), so the first loss
emptied the chromosome — ``size() -> 0`` zeroes the size-proportional event rate and freezes the
lineage. The floor now applies to LOSS too, matching ``_apply_deletion``.
"""

import numpy as np

from zombi2.genomes import simulate_nucleotide_genomes
from zombi2.genomes.nucleotide_genome import MIN_GENOME_LENGTH
from zombi2.tree import read_newick


def test_loss_never_empties_the_genome_even_at_extension_one():
    tree = read_newick("((a:1,b:1):1,(c:1,d:1):1)R:0;")
    for seed in range(25):
        res = simulate_nucleotide_genomes(tree, loss=0.4, extension=1.0, root_length=15, seed=seed)
        for genome in res.leaf_genomes.values():
            assert genome.total_length() >= MIN_GENOME_LENGTH, "a LOSS emptied the genome"


def test_loss_still_shrinks_genomes():
    # the floor must not disable loss — genomes should still lose material below root_length
    tree = read_newick("(A:1,B:1)root:0;")
    shrank = False
    for seed in range(20):
        res = simulate_nucleotide_genomes(tree, loss=0.6, extension=0.7, root_length=40, seed=seed)
        if any(g.total_length() < 40 for g in res.leaf_genomes.values()):
            shrank = True
            break
    assert shrank, "loss never shrank any genome — the floor is over-clamping"


def test_originate_does_not_reseed_a_full_chromosome():
    # a later origination inserts one novel (short) gene, it does not lay down another root_length
    # chromosome (the seed happens once, gated on the _seeded flag)
    import types
    from zombi2.genomes.nucleotide_genome import NucleotideGenome, SegmentRegistry
    from zombi2.genomes.genome import IdManager
    params = types.SimpleNamespace(extension=None)   # novel-gene length falls back to self.extension
    g = NucleotideGenome(IdManager(), root_length=100, extension=0.5, registry=SegmentRegistry())
    g.originate(np.random.default_rng(1), params)    # seed: the root chromosome (gated on _seeded)
    seeded_len = g.total_length()
    assert seeded_len == 100 and g._seeded
    g.originate(np.random.default_rng(2), params)    # a novel gene, NOT a second 100-nt chromosome
    assert g.total_length() < seeded_len + 100
