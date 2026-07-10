"""Stage 1 multiple chromosomes for the OrderedGenome model.

An ordered genome can carry more than one chromosome (``n_chromosomes``), circular or
linear. Every event stays *within* a chromosome (no translocation/fission/fusion yet):
``draw_target`` picks a chromosome (size-weighted) then a segment within it, and ``apply``
indexes into that chromosome. The single-chromosome default (``n_chromosomes=1,
circular=True``) must stay byte-identical to the pre-multichromosome engine.
"""

import numpy as np

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    SharedRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes.genome import IdManager, OrderedGene
from zombi2.genomes.events import EventType, TargetParams


def _ordered(**kw):
    return lambda ids: OrderedGenome(ids, extension=0.6, **kw)


def _chromosomes(genomes):
    """Leaf-name -> flattened list of (family, orientation), a stable fingerprint."""
    return {
        leaf.name: [(g.family, g.orientation) for g in genome.chromosome]
        for leaf, genome in sorted(genomes.leaf_genomes.items(), key=lambda kv: kv[0].name)
    }


def _rates():
    return SharedRates(duplication=0.4, loss=0.3, transfer=0.2, origination=0.1,
                       inversion=0.3, transposition=0.3)


# --- 1. byte-identity: n_chromosomes=1, circular=True == the single-chromosome engine ----

def test_default_matches_explicit_single_chromosome():
    """Passing the new defaults explicitly must change nothing (same seed -> identical leaves).

    (A stronger cross-commit byte-identity check against the pre-change base is done outside
    the test suite; here we prove the new parameters are inert at their defaults.)"""
    rates = _rates()
    for seed in (1, 2, 3, 7, 42):
        tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=seed)
        a = simulate_genomes(tree, rates, initial_families=18, seed=seed,
                             genome_factory=lambda ids: OrderedGenome(ids, extension=0.6))
        b = simulate_genomes(tree, rates, initial_families=18, seed=seed,
                             genome_factory=_ordered(n_chromosomes=1, circular=True))
        assert _chromosomes(a) == _chromosomes(b)


# --- 2. multi-chromosome structure ----------------------------------------------------

def test_four_chromosomes_structure_and_distribution():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=4)
    rates = _rates()
    g = simulate_genomes(tree, rates, initial_families=24, seed=4,
                         genome_factory=_ordered(n_chromosomes=4, circular=True))
    populated = set()
    for leaf in g.leaf_genomes.values():
        assert isinstance(leaf, OrderedGenome)
        assert len(leaf.chromosomes) == 4                       # structure preserved to the leaves
        assert leaf.size() == sum(len(c) for c in leaf.chromosomes)  # size flattens
        assert leaf.size() == len(leaf.chromosome)              # back-compat flattened view agrees
        populated |= {i for i, c in enumerate(leaf.chromosomes) if c}
    assert populated == {0, 1, 2, 3}                            # genes live on every chromosome


def test_root_originations_spread_across_chromosomes():
    """The root's repeated originate() calls distribute the initial families across chromosomes."""
    ids = IdManager()
    genome = OrderedGenome(ids, n_chromosomes=5, circular=True)
    rng = np.random.default_rng(0)
    params = TargetParams()
    for _ in range(30):
        genome.originate(rng, params)
    assert genome.size() == 30
    assert all(len(c) >= 1 for c in genome.chromosomes)         # every chromosome got some families


# --- 3. events stay within a chromosome -----------------------------------------------

def test_rearrangements_do_not_leak_genes_across_chromosomes():
    """After many duplications/inversions/transpositions each gene stays on the chromosome its
    family was seeded on (membership changes only via transfer/origination, not fired here)."""
    ids = IdManager()
    genome = OrderedGenome(ids, extension=0.6, n_chromosomes=3, circular=True)
    labels = {0: "A", 1: "B", 2: "C"}
    for cidx, lab in labels.items():
        for k in range(6):
            genome.chromosomes[cidx].append(OrderedGene(ids.new_gene(), f"{lab}{k}", 1))
    home = {g.family: lab for cidx, lab in labels.items() for g in genome.chromosomes[cidx]}

    rng = np.random.default_rng(0)
    params = TargetParams(extension=0.6)
    events = (EventType.DUPLICATION, EventType.INVERSION, EventType.TRANSPOSITION)
    for _ in range(600):
        event = events[int(rng.integers(3))]
        sel = genome.draw_target(event, rng, params)
        genome.apply(event, sel, rng, params)
        for cidx, lab in labels.items():
            assert all(home[g.family] == lab for g in genome.chromosomes[cidx]), \
                f"a gene leaked onto chromosome {cidx} ({lab}) after {event}"
    # duplications actually happened somewhere (content grew) but every chromosome kept its families
    assert genome.size() > 18


# --- 4. linear vs circular ends -------------------------------------------------------

def _seed_single(circular, extension, n=10):
    ids = IdManager()
    genome = OrderedGenome(ids, extension=extension, n_chromosomes=1, circular=circular)
    for k in range(n):
        genome.chromosomes[0].append(OrderedGene(ids.new_gene(), f"f{k}", 1))
    return genome


def test_linear_segments_never_wrap_the_origin():
    genome = _seed_single(circular=False, extension=0.9, n=10)
    rng = np.random.default_rng(1)
    params = TargetParams(extension=0.9)
    for _ in range(400):
        sel = genome.draw_target(EventType.INVERSION, rng, params)
        r = sel.region
        assert r.start + r.length <= 10                          # clamped, never crosses the end
        assert sel.genes == tuple(genome.chromosomes[0][r.start:r.start + r.length])


def test_circular_segments_may_wrap_the_origin():
    genome = _seed_single(circular=True, extension=0.9, n=10)
    rng = np.random.default_rng(1)
    params = TargetParams(extension=0.9)
    wrapped = sum(
        (lambda r: r.start + r.length > 10)(
            genome.draw_target(EventType.INVERSION, rng, params).region)
        for _ in range(400)
    )
    assert wrapped > 0                                           # a circular ring does wrap


# --- 5. speciation preserves the karyotype -------------------------------------------

def test_clone_reminting_preserves_n_chromosomes_and_circular():
    ids = IdManager()
    parent = OrderedGenome(ids, extension=0.6, transposition_flip=0.3,
                           n_chromosomes=4, circular=False)
    for cidx in range(4):
        for k in range(3):
            parent.chromosomes[cidx].append(OrderedGene(ids.new_gene(), f"f{cidx}_{k}", 1))
    child, mapping = parent.clone_reminting()

    assert isinstance(child, OrderedGenome)
    assert len(child.chromosomes) == 4
    assert child.circular is False
    assert child.transposition_flip == 0.3
    assert child.extension == 0.6
    for pc, cc in zip(parent.chromosomes, child.chromosomes):
        assert [g.family for g in pc] == [g.family for g in cc]   # per-chromosome content copied
        assert all(a.gid != b.gid for a, b in zip(pc, cc))        # ids re-minted (fresh lineages)
    assert len(mapping) == parent.size()


def test_multichromosome_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=5)
    rates = _rates()
    a = simulate_genomes(tree, rates, initial_families=16, seed=6,
                         genome_factory=_ordered(n_chromosomes=3, circular=False))
    b = simulate_genomes(tree, rates, initial_families=16, seed=6,
                         genome_factory=_ordered(n_chromosomes=3, circular=False))
    assert _chromosomes(a) == _chromosomes(b)
    assert all(len(leaf.chromosomes) == 3 for leaf in a.leaf_genomes.values())
