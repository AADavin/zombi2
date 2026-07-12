"""Multiple chromosomes for the OrderedGenome model.

An ordered genome can carry more than one chromosome (``n_chromosomes``), circular or linear. The
gene events covered here stay *within* a chromosome (whole-chromosome fission / fusion / plasmid /
loss are the separate chromosome tier, tested in test_chromosome_events.py): ``draw_target`` picks a
chromosome (size-weighted) then a segment within it, and ``apply`` looks that chromosome up by
chrom_id. The single-chromosome default (``n_chromosomes=1, circular=True``) must stay
byte-identical to the pre-multichromosome engine.
"""

import numpy as np
import pytest

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    SharedRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.cli import main
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

    (The cross-commit byte-identity check against the frozen pre-refactor baseline lives in
    tests/test_chromosome_fingerprints.py; here we prove the new parameters are inert at defaults.)"""
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
        assert leaf.size() == sum(len(c) for c in leaf.chromosomes.values())  # size flattens
        assert leaf.size() == len(leaf.chromosome)              # back-compat flattened view agrees
        populated |= {i for i, c in enumerate(leaf.chromosomes.values()) if c}
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
    assert all(len(c) >= 1 for c in genome.chromosomes.values())  # every chromosome got some families


# --- 3. events stay within a chromosome -----------------------------------------------

def test_rearrangements_do_not_leak_genes_across_chromosomes():
    """After many duplications/inversions/transpositions each gene stays on the chromosome its
    family was seeded on (membership changes only via transfer/origination, not fired here)."""
    ids = IdManager()
    genome = OrderedGenome(ids, extension=0.6, n_chromosomes=3, circular=True)
    labels = ["A", "B", "C"]
    chroms = list(genome.chromosomes.values())               # stable positional handles
    for chrom, lab in zip(chroms, labels):
        for k in range(6):
            chrom.genes.append(OrderedGene(ids.new_gene(), f"{lab}{k}", 1))
    home = {g.family: lab for chrom, lab in zip(chroms, labels) for g in chrom.genes}

    rng = np.random.default_rng(0)
    params = TargetParams(extension=0.6)
    events = (EventType.DUPLICATION, EventType.INVERSION, EventType.TRANSPOSITION)
    for _ in range(600):
        event = events[int(rng.integers(3))]
        sel = genome.draw_target(event, rng, params)
        genome.apply(event, sel, rng, params)
        for chrom, lab in zip(chroms, labels):
            assert all(home[g.family] == lab for g in chrom.genes), \
                f"a gene leaked onto chromosome {lab} after {event}"
    # duplications actually happened somewhere (content grew) but every chromosome kept its families
    assert genome.size() > 18


# --- 4. linear vs circular ends -------------------------------------------------------

def _seed_single(circular, extension, n=10):
    ids = IdManager()
    genome = OrderedGenome(ids, extension=extension, n_chromosomes=1, circular=circular)
    chrom = next(iter(genome.chromosomes.values()))
    for k in range(n):
        chrom.genes.append(OrderedGene(ids.new_gene(), f"f{k}", 1))
    return genome


def test_linear_segments_never_wrap_the_origin():
    genome = _seed_single(circular=False, extension=0.9, n=10)
    rng = np.random.default_rng(1)
    params = TargetParams(extension=0.9)
    for _ in range(400):
        sel = genome.draw_target(EventType.INVERSION, rng, params)
        r = sel.region
        assert r.start + r.length <= 10                          # clamped, never crosses the end
        genes = next(iter(genome.chromosomes.values())).genes
        assert sel.genes == tuple(genes[r.start:r.start + r.length])


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
    for cidx, chrom in enumerate(parent.chromosomes.values()):
        for k in range(3):
            chrom.genes.append(OrderedGene(ids.new_gene(), f"f{cidx}_{k}", 1))
    child, mapping = parent.clone_reminting()

    assert isinstance(child, OrderedGenome)
    assert len(child.chromosomes) == 4
    assert child.circular is False
    assert all(c.circular is False for c in child.chromosomes.values())  # per-chromosome topology
    assert child.transposition_flip == 0.3
    assert child.extension == 0.6
    for pc, cc in zip(parent.chromosomes.values(), child.chromosomes.values()):
        assert [g.family for g in pc.genes] == [g.family for g in cc.genes]  # content copied
        assert all(a.gid != b.gid for a, b in zip(pc.genes, cc.genes))       # ids re-minted
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


# --- 6. transfer is the one event that moves a family across chromosomes --------------
#
# Rearrangements stay within a chromosome (section 3); a *transfer* inserts the segment into a
# recipient chromosome chosen uniformly, so it is the one *gene* event by which a family reaches a
# new chromosome (fission/fusion aside). draw_target -> extract_segment -> choose_insertion_point ->
# insert_segment.

def _one_gene_transfer_segment(ids, rng, family="X"):
    """Extract a single-gene transfer segment of ``family`` from a throwaway donor."""
    donor = OrderedGenome(ids, n_chromosomes=1, circular=True)
    next(iter(donor.chromosomes.values())).genes.append(OrderedGene(ids.new_gene(), family, 1))
    sel = donor.draw_target(EventType.TRANSFER, rng, TargetParams(), family=family)
    return donor.extract_segment(sel, rng)


def test_transfer_inserts_into_the_chosen_recipient_chromosome():
    """choose_insertion_point returns a (chromosome, position) tuple and insert_segment places the
    transferred family on exactly that chromosome — the cross-chromosome transfer path end to end."""
    ids = IdManager()
    rng = np.random.default_rng(0)
    segment = _one_gene_transfer_segment(ids, rng, family="X")

    recipient = OrderedGenome(ids, n_chromosomes=3, circular=True)
    for cidx, chrom in enumerate(recipient.chromosomes.values()):
        chrom.genes.append(OrderedGene(ids.new_gene(), f"seed{cidx}", 1))

    at = recipient.choose_insertion_point(segment, rng)
    assert isinstance(at, tuple) and len(at) == 2          # (chrom_id, position)
    cid, _pos = at
    recipient.insert_segment(segment, at, rng)

    on = {c for c, chrom in recipient.chromosomes.items()
          for g in chrom.genes if g.family == "X"}
    assert on == {cid}                                     # landed on exactly the chosen chromosome


def test_transfer_recipient_chromosome_is_uniform_not_size_weighted():
    """The recipient chromosome is chosen *uniformly*: every chromosome is reachable and a tiny
    chromosome receives transfers about as often as a large one (unlike the size-weighted choice
    that targets rearrangements). This pins down the Q4 behaviour."""
    ids = IdManager()
    rng = np.random.default_rng(0)
    segment = _one_gene_transfer_segment(ids, rng, family="X")

    recipient = OrderedGenome(ids, n_chromosomes=4, circular=True)
    rchroms = list(recipient.chromosomes.values())
    for chrom, n in zip(rchroms, (1, 2, 3, 40)):           # deliberately very uneven sizes
        for _ in range(n):
            chrom.genes.append(OrderedGene(ids.new_gene(), "f", 1))

    hits = [recipient.choose_insertion_point(segment, rng)[0] for _ in range(400)]
    assert set(hits) == set(recipient.chromosomes)         # every chromosome reachable (by chrom_id)
    # under uniform(4) each is expected ~100 times; the 1-gene and 40-gene chromosomes are hit at
    # comparable rates, which a size-weighted choice (0.02 vs 0.87) never would.
    assert hits.count(rchroms[0].chrom_id) > 40 and hits.count(rchroms[3].chrom_id) > 40


def test_transfer_moves_families_across_chromosomes_end_to_end():
    """A full simulation with heavy transfer: some family ends up on more than one chromosome
    within a single leaf, which (with no chromosome-tier events) can only happen via a transfer."""
    rates = SharedRates(duplication=0.2, loss=0.2, transfer=0.6, origination=0.1,
                        inversion=0.1, transposition=0.1)
    spanned = False
    for seed in (1, 2, 3, 4, 5):
        tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=4.0, seed=seed)
        g = simulate_genomes(tree, rates, initial_families=20, seed=seed,
                             genome_factory=_ordered(n_chromosomes=3, circular=True))
        for leaf in g.leaf_genomes.values():
            fam_chroms: dict = {}
            for cidx, chrom in enumerate(leaf.chromosomes.values()):
                for gene in chrom.genes:
                    fam_chroms.setdefault(gene.family, set()).add(cidx)
            if any(len(cs) >= 2 for cs in fam_chroms.values()):
                spanned = True
                break
        if spanned:
            break
    assert spanned, "no family spanned >1 chromosome despite heavy transfer across five seeds"


# --- 7. CLI wiring: --n-chromosomes / --linear-chromosomes ----------------------------

def _species_tree_file(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3", "--tips", "10", "--age", "3",
          "--seed", "1", "-o", str(sp)])
    return str(sp / "species_tree.nwk")


def test_cli_ordered_multichromosome_runs(tmp_path):
    """`genomes --genome-model ordered --n-chromosomes N --linear-chromosomes` runs and writes."""
    tree = _species_tree_file(tmp_path)
    gen = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "ordered",
               "--dup", "0.2", "--trans", "0.3", "--loss", "0.2", "--orig", "0.4",
               "--inversion", "0.2", "--transposition", "0.2", "--mean-length", "2",
               "--n-chromosomes", "4", "--linear-chromosomes",
               "--initial-families", "20", "--seed", "3", "--write", "profiles", "-o", str(gen)])
    assert rc == 0
    assert (gen / "Profiles.tsv").exists()


def test_cli_rejects_n_chromosomes_without_ordered(tmp_path):
    """--n-chromosomes only applies to the ordered model; other models must error cleanly."""
    tree = _species_tree_file(tmp_path)
    with pytest.raises(SystemExit):
        main(["genomes", "--tree", tree, "--n-chromosomes", "3",
              "--seed", "1", "-o", str(tmp_path / "x")])


def test_cli_rejects_linear_chromosomes_without_ordered(tmp_path):
    tree = _species_tree_file(tmp_path)
    with pytest.raises(SystemExit):
        main(["genomes", "--tree", tree, "--linear-chromosomes",
              "--seed", "1", "-o", str(tmp_path / "x")])


def test_cli_rejects_non_positive_n_chromosomes(tmp_path):
    tree = _species_tree_file(tmp_path)
    with pytest.raises(SystemExit):
        main(["genomes", "--tree", tree, "--genome-model", "ordered", "--n-chromosomes", "0",
              "--seed", "1", "-o", str(tmp_path / "x")])
