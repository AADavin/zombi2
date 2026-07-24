"""The opt-in parallel engines (``parallel=``) at the genomes and sequences levels.

The contract, level-independent (SPEC: serial by default, parallel a separate opt-in engine):

- **worker-count invariance** — a parallel run is bit-identical for any worker count, because each unit
  (a gene family, a gene tree) draws from its own spawned RNG stream. So ``parallel=1`` (inline) equals
  ``parallel=2`` (a process pool) to the byte;
- **a separate engine** — that shared realisation *differs* from the serial reference for a given seed
  (decision A: the default path is left untouched, which the rest of the suite pins down);
- **still a valid run** — the per-family genome engine reproduces the strong invariants of the global
  one (its statistical equivalence is checked separately, offline, by a KS panel over many seeds);
- **loud, never silent** — an unsupported configuration falls back to the serial engine with a note,
  and a resolution with no per-family engine rejects the flag.
"""

import pytest

from zombi2._parallel import flatten_gene_tree, rebuild_gene_tree, resolve_workers
from zombi2.genomes import simulate_genomes_nucleotide, simulate_genomes_unordered
from zombi2.genomes.events import node_label
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.rates import modifiers as mod
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, jc69
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete


# --- shared scaffolding ---------------------------------------------------------------------------

def test_resolve_workers():
    import os
    assert resolve_workers(False) == 1 and resolve_workers(None) == 1
    assert resolve_workers(True) == (os.cpu_count() or 1)
    assert resolve_workers(4) == 4
    for bad in (0, -1, 2.5, "x"):
        with pytest.raises(ValueError):
            resolve_workers(bad)


def test_flatten_rebuild_preserves_structure_and_order():
    root = GeneNode("speciation", 0, 1.0, 0)
    root.children = [GeneNode("duplication", 1, 2.0, 1), GeneNode("extant", 2, 3.0, 2)]
    root.children[0].children = [GeneNode("extant", 1, 2.5, 3), GeneNode("loss", 1, 2.6, 4)]
    gt = rebuild_gene_tree(flatten_gene_tree(GeneTree(7, root, 0.5)))
    assert gt.family == 7 and gt.origination == 0.5
    assert [c.copy for c in gt.complete.children] == [1, 2]                 # sibling order kept
    assert [c.copy for c in gt.complete.children[0].children] == [3, 4]
    assert gt.complete.children[0].kind == "duplication"


def test_flatten_rebuild_is_recursion_safe():
    # a deep chain overflows the pickle recursion limit as an object graph; the flat form must not.
    root = GeneNode("extant", 0, 0.0, 0)
    cur = root
    for i in range(1, 8000):
        cur.children.append(GeneNode("extant", 0, float(i), i))
        cur = cur.children[0]
    rebuilt = rebuild_gene_tree(flatten_gene_tree(GeneTree(0, root, 0.0)))
    depth = 0
    n = rebuilt.complete
    while n.children:
        n = n.children[0]
        depth += 1
    assert depth == 7999


# --- sequences ------------------------------------------------------------------------------------

def _seq_fingerprint(r):
    return (tuple(sorted((f, tuple(sorted(a.items()))) for f, a in r.alignments.items())),
            tuple(sorted((f, tuple(sorted(a.items()))) for f, a in r.ancestral.items())),
            tuple(sorted(r.phylograms.items())), tuple(sorted(r.founding.items())))


@pytest.fixture(scope="module")
def genome_run():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=18, seed=1)
    return simulate_genomes_unordered(sp, duplication=0.4, transfer=0.2, loss=0.3,
                                      origination=0.1, initial_families=25, seed=2)


def test_sequences_parallel_is_worker_count_invariant(genome_run):
    def run(p):
        return simulate_sequences(genome_run, model=hky85(kappa=2.0), length=200, seed=7, parallel=p)
    inline = _seq_fingerprint(run(1))
    pooled = _seq_fingerprint(run(2))
    assert inline == pooled                                   # 1 (inline) == 2 (pool), to the byte


def test_sequences_parallel_differs_from_serial_but_matches_family_set(genome_run):
    serial = simulate_sequences(genome_run, model=hky85(kappa=2.0), length=200, seed=7)
    par = simulate_sequences(genome_run, model=hky85(kappa=2.0), length=200, seed=7, parallel=2)
    assert set(serial.alignments) == set(par.alignments)      # same families, valid run
    assert _seq_fingerprint(serial) != _seq_fingerprint(par)  # different (valid) draw — decision A


def test_sequences_parallel_nucleotide_assembles_every_node():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=3)
    gen = simulate_genomes_nucleotide(sp, loss=0.5, loss_length=40, duplication=0.5,
                                      duplication_length=40, root_length=400, genes=2,
                                      gene_length=90, seed=3)
    a = simulate_sequences(gen, model=hky85(kappa=2.0), intergene_model=jc69(), seed=9, parallel=1)
    b = simulate_sequences(gen, model=hky85(kappa=2.0), intergene_model=jc69(), seed=9, parallel=2)
    assert a.genomes == b.genomes and a.initial_genome == b.initial_genome
    assert set(a.genomes) == {node_label(i) for i in sp.complete_tree.nodes}


# --- genomes: the per-family engine ---------------------------------------------------------------

def _extant_leaves(node):
    if node is None:
        return 0
    if node.is_leaf:
        return 1 if node.kind == "extant" else 0
    return sum(_extant_leaves(c) for c in node.children)


def _gen_fingerprint(r):
    ev = sorted((round(e.time, 9), e.kind, e.lineage, e.family, e.copy, e.parent,
                 e.recipient, e.donor) for e in r.events)
    gen = {i: sorted((c.family, c.id) for c in r.genomes[i]) for i in sorted(r.genomes)}
    return ev, tuple(sorted(gen.items()))


@pytest.fixture(scope="module")
def species_for_genomes():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)


def test_genomes_parallel_is_worker_count_invariant(species_for_genomes):
    kw = dict(duplication=0.5, transfer=0.3, loss=0.4, origination=0.2, initial_families=30, seed=5)
    inline = _gen_fingerprint(simulate_genomes_unordered(species_for_genomes, parallel=1, **kw))
    pooled = _gen_fingerprint(simulate_genomes_unordered(species_for_genomes, parallel=2, **kw))
    assert inline == pooled


def test_genomes_parallel_differs_from_serial(species_for_genomes):
    kw = dict(duplication=0.5, transfer=0.3, loss=0.4, origination=0.2, initial_families=30, seed=5)
    serial = _gen_fingerprint(simulate_genomes_unordered(species_for_genomes, **kw))
    par = _gen_fingerprint(simulate_genomes_unordered(species_for_genomes, parallel=2, **kw))
    assert serial != par


def test_genomes_parallel_is_a_valid_run(species_for_genomes):
    # the merge produces a coherent run: every node has a genome, copy ids are globally unique, the
    # gene trees build, and the strong invariant holds — surviving gene-tree leaves == extant copies.
    sp = species_for_genomes
    r = simulate_genomes_unordered(sp, duplication=0.5, transfer=0.3, loss=0.4, origination=0.2,
                                   initial_families=30, seed=5, parallel=2)
    assert set(r.genomes) == set(sp.complete_tree.nodes)
    born = [e.copy for e in r.events
            if e.kind in ("origination", "duplication", "transfer", "speciation")]
    assert len(born) == len(set(born))                       # copy ids globally unique
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    for fam, gt in r.gene_trees.items():                     # builds without collision + strong invariant
        copies = sum(r.profiles.counts.get((fam, s), 0) for s in extant_sp)
        assert _extant_leaves(gt.extant) == copies


def test_genomes_parallel_named_families_survive(species_for_genomes):
    r = simulate_genomes_unordered(species_for_genomes, origination=0.2, loss=0.1,
                                   families=["toxin", "operon"], initial_families=10,
                                   seed=5, parallel=2)
    assert set(r.family_names) == {"toxin", "operon"}
    assert all(fid in {c.family for g in r.genomes.values() for c in g} for fid in r.family_names.values())


def test_genomes_parallel_falls_back_on_driven_rate(species_for_genomes, capsys):
    # a driven rate has no per-family engine yet: the run must announce the fallback and return the
    # *serial* result unchanged, not a quietly-degraded parallel one.
    sp = species_for_genomes
    habitat = simulate_discrete(sp, states=["a", "b"], switch=0.8, seed=2)
    kw = dict(duplication=0.5, loss=0.25 * mod.DrivenBy(habitat, {"a": 2.0, "b": 1.0}),
              origination=0.2, initial_families=25, seed=5)
    par = simulate_genomes_unordered(sp, parallel=4, **kw)
    note = capsys.readouterr().out
    assert "not applied" in note and "driven" in note
    serial = simulate_genomes_unordered(sp, **kw)
    assert _gen_fingerprint(par) == _gen_fingerprint(serial)


def test_genomes_parallel_true_uses_all_cores(species_for_genomes):
    # parallel=True (every core) must agree with an explicit worker count — same spawned streams.
    kw = dict(duplication=0.4, loss=0.3, origination=0.2, initial_families=15, seed=5)
    assert _gen_fingerprint(simulate_genomes_unordered(species_for_genomes, parallel=True, **kw)) == \
           _gen_fingerprint(simulate_genomes_unordered(species_for_genomes, parallel=2, **kw))
