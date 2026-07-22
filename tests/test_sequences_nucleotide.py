"""Sequences on a **nucleotide** genome: every root block evolves — spacer as well as genes — and the
evolved blocks are put back together into one genome per extant lineage.

The load-bearing test is :func:`test_assembled_genome_is_the_root_sequence_permuted`. With the
substitution rate at zero nothing mutates, so every copy of a block still carries the block's founding
sequence, and a leaf's genome must equal the **root sequence read through its own trace-back** — the
per-nucleotide ancestry the genome level already records, computed here without touching the assembly.
That pins the two things an assembly gets silently wrong: which stretch of which block each piece is,
and which strand it is read on. A genome with either wrong still looks exactly like a genome.
"""

import pytest

from zombi2.genomes import simulate_genomes_nucleotide, simulate_genomes_unordered
from zombi2.genomes.events import node_from_label, node_label
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, jc69, lg
from zombi2.species import simulate_species_tree

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _run(*, seed=1, n_extant=6, **kw):
    """A small nucleotide genome run: three genes on a 600 bp circular chromosome, rearranged and with
    copy number changing, so the leaves' blocks are neither aligned with the root partition nor all on
    the forward strand."""
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=seed)
    params = dict(inversion=2.0, inversion_length=80, loss=0.5, loss_length=40, duplication=0.5,
                  duplication_length=40, root_length=600, genes=3, gene_length=90, seed=seed)
    params.update(kw)
    return simulate_genomes_nucleotide(sp, **params)


def _reference(genomes, sequences):
    """``{(source, position): base}`` over the whole root coordinate space, read off the evolved
    blocks. Only meaningful at a **zero** substitution rate, where every copy of a block is still its
    founding sequence — which is what makes this an oracle rather than a second assembly."""
    ref = {}
    for i, (src, a, b) in enumerate(genomes.root_blocks):
        seq = sequences.founding[i]
        assert len(seq) == b - a
        ref.update({(src, a + k): seq[k] for k in range(b - a)})
    return ref


def _traced(genomes, sequences, node_id):
    """The node's chromosomes spelled out from the genome level's per-nucleotide trace-back and the
    reference bases — the assembly's answer, arrived at the long way round."""
    ref = _reference(genomes, sequences)
    return {cid: "".join(ref[(src, p)] if strand == 1 else ref[(src, p)].translate(COMPLEMENT)
                         for (src, p, strand) in trace)
            for cid, trace in genomes.trace_back(node_id).items()}


# --- the assembly ----------------------------------------------------------------------------------

def test_assembled_genome_is_the_root_sequence_permuted():
    genomes = _run(seed=3)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.0, seed=3)
    for leaf in genomes.complete_tree.extant():
        assert r.genomes[node_label(leaf.id)] == _traced(genomes, r, leaf.id)


def test_inverted_stretches_come_back_reverse_complemented():
    # the same identity on a run built so the -1 branch is certainly exercised: with no reversed block
    # anywhere, the test above would pass with the strand ignored altogether
    genomes = _run(seed=7, inversion=6.0, n_extant=8)
    leaves = [n.id for n in genomes.complete_tree.extant()]
    assert any(b.strand == -1 for lid in leaves
               for c in genomes.genomes[lid].chromosomes for b in c.blocks)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.0, seed=7)
    for lid in leaves:
        assert r.genomes[node_label(lid)] == _traced(genomes, r, lid)


def test_every_extant_lineage_gets_a_genome_of_the_right_length():
    genomes = _run(seed=5)
    r = simulate_sequences(genomes, model=hky85(kappa=3.0), substitution=0.4, seed=5)
    assert set(r.genomes) == {node_label(n.id) for n in genomes.complete_tree.extant()}
    for leaf in genomes.complete_tree.extant():
        for chrom in genomes.genomes[leaf.id].chromosomes:
            assert len(r.genomes[node_label(leaf.id)][chrom.id]) == chrom.length


def test_a_genes_own_sequence_is_in_the_genome_it_sits_in():
    # the join the assembly exists for: the sequence evolved down a gene's tree is the sequence that
    # genome carries at that locus, oriented as the genome carries it
    genomes = _run(seed=11)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.3, seed=11)
    genic = {i for i, span in enumerate(genomes.root_blocks) if span in set(genomes.gene_spans.values())}
    checked = 0
    for leaf in genomes.complete_tree.extant():
        genome = r.genomes[node_label(leaf.id)]
        for cid, pieces in genomes.assembly(leaf.id).items():
            for (block, gene, start, end, _strand) in pieces:
                seq = r.alignments[block][f"g{gene}"]
                if block not in genic or (start, end) != (0, len(seq)):
                    continue                              # a gene taken whole, which is the only way
                assert seq in genome[cid] or seq.translate(COMPLEMENT)[::-1] in genome[cid]
                checked += 1
    assert checked, "the run left no gene in any leaf — pick another seed"


def test_material_no_extant_leaf_kept_cannot_be_assembled():
    # the root partition is cut from the extant leaves, so an ancestor holding material that later died
    # out everywhere has no block for it. Say so rather than return a genome with a hole in it.
    genomes = _run(seed=4, loss=3.0, n_extant=8)
    kept = {sp for n in genomes.complete_tree.extant() for sp in genomes.ancestry(n.id)}
    doomed = [nid for nid in sorted(genomes.genomes) if set(genomes.ancestry(nid)) - kept]
    assert doomed, "nothing died out for good in this run — pick another seed"
    with pytest.raises(ValueError, match="survives in no extant leaf"):
        genomes.assembly(doomed[0])
    for leaf in genomes.complete_tree.extant():           # an extant leaf always can be assembled
        genomes.assembly(leaf.id)


# --- every block evolves, at its own length --------------------------------------------------------

def test_each_block_evolves_at_its_own_length_in_bp():
    genomes = _run(seed=2)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.2, seed=2)
    assert set(r.alignments) == set(range(len(genomes.root_blocks)))
    for i, (_src, a, b) in enumerate(genomes.root_blocks):
        assert len(r.founding[i]) == b - a
        assert all(len(s) == b - a for s in r.alignments[i].values())


def test_length_is_rejected_and_a_protein_model_too():
    genomes = _run(seed=2, n_extant=4)
    with pytest.raises(ValueError, match="length does not apply"):
        simulate_sequences(genomes, model=jc69(), length=100, seed=1)
    with pytest.raises(ValueError, match="protein model"):
        simulate_sequences(genomes, model=lg(), seed=1)
    with pytest.raises(ValueError, match="protein model"):
        simulate_sequences(genomes, model=jc69(), intergene_model=lg(), seed=1)


# --- writing ---------------------------------------------------------------------------------------

def test_write_emits_one_fasta_per_extant_lineage(tmp_path):
    genomes = _run(seed=6, n_extant=4)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.2, seed=6)
    r.write(tmp_path, outputs=("genomes",))
    for leaf in genomes.complete_tree.extant():
        label = node_label(leaf.id)
        lines = (tmp_path / f"genome_{label}.fasta").read_text().splitlines()
        headers = [ln for ln in lines if ln.startswith(">")]
        assert headers == [f">{label}_chr{c.id}" for c in genomes.genomes[leaf.id].chromosomes]
        assert sum(len(ln) for ln in lines if not ln.startswith(">")) == genomes.genomes[leaf.id].length


def test_an_unordered_run_assembles_nothing(tmp_path):
    sp = simulate_species_tree(birth=1.0, n_extant=4, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.2, loss=0.2, initial_families=3, seed=1)
    r = simulate_sequences(g, model=jc69(), length=30, seed=1)
    assert r.genomes == {}
    r.write(tmp_path, outputs=("genomes",))
    assert not list(tmp_path.glob("genome_*.fasta"))


# --- ancestral genomes -----------------------------------------------------------------------------

def test_ancestral_genomes_are_the_genomes_that_were_really_there():
    # every internal node, against the ancestry the run recorded at it — the same oracle as for the
    # leaves, which is the claim: an ancestor you reconstruct is the one that was there.
    genomes = _run(seed=3)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.0, seed=3)
    internal = {node_label(i) for i, nd in genomes.complete_tree.nodes.items()
                if nd.children is not None}
    assert r.ancestral_genomes and set(r.ancestral_genomes) <= internal
    assert set(r.ancestral_genomes) & set(r.genomes) == set()      # the two halves never overlap
    for label, chroms in r.ancestral_genomes.items():
        assert chroms == _traced(genomes, r, node_from_label(label))


def test_the_rebuilt_root_is_the_genome_the_gff_seeded(tmp_path):
    # the round trip a reader cares about: declare a genome, evolve it, rebuild the ancestor from its
    # descendants, and get back what you declared — base for base, genes where the GFF put them.
    gff = tmp_path / "seed.gff"
    gff.write_text("##gff-version 3\n##sequence-region c 1 3000\n"
                   "c\tt\tgene\t201\t500\t.\t+\t.\tID=dnaA\n"
                   "c\tt\tgene\t1001\t1300\t.\t-\t.\tID=recA\n"
                   "c\tt\tgene\t2001\t2200\t.\t+\t.\tID=gyrB\n")
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    g = simulate_genomes_nucleotide(sp, gff=gff, inversion=2.0, inversion_length=400,
                                    duplication=0.4, duplication_length=300, loss=0.4,
                                    loss_length=300, seed=1)
    root = g.complete_tree.root
    assert not [x for x in g.rearrangements if x.lineage == root], \
        "this seed rearranges on the root branch, so n0 is no longer the seeded genome"
    r = simulate_sequences(g, model=jc69(), substitution=0.0, seed=1)

    # what we started with: the founding blocks in root coordinate order. Drawn before any event and
    # never seen by the recovery, so it is a reference and not a second assembly.
    started = "".join(r.founding[i] for i in range(len(g.root_blocks)))
    chrom = g.genomes[root].chromosomes[0]
    assert len(started) == chrom.length
    assert r.ancestral_genomes[node_label(root)][chrom.id] == started
    # and the leaves have moved on, so the match above is not a run in which nothing happened
    assert any(seq != started for chroms in r.genomes.values() for seq in chroms.values())


def test_an_ancestor_can_be_rebuilt_even_where_the_root_branch_moved_things():
    # n0 sits at the *end* of the root branch, so an event there leaves it already rearranged away
    # from the seeded layout. It is still exactly what the run held, which is what is checked.
    genomes = _run(seed=2, n_extant=8)
    root = genomes.complete_tree.root
    assert [x for x in genomes.rearrangements if x.lineage == root], "no root-branch event in this run"
    r = simulate_sequences(genomes, model=jc69(), substitution=0.0, seed=2)
    started = "".join(r.founding[i] for i in range(len(genomes.root_blocks)))
    rebuilt = r.ancestral_genomes[node_label(root)]
    assert rebuilt == _traced(genomes, r, root)
    assert list(rebuilt.values()) != [started]                     # the root branch did move things


def test_write_emits_the_ancestral_genomes_on_request(tmp_path):
    genomes = _run(seed=6, n_extant=4)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.2, seed=6)
    r.write(tmp_path)                                              # not in the default set
    assert not list(tmp_path.glob("genome_ancestral_*.fasta"))
    r.write(tmp_path, outputs=("ancestral_genomes",))
    for label, chroms in r.ancestral_genomes.items():
        lines = (tmp_path / f"genome_ancestral_{label}.fasta").read_text().splitlines()
        assert [ln for ln in lines if ln.startswith(">")] == [f">{label}_chr{c}" for c in chroms]
