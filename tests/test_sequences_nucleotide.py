"""Sequences on a **nucleotide** genome: every root block evolves — spacer as well as genes — and the
evolved blocks are put back together into one genome per extant lineage.

The load-bearing test is :func:`test_assembled_genome_is_the_root_sequence_permuted`. With the
substitution rate at zero nothing mutates, so every copy of a block still carries the block's founding
sequence, and a leaf's genome must equal the **root sequence read through its own trace-back** — the
per-nucleotide ancestry the genome level already records, computed here without touching the assembly.
That pins the two things an assembly gets silently wrong: which stretch of which block each piece is,
and which strand it is read on. A genome with either wrong still looks exactly like a genome.
"""

import collections

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


def test_every_lineage_gets_a_genome_of_the_right_length():
    genomes = _run(seed=5)
    r = simulate_sequences(genomes, model=hky85(kappa=3.0), substitution=0.4, seed=5)
    assert set(r.genomes) == {node_label(i) for i in genomes.complete_tree.nodes}
    for node_id, genome in genomes.genomes.items():
        for chrom in genome.chromosomes:
            assert len(r.genomes[node_label(node_id)][chrom.id]) == chrom.length


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
            for (block, gene, _strand) in pieces:
                if block not in genic:
                    continue
                seq = r.alignments[block][f"g{gene}"]
                assert seq in genome[cid] or seq.translate(COMPLEMENT)[::-1] in genome[cid]
                checked += 1
    assert checked, "the run left no gene in any leaf — pick another seed"


def test_material_no_extant_leaf_kept_is_still_reconstructed():
    # the partition is cut at every node, not at the survivors, so a lineage holding material that
    # later died out everywhere still has blocks for it — and a genome, exact like any other
    genomes = _run(seed=4, loss=3.0, n_extant=8)
    kept = {sp for n in genomes.complete_tree.extant() for sp in genomes.ancestry(n.id)}
    doomed = [nid for nid in sorted(genomes.genomes) if set(genomes.ancestry(nid)) - kept]
    assert doomed, "nothing died out for good in this run — pick another seed"
    r = simulate_sequences(genomes, model=jc69(), substitution=0.0, seed=4)
    rebuilt = r.genomes
    for nid in doomed:
        assert rebuilt[node_label(nid)] == _traced(genomes, r, nid)


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
    assert internal <= set(r.genomes)                              # one map, covering every node
    assert set(r.genomes) == {node_label(i) for i in genomes.complete_tree.nodes}
    for label, chroms in r.genomes.items():
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
    assert r.genomes[node_label(root)][chrom.id] == started
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
    rebuilt = r.genomes[node_label(root)]
    assert rebuilt == _traced(genomes, r, root)
    assert list(rebuilt.values()) != [started]                     # the root branch did move things


def test_write_emits_a_genome_for_every_node_by_default(tmp_path):
    # no node is a special case: an ancestor and an extinct lineage are written like a surviving tip,
    # each named by whose genome it is
    genomes = _run(seed=6, n_extant=4)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.2, seed=6)
    r.write(tmp_path)
    assert not list(tmp_path.glob("*ancestral*"))
    assert {p.name for p in tmp_path.glob("genome_n*.fasta")} == \
           {f"genome_{node_label(i)}.fasta" for i in genomes.complete_tree.nodes}
    for label, chroms in r.genomes.items():
        lines = (tmp_path / f"genome_{label}.fasta").read_text().splitlines()
        assert [ln for ln in lines if ln.startswith(">")] == [f">{label}_chr{c}" for c in chroms]


# --- the two numbering schemes ---------------------------------------------------------------------

def test_a_nucleotide_run_names_its_files_blocks_not_families(tmp_path):
    # the keys here are block indices, not gene family ids, and the files have to say so: a gene
    # family id is a valid-looking int that names a different locus
    genomes = _run(seed=6, n_extant=4)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.2, seed=6)
    assert r.unit == "block"
    r.write(tmp_path, outputs=("alignments", "phylograms", "ancestral", "founding"))
    names = {p.name for p in tmp_path.iterdir()}
    assert not [n for n in names if "fam" in n]
    assert "block0.fasta" in names and "phylogram_block0_complete.nwk" in names
    assert any(n.startswith("sequences_ancestral_block") for n in names)
    assert ">block0\n" in (tmp_path / "sequences_founding.fasta").read_text()


def test_an_unordered_run_still_names_them_families(tmp_path):
    sp = simulate_species_tree(birth=1.0, n_extant=4, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.2, loss=0.2, initial_families=3, seed=1)
    r = simulate_sequences(g, model=jc69(), length=30, seed=1)
    assert r.unit == "family"
    r.write(tmp_path, outputs=("alignments", "phylograms"))
    names = {p.name for p in tmp_path.iterdir()}
    assert "fam1.fasta" in names and "phylogram_fam0_complete.nwk" in names
    assert not [n for n in names if n.startswith("block")]


def test_block_of_joins_a_gene_family_to_its_sequences():
    genomes = _run(seed=11)
    r = simulate_sequences(genomes, model=jc69(), substitution=0.3, seed=11)
    fam = min(genomes.gene_trees)
    block = genomes.block_of(fam)
    assert genomes.root_blocks[block] == genomes.gene_spans[fam]
    src, a, b = genomes.gene_spans[fam]
    assert all(len(s) == b - a for s in r.alignments[block].values())
    with pytest.raises(KeyError):
        genomes.block_of(99999)                            # never declared


def test_a_seeded_gene_always_has_a_block_however_hard_the_run():
    # the initial genome carries every seeded gene and votes on the partition, so a seeded gene has a
    # block even when it is lost in every lineage — which is the point of keeping it
    genomes = _run(seed=4, loss=3.0, n_extant=8)
    seeded = {e.source for e in genomes.events
              if type(e).__name__ == "Origination" and e.time == 0.0}
    declared = [f for f, (src, _a, _b) in genomes.gene_spans.items() if src in seeded]
    assert declared
    for fam in declared:
        assert genomes.root_blocks[genomes.block_of(fam)] == genomes.gene_spans[fam]


def test_block_of_says_so_when_a_gene_left_nothing_behind():
    # a de-novo gene born and lost on the same branch is carried by no genome in the run at all, so
    # nothing was reconstructed for it. That is the one way a declared family has no block.
    genomes = _run(seed=2, loss=4.0, loss_length=60, origination=3.0, origination_length=50,
                   n_extant=8)
    gone = sorted(set(genomes.gene_spans) - set(genomes.gene_trees))
    assert gone, "every declared gene left a block in this run — pick another seed"
    with pytest.raises(LookupError, match="no recovered root block"):
        genomes.block_of(gone[0])


# --- everything at once ----------------------------------------------------------------------------

def _busy(seed=2, n_extant=10):
    """A multi-chromosome run with every event kind on — the scaled-down twin of the stress run: two
    replicons of different shapes, rearrangement, copy-number change, transfer between lineages, de
    novo genes, and the chromosome tier splitting and merging replicons."""
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=n_extant, seed=seed)
    return simulate_genomes_nucleotide(
        sp.complete_tree, chromosomes=[(900, "circular"), (400, "linear")], genes=6, gene_length=60,
        inversion=3.0, inversion_length=120, translocation=1.0, transposition=1.0,
        duplication=1.5, duplication_length=100, loss=1.5, loss_length=100,
        transfer=1.0, transfer_length=100, origination=0.3, fission=0.1, fusion=0.1, seed=seed)


def test_every_rebuildable_genome_is_exact_under_every_event_kind():
    g = _busy()
    assert {type(e).__name__ for e in g.events} >= {"Duplication", "Loss", "Transfer", "Origination"}
    assert {type(x).__name__ for x in g.rearrangements} == {"Inversion", "Translocation",
                                                            "Transposition"}
    assert {e.kind for e in g.chromosome_events} >= {"fission", "fusion"}
    r = simulate_sequences(g, model=jc69(), substitution=0.0, seed=2)
    rebuilt = r.genomes
    assert len(rebuilt) > 10
    for label, chroms in rebuilt.items():
        assert chroms == _traced(g, r, node_from_label(label)), label
    # and it is not a run where every genome came out the same
    assert len({tuple(sorted(v.items())) for v in rebuilt.values()}) == len(rebuilt)


def test_a_partial_loss_leaves_no_fragment_to_go_wrong():
    # A loss need only *overlap* a block to end that copy's lineage for it. When the partition was cut
    # from the survivors alone, that left ancestors holding a fragment of a block whose genealogy had
    # no lineage for it — no sequence to read, and a KeyError deep in the assembly. Cutting at every
    # node removes the fragment itself: a loss's own endpoints are breakpoints, so it always covers a
    # whole number of blocks, and a copy is only ever ended for blocks it really lost.
    g = _busy()
    tips = g._recover_blocks()[2]
    dead = {key for key, gene in tips.items() if gene is None}
    assert dead, "no copy was ended by a loss in this run — pick another seed"
    # the invariant: no node still carrying a block is one of those. A copy's blocks are disjoint in
    # source coordinates, so if a node still holds all of a block under copy c, no loss of c can have
    # taken any of it — and it holds *all* of it, because the node voted on the cuts.
    for node_id, genome in g.genomes.items():
        for chrom in genome.chromosomes:
            for blk in chrom.blocks:
                for i, (s, a, b) in enumerate(g.root_blocks):
                    if s == blk.source and blk.start <= a and b <= blk.end:
                        assert (i, blk.copy) not in dead
    for node_id in sorted(g.genomes):                               # so every node assembles
        g.assembly(node_id)


# --- the initial genome ----------------------------------------------------------------------------

def test_the_initial_genome_is_what_the_run_was_seeded_with(tmp_path):
    """The genome at the START of the root branch, which no node holds: the root node sits at the END
    of it, and the stem is real simulated time. It votes on the partition like every other genome, so
    it reconstructs exactly — no matter what the stem did."""
    genomes = _run(seed=3, inversion=8.0, loss=2.0, n_extant=8)
    root = genomes.complete_tree.root
    on_stem = [x for x in genomes.rearrangements if x.lineage == root]
    assert on_stem, "the stem was quiet in this run — pick another seed"

    r = simulate_sequences(genomes, model=jc69(), substitution=0.0, seed=3)
    chrom = genomes.initial_genome.chromosomes[0]
    # seeded: the genes evenly spaced on the + strand, in coordinate order, nothing yet rearranged
    assert [b.strand for b in chrom.blocks] == [1] * len(chrom.blocks)
    assert all(chrom.blocks[i].end == chrom.blocks[i + 1].start for i in range(len(chrom.blocks) - 1))
    seq = r.initial_genome[chrom.id]
    assert len(seq) == chrom.length
    assert seq != r.genomes[node_label(root)][chrom.id]     # the stem moved things

    # every declared gene sits at exactly its declared coordinates, carrying its own founding sequence
    for fam, (_src, a, b) in genomes.gene_spans.items():
        assert seq[a:b] == r.founding[genomes.block_of(fam)]

    r.write(tmp_path, outputs=("initial_genome",))
    text = (tmp_path / "genome_initial.fasta").read_text()
    assert f">initial_chr{chrom.id}\n" in text
    assert "".join(text.split(">")[1].splitlines()[1:]) == seq


def test_the_genome_level_writes_the_initial_mosaic_in_its_own_file(tmp_path):
    genomes = _run(seed=3, n_extant=4)
    genomes.write(tmp_path, outputs=("initial_genome", "blocks"))
    rows = (tmp_path / "initial_genome.tsv").read_text().splitlines()
    assert rows[0] == "chromosome\tposition\tsource\tstart\tend\tstrand\tcopy\tgene"
    assert len(rows) - 1 == sum(len(c.blocks) for c in genomes.initial_genome.chromosomes)
    # and it is not smuggled into blocks.tsv, whose every lineage is a real node
    labels = {ln.split("\t")[0] for ln in (tmp_path / "blocks.tsv").read_text().splitlines()[1:]}
    assert labels == {node_label(i) for i in genomes.genomes}


# --- reading a written run back ---------------------------------------------------------------------

def test_a_written_run_reads_back_to_the_same_recovery(tmp_path):
    """The handoff has to carry the blocks, not just the event log: at this resolution the sequences
    evolve down a tree per block, and the blocks come from the genomes themselves."""
    from zombi2.genomes.nucleotide import read_nucleotide_genomes

    g = _busy(seed=2)
    g.write(tmp_path)
    back = read_nucleotide_genomes(tmp_path, g.complete_tree)
    assert back.events == g.events
    assert back.root_blocks == g.root_blocks
    assert back.gene_spans == g.gene_spans and back.gene_strands == g.gene_strands
    assert {i: t.to_newick("complete") for i, t in back.block_trees.items()} == \
           {i: t.to_newick("complete") for i, t in g.block_trees.items()}
    # and what it writes is byte-identical to what it read
    out = tmp_path / "again"
    back.write(out, outputs=("events", "blocks", "genes", "initial_genome"))
    for name in ("genome_events.tsv", "blocks.tsv", "genes.tsv", "initial_genome.tsv"):
        assert (out / name).read_text() == (tmp_path / name).read_text(), name


def test_reading_back_regroups_multi_row_events_correctly(tmp_path):
    """An event spanning several ancestral intervals is several rows, and every copy re-minted at one
    node shares a time and lineage. Both have to regroup to exactly what was written."""
    from zombi2.genomes.nucleotide import Duplication, Loss, Speciation, read_nucleotide_genomes

    g = _busy(seed=6, n_extant=12)
    kinds = collections.Counter(type(e).__name__ for e in g.events)
    assert kinds["Speciation"] > kinds["Duplication"], "no speciation crowd in this run"
    assert any(len(e.lost) > 1 for e in g.events if isinstance(e, Loss)) or \
           any(len(e.copied) > 1 for e in g.events if isinstance(e, Duplication)), \
        "no multi-interval event in this run — pick another seed"
    at_one_node = collections.Counter(
        (e.time, e.lineage) for e in g.events if isinstance(e, Speciation))
    assert max(at_one_node.values()) > 1, "no node re-minted several copies — pick another seed"
    g.write(tmp_path)
    assert read_nucleotide_genomes(tmp_path, g.complete_tree).events == g.events
