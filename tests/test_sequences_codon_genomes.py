"""Codon models on a nucleotide genome — ``model=gy94(...)`` where the genes sit on a genome.

On gene families a codon model has one sequence in one frame. On a genome a gene has a strand and a
place, may share bases with the gene beside it, and may not be a protein-coding gene at all: a
``gene`` row of a GFF is also what an RNA gene and a pseudogene are. So each gene block is cut into
parts — codons in each coding gene's own frame and strand, the bases no frame owns held, plain DNA
where a gene is not coding — and these check that every one of those does what it says.
"""

import random

import pytest

from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.sequences import gy94, hky85, simulate_sequences
from zombi2.sequences.substitution_models import GENETIC_CODES
from zombi2.species import simulate_species_tree

_TCAG = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AMINO = {code: dict(zip(_TCAG, table)) for code, table in GENETIC_CODES.items()}
_RC = str.maketrans("ACGT", "TGCA")


def _rc(seq):
    return seq.translate(_RC)[::-1]


def _codons(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def _frame_is_whole(seq, code=11):
    """A whole number of codons, ending in a stop, with no stop before it."""
    amino = AMINO[code]
    cs = _codons(seq)
    return len(seq) % 3 == 0 and amino[cs[-1]] == "*" and not any(amino[c] == "*" for c in cs[:-1])


# --- a small genome written by hand -------------------------------------------------------------------
#
# Eight genes in the ways a real annotation holds them: on either strand, sharing four bases at an
# operon junction, nested, and not coding at all.

def _build():
    rng = random.Random(11)
    sense = [c for c in _TCAG if AMINO[11][c] != "*"]

    def spacer(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def reading_frame(n):                        # ATG, n codons that are not stops, a stop
        return "ATG" + "".join(rng.choice(sense) for _ in range(n)) + "TAA"

    dna, genes = "", {}

    def lay(name, seq, strand="+", back=0):
        nonlocal dna
        start = len(dna) - back
        dna += seq[back:]
        genes[name] = (start, start + len(seq), strand)

    dna += spacer(60)
    lay("A", reading_frame(40))                                  # a gene on the plus strand
    dna += spacer(30)
    lay("B", _rc(reading_frame(50)), "-")                        # one on the minus strand
    dna += spacer(30)
    # C and D share four bases, ATGA: the end of C's last amino acid and its stop TGA are D's ATG
    c = "ATG" + "".join(rng.choice(sense) for _ in range(30)) + "GCA" + "TGA"
    lay("C", c)
    lay("D", "ATGA" + "AA" + "".join(rng.choice(sense) for _ in range(35)) + "TAA", back=4)
    dna += spacer(30)
    lay("E", spacer(100))                                        # not a whole number of codons
    dna += spacer(30)
    lay("F", "ATG" + "".join(rng.choice(sense) for _ in range(10)) + "TAG"   # a stop in the middle
        + "".join(rng.choice(sense) for _ in range(10)) + "TAA")
    dna += spacer(30)
    g = reading_frame(60)
    lay("G", g)
    # nested in G, on the other strand. 31 bases, so that it cannot read as a frame: 30 random bases
    # are a whole number of codons, and end in a stop often enough to pass for a gene by chance
    genes["H"] = (genes["G"][0] + 60, genes["G"][0] + 91, "-")
    dna += spacer(60)
    return dna, genes


DNA, GENES = _build()
CODING = {"A", "B", "C", "D", "G"}


@pytest.fixture
def files(tmp_path):
    gff = tmp_path / "genome.gff"
    gff.write_text(f"##gff-version 3\n##sequence-region c 1 {len(DNA)}\n" + "".join(
        f"c\tt\tgene\t{a + 1}\t{b}\t.\t{strand}\t.\tID={name}\n"
        for name, (a, b, strand) in sorted(GENES.items(), key=lambda kv: kv[1][0])), encoding="utf-8")
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(f">c\n{DNA}\n", encoding="utf-8")
    return gff, fasta


#: events that move, copy and remove whole stretches, on a tree deep enough to change every gene
BUSY = dict(inversion=3.0, inversion_extent=400, duplication=0.5, duplication_extent=300,
            loss=0.3, loss_extent=300, transfer=1.0, transfer_extent=300, transfer_to="distance",
            replacement=0.5)


def _genome(files, *, seed=2, n_extant=10, **params):
    gff, fasta = files
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=seed)
    return simulate_genomes_nucleotide(sp, gff=gff, fasta=fasta, seed=seed, **params)


def _gene_copies(g, r):
    """``(name, sequence on the gene's own strand)`` for every gene of every node's genome, read
    through that genome's own GFF — the coordinates a reader of the files would use."""
    names = g.complete_tree.labels()
    name_of = {fam: name for name, fam in g.gene_names.items()}
    for node, genome in g.node_genomes.items():
        label = names[node]
        for line in g._gff(label, genome).splitlines():
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            attrs = dict(x.split("=", 1) for x in cols[8].split(";"))
            chrom = int(cols[0].rsplit("_chr", 1)[1])
            seq = r.node_genomes[label][chrom][int(cols[3]) - 1:int(cols[4])]
            yield name_of[int(attrs["family"])], seq if cols[6] == "+" else _rc(seq)


def _original(name):
    a, b, strand = GENES[name]
    return DNA[a:b] if strand == "+" else _rc(DNA[a:b])


# --- which genes are coding ---------------------------------------------------------------------------

def test_a_gene_is_coding_when_its_sequence_reads_as_one(files):
    g = _genome(files)
    r = simulate_sequences(g, model=gy94(omega=0.1), divergence=0.1, seed=1)
    called = {name: r.coding_genes[fam] for name, fam in g.gene_names.items()}
    assert {name for name, yes in called.items() if yes} == CODING
    assert set(called) == set(GENES)
    assert r.summary()["genes"] == {"coding": 5, "not_coding": 3}
    # a run under a nucleotide model reads no frame, and says nothing about one
    plain = simulate_sequences(g, model=hky85(), divergence=0.1, seed=1)
    assert plain.coding_genes == {} and "genes" not in plain.summary()


def test_the_genetic_code_decides_what_a_stop_is(tmp_path):
    # TGA in the frame: a stop in code 11, tryptophan in code 4
    dna = "ACGT" * 5 + "ATG" + "GCA" * 5 + "TGA" + "GCT" * 5 + "TAA" + "ACGT" * 5
    (tmp_path / "g.gff").write_text(f"##sequence-region c 1 {len(dna)}\n"
                                    f"c\tt\tgene\t21\t{len(dna) - 20}\t.\t+\t.\tID=x\n")
    (tmp_path / "g.fasta").write_text(f">c\n{dna}\n")
    g = _genome((tmp_path / "g.gff", tmp_path / "g.fasta"))
    fam = g.gene_names["x"]
    assert simulate_sequences(g, model=gy94(omega=0.1), seed=1).coding_genes[fam] is False
    by_four = simulate_sequences(g, model=gy94(omega=0.1, code=4), divergence=0.3, seed=1)
    assert by_four.coding_genes[fam] is True
    assert all(_frame_is_whole(seq, code=4) for _name, seq in _gene_copies(g, by_four))


# --- a coding gene keeps its frame --------------------------------------------------------------------

def test_no_coding_gene_loses_its_reading_frame(files):
    g = _genome(files, **BUSY)
    assert len(g.rearrangements) > 10 and len(g.events) > 30
    codon = simulate_sequences(g, model=gy94(omega=0.1), intergene_model=hky85(), divergence=0.3,
                               seed=1)
    copies = [(name, seq) for name, seq in _gene_copies(g, codon) if name in CODING]
    assert len(copies) > 60 and {name for name, _s in copies} == CODING
    assert all(_frame_is_whole(seq) for _name, seq in copies)
    assert any(seq != _original(name) for name, seq in copies)          # and they did evolve
    # the same genomes under a nucleotide model: most of those frames are gone
    plain = simulate_sequences(g, model=hky85(), intergene_model=hky85(), divergence=0.3, seed=1)
    broken = [not _frame_is_whole(seq) for name, seq in _gene_copies(g, plain) if name in CODING]
    assert sum(broken) > len(broken) / 2


def test_most_changes_keep_the_amino_acid(files):
    # the mark of the right frame on the right strand: read in any other, omega would not show
    g = _genome(files, **BUSY)
    r = simulate_sequences(g, model=gy94(omega=0.05), divergence=0.2, seed=1)
    for strand_genes in ({"A", "C", "D", "G"}, {"B"}):                   # plus strand, then minus
        same = changed = 0
        for name, seq in _gene_copies(g, r):
            if name in strand_genes:
                for a, b in zip(_codons(_original(name)), _codons(seq)):
                    changed += a != b
                    same += a != b and AMINO[11][a] == AMINO[11][b]
        assert changed > 50 and same / changed > 0.6


def test_a_stop_codon_and_the_bases_two_genes_share_are_held(files):
    g = _genome(files, **BUSY)
    r = simulate_sequences(g, model=gy94(omega=1.0), intergene_model=hky85(), divergence=0.5, seed=1)
    shared_cd = GENES["C"][1] - GENES["D"][0]
    assert shared_cd == 4
    seen = 0
    for name, seq in _gene_copies(g, r):
        original = _original(name)
        if name in CODING:
            assert seq[-3:] == original[-3:]                     # its own stop codon
        if name == "C":                                          # the last two codons touch D
            assert seq[-6:] == original[-6:]
        if name == "D":                                          # the first two codons touch C
            assert seq[:6] == original[:6]
        if name == "H":                                          # nested in G: all of it is shared
            assert seq == original
        if name == "G":                                          # and G's codons under H with it
            assert seq[60:93] == original[60:93] and seq != original
        seen += 1
    assert seen > 80


def test_a_gene_that_is_not_coding_and_the_spacer_still_evolve(files):
    g = _genome(files)
    r = simulate_sequences(g, model=gy94(omega=0.1), intergene_model=hky85(), divergence=0.3, seed=1)
    for name in ("E", "F"):
        versions = {seq for n, seq in _gene_copies(g, r) if n == name}
        assert len(versions) > 3                                  # plain DNA, at the gene's speed
    tips = [r.node_genomes[label] for label in r.genomes]
    assert len({chrom[:60] for genome in tips for chrom in genome.values()}) > 3     # leading spacer


def test_at_rate_zero_every_genome_is_the_dna_that_was_given(files):
    g = _genome(files, **BUSY)
    r = simulate_sequences(g, model=gy94(omega=0.1), substitution=0.0, seed=1)
    assert r.initial_genome[g.initial_genome.chromosomes[0].id] == DNA
    for name, seq in _gene_copies(g, r):
        assert seq == _original(name)


def test_the_parallel_engine_lays_the_blocks_out_the_same_way(files):
    g = _genome(files, **BUSY)
    kw = dict(model=gy94(omega=0.1), intergene_model=hky85(), divergence=0.2, seed=3)
    one, two = simulate_sequences(g, parallel=1, **kw), simulate_sequences(g, parallel=2, **kw)
    assert one.alignments == two.alignments and one.coding_genes == two.coding_genes
    assert all(_frame_is_whole(seq) for name, seq in _gene_copies(g, one) if name in CODING)
    # and with the spacer's model left to its default, which this engine used to fail on at any
    # model: the workers were handed the argument, None, in place of the model it resolves to
    assert simulate_sequences(g, model=hky85(), parallel=1, divergence=0.1, seed=3).alignments
    assert simulate_sequences(g, model=gy94(omega=0.1), parallel=1, divergence=0.1, seed=3).coding_genes


# --- a genome that was drawn rather than given ---------------------------------------------------------

def test_a_drawn_gene_of_whole_codons_is_coding_and_ends_in_a_stop():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1)
    g = simulate_genomes_nucleotide(sp, root_length=3000, genes=6, gene_length=300, inversion=2.0,
                                    inversion_extent=500, seed=1)
    r = simulate_sequences(g, model=gy94(omega=0.1), divergence=0.3, seed=1)
    assert r.summary()["genes"] == {"coding": 6, "not_coding": 0}
    names = g.complete_tree.labels()
    checked = 0
    for node, genome in g.node_genomes.items():
        for line in g._gff(names[node], genome).splitlines():
            cols = line.split("\t")
            if len(cols) == 9:
                chrom = int(cols[0].rsplit("_chr", 1)[1])
                seq = r.node_genomes[names[node]][chrom][int(cols[3]) - 1:int(cols[4])]
                seq = seq if cols[6] == "+" else _rc(seq)
                assert len(seq) == 300 and seq.endswith("TAA") and _frame_is_whole(seq)
                checked += 1
    assert checked > 50
    # 301 bases is not a whole number of codons, so there is no frame to draw
    odd = simulate_genomes_nucleotide(sp, root_length=3000, genes=6, gene_length=301, seed=1)
    assert simulate_sequences(odd, model=gy94(omega=0.1), seed=1).summary()["genes"] == \
        {"coding": 0, "not_coding": 6}


# --- refusals, and the command line --------------------------------------------------------------------

def test_a_genome_with_indels_is_refused(files):
    g = _genome(files, deletion=20.0, deletion_extent=5)
    assert g.indels > 0
    with pytest.raises(ValueError, match="indels"):
        simulate_sequences(g, model=gy94(omega=0.1), seed=1)
    simulate_sequences(g, model=hky85(), divergence=0.1, seed=1)        # a nucleotide model still runs


def test_the_command_line_runs_a_codon_model_on_a_nucleotide_genome(files, tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--death", "0", "--n-extant", "6",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run), "--resolution", "nucleotide", "--gff", str(files[0]),
                 "--fasta", str(files[1]), "--inversion", "2", "--inversion-extent", "300",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["sequences", str(run), "--model", "gy94", "--omega", "0.1",
                 "--intergene-model", "hky85", "--divergence", "0.2", "--seed", "1",
                 "--quiet"]) == 0
    summary = (run / "sequences" / "sequences_summary.json").read_text(encoding="utf-8")
    assert '"coding": 5' in summary and '"not_coding": 3' in summary
