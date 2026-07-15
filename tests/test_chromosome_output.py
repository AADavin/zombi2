"""Output of the chromosome tier (Stage 5): gene_order.tsv (the per-leaf layout — which chromosome
each gene sits on, and in what order) and karyotype_trace.tsv (the fission / fusion / origination /
loss genealogy).

Both are opt-in ordered-genome parts: they never appear for a single-chromosome run's default
output, so existing output folders are unchanged.
"""
from collections import defaultdict

import pytest

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    SharedRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.cli import main


def _multichrom(seed=4, n_chromosomes=3, transfer=0.3, **rate_kw):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=seed)
    rates = SharedRates(duplication=0.3, loss=0.2, transfer=transfer, origination=0.2,
                        inversion=0.2, transposition=0.2, **rate_kw)
    return simulate_genomes(tree, rates, initial_families=18, seed=seed,
                            genome_factory=lambda i: OrderedGenome(
                                i, extension=0.6, n_chromosomes=n_chromosomes, circular=True))


def _layout_rows(path):
    lines = path.read_text().splitlines()
    assert lines[0] == "species\tchromosome\tposition\tfamily\tgid\torientation"
    return [line.split("\t") for line in lines[1:]]


def test_layout_round_trips_the_karyotype(tmp_path):
    """gene_order.tsv reproduces every leaf's chromosomes exactly (family order + gids), keyed by
    chrom_id — the information the profiles alone do not carry."""
    g = _multichrom()
    out = tmp_path / "gen"
    g.write(out, include=["layout"])
    recon = defaultdict(lambda: defaultdict(list))
    for species, chrom, pos, family, gid, _orient in _layout_rows(out / "gene_order.tsv"):
        recon[species][int(chrom)].append((int(pos), family, gid))
    for leaf, genome in g.leaf_genomes.items():
        for chrom in genome.chromosomes.values():
            got = [(fam, gid) for _pos, fam, gid in sorted(recon[leaf.name][chrom.chrom_id])]
            want = [(gene.family, gene.gid) for gene in chrom.genes]
            assert got == want


def test_layout_shows_a_family_spanning_chromosomes_after_transfer(tmp_path):
    """A transferred copy can land on a different chromosome, so a family appears under >1 chrom_id
    in the layout for a single species — karyotype structure the profiles cannot show."""
    for seed in range(8):
        g = _multichrom(seed=seed, transfer=0.6)
        out = tmp_path / f"g{seed}"
        g.write(out, include=["layout"])
        fam_chroms = defaultdict(set)
        for species, chrom, _pos, family, _gid, _orient in _layout_rows(out / "gene_order.tsv"):
            fam_chroms[(species, family)].add(chrom)
        if any(len(cs) >= 2 for cs in fam_chroms.values()):
            return
    raise AssertionError("no family spanned >1 chromosome in the layout across eight seeds")


def test_karyotype_trace_records_the_genealogy(tmp_path):
    g = _multichrom(seed=3, chromosome_origination=0.3, chromosome_loss=0.2, fission=0.3, fusion=0.2)
    out = tmp_path / "gen"
    g.write(out, include=["karyotype"])
    lines = (out / "karyotype_trace.tsv").read_text().splitlines()
    assert lines[0] == "time\tevent\tbranch\tparents\tchildren"
    events = [line.split("\t")[1] for line in lines[1:]]
    assert "FI" in events and "CL" in events                       # fission + chromosome loss appear
    assert len(lines) - 1 == len(g.event_log.chromosome_records)    # one row per recorded event


def test_layout_and_karyotype_are_opt_in(tmp_path):
    """include=None (the library default) never writes the opt-in ordered parts."""
    g = _multichrom()
    out = tmp_path / "gen"
    g.write(out)                                            # include=None
    assert not (out / "gene_order.tsv").exists()
    assert not (out / "karyotype_trace.tsv").exists()
    assert (out / "profiles.tsv").exists()                  # the usual output is unaffected


def test_single_chromosome_output_folder_is_unchanged(tmp_path):
    """A single-chromosome ordered run writes no new files; layout is still available on request."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    rates = SharedRates(duplication=0.3, loss=0.2, origination=0.3, inversion=0.2)
    g = simulate_genomes(tree, rates, initial_families=12, seed=1,
                         genome_factory=lambda i: OrderedGenome(i, extension=0.6))
    out = tmp_path / "gen"
    g.write(out, include=["profiles", "trace"])
    assert not (out / "gene_order.tsv").exists()
    assert not (out / "karyotype_trace.tsv").exists()
    g.write(out, include=["layout"])                        # exposes gene order even at one chromosome
    assert (out / "gene_order.tsv").exists()


# --- CLI: automatic inclusion + the chromosome-tier rate flags ------------------------

def _species(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--birth", "1", "--death", "0.3", "--tips", "8", "--age", "3",
          "--seed", "1", "-o", str(sp)])
    return str(sp / "species_tree.nwk")


def test_cli_multichromosome_auto_writes_layout(tmp_path):
    """`--n-chromosomes > 1` auto-adds gene_order.tsv; with no chromosome-tier events, no trace."""
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "ordered", "--n-chromosomes", "4",
               "--dup", "0.2", "--loss", "0.2", "--orig", "0.3", "--inversion", "0.2",
               "--initial-families", "15", "--seed", "3", "--write", "profiles", "-o", str(out)])
    assert rc == 0
    assert (out / "gene_order.tsv").exists()
    assert not (out / "karyotype_trace.tsv").exists()


def test_cli_fission_auto_writes_karyotype_trace(tmp_path):
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "ordered", "--n-chromosomes", "2",
               "--dup", "0.2", "--loss", "0.2", "--orig", "0.3",
               "--fission", "0.4", "--chromosome-origination", "0.3", "--chromosome-loss", "0.2",
               "--initial-families", "15", "--seed", "3", "--write", "profiles", "-o", str(out)])
    assert rc == 0
    assert (out / "gene_order.tsv").exists() and (out / "karyotype_trace.tsv").exists()
    assert len((out / "karyotype_trace.tsv").read_text().splitlines()) > 1  # events recorded


def test_cli_single_chromosome_writes_no_karyotype_files(tmp_path):
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "ordered",
               "--dup", "0.2", "--loss", "0.2", "--orig", "0.3", "--inversion", "0.2",
               "--initial-families", "15", "--seed", "3", "-o", str(out)])
    assert rc == 0
    assert not (out / "gene_order.tsv").exists()
    assert not (out / "karyotype_trace.tsv").exists()


def test_cli_rejects_chromosome_tier_rates_without_ordered(tmp_path):
    tree = _species(tmp_path)
    with pytest.raises(SystemExit):  # default genome model is unordered -> no chromosome tier
        main(["genomes", "--tree", tree, "--fission", "0.1", "--seed", "1", "-o", str(tmp_path / "x")])


# --- CLI: the nucleotide model shares the unified chromosome tier ---------------------

def _nuc_chrom_counts(path):
    """species -> number of distinct chromosomes, from a chromosomes.tsv layout file."""
    lines = path.read_text().splitlines()
    assert lines[0] == "species\tchromosome\ttopology\tposition\tsource\tstart\tend\tstrand"
    by = defaultdict(set)
    for line in lines[1:]:
        species, chrom = line.split("\t")[:2]
        by[species].add(chrom)
    return {k: len(v) for k, v in by.items()}


def test_cli_nucleotide_n_chromosomes_writes_layout(tmp_path):
    """`--genome-model nucleotide --n-chromosomes N` seeds N chromosomes and writes chromosomes.tsv
    (the unified flag; no chromosome-tier events -> no trace)."""
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--n-chromosomes", "3",
               "--inversion", "0.01", "--root-length", "200", "--seed", "3",
               "--write", "trees", "-o", str(out)])
    assert rc == 0
    assert set(_nuc_chrom_counts(out / "chromosomes.tsv").values()) == {3}  # every leaf has 3
    assert not (out / "karyotype_trace.tsv").exists() or \
        len((out / "karyotype_trace.tsv").read_text().splitlines()) == 1     # header only


def test_cli_nucleotide_fission_writes_karyotype_trace(tmp_path):
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide",
               "--inversion", "0.01", "--loss", "0.005", "--fission", "0.4",
               "--chromosome-origination", "0.2", "--chromosome-loss", "0.1",
               "--root-length", "200", "--seed", "3", "--write", "profiles", "-o", str(out)])
    assert rc == 0                                       # tier rates force the Python engine
    assert (out / "chromosomes.tsv").exists() and (out / "karyotype_trace.tsv").exists()
    events = [ln.split("\t")[1] for ln in (out / "karyotype_trace.tsv").read_text().splitlines()[1:]]
    assert "FI" in events                                # fission recorded


def test_cli_nucleotide_single_chromosome_writes_no_karyotype_files(tmp_path):
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide",
               "--inversion", "0.01", "--root-length", "200", "--seed", "3",
               "--write", "trees", "-o", str(out)])
    assert rc == 0
    assert not (out / "chromosomes.tsv").exists()
    assert not (out / "karyotype_trace.tsv").exists()


def test_cli_nucleotide_initial_chromosomes_is_a_deprecated_alias(tmp_path):
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--initial-chromosomes", "2",
               "--inversion", "0.01", "--root-length", "200", "--seed", "3",
               "--write", "trees", "-o", str(out)])
    assert rc == 0
    assert set(_nuc_chrom_counts(out / "chromosomes.tsv").values()) == {2}


def _nuc_topologies(path):
    """set of topologies present in a chromosomes.tsv layout file (asserting per-chromosome
    consistency: every chromosome's rows share one topology)."""
    rows = [ln.split("\t") for ln in path.read_text().splitlines()[1:]]
    per_chrom = defaultdict(set)
    for r in rows:
        per_chrom[(r[0], r[1])].add(r[2])                    # (species, chromosome) -> {topology}
    assert all(len(v) == 1 for v in per_chrom.values())      # a chromosome has one topology
    return {t for v in per_chrom.values() for t in v}


def test_cli_nucleotide_linear_chromosomes(tmp_path):
    """`--linear-chromosomes` makes every nucleotide chromosome linear (the topology column says so)."""
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--n-chromosomes", "2",
               "--linear-chromosomes", "--inversion", "0.005", "--root-length", "250",
               "--write", "trees", "--seed", "3", "-o", str(out)])
    assert rc == 0
    assert _nuc_topologies(out / "chromosomes.tsv") == {"linear"}


def test_cli_multisequence_gff_keeps_per_replicon_topology(tmp_path):
    """A multi-sequence GFF seeds each replicon with its own `Is_circular` topology — a linear
    chromosome plus a circular plasmid become a mixed-topology genome."""
    tree = _species(tmp_path)
    gff = tmp_path / "genome.gff"
    gff.write_text("##gff-version 3\n"
                   "##sequence-region chrom 1 400\n"
                   "chrom\tx\tregion\t1\t400\t.\t+\t.\tID=chrom;Is_circular=false\n"
                   "chrom\tx\tgene\t30\t90\t.\t+\t.\tID=gene-c1;locus_tag=c1\n"
                   "chrom\tx\tgene\t200\t280\t.\t-\t.\tID=gene-c2;locus_tag=c2\n"
                   "##sequence-region plasmid 1 200\n"
                   "plasmid\tx\tregion\t1\t200\t.\t+\t.\tID=plasmid;Is_circular=true\n"
                   "plasmid\tx\tgene\t20\t80\t.\t+\t.\tID=gene-p1;locus_tag=p1\n")
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--gff", str(gff),
               "--inversion", "0.003", "--write", "trees", "--seed", "3", "-o", str(out)])
    assert rc == 0
    assert _nuc_topologies(out / "chromosomes.tsv") == {"linear", "circular"}   # mixed genome


_MULTI_GFF = """\
##gff-version 3
##sequence-region chr1 1 600
chr1\tx\tgene\t50\t140\t.\t+\t.\tID=gene-A1;locus_tag=A1
chr1\tx\tgene\t200\t320\t.\t-\t.\tID=gene-A2;locus_tag=A2
chr1\tx\tgene\t400\t500\t.\t+\t.\tID=gene-A3;locus_tag=A3
##sequence-region plasmid 1 400
plasmid\tx\tgene\t30\t120\t.\t+\t.\tID=gene-B1;locus_tag=B1
plasmid\tx\tgene\t200\t300\t.\t-\t.\tID=gene-B2;locus_tag=B2
"""


def test_cli_multisequence_gff_seeds_one_chromosome_per_sequence(tmp_path):
    """A GFF with several sequences (a chromosome + a plasmid) seeds one chromosome per sequence —
    every gene is kept (the old single-sequence path silently dropped the plasmid)."""
    tree = _species(tmp_path)
    gff = tmp_path / "genome.gff"
    gff.write_text(_MULTI_GFF)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--gff", str(gff),
               "--inversion", "0.002", "--loss", "0.001", "--seed", "3",
               "--write", "trees", "-o", str(out)])
    assert rc == 0
    assert set(_nuc_chrom_counts(out / "chromosomes.tsv").values()) == {2}   # both replicons present
    genes = {ln.split("\t")[0] for ln in (out / "genes.tsv").read_text().splitlines()[1:]}
    assert {"A1", "A2", "A3", "B1", "B2"} <= genes                          # nothing dropped


def test_cli_gff_seqid_still_picks_one_sequence(tmp_path):
    tree = _species(tmp_path)
    gff = tmp_path / "genome.gff"
    gff.write_text(_MULTI_GFF)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--gff", str(gff),
               "--gff-seqid", "plasmid", "--inversion", "0.002", "--seed", "3",
               "--write", "trees", "-o", str(out)])
    assert rc == 0
    assert not (out / "chromosomes.tsv").exists()                            # single chromosome
    genes = {ln.split("\t")[0] for ln in (out / "genes.tsv").read_text().splitlines()[1:]}
    assert genes == {"B1", "B2"}                                             # only the picked sequence


def test_cli_multichromosome_ancestral_reproduces_each_replicon(tmp_path):
    """`--write ancestral` on a multi-sequence GFF writes one FASTA record per chromosome at every
    node; with a matching multi-record --genome-fasta and zero divergence, the root reproduces each
    replicon exactly."""
    import gzip
    tree = _species(tmp_path)
    gff = tmp_path / "genome.gff"
    gff.write_text(_MULTI_GFF)                                # chr1 (600 bp) + plasmid (400 bp)
    chr1, plasmid = ("ACGT" * 150), ("TGCA" * 100)           # 600 / 400 bp, matching the GFF
    (tmp_path / "g.fasta").write_text(f">chr1\n{chr1}\n>plasmid\n{plasmid}\n")
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--gff", str(gff),
               "--genome-fasta", str(tmp_path / "g.fasta"), "--inversion", "0.002", "--loss", "0.001",
               "--subst-rate", "0.0", "--write", "ancestral", "--seed", "3", "-o", str(out)])
    assert rc == 0
    recs = {}
    with gzip.open(out / "genomes" / "root.fasta.gz", "rt") as fh:
        name = None
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]; recs[name] = ""
            else:
                recs[name] += line
    assert len(recs) == 2                                     # one record per replicon
    assert set(recs.values()) == {chr1, plasmid}             # each replicon reproduced exactly
    header = (out / "architecture" / "root.tsv").read_text().splitlines()[0]
    assert header.split("\t")[0] == "chromosome"             # architecture keeps replicons apart


def test_cli_multichromosome_bed_is_per_replicon(tmp_path):
    """`--write bed` on a multi-sequence GFF writes one BED contig per replicon, coordinates
    restarting per chromosome, named to line up with each chromosome's FASTA record."""
    tree = _species(tmp_path)
    gff = tmp_path / "genome.gff"
    gff.write_text(_MULTI_GFF)                                # chr1 (A1/A2/A3) + plasmid (B1/B2)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--gff", str(gff),
               "--inversion", "0.002", "--write", "bed", "--seed", "3", "-o", str(out)])
    assert rc == 0
    # genes.bed: the seed annotation, one contig per replicon named by its GFF seqid
    rows = [ln.split("\t") for ln in (out / "genes.bed").read_text().splitlines()]
    by_contig = {}
    for chrom, start, _end, name, *_ in rows:
        by_contig.setdefault(chrom, []).append((name, int(start)))
    assert set(by_contig) == {"chr1", "plasmid"}             # the input sequence names
    assert {n for n, _ in by_contig["chr1"]} == {"A1", "A2", "A3"}
    assert {n for n, _ in by_contig["plasmid"]} == {"B1", "B2"}
    # the plasmid's coordinates restart at 0 (B1 at 29 on its own contig, not offset past chr1)
    assert min(s for _, s in by_contig["plasmid"]) == 29
    # BED/root.bed contigs match the ancestral FASTA record naming (<node>_chr<id>)
    node_contigs = {ln.split("\t")[0] for ln in (out / "bed" / "root.bed").read_text().splitlines()}
    assert node_contigs == {"root_chr0", "root_chr1"}


# --- CLI: translocation (nucleotide) ---------------------------------------------------

def test_cli_translocation_moves_material_across_chromosomes(tmp_path):
    """`--translocation` on a multi-chromosome nucleotide genome carries material onto another
    chromosome, so some source occupies more than one chromosome in the layout."""
    from collections import defaultdict
    tree = _species(tmp_path)
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", tree, "--genome-model", "nucleotide", "--n-chromosomes", "2",
               "--inversion", "0.005", "--translocation", "0.03", "--root-length", "250",
               "--write", "trees", "--seed", "5", "-o", str(out)])
    assert rc == 0
    # chromosomes.tsv columns: species chromosome topology position source start end strand
    by = defaultdict(lambda: defaultdict(set))          # species -> source -> {chromosome ids}
    for line in (out / "chromosomes.tsv").read_text().splitlines()[1:]:
        parts = line.split("\t")
        species, chrom, source = parts[0], parts[1], parts[4]
        by[species][source].add(chrom)
    # a source appearing on >1 chromosome in some species is the signature of a translocation
    crossed = any(len(cids) > 1 for sp in by.values() for cids in sp.values())
    assert crossed


def test_cli_translocation_accepted_for_ordered(tmp_path):
    """The ordered model now supports translocation too (cross-chromosome move); accepting it with a
    single chromosome is a harmless no-op, exactly as the nucleotide model treats it."""
    tree = _species(tmp_path)
    rc = main(["genomes", "--tree", tree, "--genome-model", "ordered", "--translocation", "0.1",
               "--seed", "1", "-o", str(tmp_path / "x")])
    assert rc == 0


def test_cli_translocation_rejected_for_unordered(tmp_path):
    tree = _species(tmp_path)
    with pytest.raises(SystemExit):  # the unordered model has no chromosomes to move between
        main(["genomes", "--tree", tree, "--genome-model", "unordered", "--translocation", "0.1",
              "--seed", "1", "-o", str(tmp_path / "x")])
