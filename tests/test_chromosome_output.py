"""Output of the chromosome tier (Stage 5): Gene_order.tsv (the per-leaf layout — which chromosome
each gene sits on, and in what order) and Karyotype_trace.tsv (the fission / fusion / origination /
loss genealogy).

Both are opt-in ordered-genome parts: they never appear for a single-chromosome run's default
output, so existing output folders are unchanged.
"""
from collections import defaultdict

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    SharedRates,
    simulate_genomes,
    simulate_species_tree,
)


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
    """Gene_order.tsv reproduces every leaf's chromosomes exactly (family order + gids), keyed by
    chrom_id — the information the profiles alone do not carry."""
    g = _multichrom()
    out = tmp_path / "gen"
    g.write(out, include=["layout"])
    recon = defaultdict(lambda: defaultdict(list))
    for species, chrom, pos, family, gid, _orient in _layout_rows(out / "Gene_order.tsv"):
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
        for species, chrom, _pos, family, _gid, _orient in _layout_rows(out / "Gene_order.tsv"):
            fam_chroms[(species, family)].add(chrom)
        if any(len(cs) >= 2 for cs in fam_chroms.values()):
            return
    raise AssertionError("no family spanned >1 chromosome in the layout across eight seeds")


def test_karyotype_trace_records_the_genealogy(tmp_path):
    g = _multichrom(seed=3, chromosome_origination=0.3, chromosome_loss=0.2, fission=0.3, fusion=0.2)
    out = tmp_path / "gen"
    g.write(out, include=["karyotype"])
    lines = (out / "Karyotype_trace.tsv").read_text().splitlines()
    assert lines[0] == "time\tevent\tbranch\tparents\tchildren"
    events = [line.split("\t")[1] for line in lines[1:]]
    assert "FI" in events and "CL" in events                       # fission + chromosome loss appear
    assert len(lines) - 1 == len(g.event_log.chromosome_records)    # one row per recorded event


def test_layout_and_karyotype_are_opt_in(tmp_path):
    """include=None (the library default) never writes the opt-in ordered parts."""
    g = _multichrom()
    out = tmp_path / "gen"
    g.write(out)                                            # include=None
    assert not (out / "Gene_order.tsv").exists()
    assert not (out / "Karyotype_trace.tsv").exists()
    assert (out / "Profiles.tsv").exists()                  # the usual output is unaffected


def test_single_chromosome_output_folder_is_unchanged(tmp_path):
    """A single-chromosome ordered run writes no new files; layout is still available on request."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    rates = SharedRates(duplication=0.3, loss=0.2, origination=0.3, inversion=0.2)
    g = simulate_genomes(tree, rates, initial_families=12, seed=1,
                         genome_factory=lambda i: OrderedGenome(i, extension=0.6))
    out = tmp_path / "gen"
    g.write(out, include=["profiles", "trace"])
    assert not (out / "Gene_order.tsv").exists()
    assert not (out / "Karyotype_trace.tsv").exists()
    g.write(out, include=["layout"])                        # exposes gene order even at one chromosome
    assert (out / "Gene_order.tsv").exists()
