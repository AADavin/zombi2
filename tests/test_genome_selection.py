"""P2.5 GFF-driven genome selection: partition a real annotated genome into CDS (codon selection)
and non-coding (neutral), respecting strand + frame, and reassemble down a tree. Torch-free."""
from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

import zombi2.experimental as ex
from zombi2.experimental.codon_selection import translate
from zombi2.experimental.genome_selection import CDS, GenomeSelection, read_cds_gff
from zombi2.experimental.selection import Critic
from zombi2.sequences.models import AMINO_ACIDS, reverse_complement

_AA = {a: i for i, a in enumerate(AMINO_ACIDS)}


class _N:
    def __init__(self, gid, children=()):
        self.gid = gid
        self.children = list(children)


class _PreferCritic(Critic):
    """Length-adaptive: every position prefers residue ``aa`` (so any CDS is driven toward all-``aa``)."""

    def __init__(self, aa="W", hi=0.95):
        self.aa, self.hi = aa, hi

    def profile(self, seq):
        p = np.full((len(seq), 20), (1.0 - self.hi) / 19.0)
        p[:, _AA[self.aa]] = self.hi
        return p


def _lineage(branch):
    root = _N("root", [_N("tip")])
    return root, {root: 0.0, root.children[0]: branch}


# --------------------------------------------------------------------------- #
# GFF reader + lifecycle
# --------------------------------------------------------------------------- #
def test_read_cds_gff(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n"
        "chr1\t.\tCDS\t7\t24\t.\t+\t0\tID=cds1;locus_tag=g1\n"
        "chr1\t.\tCDS\t31\t48\t.\t-\t0\tID=cds2\n"
        "chr2\t.\tCDS\t1\t9\t.\t+\t0\tID=other\n"
    )
    got = read_cds_gff(gff, seqid="chr1")
    assert got == [CDS(6, 24, 1, "g1"), CDS(30, 48, -1, "cds2")]     # 1-based incl -> 0-based half-open
    with pytest.raises(ValueError, match="multiple sequences"):      # multi-contig GFF needs a seqid
        read_cds_gff(gff)


def test_read_cds_gff_rejects_multiexon_and_bad_strand(tmp_path):
    multi = tmp_path / "m.gff"
    multi.write_text("chr1\t.\tCDS\t1\t9\t.\t+\t0\tParent=g1\n"
                     "chr1\t.\tCDS\t20\t28\t.\t+\t0\tParent=g1\n")    # two exons of one gene
    with pytest.raises(ValueError, match="multi-exon"):
        read_cds_gff(multi)
    bad = tmp_path / "b.gff"
    bad.write_text("chr1\t.\tCDS\t1\t9\t.\t.\t0\tID=c\n")             # strand '.'
    with pytest.raises(ValueError, match="strand"):
        read_cds_gff(bad)


def test_genome_selection_is_experimental():
    ex._warned.discard("GenomeSelection")
    with pytest.warns(ex.ExperimentalWarning, match="GenomeSelection"):
        GenomeSelection(_PreferCritic())


# --------------------------------------------------------------------------- #
# structure: length preserved, root reproduces input, partition correct
# --------------------------------------------------------------------------- #
def _toy_genome():
    #  [0,6) nc | [6,24) CDS+ | [24,30) nc | [30,48) CDS- | [48,60) nc
    fwd, rev = "GCT" * 6, reverse_complement("GCT" * 6)              # alanine codons
    genome = "ACGTAC" + fwd + "TTAACC" + rev + "GGGCCCAAATTT"
    assert len(genome) == 60
    return genome, [CDS(6, 24, 1, "f"), CDS(30, 48, -1, "r")]


def test_length_preserved_and_root_reproduces_input():
    genome, cds = _toy_genome()
    root, subst = _lineage(4.0)
    out = GenomeSelection(_PreferCritic(), beta=2.0).evolve_genome(
        root, subst, genome, cds, rng=np.random.default_rng(0))
    assert set(out) == {"root", "tip"}
    assert all(len(g) == len(genome) for g in out.values())          # length preserved everywhere
    assert out["root"] == genome                                     # zero-length root branch = input
    assert translate(reverse_complement(out["root"][30:48])) == "AAAAAA"  # minus-strand CDS read correctly


def test_forward_and_reverse_cds_recover_under_selection():
    genome, cds = _toy_genome()
    root, subst = _lineage(120.0)                                   # long branch (Trp has one isolated codon)
    tip = GenomeSelection(_PreferCritic("W"), beta=6.0).evolve_genome(
        root, subst, genome, cds, rng=np.random.default_rng(0))["tip"]
    fwd = translate(tip[6:24])
    rev = translate(reverse_complement(tip[30:48]))                 # reverse CDS read off the - strand
    assert fwd.count("W") >= 5 and rev.count("W") >= 5              # both strands driven to Trp
    assert "*" not in fwd                                           # no stop introduced in a CDS


def test_terminal_stop_codon_is_kept_fixed():
    coding = "GCT" * 4 + "TAA"                                       # 4 alanine codons + a stop
    genome = "AA" + coding + "AA"
    root, subst = _lineage(120.0)
    tip = GenomeSelection(_PreferCritic("W"), beta=6.0).evolve_genome(
        root, subst, genome, [CDS(2, 2 + len(coding), 1, "c")], rng=np.random.default_rng(1))["tip"]
    prot = translate(tip[2:2 + len(coding)])
    assert prot[:4].count("W") >= 3 and prot[4] == "*"              # sense codons evolve, stop stays a stop
    assert tip[2 + len(coding) - 3:2 + len(coding)] == "TAA"         # the stop codon is untouched


# --------------------------------------------------------------------------- #
# non-coding, edges, determinism, validation
# --------------------------------------------------------------------------- #
def test_all_intergenic_and_all_cds_edges():
    root, subst = _lineage(2.0)
    gsel = GenomeSelection(_PreferCritic(), beta=2.0)
    nc = gsel.evolve_genome(root, subst, "ACGTACGTAC", [], rng=np.random.default_rng(0))
    assert len(nc["tip"]) == 10                                      # no CDS -> pure neutral, length kept
    allc = gsel.evolve_genome(root, subst, "GCT" * 4, [CDS(0, 12, 1, "c")], rng=np.random.default_rng(0))
    assert len(allc["tip"]) == 12                                    # no intergenic


def test_deterministic_given_seed():
    genome, cds = _toy_genome()
    root, subst = _lineage(5.0)
    gsel = GenomeSelection(_PreferCritic(), beta=3.0)
    a = gsel.evolve_genome(root, subst, genome, cds, rng=np.random.default_rng(7))
    b = gsel.evolve_genome(root, subst, genome, cds, rng=np.random.default_rng(7))
    assert a == b


@pytest.mark.parametrize("bad,match", [
    ([CDS(6, 24, 1, "a"), CDS(20, 30, 1, "b")], "overlaps"),         # overlapping
    ([CDS(6, 25, 1, "a")], "multiple of 3"),                        # frame
    ([CDS(6, 90, 1, "a")], "out of"),                               # out of bounds
    ([CDS(6, 24, 0, "a")], "strand"),                               # bad strand
])
def test_invalid_cds_are_rejected(bad, match):
    genome, _ = _toy_genome()
    root, subst = _lineage(1.0)
    with pytest.raises(ValueError, match=match):
        GenomeSelection(_PreferCritic()).evolve_genome(root, subst, genome, bad,
                                                       rng=np.random.default_rng(0))


def test_bifurcating_tree_intergenic_drifts_and_records_ancestors():
    # the neutral intergenic path must ACTUALLY evolve (not identity) and record every internal node
    t1, t2 = _N("t1"), _N("t2")
    root = _N("root", [_N("internal", [t1, t2])])
    internal = root.children[0]
    subst = {root: 0.0, internal: 2.0, t1: 2.0, t2: 2.0}
    genome = "ACGT" * 6                                              # 24 nt, all intergenic
    out = GenomeSelection(_PreferCritic()).evolve_genome(
        root, subst, genome, [], rng=np.random.default_rng(0))
    assert set(out) == {"root", "internal", "t1", "t2"}             # every node recorded
    assert all(len(g) == 24 for g in out.values())                  # length preserved everywhere
    assert out["root"] == genome                                    # root reproduces input
    assert out["t1"] != genome and out["t2"] != genome             # neutral drift actually happened


def test_ambiguity_preserved_in_non_coding_and_rejected_in_cds():
    # N/IUPAC in non-coding is frozen (root reproduces it, every node keeps it); inside a CDS it errors
    genome = "ACNNGT" + "GCT" * 4 + "ACRGT"                          # N in leading nc, R in trailing nc
    root, subst = _lineage(30.0)
    out = GenomeSelection(_PreferCritic("W"), beta=6.0).evolve_genome(
        root, subst, genome, [CDS(6, 18, 1, "c")], rng=np.random.default_rng(0))
    assert out["root"] == genome                                    # ambiguity codes preserved on root
    assert out["tip"][2:4] == "NN" and out["tip"][20] == "R"        # frozen at the tip too
    with pytest.raises(ValueError, match="not ACGT"):
        GenomeSelection(_PreferCritic()).evolve_genome(
            root, subst, "AAGCTGNTGCTAA", [CDS(2, 11, 1, "c")], rng=np.random.default_rng(0))


def test_empty_genome_and_phased_cds_are_rejected():
    genome, _ = _toy_genome()
    root, subst = _lineage(1.0)
    gsel = GenomeSelection(_PreferCritic())
    with pytest.raises(ValueError, match="empty"):
        gsel.evolve_genome(root, subst, "", [], rng=np.random.default_rng(0))
    with pytest.raises(ValueError, match="phase"):
        gsel.evolve_genome(root, subst, genome, [CDS(6, 24, 1, "c", phase=1)],
                           rng=np.random.default_rng(0))


def test_minus_strand_terminal_stop_is_kept_fixed():
    coding = "GCT" * 4 + "TAA"                                       # coding ends in a stop
    interval = reverse_complement(coding)                           # its bytes on the - strand genome
    genome = "AA" + interval + "AA"
    root, subst = _lineage(120.0)
    tip = GenomeSelection(_PreferCritic("W"), beta=6.0).evolve_genome(
        root, subst, genome, [CDS(2, 2 + len(interval), -1, "c")], rng=np.random.default_rng(2))["tip"]
    prot = translate(reverse_complement(tip[2:2 + len(interval)]))
    assert prot[-1] == "*" and prot[:4].count("W") >= 3            # stop fixed, sense codons evolve
    assert tip[2:5] == interval[:3]                                 # the stop's bytes (interval start) unchanged


def test_module_has_no_top_level_ml_imports():
    from zombi2.experimental import genome_selection
    tree = ast.parse(inspect.getsource(genome_selection))
    top: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            top.add(node.module.split(".")[0])
    assert "torch" not in top and "esm" not in top, top


def test_exports_stay_in_the_experimental_namespace():
    import zombi2
    for name in ("CDS", "GenomeSelection", "read_cds_gff"):
        assert name in ex.__all__ and hasattr(ex, name)
        assert not hasattr(zombi2, name), f"{name} leaked into the top-level zombi2 namespace"
