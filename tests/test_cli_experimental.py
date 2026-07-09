"""The `zombi2 experimental selection` CLI. Argument-parsing / validation tests are torch-free; the
end-to-end run needs the optional ESM2 deps (guarded)."""
from __future__ import annotations

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.cli import main


def _write_inputs(tmp_path, *, coding_len_codons=19, cds_start0=9):
    """A tiny species tree + genome FASTA + single-CDS GFF3 in tmp_path; returns their paths."""
    from zombi2.experimental.codon_selection import SENSE_CODONS
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=4, age=1.0, seed=1)
    tpath = tmp_path / "species_tree.nwk"
    tpath.write_text(tree.to_newick())

    rng = np.random.default_rng(0)
    length = 120
    g = list("".join("ACGT"[i] for i in rng.integers(4, size=length)))
    coding = "".join(SENSE_CODONS[i] for i in rng.integers(len(SENSE_CODONS), size=coding_len_codons))
    end0 = cds_start0 + 3 * coding_len_codons
    g[cds_start0:end0] = list(coding)
    genome = "".join(g)
    fpath = tmp_path / "genome.fna"
    fpath.write_text(f">seq1\n{genome}\n")

    gff = tmp_path / "genome.gff"                             # GFF is 1-based inclusive
    gff.write_text("##gff-version 3\n"
                   f"seq1\t.\tCDS\t{cds_start0 + 1}\t{end0}\t.\t+\t0\tID=cds1;locus_tag=g1\n")
    return tpath, fpath, gff, genome


# --------------------------------------------------------------------------- #
# argument parsing / validation (torch-free)
# --------------------------------------------------------------------------- #
def test_help_exits_zero():
    with pytest.raises(SystemExit) as e:
        main(["experimental", "selection", "-h"])
    assert e.value.code == 0


def test_beta_and_target_dnds_are_mutually_exclusive(tmp_path):
    t, f, gff, _ = _write_inputs(tmp_path)
    with pytest.raises(SystemExit):                            # argparse mutual-exclusion -> exit 2
        main(["experimental", "selection", "-t", str(t), "--gff", str(gff),
              "--genome-fasta", str(f), "-o", str(tmp_path / "o"),
              "--beta", "1.0", "--target-dnds", "0.3"])


def test_missing_required_gff_errors(tmp_path):
    t, f, _, _ = _write_inputs(tmp_path)
    with pytest.raises(SystemExit):
        main(["experimental", "selection", "-t", str(t), "--genome-fasta", str(f),
              "-o", str(tmp_path / "o")])


def test_bad_subst_model_is_rejected(tmp_path):
    t, f, gff, _ = _write_inputs(tmp_path)
    with pytest.raises(SystemExit):                            # lg is a protein model, not a choice
        main(["experimental", "selection", "-t", str(t), "--gff", str(gff),
              "--genome-fasta", str(f), "-o", str(tmp_path / "o"), "--subst-model", "lg"])


def test_gff_fasta_seqid_mismatch_errors_cleanly(tmp_path, capsys):
    # a FASTA whose sequence is named differently from the GFF's CDS contig must NOT be silently paired.
    # needs the optional deps present, else the earlier dependency probe reports first (as on a bare CI).
    pytest.importorskip("torch")
    pytest.importorskip("esm")
    t, f, gff, genome = _write_inputs(tmp_path)
    f.write_text(f">chrOTHER\n{genome}\n")                     # GFF is on 'seq1', FASTA is 'chrOTHER'
    rc = main(["experimental", "selection", "-t", str(t), "--gff", str(gff),
               "--genome-fasta", str(f), "--beta", "1.0", "-o", str(tmp_path / "o")])
    assert rc == 1                                             # clean one-line error, not a traceback
    assert "seq1" in capsys.readouterr().err


class _PeakedCritic:
    """Torch-free stub critic: a profile peaked on each protein's own residues (any length)."""

    def profile(self, protein):
        from zombi2.sequences.models import AMINO_ACIDS
        idx = {a: i for i, a in enumerate(AMINO_ACIDS)}
        P = np.full((len(protein), 20), (1 - 0.9) / 19)
        for i, a in enumerate(protein):
            P[i, idx[a]] = 0.9
        return P


def test_calibrate_beta_genomewide_hits_the_target(tmp_path):
    # the fix for the OOM concat bug: calibrate one beta over MANY proteins (profiled independently),
    # matching a length-weighted-mean dN/dS -- verified against the same weighted mean.
    from zombi2.cli import _calibrate_beta_genomewide
    from zombi2.experimental.codon_selection import CodonSelection
    from zombi2.experimental.selection import FixedProfileCritic
    from zombi2.experimental.codon_selection import translate
    from zombi2.experimental.codon_selection import SENSE_CODONS

    rng = np.random.default_rng(3)
    proteins = [translate("".join(SENSE_CODONS[i] for i in rng.integers(len(SENSE_CODONS), size=n)))
                for n in (18, 25, 12)]
    crit = _PeakedCritic()
    beta = _calibrate_beta_genomewide(crit, proteins, 0.3, None)
    total = sum(len(p) for p in proteins)
    got = sum(len(p) * CodonSelection(FixedProfileCritic(crit.profile(p)), beta=beta).dnds(p)
              for p in proteins) / total
    assert abs(got - 0.3) < 5e-3


# --------------------------------------------------------------------------- #
# end-to-end (needs the optional ESM2 deps)
# --------------------------------------------------------------------------- #
def test_end_to_end_writes_genomes_and_reproduces_root(tmp_path):
    pytest.importorskip("torch")
    pytest.importorskip("esm")
    from zombi2.sequences.models import read_fasta
    t, f, gff, genome = _write_inputs(tmp_path)
    out = tmp_path / "out"
    rc = main(["experimental", "selection", "-t", str(t), "--gff", str(gff),
               "--genome-fasta", str(f), "--beta", "2.0", "--subst-rate", "0.5",
               "--seed", "7", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists()
    assert (out / "Genomes" / "root.fasta.gz").exists()
    assert (out / "Gene_alignments" / "g1.fasta").exists()
    # the root node reproduces the input genome (no incoming branch)
    root_seq = next(iter(read_fasta(str(out / "Genomes" / "root.fasta.gz")).values()))
    assert root_seq == genome
    # the report records the single gene block as selected
    report = (out / "Selection_report.tsv").read_text()
    assert "n_selected\t1" in report and "n_gene_blocks\t1" in report
