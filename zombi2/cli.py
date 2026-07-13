"""Command-line interface for ZOMBI2 (``species`` / ``genomes`` / ``trait`` / ``sequence`` /
``coevolve`` / ``tools`` subcommands)."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time

import numpy as np

from . import __version__

from zombi2.traits.biogeography import DEC, simulate_biogeography
from zombi2.species.ghosts import add_ghost_lineages
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes
from zombi2.genomes.profiles import ProfileMatrix
from zombi2.distributions import LogNormal
from zombi2.sequences.clocks import (
    AutocorrelatedLogNormalClock, CIRClock, RateVariation, StrictClock,
    UncorrelatedGammaClock, UncorrelatedLogNormalClock, WhiteNoiseClock,
)
from zombi2.genomes.genome import OrderedGenome
from zombi2.genomes.rates import BranchRates, FamilySampledRates, PerGenomeRates, SharedRates
from zombi2.genomes.conversion import ConversionModel
from zombi2.genomes.read_rates import read_branch_rates, read_family_rates
from zombi2.sequences.evolution import SequenceEvolution
from zombi2.genomes.simulation import Genomes, simulate_genomes
from zombi2.species.model import (
    BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent, EpisodicBirthDeath,
)
from zombi2.species.sim import simulate_species_tree
from zombi2.coevolve.sse import BiSSE, MuSSE, QuaSSE, HiSSE, simulate_sse
from zombi2.coevolve.gene_diversification import (
    GeneDiversification, simulate_gene_diversification, simulate_co_diversification,
)
from zombi2.coevolve.cladogenetic_genome import (
    CladogeneticGenome, simulate_cladogenetic_genome, _branch_count_and_length,
)
from zombi2.coevolve.gene_conditioned_trait import GeneConditionedTrait, simulate_gene_conditioned_trait
from zombi2.coevolve.trait_coupling import TraitGeneCoupling, simulate_trait_linked_genomes
from zombi2.coevolve.trait_gene_feedback import TraitGeneFeedback, simulate_trait_gene_feedback
from zombi2.traits.models import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst, Mk, ThresholdModel, TraitResult,
    Cladogenesis, simulate_traits,
)
from zombi2.genomes.transfers import TransferModel
from zombi2.tree import Tree, prune, read_newick

_DESCRIPTION = """\
Simulate each level on its own, or couple them into joint models. Run
'zombi2 <command> -h' for a command's options, grouped by model.

Species trees
  species              simulate a dated species tree

Gene families & sequences
  genomes              evolve gene families along a species tree (Newick)
  sequence             simulate DNA/protein alignments along a genomes run's gene trees

Traits & coevolution
  trait                evolve a phenotypic trait along a given species tree
  coevolve             co-evolve coupled processes (--couple driver:target)

Analysis tools
  tools                compute on ZOMBI2 outputs (reconcile, treedist, recon-accuracy, red, parse, export)

Experimental (unstable, opt-in)
  experimental         not-yet-validated models (selection: ESM2 dN/dS; ils: multispecies coalescent)
"""


# ── house style: an IQ-TREE-like grouped, sectioned help ────────────────────────────
_BOLD, _RESET = "\033[1m", "\033[0m"


def _use_color() -> bool:
    """Bold section headers only for an interactive terminal (never when piped/redirected,
    under NO_COLOR, or a dumb terminal) — so redirected help stays plain text."""
    if os.environ.get("NO_COLOR") or os.environ.get("TERM") == "dumb":
        return False
    return sys.stdout.isatty()


def _banner() -> str:
    return f"ZOMBI2 {__version__} — a simulator of species trees, genomes, traits and sequences"


class ZombiHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Grouped help in the IQ-TREE house style: argument-group titles become UPPERCASE
    section headers (bold on a terminal), with a wide, aligned help column. The auto usage
    line is kept short by giving each command a hand-written ``usage=``."""

    def __init__(self, prog: str) -> None:
        width = min(shutil.get_terminal_size((90, 24)).columns - 2, 92)
        super().__init__(prog, max_help_position=32, width=width)

    def start_section(self, heading: str | None) -> None:
        if heading and heading not in ("positional arguments", "options", "optional arguments"):
            heading = heading.upper()
            if _use_color():
                heading = _BOLD + heading + _RESET
        super().start_section(heading)

    def _format_action(self, action: argparse.Action) -> str:
        # Hide the auto subcommand list from the top-level help — the commands are curated,
        # grouped by theme, in the description instead (avoids a duplicate, ungrouped dump).
        if isinstance(action, argparse._SubParsersAction):
            return ""
        return super()._format_action(action)


def _examples(*lines: str) -> str:
    """Build an ``EXAMPLES`` epilog block in the house style: a bold header on a TTY
    (plain when the output is piped), followed by the given lines verbatim. Safe because
    the parser's formatter is ``RawDescription``-based, so these line breaks are kept."""
    header = _BOLD + "EXAMPLES" + _RESET if _use_color() else "EXAMPLES"
    return header + "\n" + "\n".join(lines)


def _int_or_float(text: str) -> int | float:
    """Parse ``--max-family-size``: a plain integer is an absolute cap, a value with a
    decimal point is a fraction of the number of species (e.g. ``0.5`` -> half of N)."""
    return float(text) if "." in text else int(text)


def _add_species_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("--mode", dest="model", choices=("backward", "forward"), default="backward",
                   metavar="MODE",
                   help="backward: reconstructed tree conditioned on --tips extant species "
                        "(default); forward: complete tree grown in time, keeping extinct "
                        "lineages (and fossils)")
    g.add_argument("--tips", type=int, default=None, metavar="N",
                   help="number of extant species (backward default 50; forward: --tips OR --age)")
    g.add_argument("--age", type=float, default=None, metavar="T",
                   help="tree age / timescale, in the same time units as the rates "
                        "(backward default 1.0; forward: --tips OR --age)")
    g.add_argument("--age-type", choices=("crown", "stem"), default="crown", metavar="KIND",
                   help="[backward] interpret --age as crown (default) or stem age")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("diversification model",
                             "the rate process, chosen by --diversification (forward only)")
    g.add_argument("--diversification", choices=("constant", "clads", "diversity-dependent"),
                   default="constant", metavar="PROCESS",
                   help="constant-rate birth-death (default); clads = per-lineage rates that "
                        "shift at each speciation (ClaDS); diversity-dependent = rates decline "
                        "toward a carrying capacity")
    g.add_argument("--birth", type=float, nargs="+", default=[1.0], metavar="RATE",
                   help="speciation rate (default 1.0); several values with --shifts give an "
                        "episodic (skyline) model. For clads/diversity-dependent it is the "
                        "root/intrinsic rate λ₀ (a single value)")
    g.add_argument("--death", type=float, nargs="+", default=[0.3], metavar="RATE",
                   help="extinction rate (default 0.3); several values with --shifts give an "
                        "episodic (skyline) model. The constant μ for --diversification "
                        "diversity-dependent (clads uses --turnover instead)")
    g.add_argument("--shifts", type=float, nargs="+", default=None, metavar="AGE",
                   help="[episodic] rate-shift ages, present -> past (K-1 ages for K rate values)")

    g = p.add_argument_group("clads model", "with --diversification clads")
    g.add_argument("--clads-alpha", type=float, default=0.9, metavar="ALPHA",
                   help="speciation-rate trend per branch; α<1 = rates slow toward the present "
                        "(default 0.9)")
    g.add_argument("--clads-sigma", type=float, default=0.1, metavar="SIGMA",
                   help="lognormal spread of the per-branch rate shift (default 0.1)")
    g.add_argument("--turnover", type=float, default=0.0, metavar="EPS",
                   help="extinction/speciation ratio ε=μ/λ, in [0,1) (0 = pure birth; default 0.0)")

    g = p.add_argument_group("diversity-dependent model",
                             "with --diversification diversity-dependent")
    g.add_argument("--carrying-capacity", "-K", type=float, default=None, metavar="K",
                   help="carrying capacity K; the speciation rate is λ₀·(1−n/K) (required for "
                        "this model)")

    g = p.add_argument_group("clade-specific shifts", "forward only")
    g.add_argument("--clade-shift", action="append", nargs=3, type=float,
                   metavar=("AGE", "BIRTH", "DEATH"), default=None, dest="clade_shift",
                   help="a clade-specific rate shift: at AGE before the present, one random "
                        "lineage then alive (and its descendants) switches to speciation BIRTH / "
                        "extinction DEATH. Repeat for several shifting clades, e.g. "
                        "--clade-shift 3.0 2.5 0.1")

    g = p.add_argument_group("forward sampling & fossils", "only with --mode forward")
    g.add_argument("--sampling-fraction", type=float, default=1.0, metavar="RHO",
                   help="fraction of extant species sampled, 0<rho<=1 (default 1.0)")
    g.add_argument("--fossilization", type=float, default=0.0, metavar="PSI",
                   help="fossil (serial) sampling rate psi — fossilized birth-death "
                        "(default 0 = no fossils)")
    g.add_argument("--removal", type=float, default=1.0, metavar="R",
                   help="removal probability on sampling, 0<=r<=1 (r<1 keeps sampled ancestors; "
                        "default 1.0)")
    g.add_argument("--mass-extinction", action="append", nargs=2, type=float,
                   metavar=("AGE", "FRACTION"), default=None, dest="mass_extinction",
                   help="a mass extinction: at AGE before the present, each lineage dies with "
                        "probability FRACTION (0<FRACTION<=1). Repeat for several pulses, e.g. "
                        "--mass-extinction 1.0 0.75 --mass-extinction 2.5 0.5")

    g = p.add_argument_group("ghost lineages", "backward only")
    g.add_argument("--ghosts", action="store_true",
                   help="graft the extinct/unsampled 'ghost' lineages back onto the tree")
    g.add_argument("--ghost-method", choices=("rejection", "htransform"), default="rejection",
                   metavar="METHOD",
                   help="ghost-subtree sampler used with --ghosts (default rejection)")

    g = p.add_argument_group("run limits", "forward only")
    g.add_argument("--max-attempts", type=int, default=10000, metavar="N",
                   help="retries before giving up when the process goes extinct (default 10000)")
    g.add_argument("--max-lineages", type=int, default=1_000_000, metavar="N",
                   help="abort a run exceeding this many live lineages (default 1000000)")


def _add_rate_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    g.add_argument("--genome-model", dest="genome_model",
                   choices=("unordered", "ordered", "nucleotide"), default="unordered",
                   metavar="LEVEL",
                   help="genome level: unordered (default) evolves gene families with no "
                        "positional structure; ordered places genes on a chromosome where order "
                        "matters (adds inversion/transposition on gene segments; distance counted "
                        "in genes, not nucleotides); nucleotide evolves nucleotide-resolution "
                        "genomes by variable-length structural events, genes emerge as 'blocks' "
                        "(see the nucleotide sections)")
    g.add_argument("--rate-model", choices=("shared", "per-genome", "family"),
                   default="shared", metavar="MODEL",
                   help="rate heterogeneity within the unordered/ordered genome levels: shared: "
                        "same per-copy rates for every family (Rust for unordered; default); "
                        "per-genome: constant per-genome rates, linear growth (Python); "
                        "family: each family has its own rates, from --family-rates (Python). "
                        "Rearrangements (--inversion/--transposition on ordered genomes) need "
                        "shared")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group(
        "gene-family rates",
        "per copy for --rate-model shared/per-genome; "
        "per nucleotide for --genome-model nucleotide")
    g.add_argument("--dup", type=float, default=0.0, metavar="RATE", help="duplication rate")
    g.add_argument("--trans", type=float, default=0.0, metavar="RATE", help="transfer (HGT) rate")
    g.add_argument("--loss", type=float, default=0.0, metavar="RATE", help="loss/deletion rate")
    g.add_argument("--orig", type=float, default=0.0, metavar="RATE",
                   help="origination rate (per branch)")
    g.add_argument("--initial-families", type=int, default=None, metavar="N",
                   dest="initial_families",
                   help="number of gene families seeded at the root, for the unordered and ordered "
                        "genome levels (--genome-model unordered/ordered) (default: 20)")
    g.add_argument("--max-family-size", type=_int_or_float, default=None, metavar="CAP",
                   help="bound family growth: integer = absolute cap, decimal = fraction of the "
                        "number of species (e.g. 0.5) [not used by --genome-model nucleotide]")

    g = p.add_argument_group(
        "gene conversion",
        "intra-genome (ectopic) gene conversion; --rate-model shared on unordered genomes")
    g.add_argument("--conversion", type=float, default=0.0, metavar="RATE",
                   help="per-copy intra-genome gene-conversion rate: one copy of a family overwrites "
                        "another copy of the SAME family (concerted evolution), pulling their "
                        "coalescence toward the present. Fires only on families with >= 2 copies")
    g.add_argument("--conversion-bias", type=float, default=0.0, metavar="B",
                   dest="conversion_bias",
                   help="donor directionality in [0, 1]: 0 = donor drawn uniformly (default); "
                        "1 = always the family's oldest copy (homogenise toward the founder)")

    g = p.add_argument_group(
        "custom rate tables",
        "user-supplied per-family and per-branch rates (unordered genomes; Python engine "
        "except a receptivity-only --branch-rates, which stays on Rust)")
    g.add_argument("--family-rates", metavar="FILE", dest="family_rates",
                   help="TSV of explicit per-family duplication/transfer/loss rates "
                        "(columns: family duplication transfer loss). Selects --rate-model family; "
                        "families not listed fall back to --dup/--trans/--loss")
    g.add_argument("--branch-rates", metavar="FILE", dest="branch_rates",
                   help="TSV of per-branch transfer emission (donation-rate factor) and/or "
                        "receptivity (absorption weight) (columns: branch emission receptivity). "
                        "Emission scales that branch's transfer rate; receptivity biases which "
                        "branch receives")

    g = p.add_argument_group("output")
    g.add_argument("--write", dest="write", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "ancestral", "bed", "geneorder", "all"),
                   default=["profiles", "trees"],
                   help="which output files to write — any of {profiles, trace, trees, events, "
                        "transfers, summary, branch_events, reconciliations, layout, karyotype} or "
                        "'all' (default: profiles trees). 'reconciliations' writes "
                        "Reconciled_complete/extant.nwk (tips <species>|<gid>) — the truth input for "
                        "'tools recon-accuracy' and 'tools reconcile'. "
                        "species_tree.nwk is always written; 'profiles' alone takes the fast Rust "
                        "counts-only path; 'trace' (optionally with 'profiles') writes the compact "
                        "single-file event log Events_trace.tsv near counts-only speed, from which "
                        "gene trees can be reconstructed later on demand; 'branch_events' writes "
                        "Branch_events.tsv, the per-species-branch event counts (with an is_extant "
                        "flag). [ordered] 'layout' writes Gene_order.tsv (which chromosome each gene "
                        "sits on) and 'karyotype' writes Karyotype_trace.tsv (fission/fusion/"
                        "origination/loss) — both added automatically for a multi-chromosome or "
                        "fission/fusion run. [nucleotide] 'ancestral' simulates DNA and reconstructs the genome "
                        "(architecture + gzipped FASTA) at every node; 'bed' writes BED gene "
                        "annotations — genes.bed for the root genome and BED/<node>.bed per node "
                        "(needs --genes/--gff); 'geneorder' writes Geneorder_events.tsv, the "
                        "structural-event log with physical breakpoints (chrom/start/length/strand) "
                        "per branch — the input for gene-order / breakpoint export)")
    g.add_argument("--sparse", action="store_true",
                   help="write the profile as a sparse long table (Profiles_sparse.tsv: "
                        "family/species/copies, present cells only) instead of the dense matrix — "
                        "the scalable output for huge trees (needs 'profiles' in --write)")
    g.add_argument("--threads", type=int, default=1, metavar="N",
                   help="parallelise the counts-only profile simulation across N cores (only with "
                        "--write profiles; Poisson-thins the gene families into N independent "
                        "copies and sums them — a different but statistically identical "
                        "realisation, whose output depends on N). Default 1 (serial)")
    g.add_argument("--annotate-species", action="store_true",
                   help="label internal gene-tree nodes <gid>|<species-branch> (e.g. g570|i5)")

    g = p.add_argument_group(
        "reconciliation likelihoods (ALE)",
        "score every extant gene tree (forces the full gene-family path)")
    g.add_argument("--score-likelihoods", action="store_true",
                   help="also write Reconciliation_likelihoods.tsv: the ALE reconciliation "
                        "log-likelihood of every extant family's gene tree, under each "
                        "--score-model, at the --dup/--trans/--loss rates (zombi2.tools ALElite)")
    g.add_argument("--score-model", nargs="+", metavar="MODEL",
                   choices=("dated", "undated", "reldated"), default=["dated", "undated"],
                   help="ALE model(s) to score with (default: dated undated). dated = faithful "
                        "time-sliced likelihood (rates are per-unit-time); undated/reldated use "
                        "per-branch odds")
    g.add_argument("--score-nsteps", type=int, default=100, metavar="N",
                   help="dated model time-grid resolution (sub-steps per slice; default 100)")
    g.add_argument("--score-origination", choices=("root", "uniform"), default="root",
                   metavar="WHERE",
                   help="where the family enters the tree: 'root' (default; exact for root-seeded "
                        "families) or 'uniform' over branches")

    g = p.add_argument_group("structural events (rearrangements)",
                             "--genome-model ordered/nucleotide")
    g.add_argument("--inversion", type=float, default=None, metavar="RATE",
                   help="inversion rate — per gene copy for --genome-model ordered (default 0), "
                        "per nucleotide for nucleotide (default 0.001)")
    g.add_argument("--transposition", type=float, default=None, metavar="RATE",
                   help="transposition rate — per gene copy for --genome-model ordered, per "
                        "nucleotide for nucleotide (default 0)")
    g.add_argument("--mean-length", type=float, default=None, metavar="L", dest="mean_length",
                   help="mean length of an inversion/transposition segment (geometric): in genes "
                        "for --genome-model ordered (default 1 = single-gene events), in "
                        "nucleotides for nucleotide (default 100)")
    g.add_argument("--transposition-flip", type=float, default=0.0, metavar="P",
                   dest="transposition_flip",
                   help="probability a transposed segment reinserts reverse-complemented "
                        "(gene order reversed and strands flipped), for --genome-model ordered "
                        "(default 0 = always keep orientation)")
    g.add_argument("--translocation", type=float, default=0.0, metavar="RATE",
                   help="[nucleotide] translocation rate, per nucleotide: an arc moves to a "
                        "different chromosome of the same genome (needs >1 chromosome; distinct from "
                        "transposition, which stays on one chromosome) (default 0 = off)")
    g.add_argument("--n-chromosomes", type=int, default=1, metavar="N", dest="n_chromosomes",
                   help="number of chromosomes seeded at the root, for --genome-model "
                        "ordered/nucleotide (default 1). [ordered] the root's initial families are "
                        "spread across them; [nucleotide] each is an independent full-length copy of "
                        "the root chromosome. Rearrangements stay within a chromosome (see "
                        "--fission/--fusion for chromosome-level changes)")
    g.add_argument("--linear-chromosomes", action="store_true", dest="linear_chromosomes",
                   help="ordered chromosomes are linear (segments never wrap the origin), for "
                        "--genome-model ordered (default: circular, as for bacteria). Nucleotide "
                        "chromosomes are always circular")
    # chromosome-tier events (ordered + nucleotide genomes; off by default). When any is set — or
    # with more than one chromosome — the run also writes the karyotype (Gene_order.tsv /
    # Chromosomes.tsv layout) + Karyotype_trace.tsv.
    g.add_argument("--fission", type=float, default=0.0, metavar="RATE",
                   help="[ordered/nucleotide] chromosome fission rate, per chromosome: a chromosome "
                        "splits in two (linear: one breakpoint; circular: two) (default 0 = off)")
    g.add_argument("--fusion", type=float, default=0.0, metavar="RATE",
                   help="[ordered/nucleotide] chromosome fusion rate, per chromosome: two "
                        "chromosomes merge into one (default 0 = off)")
    g.add_argument("--chromosome-origination", type=float, default=0.0, metavar="RATE",
                   dest="chromosome_origination",
                   help="[ordered/nucleotide] chromosome origination rate, per genome: a de-novo "
                        "replicon (a plasmid) appears (default 0 = off)")
    g.add_argument("--chromosome-loss", type=float, default=0.0, metavar="RATE",
                   dest="chromosome_loss",
                   help="[ordered/nucleotide] chromosome loss rate, per chromosome: a whole "
                        "chromosome and its genes are lost (default 0 = off)")

    g = p.add_argument_group("nucleotide model", "with --genome-model nucleotide")
    g.add_argument("--initial-chromosomes", type=int, default=None, metavar="N",
                   dest="initial_chromosomes",
                   help="deprecated alias for --n-chromosomes (--genome-model nucleotide)")
    g.add_argument("--root-length", type=int, default=1000, metavar="BP",
                   help="length of the root chromosome, in nucleotides (default 1000)")
    g.add_argument("--insertion", type=float, default=0.0, metavar="RATE",
                   help="per-nucleotide intergenic insertion rate: lay down a run of novel "
                        "nucleotides (a fresh block) inside an intergene (default 0)")
    g.add_argument("--deletion", type=float, default=0.0, metavar="RATE",
                   help="per-nucleotide intergenic deletion rate: remove a run from within a single "
                        "intergene, never touching a gene, never below the min-genome floor "
                        "(default 0)")
    g.add_argument("--indel-mean-length", type=float, default=10.0, metavar="L",
                   dest="indel_mean_length",
                   help="mean length (in nucleotides) of an insertion/deletion run — geometric, a "
                        "separate knob from --mean-length (default 10)")

    g = p.add_argument_group("genes & intergenes",
                             "--genome-model nucleotide; declare genes to enable genic mode")
    g.add_argument("--gff", metavar="FILE", default=None,
                   help="a GFF3 annotation (optionally .gz) to start from a real genome: sets the "
                        "chromosome length and the gene coordinates (intergenes are the gaps). "
                        "Overlapping genes are trimmed to be disjoint. Supersedes --genes/--root-length")
    g.add_argument("--gff-seqid", metavar="ID", default=None,
                   help="which GFF sequence to read (default: the most-annotated one — the "
                        "chromosome of a single-chromosome bacterium)")
    g.add_argument("--genes", metavar="FILE", default=None,
                   help="BED/TSV of gene intervals on the root chromosome (columns: start end "
                        "[name], 0-based half-open) — an alternative to --gff. Event breakpoints "
                        "fall only in intergene positions so genes are never split; genes and "
                        "intergenes are recovered as separate tree sets")
    g.add_argument("--pseudogenization", type=float, default=0.0, metavar="P",
                   help="[genic] probability a loss hitting a gene demotes it to intergene "
                        "(sequence retained, a state change) instead of deleting it (default 0)")
    g.add_argument("--replacement", type=float, default=0.0, metavar="P",
                   help="[genic] probability a transfer is a homologous replacement (the copy "
                        "replaces the recipient's syntenic locus, located via flanking genes; "
                        "additive when no homolog) (default 0)")

    g = p.add_argument_group("sequences & ancestral genomes",
                             "--genome-model nucleotide, with --write ancestral")
    g.add_argument("--subst-model", choices=("jc69", "k80", "hky85", "gtr"), default="hky85",
                   metavar="MODEL",
                   help="nucleotide substitution model for the sequences (default hky85)")
    g.add_argument("--kappa", type=float, default=2.0, metavar="K",
                   help="transition/transversion ratio for k80/hky85 (default 2.0)")
    g.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="equilibrium base frequencies for hky85/gtr (default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="the 6 GTR exchangeabilities (default all 1)")
    g.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="discrete-Gamma across-site rate heterogeneity shape (default: none)")
    g.add_argument("--subst-rate", type=float, default=1.0, metavar="RATE",
                   help="overall substitutions/site per unit time — scales sequence divergence "
                        "(default 1.0)")
    g.add_argument("--genome-fasta", metavar="FILE", default=None,
                   help="the input genome's DNA (FASTA, optionally .gz) to seed the root sequence; "
                        "without it the root is drawn at random")


def _build_species_model(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Construct a species-tree model (BirthDeath / EpisodicBirthDeath / ClaDS /
    DiversityDependent) from the CLI args (validated)."""
    if args.model == "backward" and (args.fossilization or args.removal != 1.0
                                     or args.sampling_fraction != 1.0):
        parser.error("--fossilization / --removal / --sampling-fraction require --mode forward "
                     "(the backward reconstructed sampler assumes complete sampling)")
    if args.model == "backward" and args.mass_extinction:
        parser.error("--mass-extinction requires --mode forward (mass extinctions kill real "
                     "lineages forward in time; the backward reconstructed sampler never sees them)")
    # [(age, fraction), ...] pulses, or None; carried by whichever model is built
    mass_ext = args.mass_extinction

    if args.clade_shift and args.diversification != "constant":
        parser.error("--clade-shift is its own constant-background model; it does not combine "
                     "with --diversification clads/diversity-dependent")
    if args.diversification != "constant":
        return _build_heterogeneous_model(args, parser, mass_ext)
    if args.clade_shift:
        return _build_clade_shift_model(args, parser, mass_ext)

    episodic = args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1
    if not episodic:
        return BirthDeath(args.birth[0], args.death[0], fossilization=args.fossilization,
                          sampling_fraction=args.sampling_fraction, removal=args.removal,
                          mass_extinctions=mass_ext)
    shifts = args.shifts or []
    if len(args.birth) != len(args.death) or len(shifts) != len(args.birth) - 1:
        parser.error("episodic model needs len(--birth) == len(--death) == len(--shifts)+1 "
                     f"(got {len(args.birth)} birth, {len(args.death)} death, {len(shifts)} shifts)")
    return EpisodicBirthDeath(birth=args.birth, death=args.death, shifts=shifts,
                              fossilization=(args.fossilization or None),
                              sampling_fraction=args.sampling_fraction, removal=args.removal,
                              mass_extinctions=mass_ext)


def _build_heterogeneous_model(args: argparse.Namespace, parser: argparse.ArgumentParser,
                               mass_ext):
    """Build a ClaDS or DiversityDependent model — both forward-only, per-lineage/diversity-
    dependent rate processes selected by ``--diversification``."""
    if args.model != "forward":
        parser.error(f"--diversification {args.diversification} is a forward-in-time process; "
                     "add --mode forward")
    if args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1:
        parser.error(f"--diversification {args.diversification} takes a single --birth/--death "
                     "(no --shifts / multiple rates — those are the episodic model)")
    if args.fossilization or args.removal != 1.0:
        parser.error(f"--fossilization / --removal are not supported by --diversification "
                     f"{args.diversification}")
    if args.diversification == "clads":
        return ClaDS(args.birth[0], alpha=args.clads_alpha, sigma=args.clads_sigma,
                     turnover=args.turnover, sampling_fraction=args.sampling_fraction,
                     mass_extinctions=mass_ext)
    # diversity-dependent
    if args.carrying_capacity is None:
        parser.error("--diversification diversity-dependent needs --carrying-capacity/-K")
    return DiversityDependent(args.birth[0], args.death[0],
                              carrying_capacity=args.carrying_capacity,
                              sampling_fraction=args.sampling_fraction,
                              mass_extinctions=mass_ext)


def _build_clade_shift_model(args: argparse.Namespace, parser: argparse.ArgumentParser,
                             mass_ext):
    """Build a CladeShiftBirthDeath — constant background plus scheduled clade-specific rate
    shifts (forward-only, age mode)."""
    if args.model != "forward":
        parser.error("--clade-shift requires --mode forward (the shifts play out forward in time)")
    if args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1:
        parser.error("--clade-shift takes a single background --birth/--death (no --shifts; those "
                     "are the episodic model)")
    if args.fossilization or args.removal != 1.0:
        parser.error("--fossilization / --removal are not supported with --clade-shift")
    shifts = [(a, b, d) for a, b, d in args.clade_shift]
    return CladeShiftBirthDeath(args.birth[0], args.death[0], clade_shifts=shifts,
                                sampling_fraction=args.sampling_fraction,
                                mass_extinctions=mass_ext)


def _write_params_log(path: str, args: argparse.Namespace, summary: str) -> None:
    """Write the full set of run parameters to ``path`` — always, for reproducibility."""
    import datetime

    from . import __version__
    lines = ["# ZOMBI2 run parameters",
             f"zombi2_version\t{__version__}",
             f"timestamp\t{datetime.datetime.now().isoformat(timespec='seconds')}",
             f"command_line\t{' '.join(sys.argv)}"]
    for key, value in sorted(vars(args).items()):
        lines.append(f"{key}\t{value}")
    lines.append(f"result\t{summary}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _add_trait_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    g.add_argument("--model", choices=("bm", "ou", "eb", "mk", "threshold", "dec"), default="bm",
                   metavar="MODEL",
                   help="trait model: bm=Brownian motion, ou=Ornstein-Uhlenbeck, "
                        "eb=early burst/ACDC, mk=discrete k-state, threshold, "
                        "dec=geographic-range evolution (default: bm)")
    g.add_argument("--replicates", type=int, default=1, metavar="N",
                   help="simulate the trait this many times with the same parameters; writes "
                        "traits.tsv with one column per replicate (default: 1)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("continuous traits", "bm / ou / eb / threshold")
    g.add_argument("--sigma2", type=float, default=1.0, metavar="S2",
                   help="diffusion rate (default: 1.0)")
    g.add_argument("--x0", type=float, default=None, metavar="X0",
                   help="root value [bm/eb/threshold]; OU defaults it to --theta")
    g.add_argument("--trend", type=float, default=0.0, metavar="MU", help="directional drift [bm/eb]")

    g = p.add_argument_group("ornstein-uhlenbeck", "--model ou")
    g.add_argument("--alpha", type=float, default=1.0, metavar="A",
                   help="mean-reversion strength (default: 1.0)")
    g.add_argument("--theta", type=float, default=0.0, metavar="T", help="optimum (default: 0.0)")

    g = p.add_argument_group("early burst & Mk rate", "--model eb / --model mk")
    g.add_argument("--rate", type=float, default=1.0, metavar="R",
                   help="EB rate-of-change (negative = early burst) [eb], or the per-transition "
                        "rate [mk] (default: 1.0)")

    g = p.add_argument_group("discrete Mk", "--model mk")
    g.add_argument("--states", type=int, default=2, metavar="K",
                   help="number of states for the mk model (default: 2)")
    g.add_argument("--ordered", action="store_true",
                   help="only allow transitions between adjacent states (i <-> i±1)")
    g.add_argument("--q-matrix", default=None, metavar="FILE",
                   help="path to a whitespace/comma-separated k x k rate matrix (an arbitrary "
                        "Markov chain); overrides --states/--rate/--ordered")

    g = p.add_argument_group("threshold", "--model threshold")
    g.add_argument("--thresholds", default="0.0", metavar="CUTS",
                   help="comma-separated liability cut points (default: 0.0)")

    g = p.add_argument_group("DEC biogeography", "--model dec")
    g.add_argument("--areas", default="3", metavar="SPEC",
                   help="number of areas (e.g. 3) or comma-separated area labels (e.g. A,B,C) "
                        "(default: 3)")
    g.add_argument("--dispersal", type=float, default=0.1, metavar="RATE",
                   help="rate of gaining an area (dispersal) (default: 0.1)")
    g.add_argument("--extinction", type=float, default=0.1, metavar="RATE",
                   help="rate of losing an area (local extinction) (default: 0.1)")
    g.add_argument("--max-range-size", type=int, default=None, metavar="N",
                   help="maximum number of areas a range may span (default: all)")
    g.add_argument("--root-range", default=None, metavar="AREAS",
                   help="comma-separated area labels for the root range (e.g. A); default: random")


def _read_q_matrix(path: str):
    """Read a ``k x k`` rate matrix from a whitespace/comma-separated file.

    Blank lines and ``#`` comments are skipped; the diagonal is ignored (recomputed by
    :class:`~zombi2.Mk`). Each row is the *from*-state, each column the *to*-state.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.replace(",", " ").split()])
    if rows and any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("q-matrix rows must all have the same length (a square k x k matrix)")
    return rows


def _build_trait_model(args):
    x0 = args.x0
    if args.model == "bm":
        return BrownianMotion(sigma2=args.sigma2, x0=(0.0 if x0 is None else x0), trend=args.trend)
    if args.model == "ou":
        return OrnsteinUhlenbeck(sigma2=args.sigma2, alpha=args.alpha, theta=args.theta, x0=x0)
    if args.model == "eb":
        return EarlyBurst(sigma2=args.sigma2, rate=args.rate,
                          x0=(0.0 if x0 is None else x0), trend=args.trend)
    if args.model == "mk":
        if args.q_matrix:                                    # arbitrary user-supplied Markov chain
            return Mk(_read_q_matrix(args.q_matrix))
        if args.ordered:                                     # adjacent-only (meristic) character
            return Mk.ordered(args.states, args.rate)
        return Mk.equal_rates(args.states, args.rate)        # equal rates (all-to-all)
    # threshold
    thresholds = [float(t) for t in str(args.thresholds).split(",")]
    return ThresholdModel(thresholds=thresholds, sigma2=args.sigma2,
                          x0=(0.0 if x0 is None else x0))


def _parse_areas(text: str):
    """``--areas``: an integer count, or a comma-separated list of area labels."""
    text = str(text)
    if "," in text:
        return [a.strip() for a in text.split(",")]
    try:
        return int(text)
    except ValueError:
        return [text]


def _build_dec_model(args) -> DEC:
    return DEC(areas=_parse_areas(args.areas), dispersal=args.dispersal,
               extinction=args.extinction, max_range_size=args.max_range_size)


def _dec_root(args):
    """The root range from ``--root-range`` (a set of area labels), or ``None``."""
    if args.root_range is None:
        return None
    return {a.strip() for a in str(args.root_range).split(",")}


def _fmt_cell(value) -> str:
    """Format one trait value for a TSV cell (range tuple -> {A,B}, float -> 6 sig figs)."""
    if isinstance(value, tuple):
        return "{" + ",".join(str(v) for v in value) + "}"
    return f"{value:.6g}" if isinstance(value, float) else str(value)


def _replicate_table(results) -> str:
    """A wide ``node``-by-replicate table: one column (rep_1..rep_N) per simulation, one row per
    node (tips and ancestral)."""
    tree = results[0].tree
    header = "node\t" + "\t".join(f"rep_{i + 1}" for i in range(len(results)))
    lines = [header]
    for node in tree.nodes():
        cells = [_fmt_cell(res.label(res.node_values[node])) for res in results]
        lines.append(node.name + "\t" + "\t".join(cells))
    return "\n".join(lines) + "\n"


def _run_trait(args) -> str:
    """Evolve a trait along the supplied species tree and write the output folder.

    DEC (geographic ranges) runs the biogeography driver (it splits ranges at speciations);
    every other model overlays a trait with the standard driver. Both return a ``TraitResult``,
    so the output writing is shared.
    """
    if args.replicates < 1:
        raise ValueError("--replicates must be >= 1")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    os.makedirs(args.out, exist_ok=True)

    if args.model == "dec":
        dec = _build_dec_model(args)
        root = _dec_root(args)
        simulate = lambda rng: simulate_biogeography(tree, dec, root_state=root, rng=rng)  # noqa: E731
    else:
        model = _build_trait_model(args)
        simulate = lambda rng: simulate_traits(tree, model, rng=rng)  # noqa: E731

    if args.replicates > 1:
        rng = np.random.default_rng(args.seed)
        results = [simulate(rng) for _ in range(args.replicates)]
        with open(os.path.join(args.out, "traits.tsv"), "w") as f:
            f.write(_replicate_table(results))       # one column per replicate, all nodes
        return (f"wrote {args.replicates} trait replicates to {args.out}/traits.tsv "
                f"(model={args.model}; {len(tree.nodes())} nodes x {args.replicates} columns)")

    res = simulate(np.random.default_rng(args.seed))
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))             # every node: tips AND ancestral states
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")              # values annotated on every node too
    return (f"wrote traits to {args.out}/ (model={args.model}; "
            f"{len(tree.extant_leaves())} tip + {len(tree.internal_nodes())} ancestral values)")


# ═══════════════════════════════════════════════════════════════════════════════
# coevolve --couple traits:genes: trait-conditioned gene-family dynamics
# ═══════════════════════════════════════════════════════════════════════════════
def _add_trait_model_args(g) -> None:
    """Scalar trait-model flags (DEC ranges do not apply), added to the argument group ``g``.
    Used by the ``coevolve --couple traits:genes`` edge for the trait it simulates.

    ``--trait-model`` stores into ``args.model`` so :func:`_build_trait_model` is reused as-is.
    The Mk rate matrix reuses the command's shared ``--q-matrix`` (not re-added here).
    """
    g.add_argument("--trait-model", dest="model", default="bm", metavar="MODEL",
                   choices=("bm", "ou", "eb", "mk", "threshold"),
                   help="trait to evolve then couple to gene families: bm=Brownian motion, "
                        "ou=Ornstein-Uhlenbeck, eb=early burst, mk=discrete k-state, threshold "
                        "(default: bm). Use --trait-file to supply a precomputed trait instead")
    g.add_argument("--sigma2", type=float, default=1.0, metavar="S2",
                   help="diffusion rate [bm/ou/eb/threshold] (default: 1.0)")
    g.add_argument("--x0", type=float, default=None, metavar="X0",
                   help="root value [bm/eb/threshold]; OU defaults it to --theta")
    g.add_argument("--trend", type=float, default=0.0, metavar="MU", help="directional drift [bm/eb]")
    g.add_argument("--alpha", type=float, default=1.0, metavar="A",
                   help="OU mean-reversion strength [ou]")
    g.add_argument("--theta", type=float, default=0.0, metavar="T", help="OU optimum [ou]")
    g.add_argument("--rate", type=float, default=1.0, metavar="R",
                   help="EB rate-of-change (negative = early burst) [eb], or per-transition rate [mk]")
    g.add_argument("--states", type=int, default=2, metavar="K",
                   help="number of states for the mk model (default: 2)")
    g.add_argument("--ordered", action="store_true",
                   help="[mk] only allow transitions between adjacent states (i <-> i±1)")
    g.add_argument("--thresholds", default="0.0", metavar="CUTS",
                   help="comma-separated liability cut points [threshold] (default: 0.0)")


def _add_traits_genes_args(p: argparse.ArgumentParser) -> None:
    """The ``coevolve --couple traits:genes`` flags — a trait conditions a gene-family panel's
    loss/gain (formerly the standalone ``coevolve-genetrait`` command). Runs on a GIVEN -t tree.
    The trait's Mk matrix reuses the shared ``--q-matrix``; the panel writes with ``--write``."""
    g = p.add_argument_group("traits:genes trait model",
                             "--couple traits:genes: the trait to simulate (on a GIVEN tree; needs "
                             "-t). mk uses the shared --q-matrix")
    _add_trait_model_args(g)
    g.add_argument("--trait-file", default=None, metavar="TSV",
                   help="use a precomputed trait instead of simulating one: a node<TAB>value table "
                        "over ALL nodes (tips and ancestors), as 'zombi2 trait' writes with "
                        "nodes=all; values must be numeric (encode discrete states as numbers). "
                        "Overrides --trait-model")
    g.add_argument("--trait-center", action="store_true",
                   help="[discrete trait] center the state values around their mean so the trait "
                        "pushes retention both up and down — recommended for a binary "
                        "aerobic/anaerobic trait; by default states are 0,1,..,k-1")
    g.add_argument("--trait-steps", type=int, default=16, metavar="K",
                   help="[continuous trait] within-branch resolution: sub-segment each branch into "
                        "K pieces (default 16; ignored for discrete traits, which use their exact "
                        "stochastic map)")

    g = p.add_argument_group("traits:genes gene panel & coupling",
                             "the panel, its trait-neutral base rates, and which families respond")
    g.add_argument("--panel", type=int, default=50, metavar="N",
                   help="number of gene families in the panel (default 50)")
    g.add_argument("--loss", type=float, default=0.5, metavar="RATE",
                   help="baseline per-copy loss rate — the loss where the trait is neutral (default 0.5)")
    g.add_argument("--trans", type=float, default=1.0, metavar="RATE",
                   help="per-copy transfer (HGT) rate — the field-blind gain channel (default 1.0)")
    g.add_argument("--dup", type=float, default=0.0, metavar="RATE",
                   help="per-copy duplication rate, trait-independent (default 0)")
    g.add_argument("--orig", type=float, default=0.0, metavar="RATE",
                   help="background origination rate of brand-new, uncoupled families (default 0)")
    g.add_argument("--responsive", default="0.3", metavar="SPEC",
                   help="which families respond to the trait: an integer count, a fraction "
                        "(e.g. 0.3), a comma-separated id/index list (e.g. F3,F7,12), or @FILE of "
                        "ids/indices (default: 0.3 = 30%% of the panel, chosen at random)")
    g.add_argument("--weight", type=float, default=1.0, metavar="W",
                   help="coupling weight of each responsive family (default 1.0)")
    g.add_argument("--signed", action="store_true",
                   help="randomise the sign of each responsive weight (some families favoured by a "
                        "high trait value, some by a low one); by default all favour a high value")
    g.add_argument("--effect-loss", type=float, default=2.0, metavar="B",
                   help="retention coupling strength: a responsive family's loss scales by "
                        "exp(-effect_loss * weight * trait) (default 2.0; 0 = no coupling)")
    g.add_argument("--effect-gain", type=float, default=0.0, metavar="B",
                   help="optional HGT-activity coupling: a lineage's transfer rate scales by "
                        "exp(effect_gain * trait) (default 0 = trait-blind gain)")
    g.add_argument("--panel-root-fraction", dest="panel_root_fraction", type=float, default=0.5,
                   metavar="F",
                   help="[traits:genes + genes:traits JOINT model only] fraction of the panel present "
                        "at the root (default 0.5). In the joint model --theta-absent/--theta-present "
                        "are the trait's optima at an empty/full panel, and --trait-alpha/--trait-sigma2 "
                        "its OU dynamics")

    g = p.add_argument_group("traits:genes output")
    g.add_argument("--write", dest="write", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "all"), default=["profiles", "trees"],
                   help="which gene-family outputs to write — any of {profiles, trace, trees, "
                        "events, transfers, summary} or 'all' (default: profiles trees). "
                        "traits.tsv, trait_tree.nwk and coupling.tsv (the responsive-family "
                        "manifest) are always written alongside")
    g.add_argument("--sparse", action="store_true",
                   help="write the profile as a sparse long table (needs 'profiles' in --write)")
    g.add_argument("--annotate-species", action="store_true",
                   help="label internal gene-tree nodes <gid>|<species-branch> (e.g. g570|i5)")


def _parse_responsive(text: str):
    """``--responsive``: a count (int), a fraction (float), or an id/index list (``@FILE`` or CSV)."""
    text = str(text).strip()
    if text.startswith("@"):
        with open(text[1:]) as f:
            content = f.read()
        return [tok for tok in content.replace(",", " ").split() if tok]
    if "," in text:
        return [tok.strip() for tok in text.split(",") if tok.strip()]
    if "." in text:
        return float(text)
    return int(text)


def _load_trait_result(tree: Tree, path: str) -> TraitResult:
    """Load a precomputed trait: a ``node<TAB>value`` table over every node (numeric values).

    Returns a continuous-kind :class:`~zombi2.traits.TraitResult` (values used at node
    resolution; a supplied trait carries no within-branch stochastic map). Every node — tips and
    ancestors — must be present, since the gene simulation reads the trait on every branch.
    """
    name2node = {n.name: n for n in tree.nodes()}
    values: dict = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or parts[0] == "node":  # skip a header row
                continue
            node = name2node.get(parts[0])
            if node is None:
                continue
            try:
                values[node] = float(parts[1])
            except ValueError:
                raise ValueError(
                    f"--trait-file needs numeric trait values; got {parts[1]!r} for node "
                    f"{parts[0]!r}. Encode discrete states as numbers (e.g. 0 / 1).") from None
    missing = [n.name for n in tree.nodes() if n not in values]
    if missing:
        raise ValueError(
            f"--trait-file is missing values for {len(missing)} node(s) (e.g. {missing[:3]}); it "
            "must cover every node — tips AND ancestors — like the traits.tsv that 'zombi2 trait' "
            "writes (its nodes=all output)")
    return TraitResult(tree=tree, model=None, node_values=values, history=None, kind="continuous")


def _write_coupling_manifest(out: str, coupling: TraitGeneCoupling) -> None:
    """Write ``coupling.tsv`` — the per-family coupling weights plus the effect sizes, so the
    trait↔gene linkage that generated the profiles is recorded for downstream inference."""
    lines = [f"# effect_loss\t{coupling.effect_loss}",
             f"# effect_gain\t{coupling.effect_gain}",
             f"# base_loss\t{coupling.base_loss}",
             f"# transfer\t{coupling.transfer}",
             f"# duplication\t{coupling.duplication}",
             f"# origination\t{coupling.origination}",
             "family\tweight"]
    for i, fam in enumerate(coupling.panel_ids):
        lines.append(f"{fam}\t{coupling.weights[i]:.6g}")
    with open(os.path.join(out, "coupling.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _run_traits_genes(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """``coevolve --couple traits:genes``: simulate a trait, then evolve a gene-family panel whose
    loss/gain is conditioned on it. An overlay edge — runs along a GIVEN -t tree."""
    if not args.tree:
        parser.error("traits:genes runs on a GIVEN tree — pass -t/--tree (a trait conditions gene "
                     "content along it; there is nothing to grow)")
    if args.age is not None or args.tips is not None:
        parser.error("traits:genes uses the given -t tree; --age/--tips only apply to the "
                     "into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    parts = set(Genomes.WRITE_PARTS) if "all" in args.write else set(args.write)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --write")
    if args.panel < 1:
        raise ValueError("--panel must be >= 1")
    rng = np.random.default_rng(args.seed)

    # 1) the trait: simulate one (--trait-model) or load a precomputed one (--trait-file)
    if args.trait_file:
        result = _load_trait_result(tree, args.trait_file)
        trait_desc = f"file:{os.path.basename(args.trait_file)}"
    else:
        result = simulate_traits(tree, _build_trait_model(args), rng=rng)
        trait_desc = args.model

    # optional: center discrete states so the coupling is two-sided (recommended for binary)
    state_values = None
    if args.trait_center and result.kind == "discrete":
        k = len(result.model.states)
        state_values = [i - (k - 1) / 2.0 for i in range(k)]

    # 2) the coupling: choose the responsive families and the effect sizes
    coupling = TraitGeneCoupling.build(
        args.panel, _parse_responsive(args.responsive), weight=args.weight, signed=args.signed,
        effect_loss=args.effect_loss, effect_gain=args.effect_gain, base_loss=args.loss,
        transfer=args.trans, duplication=args.dup, origination=args.orig,
        state_values=state_values, rng=rng)

    if args.null == "cid":
        return _run_traits_genes_cid_null(args, parser)
    if args.null == "timing":
        parser.error("traits:genes has no 'timing' null; use --null neutral")
    if args.null == "neutral":
        coupling = coupling.null("neutral")

    # 3) run
    t0 = time.perf_counter()
    res = simulate_trait_linked_genomes(tree, result, coupling, trait_steps=args.trait_steps, rng=rng)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    res.genomes().write(args.out, include=parts, sparse=args.sparse,
                        annotate_species=args.annotate_species)
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.trait.to_tsv(nodes="all"))
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.trait.to_newick() + "\n")
    _write_coupling_manifest(args.out, coupling)
    if args.null == "neutral":
        _write_null_manifest(args.out, "traits:genes", "neutral", cut="trait -> panel loss/gain",
                             preserved="removed (effect_loss = effect_gain = 0; uncoupled panel)")
    tag = " [neutral null]" if args.null == "neutral" else ""
    return (f"wrote{tag} [{' '.join(sorted(parts))}] + trait to {args.out}/ (trait={trait_desc}, "
            f"panel {coupling.n_families} families, {coupling.n_responsive} responsive, "
            f"{len(tree.extant_leaves())} tips) in {dt:.3g} s")


# ═══════════════════════════════════════════════════════════════════════════════
# coevolve: the directed-coupling umbrella over {species, traits, genes}
# ═══════════════════════════════════════════════════════════════════════════════
_COEVOLVE_NODES = ("species", "traits", "genes")
# every directed edge in the coevolution design (docs/guide/coevolution.md); all six directed
# edges and the joint (both-arrow) models are implemented.
_COEVOLVE_EDGES = {
    "traits:species", "genes:species", "species:traits",
    "species:genes", "traits:genes", "genes:traits",
}


def _add_coevolve_mode_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("--couple", action="append", nargs="+", metavar="DRIVER:TARGET", default=None,
                   help="a directed coupling edge 'driver:target' over {species, traits, genes} — "
                        "the driver's state modulates the target's rates. Implemented: "
                        "'traits:species' (SSE), 'species:traits' (cladogenetic), their "
                        "combination = ClaSSE, 'genes:species' (key innovations), 'species:genes' "
                        "(cladogenetic genome), 'genes:traits' (a modifier gene switches a trait "
                        "optimum) and 'traits:genes' (a trait conditions a gene-family panel). "
                        "Repeatable; default traits:species. See docs/guide/coevolution.md for "
                        "the full edge set")
    g.add_argument("-t", "--tree", default=None, metavar="FILE",
                   help="input species tree (Newick) — required for the on-a-given-tree edges "
                        "(species:traits, species:genes, genes:traits, traits:genes). Omit for the "
                        "into-species edges (traits:species / ClaSSE / genes:species), which GROW "
                        "the tree via --age/--tips")
    g.add_argument("--age", type=float, default=None, metavar="T",
                   help="[into-species] crown age to grow for (the extant tip count is random)")
    g.add_argument("--tips", type=int, default=None, metavar="N",
                   help="[into-species] stop when this many extant tips first coexist (age random)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("--null", choices=("none", "neutral", "cid", "timing"), default="none",
                   help="generate the matched DECOUPLED null instead of the coupled model — cut the "
                        "driver→target arrow while keeping the target's variance (for calibrating a "
                        "detector's false-positive rate). 'neutral' (all edges): the driver stops "
                        "setting the rates. 'cid': the variance comes from a HIDDEN, uncorrelated "
                        "driver (traits:species natively; the gene/trait edges via a neutral "
                        "observed channel). 'timing': an at-speciation burst is spread along "
                        "branches (species:traits, species:genes). See "
                        "docs/guide/coevolution_nulls.md")
    g.add_argument("--hidden", type=int, default=2, metavar="H",
                   help="[--null cid, traits:species] number of hidden rate classes (2 = CID-2, "
                        "4 = CID-4)")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("SSE model", "--couple traits:species (trait drives diversification)")
    g.add_argument("--sse-model", dest="sse_model", choices=("bisse", "musse", "quasse", "hisse"),
                   default="bisse", metavar="MODEL",
                   help="which state-dependent model drives diversification: bisse (binary trait), "
                        "musse (k-state), quasse (continuous trait), hisse (binary trait + hidden "
                        "diversification classes) (default: bisse)")
    g.add_argument("--root-state", type=int, default=None, metavar="I",
                   help="[bisse/musse] root state index (default: the character's stationary "
                        "distribution)")

    g = p.add_argument_group("BiSSE", "--sse-model bisse (binary trait)")
    g.add_argument("--lambda0", type=float, default=1.0, metavar="RATE", help="speciation in state 0")
    g.add_argument("--lambda1", type=float, default=2.0, metavar="RATE", help="speciation in state 1")
    g.add_argument("--mu0", type=float, default=0.3, metavar="RATE", help="extinction in state 0")
    g.add_argument("--mu1", type=float, default=0.3, metavar="RATE", help="extinction in state 1")
    g.add_argument("--q01", type=float, default=0.1, metavar="RATE", help="transition rate 0 -> 1")
    g.add_argument("--q10", type=float, default=0.1, metavar="RATE", help="transition rate 1 -> 0")

    g = p.add_argument_group("MuSSE", "--sse-model musse (k-state trait)")
    g.add_argument("--birth", type=float, nargs="+", default=None, metavar="RATE",
                   help="per-state speciation rates (k values)")
    g.add_argument("--death", type=float, nargs="+", default=None, metavar="RATE",
                   help="per-state extinction rates (k values)")
    g.add_argument("--q-matrix", default=None, metavar="FILE",
                   help="path to a k x k anagenetic transition-rate matrix (same format as "
                        "'zombi2 trait --q-matrix')")

    g = p.add_argument_group("QuaSSE", "--sse-model quasse (continuous trait)")
    g.add_argument("--spec-low", type=float, default=0.5, metavar="RATE",
                   help="speciation rate at low trait values")
    g.add_argument("--spec-high", type=float, default=2.0, metavar="RATE",
                   help="speciation rate at high trait values")
    g.add_argument("--spec-center", type=float, default=0.0, metavar="X",
                   help="trait value at the middle of the speciation sigmoid")
    g.add_argument("--spec-slope", type=float, default=1.0, metavar="S",
                   help="steepness of the speciation sigmoid")
    g.add_argument("--qmu", type=float, default=0.1, metavar="RATE", help="constant extinction rate")
    g.add_argument("--diffusion", type=float, default=1.0, metavar="S2",
                   help="trait diffusion rate sigma^2 (Brownian motion)")
    g.add_argument("--root-value", type=float, default=0.0, metavar="X0", help="root trait value x0")

    g = p.add_argument_group("HiSSE", "--sse-model hisse (binary trait + hidden diversification "
                                      "classes)")
    g.add_argument("--hidden-classes", dest="hidden_classes", type=int, default=2, metavar="H",
                   help="number of hidden diversification classes (>= 2; default 2)")
    g.add_argument("--hidden-scale", dest="hidden_scale", type=float, default=3.0, metavar="S",
                   help="speciation spread across the hidden classes: the classes span the base "
                        "--lambda0/--lambda1 up to S times faster (geometric; default 3.0). The "
                        "observed transitions --q01/--q10 and extinction --mu0/--mu1 are shared")
    g.add_argument("--hidden-switch", dest="hidden_switch", type=float, default=0.1, metavar="RATE",
                   help="rate of switching between hidden classes (symmetric; default 0.1)")

    g = p.add_argument_group("cladogenetic kernel",
                             "--couple species:traits (speciation jumps the trait)")
    g.add_argument("--clado-shift", dest="clado_shift", type=float, default=0.3, metavar="P",
                   help="[discrete trait] probability a daughter hops to another state AT each "
                        "speciation (cladogenetic change; default 0.3)")
    g.add_argument("--clado-jump", dest="clado_jump", type=float, default=1.0, metavar="S2",
                   help="[continuous trait] variance of the Gaussian jump added to each daughter's "
                        "value AT each speciation (default 1.0)")

    g = p.add_argument_group(
        "gene-driven diversification",
        "--couple genes:species (key-innovation gene families; base rates reuse --lambda0/--mu0)")
    g.add_argument("--drivers", type=int, default=2, metavar="N",
                   help="number of binary 'driver' (key-innovation) gene families")
    g.add_argument("--driver-speciation", dest="driver_speciation", type=float, default=1.0,
                   metavar="B",
                   help="per-driver effect on log speciation: a present driver scales lambda by "
                        "exp(this) (>0 = a key innovation; default 1.0)")
    g.add_argument("--driver-extinction", dest="driver_extinction", type=float, default=0.0,
                   metavar="B",
                   help="per-driver effect on log extinction: a present driver scales mu by "
                        "exp(this) (default 0)")
    g.add_argument("--driver-loss", dest="driver_loss", type=float, default=0.1, metavar="RATE",
                   help="rate a present driver is lost/deleted (default 0.1)")
    g.add_argument("--driver-origination", dest="driver_origination", type=float, default=0.05,
                   metavar="RATE", help="rate an absent driver appears de novo (default 0.05)")
    g.add_argument("--driver-transfer", dest="driver_transfer", type=float, default=0.5,
                   metavar="RATE",
                   help="per-donor HGT rate of a driver — frequency-dependent gain: a driver in "
                        "more live genomes spreads faster (default 0.5)")
    g.add_argument("--root-drivers", dest="root_drivers", type=int, default=0, metavar="M",
                   help="number of drivers present at the root (the first m; default 0 = drivers "
                        "enter by origination)")
    g.add_argument("--driver-clado-loss", dest="driver_clado_loss", type=float, default=0.0,
                   metavar="P",
                   help="probability a daughter drops each driver it carries AT each speciation. "
                        ">0 adds the species:genes burst, making this the species<->genes JOINT "
                        "model (--couple genes:species --couple species:genes; default 0)")
    g.add_argument("--driver-clado-gain", dest="driver_clado_gain", type=float, default=0.0,
                   metavar="P",
                   help="probability a daughter gains each absent driver AT each speciation "
                        "(the species:genes burst; default 0)")

    g = p.add_argument_group("cladogenetic genome",
                             "--couple species:genes (on a GIVEN tree; needs -t)")
    g.add_argument("--genome-size", dest="genome_size", type=int, default=30, metavar="N",
                   help="number of families in the root genome (default 30)")
    g.add_argument("--gene-loss", dest="gene_loss", type=float, default=0.0, metavar="RATE",
                   help="anagenetic per-family loss rate along a branch (default 0)")
    g.add_argument("--gene-origination", dest="gene_origination", type=float, default=0.0,
                   metavar="RATE",
                   help="anagenetic origination rate of new families, per lineage (default 0). "
                        "With both anagenetic rates 0 the change is purely cladogenetic")
    g.add_argument("--clado-gene-loss", dest="clado_gene_loss", type=float, default=0.1, metavar="P",
                   help="probability a daughter drops each family AT each speciation (the "
                        "founder-effect burst; default 0.1)")
    g.add_argument("--clado-gene-gain", dest="clado_gene_gain", type=float, default=2.0,
                   metavar="MEAN",
                   help="mean number of new families a daughter gains AT each speciation "
                        "(Poisson; default 2.0)")

    g = p.add_argument_group("gene-conditioned trait",
                             "--couple genes:traits (on a GIVEN tree; needs -t)")
    g.add_argument("--modifier-gain", dest="modifier_gain", type=float, default=0.5, metavar="RATE",
                   help="rate the modifier gene is gained (absent -> present; default 0.5)")
    g.add_argument("--modifier-loss", dest="modifier_loss", type=float, default=0.5, metavar="RATE",
                   help="rate the modifier gene is lost (present -> absent; default 0.5)")
    g.add_argument("--root-modifier", dest="root_modifier", action="store_true",
                   help="start with the modifier gene present at the root")
    g.add_argument("--theta-absent", dest="theta_absent", type=float, default=0.0, metavar="T",
                   help="the trait's OU optimum while the modifier is absent (default 0)")
    g.add_argument("--theta-present", dest="theta_present", type=float, default=5.0, metavar="T",
                   help="the trait's OU optimum while the modifier is present (default 5) — "
                        "acquiring the gene pulls the trait toward this peak")
    g.add_argument("--trait-alpha", dest="trait_alpha", type=float, default=1.0, metavar="A",
                   help="OU mean-reversion strength of the trait (0 = Brownian; default 1.0)")
    g.add_argument("--trait-sigma2", dest="trait_sigma2", type=float, default=1.0, metavar="S2",
                   help="trait diffusion rate sigma^2 (default 1.0)")
    g.add_argument("--trait-x0", dest="trait_x0", type=float, default=None, metavar="X0",
                   help="root trait value (default: the optimum of the root modifier state)")

    _add_traits_genes_args(p)   # --couple traits:genes (a trait conditions a gene-family panel)


def _build_anagenetic_trait(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """The along-branch trait model for the ``species:traits`` edge on a **given** tree, taken from
    ``--sse-model``. No diversification happens on a fixed tree, so only the transition/diffusion
    structure is used (bisse/musse -> the Q as an :class:`Mk`; quasse -> Brownian ``--diffusion``);
    the speciation/extinction rates are inactive here. Returns ``(model, kind_label)``."""
    if args.sse_model == "quasse":
        return BrownianMotion(sigma2=args.diffusion, x0=args.root_value), "continuous"
    if args.sse_model == "musse":
        if args.q_matrix is None:
            parser.error("species:traits with --sse-model musse needs --q-matrix (the k-state "
                         "anagenetic transition matrix)")
        return Mk(_read_q_matrix(args.q_matrix)), "discrete"
    # bisse -> a binary Mk from the q01/q10 rates
    return Mk([[0.0, args.q01], [args.q10, 0.0]]), "discrete"


def _build_sse_model(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Construct the traits:species (SSE) model selected by ``--sse-model`` from the CLI args."""
    if args.sse_model == "bisse":
        return BiSSE(args.lambda0, args.lambda1, args.mu0, args.mu1, args.q01, args.q10)
    if args.sse_model == "musse":
        if args.birth is None or args.death is None or args.q_matrix is None:
            parser.error("--sse-model musse needs --birth and --death (k rates each) plus "
                         "--q-matrix (a k x k transition-rate matrix file)")
        return MuSSE(birth=args.birth, death=args.death, Q=_read_q_matrix(args.q_matrix))
    if args.sse_model == "hisse":
        H = args.hidden_classes
        if H < 2:
            parser.error("--sse-model hisse needs --hidden-classes >= 2")
        if args.hidden_scale <= 0:
            parser.error("--hidden-scale must be > 0")
        # H hidden classes spanning the base rates up to hidden_scale x faster (geometric); the
        # observed trait (q01/q10) and extinction (mu0/mu1) are shared across classes.
        factors = np.geomspace(1.0, args.hidden_scale, H)
        classes = [BiSSE(args.lambda0 * f, args.lambda1 * f, args.mu0, args.mu1, args.q01, args.q10)
                   for f in factors]
        return HiSSE(classes, hidden_transition=args.hidden_switch)
    # quasse: sigmoidal speciation in the trait + constant extinction (bounded for exact thinning)
    spec = QuaSSE.sigmoid(args.spec_low, args.spec_high, args.spec_center, args.spec_slope)
    bound = max(args.spec_low, args.spec_high) + args.qmu
    return QuaSSE(spec, lambda x: args.qmu, sigma2=args.diffusion,
                  rate_bound=bound, x0=args.root_value)


def _sse_tip_signal(res: TraitResult) -> str:
    """A short summary of the tip-state distribution — the diversification signal, for the log."""
    vals = list(res.labeled_values().values())
    if not vals:
        return ""
    if res.kind == "continuous":
        return f", tip trait mean {sum(vals) / len(vals):.3g}"
    from collections import Counter
    counts = Counter(vals)
    total = len(vals)
    frac = " ".join(f"{k}:{100 * n / total:.0f}%"
                    for k, n in sorted(counts.items(), key=lambda kv: str(kv[0])))
    return f", tip states {frac}"


def _write_coevolve_outputs(out: str, tree: Tree, res: TraitResult) -> None:
    """Write the shared coevolve outputs: the tree, the trait at every node, and the annotated
    trait tree."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # every node: tips AND ancestral states
    with open(os.path.join(out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")               # trait annotated on every node


# --------------------------------------------------------------------------- null models (--null)
def _write_null_manifest(out: str, edge: str, kind: str, *, cut: str, preserved: str,
                         extra: dict | None = None) -> None:
    """Record the null's provenance — which arrow was cut and how the target's variance was kept —
    so a downstream calibration is self-documenting. See docs/guide/coevolution_nulls.md."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "null_manifest.tsv"), "w") as f:
        f.write("field\tvalue\n")
        f.write(f"null\t{kind}\n")
        f.write(f"edge\t{edge}\n")
        f.write(f"cut_arrow\t{cut}\n")
        f.write(f"variance_preserved_by\t{preserved}\n")
        for k, v in (extra or {}).items():
            f.write(f"{k}\t{v}\n")


def _null_sse_model(model, args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Apply ``--null`` to a traits:species SSE model; returns ``(null_model, manifest_kwargs)``."""
    if args.null == "neutral":
        return model.null("neutral"), dict(cut="trait -> (lambda, mu)",
                                            preserved="removed (rates no longer depend on the trait)")
    if args.null == "cid":
        if isinstance(model, QuaSSE):
            parser.error("--null cid needs a discrete character (bisse/musse); QuaSSE is "
                         "continuous — use --null neutral")
        if isinstance(model, HiSSE):
            parser.error("--sse-model hisse is already a hidden-state model; its null is "
                         "--null neutral (a constant-rate tree)")
        # bisse -> the binary CID; musse -> the k-state CID (both hide the hidden class in output)
        return (model.null("cid", n_hidden=args.hidden),
                dict(cut="trait -> (lambda, mu)",
                     preserved=f"a hidden CID class (CID-{args.hidden})",
                     extra={"hidden_classes": args.hidden}))
    parser.error("traits:species has no 'timing' null; use --null neutral or --null cid")


def _null_species_traits_timing(args: argparse.Namespace, tree: Tree,
                                parser: argparse.ArgumentParser):
    """The species:traits ``timing`` null: drop cladogenesis, add a matched anagenetic rate so the
    same expected change is spread **along branches** (analytic, from the tree's branch stats).
    Returns ``(anagenetic_model, manifest_extra)``."""
    n_branches, total_len = _branch_count_and_length(tree)
    if total_len <= 0.0:
        parser.error("tree has zero total branch length; cannot spread the cladogenetic change")
    per_len = n_branches / total_len
    if args.sse_model == "quasse":                         # continuous: boost Brownian diffusion
        extra = args.clado_jump * per_len
        return (BrownianMotion(sigma2=args.diffusion + extra, x0=args.root_value),
                {"matched_sigma2_extra": f"{extra:.6g}"})
    if args.sse_model == "bisse":                          # binary: boost the Mk transition rates
        s = args.clado_shift * per_len
        return Mk([[0.0, args.q01 + s], [args.q10 + s, 0.0]]), {"matched_rate_extra": f"{s:.6g}"}
    parser.error("--null timing for species:traits supports --sse-model bisse or quasse (a matched "
                 "anagenetic rate); musse timing is not implemented — use --null neutral")


# The gene/trait CID nulls are the neutral-CHANNEL workflow: run the coupled model (so the target's
# heterogeneity is real), then hand the analyst a *neutral channel of the observed type on the same
# tree* — a neutral overlay genome, or a second independent neutral trait — while withholding the
# driver as ground-truth. The neutral overlay uses a plain transfer/loss panel (recorded in the
# manifest); re-run zombi2 genomes for other rates. See docs/guide/coevolution_nulls.md.
_NULL_OVERLAY = dict(duplication=0.0, transfer=1.0, loss=0.5, origination=0.0)   # neutral genome


def _run_genes_species_cid_null(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """genes:species CID null: the drivers shape a genuinely heterogeneous tree; a NEUTRAL overlay
    genome is the decoupled observed channel; the drivers are withheld as ground-truth."""
    if args.tree:
        parser.error("genes:species grows the tree (it is an OUTPUT); give --age/--tips, not -t")
    if (args.age is None) == (args.tips is None):
        parser.error("genes:species grows the tree — give exactly one of --age or --tips")
    model = GeneDiversification(
        args.drivers, lambda0=args.lambda0, mu0=args.mu0,
        driver_speciation=args.driver_speciation, driver_extinction=args.driver_extinction,
        loss=args.driver_loss, origination=args.driver_origination,
        transfer=args.driver_transfer, root_drivers=args.root_drivers)
    t0 = time.perf_counter()
    res = simulate_gene_diversification(model, age=args.age, n_tips=args.tips, seed=args.seed)
    tree = prune(res.tree)                                  # extant-only, for a clean null dataset
    seed2 = None if args.seed is None else args.seed + 1
    profiles = simulate_genomes(tree, output="profiles", initial_families=args.genome_size,
                                seed=seed2, **_NULL_OVERLAY)
    dt = time.perf_counter() - t0

    _write_profiles_only(args.out, tree, profiles)         # tree + OBSERVED neutral genome
    with open(os.path.join(args.out, "drivers_ground_truth.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))                   # the HIDDEN drivers (withheld from analysis)
    _write_drivers_manifest(args.out, model)
    _write_null_manifest(
        args.out, "genes:species", "cid", cut="gene content -> (lambda, mu)",
        preserved="a hidden driver panel (drivers_ground_truth.tsv) shaped the tree; Profiles.tsv is "
                  "a neutral overlay genome, decoupled from diversification",
        extra={"observed": f"neutral overlay: {args.genome_size} families, transfer=1.0 loss=0.5",
               "hidden_drivers": model.n_drivers})
    return (f"wrote genes:species [cid null] to {args.out}/ ({len(tree.leaves())} tips; "
            f"{len(profiles.families)} neutral observed families, {model.n_drivers} hidden drivers) "
            f"in {dt:.3g} s")


def _run_genes_traits_cid_null(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """genes:traits CID null: a modifier gene shapes a trait with real optimum shifts; a NEUTRAL
    overlay genome is the observed gene content; the modifier is withheld as ground-truth."""
    if not args.tree:
        parser.error("genes:traits runs on a GIVEN tree — pass -t/--tree")
    if args.age is not None or args.tips is not None:
        parser.error("genes:traits uses the given -t tree; --age/--tips only apply to into-species edges")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = GeneConditionedTrait(
        gene_gain=args.modifier_gain, gene_loss=args.modifier_loss, root_gene=args.root_modifier,
        theta_absent=args.theta_absent, theta_present=args.theta_present,
        alpha=args.trait_alpha, sigma2=args.trait_sigma2, x0=args.trait_x0)
    t0 = time.perf_counter()
    res = simulate_gene_conditioned_trait(tree, model, seed=args.seed)     # modifier shapes the trait
    seed2 = None if args.seed is None else args.seed + 1
    profiles = simulate_genomes(tree, output="profiles", initial_families=args.genome_size,
                                seed=seed2, **_NULL_OVERLAY)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())                         # OBSERVED neutral genome
    with open(os.path.join(args.out, "Presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:  # the trait (real optimum shifts), observed
        f.write("node\ttrait\n")
        for n in tree.nodes():
            f.write(f"{n.name}\t{res.node_trait[n]:.6g}\n")
    with open(os.path.join(args.out, "modifier_ground_truth.tsv"), "w") as f:   # the HIDDEN modifier
        f.write("node\tmodifier\n")
        for n in tree.nodes():
            f.write(f"{n.name}\t{int(res.gene.node_values[n])}\n")
    _write_null_manifest(
        args.out, "genes:traits", "cid", cut="gene presence -> trait optimum",
        preserved="a hidden modifier (modifier_ground_truth.tsv) shifted the trait's optimum; "
                  "Profiles.tsv is a neutral overlay genome, decoupled from the trait",
        extra={"observed": f"neutral overlay: {args.genome_size} families, transfer=1.0 loss=0.5"})
    return (f"wrote genes:traits [cid null] to {args.out}/ ({len(tree.extant_leaves())} tips; "
            f"{len(profiles.families)} neutral observed families, hidden modifier) in {dt:.3g} s")


def _run_traits_genes_cid_null(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """traits:genes CID null: a HIDDEN trait drives the panel's retention; a SECOND, independent
    neutral trait is the observed trait; the driving trait is withheld as ground-truth."""
    if not args.tree:
        parser.error("traits:genes runs on a GIVEN tree — pass -t/--tree")
    if args.age is not None or args.tips is not None:
        parser.error("traits:genes uses the given -t tree; --age/--tips only apply to into-species edges")
    if args.trait_file:
        parser.error("--null cid for traits:genes simulates its own hidden + observed traits; "
                     "--trait-file is not supported here")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    parts = set(Genomes.WRITE_PARTS) if "all" in args.write else set(args.write)
    if args.panel < 1:
        raise ValueError("--panel must be >= 1")
    rng = np.random.default_rng(args.seed)

    hidden_trait = simulate_traits(tree, _build_trait_model(args), rng=rng)   # drives the panel
    state_values = None
    if args.trait_center and hidden_trait.kind == "discrete":
        k = len(hidden_trait.model.states)
        state_values = [i - (k - 1) / 2.0 for i in range(k)]
    coupling = TraitGeneCoupling.build(
        args.panel, _parse_responsive(args.responsive), weight=args.weight, signed=args.signed,
        effect_loss=args.effect_loss, effect_gain=args.effect_gain, base_loss=args.loss,
        transfer=args.trans, duplication=args.dup, origination=args.orig,
        state_values=state_values, rng=rng)
    t0 = time.perf_counter()
    res = simulate_trait_linked_genomes(tree, hidden_trait, coupling, trait_steps=args.trait_steps,
                                        rng=rng)
    observed_trait = simulate_traits(tree, _build_trait_model(args), rng=rng)   # independent, decoupled
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    res.genomes().write(args.out, include=parts, sparse=args.sparse,
                        annotate_species=args.annotate_species)   # the panel (shaped by hidden trait)
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(observed_trait.to_tsv(nodes="all"))               # OBSERVED, decoupled trait
    with open(os.path.join(args.out, "trait_ground_truth.tsv"), "w") as f:
        f.write(hidden_trait.to_tsv(nodes="all"))                 # the HIDDEN driving trait
    _write_coupling_manifest(args.out, coupling)
    _write_null_manifest(
        args.out, "traits:genes", "cid", cut="trait -> panel loss/gain",
        preserved="a hidden trait (trait_ground_truth.tsv) shaped the panel; traits.tsv is a second, "
                  "independent neutral trait, decoupled from the panel",
        extra={"observed_trait": args.model, "responsive_families": coupling.n_responsive})
    return (f"wrote traits:genes [cid null] to {args.out}/ (observed trait decoupled from a "
            f"{coupling.n_families}-family panel shaped by a hidden trait, "
            f"{len(tree.extant_leaves())} tips) in {dt:.3g} s")


def _run_coevolve_mode(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Run the ``coevolve`` umbrella over the six directed edges (``--couple``): ``traits:species``
    (SSE), ``species:traits`` (cladogenetic) and their combination = **ClaSSE**; ``genes:species``
    (key innovations); ``species:genes`` (cladogenetic genome); ``genes:traits`` (gene-conditioned
    trait); and ``traits:genes`` (trait-conditioned genes). Each node-pair's two edges can also be
    combined: ``genes:species``+``species:genes`` = **co-diversification** and
    ``traits:genes``+``genes:traits`` = **trait-gene feedback** (as ``traits:species``+
    ``species:traits`` = ClaSSE). Whether the tree is grown (an arrow into species) or read from
    ``-t`` follows the arrows-into-S rule."""
    # --couple accepts both repeated flags and space-separated lists (append + nargs); flatten
    raw = args.couple or [["traits:species"]]
    edges = [e.strip().lower() for group in raw for e in group]
    for e in edges:
        if e not in _COEVOLVE_EDGES:
            parser.error(f"unknown --couple edge {e!r}: expected 'driver:target' over "
                         f"{{{', '.join(_COEVOLVE_NODES)}}} (e.g. traits:species); see "
                         "docs/guide/coevolution.md for the full edge set")
    eset = set(edges)

    # --null cuts a single directed arrow; a joint/both-arrows model has two, so decline it there
    # (run the null for one arrow at a time). Each edge below applies its own archetype.
    if args.null != "none" and len(eset) > 1:
        parser.error("--null cuts a single directed edge; a joint (both-arrows) model has two — "
                     "run the null for one arrow at a time (e.g. --couple traits:species --null cid)")

    # ---- joint (both-arrow) models: a node-pair with BOTH its directed edges on ----
    # species<->genes: driver gene content drives diversification AND bursts at each speciation
    # (one arrow into S -> the tree is an OUTPUT), the genomic analogue of ClaSSE.
    if eset == {"genes:species", "species:genes"}:
        return _run_co_diversification(args, parser)

    # traits<->genes: the trait and a coupled panel modulate each other (no arrow into S -> an
    # overlay on a given tree). The closed feedback loop writes a trait<->gene tip association.
    if eset == {"traits:genes", "genes:traits"}:
        return _run_trait_gene_feedback(args, parser)

    # traits:genes — a trait conditions a gene-family panel's loss/gain (formerly the standalone
    # 'coevolve-genetrait' command). An overlay edge (no arrow into S), so the tree is an INPUT.
    if "traits:genes" in eset:
        if eset != {"traits:genes"}:
            parser.error("traits:genes combines only with genes:traits (the trait-gene feedback "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_traits_genes(args, parser)

    # genes:species — gene content drives diversification (a forward joint loop). Combines with
    # species:genes above (co-diversification); other combinations are still on the roadmap.
    if "genes:species" in eset:
        if eset != {"genes:species"}:
            parser.error("genes:species combines only with species:genes (the co-diversification "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_genes_species(args, parser)

    # genes:traits — gene content conditions a trait (a modifier gene switches the trait's OU
    # optimum). An overlay edge (no arrow into S), so the tree is an INPUT; runs on a given -t tree.
    if "genes:traits" in eset:
        if eset != {"genes:traits"}:
            parser.error("genes:traits combines only with traits:genes (the trait-gene feedback "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_genes_traits(args, parser)

    # species:genes — speciation drives gene content (cladogenetic genome). An overlay edge (no
    # arrow into S), so the tree is an INPUT; runs on a given -t tree.
    if "species:genes" in eset:
        if eset != {"species:genes"}:
            parser.error("species:genes combines only with genes:species (the co-diversification "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_species_genes(args, parser)

    traits_species = "traits:species" in eset      # SSE arrow (trait -> diversification), into S
    species_traits = "species:traits" in eset      # cladogenetic arrow (speciation -> trait)
    clado = (Cladogenesis(jump_sigma2=args.clado_jump, shift=args.clado_shift)
             if species_traits else None)

    # species:traits ALONE — no arrow into S, so the tree is an INPUT (nothing to grow): evolve the
    # trait along the given tree with cladogenetic jumps at its speciation nodes.
    if species_traits and not traits_species:
        if not args.tree:
            parser.error("species:traits alone runs on a GIVEN tree — pass -t/--tree (no "
                         "diversification happens on this edge, so there is nothing to grow)")
        if args.age is not None or args.tips is not None:
            parser.error("species:traits alone uses the given -t tree; --age/--tips only apply to "
                         "the into-species edges that grow a tree")
        with open(args.tree) as f:
            tree = read_newick(f.read())
        model, kind = _build_anagenetic_trait(args, parser)
        nman = None
        if args.null == "neutral":                         # drop the at-speciation jump entirely
            clado = None
            nman = dict(cut="speciation -> trait jump", preserved="removed (anagenetic change only)")
        elif args.null == "timing":                        # spread the jump along branches (matched)
            clado = None
            model, extra = _null_species_traits_timing(args, tree, parser)
            nman = dict(cut="speciation -> trait jump",
                        preserved="spread along branches (matched anagenetic rate)", extra=extra)
        elif args.null == "cid":
            parser.error("species:traits has no 'cid' null (its driver is the speciation process, "
                         "not a state); use --null neutral or --null timing")
        t0 = time.perf_counter()
        res = simulate_traits(tree, model, cladogenesis=clado,
                              root_state=args.root_state, seed=args.seed)
        dt = time.perf_counter() - t0
        _write_coevolve_outputs(args.out, tree, res)
        if nman is not None:
            _write_null_manifest(args.out, "species:traits", args.null, **nman)
        tag = f" [{args.null} null]" if args.null != "none" else ""
        mode = "anagenetic" if args.null != "none" else "cladogenetic"
        return (f"wrote species:traits{tag} ({mode} {kind}) to {args.out}/ "
                f"({len(tree.extant_leaves())} tips{_sse_tip_signal(res)}) in {dt:.3g} s")

    # traits:species (SSE) or traits:species + species:traits (ClaSSE): an arrow INTO S, so the
    # tree is an OUTPUT — grow it forward with exactly one stopping condition (no input -t tree).
    if args.tree:
        parser.error("traits:species grows the tree (it is an OUTPUT); give --age/--tips, not an "
                     "input -t tree (that is the species:traits-alone edge)")
    if (args.age is None) == (args.tips is None):
        parser.error("traits:species grows the tree — give exactly one of --age or --tips")

    model = _build_sse_model(args, parser)
    nman = None
    if args.null != "none":                                # cut trait -> diversification
        model, nman = _null_sse_model(model, args, parser)
    t0 = time.perf_counter()
    res = simulate_sse(model, age=args.age, n_tips=args.tips, root_state=args.root_state,
                       cladogenesis=clado, seed=args.seed)
    dt = time.perf_counter() - t0
    _write_coevolve_outputs(args.out, res.tree, res)
    if nman is not None:
        _write_null_manifest(args.out, "traits:species", args.null, **nman)
    n_extant = len(res.tree.extant_leaves())
    tag = f" [{args.null} null]" if args.null != "none" else ""
    edge_label = ("traits:species+species:traits" if clado is not None else "traits:species") + tag
    model_label = f"ClaSSE {args.sse_model}" if clado is not None else f"SSE {args.sse_model}"
    return (f"wrote {edge_label} ({model_label}) to {args.out}/ "
            f"({n_extant} extant tips{_sse_tip_signal(res)}) in {dt:.3g} s")


def _write_drivers_manifest(out: str, model: GeneDiversification) -> None:
    """Write ``drivers_manifest.tsv`` — the per-driver effect sizes and rates behind the tree, so
    the gene↔diversification linkage that shaped the profiles is on record for inference."""
    root = ",".join(f"D{i}" for i in sorted(model.root_set)) or "-"
    lines = [f"# lambda0\t{model.lambda0:g}", f"# mu0\t{model.mu0:g}",
             f"# loss\t{model.loss:g}", f"# origination\t{model.origination:g}",
             f"# transfer\t{model.transfer:g}", f"# root_drivers\t{root}",
             "driver\tbeta_speciation\tbeta_extinction"]
    for i in range(model.n_drivers):
        lines.append(f"D{i}\t{model.beta_lambda[i]:.6g}\t{model.beta_mu[i]:.6g}")
    with open(os.path.join(out, "drivers_manifest.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _run_genes_species(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Grow a tree whose diversification is driven by a panel of binary key-innovation gene
    families (``genes:species``); the neutral genome is overlaid afterward with ``zombi2 genomes``
    on the resulting tree (exact under independent families)."""
    if args.tree:
        parser.error("genes:species grows the tree (it is an OUTPUT); give --age/--tips, not an "
                     "input -t tree")
    if (args.age is None) == (args.tips is None):
        parser.error("genes:species grows the tree — give exactly one of --age or --tips")

    model = GeneDiversification(
        args.drivers, lambda0=args.lambda0, mu0=args.mu0,
        driver_speciation=args.driver_speciation, driver_extinction=args.driver_extinction,
        loss=args.driver_loss, origination=args.driver_origination,
        transfer=args.driver_transfer, root_drivers=args.root_drivers)
    if args.null == "cid":
        return _run_genes_species_cid_null(args, parser)
    if args.null == "timing":
        parser.error("genes:species has no 'timing' null; use --null neutral or --null cid")
    if args.null == "neutral":
        model = model.null("neutral")
    t0 = time.perf_counter()
    res = simulate_gene_diversification(model, age=args.age, n_tips=args.tips, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(res.tree.to_newick() + "\n")          # the tree the drivers shaped
    with open(os.path.join(args.out, "drivers.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # per-node driver presence (0/1 columns)
    _write_drivers_manifest(args.out, model)
    n_extant = len(res.tree.extant_leaves())
    prev = " ".join(f"D{i}:{100 * p:.0f}%" for i, p in enumerate(res.tip_prevalence()))
    print(f"  overlay the neutral genome with: zombi2 genomes -t {args.out}/species_tree.nwk "
          f"--trans 1 --loss 0.5 --write profiles trees -o {args.out}")
    if args.null == "neutral":
        _write_null_manifest(args.out, "genes:species", "neutral", cut="gene content -> (lambda, mu)",
                             preserved="removed (drivers no longer set the rates)")
    tag = " [neutral null]" if args.null == "neutral" else ""
    return (f"wrote genes:species{tag} (key innovations) to {args.out}/ "
            f"({n_extant} extant tips, {model.n_drivers} drivers, tip prevalence {prev}) "
            f"in {dt:.3g} s")


def _run_co_diversification(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """The species<->genes JOINT model (``genes:species`` + ``species:genes``): the same driver
    families both drive diversification AND are reshuffled by a cladogenetic burst at every
    speciation. One arrow points into S, so the tree is an OUTPUT (grown jointly)."""
    if args.tree:
        parser.error("genes:species+species:genes grows the tree (it is an OUTPUT); give "
                     "--age/--tips, not an input -t tree")
    if (args.age is None) == (args.tips is None):
        parser.error("the species<->genes joint model grows the tree — give exactly one of "
                     "--age or --tips")
    if args.driver_clado_loss <= 0.0 and args.driver_clado_gain <= 0.0:
        parser.error("the species:genes arrow needs a cladogenetic burst on the drivers: set "
                     "--driver-clado-loss and/or --driver-clado-gain > 0 (with both 0 there is no "
                     "species:genes coupling — that is plain genes:species)")

    model = GeneDiversification(
        args.drivers, lambda0=args.lambda0, mu0=args.mu0,
        driver_speciation=args.driver_speciation, driver_extinction=args.driver_extinction,
        loss=args.driver_loss, origination=args.driver_origination,
        transfer=args.driver_transfer, root_drivers=args.root_drivers,
        cladogenetic_loss=args.driver_clado_loss, cladogenetic_gain=args.driver_clado_gain)
    t0 = time.perf_counter()
    res = simulate_co_diversification(model, age=args.age, n_tips=args.tips, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(res.tree.to_newick() + "\n")          # the tree the drivers shaped and burst
    with open(os.path.join(args.out, "drivers.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # per-node driver presence (0/1 columns)
    _write_drivers_manifest(args.out, model)
    n_extant = len(res.tree.extant_leaves())
    prev = " ".join(f"D{i}:{100 * p:.0f}%" for i, p in enumerate(res.tip_prevalence()))
    print(f"  overlay the neutral genome with: zombi2 genomes -t {args.out}/species_tree.nwk "
          f"--trans 1 --loss 0.5 --write profiles trees -o {args.out}")
    return (f"wrote genes:species+species:genes (co-diversification) to {args.out}/ "
            f"({n_extant} extant tips, {model.n_drivers} drivers, tip prevalence {prev}) "
            f"in {dt:.3g} s")


def _run_species_genes(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Evolve a genome down a GIVEN tree with a cladogenetic ('punctuational') burst of gene loss
    and gain at every speciation (the ``species:genes`` edge — speciation drives gene content)."""
    if not args.tree:
        parser.error("species:genes runs on a GIVEN tree — pass -t/--tree (speciation drives the "
                     "genome; there is no diversification to grow here)")
    if args.age is not None or args.tips is not None:
        parser.error("species:genes uses the given -t tree; --age/--tips only apply to the "
                     "into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = CladogeneticGenome(
        initial_families=args.genome_size, loss=args.gene_loss, origination=args.gene_origination,
        cladogenetic_loss=args.clado_gene_loss, cladogenetic_gain=args.clado_gene_gain)
    nman = None
    if args.null == "cid":
        parser.error("species:genes has no 'cid' null (its driver is the speciation process, not a "
                     "state); use --null neutral or --null timing")
    if args.null == "neutral":
        model = model.null("neutral")
        nman = dict(cut="speciation -> gene burst", preserved="removed (anagenetic turnover only)")
    elif args.null == "timing":
        model = model.null("timing", tree=tree)
        nman = dict(cut="speciation -> gene burst",
                    preserved="spread along branches (matched anagenetic rate)",
                    extra={"loss": f"{model.loss:.6g}", "origination": f"{model.origination:.6g}"})
    t0 = time.perf_counter()
    res = simulate_cladogenetic_genome(tree, model, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")              # the given tree, for provenance
    pm = res.profile_matrix()
    with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
        f.write(pm.to_tsv())
    with open(os.path.join(args.out, "Presence.tsv"), "w") as f:
        f.write(pm.to_tsv(presence=True))
    sizes = res.genome_sizes()
    with open(os.path.join(args.out, "genome_sizes.tsv"), "w") as f:
        f.write("node\tgenome_size\n")
        for node in tree.nodes():
            f.write(f"{node.name}\t{sizes[node]}\n")
    tips = tree.extant_leaves()
    mean_size = sum(sizes[t] for t in tips) / len(tips) if tips else 0
    if nman is not None:
        _write_null_manifest(args.out, "species:genes", args.null, **nman)
    tag = f" [{args.null} null]" if args.null != "none" else ""
    return (f"wrote species:genes{tag} (cladogenetic genome) to {args.out}/ "
            f"({len(tips)} tips, {len(pm.families)} families, mean genome {mean_size:.0f}) "
            f"in {dt:.3g} s")


def _run_genes_traits(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Evolve a trait down a GIVEN tree whose OU optimum is switched by a modifier gene's presence
    (the ``genes:traits`` edge — gene content conditions the trait)."""
    if not args.tree:
        parser.error("genes:traits runs on a GIVEN tree — pass -t/--tree (gene content conditions "
                     "the trait; there is no diversification to grow here)")
    if args.age is not None or args.tips is not None:
        parser.error("genes:traits uses the given -t tree; --age/--tips only apply to the "
                     "into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = GeneConditionedTrait(
        gene_gain=args.modifier_gain, gene_loss=args.modifier_loss, root_gene=args.root_modifier,
        theta_absent=args.theta_absent, theta_present=args.theta_present,
        alpha=args.trait_alpha, sigma2=args.trait_sigma2, x0=args.trait_x0)
    if args.null == "cid":
        return _run_genes_traits_cid_null(args, parser)
    if args.null == "timing":
        parser.error("genes:traits has no 'timing' null; use --null neutral")
    if args.null == "neutral":
        model = model.null("neutral")
    t0 = time.perf_counter()
    res = simulate_gene_conditioned_trait(tree, model, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")              # the given tree, for provenance
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # per-node modifier presence + trait value
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")               # trait annotated on every node

    tips = tree.extant_leaves()
    tv, gp = res.trait_values(), res.gene_presence()
    car = [tv[t] for t in tips if gp[t]]
    non = [tv[t] for t in tips if not gp[t]]
    car_m = f"{sum(car) / len(car):.2g}" if car else "-"
    non_m = f"{sum(non) / len(non):.2g}" if non else "-"
    if args.null == "neutral":
        _write_null_manifest(args.out, "genes:traits", "neutral", cut="gene presence -> trait optimum",
                             preserved="removed (theta_present = theta_absent; plain OU)")
    tag = " [neutral null]" if args.null == "neutral" else ""
    return (f"wrote genes:traits{tag} (gene-conditioned trait) to {args.out}/ "
            f"({len(tips)} tips; carrier trait mean {car_m} vs non-carrier {non_m}) in {dt:.3g} s")


def _run_trait_gene_feedback(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """The traits<->genes JOINT model (``traits:genes`` + ``genes:traits``): the trait and a coupled
    panel evolve together, each modulating the other. An overlay (no arrow into S), so it needs a
    given ``-t`` tree. Reuses --panel/--effect-loss/--loss/--trans for the panel and
    --theta-absent/--theta-present/--trait-alpha/--trait-sigma2 for the trait (see --panel-root-fraction)."""
    if not args.tree:
        parser.error("traits:genes+genes:traits runs on a GIVEN tree — pass -t/--tree (neither "
                     "arrow points into S, so there is nothing to grow)")
    if args.age is not None or args.tips is not None:
        parser.error("the traits<->genes joint model uses the given -t tree; --age/--tips only "
                     "apply to the into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = TraitGeneFeedback(
        n_families=args.panel, effect_loss=args.effect_loss, base_loss=args.loss, gain=args.trans,
        theta_low=args.theta_absent, theta_high=args.theta_present,
        alpha=args.trait_alpha, sigma2=args.trait_sigma2, x0=args.trait_x0,
        root_fraction=args.panel_root_fraction, steps=args.trait_steps)
    t0 = time.perf_counter()
    res = simulate_trait_gene_feedback(tree, model, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")                  # the given tree, for provenance
    with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
        f.write(res.profiles.to_tsv(presence=True))       # panel presence at the extant tips
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write("node\ttrait\tpanel_occupancy\n")         # the coupled trait + panel at every node
        for n in tree.nodes():
            f.write(f"{n.name}\t{res.node_trait[n]:.6g}\t{res.node_presence[n].mean():.6g}\n")

    corr = res.trait_gene_correlation()
    corr_s = f"{corr:.2f}" if corr == corr else "n/a"
    return (f"wrote traits:genes+genes:traits (trait-gene feedback) to {args.out}/ "
            f"({len(tree.extant_leaves())} tips; tip trait-panel corr {corr_s}) in {dt:.3g} s")


def _write_profiles_only(out: str, tree: Tree, profiles, sparse: bool = False) -> None:
    """Emit the reduced profiles-only output: tree + copy-number/presence matrices.

    With ``sparse=True`` the profile is written as a single COO long table
    (``Profiles_sparse.tsv``) that is O(present cells), so the output scales to trees
    where the dense families x species matrix would be astronomically large.
    """
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    if sparse:
        with open(os.path.join(out, "Profiles_sparse.tsv"), "w") as f:
            f.write(profiles.to_coo_tsv())
        return
    with open(os.path.join(out, "Profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())
    with open(os.path.join(out, "Presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))


def _extension_from_mean_length(mean_length: float | None) -> float | None:
    """User-facing knob → engine parameter. The user gives the *mean* segment length L (genes or
    nucleotides); the engine wants the geometric continuation probability. ``None`` keeps the
    per-level default; otherwise ``extension = 1 - 1/L`` (L=1 → single-element events)."""
    if mean_length is None:
        return None
    if mean_length < 1.0:
        raise ValueError(f"--mean-length must be >= 1 (a segment spans at least one unit), got {mean_length}")
    return 1.0 - 1.0 / mean_length


def _run_genomes(tree: Tree, args: argparse.Namespace,
                 parser: argparse.ArgumentParser) -> str:
    """Simulate gene families along ``tree``, write output, and return a one-line summary.

    The default ``shared`` rate model runs on the Rust engine automatically (``simulate_genomes``
    raises a build hint if the extension is missing); ``per-genome`` runs on Python.
    """
    args.extension = _extension_from_mean_length(args.mean_length)   # mean-length knob → engine p
    parts = set(Genomes.WRITE_PARTS) if "all" in args.write else set(args.write)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --write")
    if args.threads > 1:
        # --threads parallelises ONLY the counts-only Rust fast path: built-in shared rates on an
        # unordered genome, profiles-only, no conversion / branch-rates / scoring. Reject every other
        # combination up front with a flag-level message — otherwise per-genome/family rates crash
        # deep in the engine and --conversion / ordered silently run serial (the threads ignored).
        reason = None
        if parts != {"profiles"}:
            reason = "use it with exactly --write profiles"
        elif args.genome_model != "unordered":
            reason = f"--genome-model {args.genome_model} runs serially; use --genome-model unordered"
        elif args.rate_model != "shared" or args.family_rates:
            reason = "the built-in --rate-model shared is required (per-genome/family run on Python)"
        elif args.conversion:
            reason = "--conversion runs on the full (serial) path"
        elif args.branch_rates:
            reason = "--branch-rates runs on the full (serial) path"
        elif getattr(args, "score_likelihoods", False):
            reason = "--score-likelihoods forces the full (serial) gene-tree path"
        if reason is not None:
            parser.error(f"--threads > 1 parallelises only the counts-only path: {reason}")

    if args.genome_model == "nucleotide":
        if args.initial_families is not None:
            parser.error("--initial-families is for the unordered genome level "
                         "(--genome-model unordered); the nucleotide model uses --n-chromosomes")
        if getattr(args, "score_likelihoods", False):
            parser.error("--score-likelihoods scores reconstructed gene-family trees, which the "
                         "nucleotide genome model does not produce; use --genome-model "
                         "unordered/ordered to score reconciliation likelihoods")
        if args.linear_chromosomes:
            parser.error("--linear-chromosomes applies to --genome-model ordered; nucleotide "
                         "chromosomes are always circular")
        return _run_nucleotides(tree, args, parts)

    if args.initial_chromosomes is not None:
        parser.error("--initial-chromosomes is only for --genome-model nucleotide; the "
                     "unordered and ordered genome levels use --initial-families")
    if args.translocation:
        parser.error("--translocation is only for --genome-model nucleotide (it moves an arc "
                     "between chromosomes of a multi-chromosome nucleotide genome)")

    ordered = args.genome_model == "ordered"
    initial_families = 20 if args.initial_families is None else args.initial_families
    args.initial_families = initial_families  # record the effective value in the params log
    # user-supplied custom rate tables are unordered-only (Python engine, except a
    # receptivity-only --branch-rates on a plain SharedRates, which stays on Rust)
    if (args.family_rates or args.branch_rates) and args.genome_model != "unordered":
        parser.error("--family-rates / --branch-rates are only for --genome-model unordered")
    if not (0.0 <= args.transposition_flip <= 1.0):
        parser.error("--transposition-flip must be a probability in [0, 1]")
    if args.transposition_flip and not ordered:
        parser.error("--transposition-flip applies to transpositions on an ordered chromosome; "
                     "use --genome-model ordered")
    if args.n_chromosomes < 1:
        parser.error("--n-chromosomes must be >= 1")
    if args.n_chromosomes != 1 and not ordered:  # nucleotide is handled in its own branch above
        parser.error("--n-chromosomes applies to --genome-model ordered or nucleotide")
    if args.linear_chromosomes and not ordered:
        parser.error("--linear-chromosomes applies to --genome-model ordered")
    chrom_tier = bool(args.fission or args.fusion or args.chromosome_origination
                      or args.chromosome_loss)
    if chrom_tier and not ordered:
        parser.error("--fission / --fusion / --chromosome-origination / --chromosome-loss apply to "
                     "--genome-model ordered or nucleotide")
    if ordered and (args.n_chromosomes > 1 or chrom_tier):
        # auto-surface the karyotype outputs when non-trivial, so a multi-chromosome or fission /
        # fusion run captures its layout (and genealogy) without the user asking; single-chromosome
        # runs are untouched.
        parts.add("layout")
        if chrom_tier:
            parts.add("karyotype")
    family_mode = args.rate_model == "family" or args.family_rates is not None
    if args.rate_model == "family" and args.family_rates is None:
        parser.error("--rate-model family needs a --family-rates FILE")
    if family_mode and args.rate_model == "per-genome":
        parser.error("--family-rates is itself a per-family model; do not also pass "
                     "--rate-model per-genome")
    if args.conversion and family_mode:
        parser.error("gene conversion (--conversion) needs --rate-model shared; the per-family "
                     "table carries no conversion rate")

    rates = None  # None => use the D/T/L/O shorthand (plain shared, unordered — the Rust fast path)
    if args.rate_model == "per-genome":
        if args.conversion:
            parser.error("gene conversion (--conversion) needs --rate-model shared; "
                         "per-genome rates do not carry it")
        if ordered and (args.inversion is not None or args.transposition is not None):
            parser.error("rearrangements (--inversion/--transposition) need --rate-model shared; "
                         "per-genome rates do not carry them")
        rates = PerGenomeRates(args.dup, args.trans, args.loss, args.orig)
    elif ordered:  # shared per-copy rates + rearrangements on an ordered chromosome
        if args.conversion:
            parser.error("gene conversion (--conversion) is only supported on unordered genomes "
                         "(--genome-model unordered)")
        inv = 0.0 if args.inversion is None else args.inversion
        tps = 0.0 if args.transposition is None else args.transposition
        args.inversion, args.transposition = inv, tps  # record effective values in the params log
        rates = SharedRates(args.dup, args.trans, args.loss, args.orig,
                            inversion=inv, transposition=tps,
                            chromosome_origination=args.chromosome_origination,
                            chromosome_loss=args.chromosome_loss,
                            fission=args.fission, fusion=args.fusion)
    elif family_mode:  # each family its own rates, from the table (unlisted -> --dup/--trans/--loss)
        rates = FamilySampledRates(duplication=args.dup, transfer=args.trans, loss=args.loss,
                                   origination=args.orig,
                                   rates=read_family_rates(args.family_rates))

    # --branch-rates overlay: per-branch transfer emission (BranchRates) + receptivity (TransferModel)
    transfers = None
    if args.branch_rates is not None:
        emission, receptivity = read_branch_rates(args.branch_rates)
        if rates is None:  # plain shared base — carry the conversion rate through the overlay
            rates = SharedRates(args.dup, args.trans, args.loss, args.orig,
                                conversion=args.conversion)
        if emission:
            rates = BranchRates(rates, factors=emission, events=("transfer",))
        if receptivity:
            transfers = TransferModel(receptivity=receptivity)

    if rates is None:  # plain shared (unordered): D/T/L/O shorthand + optional gene conversion
        model_kw = dict(duplication=args.dup, transfer=args.trans, loss=args.loss,
                        origination=args.orig, conversion=args.conversion)
    else:
        model_kw = dict(rates=rates)
    if args.conversion:
        model_kw["conversions"] = ConversionModel(bias=args.conversion_bias)
    if transfers is not None:
        model_kw["transfers"] = transfers
    rate_kw = dict(**model_kw, initial_families=initial_families,
                   max_family_size=args.max_family_size, seed=args.seed)
    if ordered:
        ext = args.extension  # ordered event length is counted in genes; None -> single-gene events
        flip = args.transposition_flip  # P(a transposed segment reinserts reverse-complemented)
        n_chrom = args.n_chromosomes     # number of chromosomes to seed at the root
        circular = not args.linear_chromosomes  # linear = eukaryotic ends (no wrap); circular default
        rate_kw["genome_factory"] = (
            lambda ids, _e=ext, _f=flip, _n=n_chrom, _c=circular: OrderedGenome(
                ids, extension=_e, transposition_flip=_f, n_chromosomes=_n, circular=_c)
        )

    # scoring reconciliation likelihoods needs the full gene-family genealogy, so it forces the
    # full path (the fast counts-only / trace paths don't reconstruct gene trees).
    score = getattr(args, "score_likelihoods", False)

    t0 = time.perf_counter()
    if parts == {"profiles"} and not score and not ordered and not args.conversion:
        # counts-only Rust fast path: no genealogy reconstructed (parallel when --threads > 1).
        # Conversion only reshapes gene trees (copy numbers are unchanged), so it has no bearing on
        # a profiles-only run and stays on the full path below.
        profiles = simulate_genomes(tree, output="profiles", threads=args.threads, **rate_kw)
        dt = time.perf_counter() - t0
        _write_profiles_only(args.out, tree, profiles, sparse=args.sparse)
        n_families = len(profiles.families)
    elif ("trace" in parts and parts <= {"trace", "profiles"} and not score and not ordered
          and not args.conversion):
        # event-trace fast path: compact Events_trace.tsv (+ profile), no per-event objects,
        # no gene-tree reconstruction — near counts-only speed, trees reconstructable later
        trace = simulate_genomes(tree, output="trace", **rate_kw)
        dt = time.perf_counter() - t0
        trace.write(args.out, include=parts, sparse=args.sparse)
        n_families = len(trace.profiles.families)
    else:
        genomes = simulate_genomes(tree, **rate_kw)
        dt = time.perf_counter() - t0
        genomes.write(args.out, include=parts, sparse=args.sparse,
                      annotate_species=args.annotate_species)
        n_families = len(genomes.profiles.families)
        if score:
            _write_reconciliation_likelihoods(genomes, args)
    suffix = " + Reconciliation_likelihoods.tsv" if score else ""
    return (f"wrote [{' '.join(sorted(parts))}]{suffix} to {args.out}/ "
            f"({len(tree.leaves())} tips, {n_families} gene families) in {dt:.3g} s")


def _write_reconciliation_likelihoods(genomes, args: argparse.Namespace) -> None:
    """Score every extant family's gene tree (ALElite) and write Reconciliation_likelihoods.tsv."""
    from .tools.reconciliation import write_scores_tsv

    models = list(dict.fromkeys(args.score_model))  # de-dupe, keep order
    rows = genomes.reconciliation_likelihoods(
        args.dup, args.trans, args.loss, models=models,
        origination=args.score_origination, n_steps=args.score_nsteps,
    )
    write_scores_tsv(rows, os.path.join(args.out, "Reconciliation_likelihoods.tsv"), models=models)


def _run_tools_reconcile(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools reconcile`` — the ALE reconciliation log-likelihood (ALElite) of one or
    more given gene trees, *evaluated* at fixed DTL rates (no rate fitting)."""
    from .tools import (GeneTree, SpeciesTree, FamilyScore, reconciliation_likelihood,
                        write_scores_tsv)

    with open(args.species_tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.species_tree} is not a usable species tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    species_names = {n.name for n in tree.leaves()}
    sp = SpeciesTree.from_tree(tree)               # build the dated species index once

    with open(args.gene_tree) as f:                # one Newick per non-blank, non-comment line
        newicks = [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    if not newicks:
        raise ValueError(f"no gene trees found in {args.gene_tree}")

    models = list(dict.fromkeys(args.model))       # de-dupe, keep order
    rows = []
    for i, nwk in enumerate(newicks, 1):
        gt = GeneTree.from_newick(nwk)
        unknown = gt.species_set() - species_names
        if unknown:
            raise ValueError(
                f"gene tree {i} references species absent from the species tree: "
                f"{', '.join(sorted(unknown))} — tip labels must be '<species>|<gid>' with "
                "<species> a species-tree leaf.")
        tips = sum(g.is_leaf for g in gt.nodes)
        logliks = {m: reconciliation_likelihood(
                        gene_tree=gt, species_tree=sp,
                        duplication=args.dup, transfer=args.trans, loss=args.loss,
                        model=m, origination=args.origination, n_steps=args.n_steps)
                   for m in models}
        rows.append(FamilyScore(family=str(i), extant_tips=tips, logliks=logliks))

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "Reconciliation_likelihoods.tsv")
        write_scores_tsv(rows, path, models=tuple(models))
        print(f"wrote {path} ({len(rows)} gene tree(s) x {len(models)} model(s))")
    elif len(rows) == 1 and len(models) == 1:
        print(f"{rows[0].logliks[models[0]]:.6f}")     # bare number — scripting-friendly
    else:                                              # same columns as write_scores_tsv
        print("family\textant_copies\t" + "\t".join(f"{m}_loglik" for m in models))
        for r in rows:
            print(f"{r.family}\t{r.extant_tips}\t"
                  + "\t".join(f"{r.logliks[m]:.6f}" for m in models))
    return 0


def _run_tools_simulate(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools simulate`` — forward-simulate gene families under the ALE undated/reldated
    model (the generative twin of 'reconcile') and write per-family ground-truth reconciliations."""
    from .tools.reconciliation.undated import UndatedDTL
    from .tools.reconciliation.undated_sim import simulate_undated

    with open(args.species_tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.species_tree} is not a usable species tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    transfers = "dated" if args.model == "reldated" else "global"
    model = UndatedDTL(args.dup, args.trans, args.loss)
    try:
        res = simulate_undated(tree, model, n_families=args.families,
                               origination=args.origination, transfers=transfers,
                               seed=args.seed, max_events=args.max_events)
    except (ValueError, RuntimeError) as e:
        parser.error(str(e))

    c = res.event_counts
    print(f"simulated {res.n_families} families under {args.model} "
          f"(d={args.dup}, t={args.trans}, l={args.loss}): "
          f"{res.n_surviving} survived, {res.n_extinct} went extinct")
    print(f"events: {c.get('D', 0)} D, {c.get('T', 0)} T, {c.get('L', 0)} L, {c.get('S', 0)} S")

    if args.score:
        from .tools.reconciliation.undated import undated_joint_loglik
        ll = undated_joint_loglik(res.gene_trees(), res.species_tree, model,
                                  origination=args.origination, transfers=transfers,
                                  n_extinct=res.n_extinct)
        print(f"joint undated log-likelihood of the {res.n_surviving} survivors "
              f"(+{res.n_extinct} extinct) under the generating odds: {ll:.6f}")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        n_extant = _write_undated_sim(res, args.out)
        print(f"wrote Reconciled_complete.nwk, Reconciled_extant.nwk ({n_extant} survivors), "
              f"Reconciliation_events.tsv and Gene_family_profiles.tsv into {args.out}/")
    return 0


def _write_undated_sim(res, out: str) -> int:
    """Write a simulated result's ground-truth reconciliations: bare annotated Newicks (one family
    per line — the format 'recon-accuracy' reads) plus a flat S/D/T/L event table. Returns the
    number of surviving (extant) families written."""
    complete_lines, extant_lines = [], []
    ev_lines = ["family\tevent\tspecies\trecipient\ttime\tgene"]
    for i, recon in enumerate(res.reconciliations, 1):
        if recon.complete is not None:
            complete_lines.append(recon.complete)
        if recon.extant is not None:
            extant_lines.append(recon.extant)
        for e in recon.events:
            ev_lines.append(f"{i}\t{e.event}\t{e.species}\t{e.recipient or ''}\t"
                            f"{e.time:.10g}\t{e.gene or ''}")
    prof_lines = ["family\t" + "\t".join(res.leaf_names)]
    for fam, counts in res.profile_rows():
        prof_lines.append(fam + "\t" + "\t".join(str(c) for c in counts))
    for name, lines in (("Reconciled_complete.nwk", complete_lines),
                        ("Reconciled_extant.nwk", extant_lines),
                        ("Reconciliation_events.tsv", ev_lines),
                        ("Gene_family_profiles.tsv", prof_lines)):
        with open(os.path.join(out, name), "w") as f:
            f.write("\n".join(lines) + ("\n" if lines else ""))
    return len(extant_lines)


_TREEDIST_COLS = ("tree", "n_leaves", "rf", "rf_norm", "rf_unrooted",
                  "branch_score", "quartet", "quartet_norm", "matching", "matching_norm")


def _treedist_row(label: str, c) -> str:
    """One TSV line for a :class:`~zombi2.tools.treedist.TreeComparison` (blank when a metric was
    skipped: quartet over max-leaves, or matching over max-leaves / SciPy missing)."""
    quartet = "" if c.quartet is None else str(c.quartet)
    quartet_norm = "" if c.quartet_normalized is None else f"{c.quartet_normalized:.6f}"
    matching = "" if c.matching is None else str(c.matching)
    matching_norm = "" if c.matching_normalized is None else f"{c.matching_normalized:.6f}"
    return "\t".join((
        label, str(c.n_leaves), str(c.rf), f"{c.rf_normalized:.6f}", str(c.rf_unrooted),
        f"{c.branch_score:.6f}", quartet, quartet_norm, matching, matching_norm,
    ))


def _run_tools_treedist(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools treedist`` — RF, branch-score, and quartet distances between a reference
    tree and one or more comparison trees (e.g. a simulated truth vs. inferred trees)."""
    from .tools.treedist import compare_trees

    def _read_trees(path):
        with open(path) as f:
            return [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]

    ref = _read_trees(args.reference)
    if len(ref) != 1:
        parser.error(f"{args.reference} must contain exactly one reference tree (found {len(ref)})")
    reference = read_newick(ref[0])

    estimates = _read_trees(args.estimate)
    if not estimates:
        raise ValueError(f"no trees found in {args.estimate}")

    rows = []
    for i, nwk in enumerate(estimates, 1):
        try:
            c = compare_trees(reference, read_newick(nwk), quartet=not args.no_quartet,
                              max_leaves=args.max_leaves, branch_score_order=args.branch_order)
        except ValueError as e:
            parser.error(f"tree {i} in {args.estimate}: {e}")
        rows.append(_treedist_row(str(i) if len(estimates) > 1 else "1", c))

    header = "\t".join(_TREEDIST_COLS)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "Tree_distances.tsv")
        with open(path, "w") as f:
            f.write(header + "\n" + "\n".join(rows) + "\n")
        print(f"wrote {path} ({len(rows)} comparison(s))")
    else:
        print(header)
        for r in rows:
            print(r)
    return 0


_RECONACC_COLS = ("family", "n_nodes", "event_acc", "mapping_acc", "joint_acc",
                  "transfers", "transfers_recovered")


def _run_tools_recon_accuracy(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools recon-accuracy`` — node-by-node accuracy of an inferred reconciliation
    against a true one, per family (paired by line) and pooled over all families."""
    from .tools.recon_accuracy import reconciliation_accuracy

    def _read(path):
        with open(path) as f:
            return [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]

    truth, inferred = _read(args.truth), _read(args.inferred)
    if len(truth) != len(inferred):
        parser.error(f"--truth has {len(truth)} tree(s) but --inferred has {len(inferred)}; "
                     "they are paired line by line and must match")
    if not truth:
        raise ValueError(f"no reconciled trees found in {args.truth}")

    accs = []
    for i, (t, e) in enumerate(zip(truth, inferred), 1):
        try:
            accs.append(reconciliation_accuracy(t, e))
        except ValueError as err:
            parser.error(f"family {i}: {err}")

    rows = []
    for i, a in enumerate(accs, 1):
        rows.append("\t".join((
            str(i), str(a.n_nodes), f"{a.event_accuracy:.6f}", f"{a.mapping_accuracy:.6f}",
            f"{a.joint_accuracy:.6f}", str(a.transfer.n_true), str(a.transfer.both_correct),
        )))

    # pooled (micro-averaged over all nodes) summary
    N = sum(a.n_nodes for a in accs)
    ev = sum(round(a.event_accuracy * a.n_nodes) for a in accs)
    mp = sum(round(a.mapping_accuracy * a.n_nodes) for a in accs)
    jt = sum(round(a.joint_accuracy * a.n_nodes) for a in accs)
    nT = sum(a.transfer.n_true for a in accs)
    det = sum(a.transfer.detected for a in accs)
    both = sum(a.transfer.both_correct for a in accs)
    pooled = (
        f"# pooled over {len(accs)} family(ies), {N} node(s): "
        f"event_acc={ev / N:.4f} mapping_acc={mp / N:.4f} joint_acc={jt / N:.4f}"
        if N else "# pooled: no internal nodes to score"
    )
    if nT:
        pooled += f" | transfers: {det}/{nT} detected, {both}/{nT} donor+recipient recovered"

    header = "\t".join(_RECONACC_COLS)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "Reconciliation_accuracy.tsv")
        with open(path, "w") as f:
            f.write(header + "\n" + "\n".join(rows) + "\n" + pooled + "\n")
        print(f"wrote {path} ({len(rows)} family(ies))")
        print(pooled)
    else:
        print(header)
        for r in rows:
            print(r)
        print(pooled)
    return 0


def _run_tools_parse(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools parse`` — read an external reconciliation run (ALE or AleRax) and print a
    summary (rates, log-likelihood, top transfers); with -o, also write the tables as TSV."""
    try:                                       # optional 'reconparser' extra (pandas)
        from .tools.reconparser import ALEParser, AleRaxRun
    except ImportError as e:
        raise RuntimeError(str(e)) from e

    tool = args.tool
    if tool == "auto":
        tool = "alerax" if os.path.isdir(args.path) else "ale"

    if args.out:
        os.makedirs(args.out, exist_ok=True)

    if tool == "ale":
        p = ALEParser(args.path)
        present = [k for k, ok in p.files_exist().items() if ok]
        if not present:
            raise FileNotFoundError(
                f"no ALE output files found for base path {args.path!r} "
                "(expected .ucons_tree / .uTs / .uml_rec) — is --tool right?")
        print(f"ALE reconciliation: {p.base_path}")
        print(f"  files present: {', '.join(present)}")
        try:
            print(f"  log-likelihood: {p.get_log_likelihood():.6f}")
        except (FileNotFoundError, ValueError):
            pass
        try:
            r = p.get_ml_rates()
            print(f"  ML rates:  D={r['duplications']:.4g}  "
                  f"T={r['transfers']:.4g}  L={r['losses']:.4g}")
        except (FileNotFoundError, ValueError):
            pass
        try:
            s = p.get_summary_statistics()
            print(f"  total events:  D={s['total_duplications']:g}  "
                  f"T={s['total_transfers']:g}  L={s['total_losses']:g}  "
                  f"S={s['total_speciations']:g}")
        except (FileNotFoundError, ValueError):
            pass
        try:
            tr = p.get_transfers()
            print(f"  transfers: {len(tr)} edge(s)"
                  + (f" (top {args.top} by frequency)" if len(tr) else ""))
            for _, row in tr.nlargest(args.top, "freq").iterrows():
                print(f"     {row['from']} -> {row['to']}   {row['freq']:.3f}")
            if args.out:
                tp = os.path.join(args.out, "ale_transfers.tsv")
                tr.to_csv(tp, sep="\t", index=False)
                print(f"  wrote {tp}")
        except (FileNotFoundError, ValueError):
            pass
        if args.out:
            try:
                bs = p.get_branch_statistics()
                bp = os.path.join(args.out, "ale_branch_statistics.tsv")
                bs.to_csv(bp, sep="\t", index=False)
                print(f"  wrote {bp}")
            except (FileNotFoundError, ValueError):
                pass
        return 0

    # tool == "alerax"
    run = AleRaxRun(args.path)                  # raises NotADirectoryError on a non-dir path
    print(f"AleRax run: {run.output_dir}")
    try:
        info = run.get_run_info()
        bits = [f"version: {info.get('version', '?')}"]
        if "num_families" in info:
            bits.append(f"families: {info['num_families']}")
        if "num_species" in info:
            bits.append(f"species: {info['num_species']}")
        print("  " + "   ".join(bits))
    except (FileNotFoundError, ValueError):
        pass
    try:
        print(f"  total log-likelihood: {run.get_total_log_likelihood():.6f}")
    except (FileNotFoundError, ValueError):
        pass
    try:
        tr = run.get_transfers()
        score = "score" if "score" in tr.columns else tr.columns[-1]
        print(f"  global transfers: {len(tr)} edge(s)"
              + (f" (top {args.top} by {score})" if len(tr) else ""))
        for _, row in tr.nlargest(args.top, score).iterrows():
            print(f"     {row['from']} -> {row['to']}   {row[score]:.3f}")
        if args.out:
            tp = os.path.join(args.out, "alerax_transfers.tsv")
            tr.to_csv(tp, sep="\t", index=False)
            print(f"  wrote {tp}")
    except (FileNotFoundError, ValueError):
        pass
    if args.out:
        try:
            lk = run.get_per_family_likelihoods()
            lp = os.path.join(args.out, "alerax_per_family_likelihoods.tsv")
            lk.to_csv(lp, sep="\t", index=False)
            print(f"  wrote {lp}")
        except (FileNotFoundError, ValueError):
            pass
    return 0


def _run_tools_red(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools red`` — Relative Evolutionary Divergence of every node of a tree."""
    from .tools import relative_evolutionary_divergence

    with open(args.tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.tree} is not a usable tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    red = relative_evolutionary_divergence(tree)
    rows = [(n.name, n.is_leaf(), red[n]) for n in tree.nodes_preorder()]

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "RED.tsv")
        with open(path, "w") as f:
            f.write("node\tis_leaf\tred\n")
            for name, leaf, r in rows:
                f.write(f"{name}\t{leaf}\t{r:.6f}\n")
        print(f"wrote {path} ({len(rows)} node(s))")
    else:
        print("node\tis_leaf\tred")
        for name, leaf, r in rows:
            print(f"{name}\t{leaf}\t{r:.6f}")
    return 0


def _run_tools_export(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools export`` — gene-order study formats from a nucleotide genomes run."""
    from .tools.geneorder_export import breakpoints_tsv, gff_text, posortho_tsv

    if not os.path.isdir(args.genomes_dir):
        parser.error(f"{args.genomes_dir} is not a directory")
    builders = {"breakpoints": ("Breakpoints.tsv", breakpoints_tsv),
                "gff": ("Genes.gff", gff_text),
                "posortho": ("Positional_orthologs.tsv", posortho_tsv)}
    if args.out:
        os.makedirs(args.out, exist_ok=True)
    for fmt in args.formats:
        filename, build = builders[fmt]
        try:
            text = build(args.genomes_dir)
        except FileNotFoundError as e:
            parser.error(str(e))
        if args.out:
            path = os.path.join(args.out, filename)
            with open(path, "w") as f:
                f.write(text)
            print(f"wrote {path} ({max(0, len(text.splitlines()) - 1)} row(s))")
        else:
            print(text, end="")
    return 0


def _run_nucleotides(tree: Tree, args: argparse.Namespace, parts: set) -> str:
    """Simulate nucleotide-resolution genomes (variable-length structural events) along ``tree``.

    Genes are not atomic here — they emerge as **blocks** (maximal intervals with one shared
    history). ``profiles`` writes the emergent block-by-species profile (plus ``blocks.tsv`` and
    the per-leaf ``Mosaics.tsv``); ``trees`` writes the per-block gene trees and their
    reconciliations. Only ``profiles``/``trees`` apply here (the family-model ``events`` /
    ``transfers`` / ``summary`` do not). ``profiles`` alone takes the fast Rust path.
    """
    want = parts & {"profiles", "trees", "reconciliations", "ancestral", "bed", "geneorder"}
    if not want:
        raise ValueError("the nucleotide model writes 'profiles', 'trees', 'reconciliations', "
                         "'ancestral', 'bed' and/or 'geneorder'; --write events/transfers/summary/"
                         "branch_events do not apply to it")
    ancestral = "ancestral" in want
    bed = "bed" in want
    # --n-chromosomes is the unified flag (both models); --initial-chromosomes is a deprecated alias
    # that takes precedence when explicitly given.
    initial_chromosomes = (args.initial_chromosomes if args.initial_chromosomes is not None
                           else args.n_chromosomes)
    if initial_chromosomes < 1:
        raise ValueError("--n-chromosomes must be >= 1")
    args.n_chromosomes = args.initial_chromosomes = initial_chromosomes  # effective value in the log
    # the structural knobs are shared with the ordered level, so their defaults are resolved here
    args.inversion = 0.001 if args.inversion is None else args.inversion
    args.transposition = 0.0 if args.transposition is None else args.transposition
    args.extension = 0.99 if args.extension is None else args.extension
    if args.gff and args.genes:
        raise ValueError("give either --gff or --genes (not both) to set the gene coordinates")
    gff_info = None
    gff_all = None
    root_chromosomes = None
    if args.gff:                              # start from a real genome: length + gene coordinates
        from zombi2.genomes.gff import read_gff, read_gff_all
        if args.gff_seqid:                    # one named sequence -> a single chromosome
            gff_info = read_gff(args.gff, seqid=args.gff_seqid)
            genes = gff_info.genes
            args.root_length = gff_info.length
        else:
            gff_all = read_gff_all(args.gff)  # every sequence becomes its own chromosome
            if len(gff_all) == 1:
                gff_info = gff_all[0]
                genes = gff_info.genes
                args.root_length = gff_info.length
            else:                             # a multi-replicon genome -> heterogeneous chromosomes
                if initial_chromosomes != 1:
                    raise ValueError("--gff with several sequences already sets the chromosomes; "
                                     "drop --n-chromosomes (or --gff-seqid to pick one sequence)")
                root_chromosomes = [(g.length, g.genes) for g in gff_all]
                genes = None
    elif args.genes:
        genes = _read_gene_intervals(args.genes)
    else:
        genes = None
    genic = bool(genes) or root_chromosomes is not None
    if bed and not genic:
        raise ValueError("--write bed annotates genes on each genome, so it needs gene "
                         "coordinates: supply --genes or --gff")
    transfers = TransferModel(replacement=0.0) if genic else None  # homologous repl. is genome-side
    indels = bool(args.insertion or args.deletion)
    chrom_tier = bool(args.fission or args.fusion or args.chromosome_origination
                      or args.chromosome_loss)
    sim_kw = dict(inversion=args.inversion, loss=args.loss, duplication=args.dup,
                  transfer=args.trans, transposition=args.transposition,
                  translocation=args.translocation,
                  origination=args.orig, insertion=args.insertion, deletion=args.deletion,
                  indel_mean_length=args.indel_mean_length, root_length=args.root_length,
                  extension=args.extension, initial_chromosomes=initial_chromosomes, seed=args.seed,
                  fission=args.fission, fusion=args.fusion,
                  chromosome_origination=args.chromosome_origination,
                  chromosome_loss=args.chromosome_loss,
                  gene_intervals=genes, root_chromosomes=root_chromosomes,
                  pseudogenization=args.pseudogenization,
                  replacement=args.replacement, transfers=transfers,
                  retain_internal=ancestral or bed)  # bed annotates every node's genome

    t0 = time.perf_counter()
    # the Python engine is needed for the event log, the genic model, indels, the chromosome tier,
    # translocation, and explicit heterogeneous replicons (the Rust profiles path is single-chromosome).
    if ("trees" in want or "reconciliations" in want or "geneorder" in want or genic or ancestral
            or bed or indels or chrom_tier or args.translocation or root_chromosomes is not None):
        result = simulate_nucleotide_genomes(tree, output="genomes", **sim_kw)
    else:                                     # profiles only -> Rust fast path (Python fallback)
        try:
            result = simulate_nucleotide_genomes(tree, output="profiles", **sim_kw)
        except (ImportError, RuntimeError):
            result = simulate_nucleotide_genomes(tree, output="genomes", **sim_kw)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    if genic:
        _write_genes_table(args.out, result.registry)

    if "profiles" in want:
        _write_blocks_table(args.out, result.blocks)
        block_ids, species, matrix = result.profile_matrix()
        pm = ProfileMatrix([f"block{a}" for a in block_ids], species, matrix)
        if args.sparse:
            with open(os.path.join(args.out, "Profiles_sparse.tsv"), "w") as f:
                f.write(pm.to_coo_tsv())
        else:
            with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
                f.write(pm.to_tsv())
            with open(os.path.join(args.out, "Presence.tsv"), "w") as f:
                f.write(pm.to_tsv(presence=True))
        _write_mosaics(args.out, result)
    if "trees" in want or "reconciliations" in want:
        result.write_reconciliations(args.out)   # Reconciled_complete/extant.nwk + events.tsv
    if "trees" in want:
        _write_block_gene_trees(args.out, result, genic=genic)
        if genic:
            _write_pseudogenizations(args.out, result)
    if ancestral:
        _write_ancestral(args.out, result, tree, args, gff_info, gff_all)
    if bed:
        _write_bed(args.out, result, tree, gff_info, gff_all)
    if "geneorder" in want:
        from zombi2.genomes.simulation import geneorder_events_from_log
        with open(os.path.join(args.out, "Geneorder_events.tsv"), "w") as f:
            f.write(geneorder_events_from_log(result.event_log))
    # karyotype: when the run is multi-chromosome or uses the chromosome tier, surface the layout
    # (Chromosomes.tsv) and the fission/fusion/origination/loss genealogy (Karyotype_trace.tsv).
    # Needs the Python engine (the Rust profiles path is single-chromosome and event-log-free).
    py_engine = any(isinstance(getattr(g, "chromosomes", None), dict)
                    for g in result.leaf_genomes.values())
    multi_chrom = initial_chromosomes > 1 or (root_chromosomes is not None
                                              and len(root_chromosomes) > 1)
    if py_engine and (chrom_tier or multi_chrom):
        _write_nucleotide_karyotype(args.out, result)

    if gff_all is not None and len(gff_all) > 1:
        print(f"  GFF: {len(gff_all)} sequences -> {len(gff_all)} chromosomes "
              f"({', '.join(f'{g.seqid}:{g.length}bp/{len(g.genes)}g' for g in gff_all)})")
    elif gff_info is not None:
        print(f"  GFF {gff_info.seqid}: {gff_info.length} bp, {gff_info.n_features} genes "
              f"-> {len(gff_info.genes)} after trimming ({gff_info.n_trimmed} trimmed, "
              f"{gff_info.n_dropped} dropped as overlapping)")
    extra = f", {len(result.gene_blocks())} genes" if genic else ""
    return (f"wrote [{' '.join(sorted(want))}] (nucleotide{'/genic' if genic else ''}) to "
            f"{args.out}/ ({len(result.leaf_genomes)} tips, {len(result.blocks)} blocks{extra}) "
            f"in {dt:.3g} s")


def _write_nucleotide_karyotype(out: str, result) -> None:
    """Write the karyotype of a multi-chromosome / chromosome-tier nucleotide run.

    ``Chromosomes.tsv`` — per extant leaf, which chromosome each segment sits on and in what order
    (``species chromosome position source start end strand``); ``Karyotype_trace.tsv`` — the
    fission / fusion / origination / loss genealogy (``time event branch parents children``), one
    row per chromosome-tier event (header-only if the karyotype never changed).
    """
    lay = ["species\tchromosome\tposition\tsource\tstart\tend\tstrand"]
    for leaf, genome in sorted(result.leaf_genomes.items(), key=lambda kv: kv[0].name):
        chroms = getattr(genome, "chromosomes", None)
        if not isinstance(chroms, dict):
            continue
        for chrom in chroms.values():
            for pos, s in enumerate(chrom.elements):
                strand = "+" if s.strand >= 0 else "-"
                lay.append(f"{leaf.name}\t{chrom.chrom_id}\t{pos}\t{s.source}\t"
                           f"{s.src_start}\t{s.src_end}\t{strand}")
    with open(os.path.join(out, "Chromosomes.tsv"), "w") as f:
        f.write("\n".join(lay) + "\n")

    kar = ["time\tevent\tbranch\tparents\tchildren"]
    for r in result.event_log.chromosome_records:
        parents = ";".join(str(p) for p in r.parents)
        children = ";".join(str(c) for c in r.children)
        kar.append(f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t{parents}\t{children}")
    with open(os.path.join(out, "Karyotype_trace.tsv"), "w") as f:
        f.write("\n".join(kar) + "\n")


def _write_ancestral(out: str, result, tree, args, gff_info, gff_all=None) -> None:
    """Simulate sequences and write the genome (architecture + gzipped DNA) at every node.

    ``Architecture/<node>.tsv`` — the ordered, oriented gene/intergene mosaic of the node's genome
    (a ``chromosome`` column keeps replicons apart); ``Genomes/<node>.fasta.gz`` — its assembled DNA,
    one FASTA record per chromosome for a multi-chromosome genome (else one record, the whole
    genome); ``Gene_alignments/<gene>.fasta`` — the extant per-gene alignments. The root sequence is
    seeded from ``--genome-fasta`` when given (a multi-record FASTA for a multi-sequence GFF, matched
    to each replicon by sequence name), else drawn at random.
    """
    from zombi2.sequences.models import make_model, GammaRates, read_fasta, write_fasta
    model = make_model(args.subst_model, kappa=args.kappa,
                       freqs=args.base_freqs, rates=args.gtr_rates)
    gamma = GammaRates(args.gamma_shape) if args.gamma_shape else None
    multi = gff_all is not None and len(gff_all) > 1
    root_fasta = None
    if args.genome_fasta:
        fa = read_fasta(args.genome_fasta)
        if multi:
            # a multi-record FASTA: match each replicon's record (by GFF sequence name) to its source
            root_chroms = list(result.node_genomes[tree.root].chromosomes.values())
            root_fasta = {}
            for g, chrom in zip(gff_all, root_chroms):
                if g.seqid not in fa:
                    raise ValueError(f"--genome-fasta has no record {g.seqid!r} "
                                     f"(have: {', '.join(sorted(fa))})")
                dna = fa[g.seqid]
                if len(dna) != g.length:
                    raise ValueError(f"--genome-fasta record {g.seqid} is {len(dna)} bp but the GFF "
                                     f"says {g.length} bp")
                root_fasta[chrom.elements[0].source] = dna
        else:
            seqid = gff_info.seqid if gff_info is not None else None
            root_fasta = fa[seqid] if seqid in fa else next(iter(fa.values()))
            if len(root_fasta) != args.root_length:
                raise ValueError(f"--genome-fasta sequence is {len(root_fasta)} bp but the genome is "
                                 f"{args.root_length} bp; supply the matching chromosome FASTA")
    result.simulate_sequences(model, gamma=gamma, root_fasta=root_fasta,
                              subst_rate=args.subst_rate, seed=args.seed)

    adir = os.path.join(out, "Architecture")
    gdir = os.path.join(out, "Genomes")
    os.makedirs(adir, exist_ok=True)
    os.makedirs(gdir, exist_ok=True)
    for node in tree.nodes_preorder():
        name = node.name
        mosaics = result.node_mosaics(node)      # {chrom_id: [(block_id, strand), ...]}
        seqs = result.node_sequences(node)        # {chrom_id: dna}
        lines = ["chromosome\torder\tblock\tkind\tgene_id\tstrand\tlength"]
        for cid, mosaic in mosaics.items():
            for i, (aid, strand) in enumerate(mosaic):
                a = result._block_by_id[aid]
                lines.append(f"{cid}\t{i}\tblock{aid}\t{a.kind}\t{a.gene_id or '-'}\t"
                             f"{'+' if strand > 0 else '-'}\t{a.length}")
        with open(os.path.join(adir, f"{name}.tsv"), "w") as f:
            f.write("\n".join(lines) + "\n")
        # one FASTA record per chromosome (multi-chromosome), else one record for the whole genome
        records = ({name: next(iter(seqs.values()))} if len(seqs) == 1
                   else {f"{name}_chr{cid}": dna for cid, dna in seqs.items()})
        write_fasta(os.path.join(gdir, f"{name}.fasta.gz"), records, gzip_out=True)

    aln_dir = os.path.join(out, "Gene_alignments")
    os.makedirs(aln_dir, exist_ok=True)
    for gene, aln in result.gene_alignments().items():
        write_fasta(os.path.join(aln_dir, f"{gene}.fasta"), aln)


def _read_gene_intervals(path: str) -> list[tuple]:
    """Read a BED/TSV of gene intervals: ``start end [name]`` per line (0-based half-open).

    Blank lines and ``#`` comments are skipped; a leading ``track``/``chrom``-style header is
    tolerated (any line whose first field is not an integer).
    """
    out: list[tuple] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            try:
                a, b = int(fields[0]), int(fields[1])
            except (ValueError, IndexError):
                continue  # header / non-numeric row
            name = fields[2] if len(fields) > 2 else None
            out.append((a, b, name))
    if not out:
        raise ValueError(f"no gene intervals found in {path!r} (expected 'start end [name]' lines)")
    return out


def _write_bed(out: str, result, tree, gff_info, gff_all=None) -> None:
    """Write BED gene annotations — one BED6 feature per gene on each genome.

    ``genes.bed`` is the root (seed) genome's annotation, using the input sequence name as the
    chromosome (the GFF/FASTA seqid when a real genome was supplied) so it loads against the
    original genome. ``BED/<node>.bed`` is the annotation of every node's genome *after*
    rearrangements — genes at their coordinates on that node's chromosome, whose chromosome name
    matches the corresponding ``Genomes/<node>.fasta.gz`` record (written by ``--write ancestral``).

    Each chromosome is annotated separately, with coordinates that restart at 0 per chromosome — so
    a multi-chromosome genome (a chromosome plus its plasmids) gets one BED contig per replicon,
    named to line up with its FASTA record: ``<seqid>`` at the root (the input names) and
    ``<node>_chr<id>`` at every node (single-chromosome runs keep the plain ``<seqid>`` / ``<node>``).

    Columns are standard BED6: ``chrom  chromStart  chromEnd  name  score  strand`` — 0-based,
    half-open, the same coordinate convention ZOMBI2 uses internally, so no conversion is needed.
    Only gene blocks are emitted (intergenes are the gaps); the score field is a constant 0.

    ``strand`` is the gene's orientation *relative to the root genome* — every gene is ``+`` at the
    root, and an inversion during evolution flips it. It is not a GFF-annotated coding strand: the
    genic model does not track that (``read_gff`` reads only coordinates), so ``genes.bed`` is
    always all ``+``.
    """
    multi = gff_all is not None and len(gff_all) > 1
    # root chromosome id -> input sequence name (for genes.bed, so it loads on the real genome)
    seqid_by_cid = {}
    if multi:
        root_chroms = list(result.node_genomes[tree.root].chromosomes.values())
        seqid_by_cid = {c.chrom_id: g.seqid for g, c in zip(gff_all, root_chroms)}

    def bed_rows(node, chrom_name) -> list[str]:
        """One BED6 row per gene, coordinates restarting at 0 on each chromosome. ``chrom_name(cid)``
        gives the contig name to match this node's FASTA record(s)."""
        rows: list[str] = []
        for cid, mosaic in result.node_mosaics(node).items():
            chrom, offset = chrom_name(cid), 0
            for block_id, strand in mosaic:
                block = result._block_by_id[block_id]
                if block.kind == "gene":
                    rows.append(f"{chrom}\t{offset}\t{offset + block.length}\t"
                                f"{block.gene_id}\t0\t{'+' if strand > 0 else '-'}")
                offset += block.length
        return rows

    def write_bed(path: str, rows: list[str]) -> None:
        with open(path, "w") as f:
            f.write("\n".join(rows) + ("\n" if rows else ""))

    # root (seed) annotation — chromosome named after the input sequence(s)
    root_seqid = gff_info.seqid if gff_info is not None else "root_chromosome"
    root_name = (lambda cid: seqid_by_cid.get(cid, f"chr{cid}")) if multi else (lambda cid: root_seqid)
    write_bed(os.path.join(out, "genes.bed"), bed_rows(tree.root, root_name))

    # every node's genome (ancestral + extant), each contig keyed to its FASTA record id
    bdir = os.path.join(out, "BED")
    os.makedirs(bdir, exist_ok=True)
    for node in tree.nodes_preorder():
        node_name = ((lambda cid, n=node.name: f"{n}_chr{cid}") if multi
                     else (lambda cid, n=node.name: n))
        write_bed(os.path.join(bdir, f"{node.name}.bed"), bed_rows(node, node_name))


def _write_genes_table(out: str, registry) -> None:
    """Write ``genes.tsv`` — the gene annotation (seed genes + any originated novel genes)."""
    lines = ["gene\tsource\tstart\tend\tlength"]
    for source in sorted(registry.genes):
        for gi in registry.genes[source]:
            lines.append(f"{gi.gene_id}\t{gi.source}\t{gi.start}\t{gi.end}\t{gi.length}")
    with open(os.path.join(out, "genes.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_pseudogenizations(out: str, result) -> None:
    """Write ``Pseudogenizations.tsv`` — every gene->intergene state flip (branch, time, lineage)."""
    lines = ["block\tgene\tspecies_branch\ttime\tgene_lineage"]
    for block_id, gene_id, species, t, gid in result.pseudogenizations():
        lines.append(f"block{block_id}\t{gene_id}\t{species}\t{t:.10g}\t{gid}")
    with open(os.path.join(out, "Pseudogenizations.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_blocks_table(out: str, blocks) -> None:
    """Write ``blocks.tsv`` — the emergent gene families (uncut ancestral intervals).

    Carries the ``kind`` (gene/intergene) and ``gene_id`` classification (``-`` for intergene).
    """
    lines = ["block\tsource\tstart\tend\tlength\tkind\tgene_id"]
    for a in blocks:
        lines.append(f"block{a.block_id}\t{a.source}\t{a.start}\t{a.end}\t{a.length}\t"
                     f"{a.kind}\t{a.gene_id or '-'}")
    with open(os.path.join(out, "blocks.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_mosaics(out: str, result) -> None:
    """Write ``Mosaics.tsv`` — each extant genome as an ordered, signed sequence of blocks."""
    lines = ["leaf\tmosaic"]
    for leaf in sorted(result.leaf_genomes, key=lambda n: n.name):
        seq = " ".join(("+" if s > 0 else "-") + f"block{aid}"
                       for aid, s in result.leaf_mosaic(leaf))
        lines.append(f"{leaf.name}\t{seq}")
    with open(os.path.join(out, "Mosaics.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_block_gene_trees(out: str, result, genic: bool = False) -> None:
    """Write per-block trees to ``block<id>_complete.nwk`` / ``_extant.nwk``.

    Plain nucleotide model: everything under ``gene_trees/``. Genic mode: gene blocks under
    ``Gene_trees/`` and intergene blocks under ``Intergene_trees/`` (both tree sets recovered).
    """
    def dump(tdir: str, trees: dict) -> None:
        os.makedirs(tdir, exist_ok=True)
        for block_id, (complete, extant) in trees.items():
            if complete:
                with open(os.path.join(tdir, f"block{block_id}_complete.nwk"), "w") as f:
                    f.write(complete + "\n")
            if extant:
                with open(os.path.join(tdir, f"block{block_id}_extant.nwk"), "w") as f:
                    f.write(extant + "\n")

    if genic:
        dump(os.path.join(out, "Gene_trees"), result.gene_trees())
        dump(os.path.join(out, "Intergene_trees"), result.intergene_trees())
    else:
        dump(os.path.join(out, "gene_trees"), result.block_gene_trees())


def _add_sequence_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("--genomes", required=True, metavar="DIR",
                   help="a prior 'zombi2 genomes' output directory — reads its species_tree.nwk "
                        "and Events_trace.tsv (run genomes with 'trace' in --write)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("per-family speed")
    g.add_argument("--family-speed", type=float, default=0.0, metavar="SIGMA",
                   help="per-gene-family intrinsic substitution speed: each family draws a "
                        "constant multiplier ~ LogNormal(0, SIGMA) (0 = every family the same)")

    g = p.add_argument_group("lineage clock",
                             "the relaxed molecular clock shared across families (chronogram "
                             "-> phylogram). Pick one with --clock; its parameter knobs below")
    g.add_argument("--clock", default=None, metavar="MODEL",
                   choices=["strict", "autocorrelated-lognormal", "uncorrelated-lognormal",
                            "uncorrelated-gamma", "white-noise", "cir", "discrete-bin"],
                   help="relaxed clock model: strict | autocorrelated-lognormal | "
                        "uncorrelated-lognormal | uncorrelated-gamma | white-noise | cir | "
                        "discrete-bin (default: autocorrelated-lognormal via --branch-speed, "
                        "or discrete-bin if --branch-bins is given)")
    g.add_argument("--clock-mean", type=float, default=1.0, metavar="M",
                   help="mean / strict / root rate of the clock (default 1.0)")
    g.add_argument("--clock-sigma", type=float, default=0.5, metavar="SIGMA",
                   help="rate spread for lognormal / white-noise / cir clocks (the volatility; "
                        "default 0.5). 0 = strict for the uncorrelated & autocorrelated lognormal")
    g.add_argument("--clock-shape", type=float, default=3.0, metavar="ALPHA",
                   help="[--clock uncorrelated-gamma] gamma shape; larger = tighter around the "
                        "mean (default 3.0)")
    g.add_argument("--clock-theta", type=float, default=1.0, metavar="THETA",
                   help="[--clock cir] CIR mean-reversion speed (default 1.0)")
    g.add_argument("--branch-speed", type=float, default=0.0, metavar="SIGMA",
                   help="shorthand for the autocorrelated-lognormal clock: drift SIGMA per "
                        "sqrt(time) (0 = strict). Used when --clock is not given")
    g.add_argument("--branch-bins", default=None, metavar="R1,R2,...",
                   help="[--clock discrete-bin] comma-separated ORDERED rate multipliers "
                        "(e.g. 0.25,0.5,1,2,4), a Markov walk between adjacent bins")
    g.add_argument("--branch-switch-rate", type=float, default=1.0, metavar="RATE",
                   help="[--clock discrete-bin] rate of stepping to a neighbouring bin (default 1.0)")
    g.add_argument("--branch-up-bias", type=float, default=0.5, metavar="P",
                   help="[--clock discrete-bin] probability a step goes to the faster neighbour "
                        "(default 0.5 = symmetric walk)")

    g = p.add_argument_group("sequence alignments",
                             "give --subst-model to evolve DNA/protein along the rescaled trees")
    g.add_argument("--subst-model", default=None, metavar="MODEL",
                   help="substitution model to simulate an alignment per family: DNA "
                        "(jc69/k80/hky85/gtr) or protein (lg/wag/jtt/dayhoff/poisson). DNA vs "
                        "protein is auto-detected from the name. Omit to only rescale the trees "
                        "(no sequences)")
    g.add_argument("--seq-length", type=int, default=300, metavar="N",
                   help="alignment length in sites (nt for DNA, aa for protein; default 300); "
                        "ignored when --root-fasta seeds each family's root")
    g.add_argument("--root-fasta", metavar="FILE", default=None,
                   help="FASTA (optionally .gz) of per-family root sequences keyed by family id "
                        "(header = family id); seeds each family's root instead of a random draw. "
                        "Its length overrides --seq-length per family")
    g.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="discrete-Gamma across-site rate heterogeneity shape (default: none)")
    g.add_argument("--kappa", type=float, default=2.0, metavar="K",
                   help="[DNA k80/hky85] transition/transversion ratio (default 2.0)")
    g.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[DNA hky85/gtr] equilibrium base frequencies (default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[DNA gtr] the 6 exchangeabilities (default all 1)")


def _build_lineage_clock(args: argparse.Namespace):
    """Resolve the sequence command's lineage-clock flags to a (Clock, description) pair.

    ``--clock`` selects the model explicitly; without it we fall back to the historical flags
    (``--branch-bins`` -> discrete-bin, otherwise the autocorrelated lognormal at
    ``--branch-speed``), so old command lines keep working.
    """
    mean = args.clock_mean
    model = args.clock
    if model is None:
        model = "discrete-bin" if args.branch_bins else "autocorrelated-lognormal-legacy"

    if model == "strict":
        return StrictClock(mean), f"strict {mean:g}"
    if model == "autocorrelated-lognormal":
        return AutocorrelatedLogNormalClock(args.clock_sigma, root_rate=mean), \
            f"autocorrelated-lognormal sigma={args.clock_sigma:g}"
    if model == "autocorrelated-lognormal-legacy":
        return AutocorrelatedLogNormalClock(args.branch_speed, root_rate=mean), \
            f"autocorrelated-lognormal sigma={args.branch_speed:g}"
    if model == "uncorrelated-lognormal":
        return UncorrelatedLogNormalClock(args.clock_sigma, mean=mean), \
            f"uncorrelated-lognormal sigma={args.clock_sigma:g}"
    if model == "uncorrelated-gamma":
        return UncorrelatedGammaClock(args.clock_shape, mean=mean), \
            f"uncorrelated-gamma shape={args.clock_shape:g}"
    if model == "white-noise":
        return WhiteNoiseClock(args.clock_sigma, mean=mean), \
            f"white-noise sigma={args.clock_sigma:g}"
    if model == "cir":
        return CIRClock(theta=args.clock_theta, sigma=args.clock_sigma, mean=mean), \
            f"cir theta={args.clock_theta:g} sigma={args.clock_sigma:g}"
    if model == "discrete-bin":
        if not args.branch_bins:
            raise ValueError("--clock discrete-bin needs --branch-bins R1,R2,... (ordered rates)")
        bins = [float(x) for x in args.branch_bins.split(",") if x.strip() != ""]
        return RateVariation(bins=bins, switch_rate=args.branch_switch_rate,
                             up_bias=args.branch_up_bias), f"discrete-bin [{args.branch_bins}]"
    raise ValueError(f"unknown --clock {model!r}")


def _run_sequence(args: argparse.Namespace) -> str:
    """Overlay the gene x lineage substitution clock on a prior genomes run's gene trees, and —
    with ``--subst-model`` — simulate a DNA or protein alignment down each rescaled tree.

    Replays the compact ``Events_trace.tsv`` (no re-simulation of gene content), rescales every
    reconciled gene tree from time into substitutions/site, and writes the phylograms plus the
    drawn per-family speeds and per-branch rates. The lineage clock is shared across families
    (``--branch-speed`` lognormal or ``--branch-bins`` discrete-bin); each family draws one
    constant speed (``--family-speed``). When ``--subst-model`` is given, a sequence is evolved
    along each rescaled **extant** gene tree (the rescaled branch lengths ARE the substitutions/
    site) and the leaf alignment is written as ``alignments/<family>.fasta``.
    """
    from zombi2.genomes.profiles import _natkey
    from zombi2.genomes.reconciliation import extant_species_from_records
    from zombi2.sequences.models import (GammaRates, evolve_on_tree, is_protein_model, make_model,
                               read_fasta, write_fasta)
    from zombi2.genomes.simulation import read_events_trace

    if args.family_speed < 0 or args.branch_speed < 0:
        raise ValueError("--family-speed / --branch-speed must be >= 0")
    if args.clock is None and args.branch_speed > 0 and args.branch_bins:
        raise ValueError("--branch-speed (lognormal clock) and --branch-bins (discrete-bin "
                         "clock) are two lineage clocks; give at most one, or select one "
                         "explicitly with --clock")
    lineage_clock, clock_desc = _build_lineage_clock(args)
    model = None
    if args.subst_model:
        model = make_model(args.subst_model, kappa=args.kappa,
                           freqs=args.base_freqs, rates=args.gtr_rates)
    elif args.gamma_shape or args.root_fasta:
        raise ValueError("--gamma-shape / --root-fasta only apply with --subst-model "
                         "(which turns on sequence simulation)")
    gamma = GammaRates(args.gamma_shape) if args.gamma_shape else None

    tree_path = os.path.join(args.genomes, "species_tree.nwk")
    trace_path = os.path.join(args.genomes, "Events_trace.tsv")
    if not os.path.exists(trace_path):
        raise FileNotFoundError(
            f"{trace_path} not found — re-run 'zombi2 genomes' on that tree with 'trace' in "
            f"--write (e.g. --write trace profiles) so the genealogy can be replayed")
    with open(tree_path) as f:
        tree = read_newick(f.read())
    with open(trace_path) as f:
        # pass the tree so a compact (speciation-free) trace is replayed back to a full genealogy
        families = read_events_trace(f.read(), tree)
    gid2species = extant_species_from_records(families, tree)

    family_speed = LogNormal(0.0, args.family_speed) if args.family_speed > 0 else 1.0
    se = SequenceEvolution(lineage=lineage_clock, family_speed=family_speed)

    t0 = time.perf_counter()
    phylo, node_trees = se.scale_families_trees(tree, families, gid2species, seed=args.seed)
    dt = time.perf_counter() - t0

    tdir = os.path.join(args.out, "gene_trees")
    os.makedirs(tdir, exist_ok=True)
    n = 0
    for fam, extant in phylo.extant.items():
        if extant:
            with open(os.path.join(tdir, f"{fam}_extant_subst.nwk"), "w") as f:
                f.write(extant + "\n")
            n += 1
    for fam, complete in phylo.complete.items():
        if complete:
            with open(os.path.join(tdir, f"{fam}_complete_subst.nwk"), "w") as f:
                f.write(complete + "\n")
    with open(os.path.join(args.out, "gene_family_speeds.tsv"), "w") as f:
        f.write("family\tspeed\n")
        for fam, s in sorted(phylo.family_speed.items()):
            f.write(f"{fam}\t{s:.10g}\n")
    with open(os.path.join(args.out, "branch_rates.tsv"), "w") as f:
        f.write("species_branch\trate\n")
        for name, r in phylo.branch_rate.items():
            f.write(f"{name}\t{r:.10g}\n")

    msg = (f"wrote substitution-unit gene trees for {n} families to {args.out}/gene_trees/ "
           f"(clock: {clock_desc}, family-speed {args.family_speed}) in {dt:.3g} s")

    if model is None:
        return msg

    # Evolve a sequence down each rescaled extant gene tree and write the leaf alignment.
    root_seqs: dict = {}
    if args.root_fasta:
        root_seqs = read_fasta(args.root_fasta)
    rng = np.random.default_rng(args.seed)
    kind = "protein" if is_protein_model(args.subst_model) else "DNA"
    aln_dir = os.path.join(args.out, "alignments")
    os.makedirs(aln_dir, exist_ok=True)

    n_aln = 0
    for fam in sorted(node_trees, key=_natkey):
        entry = node_trees[fam]["extant"]
        if entry is None:                     # no survivors -> no alignment
            continue
        root_node, subst = entry
        kw = {}
        if fam in root_seqs:
            kw["root_seq"] = root_seqs[fam]
        else:
            kw["length"] = args.seq_length
        seqs = evolve_on_tree(root_node, subst, model, rng, gamma=gamma, **kw)
        records = {f"{leaf.species}_{leaf.gid}": seqs[leaf.gid]
                   for leaf in _iter_leaves(root_node)}
        if records:
            write_fasta(os.path.join(aln_dir, f"{fam}.fasta"), records)
            n_aln += 1

    return (f"{msg}; simulated {kind} alignments ({model.name}) for {n_aln} families "
            f"to {args.out}/alignments/")


def _iter_leaves(node):
    """Yield the leaves (childless nodes) of a reconciliation ``_Node`` tree, left to right."""
    if not node.children:
        yield node
        return
    for child in node.children:
        yield from _iter_leaves(child)


def _add_tools_args(p: argparse.ArgumentParser) -> None:
    """The ``tools`` command groups analyses that compute on ZOMBI2 outputs (the ``zombi2.tools``
    layer). Each tool is its own sub-subcommand: ``reconcile`` (ALElite likelihood),
    ``simulate`` (its generative twin — sample gene families under the undated model),
    ``treedist`` (tree distances), ``recon-accuracy`` (reconciliation accuracy), ``red`` (RED)
    and ``parse`` (read external ALE / AleRax reconciliation output)."""
    tsub = p.add_subparsers(dest="tools_command", metavar="<tool>", required=True)
    rp = tsub.add_parser(
        "reconcile",
        help="ALE reconciliation log-likelihood of a gene tree given a species tree",
        description=(
            "Compute the ALE reconciliation log-likelihood P(gene tree | species tree, DTL "
            "rates) of one or more gene trees, EVALUATED at the given --dup/--trans/--loss "
            "(ALElite). This is not inference: it scores fixed rates, it does not fit them."
        ),
        usage="zombi2 tools reconcile -g FILE -t FILE --dup D --trans T --loss L [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # log-likelihood of a reconciled gene tree under the faithful dated model",
            "  zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15",
            "",
            "  # compare all three ALE models and save the table into out/",
            "  zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15 --model dated undated reldated -o out/",
        ),
    )
    _add_tools_reconcile_args(rp)

    smp = tsub.add_parser(
        "simulate",
        help="simulate gene families under the ALE undated/reldated model (generative twin of 'reconcile')",
        description=(
            "Forward-simulate gene families under the ALEml_undated / GeneRax UndatedDTL model — "
            "the exact generative twin of 'zombi2 tools reconcile' (simulate here, then score there, "
            "and the rates round-trip). --dup/--trans/--loss are per-branch ODDS (dimensionless, "
            "relative to a speciation), NOT per-unit-time rates, and the species tree needs no dates "
            "(a cladogram is fine; unit branches are assumed). Writes a ground-truth reconciliation "
            "per family (complete + extant annotated Newicks and an S/D/T/L event table) — the same "
            "format 'zombi2 tools recon-accuracy' scores. For a dated, contemporaneous-transfer "
            "forward simulation, use 'zombi2 genomes' instead."
        ),
        usage="zombi2 tools simulate -t FILE --dup D --trans T --loss L [-n N] [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # sample 200 families on a cladogram under undated odds, write ground-truth reconciliations",
            "  zombi2 tools simulate -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.3 -n 200 -o truth/",
            "",
            "  # then score an inferred reconciliation against that truth",
            "  zombi2 tools recon-accuracy -t truth/Reconciled_extant.nwk -i inferred_recon.nwk",
        ),
    )
    _add_tools_simulate_args(smp)

    tp = tsub.add_parser(
        "treedist",
        help="RF, branch-score, and quartet distances between two trees",
        description=(
            "Tree distances between a REFERENCE tree (e.g. a simulated truth) and one or more "
            "comparison trees (e.g. inferred estimates), over their shared leaf set: rooted and "
            "unrooted Robinson-Foulds, the Kuhner-Felsenstein branch score, and the quartet "
            "distance. One output row per comparison tree."
        ),
        usage="zombi2 tools treedist -r FILE -e FILE [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # distances between a true species tree and an inferred one",
            "  zombi2 tools treedist -r true_tree.nwk -e inferred_tree.nwk",
            "",
            "  # score many bootstrap/replicate trees against one reference, saved to out/",
            "  zombi2 tools treedist -r true_tree.nwk -e replicates.nwk -o out/",
        ),
    )
    _add_tools_treedist_args(tp)

    ap = tsub.add_parser(
        "recon-accuracy",
        help="accuracy of an inferred reconciliation against a known (simulated) one",
        description=(
            "Node-by-node accuracy of an INFERRED reconciliation against the TRUE one for the "
            "same gene tree: event-type accuracy and per-class precision/recall, species "
            "(MRCA) mapping accuracy, and transfer donor/recipient recovery. Inputs are ZOMBI2 "
            "annotated reconciled Newicks (as written by 'zombi2 tools simulate'), "
            "one family per line, --truth and --inferred paired by line."
        ),
        usage="zombi2 tools recon-accuracy -t FILE -i FILE [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # score inferred reconciliations against the simulated truth",
            "  zombi2 tools recon-accuracy -t true_recon.nwk -i inferred_recon.nwk",
        ),
    )
    _add_tools_recon_accuracy_args(ap)

    rp = tsub.add_parser(
        "red",
        help="Relative Evolutionary Divergence (RED) of every node of a tree",
        description=(
            "Compute the Relative Evolutionary Divergence (RED, Parks et al. 2018) of every node "
            "of a rooted tree: the root is 0, every leaf is 1, and each internal node sits at its "
            "relative position along the root-to-tip path. RED is invariant to a global rate "
            "rescaling, so on a phylogram it approximates each node's relative age without a "
            "clock — GTDB's rank-normalisation quantity."
        ),
        usage="zombi2 tools red -t FILE [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # RED of every node (a phylogram recovers relative ages; a dated tree gives them exactly)",
            "  zombi2 tools red -t species_tree.nwk -o out/",
        ),
    )
    _add_tools_red_args(rp)

    pp = tsub.add_parser(
        "parse",
        help="parse external reconciliation output (ALE, AleRax) and summarize it",
        description=(
            "Read the output of an established reconciliation program and print a summary — the "
            "ML DTL rates, the log-likelihood, and the top transfers. Understands classic ALE "
            "(.ucons_tree / .uTs / .uml_rec, v0.4 and v1.0) and AleRax run directories (v1.2+); "
            "the tool is auto-detected from the path (a directory is an AleRax run). With -o it "
            "also writes the transfer / per-branch tables as TSV. This is the reconparser interop "
            "bridge — needs the optional extra:  pip install 'zombi2[reconparser]'."
        ),
        usage="zombi2 tools parse PATH [--tool auto|ale|alerax] [--top N] [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # summarize a classic ALE result (base path, without the .uml_rec extension)",
            "  zombi2 tools parse results.ale",
            "",
            "  # summarize an AleRax run directory and save its transfer tables into out/",
            "  zombi2 tools parse alerax_output/ --top 20 -o out/",
        ),
    )
    _add_tools_parse_args(pp)

    xp = tsub.add_parser(
        "export",
        help="export gene-order study formats from a nucleotide 'zombi2 genomes' run",
        description=(
            "Derive gene-order study formats from a nucleotide genomes output directory (the "
            "complement of the fork's zombiExporter). 'breakpoints' (adjacencies broken per tree "
            "edge), 'gff' (every node's genes as one GFF3) and 'posortho' (positional ortholog "
            "sets over the leaves) come from the per-node gene orders in BED/, so the run needs "
            "'bed' in --write. breakpoints / posortho are exact for content-conserving "
            "rearrangements (inversion / transposition); under duplication / loss gene content "
            "changes, so interpret those accordingly. ('dupinfo' and 'ffgc' are planned.)"
        ),
        usage="zombi2 tools export GENOMES_DIR --format {breakpoints,gff,posortho} [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # simulate with the gene-order outputs, then export the broken adjacencies + GFF",
            "  zombi2 genomes -t species_tree.nwk --genome-model nucleotide --genes genes.tsv \\",
            "      --root-length 3000 --inversion 0.01 --transposition 0.005 --write bed geneorder -o run/",
            "  zombi2 tools export run/ --format breakpoints gff posortho -o export/",
        ),
    )
    _add_tools_export_args(xp)


def _add_tools_export_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("genomes_dir", metavar="GENOMES_DIR",
                   help="a 'zombi2 genomes --genome-model nucleotide' output directory")
    g.add_argument("--format", dest="formats", nargs="+", required=True,
                   choices=("breakpoints", "gff", "posortho"), metavar="FORMAT",
                   help="which format(s) to export: breakpoints / gff / posortho (all need "
                        "--write bed)")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write the export file(s) into DIR (default: print to stdout)")


def _add_tools_reconcile_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-g", "--gene-tree", required=True, metavar="FILE",
                   help="Newick file of one or more reconciled gene trees (one per line); tip "
                        "labels '<species>|<gid>', <species> matching a species-tree leaf")
    g.add_argument("-t", "--species-tree", required=True, metavar="FILE",
                   help="dated species-tree Newick (as written by 'zombi2 species')")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Reconciliation_likelihoods.tsv into DIR (default: print to "
                        "stdout — a bare number for one tree and one model, else a table)")

    g = p.add_argument_group("DTL rates")
    g.add_argument("--dup", type=float, default=0.0, metavar="RATE",
                   help="duplication rate (per-unit-time for dated; per-branch odds for undated/reldated)")
    g.add_argument("--trans", type=float, default=0.0, metavar="RATE", help="transfer rate")
    g.add_argument("--loss", type=float, default=0.0, metavar="RATE", help="loss rate")

    g = p.add_argument_group("model")
    g.add_argument("--model", nargs="+", default=["dated"], metavar="MODEL",
                   choices=("dated", "undated", "reldated"),
                   help="ALE model(s) to score with (default: dated). dated = faithful "
                        "time-sliced likelihood (rates per-unit-time); undated = GeneRax "
                        "UndatedDTL (per-branch odds); reldated = time-overlap-constrained undated")
    g.add_argument("--n-steps", type=int, default=100, metavar="N",
                   help="dated model time-grid resolution (sub-steps per slice; default 100)")
    g.add_argument("--origination", choices=("root", "uniform"), default="root", metavar="WHERE",
                   help="where the family originates: 'root' (default; exact for root-seeded "
                        "families) or 'uniform' over branches")


def _add_tools_simulate_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--species-tree", required=True, metavar="FILE",
                   help="species-tree Newick; a cladogram with no branch lengths is fine for the "
                        "undated model (unit branches are assumed) — reldated needs real dates")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Reconciled_complete.nwk, Reconciled_extant.nwk and "
                        "Reconciliation_events.tsv into DIR (default: print a summary to stdout)")

    g = p.add_argument_group("DTL odds (per-branch, relative to a speciation — NOT per-unit-time)")
    g.add_argument("--dup", type=float, default=0.0, metavar="ODDS", help="duplication odds d")
    g.add_argument("--trans", type=float, default=0.0, metavar="ODDS", help="transfer odds t")
    g.add_argument("--loss", type=float, default=0.0, metavar="ODDS", help="loss odds l")

    g = p.add_argument_group("model")
    g.add_argument("--model", choices=("undated", "reldated"), default="undated", metavar="MODEL",
                   help="undated = a transfer may land on any branch (default); reldated = only on "
                        "a branch that overlaps the donor in time (needs a dated tree)")
    g.add_argument("--origination", choices=("root", "uniform"), default="root", metavar="WHERE",
                   help="where each family originates: 'root' (default) or 'uniform' over branches")
    g.add_argument("-n", "--families", type=int, default=100, metavar="N",
                   help="number of families to simulate (default 100)")
    g.add_argument("--seed", type=int, default=None, metavar="INT",
                   help="random seed for a reproducible draw")
    g.add_argument("--max-events", type=int, default=1_000_000, metavar="N",
                   help="per-family event cap; guards against runaway families at supercritical "
                        "odds (default 1000000)")
    g.add_argument("--score", action="store_true",
                   help="also report the joint undated log-likelihood of the simulated survivors "
                        "(with the true extinct count) under the generating odds — a round-trip check")


def _add_tools_treedist_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-r", "--reference", required=True, metavar="FILE",
                   help="Newick file with exactly one reference tree (e.g. a simulated truth)")
    g.add_argument("-e", "--estimate", required=True, metavar="FILE",
                   help="Newick file of one or more comparison trees (one per line); each is "
                        "compared to the reference and must share its leaf label set")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Tree_distances.tsv into DIR (default: print the table to stdout)")

    g = p.add_argument_group("metrics")
    g.add_argument("--no-quartet", action="store_true",
                   help="skip the quartet distance (it is O(n^4) in the number of leaves)")
    g.add_argument("--max-leaves", type=int, default=100, metavar="N",
                   help="quartet-distance guard: skip it above N leaves (default 100); raise to "
                        "force it on larger trees")
    g.add_argument("--branch-order", type=int, choices=(1, 2), default=2, metavar="P",
                   help="branch-score norm: 2 = L2 / Kuhner-Felsenstein (default), 1 = L1")


def _add_tools_recon_accuracy_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--truth", required=True, metavar="FILE",
                   help="annotated reconciled Newick(s) of the TRUE reconciliation, one family "
                        "per line (labels '<species>|<EVENT>', '<donor>|T>recipient', tips "
                        "'<species>|<gid>')")
    g.add_argument("-i", "--inferred", required=True, metavar="FILE",
                   help="annotated reconciled Newick(s) of the INFERRED reconciliation, paired "
                        "with --truth line by line (same gene-tree topology and tip labels)")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Reconciliation_accuracy.tsv into DIR (default: print to stdout)")


def _add_tools_red_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="Newick tree (one tree). Branch lengths are read as-is: pass a phylogram "
                        "(substitutions) to recover relative ages, or a dated tree for exact "
                        "relative ages. Works with the trees 'zombi2 species'/'sequence' write.")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write RED.tsv (node, is_leaf, red) into DIR (default: print the table to stdout)")


def _add_tools_parse_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("path", metavar="PATH",
                   help="the reconciliation output: an ALE base path (e.g. results.ale, or any "
                        "of its .ucons_tree/.uTs/.uml_rec files) or an AleRax run directory")
    g.add_argument("--tool", choices=("auto", "ale", "alerax"), default="auto", metavar="NAME",
                   help="which reconciliation tool produced PATH (default: auto — a directory "
                        "is treated as an AleRax run, anything else as classic ALE)")
    g.add_argument("--top", type=int, default=10, metavar="N",
                   help="how many top transfers (by frequency/score) to print (default: 10)")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="also write the transfer and per-branch/per-family tables as TSV into DIR")


# --------------------------------------------------------------------------- #
# experimental: the zombi2.experimental layer (unstable, opt-in)
# --------------------------------------------------------------------------- #
def _add_experimental_args(p: argparse.ArgumentParser) -> None:
    """The ``experimental`` command groups unstable, not-yet-validated models (the
    ``zombi2.experimental`` layer). Each is a sub-subcommand: ``selection`` (ESM2 codon dN/dS) and
    ``ils`` (multispecies coalescent)."""
    esub = p.add_subparsers(dest="experimental_command", metavar="<model>", required=True)
    sp = esub.add_parser(
        "selection",
        help="language-model (ESM2) codon selection on a real annotated genome (emergent dN/dS)",
        description=(
            "Evolve a real annotated genome down a species tree with protein-language-model "
            "selection on its coding genes. The nucleotide genome model runs the structural "
            "simulation (inversion/duplication/loss/transfer/...); each gene evolves as coding DNA "
            "along its own gene tree under a codon mutation-selection process whose selection comes "
            "from an ESM2 critic (mutation on DNA, selection on the encoded protein -> emergent "
            "dN/dS), while intergenic DNA drifts neutrally. Genomes are reconstructed at every node.\n\n"
            "EXPERIMENTAL: APIs and outputs may change; needs the optional deps "
            "(pip install 'zombi2[selection]': torch, fair-esm, scipy)."
        ),
        usage="zombi2 experimental selection -t FILE --gff FILE --genome-fasta FILE -o DIR [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # evolve a real genome with ESM2 purifying selection on its genes",
            "  zombi2 experimental selection -t species_tree.nwk --gff genome.gff --genome-fasta genome.fna \\",
            "      --beta 1.0 --dup 0.01 --loss 0.01 --seed 1 -o out/",
            "",
            "  # ...or calibrate the selection strength to a target genome-wide dN/dS, with the big ESM2",
            "  zombi2 experimental selection -t species_tree.nwk --gff genome.gff --genome-fasta genome.fna \\",
            "      --target-dnds 0.2 --esm-model esm2_t33_650M_UR50D -o out/",
        ),
    )
    _add_experimental_selection_args(sp)

    sp_ils = esub.add_parser(
        "ils",
        help="incomplete lineage sorting: gene trees under the multispecies coalescent",
        description=(
            "Simulate gene trees under the multispecies coalescent, so gene lineages need not "
            "coalesce at the nodes they pass through; deep coalescence makes gene trees disagree with "
            "the containing tree -- incomplete lineage sorting (ILS). The amount of ILS is set by "
            "--population-size N, in the tree's own time units: it grows with branch_length / N.\n\n"
            "Two modes. Plain: the coalescent inside the species tree -t (single-copy orthologs). "
            "DTL + ILS: add --events-trace from a 'zombi2 genomes' run and the coalescent runs inside "
            "each gene family's locus tree (duplications/transfers/losses), one gene tree per family; "
            "a duplication's new copy, a transferred copy and the family origination are single-copy "
            "foundings (bounded coalescent), while speciations allow deep coalescence.\n\n"
            "EXPERIMENTAL: APIs and outputs may change. Pure numpy (no optional dependencies)."
        ),
        usage="zombi2 experimental ils -t FILE -N POP [--events-trace FILE] [-n R] [-k C] -o DIR",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # plain ILS: 1000 gene trees under the MSC on a species tree, one allele per species",
            "  zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 --seed 1 -o out/",
            "",
            "  # DTL + ILS: a coalescent gene tree per family from a 'genomes' run (write it with --write trace)",
            "  zombi2 experimental ils -t species_tree.nwk --events-trace run/Events_trace.tsv -N 0.5 -o out/",
        ),
    )
    _add_experimental_ils_args(sp_ils)


def _add_experimental_selection_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="dated species-tree Newick (as written by 'zombi2 species')")
    g.add_argument("--gff", required=True, metavar="FILE",
                   help="GFF3 with CDS features (keeps strand + reading frame); defines the coding "
                        "genes. Single-exon CDS only")
    g.add_argument("--gff-seqid", default=None, metavar="ID",
                   help="which sequence/contig of the GFF (and FASTA) to use (default: the GFF's "
                        "sole sequence; required if it has several)")
    g.add_argument("--genome-fasta", required=True, metavar="FILE",
                   help="the root genome FASTA (optionally .gz); its length sets the chromosome "
                        "length and it seeds the root that then evolves")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")
    g.add_argument("--seed", type=int, default=None, metavar="N", help="RNG seed for reproducibility")

    g = p.add_argument_group("selection (the language-model critic)")
    g.add_argument("--critic", choices=["esm2"], default="esm2", metavar="NAME",
                   help="the protein-language-model critic (default: esm2; the Critic API is "
                        "pluggable from Python for other models)")
    g.add_argument("--esm-model", default="esm2_t6_8M_UR50D", metavar="NAME",
                   help="ESM2 model: small default esm2_t6_8M_UR50D (8M params); go big with e.g. "
                        "esm2_t33_650M_UR50D (650M, GPU recommended)")
    sel = g.add_mutually_exclusive_group()
    sel.add_argument("--beta", type=float, default=None, metavar="B",
                     help="selection strength (>= 0; 0 = neutral). Default 1.0 unless --target-dnds "
                          "is given. Larger = stronger purifying selection (lower dN/dS)")
    sel.add_argument("--target-dnds", type=float, default=None, metavar="W",
                     help="instead of --beta, calibrate beta so the genome-wide expected dN/dS is "
                          "about W (in (0, 1)); measured on the root proteins")

    g = p.add_argument_group("mutation model (nucleotide; codon backbone + intergene)")
    g.add_argument("--subst-model", default="hky85", metavar="MODEL",
                   choices=["jc69", "k80", "hky85", "gtr"],
                   help="nucleotide substitution model: jc69 | k80 | hky85 | gtr (default hky85)")
    g.add_argument("--kappa", type=float, default=2.0, metavar="K",
                   help="[k80/hky85] transition/transversion ratio (default 2.0)")
    g.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[hky85/gtr] equilibrium base frequencies (default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[gtr] the 6 exchangeabilities (default all 1)")
    g.add_argument("--subst-rate", type=float, default=1.0, metavar="R",
                   help="overall divergence scale: neutral substitutions/site at the root (default 1.0)")
    g.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="discrete-Gamma across-site rate heterogeneity for INTERGENE blocks "
                        "(coding-block heterogeneity is emergent from selection; default: none)")

    g = p.add_argument_group("genome structural events (per-nucleotide rates)")
    g.add_argument("--inversion", type=float, default=0.001, metavar="R",
                   help="inversion rate (default 0.001)")
    g.add_argument("--duplication", "--dup", type=float, default=0.0, dest="duplication",
                   metavar="R", help="segmental duplication rate (default 0)")
    g.add_argument("--loss", type=float, default=0.0, metavar="R",
                   help="loss / deletion rate (default 0)")
    g.add_argument("--transfer", "--trans", type=float, default=0.0, dest="transfer",
                   metavar="R", help="transfer rate (default 0)")
    g.add_argument("--transposition", type=float, default=0.0, metavar="R",
                   help="transposition rate (default 0)")
    g.add_argument("--origination", "--orig", type=float, default=0.0, dest="origination",
                   metavar="R", help="per-branch novel-gene origination rate (default 0)")
    g.add_argument("--pseudogenization", type=float, default=0.0, metavar="P",
                   help="probability a loss demotes a gene to intergene rather than deleting it "
                        "(default 0). Note: pseudogenized lineages currently stay under selection")


def _add_experimental_ils_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="species-tree Newick (as written by 'zombi2 species') -- the coalescent "
                        "container, and (with --events-trace) the frame the locus trees live in")
    g.add_argument("--events-trace", default=None, metavar="FILE", dest="events_trace",
                   help="a 'genomes' run's Events_trace.tsv (write it with 'zombi2 genomes ... "
                        "--write trace'). When given, run DTL + ILS: the coalescent within each gene "
                        "family's locus tree, one gene tree per family. Without it, plain species-tree ILS")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")
    g.add_argument("--seed", type=int, default=None, metavar="N", help="RNG seed for reproducibility")

    g = p.add_argument_group("coalescent")
    g.add_argument("-N", "--population-size", type=float, required=True, metavar="POP",
                   dest="population_size",
                   help="effective population size in the tree's time units: pairwise coalescence "
                        "rate 1/POP per unit time. Larger POP => more ILS (governed by branch / POP)")
    g.add_argument("-n", "--replicates", type=int, default=1, metavar="R",
                   help="independent gene trees to draw (default 1); with --events-trace, per family")
    g.add_argument("-k", "--samples", type=int, default=1, metavar="C",
                   help="alleles sampled per species tip / per extant gene copy (default 1)")


def _selection_cds_protein(genome: str, c, translate, reverse_complement):
    """The clean 5'->3' protein of one CDS, or ``None`` if it is out-of-frame / has an internal stop /
    is not ACGT (so calibration only sees translatable coding sequence)."""
    sub = genome[c.start:c.end]
    coding = sub if c.strand == 1 else reverse_complement(sub)
    if c.phase != 0 or len(coding) % 3 or any(ch not in "ACGT" for ch in coding):
        return None
    if len(coding) >= 3 and translate(coding[-3:]) == "*":
        coding = coding[:-3]
    if not coding:
        return None
    prot = translate(coding)
    return None if "*" in prot else prot


def _gff_contigs(path: str) -> set:
    """The set of sequence ids carrying a CDS feature in a GFF3 (to pair the GFF with the FASTA)."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs: set = set()
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) >= 3 and f[2] == "CDS":
                contigs.add(f[0])
    return contigs


def _calibrate_beta_genomewide(critic, proteins, target_dnds, nuc_model, *,
                               hi: float = 64.0, tol: float = 1e-3, max_iter: int = 48) -> float:
    """Find one ``beta`` whose LENGTH-WEIGHTED-MEAN expected dN/dS over ``proteins`` is ~ ``target_dnds``.

    Each protein is profiled by the critic **once** (bounded by the CDS length, so tractable), then dN/dS
    is analytic per beta — unlike a single genome-length concatenation, which would blow up the critic's
    O(L^2) attention. dN/dS is monotone decreasing in beta, so a bisection on ``[0, hi]`` converges.
    """
    from zombi2.experimental.codon_selection import CodonSelection
    from zombi2.experimental.selection import FixedProfileCritic
    if not 0.0 < target_dnds < 1.0:
        raise ValueError(f"--target-dnds must be in (0, 1), got {target_dnds}")
    # one critic call per CDS, then a reusable analytic model per CDS (beta is set live in the loop)
    models = [(CodonSelection(FixedProfileCritic(critic.profile(p)), beta=1.0, nuc_model=nuc_model), p)
              for p in proteins]
    total = float(sum(len(p) for _, p in models))

    def omega(b: float) -> float:
        acc = 0.0
        for sel, p in models:
            sel.beta = b
            acc += len(p) * sel.dnds(p)
        return acc / total

    if omega(hi) > target_dnds:
        raise ValueError(f"--target-dnds {target_dnds} needs beta > {hi}; choose a less extreme target")
    lo, high = 0.0, float(hi)
    for _ in range(max_iter):
        mid = 0.5 * (lo + high)
        w = omega(mid)
        if w > 0.0 and abs(w - target_dnds) <= tol:
            return mid
        lo, high = (mid, high) if w > target_dnds else (lo, mid)
    mid = 0.5 * (lo + high)
    if abs(omega(mid) - target_dnds) > tol:
        raise ValueError(f"--target-dnds calibration did not converge within {max_iter} iterations; "
                         "try a less extreme target")
    return mid


def _run_experimental_selection(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    print("zombi2: 'experimental selection' is unstable — APIs and outputs may change "
          "(zombi2.experimental).", file=sys.stderr)
    import importlib.util
    missing = [m for m in ("torch", "esm", "scipy") if importlib.util.find_spec(m) is None]
    if missing:                                        # probe the ACTUAL optional deps (all lazy inside)
        raise RuntimeError(f"experimental selection needs optional dependencies {missing}; install them "
                           "with:  pip install 'zombi2[selection]'  (torch, fair-esm, scipy)")
    from zombi2.experimental import read_cds_gff, simulate_nucleotide_selection
    from zombi2.experimental.codon_selection import translate
    from zombi2.experimental.selection import ESM2Critic
    from zombi2.sequences.models import GammaRates, make_model, read_fasta, reverse_complement

    with open(args.tree) as f:
        tree = read_newick(f.read())
    fa = read_fasta(args.genome_fasta)
    # pick ONE sequence and use the same id for both the GFF and the FASTA, so CDS coordinates can
    # never be silently applied to the wrong contig
    if args.gff_seqid is not None:
        seqid = args.gff_seqid
    else:
        contigs = _gff_contigs(args.gff)
        if len(contigs) != 1:
            raise ValueError(f"the GFF spans {len(contigs)} sequences {sorted(contigs)}; "
                             "pass --gff-seqid to pick one")
        seqid = next(iter(contigs))
    if seqid not in fa:
        raise ValueError(f"the GFF is annotated on {seqid!r} but --genome-fasta has no such sequence "
                         f"(have: {', '.join(fa)}); supply the matching FASTA or --gff-seqid")
    genome = fa[seqid].upper()
    cds = read_cds_gff(args.gff, seqid=seqid)
    if not cds:
        raise ValueError(f"no CDS features found for {seqid!r} in {args.gff!r}")

    model = make_model(args.subst_model, kappa=args.kappa, freqs=args.base_freqs, rates=args.gtr_rates)
    gamma = GammaRates(args.gamma_shape) if args.gamma_shape else None
    critic = ESM2Critic(args.esm_model)                # args.critic == "esm2" (the only choice)

    if args.target_dnds is not None:
        proteins = [p for p in (_selection_cds_protein(genome, c, translate, reverse_complement)
                                for c in cds) if p]
        if not proteins:
            raise ValueError("no cleanly-translatable CDS to calibrate --target-dnds on")
        beta = _calibrate_beta_genomewide(critic, proteins, args.target_dnds, model)
        print(f"zombi2: calibrated beta = {beta:.4g} for a genome-wide dN/dS ~ {args.target_dnds} "
              f"(length-weighted over {len(proteins)} CDS)", file=sys.stderr)
    else:
        beta = 1.0 if args.beta is None else args.beta

    result, report = simulate_nucleotide_selection(
        tree, genome, cds, critic=critic, beta=beta, nuc_model=model, gamma=gamma,
        subst_rate=args.subst_rate, seed=args.seed,
        inversion=args.inversion, duplication=args.duplication, loss=args.loss,
        transfer=args.transfer, transposition=args.transposition, origination=args.origination,
        pseudogenization=args.pseudogenization)

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick())
    _write_selection_outputs(args.out, result, tree, report, beta)

    n_nodes = sum(1 for _ in tree.nodes_preorder())
    summary = (f"experimental selection: {report.n_selected}/{report.n_gene_blocks} gene blocks under "
               f"selection ({report.n_neutral_fallback} fell back to neutral), "
               f"{report.n_intergene} intergene blocks; beta={beta:.4g}. "
               f"Genomes for {n_nodes} nodes -> {args.out}/")
    print(summary)
    _write_params_log(os.path.join(args.out, "selection.log"), args, summary)
    return 0


def _run_experimental_ils(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    print("zombi2: 'experimental ils' is unstable — APIs and outputs may change "
          "(zombi2.experimental).", file=sys.stderr)
    from zombi2.experimental.ils import MultispeciesCoalescent, is_concordant

    if args.replicates < 1:
        raise ValueError("--replicates must be >= 1")
    if args.samples < 1:
        raise ValueError("--samples must be >= 1")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    msc = MultispeciesCoalescent(population_size=args.population_size)
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick())

    if args.events_trace:                          # DTL + ILS: the coalescent within each locus tree
        from zombi2.genomes.reconciliation import extant_species_from_records
        from zombi2.genomes.simulation import read_events_trace
        with open(args.events_trace) as f:
            families = read_events_trace(f.read(), tree)
        if not families:
            raise ValueError(f"no gene-family events in {args.events_trace!r}")
        gid2species = extant_species_from_records(families, tree)
        fam_trees = msc._family_trees(families, gid2species, tree.total_age,
                                      args.samples, args.replicates, rng)
        gdir = os.path.join(args.out, "gene_trees")
        os.makedirs(gdir, exist_ok=True)
        for family, trees in sorted(fam_trees.items()):
            with open(os.path.join(gdir, f"{family}.nwk"), "w") as fh:
                for t in trees:
                    fh.write(t.to_newick(include_internal_names=False) + "\n")
        reps = "" if args.replicates == 1 else f", x{args.replicates} reps"
        summary = (f"experimental ils (DTL + ILS): coalescent gene trees for {len(fam_trees)} of "
                   f"{len(families)} families -> {args.out}/gene_trees/ (N={args.population_size:g}{reps})")
        print(summary)
        _write_params_log(os.path.join(args.out, "ils.log"), args, summary)
        return 0

    genes = msc.sample_gene_trees(tree, args.replicates, samples=args.samples, rng=rng)
    with open(os.path.join(args.out, "gene_trees.nwk"), "w") as f:
        for g in genes:
            f.write(g.to_newick(include_internal_names=False) + "\n")
    copies = "1 copy/species" if args.samples == 1 else f"{args.samples} copies/species"
    summary = (f"experimental ils: {len(genes)} gene tree(s) under the multispecies coalescent "
               f"(N={args.population_size:g}, {copies}) -> {args.out}/gene_trees.nwk")
    if args.samples == 1:
        conc = sum(is_concordant(g, tree) for g in genes) / len(genes)
        summary += f"; {conc:.1%} match the species-tree topology"
    print(summary)
    _write_params_log(os.path.join(args.out, "ils.log"), args, summary)
    return 0


def _write_selection_outputs(out: str, result, tree, report, beta: float) -> None:
    """Write the per-node genomes + architecture, the extant gene alignments, and the selection report."""
    from zombi2.sequences.models import write_fasta
    adir = os.path.join(out, "Architecture")
    gdir = os.path.join(out, "Genomes")
    os.makedirs(adir, exist_ok=True)
    os.makedirs(gdir, exist_ok=True)
    for node in tree.nodes_preorder():
        name = node.name
        lines = ["order\tblock\tkind\tgene_id\tstrand\tlength"]
        for i, (aid, strand) in enumerate(result.node_mosaic(node)):
            a = result._block_by_id[aid]
            lines.append(f"{i}\tblock{aid}\t{a.kind}\t{a.gene_id or '-'}\t"
                         f"{'+' if strand > 0 else '-'}\t{a.length}")
        with open(os.path.join(adir, f"{name}.tsv"), "w") as f:
            f.write("\n".join(lines) + "\n")
        write_fasta(os.path.join(gdir, f"{name}.fasta.gz"), {name: result.node_sequence(node)},
                    gzip_out=True)

    aln_dir = os.path.join(out, "Gene_alignments")
    os.makedirs(aln_dir, exist_ok=True)
    for gene, aln in result.gene_alignments().items():
        write_fasta(os.path.join(aln_dir, f"{gene}.fasta"), aln)

    with open(os.path.join(out, "Selection_report.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        for k in ("n_blocks", "n_gene_blocks", "n_selected", "n_neutral_fallback",
                  "n_intergene", "n_empty"):
            f.write(f"{k}\t{getattr(report, k)}\n")
        f.write(f"beta\t{beta:.6g}\n")
    if report.fallbacks:
        with open(os.path.join(out, "Selection_fallbacks.tsv"), "w") as f:
            f.write("block_id\tgene_id\treason\n")
            for block_id, gene_id, reason in report.fallbacks:
                f.write(f"{block_id}\t{gene_id or '-'}\t{reason}\n")


def _add_subcommand(sub, name: str, help: str, description: str, usage: str, adder,
                    epilog: str | None = None):
    """Register a subcommand with the house-style formatter and a hand-written compact usage.

    The command list itself is curated (grouped by theme) in the top-level description, so the
    per-command ``help`` is suppressed from argparse's auto listing to avoid a duplicate dump.
    ``epilog`` (built with :func:`_examples`) adds a worked-example block below the options.
    """
    p = sub.add_parser(name, help=help, description=description, usage=usage,
                       epilog=epilog, formatter_class=ZombiHelpFormatter)
    adder(p)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="zombi2", description=_banner() + "\n\n" + _DESCRIPTION,
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # 1. a dated species tree (20 extant tips)",
            "  zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o out/",
            "",
            "  # 2. gene families along it",
            "  zombi2 genomes -t out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/",
            "",
            "  # 3. a trait (or sequences, or coupled processes) on the same tree",
            "  zombi2 trait -t out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/",
            "",
            "Run 'zombi2 <command> -h' for a command's options and its own examples.",
        ),
    )
    parser.add_argument("--version", action="version", version=f"ZOMBI2 {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)

    _add_subcommand(
        sub, "species", "simulate a dated species tree",
        "Simulate a dated species tree.",
        "zombi2 species -o DIR [--mode MODE] [--diversification PROCESS] [options]",
        _add_species_args,
        epilog=_examples(
            "  # backward birth-death tree, 50 extant tips at age 5",
            "  zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o out/",
            "",
            "  # forward simulation to a fixed age (complete tree, keeps extinct lineages)",
            "  zombi2 species --mode forward --birth 1 --death 0.4 --age 5 --seed 1 -o out/",
        ))

    _add_subcommand(
        sub, "genomes", "evolve gene families along a species tree",
        "Evolve gene families along a species tree.",
        "zombi2 genomes -t FILE -o DIR [--genome-model LEVEL] [--rate-model MODEL] "
        "[--write PART ...] [options]",
        _add_rate_args,
        epilog=_examples(
            "  # DTL gene families with a full event log and gene trees",
            "  zombi2 genomes -t out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/",
            "",
            "  # counts-only profiles (scales to very large trees)",
            "  zombi2 genomes -t out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --write profiles --seed 42 -o out/",
        ))

    _add_subcommand(
        sub, "trait", "evolve a phenotypic trait along a given species tree",
        "Evolve a phenotypic trait along a species tree, writing tip and ancestral values.",
        "zombi2 trait -t FILE -o DIR [--model MODEL] [options]", _add_trait_args,
        epilog=_examples(
            "  # Ornstein-Uhlenbeck continuous trait",
            "  zombi2 trait -t out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/",
            "",
            "  # 3-state discrete Mk trait, 20 replicates",
            "  zombi2 trait -t out/species_tree.nwk --model mk --states 3 --replicates 20 --seed 1 -o out/",
        ))

    _add_subcommand(
        sub, "coevolve", "co-evolve coupled processes (--couple driver:target)",
        "Co-evolve coupled processes over {species, traits, genes} — pick directed edges with "
        "--couple (e.g. traits:species = SSE, traits:genes = trait-conditioned gene families).",
        "zombi2 coevolve -o DIR --couple DRIVER:TARGET [-t FILE] [--age T|--tips N] [options]",
        _add_coevolve_mode_args,
        epilog=_examples(
            "  # trait-conditioned gene families (loss/gain depends on a simulated trait)",
            "  zombi2 coevolve --couple traits:genes -t out/species_tree.nwk --trait-model mk --states 2 --trait-center --responsive 0.3 --effect-loss 3 --seed 1 -o out/",
            "",
            "  # trait-dependent diversification (BiSSE), grows the tree",
            "  zombi2 coevolve --couple traits:species --sse-model bisse --tips 50 --seed 1 -o out/",
        ))

    _add_subcommand(
        sub, "sequence", "simulate DNA/protein alignments along a genomes run's gene trees",
        "Rescale a 'genomes' run's gene trees from time into substitutions/site under a "
        "gene × lineage clock, then (with --subst-model) simulate a DNA or protein sequence "
        "alignment along each rescaled gene tree.",
        "zombi2 sequence --genomes DIR -o DIR [--subst-model MODEL] "
        "[--clock MODEL [--clock-sigma S]] [options]",
        _add_sequence_args,
        epilog=_examples(
            "  # rescale gene trees into substitutions/site (needs a 'genomes' run done with --write trace)",
            "  zombi2 sequence --genomes out/ --branch-speed 0.4 --family-speed 0.5 --seed 7 -o out/",
            "",
            "  # ...and also simulate DNA alignments under HKY85",
            "  zombi2 sequence --genomes out/ --subst-model hky85 --branch-speed 0.4 --seed 7 -o out/",
        ))

    _add_subcommand(
        sub, "tools", "compute on ZOMBI2 outputs (reconcile, treedist, recon-accuracy, red, parse, export)",
        "Analysis tools that compute on ZOMBI2 outputs — the stable analysis complement to the "
        "simulator (the zombi2.tools layer). Each tool is a sub-subcommand; run "
        "'zombi2 tools <tool> -h' for its options.\n\n"
        "Tools\n"
        "  reconcile            ALE reconciliation likelihood of a gene tree (ALElite)\n"
        "  treedist             tree distances (RF, branch-score, quartet, matching) vs a reference\n"
        "  recon-accuracy       accuracy of an inferred reconciliation vs a known one\n"
        "  red                  Relative Evolutionary Divergence of every node (Parks et al. 2018)\n"
        "  parse                read external ALE / AleRax reconciliation output (reconparser)",
        "zombi2 tools <tool> [options]",
        _add_tools_args,
        epilog=_examples(
            "  # ALE reconciliation log-likelihood of a gene tree at given DTL rates",
            "  zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15",
            "",
            "  # Relative Evolutionary Divergence of every node of a tree",
            "  zombi2 tools red -t species_tree.nwk -o out/",
            "",
            "  # summarize an existing ALE / AleRax reconciliation (needs zombi2[reconparser])",
            "  zombi2 tools parse results.ale",
        ))

    _add_subcommand(
        sub, "experimental", "unstable, opt-in models (selection: ESM2 dN/dS; ils: multispecies coalescent)",
        "Experimental, not-yet-validated models (the zombi2.experimental layer) — APIs and outputs "
        "may change. Each is a sub-subcommand; run 'zombi2 experimental <model> -h' for its options.\n\n"
        "Models\n"
        "  selection            language-model (ESM2) codon selection on a real annotated genome\n"
        "  ils                  incomplete lineage sorting (multispecies-coalescent gene trees)",
        "zombi2 experimental <model> [options]",
        _add_experimental_args,
        epilog=_examples(
            "  # evolve a real genome down a species tree with ESM2 purifying selection on its genes",
            "  zombi2 experimental selection -t species_tree.nwk --gff genome.gff --genome-fasta genome.fna --beta 1 -o out/",
            "",
            "  # draw 1000 gene trees under the multispecies coalescent (incomplete lineage sorting)",
            "  zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 -o out/",
        ))

    args = parser.parse_args(argv)              # the banner shows on --help only, not on every run
    try:
        return _dispatch(args, parser)
    except (ValueError, RuntimeError, FileNotFoundError, OSError) as e:
        # Report expected failures as a clean one-line error, never a traceback.
        print(f"zombi2: error: {e}", file=sys.stderr)
        return 1


def _dispatch(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.command == "species":
        if args.model == "forward" and args.ghosts:
            parser.error("--ghosts un-prunes a reconstructed (backward) tree; forward trees "
                         "already include extinct lineages")
        model = _build_species_model(args, parser)
        common = dict(age_type=args.age_type, max_attempts=args.max_attempts,
                      max_lineages=args.max_lineages, seed=args.seed)

        t0 = time.perf_counter()
        if args.model == "backward":
            n_tips = args.tips if args.tips is not None else 50
            age = args.age if args.age is not None else 1.0
            tree = simulate_species_tree(model, n_tips=n_tips, age=age,
                                         direction="backward", **common)
            if args.ghosts:
                add_ghost_lineages(tree, model, method=args.ghost_method, seed=args.seed)
        else:  # forward
            if (args.tips is None) == (args.age is None):
                parser.error("forward model needs exactly one of --tips or --age "
                             "(--tips to stop at that many extant species; "
                             "--age to grow for that long)")
            if args.mass_extinction and args.age is None:
                parser.error("--mass-extinction needs --age: its times are ages before a fixed "
                             "present, which --tips (random age) leaves undefined")
            if args.clade_shift and args.age is None:
                parser.error("--clade-shift needs --age: its times are ages before a fixed "
                             "present, which --tips (random age) leaves undefined")
            try:
                tree = simulate_species_tree(model, n_tips=args.tips, age=args.age,
                                             direction="forward", **common)
            except RuntimeError:
                raise RuntimeError(
                    f"forward simulation kept going extinct in {args.max_attempts} attempts. "
                    f"With --death {args.death} vs --birth {args.birth}, most runs die out — "
                    f"lower --death, raise --max-attempts, or use --mode backward.") from None
        dt = time.perf_counter() - t0

        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
            f.write(tree.to_newick() + "\n")               # the complete tree (extinct/ghost tips kept)
        leaves = tree.leaves()
        n_extant = sum(1 for n in leaves if n.is_extant)
        n_unsampled = sum(1 for n in leaves if n.name.startswith("u"))   # ghost tips (u*), from ρ<1
        n_extinct = len(leaves) - n_extant - n_unsampled
        wrote = "species_tree.nwk"
        if n_extant and n_extant < len(leaves):            # dead tips present: also the pruned tree
            with open(os.path.join(args.out, "species_tree_extant.nwk"), "w") as f:
                f.write(prune(tree, keep="extant").to_newick() + "\n")
            wrote += " + species_tree_extant.nwk"
        with open(os.path.join(args.out, "species_nodes.tsv"), "w") as f:   # node metadata table
            f.write("name\ttime\tis_leaf\tis_extant\n")
            for node in tree.nodes():
                is_leaf = not node.children
                f.write(f"{node.name}\t{node.time:.10g}\t{is_leaf}\t{bool(node.is_extant)}\n")
        wrote += " + species_nodes.tsv"
        parts = [f"{n_extant} extant"]
        if n_extinct:
            parts.append(f"{n_extinct} extinct")
        if n_unsampled:
            parts.append(f"{n_unsampled} unsampled")
        summary = " + ".join(parts) + " tips"
        print(f"wrote {args.out}/{wrote} ({summary}) in {dt:.3g} s")
        _write_params_log(os.path.join(args.out, "species_tree.log"), args, summary)
        return 0

    if args.command == "genomes":
        with open(args.tree) as f:
            tree = read_newick(f.read())
        if len(tree.leaves()) < 2:
            parser.error(f"{args.tree} is not a usable species tree — fewer than 2 tips "
                         "(is it a valid Newick file?)")
        summary = _run_genomes(tree, args, parser)
        print(summary)
        _write_params_log(os.path.join(args.out, "genomes.log"), args, summary)
        return 0

    if args.command == "trait":
        summary = _run_trait(args)
        print(summary)
        _write_params_log(os.path.join(args.out, "trait.log"), args, summary)
        return 0

    if args.command == "sequence":
        summary = _run_sequence(args)
        print(summary)
        _write_params_log(os.path.join(args.out, "sequence.log"), args, summary)
        return 0

    if args.command == "coevolve":
        summary = _run_coevolve_mode(args, parser)
        print(summary)
        _write_params_log(os.path.join(args.out, "coevolve.log"), args, summary)
        return 0

    if args.command == "tools":
        if args.tools_command == "reconcile":
            return _run_tools_reconcile(args, parser)
        if args.tools_command == "simulate":
            return _run_tools_simulate(args, parser)
        if args.tools_command == "treedist":
            return _run_tools_treedist(args, parser)
        if args.tools_command == "recon-accuracy":
            return _run_tools_recon_accuracy(args, parser)
        if args.tools_command == "red":
            return _run_tools_red(args, parser)
        if args.tools_command == "export":
            return _run_tools_export(args, parser)
        if args.tools_command == "parse":
            return _run_tools_parse(args, parser)
        parser.error(f"unknown tool {args.tools_command!r}")   # unreachable: subparsers validate

    if args.command == "experimental":
        if args.experimental_command == "selection":
            return _run_experimental_selection(args, parser)
        if args.experimental_command == "ils":
            return _run_experimental_ils(args, parser)
        parser.error(f"unknown experimental model {args.experimental_command!r}")   # unreachable

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
