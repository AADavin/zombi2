"""Command-line interface for ZOMBI2 (``zombi2 species`` / ``genomes`` / ``trait``)."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time

import numpy as np

from . import __version__

from .biogeography import DEC, simulate_biogeography
from .ghosts import add_ghost_lineages
from .matching import match_profiles, match_profiles_smc
from .nucleotide_sim import simulate_nucleotide_genomes
from .profiles import ProfileMatrix
from .distributions import LogNormal
from .rate_variation import RateVariation
from .genome import OrderedGenome
from .rates import PerGenomeRates, SharedRates
from .sequence_evolution import SequenceEvolution
from .simulation import Genomes, simulate_genomes
from .species_model import (
    BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent, EpisodicBirthDeath,
)
from .species_sim import simulate_species_tree
from .sse import BiSSE, MuSSE, QuaSSE, simulate_sse
from .gene_diversification import GeneDiversification, simulate_gene_diversification
from .cladogenetic_genome import CladogeneticGenome, simulate_cladogenetic_genome
from .gene_conditioned_trait import GeneConditionedTrait, simulate_gene_conditioned_trait
from .coupling import pathway_blocks, simulate_coupled
from .trait_coupling import TraitGeneCoupling, simulate_trait_linked_genomes
from ._traits_impl import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst, Mk, ThresholdModel, TraitResult,
    Cladogenesis, simulate_traits,
)
from .transfers import TransferModel
from .tree import Tree, read_newick

_DESCRIPTION = """\
Simulate each level on its own, or couple them into joint models; or run the inverse and
fit rates to data. Run 'zombi2 <command> -h' for a command's options, grouped by model.

Species trees
  species              simulate a dated species tree

Gene families & sequences
  genomes              evolve gene families along a species tree (Newick)
  sequence             simulate DNA/protein alignments along a genomes run's gene trees

Traits & coevolution
  trait                evolve a phenotypic trait along a given species tree
  coevolve             co-evolve coupled processes (--couple driver:target)

Inference
  abc                  fit gene-family rates to an empirical profile (ABC)
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
    g.add_argument("--rate-model", choices=("shared", "per-genome", "coupled"),
                   default="shared", metavar="MODEL",
                   help="rate heterogeneity within the unordered/ordered genome levels: shared: "
                        "same per-copy rates for every family (Rust for unordered; default); "
                        "per-genome: constant per-genome rates, linear growth (Python); coupled: "
                        "non-independent families (Potts/Ising) whose presence/absence co-varies by "
                        "pathway (see the coupled model section; unordered only). Rearrangements "
                        "(--inversion/--transposition on ordered genomes) need shared")
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

    g = p.add_argument_group("output")
    g.add_argument("--write", dest="output", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "ancestral", "all"), default=["profiles", "trees"],
                   help="which output files to write — any of {profiles, trace, trees, events, "
                        "transfers, summary} or 'all' (default: profiles trees). species_tree.nwk "
                        "is always written; 'profiles' alone takes the fast Rust counts-only path; "
                        "'trace' (optionally with 'profiles') writes the compact single-file event "
                        "log Events_trace.tsv near counts-only speed, from which gene trees can be "
                        "reconstructed later on demand. [nucleotide] 'ancestral' simulates DNA and "
                        "reconstructs the genome (architecture + gzipped FASTA) at every node")
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
                        "--score-model, at the --dup/--trans/--loss rates")
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
    g.add_argument("--extension", type=float, default=None, metavar="P",
                   help="geometric event-length parameter (mean length 1/(1-extension)): counted "
                        "in genes for --genome-model ordered (default None = single-gene events), "
                        "in nucleotides for nucleotide (default 0.99)")

    g = p.add_argument_group("nucleotide model", "with --genome-model nucleotide")
    g.add_argument("--initial-chromosomes", type=int, default=None, metavar="N",
                   dest="initial_chromosomes",
                   help="number of root chromosomes seeded at the root, for --genome-model "
                        "nucleotide (default: 1)")
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
                        "separate knob from --extension (default 10)")

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

    g = p.add_argument_group(
        "coupled model", "with --rate-model coupled (Potts/Ising gene-family non-independence); "
        "reuses --trans (HGT) and --orig; supports the full --write set (profiles, trees, "
        "events, ...) over the fixed family panel")
    g.add_argument("--pathways", metavar="SIZES", default="4,4",
                   help="comma-separated pathway block sizes defining the family panel and its "
                        "coupling blocks, e.g. 4,4,3 -> an 11-family panel of three co-occurring "
                        "pathways (default 4,4)")
    g.add_argument("--within", type=float, default=3.0, metavar="J",
                   help="coupling J between families in the SAME pathway (>0 -> co-occurrence; "
                        "default 3.0)")
    g.add_argument("--between", type=float, default=0.0, metavar="J",
                   help="coupling J between families in DIFFERENT pathways (<0 -> mutually "
                        "exclusive rival pathways; default 0 = independent blocks)")
    g.add_argument("--field", type=float, default=2.0, metavar="H",
                   help="intrinsic per-family retention field h (larger -> more universally "
                        "present 'core' families; default 2.0)")
    g.add_argument("--beta", type=float, default=1.0, metavar="B",
                   help="global coupling strength (inverse temperature) scaling the whole field "
                        "(0 -> independent; default 1.0)")
    g.add_argument("--base-loss", type=float, default=1.0, metavar="RATE", dest="base_loss",
                   help="baseline per-family loss rate, the loss at field 0 (default 1.0); "
                        "coupling enters through loss as loss_i = base-loss * exp(-beta * f_i)")
    g.add_argument("--gain-coupling", type=float, default=0.0, metavar="G", dest="gain_coupling",
                   help="gain-side coupling: field-bias HGT establishment so a transferred copy "
                        "sticks preferentially where its pathway partners are present "
                        "(0 -> field-blind gain, default; needs --trans > 0 to have any effect)")


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
                        "exp(effect_gain * trait) (default 0 = field-blind gain, as in the Potts model)")

    g = p.add_argument_group("traits:genes output")
    g.add_argument("--write", dest="output", nargs="+", metavar="PART",
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
    parts = set(Genomes.WRITE_PARTS) if "all" in args.output else set(args.output)
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
    return (f"wrote [{' '.join(sorted(parts))}] + trait to {args.out}/ (trait={trait_desc}, "
            f"panel {coupling.n_families} families, {coupling.n_responsive} responsive, "
            f"{len(tree.extant_leaves())} tips) in {dt:.3g} s")


# ═══════════════════════════════════════════════════════════════════════════════
# coevolve: the directed-coupling umbrella (Phase 1: traits:species = SSE)
# ═══════════════════════════════════════════════════════════════════════════════
_COEVOLVE_NODES = ("species", "traits", "genes")
# every directed edge in the coevolution design (docs/coevolution_models.md); Phase 1
# implements only traits:species.
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
                        "Repeatable; default traits:species. See docs/coevolution_models.md for "
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
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("SSE model", "--couple traits:species (trait drives diversification)")
    g.add_argument("--sse-model", dest="sse_model", choices=("bisse", "musse", "quasse"),
                   default="bisse", metavar="MODEL",
                   help="which state-dependent model drives diversification: bisse (binary trait), "
                        "musse (k-state), quasse (continuous trait) (default: bisse)")
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


def _run_coevolve_mode(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Run the ``coevolve`` umbrella over the six directed edges (``--couple``): ``traits:species``
    (SSE), ``species:traits`` (cladogenetic) and their combination = **ClaSSE**; ``genes:species``
    (key innovations); ``species:genes`` (cladogenetic genome); ``genes:traits`` (gene-conditioned
    trait); and ``traits:genes`` (trait-conditioned genes). Whether the tree is grown (an arrow into
    species) or read from ``-t`` follows the arrows-into-S rule."""
    # --couple accepts both repeated flags and space-separated lists (append + nargs); flatten
    raw = args.couple or [["traits:species"]]
    edges = [e.strip().lower() for group in raw for e in group]
    for e in edges:
        if e not in _COEVOLVE_EDGES:
            parser.error(f"unknown --couple edge {e!r}: expected 'driver:target' over "
                         f"{{{', '.join(_COEVOLVE_NODES)}}} (e.g. traits:species); see "
                         "docs/coevolution_models.md for the full edge set")
    eset = set(edges)
    supported = {"traits:species", "species:traits", "genes:species", "species:genes",
                 "genes:traits", "traits:genes"}
    unsupported = eset - supported
    if unsupported:
        parser.error(f"--couple {', '.join(sorted(unsupported))} is planned but not yet "
                     "implemented; the built edges are traits:species (SSE), species:traits "
                     "(cladogenetic), their combination (ClaSSE), genes:species (key innovations), "
                     "species:genes (cladogenetic genome), genes:traits (gene-conditioned trait), "
                     "and traits:genes (trait-conditioned genes). See docs/coevolution_models.md")

    # traits:genes — a trait conditions a gene-family panel's loss/gain (formerly the standalone
    # 'coevolve-genetrait' command). An overlay edge (no arrow into S), so the tree is an INPUT.
    if "traits:genes" in eset:
        if eset != {"traits:genes"}:
            parser.error("traits:genes runs on its own in this phase; combining it with other "
                         "edges is future work — see docs/coevolution_models.md")
        return _run_traits_genes(args, parser)

    # genes:species — gene content drives diversification (its own forward joint loop, v1 stands
    # alone; combining it with other edges is the full joint model, still on the roadmap).
    if "genes:species" in eset:
        if eset != {"genes:species"}:
            parser.error("genes:species runs on its own in this phase; combining it with other "
                         "edges (the fully joint model) is future work — see docs/coevolution_models.md")
        return _run_genes_species(args, parser)

    # genes:traits — gene content conditions a trait (a modifier gene switches the trait's OU
    # optimum). An overlay edge (no arrow into S), so the tree is an INPUT; runs on a given -t tree.
    if "genes:traits" in eset:
        if eset != {"genes:traits"}:
            parser.error("genes:traits runs on its own in this phase; combining it with other "
                         "edges is future work — see docs/coevolution_models.md")
        return _run_genes_traits(args, parser)

    # species:genes — speciation drives gene content (cladogenetic genome). An overlay edge (no
    # arrow into S), so the tree is an INPUT; runs on a given -t tree.
    if "species:genes" in eset:
        if eset != {"species:genes"}:
            parser.error("species:genes runs on its own in this phase; combining it with other "
                         "edges is future work — see docs/coevolution_models.md")
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
        t0 = time.perf_counter()
        res = simulate_traits(tree, model, cladogenesis=clado,
                              root_state=args.root_state, seed=args.seed)
        dt = time.perf_counter() - t0
        _write_coevolve_outputs(args.out, tree, res)
        return (f"wrote species:traits (cladogenetic {kind}) to {args.out}/ "
                f"({len(tree.extant_leaves())} tips{_sse_tip_signal(res)}) in {dt:.3g} s")

    # traits:species (SSE) or traits:species + species:traits (ClaSSE): an arrow INTO S, so the
    # tree is an OUTPUT — grow it forward with exactly one stopping condition (no input -t tree).
    if args.tree:
        parser.error("traits:species grows the tree (it is an OUTPUT); give --age/--tips, not an "
                     "input -t tree (that is the species:traits-alone edge)")
    if (args.age is None) == (args.tips is None):
        parser.error("traits:species grows the tree — give exactly one of --age or --tips")

    model = _build_sse_model(args, parser)
    t0 = time.perf_counter()
    res = simulate_sse(model, age=args.age, n_tips=args.tips, root_state=args.root_state,
                       cladogenesis=clado, seed=args.seed)
    dt = time.perf_counter() - t0
    _write_coevolve_outputs(args.out, res.tree, res)
    n_extant = len(res.tree.extant_leaves())
    edge_label = "traits:species+species:traits" if clado is not None else "traits:species"
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
    return (f"wrote genes:species (key innovations) to {args.out}/ "
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
    return (f"wrote species:genes (cladogenetic genome) to {args.out}/ "
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
    return (f"wrote genes:traits (gene-conditioned trait) to {args.out}/ "
            f"({len(tips)} tips; carrier trait mean {car_m} vs non-carrier {non_m}) in {dt:.3g} s")


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


def _run_genomes(tree: Tree, args: argparse.Namespace,
                 parser: argparse.ArgumentParser) -> str:
    """Simulate gene families along ``tree``, write output, and return a one-line summary.

    The default ``shared`` rate model runs on the Rust engine automatically (``simulate_genomes``
    raises a build hint if the extension is missing); ``per-genome`` runs on Python.
    """
    parts = set(Genomes.WRITE_PARTS) if "all" in args.output else set(args.output)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --write")
    if args.threads > 1 and parts != {"profiles"}:
        raise ValueError("--threads > 1 parallelises only the counts-only path; use it with "
                         "exactly --write profiles")

    if args.rate_model == "coupled":
        if args.genome_model != "unordered":
            parser.error("--rate-model coupled is an unordered-genome model and cannot combine "
                         f"with --genome-model {args.genome_model}")
        return _run_coupled(tree, args, parser, parts)

    if args.genome_model == "nucleotide":
        if args.initial_families is not None:
            parser.error("--initial-families is for the unordered genome level "
                         "(--genome-model unordered); the nucleotide model uses "
                         "--initial-chromosomes")
        return _run_nucleotides(tree, args, parts)

    if args.initial_chromosomes is not None:
        parser.error("--initial-chromosomes is only for --genome-model nucleotide; the "
                     "unordered and ordered genome levels use --initial-families")

    ordered = args.genome_model == "ordered"
    initial_families = 20 if args.initial_families is None else args.initial_families
    args.initial_families = initial_families  # record the effective value in the params log
    if args.rate_model == "per-genome":
        if ordered and (args.inversion is not None or args.transposition is not None):
            parser.error("rearrangements (--inversion/--transposition) need --rate-model shared; "
                         "per-genome rates do not carry them")
        model_kw = dict(rates=PerGenomeRates(args.dup, args.trans, args.loss, args.orig))
    elif ordered:  # shared per-copy rates + rearrangements on an ordered chromosome
        inv = 0.0 if args.inversion is None else args.inversion
        tps = 0.0 if args.transposition is None else args.transposition
        args.inversion, args.transposition = inv, tps  # record effective values in the params log
        model_kw = dict(rates=SharedRates(args.dup, args.trans, args.loss, args.orig,
                                          inversion=inv, transposition=tps))
    else:  # shared (unordered)
        model_kw = dict(duplication=args.dup, transfer=args.trans, loss=args.loss,
                        origination=args.orig)
    rate_kw = dict(**model_kw, initial_families=initial_families,
                   max_family_size=args.max_family_size, seed=args.seed)
    if ordered:
        ext = args.extension  # ordered event length is counted in genes; None -> single-gene events
        rate_kw["genome_factory"] = lambda ids, _e=ext: OrderedGenome(ids, extension=_e)

    # scoring reconciliation likelihoods needs the full gene-family genealogy, so it forces the
    # full path (the fast counts-only / trace paths don't reconstruct gene trees).
    score = getattr(args, "score_likelihoods", False)

    t0 = time.perf_counter()
    if parts == {"profiles"} and not score and not ordered:
        # counts-only Rust fast path: no genealogy reconstructed (parallel when --threads > 1)
        profiles = simulate_genomes(tree, output="profiles", threads=args.threads, **rate_kw)
        dt = time.perf_counter() - t0
        _write_profiles_only(args.out, tree, profiles, sparse=args.sparse)
        n_families = len(profiles.families)
    elif "trace" in parts and parts <= {"trace", "profiles"} and not score and not ordered:
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


def _run_coupled(tree: Tree, args: argparse.Namespace,
                 parser: argparse.ArgumentParser, parts: set) -> str:
    """Simulate a non-independent (Potts/Ising) gene-family panel along ``tree``.

    The coupling structure is a set of pathway blocks (``--pathways``): families within a block
    co-occur (``--within`` J), families across blocks couple by ``--between`` J. Coupling enters
    through loss (``loss_i = base-loss·exp(-beta·f_i)``); ``--gain-coupling`` additionally
    field-biases HGT establishment. The full genealogy is recorded, so every output component
    (profiles, gene trees, events, transfers, ...) is produced exactly as for the other rate
    models — gene-tree reconstruction is a model-agnostic function of the event log."""
    try:
        sizes = [int(s) for s in str(args.pathways).replace(" ", "").split(",") if s]
    except ValueError:
        parser.error("--pathways must be comma-separated integers, e.g. 4,4,3")
    if not sizes or any(s <= 0 for s in sizes):
        parser.error("--pathways must be positive integers, e.g. 4,4,3")

    spec = pathway_blocks(
        sizes, within=args.within, between=args.between, h=args.field,
        base_loss=args.base_loss, transfer=args.trans, origination=args.orig,
        beta=args.beta, gain_coupling=args.gain_coupling,
    )
    t0 = time.perf_counter()
    res = simulate_coupled(tree, spec, seed=args.seed)
    # Wrap the coupled result in the standard Genomes so it shares the full output machinery.
    # Keep res.profiles (the fixed panel: all N rows, incl. any globally-extinct family) rather
    # than letting Genomes derive a reduced one from the leaf genomes.
    genomes = Genomes(species_tree=tree, leaf_genomes=res.leaf_genomes,
                      event_log=res.event_log, profiles=res.profiles)
    genomes.write(args.out, include=parts, sparse=args.sparse,
                  annotate_species=args.annotate_species)
    dt = time.perf_counter() - t0
    return (f"wrote [{' '.join(sorted(parts))}] to {args.out}/ ({len(tree.leaves())} tips, "
            f"{spec.n_families} panel families, gain_coupling={args.gain_coupling:g}) in {dt:.3g} s")


def _write_reconciliation_likelihoods(genomes, args: argparse.Namespace) -> None:
    """Score every extant family's gene tree (ALElite) and write Reconciliation_likelihoods.tsv."""
    from .alelite import write_scores_tsv

    models = list(dict.fromkeys(args.score_model))  # de-dupe, keep order
    rows = genomes.reconciliation_likelihoods(
        args.dup, args.trans, args.loss, models=models,
        origination=args.score_origination, n_steps=args.score_nsteps,
    )
    write_scores_tsv(rows, os.path.join(args.out, "Reconciliation_likelihoods.tsv"), models=models)


def _run_nucleotides(tree: Tree, args: argparse.Namespace, parts: set) -> str:
    """Simulate nucleotide-resolution genomes (variable-length structural events) along ``tree``.

    Genes are not atomic here — they emerge as **blocks** (maximal intervals with one shared
    history). ``profiles`` writes the emergent block-by-species profile (plus ``blocks.tsv`` and
    the per-leaf ``Mosaics.tsv``); ``trees`` writes the per-block gene trees and their
    reconciliations. Only ``profiles``/``trees`` apply here (the family-model ``events`` /
    ``transfers`` / ``summary`` do not). ``profiles`` alone takes the fast Rust path.
    """
    want = parts & {"profiles", "trees", "ancestral"}
    if not want:
        raise ValueError("the nucleotide model writes 'profiles', 'trees' and/or 'ancestral'; "
                         "--write events/transfers/summary do not apply to it")
    ancestral = "ancestral" in want
    initial_chromosomes = 1 if args.initial_chromosomes is None else args.initial_chromosomes
    args.initial_chromosomes = initial_chromosomes  # record the effective value in the params log
    # the structural knobs are shared with the ordered level, so their defaults are resolved here
    args.inversion = 0.001 if args.inversion is None else args.inversion
    args.transposition = 0.0 if args.transposition is None else args.transposition
    args.extension = 0.99 if args.extension is None else args.extension
    if args.gff and args.genes:
        raise ValueError("give either --gff or --genes (not both) to set the gene coordinates")
    gff_info = None
    if args.gff:                              # start from a real genome: length + gene coordinates
        from .gff import read_gff
        gff_info = read_gff(args.gff, seqid=args.gff_seqid)
        genes = gff_info.genes
        args.root_length = gff_info.length    # GFF is authoritative for the chromosome length
    elif args.genes:
        genes = _read_gene_intervals(args.genes)
    else:
        genes = None
    genic = bool(genes)
    transfers = TransferModel(replacement=0.0) if genic else None  # homologous repl. is genome-side
    indels = bool(args.insertion or args.deletion)
    sim_kw = dict(inversion=args.inversion, loss=args.loss, duplication=args.dup,
                  transfer=args.trans, transposition=args.transposition,
                  origination=args.orig, insertion=args.insertion, deletion=args.deletion,
                  indel_mean_length=args.indel_mean_length, root_length=args.root_length,
                  extension=args.extension, initial_chromosomes=initial_chromosomes, seed=args.seed,
                  gene_intervals=genes, pseudogenization=args.pseudogenization,
                  replacement=args.replacement, transfers=transfers, retain_internal=ancestral)

    t0 = time.perf_counter()
    if "trees" in want or genic or ancestral or indels:  # these need the Python engine
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
    if "trees" in want:
        _write_block_gene_trees(args.out, result, genic=genic)
        result.write_reconciliations(args.out)   # Reconciled_complete/extant.nwk + events.tsv
        if genic:
            _write_pseudogenizations(args.out, result)
    if ancestral:
        _write_ancestral(args.out, result, tree, args, gff_info)

    if gff_info is not None:
        print(f"  GFF {gff_info.seqid}: {gff_info.length} bp, {gff_info.n_features} genes "
              f"-> {len(gff_info.genes)} after trimming ({gff_info.n_trimmed} trimmed, "
              f"{gff_info.n_dropped} dropped as overlapping)")
    extra = f", {len(result.gene_blocks())} genes" if genic else ""
    return (f"wrote [{' '.join(sorted(want))}] (nucleotide{'/genic' if genic else ''}) to "
            f"{args.out}/ ({len(result.leaf_genomes)} tips, {len(result.blocks)} blocks{extra}) "
            f"in {dt:.3g} s")


def _write_ancestral(out: str, result, tree, args, gff_info) -> None:
    """Simulate sequences and write the genome (architecture + gzipped DNA) at every node.

    ``Architecture/<node>.tsv`` — the ordered, oriented gene/intergene mosaic of the node's genome;
    ``Genomes/<node>.fasta.gz`` — its full assembled DNA (the root reproduces the input genome);
    ``Gene_alignments/<gene>.fasta`` — the extant per-gene alignments. The root sequence is seeded
    from ``--genome-fasta`` (the real genome) when given, else drawn at random.
    """
    from .sequence_sim import make_model, GammaRates, read_fasta, write_fasta
    model = make_model(args.subst_model, kappa=args.kappa,
                       freqs=args.base_freqs, rates=args.gtr_rates)
    gamma = GammaRates(args.gamma_shape) if args.gamma_shape else None
    root_fasta = None
    if args.genome_fasta:
        fa = read_fasta(args.genome_fasta)
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


def _add_abc_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="species tree (Newick) the empirical data evolved along")
    g.add_argument("--profiles", required=True, metavar="TSV",
                   help="empirical copy-number profile table (families x species TSV, like the "
                        "Profiles.tsv that 'zombi2 genomes' writes)")
    g.add_argument("--rate-model", dest="model", choices=("uniform", "family"), default="uniform",
                   metavar="MODEL",
                   help="uniform: one shared scalar rate per type (Rust; default); "
                        "family: per-family sampled rates, fitting each rate's mean (Python)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group(
        "priors",
        "the rates to fit — two values LOW HIGH (uniform) or one (fixed); omit to hold at 0. "
        "At least one must be a range")
    # priors reuse the genomes rate flags, but each takes a PRIOR (see _build_priors)
    for flag, param in (("--dup", "duplication"), ("--trans", "transfer"),
                        ("--loss", "loss"), ("--orig", "origination")):
        g.add_argument(flag, type=float, nargs="+", default=None, metavar="RATE",
                       help=f"{param} prior")

    g = p.add_argument_group("family model", "--rate-model family")
    g.add_argument("--family-shape", type=float, default=2.0, metavar="A",
                   help="Gamma shape for per-family rate dispersion (default 2.0)")

    g = p.add_argument_group("rejection ABC", "the default sampler")
    g.add_argument("--n-sims", type=int, default=1000, metavar="N",
                   help="number of prior simulations (default 1000)")
    g.add_argument("--accept", type=float, default=0.05, metavar="FRAC",
                   help="fraction of closest simulations to accept (default 0.05)")
    g.add_argument("--processes", type=int, default=None, metavar="N",
                   help="parallel worker processes (default: serial)")

    g = p.add_argument_group("ABC-SMC", "sequential sampler with shrinking tolerance")
    g.add_argument("--smc", action="store_true", help="use ABC-SMC instead of rejection")
    g.add_argument("--rounds", type=int, default=5, metavar="N",
                   help="number of SMC rounds (default 5)")
    g.add_argument("--particles", type=int, default=200, metavar="N",
                   help="particles per round (default 200)")
    g.add_argument("--quantile", type=float, default=0.5, metavar="Q",
                   help="tolerance quantile carried between rounds (default 0.5)")

    g = p.add_argument_group("posterior & simulation")
    g.add_argument("--regression-adjust", action="store_true",
                   help="also write the regression-adjusted posterior (Beaumont 2002)")
    g.add_argument("--initial-families", type=int, default=20, metavar="N",
                   dest="initial_families",
                   help="gene families seeded at the root of each simulation (default 20)")
    g.add_argument("--max-family-size", type=_int_or_float, default=None, metavar="CAP",
                   help="growth cap for each simulation — recommended with --rate-model family to "
                        "avoid runaway growth (integer = absolute, decimal = fraction of N)")


def _build_priors(args: argparse.Namespace) -> dict:
    """Turn the ``--dup/--trans/--loss/--orig`` flags into a priors dict for ``match_profiles``.

    Two values ``LOW HIGH`` -> a uniform prior on that rate; one value -> fixed; omitted ->
    the rate is held at 0. At least one rate must be given as a range (there must be something
    to fit).
    """
    priors: dict = {}
    for flag, param in (("dup", "duplication"), ("trans", "transfer"),
                        ("loss", "loss"), ("orig", "origination")):
        spec = getattr(args, flag)
        if spec is None:
            continue
        if len(spec) == 1:
            priors[param] = spec[0]                       # fixed value
        elif len(spec) == 2:
            priors[param] = (spec[0], spec[1])            # uniform (low, high)
        else:
            raise ValueError(f"--{flag} takes one value (fixed) or two (LOW HIGH), got {len(spec)}")
    if not any(isinstance(v, tuple) for v in priors.values()):
        raise ValueError("give at least one rate to fit as a range, e.g. --loss 0 1.5 (LOW HIGH)")
    return priors


def _write_abc_outputs(out: str, fit, adjusted: bool = False) -> None:
    """Write the ABC posterior, the per-parameter summary, and the spectrum diagnostic."""
    post = fit.posterior
    names = list(post)
    n_accept = len(next(iter(post.values())))
    lines = ["\t".join(names)]
    for i in range(n_accept):
        lines.append("\t".join(f"{post[nm][i]:.6g}" for nm in names))
    with open(os.path.join(out, "posterior.tsv"), "w") as f:      # accepted draws, one col/param
        f.write("\n".join(lines) + "\n")

    slines = ["parameter\tmean\tmedian\tlo95\thi95"]
    for nm, s in fit.summary().items():
        slines.append(f"{nm}\t{s['mean']:.6g}\t{s['median']:.6g}\t{s['lo95']:.6g}\t{s['hi95']:.6g}")
    if adjusted:
        slines.append("# regression-adjusted (Beaumont 2002)")
        for nm, s in fit.summary(adjusted=True).items():
            slines.append(f"{nm}_adj\t{s['mean']:.6g}\t{s['median']:.6g}\t"
                          f"{s['lo95']:.6g}\t{s['hi95']:.6g}")
    with open(os.path.join(out, "summary.tsv"), "w") as f:
        f.write("\n".join(slines) + "\n")

    if fit.uses_default_summary:                                  # posterior-predictive spectrum
        d = fit.spectra_data()
        lo, med, hi = np.percentile(d["accepted"], [2.5, 50, 97.5], axis=0)
        flines = ["k\tempirical\tacc_median\tacc_lo95\tacc_hi95"]
        for i, k in enumerate(d["k"]):
            flines.append(f"{int(k)}\t{d['empirical'][i]:.6g}\t{med[i]:.6g}\t"
                          f"{lo[i]:.6g}\t{hi[i]:.6g}")
        with open(os.path.join(out, "spectra.tsv"), "w") as f:
            f.write("\n".join(flines) + "\n")


def _run_abc(args: argparse.Namespace) -> str:
    """Fit gene-family rates to an empirical profile by ABC and write the posterior."""
    with open(args.tree) as f:
        tree = read_newick(f.read())
    empirical = ProfileMatrix.from_tsv(args.profiles)
    priors = _build_priors(args)
    common = dict(model=args.model, family_shape=args.family_shape,
                  initial_families=args.initial_families, max_family_size=args.max_family_size,
                  seed=args.seed)

    t0 = time.perf_counter()
    if args.smc:
        fit = match_profiles_smc(tree, empirical, priors, rounds=args.rounds,
                                 n_particles=args.particles, quantile=args.quantile, **common)
        effort = f"{args.rounds} SMC rounds x {args.particles} particles"
    else:
        fit = match_profiles(tree, empirical, priors, n_sims=args.n_sims, accept=args.accept,
                             processes=args.processes, **common)
        effort = f"{args.n_sims} sims"
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    _write_abc_outputs(args.out, fit, adjusted=args.regression_adjust)
    posterior = " ".join(f"{n}={s['median']:.3g}[{s['lo95']:.3g},{s['hi95']:.3g}]"
                         for n, s in fit.summary().items())
    return (f"fit {len(fit.accepted)} accepted / {effort}, tol={fit.tolerance:.3g} in {dt:.3g} s "
            f"-> {args.out}/ (median [95% CI]: {posterior})")


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
                             "shared across families; pick lognormal OR discrete-bin")
    g.add_argument("--branch-speed", type=float, default=0.0, metavar="SIGMA",
                   help="lognormal relaxed clock: autocorrelated, drift SIGMA per sqrt(time) "
                        "(0 = strict clock). Exclusive with --branch-bins")
    g.add_argument("--branch-bins", default=None, metavar="R1,R2,...",
                   help="discrete-bin within-branch GTDB clock: comma-separated ORDERED rate "
                        "multipliers (e.g. 0.25,0.5,1,2,4), a Markov walk between adjacent bins")
    g.add_argument("--branch-switch-rate", type=float, default=1.0, metavar="RATE",
                   help="[--branch-bins] rate of stepping to a neighbouring bin (default 1.0)")
    g.add_argument("--branch-up-bias", type=float, default=0.5, metavar="P",
                   help="[--branch-bins] probability a step goes to the faster neighbour "
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
    from .profiles import _natkey
    from .reconciliation import extant_species_from_records
    from .sequence_sim import (GammaRates, evolve_on_tree, is_protein_model, make_model,
                               read_fasta, write_fasta)
    from .simulation import read_events_trace

    if args.family_speed < 0 or args.branch_speed < 0:
        raise ValueError("--family-speed / --branch-speed must be >= 0")
    if args.branch_speed > 0 and args.branch_bins:
        raise ValueError("--branch-speed (lognormal clock) and --branch-bins (discrete-bin "
                         "clock) are two lineage clocks; give at most one")
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
    if args.branch_bins:
        bins = [float(x) for x in args.branch_bins.split(",") if x.strip() != ""]
        se = SequenceEvolution(
            lineage=RateVariation(bins=bins, switch_rate=args.branch_switch_rate,
                                  up_bias=args.branch_up_bias),
            family_speed=family_speed)
    else:
        se = SequenceEvolution(branch_sigma=args.branch_speed, family_speed=family_speed)

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

    clock = f"branch-bins [{args.branch_bins}]" if args.branch_bins else f"branch-speed {args.branch_speed}"
    msg = (f"wrote substitution-unit gene trees for {n} families to {args.out}/gene_trees/ "
           f"(clock: {clock}, family-speed {args.family_speed}) in {dt:.3g} s")

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


def _add_subcommand(sub, name: str, help: str, description: str, usage: str, adder):
    """Register a subcommand with the house-style formatter and a hand-written compact usage.

    The command list itself is curated (grouped by theme) in the top-level description, so the
    per-command ``help`` is suppressed from argparse's auto listing to avoid a duplicate dump.
    """
    p = sub.add_parser(name, help=help, description=description, usage=usage,
                       formatter_class=ZombiHelpFormatter)
    adder(p)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="zombi2", description=_banner() + "\n\n" + _DESCRIPTION,
        formatter_class=ZombiHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"ZOMBI2 {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)

    _add_subcommand(
        sub, "species", "simulate a dated species tree",
        "Simulate a dated species tree.",
        "zombi2 species -o DIR [--mode MODE] [--diversification PROCESS] [options]",
        _add_species_args)

    _add_subcommand(
        sub, "genomes", "evolve gene families along a species tree",
        "Evolve gene families along a species tree.",
        "zombi2 genomes -t FILE -o DIR [--genome-model LEVEL] [--rate-model MODEL] "
        "[--write PART ...] [options]",
        _add_rate_args)

    _add_subcommand(
        sub, "trait", "evolve a phenotypic trait along a given species tree",
        "Evolve a phenotypic trait along a species tree, writing tip and ancestral values.",
        "zombi2 trait -t FILE -o DIR [--model MODEL] [options]", _add_trait_args)

    _add_subcommand(
        sub, "abc", "fit gene-family rates to an empirical profile by ABC",
        "Fit gene-family rates to an empirical copy-number profile by Approximate Bayesian "
        "Computation (the inverse of 'genomes').",
        "zombi2 abc -t FILE --profiles TSV -o DIR [--dup LOW HIGH ...] [options]", _add_abc_args)

    _add_subcommand(
        sub, "coevolve", "co-evolve coupled processes (--couple driver:target)",
        "Co-evolve coupled processes over {species, traits, genes} — pick directed edges with "
        "--couple (e.g. traits:species = SSE, traits:genes = trait-conditioned gene families).",
        "zombi2 coevolve -o DIR --couple DRIVER:TARGET [-t FILE] [--age T|--tips N] [options]",
        _add_coevolve_mode_args)

    _add_subcommand(
        sub, "sequence", "simulate DNA/protein alignments along a genomes run's gene trees",
        "Rescale a 'genomes' run's gene trees from time into substitutions/site under a "
        "gene × lineage clock, then (with --subst-model) simulate a DNA or protein sequence "
        "alignment along each rescaled gene tree.",
        "zombi2 sequence --genomes DIR -o DIR [--subst-model MODEL] "
        "[--branch-speed SIGMA|--branch-bins ...] [options]",
        _add_sequence_args)

    args = parser.parse_args(argv)
    print(_banner(), file=sys.stderr)          # a banner on each run (stderr keeps stdout clean)
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
            f.write(tree.to_newick() + "\n")
        leaves = tree.leaves()
        n_extant = sum(1 for n in leaves if n.is_extant)
        dead = len(leaves) - n_extant
        extra = f" + {dead} extinct" if dead else ""
        summary = f"{n_extant} extant{extra} tips"
        print(f"wrote {args.out}/species_tree.nwk ({summary}) in {dt:.3g} s")
        _write_params_log(os.path.join(args.out, "species_tree.log"), args, summary)
        return 0

    if args.command == "genomes":
        with open(args.tree) as f:
            tree = read_newick(f.read())
        summary = _run_genomes(tree, args, parser)
        print(summary)
        _write_params_log(os.path.join(args.out, "genomes.log"), args, summary)
        return 0

    if args.command == "trait":
        summary = _run_trait(args)
        print(summary)
        _write_params_log(os.path.join(args.out, "trait.log"), args, summary)
        return 0

    if args.command == "abc":
        summary = _run_abc(args)
        print(summary)
        _write_params_log(os.path.join(args.out, "abc.log"), args, summary)
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

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
