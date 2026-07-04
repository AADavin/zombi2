"""Command-line interface for ZOMBI2 (``zombi2 species`` / ``genomes`` / ``trait``)."""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

from .biogeography import DEC, simulate_biogeography
from .ghosts import add_ghost_lineages
from .matching import match_profiles, match_profiles_smc
from .nucleotide_sim import simulate_nucleotide_genomes
from .profiles import ProfileMatrix
from .distributions import LogNormal
from .rate_variation import RateVariation
from .rates import GenomeWiseRates
from .sequence_evolution import SequenceEvolution
from .simulation import Genomes, simulate_genomes
from .species_model import (
    BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent, EpisodicBirthDeath,
)
from .species_sim import simulate_species_tree
from .sse import BiSSE, MuSSE, QuaSSE, simulate_sse
from .gene_diversification import GeneDiversification, simulate_gene_diversification
from .trait_coupling import TraitGeneCoupling, simulate_trait_linked_genomes
from .traits import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst, Mk, ThresholdModel, TraitResult,
    Cladogenesis, simulate_traits,
)
from .transfers import TransferModel
from .tree import Tree, read_newick

_DESCRIPTION = """\
ZOMBI2 — a phylogenetic simulator of species trees and gene families.

Simulate in two steps: build a species tree, then evolve gene families along it.

  zombi2 species   simulate a species tree
  zombi2 genomes   evolve gene families along a species tree (Newick)
  zombi2 trait     evolve a phenotypic trait along a given species tree
  zombi2 abc       fit gene-family rates to an empirical profile (ABC inference)
  zombi2 sequence  rescale a genomes run's gene trees into substitutions/site
  zombi2 coevolve-genetrait  evolve gene families conditioned on a trait
  zombi2 coevolve  co-evolve coupled processes (--couple traits:species = SSE)

Run 'zombi2 <command> -h' for a command's options.
"""


def _int_or_float(text: str) -> int | float:
    """Parse ``--max-family-size``: a plain integer is an absolute cap, a value with a
    decimal point is a fraction of the number of species (e.g. ``0.5`` -> half of N)."""
    return float(text) if "." in text else int(text)


def _add_species_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--model", choices=("backward", "forward"), default="backward",
                   help="backward: reconstructed tree conditioned on --tips extant species "
                        "(default); forward: complete tree grown in time, keeping extinct "
                        "lineages (and fossils)")
    p.add_argument("--diversification", choices=("constant", "clads", "diversity-dependent"),
                   default="constant",
                   help="diversification process (forward only): constant-rate birth–death "
                        "(default); clads = per-lineage rates that shift at each speciation "
                        "(ClaDS); diversity-dependent = rates decline toward a carrying capacity")
    p.add_argument("--birth", type=float, nargs="+", default=[1.0], metavar="RATE",
                   help="speciation rate (default 1.0); several values with --shifts give an "
                        "episodic (skyline) model. For clads/diversity-dependent it is the "
                        "root/intrinsic rate λ₀ (a single value)")
    p.add_argument("--death", type=float, nargs="+", default=[0.3], metavar="RATE",
                   help="extinction rate (default 0.3); several values with --shifts give an "
                        "episodic (skyline) model. The constant μ for --diversification "
                        "diversity-dependent (clads uses --turnover instead)")
    p.add_argument("--shifts", type=float, nargs="+", default=None, metavar="AGE",
                   help="episodic rate-shift ages, present -> past (K-1 ages for K rate values)")
    p.add_argument("--tips", type=int, default=None,
                   help="number of extant species (backward default 50; forward: --tips OR --age)")
    p.add_argument("--age", type=float, default=None,
                   help="tree age / timescale, in the same time units as the rates "
                        "(backward default 1.0; forward: --tips OR --age)")
    p.add_argument("--age-type", choices=("crown", "stem"), default="crown",
                   help="interpret --age as crown (default) or stem age [backward]")
    p.add_argument("--sampling-fraction", type=float, default=1.0, metavar="RHO",
                   help="[forward] fraction of extant species sampled, 0<rho<=1 (default 1.0)")
    p.add_argument("--fossilization", type=float, default=0.0, metavar="PSI",
                   help="[forward] fossil (serial) sampling rate psi — fossilized birth–death "
                        "(default 0 = no fossils)")
    p.add_argument("--removal", type=float, default=1.0, metavar="R",
                   help="[forward] removal probability on sampling, 0<=r<=1 (r<1 keeps sampled "
                        "ancestors; default 1.0)")
    p.add_argument("--mass-extinction", action="append", nargs=2, type=float,
                   metavar=("AGE", "FRACTION"), default=None, dest="mass_extinction",
                   help="[forward] a mass extinction: at AGE before the present, each lineage "
                        "dies with probability FRACTION (0<FRACTION<=1). Repeat for several "
                        "pulses, e.g. --mass-extinction 1.0 0.75 --mass-extinction 2.5 0.5")
    p.add_argument("--clads-alpha", type=float, default=0.9, metavar="ALPHA",
                   help="[--diversification clads] speciation-rate trend per branch; α<1 = "
                        "rates slow toward the present (default 0.9)")
    p.add_argument("--clads-sigma", type=float, default=0.1, metavar="SIGMA",
                   help="[--diversification clads] lognormal spread of the per-branch rate "
                        "shift (default 0.1)")
    p.add_argument("--turnover", type=float, default=0.0, metavar="EPS",
                   help="[--diversification clads] extinction/speciation ratio ε=μ/λ, in [0,1) "
                        "(0 = pure birth; default 0.0)")
    p.add_argument("--carrying-capacity", "-K", type=float, default=None, metavar="K",
                   help="[--diversification diversity-dependent] carrying capacity K; the "
                        "speciation rate is λ₀·(1−n/K) (required for this model)")
    p.add_argument("--clade-shift", action="append", nargs=3, type=float,
                   metavar=("AGE", "BIRTH", "DEATH"), default=None, dest="clade_shift",
                   help="[forward] a clade-specific rate shift: at AGE before the present, one "
                        "random lineage then alive (and its descendants) switches to speciation "
                        "BIRTH / extinction DEATH. Repeat for several shifting clades, e.g. "
                        "--clade-shift 3.0 2.5 0.1")
    p.add_argument("--ghosts", action="store_true",
                   help="[backward] graft the extinct/unsampled 'ghost' lineages back onto the tree")
    p.add_argument("--ghost-method", choices=("rejection", "htransform"), default="rejection",
                   help="ghost-subtree sampler used with --ghosts (default rejection)")
    p.add_argument("--max-attempts", type=int, default=10000,
                   help="[forward] retries before giving up when the process goes extinct "
                        "(default 10000)")
    p.add_argument("--max-lineages", type=int, default=1_000_000,
                   help="[forward] abort a run exceeding this many live lineages (default 1000000)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


def _add_rate_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--rate-model", choices=("uniform", "genome-wise", "nucleotide"),
                   default="uniform",
                   help="uniform: same per-copy rates for every family (Rust; default); "
                        "genome-wise: constant per-genome rates, linear growth (Python); "
                        "nucleotide: nucleotide-resolution genomes evolving by variable-length "
                        "structural events, genes emerge as 'atoms' (see the nucleotide options)")
    p.add_argument("--dup", type=float, default=0.0,
                   help="duplication rate (per copy; per nucleotide when --rate-model nucleotide)")
    p.add_argument("--trans", type=float, default=0.0,
                   help="transfer rate (per copy; per nucleotide when --rate-model nucleotide)")
    p.add_argument("--loss", type=float, default=0.0,
                   help="loss/deletion rate (per copy; per nucleotide when --rate-model nucleotide)")
    p.add_argument("--orig", type=float, default=0.0, help="origination rate (per branch)")
    p.add_argument("--initial-size", type=int, default=None,
                   help="genomes seeded at the root (default: 20 gene families; "
                        "1 root chromosome for --rate-model nucleotide)")
    p.add_argument("--max-family-size", type=_int_or_float, default=None,
                   help="bound family growth: integer = absolute cap, "
                        "decimal = fraction of the number of species (e.g. 0.5) "
                        "[not used by --rate-model nucleotide]")
    p.add_argument("--output", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "ancestral", "all"), default=["profiles", "trees"],
                   help="which output files to write — any of {profiles, trace, trees, events, "
                        "transfers, summary} or 'all' (default: profiles trees). "
                        "species_tree.nwk is always written; 'profiles' alone takes the fast "
                        "Rust counts-only path; 'trace' (optionally with 'profiles') writes the "
                        "compact single-file event log Events_trace.tsv near counts-only speed, "
                        "from which gene trees can be reconstructed later on demand. "
                        "[nucleotide] 'ancestral' simulates DNA sequences and reconstructs the "
                        "genome (architecture + gzipped FASTA) at every node (root = input genome)")
    p.add_argument("--sparse", action="store_true",
                   help="write the profile as a sparse long table (Profiles_sparse.tsv: "
                        "family/species/copies, present cells only) instead of the dense "
                        "matrix — the scalable output for huge trees (needs 'profiles' in --output)")
    p.add_argument("--annotate-species", action="store_true",
                   help="label internal gene-tree nodes <gid>|<species-branch> (e.g. g570|i5)")
    # --- nucleotide model only (--rate-model nucleotide) ---
    p.add_argument("--inversion", type=float, default=0.001,
                   help="[nucleotide] per-nucleotide inversion rate (default 0.001)")
    p.add_argument("--transposition", type=float, default=0.0,
                   help="[nucleotide] per-nucleotide transposition rate (default 0)")
    p.add_argument("--root-length", type=int, default=1000,
                   help="[nucleotide] length of the root chromosome, in nucleotides (default 1000)")
    p.add_argument("--extension", type=float, default=0.99,
                   help="[nucleotide] geometric event-length parameter; mean event length is "
                        "1/(1-extension) nucleotides (default 0.99)")
    # --- genes & intergenes (nucleotide model) ---
    p.add_argument("--gff", metavar="FILE", default=None,
                   help="[nucleotide] a GFF3 annotation (optionally .gz) to start from a real "
                        "genome: sets the chromosome length and the gene coordinates (intergenes "
                        "are the gaps). Overlapping genes are trimmed to be disjoint. Enables "
                        "genic mode; supersedes --genes/--root-length")
    p.add_argument("--gff-seqid", metavar="ID", default=None,
                   help="[nucleotide] which GFF sequence to read (default: the most-annotated "
                        "one — the chromosome of a single-chromosome bacterium)")
    p.add_argument("--genes", metavar="FILE", default=None,
                   help="[nucleotide] BED/TSV of gene intervals on the root chromosome (columns: "
                        "start end [name], 0-based half-open) — an alternative to --gff. Enables "
                        "genic mode: event breakpoints fall only in intergene positions so genes "
                        "are never split; genes and intergenes are recovered as separate tree sets")
    p.add_argument("--pseudogenization", type=float, default=0.0,
                   help="[nucleotide, genic] probability a loss hitting a gene demotes it to "
                        "intergene (sequence retained, a state change) instead of deleting it "
                        "(default 0)")
    p.add_argument("--replacement", type=float, default=0.0,
                   help="[nucleotide, genic] probability a transfer is a homologous replacement "
                        "(the copy replaces the recipient's syntenic locus, located via flanking "
                        "genes; additive when no homolog) (default 0)")
    # --- sequences + ancestral genomes (--output ancestral) ---
    p.add_argument("--subst-model", choices=("jc69", "k80", "hky85", "gtr"), default="hky85",
                   help="[nucleotide, --output ancestral] nucleotide substitution model for the "
                        "sequences (default hky85)")
    p.add_argument("--kappa", type=float, default=2.0,
                   help="[nucleotide] transition/transversion ratio for k80/hky85 (default 2.0)")
    p.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[nucleotide] equilibrium base frequencies for hky85/gtr (default equal)")
    p.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[nucleotide] the 6 GTR exchangeabilities (default all 1)")
    p.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="[nucleotide] discrete-Gamma across-site rate heterogeneity shape "
                        "(default: none / uniform rates)")
    p.add_argument("--subst-rate", type=float, default=1.0,
                   help="[nucleotide] overall substitutions/site per unit time — scales sequence "
                        "divergence (default 1.0)")
    p.add_argument("--genome-fasta", metavar="FILE", default=None,
                   help="[nucleotide, --output ancestral] the input genome's DNA (FASTA, optionally "
                        ".gz) to seed the root sequence; without it the root is drawn at random")


def _build_species_model(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Construct a species-tree model (BirthDeath / EpisodicBirthDeath / ClaDS /
    DiversityDependent) from the CLI args (validated)."""
    if args.model == "backward" and (args.fossilization or args.removal != 1.0
                                     or args.sampling_fraction != 1.0):
        parser.error("--fossilization / --removal / --sampling-fraction require --model forward "
                     "(the backward reconstructed sampler assumes complete sampling)")
    if args.model == "backward" and args.mass_extinction:
        parser.error("--mass-extinction requires --model forward (mass extinctions kill real "
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
                     "add --model forward")
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
        parser.error("--clade-shift requires --model forward (the shifts play out forward in time)")
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
    p.add_argument("-t", "--tree", required=True,
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    p.add_argument("--model", choices=("bm", "ou", "eb", "mk", "threshold", "dec"), default="bm",
                   help="trait model: bm=Brownian motion, ou=Ornstein-Uhlenbeck, "
                        "eb=early burst/ACDC, mk=discrete k-state, threshold, "
                        "dec=geographic-range evolution (default: bm)")
    p.add_argument("--sigma2", type=float, default=1.0,
                   help="diffusion rate [bm/ou/eb/threshold] (default: 1.0)")
    p.add_argument("--x0", type=float, default=None,
                   help="root value [bm/eb/threshold]; OU defaults it to --theta")
    p.add_argument("--trend", type=float, default=0.0, help="directional drift [bm/eb]")
    p.add_argument("--alpha", type=float, default=1.0,
                   help="OU mean-reversion strength [ou] (default: 1.0)")
    p.add_argument("--theta", type=float, default=0.0, help="OU optimum [ou] (default: 0.0)")
    p.add_argument("--rate", type=float, default=1.0,
                   help="EB rate-of-change (negative = early burst) [eb], "
                        "or the per-transition rate [mk] (default: 1.0)")
    p.add_argument("--states", type=int, default=2,
                   help="number of states for the mk model (default: 2)")
    p.add_argument("--ordered", action="store_true",
                   help="[mk] only allow transitions between adjacent states (i <-> i±1)")
    p.add_argument("--q-matrix", default=None,
                   help="[mk] path to a whitespace/comma-separated k x k rate matrix (an "
                        "arbitrary Markov chain); overrides --states/--rate/--ordered")
    p.add_argument("--thresholds", default="0.0",
                   help="comma-separated liability cut points [threshold] (default: 0.0)")
    # DEC (geographic-range evolution)
    p.add_argument("--areas", default="3",
                   help="[dec] number of areas (e.g. 3) or comma-separated area labels "
                        "(e.g. A,B,C) (default: 3)")
    p.add_argument("--dispersal", type=float, default=0.1,
                   help="[dec] rate of gaining an area (dispersal) (default: 0.1)")
    p.add_argument("--extinction", type=float, default=0.1,
                   help="[dec] rate of losing an area (local extinction) (default: 0.1)")
    p.add_argument("--max-range-size", type=int, default=None,
                   help="[dec] maximum number of areas a range may span (default: all)")
    p.add_argument("--root-range", default=None,
                   help="[dec] comma-separated area labels for the root range (e.g. A); "
                        "default: a random range")
    p.add_argument("--replicates", type=int, default=1,
                   help="simulate the trait this many times with the same parameters; writes "
                        "traits.tsv with one column per replicate (default: 1)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


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
# coevolve-genetrait: trait-conditioned gene-family dynamics
# ═══════════════════════════════════════════════════════════════════════════════
def _add_trait_model_args(p: argparse.ArgumentParser) -> None:
    """Trait-model flags for the coevolve command (scalar traits; DEC ranges do not apply).

    ``--trait-model`` stores into ``args.model`` so :func:`_build_trait_model` is reused as-is.
    """
    p.add_argument("--trait-model", dest="model", default="bm",
                   choices=("bm", "ou", "eb", "mk", "threshold"),
                   help="trait to evolve then couple to gene families: bm=Brownian motion, "
                        "ou=Ornstein-Uhlenbeck, eb=early burst, mk=discrete k-state, threshold "
                        "(default: bm). Use --trait-file to supply a precomputed trait instead")
    p.add_argument("--sigma2", type=float, default=1.0,
                   help="diffusion rate [bm/ou/eb/threshold] (default: 1.0)")
    p.add_argument("--x0", type=float, default=None,
                   help="root value [bm/eb/threshold]; OU defaults it to --theta")
    p.add_argument("--trend", type=float, default=0.0, help="directional drift [bm/eb]")
    p.add_argument("--alpha", type=float, default=1.0, help="OU mean-reversion strength [ou]")
    p.add_argument("--theta", type=float, default=0.0, help="OU optimum [ou]")
    p.add_argument("--rate", type=float, default=1.0,
                   help="EB rate-of-change (negative = early burst) [eb], or per-transition rate [mk]")
    p.add_argument("--states", type=int, default=2, help="number of states for the mk model (default: 2)")
    p.add_argument("--ordered", action="store_true",
                   help="[mk] only allow transitions between adjacent states (i <-> i±1)")
    p.add_argument("--q-matrix", default=None,
                   help="[mk] path to a k x k rate matrix; overrides --states/--rate/--ordered")
    p.add_argument("--thresholds", default="0.0",
                   help="comma-separated liability cut points [threshold] (default: 0.0)")


def _add_coevolve_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("-t", "--tree", required=True,
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    _add_trait_model_args(p)
    p.add_argument("--trait-file", default=None, metavar="TSV",
                   help="use a precomputed trait instead of simulating one: a node<TAB>value table "
                        "over ALL nodes (tips and ancestors), as 'zombi2 trait' writes with "
                        "nodes=all; values must be numeric (encode discrete states as numbers). "
                        "Overrides --trait-model")
    p.add_argument("--trait-center", action="store_true",
                   help="[discrete trait] center the state values around their mean so the trait "
                        "pushes retention both up and down — recommended for a binary "
                        "aerobic/anaerobic trait; by default states are 0,1,..,k-1")
    p.add_argument("--trait-steps", type=int, default=16, metavar="K",
                   help="[continuous trait] within-branch resolution: sub-segment each branch into "
                        "K pieces (default 16; ignored for discrete traits, which use their exact "
                        "stochastic map)")
    # --- gene-family panel and its base rates ---
    p.add_argument("--panel", type=int, default=50,
                   help="number of gene families in the panel (default 50)")
    p.add_argument("--loss", type=float, default=0.5,
                   help="baseline per-copy loss rate — the loss where the trait is neutral (default 0.5)")
    p.add_argument("--trans", type=float, default=1.0,
                   help="per-copy transfer (HGT) rate — the field-blind gain channel (default 1.0)")
    p.add_argument("--dup", type=float, default=0.0,
                   help="per-copy duplication rate, trait-independent (default 0)")
    p.add_argument("--orig", type=float, default=0.0,
                   help="background origination rate of brand-new, uncoupled families (default 0)")
    # --- the trait <-> gene-family coupling ---
    p.add_argument("--responsive", default="0.3", metavar="SPEC",
                   help="which families respond to the trait: an integer count, a fraction "
                        "(e.g. 0.3), a comma-separated id/index list (e.g. F3,F7,12), or @FILE of "
                        "ids/indices (default: 0.3 = 30%% of the panel, chosen at random)")
    p.add_argument("--weight", type=float, default=1.0,
                   help="coupling weight of each responsive family (default 1.0)")
    p.add_argument("--signed", action="store_true",
                   help="randomise the sign of each responsive weight (some families favoured by a "
                        "high trait value, some by a low one); by default all favour a high value")
    p.add_argument("--effect-loss", type=float, default=2.0,
                   help="retention coupling strength: a responsive family's loss scales by "
                        "exp(-effect_loss * weight * trait) (default 2.0; 0 = no coupling)")
    p.add_argument("--effect-gain", type=float, default=0.0,
                   help="optional HGT-activity coupling: a lineage's transfer rate scales by "
                        "exp(effect_gain * trait) (default 0 = field-blind gain, as in the Potts model)")
    p.add_argument("--output", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "all"), default=["profiles", "trees"],
                   help="which gene-family outputs to write — any of {profiles, trace, trees, "
                        "events, transfers, summary} or 'all' (default: profiles trees). "
                        "traits.tsv, trait_tree.nwk and coupling.tsv (the responsive-family "
                        "manifest) are always written alongside")
    p.add_argument("--sparse", action="store_true",
                   help="write the profile as a sparse long table (needs 'profiles' in --output)")
    p.add_argument("--annotate-species", action="store_true",
                   help="label internal gene-tree nodes <gid>|<species-branch> (e.g. g570|i5)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


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


def _run_coevolve(args: argparse.Namespace) -> str:
    """Simulate a trait, then evolve a gene-family panel whose loss/gain is conditioned on it."""
    with open(args.tree) as f:
        tree = read_newick(f.read())
    parts = set(Genomes.WRITE_PARTS) if "all" in args.output else set(args.output)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --output")
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
    p.add_argument("--couple", action="append", nargs="+", metavar="DRIVER:TARGET", default=None,
                   help="a directed coupling edge 'driver:target' over {species, traits, genes} — "
                        "the driver's state modulates the target's rates. Phase 1 implements "
                        "'traits:species' (SSE: a trait drives speciation/extinction), "
                        "'species:traits' (cladogenetic: speciation jumps the trait), and their "
                        "combination = ClaSSE. Repeatable; default traits:species. See "
                        "docs/coevolution_models.md for the full edge set")
    p.add_argument("-t", "--tree", default=None,
                   help="input species tree (Newick) — required for 'species:traits' ALONE (the "
                        "trait evolves along a GIVEN tree). Omit for the into-species edges "
                        "(traits:species / ClaSSE), which GROW the tree via --age/--tips")
    # into-species edges grow the tree, so they take a stopping condition (not an input -t tree)
    p.add_argument("--age", type=float, default=None,
                   help="[traits:species] crown age to grow for (the extant tip count is random)")
    p.add_argument("--tips", type=int, default=None,
                   help="[traits:species] stop when this many extant tips first coexist (age random)")
    p.add_argument("--sse-model", dest="sse_model", choices=("bisse", "musse", "quasse"),
                   default="bisse",
                   help="[traits:species] which state-dependent model drives diversification: "
                        "bisse (binary trait), musse (k-state), quasse (continuous trait) "
                        "(default: bisse)")
    # BiSSE — per-state speciation/extinction and asymmetric transitions
    p.add_argument("--lambda0", type=float, default=1.0, help="[bisse] speciation in state 0")
    p.add_argument("--lambda1", type=float, default=2.0, help="[bisse] speciation in state 1")
    p.add_argument("--mu0", type=float, default=0.3, help="[bisse] extinction in state 0")
    p.add_argument("--mu1", type=float, default=0.3, help="[bisse] extinction in state 1")
    p.add_argument("--q01", type=float, default=0.1, help="[bisse] transition rate 0 -> 1")
    p.add_argument("--q10", type=float, default=0.1, help="[bisse] transition rate 1 -> 0")
    # MuSSE — general k-state
    p.add_argument("--birth", type=float, nargs="+", default=None, metavar="RATE",
                   help="[musse] per-state speciation rates (k values)")
    p.add_argument("--death", type=float, nargs="+", default=None, metavar="RATE",
                   help="[musse] per-state extinction rates (k values)")
    p.add_argument("--q-matrix", default=None, metavar="FILE",
                   help="[musse] path to a k x k anagenetic transition-rate matrix (same format as "
                        "'zombi2 trait --q-matrix')")
    p.add_argument("--root-state", type=int, default=None,
                   help="[bisse/musse] root state index (default: the character's stationary "
                        "distribution)")
    # QuaSSE — continuous trait, sigmoidal speciation + constant extinction
    p.add_argument("--spec-low", type=float, default=0.5,
                   help="[quasse] speciation rate at low trait values")
    p.add_argument("--spec-high", type=float, default=2.0,
                   help="[quasse] speciation rate at high trait values")
    p.add_argument("--spec-center", type=float, default=0.0,
                   help="[quasse] trait value at the middle of the speciation sigmoid")
    p.add_argument("--spec-slope", type=float, default=1.0,
                   help="[quasse] steepness of the speciation sigmoid")
    p.add_argument("--qmu", type=float, default=0.1, help="[quasse] constant extinction rate")
    p.add_argument("--diffusion", type=float, default=1.0,
                   help="[quasse] trait diffusion rate sigma^2 (Brownian motion)")
    p.add_argument("--root-value", type=float, default=0.0, help="[quasse] root trait value x0")
    # species:traits — the cladogenetic (speciation -> trait) kernel; used when species:traits is
    # in --couple (on its own, or combined with traits:species for ClaSSE)
    p.add_argument("--clado-shift", dest="clado_shift", type=float, default=0.3,
                   help="[species:traits, discrete trait] probability a daughter hops to another "
                        "state AT each speciation (cladogenetic change; default 0.3)")
    p.add_argument("--clado-jump", dest="clado_jump", type=float, default=1.0,
                   help="[species:traits, continuous trait] variance of the Gaussian jump added to "
                        "each daughter's value AT each speciation (default 1.0)")
    # genes:species — gene-content-dependent diversification (key innovations + HGT). The base
    # (no-driver) speciation/extinction rates reuse --lambda0/--mu0.
    p.add_argument("--drivers", type=int, default=2,
                   help="[genes:species] number of binary 'driver' (key-innovation) gene families")
    p.add_argument("--driver-speciation", dest="driver_speciation", type=float, default=1.0,
                   help="[genes:species] per-driver effect on log speciation: a present driver "
                        "scales lambda by exp(this) (>0 = a key innovation; default 1.0). Base "
                        "lambda0 = --lambda0")
    p.add_argument("--driver-extinction", dest="driver_extinction", type=float, default=0.0,
                   help="[genes:species] per-driver effect on log extinction: a present driver "
                        "scales mu by exp(this) (default 0). Base mu0 = --mu0")
    p.add_argument("--driver-loss", dest="driver_loss", type=float, default=0.1,
                   help="[genes:species] rate a present driver is lost/deleted (default 0.1)")
    p.add_argument("--driver-origination", dest="driver_origination", type=float, default=0.05,
                   help="[genes:species] rate an absent driver appears de novo (default 0.05)")
    p.add_argument("--driver-transfer", dest="driver_transfer", type=float, default=0.5,
                   help="[genes:species] per-donor HGT rate of a driver — frequency-dependent gain: "
                        "a driver in more live genomes spreads faster (default 0.5)")
    p.add_argument("--root-drivers", dest="root_drivers", type=int, default=0,
                   help="[genes:species] number of drivers present at the root (the first m; "
                        "default 0 = drivers enter by origination)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


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
    """Run the ``coevolve`` umbrella. Implemented edges: ``traits:species`` (SSE), ``species:traits``
    (cladogenetic — the trait jumps at speciation), and their combination = **ClaSSE**. Whether the
    tree is grown (an arrow into species) or read from ``-t`` follows the arrows-into-S rule."""
    # --couple accepts both repeated flags and space-separated lists (append + nargs); flatten
    raw = args.couple or [["traits:species"]]
    edges = [e.strip().lower() for group in raw for e in group]
    for e in edges:
        if e not in _COEVOLVE_EDGES:
            parser.error(f"unknown --couple edge {e!r}: expected 'driver:target' over "
                         f"{{{', '.join(_COEVOLVE_NODES)}}} (e.g. traits:species); see "
                         "docs/coevolution_models.md for the full edge set")
    eset = set(edges)
    supported = {"traits:species", "species:traits", "genes:species"}
    unsupported = eset - supported
    if unsupported:
        if eset == {"traits:genes"}:
            parser.error("the traits:genes edge ships today as 'zombi2 coevolve-genetrait' — use "
                         "that command (it will be folded in as 'coevolve --couple traits:genes')")
        parser.error(f"--couple {', '.join(sorted(unsupported))} is planned but not yet "
                     "implemented; the built edges are traits:species (SSE), species:traits "
                     "(cladogenetic), their combination (ClaSSE), and genes:species (key "
                     "innovations). See docs/coevolution_models.md")

    # genes:species — gene content drives diversification (its own forward joint loop, v1 stands
    # alone; combining it with other edges is the full joint model, still on the roadmap).
    if "genes:species" in eset:
        if eset != {"genes:species"}:
            parser.error("genes:species runs on its own in this phase; combining it with other "
                         "edges (the fully joint model) is future work — see docs/coevolution_models.md")
        return _run_genes_species(args, parser)

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
          f"--trans 1 --loss 0.5 --output profiles trees -o {args.out}")
    return (f"wrote genes:species (key innovations) to {args.out}/ "
            f"({n_extant} extant tips, {model.n_drivers} drivers, tip prevalence {prev}) "
            f"in {dt:.3g} s")


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


def _run_genomes(tree: Tree, args: argparse.Namespace) -> str:
    """Simulate gene families along ``tree``, write output, and return a one-line summary.

    The default ``uniform`` model runs on the Rust engine automatically (``simulate_genomes``
    raises a build hint if the extension is missing); ``genome-wise`` runs on Python.
    """
    parts = set(Genomes.WRITE_PARTS) if "all" in args.output else set(args.output)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --output")

    if args.rate_model == "nucleotide":
        return _run_nucleotides(tree, args, parts)

    initial_size = 20 if args.initial_size is None else args.initial_size
    args.initial_size = initial_size          # record the effective value in the params log
    if args.rate_model == "genome-wise":
        model_kw = dict(rates=GenomeWiseRates(args.dup, args.trans, args.loss, args.orig))
    else:  # uniform
        model_kw = dict(duplication=args.dup, transfer=args.trans, loss=args.loss,
                        origination=args.orig)
    rate_kw = dict(**model_kw, initial_size=initial_size,
                   max_family_size=args.max_family_size, seed=args.seed)

    t0 = time.perf_counter()
    if parts == {"profiles"}:
        # counts-only Rust fast path: no genealogy reconstructed
        profiles = simulate_genomes(tree, output="profiles", **rate_kw)
        dt = time.perf_counter() - t0
        _write_profiles_only(args.out, tree, profiles, sparse=args.sparse)
        n_families = len(profiles.families)
    elif "trace" in parts and parts <= {"trace", "profiles"}:
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
    return (f"wrote [{' '.join(sorted(parts))}] to {args.out}/ "
            f"({len(tree.leaves())} tips, {n_families} gene families) in {dt:.3g} s")


def _run_nucleotides(tree: Tree, args: argparse.Namespace, parts: set) -> str:
    """Simulate nucleotide-resolution genomes (variable-length structural events) along ``tree``.

    Genes are not atomic here — they emerge as **atoms** (maximal intervals with one shared
    history). ``profiles`` writes the emergent atom-by-species profile (plus ``atoms.tsv`` and
    the per-leaf ``Mosaics.tsv``); ``trees`` writes the per-atom gene trees and their
    reconciliations. Only ``profiles``/``trees`` apply here (the family-model ``events`` /
    ``transfers`` / ``summary`` do not). ``profiles`` alone takes the fast Rust path.
    """
    want = parts & {"profiles", "trees", "ancestral"}
    if not want:
        raise ValueError("the nucleotide model writes 'profiles', 'trees' and/or 'ancestral'; "
                         "--output events/transfers/summary do not apply to it")
    ancestral = "ancestral" in want
    initial_size = 1 if args.initial_size is None else args.initial_size
    args.initial_size = initial_size          # record the effective value in the params log
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
    sim_kw = dict(inversion=args.inversion, loss=args.loss, duplication=args.dup,
                  transfer=args.trans, transposition=args.transposition,
                  origination=args.orig, root_length=args.root_length,
                  extension=args.extension, initial_size=initial_size, seed=args.seed,
                  gene_intervals=genes, pseudogenization=args.pseudogenization,
                  replacement=args.replacement, transfers=transfers, retain_internal=ancestral)

    t0 = time.perf_counter()
    if "trees" in want or genic or ancestral:  # genealogy / genic / ancestral need the Python engine
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
        _write_atoms_table(args.out, result.atoms)
        atom_ids, species, matrix = result.profile_matrix()
        pm = ProfileMatrix([f"atom{a}" for a in atom_ids], species, matrix)
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
        _write_atom_gene_trees(args.out, result, genic=genic)
        result.write_reconciliations(args.out)   # Reconciled_complete/extant.nwk + events.tsv
        if genic:
            _write_pseudogenizations(args.out, result)
    if ancestral:
        _write_ancestral(args.out, result, tree, args, gff_info)

    if gff_info is not None:
        print(f"  GFF {gff_info.seqid}: {gff_info.length} bp, {gff_info.n_features} genes "
              f"-> {len(gff_info.genes)} after trimming ({gff_info.n_trimmed} trimmed, "
              f"{gff_info.n_dropped} dropped as overlapping)")
    extra = f", {len(result.gene_atoms())} genes" if genic else ""
    return (f"wrote [{' '.join(sorted(want))}] (nucleotide{'/genic' if genic else ''}) to "
            f"{args.out}/ ({len(result.leaf_genomes)} tips, {len(result.atoms)} atoms{extra}) "
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
        lines = ["order\tatom\tkind\tgene_id\tstrand\tlength"]
        for i, (aid, strand) in enumerate(result.node_mosaic(node)):
            a = result._atom_by_id[aid]
            lines.append(f"{i}\tatom{aid}\t{a.kind}\t{a.gene_id or '-'}\t"
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
    lines = ["atom\tgene\tspecies_branch\ttime\tgene_lineage"]
    for atom_id, gene_id, species, t, gid in result.pseudogenizations():
        lines.append(f"atom{atom_id}\t{gene_id}\t{species}\t{t:.10g}\t{gid}")
    with open(os.path.join(out, "Pseudogenizations.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_atoms_table(out: str, atoms) -> None:
    """Write ``atoms.tsv`` — the emergent gene families (uncut ancestral intervals).

    Carries the ``kind`` (gene/intergene) and ``gene_id`` classification (``-`` for intergene).
    """
    lines = ["atom\tsource\tstart\tend\tlength\tkind\tgene_id"]
    for a in atoms:
        lines.append(f"atom{a.atom_id}\t{a.source}\t{a.start}\t{a.end}\t{a.length}\t"
                     f"{a.kind}\t{a.gene_id or '-'}")
    with open(os.path.join(out, "atoms.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_mosaics(out: str, result) -> None:
    """Write ``Mosaics.tsv`` — each extant genome as an ordered, signed sequence of atoms."""
    lines = ["leaf\tmosaic"]
    for leaf in sorted(result.leaf_genomes, key=lambda n: n.name):
        seq = " ".join(("+" if s > 0 else "-") + f"atom{aid}"
                       for aid, s in result.leaf_mosaic(leaf))
        lines.append(f"{leaf.name}\t{seq}")
    with open(os.path.join(out, "Mosaics.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_atom_gene_trees(out: str, result, genic: bool = False) -> None:
    """Write per-atom trees to ``atom<id>_complete.nwk`` / ``_extant.nwk``.

    Plain nucleotide model: everything under ``gene_trees/``. Genic mode: gene atoms under
    ``Gene_trees/`` and intergene atoms under ``Intergene_trees/`` (both tree sets recovered).
    """
    def dump(tdir: str, trees: dict) -> None:
        os.makedirs(tdir, exist_ok=True)
        for atom_id, (complete, extant) in trees.items():
            if complete:
                with open(os.path.join(tdir, f"atom{atom_id}_complete.nwk"), "w") as f:
                    f.write(complete + "\n")
            if extant:
                with open(os.path.join(tdir, f"atom{atom_id}_extant.nwk"), "w") as f:
                    f.write(extant + "\n")

    if genic:
        dump(os.path.join(out, "Gene_trees"), result.gene_trees())
        dump(os.path.join(out, "Intergene_trees"), result.intergene_trees())
    else:
        dump(os.path.join(out, "gene_trees"), result.atom_gene_trees())


def _add_abc_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("-t", "--tree", required=True,
                   help="species tree (Newick) the empirical data evolved along")
    p.add_argument("--profiles", required=True, metavar="TSV",
                   help="empirical copy-number profile table (families x species TSV, like the "
                        "Profiles.tsv that 'zombi2 genomes' writes)")
    # priors — reuse the genomes rate flags, but each takes a PRIOR: two values LOW HIGH
    # (a uniform prior) or one value (fixed). An omitted rate is held at 0.
    for flag, param in (("--dup", "duplication"), ("--trans", "transfer"),
                        ("--loss", "loss"), ("--orig", "origination")):
        p.add_argument(flag, type=float, nargs="+", default=None, metavar="RATE",
                       help=f"{param} prior: two values LOW HIGH (uniform) or one value (fixed); "
                            f"omit to hold {param} at 0")
    p.add_argument("--model", choices=("uniform", "family"), default="uniform",
                   help="uniform: one shared scalar rate per type (Rust; default); "
                        "family: per-family sampled rates, fitting each rate's mean (Python)")
    p.add_argument("--family-shape", type=float, default=2.0,
                   help="[--model family] Gamma shape for per-family rate dispersion (default 2.0)")
    p.add_argument("--n-sims", type=int, default=1000,
                   help="[rejection] number of prior simulations (default 1000)")
    p.add_argument("--accept", type=float, default=0.05,
                   help="[rejection] fraction of closest simulations to accept (default 0.05)")
    p.add_argument("--processes", type=int, default=None,
                   help="[rejection] parallel worker processes (default: serial)")
    p.add_argument("--smc", action="store_true",
                   help="use ABC-SMC (sequential, shrinking tolerance) instead of rejection")
    p.add_argument("--rounds", type=int, default=5, help="[--smc] number of SMC rounds (default 5)")
    p.add_argument("--particles", type=int, default=200,
                   help="[--smc] particles per round (default 200)")
    p.add_argument("--quantile", type=float, default=0.5,
                   help="[--smc] tolerance quantile carried between rounds (default 0.5)")
    p.add_argument("--regression-adjust", action="store_true",
                   help="also write the regression-adjusted posterior (Beaumont 2002)")
    p.add_argument("--initial-size", type=int, default=20,
                   help="gene families seeded at the root of each simulation (default 20)")
    p.add_argument("--max-family-size", type=_int_or_float, default=None,
                   help="growth cap for each simulation — recommended with --model family to "
                        "avoid runaway growth (integer = absolute, decimal = fraction of N)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


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
                  initial_size=args.initial_size, max_family_size=args.max_family_size,
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
    p.add_argument("--genomes", required=True, metavar="DIR",
                   help="a prior 'zombi2 genomes' output directory — reads its species_tree.nwk "
                        "and Events_trace.tsv (run genomes with 'trace' in --output)")
    p.add_argument("--family-speed", type=float, default=0.0, metavar="SIGMA",
                   help="per-gene-family intrinsic substitution speed: each family draws a "
                        "constant multiplier ~ LogNormal(0, SIGMA) (0 = every family the same)")
    p.add_argument("--branch-speed", type=float, default=0.0, metavar="SIGMA",
                   help="shared species-tree lineage clock: an autocorrelated lognormal relaxed "
                        "clock, drift SIGMA per sqrt(time) (0 = strict clock). Combine with "
                        "--family-speed for the full gene x lineage model; exclusive with --branch-bins")
    p.add_argument("--branch-bins", default=None, metavar="R1,R2,...",
                   help="alternative lineage clock — the discrete-bin within-branch GTDB model: "
                        "comma-separated ORDERED rate multipliers (e.g. 0.25,0.5,1,2,4), a Markov "
                        "walk between adjacent bins (--branch-switch-rate, --branch-up-bias)")
    p.add_argument("--branch-switch-rate", type=float, default=1.0, metavar="RATE",
                   help="[--branch-bins] rate of stepping to a neighbouring bin (default 1.0)")
    p.add_argument("--branch-up-bias", type=float, default=0.5, metavar="P",
                   help="[--branch-bins] probability a step goes to the faster neighbour "
                        "(default 0.5 = symmetric walk)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


def _run_sequence(args: argparse.Namespace) -> str:
    """Overlay the gene x lineage substitution clock on a prior genomes run's gene trees.

    Replays the compact ``Events_trace.tsv`` (no re-simulation of gene content), rescales every
    reconciled gene tree from time into substitutions/site, and writes the phylograms plus the
    drawn per-family speeds and per-branch rates. The lineage clock is shared across families
    (``--branch-speed`` lognormal or ``--branch-bins`` discrete-bin); each family draws one
    constant speed (``--family-speed``).
    """
    from .reconciliation import extant_species_from_records
    from .simulation import read_events_trace

    if args.family_speed < 0 or args.branch_speed < 0:
        raise ValueError("--family-speed / --branch-speed must be >= 0")
    if args.branch_speed > 0 and args.branch_bins:
        raise ValueError("--branch-speed (lognormal clock) and --branch-bins (discrete-bin "
                         "clock) are two lineage clocks; give at most one")

    tree_path = os.path.join(args.genomes, "species_tree.nwk")
    trace_path = os.path.join(args.genomes, "Events_trace.tsv")
    if not os.path.exists(trace_path):
        raise FileNotFoundError(
            f"{trace_path} not found — re-run 'zombi2 genomes' on that tree with 'trace' in "
            f"--output (e.g. --output trace profiles) so the genealogy can be replayed")
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
    phylo = se.scale_families(tree, families, gid2species, seed=args.seed)
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
    return (f"wrote substitution-unit gene trees for {n} families to {args.out}/gene_trees/ "
            f"(clock: {clock}, family-speed {args.family_speed}) in {dt:.3g} s")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="zombi2", description=_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)

    ps = sub.add_parser("species", help="simulate a species tree")
    _add_species_args(ps)

    pg = sub.add_parser("genomes", help="evolve gene families along a species tree")
    pg.add_argument("-t", "--tree", required=True,
                    help="input species tree in Newick format (e.g. species_tree.nwk)")
    _add_rate_args(pg)
    pg.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    pg.add_argument("-o", "--out", required=True, help="output directory")

    pt = sub.add_parser("trait", help="evolve a phenotypic trait along a given species tree")
    _add_trait_args(pt)

    pa = sub.add_parser("abc", help="fit gene-family rates to an empirical profile by ABC")
    _add_abc_args(pa)

    pce = sub.add_parser(
        "coevolve-genetrait",
        help="co-evolve a trait and gene families (trait-conditioned gene-family dynamics)")
    _add_coevolve_args(pce)

    pcv = sub.add_parser(
        "coevolve",
        help="co-evolve coupled processes (Phase 1: --couple traits:species = SSE)")
    _add_coevolve_mode_args(pcv)

    pq = sub.add_parser("sequence",
                        help="rescale a genomes run's gene trees into substitutions/site")
    _add_sequence_args(pq)

    args = parser.parse_args(argv)
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
                    f"lower --death, raise --max-attempts, or use --model backward.") from None
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
        summary = _run_genomes(tree, args)
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

    if args.command == "coevolve-genetrait":
        summary = _run_coevolve(args)
        print(summary)
        _write_params_log(os.path.join(args.out, "coevolve-genetrait.log"), args, summary)
        return 0

    if args.command == "coevolve":
        summary = _run_coevolve_mode(args, parser)
        print(summary)
        _write_params_log(os.path.join(args.out, "coevolve.log"), args, summary)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
