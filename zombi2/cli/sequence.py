"""zombi2 sequence command."""
from __future__ import annotations

import argparse
import os
import time

import numpy as np


from zombi2.distributions import LogNormal
from zombi2.sequences.clocks import (
    AutocorrelatedLogNormalClock, CIRClock, RateVariation, StrictClock,
    UncorrelatedGammaClock, UncorrelatedLogNormalClock, WhiteNoiseClock,
)
from zombi2.genomes.read_rates import read_family_speeds
from zombi2.sequences.evolution import SequenceEvolution
from zombi2.tree import read_newick

from zombi2.cli.framework import _add_params_arg, _write_params_log

def _add_sequence_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    _add_params_arg(g)
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
    g.add_argument("--family-speeds", metavar="FILE", dest="family_speeds",
                   help="TSV of explicit per-family substitution-speed multipliers (columns: "
                        "family speed) for named families — composes with (multiplies) the random "
                        "--family-speed draw; families not listed default to 1.0")

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
                             "give --subst-model to evolve DNA/protein/codon along the rescaled trees")
    g.add_argument("--subst-model", default=None, metavar="MODEL",
                   help="substitution model to simulate an alignment per family: DNA "
                        "(jc69/k80/hky85/gtr), protein (lg/wag/jtt/dayhoff/poisson) or codon "
                        "(gy94/mg94, dN/dS via --omega). DNA/protein/codon is auto-detected from "
                        "the name; codon models write in-frame coding DNA. Omit to only rescale "
                        "the trees (no sequences)")
    g.add_argument("--seq-length", type=int, default=300, metavar="N",
                   help="alignment length in sites (nt for DNA, aa for protein, CODONS for "
                        "gy94/mg94 → 3N nt; default 300); ignored when --root-fasta seeds the root")
    g.add_argument("--root-fasta", metavar="FILE", default=None,
                   help="FASTA (optionally .gz) of per-family root sequences keyed by family id "
                        "(header = family id); seeds each family's root instead of a random draw. "
                        "Its length overrides --seq-length per family")
    g.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="discrete-Gamma across-site rate heterogeneity shape (default: none)")
    g.add_argument("--kappa", type=float, default=2.0, metavar="K",
                   help="[DNA k80/hky85, codon gy94/mg94] transition/transversion ratio "
                        "(default 2.0)")
    g.add_argument("--omega", type=float, default=1.0, metavar="W",
                   help="[codon gy94/mg94] dN/dS ratio: <1 purifying, 1 neutral, >1 positive "
                        "selection (default 1.0)")
    g.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[DNA hky85/gtr, codon gy94/mg94] equilibrium base frequencies; for codon "
                        "models these build the F1×4 codon frequencies (default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[DNA gtr] the 6 exchangeabilities (default all 1)")

    g = p.add_argument_group("codon site models (dN/dS varies across sites)",
                             "with --subst-model gy94/mg94; let ω differ among codon sites. "
                             "Replaces the single --omega. Mutually exclusive with --gamma-shape")
    g.add_argument("--omega-model", default=None, metavar="MODEL",
                   choices=["m1a", "m2a", "m3", "m7", "m8"],
                   help="distribution of dN/dS across codon sites: m1a (neutral: ω0<1 + ω=1), "
                        "m2a (m1a + positive class ω>1), m3 (discrete classes via --omega-classes), "
                        "m7 (Beta(p,q) on [0,1]), m8 (m7 + positive class). Omit for a single --omega")
    g.add_argument("--omega0", type=float, default=None, metavar="W",
                   help="[m1a/m2a] the purifying class ω < 1 (default 0.1)")
    g.add_argument("--omega2", type=float, default=None, metavar="W",
                   help="[m2a] the positive-selection class ω > 1 (default 2.0)")
    g.add_argument("--omega-s", type=float, default=None, metavar="W", dest="omega_s",
                   help="[m8] the positive-selection class ω ≥ 1 (default 2.0)")
    g.add_argument("--omega-p0", type=float, default=None, metavar="P", dest="omega_p0",
                   help="[m1a/m2a] proportion of the purifying class; [m8] proportion of the Beta "
                        "classes (default 0.6 for m1a/m2a, 0.9 for m8)")
    g.add_argument("--omega-p1", type=float, default=None, metavar="P", dest="omega_p1",
                   help="[m2a] proportion of the neutral (ω=1) class (default 0.3)")
    g.add_argument("--beta-p", type=float, default=None, metavar="A", dest="beta_p",
                   help="[m7/m8] Beta(p,q) shape p for the ω-in-[0,1] classes")
    g.add_argument("--beta-q", type=float, default=None, metavar="B", dest="beta_q",
                   help="[m7/m8] Beta(p,q) shape q")
    g.add_argument("--omega-cats", type=int, default=4, metavar="N", dest="omega_cats",
                   help="[m7/m8] number of discrete Beta categories (default 4)")
    g.add_argument("--omega-classes", default=None, metavar="W0:P0,W1:P1,...", dest="omega_classes",
                   help="[m3] discrete ω classes as omega:proportion pairs (proportions renormalised)")

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

def _build_codon_site_model(args: argparse.Namespace):
    """Build a CodonSiteModel from the --omega-model flags (base = --subst-model gy94/mg94)."""
    from zombi2.sequences.codon_models import is_codon_model, make_codon_site_model

    if not is_codon_model(args.subst_model):
        raise ValueError(f"--omega-model {args.omega_model} needs a codon base model "
                         f"(--subst-model gy94 or mg94), not {args.subst_model!r}")
    if args.gamma_shape:
        raise ValueError("--omega-model (per-site ω mixture) and --gamma-shape (across-site rates) "
                         "are both per-site layers; give at most one")
    omegas = proportions = None
    if args.omega_classes:
        try:
            pairs = [tok.split(":") for tok in args.omega_classes.split(",") if tok.strip()]
            omegas = [float(w) for w, _ in pairs]
            proportions = [float(p) for _, p in pairs]
        except ValueError:
            raise ValueError("--omega-classes must be omega:proportion pairs, e.g. "
                             "'0.1:0.6,1.0:0.3,3.0:0.1'")
    return make_codon_site_model(
        args.omega_model, kappa=args.kappa, base=args.subst_model, freqs=args.base_freqs,
        p0=args.omega_p0, omega0=args.omega0, p1=args.omega_p1, omega2=args.omega2,
        beta_p=args.beta_p, beta_q=args.beta_q, omega_s=args.omega_s, ncat=args.omega_cats,
        omegas=omegas, proportions=proportions)

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
    from zombi2.sequences.models import (GammaRates, evolve_on_tree, is_codon_model, is_protein_model,
                               make_model, read_fasta, write_fasta)
    from zombi2.genomes.simulation import read_events_trace

    if args.family_speed < 0 or args.branch_speed < 0:
        raise ValueError("--family-speed / --branch-speed must be >= 0")
    if args.clock is None and args.branch_speed > 0 and args.branch_bins:
        raise ValueError("--branch-speed (lognormal clock) and --branch-bins (discrete-bin "
                         "clock) are two lineage clocks; give at most one, or select one "
                         "explicitly with --clock")
    lineage_clock, clock_desc = _build_lineage_clock(args)
    model = None
    if args.subst_model and args.omega_model:
        model = _build_codon_site_model(args)
    elif args.subst_model:
        model = make_model(args.subst_model, kappa=args.kappa, omega=args.omega,
                           freqs=args.base_freqs, rates=args.gtr_rates)
    elif args.omega_model:
        raise ValueError("--omega-model needs --subst-model gy94 or mg94 (the base codon model)")
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
    family_factors = read_family_speeds(args.family_speeds) if args.family_speeds else None
    se = SequenceEvolution(lineage=lineage_clock, family_speed=family_speed,
                           family_factors=family_factors)

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
    kind = ("codon DNA" if is_codon_model(args.subst_model)
            else "protein" if is_protein_model(args.subst_model) else "DNA")
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

    detail = model.name
    if getattr(model, "components", None) is not None:      # a codon site mixture
        detail = f"{model.name}, mean dN/dS={model.mean_omega:.3g}"
    return (f"{msg}; simulated {kind} alignments ({detail}) for {n_aln} families "
            f"to {args.out}/alignments/")

def _iter_leaves(node):
    """Yield the leaves (childless nodes) of a reconciliation ``_Node`` tree, left to right."""
    if not node.children:
        yield node
        return
    for child in node.children:
        yield from _iter_leaves(child)


def run(args, parser):
    summary = _run_sequence(args)
    print(summary)
    _write_params_log(os.path.join(args.out, "sequence.log"), args, summary)
    return 0
