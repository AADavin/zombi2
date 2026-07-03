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
from .rates import GenomeWiseRates
from .simulation import Genomes, simulate_genomes
from .species_model import BirthDeath, EpisodicBirthDeath
from .species_sim import simulate_species_tree
from .traits import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst, Mk, ThresholdModel, simulate_traits,
)
from .tree import Tree, read_newick

_DESCRIPTION = """\
ZOMBI2 — a phylogenetic simulator of species trees and gene families.

Simulate in two steps: build a species tree, then evolve gene families along it.

  zombi2 species   simulate a species tree
  zombi2 genomes   evolve gene families along a species tree (Newick)
  zombi2 trait     evolve a phenotypic trait along a given species tree
  zombi2 abc       fit gene-family rates to an empirical profile (ABC inference)

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
    p.add_argument("--birth", type=float, nargs="+", default=[1.0], metavar="RATE",
                   help="speciation rate (default 1.0); several values with --shifts give an "
                        "episodic (skyline) model")
    p.add_argument("--death", type=float, nargs="+", default=[0.3], metavar="RATE",
                   help="extinction rate (default 0.3); several values with --shifts give an "
                        "episodic (skyline) model")
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
                   choices=(*Genomes.WRITE_PARTS, "all"), default=["profiles", "trees"],
                   help="which output files to write — any of {profiles, trees, events, "
                        "transfers, summary} or 'all' (default: profiles trees). "
                        "species_tree.nwk is always written; 'profiles' alone takes the fast "
                        "Rust counts-only path")
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


def _build_species_model(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Construct a BirthDeath or EpisodicBirthDeath model from the CLI args (validated)."""
    if args.model == "backward" and (args.fossilization or args.removal != 1.0
                                     or args.sampling_fraction != 1.0):
        parser.error("--fossilization / --removal / --sampling-fraction require --model forward "
                     "(the backward reconstructed sampler assumes complete sampling)")
    episodic = args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1
    if not episodic:
        return BirthDeath(args.birth[0], args.death[0], fossilization=args.fossilization,
                          sampling_fraction=args.sampling_fraction, removal=args.removal)
    shifts = args.shifts or []
    if len(args.birth) != len(args.death) or len(shifts) != len(args.birth) - 1:
        parser.error("episodic model needs len(--birth) == len(--death) == len(--shifts)+1 "
                     f"(got {len(args.birth)} birth, {len(args.death)} death, {len(shifts)} shifts)")
    return EpisodicBirthDeath(birth=args.birth, death=args.death, shifts=shifts,
                              fossilization=(args.fossilization or None),
                              sampling_fraction=args.sampling_fraction, removal=args.removal)


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
    want = parts & {"profiles", "trees"}
    if not want:
        raise ValueError("the nucleotide model writes 'profiles' and/or 'trees'; "
                         "--output events/transfers/summary do not apply to it")
    initial_size = 1 if args.initial_size is None else args.initial_size
    args.initial_size = initial_size          # record the effective value in the params log
    sim_kw = dict(inversion=args.inversion, loss=args.loss, duplication=args.dup,
                  transfer=args.trans, transposition=args.transposition,
                  origination=args.orig, root_length=args.root_length,
                  extension=args.extension, initial_size=initial_size, seed=args.seed)

    t0 = time.perf_counter()
    if "trees" in want:                       # genealogy needs the full Python event log
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
        _write_atom_gene_trees(args.out, result)
        result.write_reconciliations(args.out)   # Reconciled_complete/extant.nwk + events.tsv

    return (f"wrote [{' '.join(sorted(want))}] (nucleotide) to {args.out}/ "
            f"({len(result.leaf_genomes)} tips, {len(result.atoms)} atoms) in {dt:.3g} s")


def _write_atoms_table(out: str, atoms) -> None:
    """Write ``atoms.tsv`` — the emergent gene families (uncut ancestral intervals)."""
    lines = ["atom\tsource\tstart\tend\tlength"]
    for a in atoms:
        lines.append(f"atom{a.atom_id}\t{a.source}\t{a.start}\t{a.end}\t{a.length}")
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


def _write_atom_gene_trees(out: str, result) -> None:
    """Write per-atom gene trees to ``gene_trees/atom<id>_complete.nwk`` / ``_extant.nwk``."""
    tdir = os.path.join(out, "gene_trees")
    os.makedirs(tdir, exist_ok=True)
    for atom_id, (complete, extant) in result.atom_gene_trees().items():
        if complete:
            with open(os.path.join(tdir, f"atom{atom_id}_complete.nwk"), "w") as f:
                f.write(complete + "\n")
        if extant:
            with open(os.path.join(tdir, f"atom{atom_id}_extant.nwk"), "w") as f:
                f.write(extant + "\n")


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

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
