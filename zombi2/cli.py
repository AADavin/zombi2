"""Command-line interface for ZOMBI2 (``zombi2 species`` and ``zombi2 genomes``)."""

from __future__ import annotations

import argparse
import os

from .simulation import simulate_genomes
from .species_model import BirthDeath
from .species_sim import simulate_species_tree
from .sse import BiSSE, QuaSSE, simulate_sse
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
  zombi2 sse       simulate a tree and a trait jointly (the trait drives diversification)

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
                        "lineages")
    p.add_argument("--birth", type=float, default=1.0, help="speciation rate (default: 1.0)")
    p.add_argument("--death", type=float, default=0.3, help="extinction rate (default: 0.3)")
    p.add_argument("--tips", type=int, default=None,
                   help="number of extant species (backward default: 50; "
                        "forward: give --tips OR --age)")
    p.add_argument("--age", type=float, default=None,
                   help="tree age / timescale, in the same time units as the rates "
                        "(backward default: 1.0; forward: give --tips OR --age)")
    p.add_argument("--age-type", choices=("crown", "stem"), default="crown",
                   help="interpret --age as crown (default) or stem age [backward]")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


def _add_rate_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--dup", type=float, default=0.0, help="duplication rate")
    p.add_argument("--trans", type=float, default=0.0, help="transfer rate")
    p.add_argument("--loss", type=float, default=0.0, help="loss rate")
    p.add_argument("--orig", type=float, default=0.0, help="origination rate")
    p.add_argument("--initial-size", type=int, default=20, help="seed gene families at root")
    p.add_argument("--max-family-size", type=_int_or_float, default=None,
                   help="bound family growth: integer = absolute cap, "
                        "decimal = fraction of the number of species (e.g. 0.5)")
    p.add_argument("--profiles-only", action="store_true",
                   help="write only species_tree.nwk + Profiles.tsv/Presence.tsv (no event "
                        "log or gene trees) — the fastest path (Rust counts-only engine).")


def _add_trait_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("-t", "--tree", required=True,
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    p.add_argument("--model", choices=("bm", "ou", "eb", "mk", "threshold"), default="bm",
                   help="trait model: bm=Brownian motion, ou=Ornstein-Uhlenbeck, "
                        "eb=early burst/ACDC, mk=discrete k-state, threshold (default: bm)")
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
    p.add_argument("--thresholds", default="0.0",
                   help="comma-separated liability cut points [threshold] (default: 0.0)")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


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
        return Mk.equal_rates(args.states, args.rate)
    # threshold
    thresholds = [float(t) for t in str(args.thresholds).split(",")]
    return ThresholdModel(thresholds=thresholds, sigma2=args.sigma2,
                          x0=(0.0 if x0 is None else x0))


def _write_trait_output(out: str, res, *, tree_file: str | None = None) -> None:
    """Write traits.tsv (tips + ancestral node values) and trait_tree.nwk (annotated Newick).

    ``tree_file`` names the file to write the tree into (``sse`` writes its own simulated
    tree); ``None`` for ``trait`` (the tree was supplied on the command line).
    """
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))            # every node: tips AND ancestral states
    with open(os.path.join(out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")             # values annotated on every node too
    if tree_file is not None:
        with open(os.path.join(out, tree_file), "w") as f:
            f.write(res.tree.to_newick() + "\n")


def _run_trait(args) -> str:
    """Evolve a trait along the supplied species tree and write the output folder."""
    with open(args.tree) as f:
        tree = read_newick(f.read())
    res = simulate_traits(tree, _build_trait_model(args), seed=args.seed)
    _write_trait_output(args.out, res)
    return (f"wrote traits to {args.out}/ (model={args.model}; "
            f"{len(tree.extant_leaves())} tip + {len(tree.internal_nodes())} ancestral values)")


def _add_sse_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--model", choices=("bisse", "quasse"), default="bisse",
                   help="bisse: a binary trait drives speciation/extinction; "
                        "quasse: a continuous trait drives it (default: bisse)")
    p.add_argument("--age", type=float, default=None,
                   help="grow for this crown age (default: 3.0 if neither --age nor --tips)")
    p.add_argument("--tips", type=int, default=None,
                   help="grow until this many extant species instead")
    # BiSSE: per-state speciation/extinction and transition rates
    p.add_argument("--lambda0", type=float, default=1.0, help="[bisse] speciation in state 0")
    p.add_argument("--lambda1", type=float, default=2.0, help="[bisse] speciation in state 1")
    p.add_argument("--mu0", type=float, default=0.3, help="[bisse] extinction in state 0")
    p.add_argument("--mu1", type=float, default=0.3, help="[bisse] extinction in state 1")
    p.add_argument("--q01", type=float, default=0.1, help="[bisse] transition rate 0 -> 1")
    p.add_argument("--q10", type=float, default=0.1, help="[bisse] transition rate 1 -> 0")
    # QuaSSE: sigmoidal speciation in the trait + constant extinction
    p.add_argument("--spec-low", type=float, default=0.2,
                   help="[quasse] speciation rate at low trait values")
    p.add_argument("--spec-high", type=float, default=2.0,
                   help="[quasse] speciation rate at high trait values")
    p.add_argument("--spec-center", type=float, default=0.0,
                   help="[quasse] trait value at the middle of the speciation sigmoid")
    p.add_argument("--spec-slope", type=float, default=1.0,
                   help="[quasse] steepness of the speciation sigmoid")
    p.add_argument("--mu", type=float, default=0.1, help="[quasse] constant extinction rate")
    p.add_argument("--sigma2", type=float, default=1.0, help="[quasse] trait diffusion rate")
    p.add_argument("--x0", type=float, default=0.0, help="[quasse] root trait value")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    p.add_argument("-o", "--out", required=True, help="output directory")


def _run_sse(args, parser) -> str:
    """Simulate a tree and a trait jointly (BiSSE/QuaSSE) and write the output folder."""
    if args.age is not None and args.tips is not None:
        parser.error("provide at most one of --age or --tips")
    age = args.age
    tips = args.tips
    if age is None and tips is None:
        age = 3.0
    if args.model == "bisse":
        model = BiSSE(args.lambda0, args.lambda1, args.mu0, args.mu1, args.q01, args.q10)
    else:  # quasse: sigmoidal speciation, constant extinction; rate_bound = sup(lambda) + mu
        spec = QuaSSE.sigmoid(args.spec_low, args.spec_high, args.spec_center, args.spec_slope)
        mu = args.mu
        model = QuaSSE(spec, lambda x: mu, sigma2=args.sigma2,
                       rate_bound=args.spec_high + mu, x0=args.x0)
    res = simulate_sse(model, age=age, n_tips=tips, seed=args.seed)
    _write_trait_output(args.out, res, tree_file="species_tree.nwk")
    n_extant = len(res.tree.extant_leaves())
    n_anc = len(res.tree.internal_nodes())
    return (f"wrote {args.model} tree + traits to {args.out}/ "
            f"({n_extant} extant tips + {n_anc} ancestral values)")


def _write_profiles_only(out: str, tree: Tree, profiles) -> None:
    """Emit the reduced profiles-only output: tree + copy-number/presence matrices."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(out, "Profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())
    with open(os.path.join(out, "Presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))


def _run_genomes(tree: Tree, args) -> str:
    """Simulate gene families along ``tree`` and write output; return a one-line summary.

    The built-in model runs on the Rust engine automatically (``simulate_genomes`` raises a
    build hint if the extension is missing).
    """
    rate_kw = dict(duplication=args.dup, transfer=args.trans, loss=args.loss,
                   origination=args.orig, initial_size=args.initial_size,
                   max_family_size=args.max_family_size, seed=args.seed)

    if args.profiles_only:
        profiles = simulate_genomes(tree, output="profiles", **rate_kw)
        _write_profiles_only(args.out, tree, profiles)
        return (f"wrote profiles to {args.out}/ "
                f"({len(tree.leaves())} tips, {len(profiles.families)} gene families, "
                f"profiles only)")

    genomes = simulate_genomes(tree, **rate_kw)
    genomes.write(args.out)
    return (f"wrote simulation to {args.out}/ "
            f"({len(tree.leaves())} tips, {len(genomes.profiles.families)} gene families)")


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
    pg.add_argument("--seed", type=int, default=None)
    pg.add_argument("-o", "--out", required=True, help="output directory")

    pt = sub.add_parser("trait", help="evolve a phenotypic trait along a given species tree")
    _add_trait_args(pt)

    pe = sub.add_parser("sse", help="simulate a tree and a trait jointly (BiSSE / QuaSSE)")
    _add_sse_args(pe)

    args = parser.parse_args(argv)

    if args.command == "species":
        model = BirthDeath(args.birth, args.death)
        if args.model == "backward":
            n_tips = args.tips if args.tips is not None else 50
            age = args.age if args.age is not None else 1.0
            tree = simulate_species_tree(model, n_tips=n_tips, age=age,
                                         age_type=args.age_type, direction="backward",
                                         seed=args.seed)
        else:  # forward
            if (args.tips is None) == (args.age is None):
                parser.error("forward model needs exactly one of --tips or --age "
                             "(--tips to stop at that many extant species; "
                             "--age to grow for that long)")
            tree = simulate_species_tree(model, n_tips=args.tips, age=args.age,
                                         direction="forward", seed=args.seed)

        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
            f.write(tree.to_newick() + "\n")
        leaves = tree.leaves()
        n_extant = sum(1 for n in leaves if n.is_extant)
        extinct = len(leaves) - n_extant
        extra = f" + {extinct} extinct" if extinct else ""
        print(f"wrote {args.out}/species_tree.nwk ({n_extant} extant{extra} tips)")
        return 0

    if args.command == "genomes":
        with open(args.tree) as f:
            tree = read_newick(f.read())
        print(_run_genomes(tree, args))
        return 0

    if args.command == "trait":
        print(_run_trait(args))
        return 0

    if args.command == "sse":
        print(_run_sse(args, parser))
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
