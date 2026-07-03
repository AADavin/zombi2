"""Command-line interface for ZOMBI2 (``zombi2 species`` and ``zombi2 genomes``)."""

from __future__ import annotations

import argparse
import os

from .simulation import simulate_genomes
from .species_model import BirthDeath
from .species_sim import simulate_species_tree
from .traits import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst, Mk, ThresholdModel,
    simulate_traits, replicate_traits,
)
from .tree import Tree, read_newick

_DESCRIPTION = """\
ZOMBI2 — a phylogenetic simulator of species trees and gene families.

Simulate in two steps: build a species tree, then evolve gene families along it.

  zombi2 species   simulate a species tree
  zombi2 genomes   evolve gene families along a species tree (Newick)
  zombi2 trait     evolve a phenotypic trait along a given species tree

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
    p.add_argument("--replicates", type=int, default=1,
                   help="simulate the trait this many times with the same parameters; writes "
                        "traits.tsv with one column per replicate (default: 1)")
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


def _fmt_cell(value) -> str:
    """Format one trait value for a TSV cell (float -> 6 sig figs, else str)."""
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
    """Evolve a trait along the supplied species tree and write the output folder."""
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = _build_trait_model(args)
    os.makedirs(args.out, exist_ok=True)

    if args.replicates > 1:
        results = replicate_traits(tree, model, args.replicates, seed=args.seed)
        with open(os.path.join(args.out, "traits.tsv"), "w") as f:
            f.write(_replicate_table(results))       # one column per replicate, all nodes
        return (f"wrote {args.replicates} trait replicates to {args.out}/traits.tsv "
                f"(model={args.model}; {len(tree.nodes())} nodes x {args.replicates} columns)")

    res = simulate_traits(tree, model, seed=args.seed)
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))             # every node: tips AND ancestral states
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")              # values annotated on every node too
    return (f"wrote traits to {args.out}/ (model={args.model}; "
            f"{len(tree.extant_leaves())} tip + {len(tree.internal_nodes())} ancestral values)")


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

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
