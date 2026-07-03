"""Command-line interface for ZOMBI2 (``zombi2 species`` and ``zombi2 genomes``)."""

from __future__ import annotations

import argparse
import os
import sys

from .simulation import simulate_genomes
from .species_model import BirthDeath
from .species_sim import simulate_species_tree
from .tree import Tree, read_newick

_DESCRIPTION = """\
ZOMBI2 — a phylogenetic simulator of species trees and gene families.

Simulate in two steps: build a species tree, then evolve gene families along it.

  zombi2 species   simulate a species tree
  zombi2 genomes   evolve gene families along a species tree (Newick)

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
    p.add_argument("--max-attempts", type=int, default=10000,
                   help="[forward] retries before giving up when the process goes extinct "
                        "(default: 10000)")
    p.add_argument("--max-lineages", type=int, default=1_000_000,
                   help="[forward] abort a run that exceeds this many live lineages "
                        "(default: 1000000)")
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

    args = parser.parse_args(argv)
    try:
        return _dispatch(args, parser)
    except (ValueError, RuntimeError, FileNotFoundError, OSError) as e:
        # Report expected failures as a clean one-line error, never a traceback.
        print(f"zombi2: error: {e}", file=sys.stderr)
        return 1


def _dispatch(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.command == "species":
        model = BirthDeath(args.birth, args.death)
        common = dict(age_type=args.age_type, max_attempts=args.max_attempts,
                      max_lineages=args.max_lineages, seed=args.seed)
        if args.model == "backward":
            n_tips = args.tips if args.tips is not None else 50
            age = args.age if args.age is not None else 1.0
            tree = simulate_species_tree(model, n_tips=n_tips, age=age,
                                         direction="backward", **common)
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

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
