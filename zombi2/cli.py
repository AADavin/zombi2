"""A thin command-line wrapper over the library.

The CLI only builds model objects and calls the simulate functions — no simulation
logic lives here. Three subcommands:

* ``species`` — simulate a species tree only (backward birth–death).
* ``genomes`` — simulate gene families along an externally supplied Newick tree.
* ``all``     — simulate a species tree, then gene families along it, in one run.
"""

from __future__ import annotations

import argparse
import os

from .fast import rust_available, simulate_and_write_fast, simulate_profiles_fast
from .simulation import simulate_genomes
from .species_model import BirthDeath
from .species_sim import simulate_species_tree
from .tree import Tree, read_newick


def _int_or_float(text: str) -> int | float:
    """Parse ``--max-family-size``: a plain integer is an absolute cap, a value with a
    decimal point is a fraction of the number of species (e.g. ``0.5`` -> half of N)."""
    return float(text) if "." in text else int(text)


def _add_species_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--birth", type=float, required=True, help="speciation rate")
    p.add_argument("--death", type=float, default=0.0, help="extinction rate")
    p.add_argument("--tips", type=int, required=True, help="number of extant species N")
    p.add_argument("--age", type=float, required=True, help="tree age (crown or stem)")
    p.add_argument("--age-type", choices=("crown", "stem"), default="crown")
    p.add_argument("--seed", type=int, default=None)
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
    p.add_argument("--fast", action="store_true",
                   help="use the Rust engine (much faster). Writes the same full ZOMBI-1 "
                        "output as the default, but simulated, reconstructed and written "
                        "entirely in Rust. Requires the compiled zombi2_core extension.")
    p.add_argument("--profiles-only", action="store_true",
                   help="with --fast, write only species_tree.nwk + Profiles.tsv/Presence.tsv "
                        "(no event log or gene trees) — the fastest path.")


def _write_profiles_only(out: str, tree: Tree, profiles) -> None:
    """Emit the reduced output of the Rust fast path: tree + copy-number/presence matrices."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(out, "Profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())
    with open(os.path.join(out, "Presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))


def _run_genomes(tree: Tree, args) -> str:
    """Simulate gene families along ``tree`` and write output; return a one-line summary."""
    if args.fast:
        if not rust_available():
            raise SystemExit(
                "--fast needs the compiled zombi2_core extension; build it with "
                "`cd rust && maturin build --release -i python3 && "
                "pip install --force-reinstall target/wheels/*.whl`, or drop --fast."
            )
        if args.profiles_only:
            profiles = simulate_profiles_fast(
                tree, duplication=args.dup, transfer=args.trans, loss=args.loss,
                origination=args.orig, initial_size=args.initial_size,
                max_family_size=args.max_family_size, seed=args.seed,
            )
            _write_profiles_only(args.out, tree, profiles)
            return (f"wrote profiles to {args.out}/ (Rust fast path: "
                    f"{len(tree.leaves())} tips, {len(profiles.families)} gene families, "
                    f"profiles only)")
        summary = simulate_and_write_fast(
            tree, args.out, duplication=args.dup, transfer=args.trans, loss=args.loss,
            origination=args.orig, initial_size=args.initial_size,
            max_family_size=args.max_family_size, seed=args.seed,
        )
        return (f"wrote simulation to {args.out}/ (Rust fast path: "
                f"{summary['n_species']} tips, {summary['n_families']} gene families, "
                f"{summary['n_events']} events)")

    genomes = simulate_genomes(
        tree, duplication=args.dup, transfer=args.trans, loss=args.loss,
        origination=args.orig, initial_size=args.initial_size,
        max_family_size=args.max_family_size, seed=args.seed,
    )
    genomes.write(args.out)
    return (f"wrote simulation to {args.out}/ "
            f"({len(tree.leaves())} tips, {len(genomes.profiles.families)} gene families)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="zombi2", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    ps = sub.add_parser("species", help="simulate a species tree only")
    _add_species_args(ps)

    pg = sub.add_parser("genomes", help="simulate gene families along a supplied Newick tree")
    pg.add_argument("-t", "--tree", required=True,
                    help="input species tree in Newick format (e.g. species_tree.nwk)")
    _add_rate_args(pg)
    pg.add_argument("--seed", type=int, default=None)
    pg.add_argument("-o", "--out", required=True, help="output directory")

    pa = sub.add_parser("all", help="simulate species tree then gene families")
    _add_species_args(pa)
    _add_rate_args(pa)

    args = parser.parse_args(argv)

    if args.command == "species":
        tree = simulate_species_tree(
            BirthDeath(args.birth, args.death),
            n_tips=args.tips, age=args.age, age_type=args.age_type, seed=args.seed,
        )
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
            f.write(tree.to_newick() + "\n")
        print(f"wrote {args.out}/species_tree.nwk ({args.tips} tips)")
        return 0

    if args.command == "genomes":
        with open(args.tree) as f:
            tree = read_newick(f.read())
        print(_run_genomes(tree, args))
        return 0

    if args.command == "all":
        tree = simulate_species_tree(
            BirthDeath(args.birth, args.death),
            n_tips=args.tips, age=args.age, age_type=args.age_type, seed=args.seed,
        )
        print(_run_genomes(tree, args))
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
