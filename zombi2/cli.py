"""A thin command-line wrapper over the library.

The CLI only builds model objects and calls the simulate functions — no simulation
logic lives here. v1 provides ``species`` (species tree only) and ``all`` (full run).
Standalone ``genomes`` (on an externally supplied Newick tree) needs a Newick reader
and is deferred.
"""

from __future__ import annotations

import argparse

from .simulation import simulate_genomes
from .species_model import BirthDeath
from .species_sim import simulate_species_tree


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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="zombi2", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    ps = sub.add_parser("species", help="simulate a species tree only")
    _add_species_args(ps)

    pa = sub.add_parser("all", help="simulate species tree then gene families")
    _add_species_args(pa)
    _add_rate_args(pa)

    sub.add_parser("genomes", help="(deferred) genomes on a supplied tree")

    args = parser.parse_args(argv)

    if args.command == "species":
        tree = simulate_species_tree(
            BirthDeath(args.birth, args.death),
            n_tips=args.tips, age=args.age, age_type=args.age_type, seed=args.seed,
        )
        import os
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
            f.write(tree.to_newick() + "\n")
        print(f"wrote {args.out}/species_tree.nwk ({args.tips} tips)")
        return 0

    if args.command == "all":
        tree = simulate_species_tree(
            BirthDeath(args.birth, args.death),
            n_tips=args.tips, age=args.age, age_type=args.age_type, seed=args.seed,
        )
        genomes = simulate_genomes(
            tree, duplication=args.dup, transfer=args.trans, loss=args.loss,
            origination=args.orig, initial_size=args.initial_size, seed=args.seed,
        )
        genomes.write(args.out)
        print(f"wrote simulation to {args.out}/ "
              f"({args.tips} tips, {len(genomes.profiles.families)} gene families)")
        return 0

    if args.command == "genomes":
        raise NotImplementedError(
            "standalone `genomes` needs a Newick reader (deferred); use `all` or the "
            "Python API for now."
        )

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
