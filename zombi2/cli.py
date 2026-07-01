"""A thin command-line wrapper over the library.

The CLI only builds model objects and calls the simulator — no simulation logic lives
here (msprime idiom). v1 provides ``species`` (species tree only) and ``all`` (full
run). Standalone ``genomes`` (running against an externally supplied Newick tree) needs
a Newick reader and is deferred.
"""

from __future__ import annotations

import argparse

import numpy as np

from .rates import EventRates, RateModel
from .simulation import Simulation
from .species_model import SpeciesTreeModel
from .species_sim import SpeciesTreeSimulator


def _add_species_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--birth", type=float, required=True, help="speciation rate λ")
    p.add_argument("--death", type=float, default=0.0, help="extinction rate μ")
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

    pg = sub.add_parser("genomes", help="(deferred) genomes on a supplied tree")

    args = parser.parse_args(argv)

    if args.command == "species":
        model = SpeciesTreeModel(
            birth=args.birth, death=args.death, n_tips=args.tips,
            age=args.age, age_type=args.age_type,
        )
        rng = np.random.default_rng(args.seed)
        tree = SpeciesTreeSimulator().simulate(model, rng)
        out = args.out
        import os
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "species_tree.nwk"), "w") as f:
            f.write(tree.to_newick() + "\n")
        print(f"wrote {out}/species_tree.nwk ({args.tips} tips)")
        return 0

    if args.command == "all":
        species = SpeciesTreeModel(
            birth=args.birth, death=args.death, n_tips=args.tips,
            age=args.age, age_type=args.age_type,
        )
        rates = RateModel(EventRates(
            duplication=args.dup, transfer=args.trans, loss=args.loss, origination=args.orig,
        ))
        result = Simulation(species, rates, seed=args.seed, initial_size=args.initial_size).run()
        result.write(args.out)
        n_fam = len(result.profiles.families)
        print(f"wrote simulation to {args.out}/ ({args.tips} tips, {n_fam} gene families)")
        return 0

    if args.command == "genomes":
        raise NotImplementedError(
            "standalone `genomes` needs a Newick reader (deferred); use `all` or the "
            "Python API for now."
        )

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
