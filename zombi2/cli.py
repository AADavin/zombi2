"""Command-line interface for ZOMBI2 (``zombi2 species`` and ``zombi2 genomes``)."""

from __future__ import annotations

import argparse
import os
import sys
import time

from .ghosts import add_ghost_lineages
from .rates import GenomeWiseRates
from .simulation import simulate_genomes
from .species_model import BirthDeath, EpisodicBirthDeath
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
    p.add_argument("--log-level", choices=("low", "medium", "high"), default="medium",
                   help="detail of the parameters saved to <out>/species_tree.log "
                        "(always written; default medium)")
    p.add_argument("-o", "--out", required=True, help="output directory")


def _add_rate_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--rate-model", choices=("uniform", "genome-wise"), default="uniform",
                   help="uniform: same per-copy rates for every family (Rust; default); "
                        "genome-wise: constant per-genome rates, linear growth (Python)")
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


def _write_params_log(path: str, args: argparse.Namespace, summary: str, level: str) -> None:
    """Write the run's parameters to ``path`` at verbosity ``level`` (low/medium/high).

    ``low`` = version + command line + seed + result (enough to reproduce); ``medium`` (the
    default) adds the core scientific parameters; ``high`` adds a timestamp and every argument
    (including engine internals).
    """
    from . import __version__
    d = vars(args)
    lines = ["# ZOMBI2 run parameters",
             f"zombi2_version\t{__version__}",
             f"command_line\t{' '.join(sys.argv)}",
             f"seed\t{d.get('seed')}"]
    if level == "high":
        import datetime
        lines.append(f"timestamp\t{datetime.datetime.now().isoformat(timespec='seconds')}")
    if level in ("medium", "high"):
        skip = {"log_level", "out", "command", "seed"}
        if level == "medium":
            skip |= {"max_attempts", "max_lineages"}   # engine internals -> high only
        for key in sorted(d):
            if key not in skip:
                lines.append(f"{key}\t{d[key]}")
    lines.append(f"result\t{summary}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_profiles_only(out: str, tree: Tree, profiles) -> None:
    """Emit the reduced profiles-only output: tree + copy-number/presence matrices."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(out, "Profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())
    with open(os.path.join(out, "Presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))


def _run_genomes(tree: Tree, args: argparse.Namespace) -> str:
    """Simulate gene families along ``tree``, write output, and return a one-line summary.

    The default ``uniform`` model runs on the Rust engine automatically (``simulate_genomes``
    raises a build hint if the extension is missing); ``genome-wise`` runs on Python.
    """
    if args.rate_model == "genome-wise":
        model_kw = dict(rates=GenomeWiseRates(args.dup, args.trans, args.loss, args.orig))
    else:  # uniform
        model_kw = dict(duplication=args.dup, transfer=args.trans, loss=args.loss,
                        origination=args.orig)
    rate_kw = dict(**model_kw, initial_size=args.initial_size,
                   max_family_size=args.max_family_size, seed=args.seed)

    t0 = time.perf_counter()
    if args.profiles_only:
        profiles = simulate_genomes(tree, output="profiles", **rate_kw)
        dt = time.perf_counter() - t0
        _write_profiles_only(args.out, tree, profiles)
        return (f"wrote profiles to {args.out}/ "
                f"({len(tree.leaves())} tips, {len(profiles.families)} gene families, "
                f"profiles only) in {dt:.3g} s")

    genomes = simulate_genomes(tree, **rate_kw)
    dt = time.perf_counter() - t0
    genomes.write(args.out)
    return (f"wrote simulation to {args.out}/ "
            f"({len(tree.leaves())} tips, {len(genomes.profiles.families)} gene families) "
            f"in {dt:.3g} s")


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
    pg.add_argument("--log-level", choices=("low", "medium", "high"), default="medium",
                    help="detail of the parameters saved to <out>/genomes.log "
                         "(always written; default medium)")
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
        _write_params_log(os.path.join(args.out, "species_tree.log"), args, summary,
                          args.log_level)
        return 0

    if args.command == "genomes":
        with open(args.tree) as f:
            tree = read_newick(f.read())
        summary = _run_genomes(tree, args)
        print(summary)
        _write_params_log(os.path.join(args.out, "genomes.log"), args, summary, args.log_level)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
