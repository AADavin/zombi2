"""``zombi2 species`` — the dated species tree (a thin shell over
`zombi2.species.simulate_species_tree()`).

The long options are the API keyword names, and ``--birth`` / ``--death`` take a rate in its written
form (SPEC §5): a bare number on its natural scope (per lineage), or the same ``scope(base) ×
modifiers`` expression the Python API takes — ``--birth "1.0 * OnTime({0: 1.0, 3: 0.3})"``."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.species import IMPLEMENTED_MODIFIERS, _WRITE_OUTPUTS, simulate_species_tree
from zombi2._runtime.report import write_run_report
from zombi2.cli.framework import (resolve_seed, _add_flat_arg, _add_force_arg, _add_quiet_arg, _add_params_arg,
                                  _add_run_arg, _rate, _rates_help, _write_params_log,
                                  check_stale_downstream, clear_stale_downstream, defaults_used, signpost,
                                  input_digests,
                                  level_dir)

#: the RATES block for ``zombi2 species -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    IMPLEMENTED_MODIFIERS, "--birth", scopes="Global(1.0)",
    note="Global(base) is one budget for the tree (linear growth); a bare number is per lineage.")


def _add_species_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "where this run's files are written")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")

    # --birth and the stop condition are filled in run(), not marked argparse-`required`, so a
    # --params file can supply them (a required argument is never satisfied by a default) and so a
    # bare run can fall back to an illustrative default (recorded in the log and the report).
    g = p.add_argument_group("diversification")
    g.add_argument("--birth", type=_rate, default=None, metavar="RATE",
                   help="speciation rate, per lineage (default 1.0)")
    g.add_argument("--death", type=_rate, default=0.0, metavar="RATE",
                   help="extinction rate, per lineage (default 0 = pure birth)")

    g = p.add_argument_group("stop condition", "at most one; default --n-extant 20")
    g.add_argument("--n-extant", type=int, default=None, metavar="N", dest="n_extant",
                   help="stop at N extant lineages (conditioned on survival)")
    g.add_argument("--total-time", type=float, default=None, metavar="T", dest="total_time",
                   help="grow forward for T time units from the origin (t=0)")
    g.add_argument("--max-lineages", type=int, default=100_000, metavar="N", dest="max_lineages",
                   help="stop with an error if standing diversity passes N "
                        "(default 100000; 0 lifts the guard)"),

    g = p.add_argument_group("sampling & fossils")
    g.add_argument("--sampling", type=float, default=1.0, metavar="RHO",
                   help="extant sampling probability ρ, 0<ρ≤1 (default 1.0 = all)")
    g.add_argument("--fossils", type=float, default=0.0, metavar="RATE",
                   help="fossil (serial) recovery rate along the tree (default 0)")
    g.add_argument("--mass-extinction", action="append", nargs=2, type=float,
                   metavar=("TIME", "FRACTION"), default=None, dest="mass_extinction",
                   help="at TIME, lose each standing lineage with probability FRACTION; "
                        "repeatable; needs --total-time")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_WRITE_OUTPUTS, default=None, metavar="PART",
                   help=f"outputs to write (default all applicable): {', '.join(_WRITE_OUTPUTS)}")
    _add_flat_arg(g)
    _add_quiet_arg(g)
    _add_force_arg(g)


def run(args, parser):
    # validated here (not as argparse `required`) so a --params file can supply either
    if args.n_extant is not None and args.total_time is not None:
        parser.error("give exactly one stop condition: --n-extant N or --total-time T")
    # A bare `zombi2 species out/` runs: refusing taught nobody the shape of the command, and the
    # numbers are the part a newcomer has no way to guess. The chosen values are recorded in the log
    # and the report, so an invented default is visible there rather than shouted on the terminal.
    defaults_used(args, birth=1.0, **({} if args.total_time is not None else {"n_extant": 20}))

    # refuse up front if re-running would orphan a later level already in the run (unless --force)
    resolve_seed(args)                      # a run must be reproducible from its own log
    check_stale_downstream(args, "species")

    # [(time, fraction), ...] pulses, or None — the API places them on the timeline and needs a
    # fixed end (--total-time); it raises a clean error if that is missing.
    mass_ext = [(t, f) for t, f in args.mass_extinction] if args.mass_extinction else None

    t0 = time.perf_counter()
    result = simulate_species_tree(
        birth=args.birth, death=args.death, n_extant=args.n_extant, total_time=args.total_time,
        mass_extinctions=mass_ext, sampling=args.sampling, fossils=args.fossils, seed=args.seed,
        progress=not args.quiet, max_lineages=args.max_lineages or None)
    dt = time.perf_counter() - t0

    clear_stale_downstream(args, "species")   # --force: drop the now-stale downstream (run succeeded)
    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "species", args.flat)
    result.write(out, outputs=args.write)

    n_extant = result.n_extant
    n_total = len(result.complete_tree.nodes)
    n_leaves = len([n for n in result.complete_tree.nodes.values() if n.children is None])
    n_extinct = len(result.complete_tree.extinct_leaves())
    n_unsampled = len(result.complete_tree.unsampled_leaves())
    parts = [f"{n_extant} extant"]
    if n_extinct:
        parts.append(f"{n_extinct} extinct")
    # under --sampling these are survivors the run did not observe; without them the parts do not
    # add up to the tip count printed beside them, which reads as arithmetic going wrong
    if n_unsampled:
        parts.append(f"{n_unsampled} unsampled")
    if result.fossils:
        parts.append(f"{len(result.fossils)} fossils")
    summary = " + ".join(parts) + f" ({n_leaves} tips, {n_total} nodes)"
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(out, "species.log"), args, summary,
                      inputs=input_digests(args.birth, args.death))
    signpost(args, write_run_report(args.run), out)   # every file it wrote, then the run report
    return 0
