"""``zombi2 species`` — the dated species tree (a thin shell over
:func:`zombi2.species.simulate_species_tree`).

The long options are the API keyword names, and ``--birth`` / ``--death`` take a rate in its written
form (SPEC §5): a bare number on its natural scope (per lineage), or the same ``scope(base) ×
modifiers`` expression the Python API takes — ``--birth "1.0 * OnTime({0: 1.0, 3: 0.3})"``."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.species import WIRED_MODIFIERS, _WRITE_OUTPUTS, simulate_species_tree
from zombi2.cli.framework import (_add_flat_arg, _add_params_arg, _add_run_arg, _rate,
                                  _rates_help, _write_params_log, level_dir)

#: the RATES block for ``zombi2 species -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--birth", scopes="Global(1.0)",
    note="Global(base) is one shared budget for the whole tree (linear, not exponential, growth); "
         "a bare number is per lineage.")


def _add_species_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "where this run's files are written")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")

    # --birth and the stop condition are validated in run(), not marked argparse-`required`, so a
    # --params file can supply them (a required argument is never satisfied by a default).
    g = p.add_argument_group("diversification", "the per-lineage birth–death rates (see RATES below)")
    g.add_argument("--birth", type=_rate, default=None, metavar="RATE",
                   help="speciation rate (per lineage) — required")
    g.add_argument("--death", type=_rate, default=0.0, metavar="RATE",
                   help="extinction rate (per lineage); 0 = a pure-birth (Yule) tree (default 0)")

    g = p.add_argument_group("stop condition", "grow the tree until exactly one of these — required")
    g.add_argument("--n-extant", type=int, default=None, metavar="N", dest="n_extant",
                   help="stop at N extant (surviving) lineages — conditioned on survival")
    g.add_argument("--total-time", type=float, default=None, metavar="T", dest="total_time",
                   help="grow forward for T time units (time runs forward from the crown)")

    g = p.add_argument_group("sampling & fossils")
    g.add_argument("--sampling", type=float, default=1.0, metavar="RHO",
                   help="incomplete extant sampling ρ, 0<ρ≤1: each survivor is observed with "
                        "probability ρ (default 1.0 = all observed)")
    g.add_argument("--fossils", type=float, default=0.0, metavar="RATE",
                   help="fossil (serial) recovery rate along the tree (default 0 = no fossils)")
    g.add_argument("--mass-extinction", action="append", nargs=2, type=float,
                   metavar=("TIME", "FRACTION"), default=None, dest="mass_extinction",
                   help="a mass-extinction pulse: at TIME (forward from the crown) each standing "
                        "lineage is lost with probability FRACTION. Repeat for several pulses, e.g. "
                        "--mass-extinction 3.0 0.75. Needs --total-time.")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_WRITE_OUTPUTS, default=None, metavar="PART",
                   help=f"which outputs to write (default all applicable): "
                        f"{', '.join(_WRITE_OUTPUTS)}. Files are prefixed 'species_'.")
    _add_flat_arg(g)


def run(args, parser):
    # validated here (not as argparse `required`) so a --params file can supply either
    if args.birth is None:
        parser.error("--birth is required (give it on the command line or in --params)")
    if (args.n_extant is None) == (args.total_time is None):
        parser.error("give exactly one stop condition: --n-extant N or --total-time T")

    # [(time, fraction), ...] pulses, or None — the API places them on the timeline and needs a
    # fixed end (--total-time); it raises a clean error if that is missing.
    mass_ext = [(t, f) for t, f in args.mass_extinction] if args.mass_extinction else None

    t0 = time.perf_counter()
    result = simulate_species_tree(
        birth=args.birth, death=args.death, n_extant=args.n_extant, total_time=args.total_time,
        mass_extinctions=mass_ext, sampling=args.sampling, fossils=args.fossils, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "species", args.flat)
    result.write(out, outputs=args.write)

    n_extant = result.n_extant
    n_total = len(result.complete_tree.nodes)
    n_leaves = len([n for n in result.complete_tree.nodes.values() if n.children is None])
    n_extinct = len(result.complete_tree.extinct())
    parts = [f"{n_extant} extant"]
    if n_extinct:
        parts.append(f"{n_extinct} extinct")
    if result.fossils:
        parts.append(f"{len(result.fossils)} fossils")
    summary = " + ".join(parts) + f" ({n_leaves} tips, {n_total} nodes)"
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(level_dir(args.run, "logs", args.flat), "species.log"),
                      args, summary)
    return 0
