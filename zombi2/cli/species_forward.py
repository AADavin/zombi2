"""``zombi2 species`` for the new forward engine — the replacement-in-waiting for ``cli/species.py``.

Mirrors :func:`zombi2.species_tree.simulate_species_tree`: ``birth``/``death`` rates with the common
modifiers (skyline, diversity-dependent, clade drift), a stop of ``--n-extant`` survivors or a
``--total-time``, plus mass extinctions, incomplete sampling, and fossils. Rich modifier
compositions and death-rate modifiers stay Python-only — the CLI exposes the common cases.

**Standalone for now, not wired into the ``species`` command.** The old ``cli/species.py`` is the
pipeline entry point for ~two dozen downstream tests (they run ``species`` then feed its tree to
``genomes``/``traits``), so this flips on together with the pipeline's move to the new engine. It is
exercised by its own tests, which build the parser directly.
"""
from __future__ import annotations

import argparse

from zombi2 import modifiers as mod
from zombi2 import scope
from zombi2.cli.framework import _add_params_arg, _write_params_log
from zombi2.species_tree import simulate_species_tree


def _skyline(text: str) -> dict:
    """Parse a ``--skyline`` value ``0:1.0,3:0.3`` into ``{0.0: 1.0, 3.0: 0.3}`` for a Time modifier."""
    schedule: dict[float, float] = {}
    for step in text.split(","):
        t, sep, f = step.partition(":")
        if not sep:
            raise argparse.ArgumentTypeError(f"--skyline step {step!r} must be time:factor")
        schedule[float(t)] = float(f)
    return schedule


def _add_species_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("rates")
    g.add_argument("--birth", type=float, default=1.0, metavar="RATE",
                   help="per-lineage speciation rate (default 1.0)")
    g.add_argument("--death", type=float, default=0.0, metavar="RATE",
                   help="per-lineage extinction rate (default 0.0 = Yule)")
    g.add_argument("--global-birth", action="store_true",
                   help="one shared birth budget for the whole tree (linear growth), not per lineage")

    g = p.add_argument_group("how the birth rate varies", "the common cases; compose freely")
    g.add_argument("--skyline", type=_skyline, default=None, metavar="T:F,...",
                   help="episodic birth: a relative factor per interval, e.g. 0:1.0,3:0.3")
    g.add_argument("--diversity-cap", type=float, default=None, metavar="K",
                   help="diversity-dependent: birth slows to 0 as diversity approaches K")
    g.add_argument("--clade-drift", type=float, default=None, metavar="SIGMA",
                   help="birth rate drifts down the tree with this spread (inherited at each split)")

    g = p.add_argument_group("when to stop", "give exactly one")
    g.add_argument("--n-extant", type=int, default=None, metavar="N",
                   help="stop at N surviving lineages (conditioned on survival)")
    g.add_argument("--total-time", type=float, default=None, metavar="T",
                   help="grow to time T (measured forward from the crown)")

    g = p.add_argument_group("interventions & observation")
    g.add_argument("--mass-extinction", action="append", nargs=2, type=float, default=None,
                   metavar=("TIME", "FRACTION"),
                   help="at TIME, cull FRACTION of the living (repeatable; needs --total-time)")
    g.add_argument("--sampling", type=float, default=1.0, metavar="RHO",
                   help="observe only fraction RHO of the survivors (default 1.0)")
    g.add_argument("--fossils", type=float, default=0.0, metavar="PSI",
                   help="recover fossils along branches at rate PSI (a side output)")

    g = p.add_argument_group("output")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")
    g.add_argument("--write", nargs="+", default=None,
                   choices=("complete", "extant", "events", "fossils"), metavar="WHAT",
                   help="which outputs to write (default: all applicable)")
    g.add_argument("--seed", type=int, default=None, metavar="N", help="random seed")
    _add_params_arg(g)


def _birth_rate(args: argparse.Namespace):
    """Assemble the birth rate spec — ``scope(base) × modifiers`` from the requested flags."""
    birth = scope.Global(args.birth) if args.global_birth else args.birth
    if args.skyline is not None:
        birth = birth * mod.Time(args.skyline)
    if args.diversity_cap is not None:
        birth = birth * mod.Diversity(cap=args.diversity_cap)
    if args.clade_drift is not None:
        birth = birth * mod.Inherited(spread=args.clade_drift)
    return birth


def run(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if (args.n_extant is None) == (args.total_time is None):
        parser.error("give exactly one of --n-extant or --total-time")
    mass_ext = [(t, f) for t, f in args.mass_extinction] if args.mass_extinction else None
    result = simulate_species_tree(
        birth=_birth_rate(args), death=args.death,
        n_extant=args.n_extant, total_time=args.total_time,
        mass_extinctions=mass_ext, sampling=args.sampling, fossils=args.fossils, seed=args.seed)
    result.write(args.out, outputs=args.write)

    parts = [f"{result.n_extant} extant"]
    for count, label in ((len(result.complete_tree.extinct()), "extinct"),
                         (len(result.complete_tree.unsampled()), "unsampled"),
                         (len(result.fossils), "fossils")):
        if count:
            parts.append(f"{count} {label}")
    summary = " + ".join(parts)
    print(f"wrote species tree to {args.out}/ ({summary})")
    _write_params_log(f"{args.out}/species.log", args, summary)
    return 0
