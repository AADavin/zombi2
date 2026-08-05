"""``zombi2 traits`` — evolve a trait along a species tree.

``--kind`` picks the state space, which is what genuinely differs between the two trait engines:
``continuous`` (a real value diffusing, `zombi2.traits.simulate_continuous()` — Brownian motion,
or Ornstein–Uhlenbeck with ``--reverts-to``/``--pull``) or ``discrete`` (a finite state switching,
`simulate_discrete()` — the Mk model with ``--switch``, or the threshold model
with ``--liability``/``--threshold``). Long options are the API keyword names, and every rate flag —
``--rate``, ``--switch``, ``--liability`` — takes the written form of a rate (SPEC §5): a bare number,
or the same ``scope(base) × modifiers`` expression Python takes. ``--rate "1.0 * OnTime({0: 4.0, 1:
1.0})"`` is an early burst, ``--rate "1.0 * FromParent(spread=0.2)"`` variable-rates BM, and
``--switch "0.2 * DrivenBy('habitat/trait_events.tsv', {'aquatic': 3.0, 'terrestrial': 1.0})"`` a
discrete trait conditioned on one grown first. ``--switch`` also reads its keyword's other two shapes
— a ``{'a->b': rate}`` dict, a ``k x k`` matrix — whose entries are rates in that same form.

They do not take the same *modifiers*: each engine declares what it reads and refuses the rest by
name, so one a switch rate cannot honour raises there rather than being dropped here.

Correlated multi-trait runs (``correlation=``) and multi-optimum OU (``regimes=``) need a Python
object, so they stay in the Python API — the CLI covers the single-trait cases.
"""
from __future__ import annotations

import argparse
import os
import time

from zombi2.cli.framework import (resolve_seed, _add_flat_arg, _add_force_arg, _add_quiet_arg, _add_from_arg,
                                  _add_params_arg, _add_run_arg, _rate, _rates_help, _read_tip_fates,
                                  _write_params_log, check_stale_downstream, clear_stale_downstream,
                                  signpost, input_digests, level_dir, resolve_tree, sibling_fates,
                                  conditioned_levels, record_conditioning)
from zombi2.tree import node_label, read_newick
from zombi2._runtime.report import write_run_report
from zombi2.traits import IMPLEMENTED_MODIFIERS, simulate_continuous, simulate_discrete
from zombi2.traits.discrete import IMPLEMENTED_MODIFIERS as _SWITCH_MODIFIERS

#: the RATES block for ``zombi2 traits -h``, built from the level's own declarations — the listed
#: modifiers are the continuous ``--rate``'s, and the note names the discrete engine's shorter list
#: from its own tuple rather than by hand, so neither can fall behind what the engine takes.
RATES_HELP = _rates_help(
    IMPLEMENTED_MODIFIERS, "--rate",
    note="--rate, --switch and --liability all take a rate in this form. The list above is --rate's; "
         "--switch takes " + ", ".join(m.__name__ for m in _SWITCH_MODIFIERS) + " and nothing else, "
         "and --liability no modifier yet. --switch also reads a {'a->b': rate} dict — only the named "
         "transitions happen — or a k x k matrix of them.")

# the write vocabularies, mirroring TraitsResult.write. The event log IS the driver file now
# (a driven run replays it against the tree), so there is no separate driver output.
_CONTINUOUS_OUTPUTS = ("values", "events", "tree", "summary")
_DISCRETE_OUTPUTS = ("values", "events", "tree", "summary")

# what each kind writes when --write is not given
_CONTINUOUS_DEFAULT = ("values", "tree", "summary")
_DISCRETE_DEFAULT = ("values", "events", "tree", "summary")

# kind-specific knobs — (attribute, default) pairs — rejected under the other kind
_CONTINUOUS_ONLY = (("rate", 1.0), ("reverts_to", None), ("pull", None))
_DISCRETE_ONLY = (("states", None), ("switch", None), ("liability", None), ("threshold", None))


def _add_traits_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "the trait rides the species tree it already holds")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    _add_from_arg(g, "the tree the trait rides", "a Newick file, or another run's directory")
    # validated in run() rather than argparse-`required`, so a --params file can supply it
    g.add_argument("--kind", choices=("continuous", "discrete"), default=None, metavar="KIND",
                   help="the state space — continuous or discrete (required)")
    g.add_argument("--name", default=None, metavar="NAME",
                   help="write this trait to traits/NAME/, so a run can hold several — and one can "
                        "drive another")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("--tip-fates", metavar="FILE", dest="tip_fates",
                   help="a TSV 'tip_name<TAB>extant|extinct|unsampled'; required when the input "
                        "tree is not ultrametric")

    g = p.add_argument_group("continuous trait", "only with --kind continuous")
    g.add_argument("--rate", type=_rate, default=1.0, metavar="RATE",
                   help="the Brownian variance-rate σ² (default 1.0; see RATES below)")
    g.add_argument("--reverts-to", type=float, default=None, metavar="VALUE", dest="reverts_to",
                   help="[OU] the optimum the value is pulled toward (needs --pull)")
    g.add_argument("--pull", type=float, default=None, metavar="STRENGTH",
                   help="[OU] the strength of that pull (needs --reverts-to)")

    g = p.add_argument_group("discrete trait", "only with --kind discrete")
    g.add_argument("--states", metavar="A,B,...", default=None,
                   help="the state space (required), e.g. marine,terrestrial")
    g.add_argument("--switch", type=_rate, default=None, metavar="RATE",
                   help="[Mk] the switching rate between states: one symmetric rate, a "
                        "{'a->b': rate} dict of the transitions to allow, or a k x k matrix "
                        "(see RATES below)")
    g.add_argument("--liability", type=_rate, default=None, metavar="RATE",
                   help="[threshold] the variance-rate of the underlying liability")
    g.add_argument("--threshold", type=float, default=None, metavar="CUT",
                   help="[threshold] the liability value the state flips at")

    g = p.add_argument_group("both kinds")
    g.add_argument("--start", default=None, metavar="VALUE",
                   help="the value at time 0 — a number when continuous (default 0), a state "
                        "label when discrete (default: uniform over --states)")
    g.add_argument("--at-speciation", type=float, default=None, metavar="X",
                   dest="at_speciation",
                   help="a change at each speciation node — jump VARIANCE: Normal(0, X), so the "
                        "width is sqrt(X), when "
                        "continuous, hop probability when discrete")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_DISCRETE_OUTPUTS, default=None, metavar="PART",
                   help="which outputs to write: values (every node), events (the root state "
                        "then every switch — what a conditioned run reads), tree (annotated "
                        "Newick), summary. Default: all but events, and events too when discrete")
    _add_flat_arg(g)
    _add_quiet_arg(g)
    _add_force_arg(g)


def _traits_slot(args, parser) -> str:
    """Which directory under the run this trait belongs in: ``traits`` or ``traits/<name>``.

    A run directory holds one slot per level, which is right for the levels a run has one of. It is
    wrong for traits: a tree can carry several, and one of them can drive another (Chapter 9), so
    they need somewhere to sit side by side. ``--name`` gives each its own.

    An unnamed run keeps the plain ``traits/`` slot and overwrites what is there, which is what
    re-running a level means — the slot is the level's, and a second run of it replaces the first."""
    if args.name is None:
        return "traits"
    name = args.name.strip()
    if not name or os.sep in name or (os.altsep and os.altsep in name) or name in (".", ".."):
        parser.error(f"--name must be a plain directory name, got {args.name!r}")
    if args.flat:
        parser.error("--name and --flat ask for opposite things: --flat writes every file into the "
                     "run directory, which is what leaves two traits nowhere to go.")
    return os.path.join("traits", name)


def run(args, parser):
    if args.kind is None:
        parser.error("--kind is required: continuous (a real value diffusing) or discrete (a "
                     "finite state switching)")
    discrete = args.kind == "discrete"

    # reject the other kind's knobs, so a silently-ignored flag can't give a misleading run
    stray_spec = _CONTINUOUS_ONLY if discrete else _DISCRETE_ONLY
    stray = [f"--{attr.replace('_', '-')}" for attr, default in stray_spec
             if getattr(args, attr) != default]
    if stray:
        other = "continuous" if discrete else "discrete"
        parser.error(f"these options need --kind {other}: {', '.join(stray)}")

    vocab = _DISCRETE_OUTPUTS if discrete else _CONTINUOUS_OUTPUTS
    if args.write:
        bad = [o for o in args.write if o not in vocab]
        if bad:
            parser.error(f"--write {' '.join(bad)} not available for --kind {args.kind}; "
                         f"choose from: {', '.join(vocab)}")

    states = None
    if discrete:
        if not args.states:
            parser.error("--kind discrete needs --states (e.g. --states marine,terrestrial)")
        states = [s.strip() for s in args.states.split(",") if s.strip()]
        if len(states) < 2:
            parser.error(f"--states needs at least two states, got {args.states!r}")
        if args.switch is None and args.liability is None and args.threshold is None:
            parser.error("--kind discrete needs --switch (the Mk model) or "
                         "--liability/--threshold (the threshold model)")
        # the rate grammar can also spell a dict or a list, and neither means anything here: a dict is
        # the keyword's several-liabilities form, which needs the correlation they share and has no
        # flag for it; a list reached `as_rate` as a TypeError, which the CLI does not catch.
        if isinstance(args.liability, (dict, list, tuple)):
            parser.error(
                "--liability is one variance-rate. Several liabilities — correlated discrete traits "
                "grown together — stay in the Python API, simulate_discrete(liability={...}, "
                "correlation={...}).")

    # refuse up front if re-running would orphan a level conditioned on this trait (unless --force)
    resolve_seed(args)                      # a run must be reproducible from its own log
    check_stale_downstream(args, "traits")

    tree_path = resolve_tree(args.source or args.run, is_run_dir=args.source is None)
    # an explicit --tip-fates wins; otherwise pick up the run's own species_fates.tsv so extinct and
    # unsampled tips are read from the record rather than guessed from tip depth
    tip_fates = _read_tip_fates(args.tip_fates) if args.tip_fates else sibling_fates(tree_path)
    try:
        with open(tree_path, encoding="utf-8") as f:
            tree, names = read_newick(f.read(), tip_fates=tip_fates)
    except FileNotFoundError:
        raise FileNotFoundError(f"tree file not found: {tree_path}") from None

    t0 = time.perf_counter()
    if discrete:
        result = simulate_discrete(tree, states=states, switch=args.switch, start=args.start,
                                   liability=args.liability, threshold=args.threshold,
                                   at_speciation=args.at_speciation, seed=args.seed,
                                   progress=not args.quiet)
    else:
        if args.start is None:
            start = 0.0
        else:
            try:
                start = float(args.start)
            except ValueError:
                parser.error(f"--start must be a number when --kind continuous, got {args.start!r}")
        result = simulate_continuous(tree, start=start, rate=args.rate,
                                     reverts_to=args.reverts_to, pull=args.pull,
                                     at_speciation=args.at_speciation, seed=args.seed,
                                     progress=not args.quiet)
    dt = time.perf_counter() - t0

    clear_stale_downstream(args, "traits")   # --force: drop the now-stale downstream (run succeeded)
    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, _traits_slot(args, parser), args.flat)
    outputs = args.write or (_DISCRETE_DEFAULT if discrete else _CONTINUOUS_DEFAULT)
    result.write(out, outputs=outputs)
    if names:  # an external tree: map ZOMBI's n<id> back to the user's labels (join on the node col)
        rows = ["node\tname"] + [f"{node_label(i)}\t{lbl}" for i, lbl in sorted(names.items())]
        with open(os.path.join(out, "names.tsv"), "w", encoding="utf-8") as f:
            f.write("\n".join(rows) + "\n")

    if not args.flat:                             # record which same-run levels this run reads (if any),
        record_conditioning(out, conditioned_levels(  # so re-running one of them knows it orphans this
            args.run, (args.rate, args.switch)))
    n_tips = len(result.values)
    detail = f"{len(states)} states" if discrete else "diffusing"
    summary = f"a {result.kind} trait ({detail}) over {n_tips} extant tips"
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    # the log is the run's parameters, not the parser's: a discrete run has no --rate and a continuous
    # one no --switch, so drop the other kind's knobs (the same set rejected as stray above) — recording
    # them at their defaults would read as though the run had them and chose those values.
    _write_params_log(os.path.join(out, "traits.log"), args, summary,
                      omit={attr for attr, _default in stray_spec},
                      inputs=input_digests(tree_path, args.tip_fates,
                                           os.path.join(os.path.dirname(tree_path),
                                                        "species_fates.tsv"),
                                           args.rate, args.switch))
    signpost(args, write_run_report(args.run), out)   # every file it wrote, then the run report
    return 0
