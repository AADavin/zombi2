"""``zombi2 traits`` — evolve a trait along a species tree.

``--kind`` picks the state space, which is what genuinely differs between the two trait engines:
``continuous`` (a real value diffusing, :func:`zombi2.traits.simulate_continuous` — Brownian motion,
or Ornstein–Uhlenbeck with ``--reverts-to``/``--pull``) or ``discrete`` (a finite state switching,
:func:`~zombi2.traits.simulate_discrete` — the Mk model with ``--switch``, or the threshold model
with ``--liability``/``--threshold``). Long options are the API keyword names, and ``--rate`` takes
the written form of a rate (SPEC §5) — a bare number, or the same ``scope(base) × modifiers``
expression the Python API takes: ``--rate "1.0 * OnTime({0: 4.0, 1: 1.0})"`` is an early burst,
``--rate "1.0 * FromParent(spread=0.2)"`` variable-rates BM.

Correlated multi-trait runs (``correlation=``) and multi-optimum OU (``regimes=``) need a Python
object, so they stay in the Python API — the CLI covers the single-trait cases.
"""
from __future__ import annotations

import argparse
import os
import time

from zombi2.cli.framework import (_add_flat_arg, _add_params_arg, _rate, _rates_help,
                                  _write_params_log, level_dir, resolve_tree)
from zombi2.cli.genomes import _read_tip_fates
from zombi2.species import read_newick
from zombi2.traits import WIRED_MODIFIERS, simulate_continuous, simulate_discrete

#: the RATES block for ``zombi2 traits -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--rate",
    note="These bend a continuous trait's variance-rate (--rate). The discrete switching rate "
         "(--switch) and the liability rate (--liability) are bare numbers this slice.")

# the write vocabularies, mirroring TraitsResult.write ("driver" is discrete-only: it replays the
# stochastic character map, which a diffusion has no equivalent of)
_CONTINUOUS_OUTPUTS = ("values", "changes", "tree")
_DISCRETE_OUTPUTS = ("values", "changes", "tree", "driver")

# what each kind writes when --write is not given ("driver" stays opt-in: it exists to feed a
# conditioned coupling run, not to describe this one)
_CONTINUOUS_DEFAULT = ("values", "tree")
_DISCRETE_DEFAULT = ("values", "changes", "tree")

# kind-specific knobs — (attribute, default) pairs — rejected under the other kind
_CONTINUOUS_ONLY = (("rate", 1.0), ("reverts_to", None), ("pull", None))
_DISCRETE_ONLY = (("states", None), ("switch", None), ("liability", None), ("threshold", None))


def _add_traits_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("-t", "--tree", required=True, metavar="FILE|DIR",
                   help="the tree the trait rides: a Newick file, or a 'zombi2 species' run "
                        "directory to take its complete tree from (-t out/ instead of -t "
                        "out/species/species_complete.nwk). Any external tree works too; the trait "
                        "evolves on the complete tree, extinct lineages included")
    g.add_argument("-o", "--output", required=True, metavar="DIR", dest="output",
                   help="output directory (created if needed)")
    g.add_argument("--kind", choices=("continuous", "discrete"), default="continuous",
                   metavar="KIND",
                   help="continuous (a real value diffusing, default) or discrete (a finite state "
                        "switching)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("--tip-fates", metavar="FILE", dest="tip_fates",
                   help="[external non-ultrametric trees] a TSV 'tip_name<TAB>extant|extinct' "
                        "declaring each tip's fate; required when the input tree is not ultrametric "
                        "(ZOMBI won't guess extinct lineages from early-sampled tips)")

    g = p.add_argument_group("continuous trait", "only with --kind continuous")
    g.add_argument("--rate", type=_rate, default=1.0, metavar="RATE",
                   help="the Brownian variance-rate σ² — how fast the value diffuses (default 1.0; "
                        "see RATES below)")
    g.add_argument("--reverts-to", type=float, default=None, metavar="VALUE", dest="reverts_to",
                   help="[OU] the optimum the value is pulled toward (needs --pull)")
    g.add_argument("--pull", type=float, default=None, metavar="STRENGTH",
                   help="[OU] how strongly the value is pulled to --reverts-to (needs --reverts-to)")

    g = p.add_argument_group("discrete trait", "only with --kind discrete")
    g.add_argument("--states", metavar="A,B,...", default=None,
                   help="the state space, comma-separated (e.g. marine,terrestrial). Required for "
                        "--kind discrete")
    g.add_argument("--switch", type=float, default=None, metavar="RATE",
                   help="[Mk] the symmetric switching rate between states — a bare number this "
                        "slice (an asymmetric rate matrix needs the Python API)")
    g.add_argument("--liability", type=float, default=None, metavar="RATE",
                   help="[threshold] the variance-rate of the underlying continuous liability — a "
                        "bare number this slice")
    g.add_argument("--threshold", type=float, default=None, metavar="CUT",
                   help="[threshold] the liability value the state flips at")

    g = p.add_argument_group("both kinds")
    g.add_argument("--start", default=None, metavar="VALUE",
                   help="the value at time 0 — a number when --kind continuous (default 0), a state "
                        "label when --kind discrete (default: the first state)")
    g.add_argument("--at-speciation", type=float, default=None, metavar="X",
                   dest="at_speciation",
                   help="add a change at each speciation node: the jump width (Normal(0, X)) when "
                        "--kind continuous, the probability of hopping to another state when "
                        "--kind discrete")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_DISCRETE_OUTPUTS, default=None, metavar="PART",
                   help="which outputs to write (default: values, tree [+ changes when discrete]). "
                        "values: the tip table. changes: the event log. tree: the trait tree "
                        "(annotated Newick). driver: the segment table that feeds a conditioned "
                        "run (discrete only).")
    _add_flat_arg(g)


def run(args, parser):
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

    tip_fates = _read_tip_fates(args.tip_fates) if args.tip_fates else None
    tree_path = resolve_tree(args.tree)
    try:
        with open(tree_path) as f:
            tree, names = read_newick(f.read(), tip_fates=tip_fates)
    except FileNotFoundError:
        raise FileNotFoundError(f"tree file not found: {tree_path}") from None

    t0 = time.perf_counter()
    if discrete:
        result = simulate_discrete(tree, states=states, switch=args.switch, start=args.start,
                                   liability=args.liability, threshold=args.threshold,
                                   at_speciation=args.at_speciation, seed=args.seed)
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
                                     at_speciation=args.at_speciation, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.output, exist_ok=True)
    out = level_dir(args.output, "traits", args.flat)
    outputs = args.write or (_DISCRETE_DEFAULT if discrete else _CONTINUOUS_DEFAULT)
    result.write(out, outputs=outputs)
    if names:  # an external tree: map ZOMBI's n<id> back to the user's labels (join on the node col)
        rows = ["node\tname"] + [f"n{i}\t{lbl}" for i, lbl in sorted(names.items())]
        with open(os.path.join(out, "names.tsv"), "w") as f:
            f.write("\n".join(rows) + "\n")

    n_tips = len(result.values)
    detail = f"{len(states)} states" if discrete else "diffusing"
    summary = f"a {result.kind} trait ({detail}) over {n_tips} extant tips"
    print(f"wrote {args.output}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(level_dir(args.output, "logs", args.flat), "traits.log"),
                      args, summary)
    return 0
