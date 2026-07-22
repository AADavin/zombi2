"""``zombi2 joint`` — grow a species tree and the level that drives it, together.

**Conditioning** runs one level and then reads it: a trait run writes ``trait_driver.tsv``, and a
later genome run makes its loss rate read that file. It works because the driver can be finished
first. **Joining** is for when it cannot — when the trait drives which lineages speciate, the tree
and the trait each depend on the other at every instant, and neither exists before the other. The two
are grown in one pass (``SPEC §2–4``), which is what this command does.

The driver is named in the rate rather than passed as a file, and that is the whole difference:
``--birth "1.0 * DrivenBy('trait', {'small': 1.0, 'large': 3.0})"`` reads a **live level**, not a
path. Give exactly one driver — a discrete trait (``--states``), which is the BiSSE/MuSSE family, or
gene content (``--origination`` and friends), where a lineage's genome decides how fast it splits.
"""
from __future__ import annotations

import argparse
import os
import time

from zombi2.cli.framework import (_add_flat_arg, _add_params_arg, _add_quiet_arg, _add_run_arg,
                                  _rate, _rates_help, _write_params_log, default_outputs,
                                  level_dir)
from zombi2.cli.traits import _DISCRETE_DEFAULT as TRAITS_DEFAULT
from zombi2.genomes import unordered
from zombi2.joint import simulate_joint
from zombi2.rates.modifiers import DrivenBy
from zombi2.traits import discrete

#: the RATES block for ``zombi2 joint -h``. DrivenBy is the point of the command, so it leads.
RATES_HELP = _rates_help(
    (DrivenBy,), "--birth",
    note="DrivenBy here names a LIVE level, not a file — that is what makes a run joint rather "
         "than conditioned. 'trait' reads the discrete trait grown alongside; 'genomes:count' reads "
         "a lineage's total gene count; 'genomes:<name>' reads whether a named family is present "
         "({'present': 3.0, 'absent': 1.0}). Drive --death too for state-dependent extinction.")

#: the driver flags, by which driver they build — used to reject the other driver's flags rather
#: than ignore them, the discipline every other command follows
_TRAIT_ONLY = (("states", None), ("switch", None), ("trait_start", None))
_GENOME_ONLY = (("duplication", 0.0), ("loss", 0.0), ("origination", 0.0),
                ("initial_families", 0), ("families", None))


def _add_joint_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "where this run's files are written — both levels land here")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")

    # validated in run() rather than argparse-`required`, so a --params file can supply them
    g = p.add_argument_group("diversification", "the per-lineage rates — drive one with DrivenBy "
                                                "(see RATES below)")
    g.add_argument("--birth", type=_rate, default=None, metavar="RATE",
                   help="speciation rate (per lineage) — required")
    g.add_argument("--death", type=_rate, default=0.0, metavar="RATE",
                   help="extinction rate (per lineage); 0 = pure birth (default 0)")

    g = p.add_argument_group("stop condition", "grow until exactly one of these — required")
    g.add_argument("--n-extant", type=int, default=None, metavar="N", dest="n_extant",
                   help="stop at N extant lineages — conditioned on survival")
    g.add_argument("--total-time", type=float, default=None, metavar="T", dest="total_time",
                   help="grow forward for T time units")

    g = p.add_argument_group("driver: a discrete trait", "state-dependent diversification "
                                                         "(BiSSE / MuSSE)")
    g.add_argument("--states", metavar="A,B,...", default=None,
                   help="the trait's state space, comma-separated (e.g. small,large). Giving this "
                        "makes the trait the driver")
    g.add_argument("--switch", type=float, default=None, metavar="RATE",
                   help="the symmetric switching rate between states")
    g.add_argument("--trait-start", metavar="STATE", default=None, dest="trait_start",
                   help="the root state (default: one of --states drawn uniformly)")

    g = p.add_argument_group("driver: gene content", "a lineage's genome decides how fast it splits")
    g.add_argument("--duplication", type=_rate, default=0.0, metavar="RATE",
                   help="gene duplication rate (per copy)")
    g.add_argument("--loss", type=_rate, default=0.0, metavar="RATE",
                   help="gene loss rate (per copy)")
    g.add_argument("--origination", type=_rate, default=0.0, metavar="RATE",
                   help="new-family origination rate (per lineage)")
    g.add_argument("--initial-families", type=int, default=0, metavar="N", dest="initial_families",
                   help="gene families the root genome starts with (default 0)")
    g.add_argument("--families", metavar="A,B,...", default=None,
                   help="named families to seed, comma-separated — a name is what "
                        "DrivenBy('genomes:<name>', …) reads")

    g = p.add_argument_group("outputs")
    _add_flat_arg(g)
    _add_quiet_arg(g)


def _stray(args, spec) -> list[str]:
    return [f"--{attr.replace('_', '-')}" for attr, default in spec if getattr(args, attr) != default]


def run(args, parser):
    if args.birth is None:
        parser.error("--birth is required (give it on the command line or in --params)")
    if (args.n_extant is None) == (args.total_time is None):
        parser.error("give exactly one stop condition: --n-extant N or --total-time T")

    # exactly one driver, and no flags from the other one — a silently-ignored driver flag would
    # give a run that looks joint and is not
    genome_flags = _stray(args, _GENOME_ONLY)
    if args.states and genome_flags:
        is_are = "does" if len(genome_flags) == 1 else "do"
        parser.error(f"give one driver: --states builds a trait driver, so "
                     f"{', '.join(genome_flags)} {is_are} not apply (drop them, or drop --states "
                     "for a gene-content driver)")
    if not args.states and not genome_flags:
        parser.error("a joint run needs a driver: --states for a discrete trait "
                     "(state-dependent diversification), or --origination / --duplication / --loss "
                     "for gene content")
    if not args.states and (trait_flags := _stray(args, _TRAIT_ONLY)):
        parser.error(f"these options need --states: {', '.join(trait_flags)}")

    if args.states:
        states = [s.strip() for s in args.states.split(",") if s.strip()]
        if len(states) < 2:
            parser.error(f"--states needs at least two, got {args.states!r}")
        driver = dict(trait=discrete(states=states, switch=args.switch, start=args.trait_start))
    else:
        names = [s.strip() for s in args.families.split(",") if s.strip()] if args.families else None
        driver = dict(genome=unordered(duplication=args.duplication, loss=args.loss,
                                       origination=args.origination,
                                       initial_families=args.initial_families, families=names))

    t0 = time.perf_counter()
    result = simulate_joint(birth=args.birth, death=args.death, n_extant=args.n_extant,
                            total_time=args.total_time, seed=args.seed, **driver)
    dt = time.perf_counter() - t0

    os.makedirs(args.run, exist_ok=True)
    # Both levels belong to one run, so each is written where — and how — its own command would
    # write it: same directory, same default outputs, gene trees in their own subdirectory. Reaching
    # for Result.write's bare default here would quietly give a joint run fewer files than the two
    # commands it stands in for.
    result.species.write(level_dir(args.run, "species", args.flat))
    if result.trait is not None:
        result.trait.write(level_dir(args.run, "traits", args.flat), outputs=TRAITS_DEFAULT)
        detail = "a discrete trait driving speciation"
    else:
        out = level_dir(args.run, "genomes", args.flat)
        wanted = default_outputs(result.genome)
        result.genome.write(out, outputs=[o for o in wanted if o != "gene_trees"])
        result.genome.write(level_dir(out, "gene_trees", args.flat), outputs=("gene_trees",))
        detail = "gene content driving speciation"

    n_extant = len(result.species.complete_tree.extant())
    summary = f"{n_extant} extant tips, {detail}"
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(level_dir(args.run, "species", args.flat), "joint.log"),
                      args, summary)
    return 0
