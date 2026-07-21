"""``zombi2`` command-line entry point — assembles the subcommand parser and dispatches.

Each subcommand lives in its own module and mirrors one level's ``simulate_*`` function; this
module wires them into one argparse parser via the shared framework and routes ``args.command`` to
the module's ``run``. Adding a command is: write a module with an ``_add_*_args`` argument builder
and a ``run(args, parser)`` handler, then add one ``_add_subcommand(...)`` call and one ``_RUN``
entry here.
"""
from __future__ import annotations

import argparse
import sys

from zombi2 import __version__
from zombi2.cli import genomes, sequences, species, traits
from zombi2.cli.framework import (
    _DESCRIPTION, ZombiHelpFormatter, _add_subcommand, _apply_params_file, _banner, _examples,
)

#: command name -> handler; the single source of dispatch
_RUN = {"species": species.run, "genomes": genomes.run, "sequences": sequences.run,
        "traits": traits.run}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="zombi2", description=_banner() + "\n\n" + _DESCRIPTION,
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # 1. a dated species tree (20 extant tips)",
            "  zombi2 species --birth 1 --death 0.3 --n-extant 20 --seed 1 -o out/",
            "",
            "  # 2. gene families along it",
            "  zombi2 genomes -t out/species_complete.nwk --duplication 0.2 --transfer 0.1 "
            "--loss 0.25 --origination 0.5 --seed 42 -o out/",
            "",
            "  # 3. sequences down each gene tree",
            "  zombi2 sequences --genomes out/ --model hky85 --length 1000 --seed 1 -o out/",
            "",
            "Run 'zombi2 <command> -h' for a command's options and its own examples.",
        ),
    )
    parser.add_argument("--version", action="version", version=f"ZOMBI2 {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)

    _add_subcommand(
        sub, "species", "simulate a dated species tree",
        "Simulate a dated species tree by a per-lineage birth–death process (time runs forward "
        "from the crown).",
        "zombi2 species -o DIR --birth RATE (--n-extant N | --total-time T) [options]",
        species._add_species_args,
        epilog=_examples(
            "  # 20 extant tips, birth–death",
            "  zombi2 species --birth 1 --death 0.3 --n-extant 20 --seed 1 -o out/",
            "",
            "  # grow for a fixed time, with a mass-extinction pulse at t=3",
            "  zombi2 species --birth 1 --death 0.4 --total-time 5 --mass-extinction 3 0.75 "
            "--seed 1 -o out/",
            "",
            "  # a skyline: speciation drops to a third at time 3 (see RATES)",
            "  zombi2 species --birth \"1.0 * OnTime({0: 1.0, 3: 0.3})\" --death 0.3 "
            "--total-time 5 --seed 1 -o out/",
        ) + "\n\n" + species.RATES_HELP)

    _add_subcommand(
        sub, "genomes", "evolve gene families along a species tree",
        "Evolve gene families along a species tree, at the unordered (gene-family counts) or "
        "ordered (genes positioned on chromosomes) resolution.",
        "zombi2 genomes -t FILE -o DIR [--resolution RESOLUTION] [options]",
        genomes._add_genomes_args,
        epilog=_examples(
            "  # unordered D/T/L/O gene families, with the event log and profiles",
            "  zombi2 genomes -t out/species_complete.nwk --duplication 0.2 --transfer 0.1 "
            "--loss 0.25 --origination 0.5 --seed 42 -o out/",
            "",
            "  # ordered genomes with inversions on 3 chromosomes",
            "  zombi2 genomes -t out/species_complete.nwk --resolution ordered --duplication 0.2 "
            "--loss 0.2 --origination 0.5 --inversion 0.3 --chromosomes 3 --seed 42 -o out/",
            "",
            "  # loss twice as fast from time 2 onward (see RATES)",
            "  zombi2 genomes -t out/species_complete.nwk --duplication 0.2 "
            "--loss \"0.25 * OnTime({0: 1.0, 2: 2.0})\" --origination 0.5 --seed 42 -o out/",
        ) + "\n\n" + genomes.RATES_HELP)

    _add_subcommand(
        sub, "sequences", "evolve sequences down each gene tree",
        "Evolve one sequence inside each gene, down its gene tree, under a nucleotide substitution "
        "model and a per-site rate. Replays a prior 'zombi2 genomes' run (--genomes DIR).",
        "zombi2 sequences --genomes DIR -o DIR --model MODEL [options]",
        sequences._add_sequence_args,
        epilog=_examples(
            "  # HKY85, 1000 sites, strict clock, along a prior genomes run",
            "  zombi2 sequences --genomes out/ --model hky85 --kappa 2 --length 1000 --seed 1 "
            "-o seqs/",
            "",
            "  # GTR with an uncorrelated (relaxed) lineage clock",
            "  zombi2 sequences --genomes out/ --model gtr --frequencies 0.3 0.2 0.2 0.3 "
            "--substitution \"1.0 * ByLineage(spread=0.3)\" --seed 1 -o seqs/",
        ) + "\n\n" + sequences.RATES_HELP)

    _add_subcommand(
        sub, "traits", "evolve a trait along a species tree",
        "Evolve a trait along a species tree, with a continuous (a real value diffusing) or "
        "discrete (a finite state switching) state space.",
        "zombi2 traits -t FILE -o DIR [--kind KIND] [options]",
        traits._add_traits_args,
        epilog=_examples(
            "  # a continuous trait diffusing by Brownian motion (variance-rate 1.0)",
            "  zombi2 traits -t out/species_complete.nwk --rate 1.0 --seed 1 -o out/",
            "",
            "  # the same value pulled toward an optimum (Ornstein-Uhlenbeck)",
            "  zombi2 traits -t out/species_complete.nwk --rate 1.0 --reverts-to 2 --pull 0.5 "
            "--seed 1 -o out/",
            "",
            "  # a discrete habitat flipping between two states (Mk)",
            "  zombi2 traits -t out/species_complete.nwk --kind discrete "
            "--states marine,terrestrial --switch 0.1 --seed 1 -o out/",
            "",
            "  # an early burst: the variance-rate starts at 4 and settles to 1 (see RATES)",
            "  zombi2 traits -t out/species_complete.nwk "
            "--rate \"1.0 * OnTime({0: 4.0, 1: 1.0})\" --seed 1 -o out/",
        ) + "\n\n" + traits.RATES_HELP)

    _apply_params_file(sub, argv)               # --params FILE seeds defaults; CLI flags override
    args = parser.parse_args(argv)              # the banner shows on --help only, not on every run
    try:
        return _RUN[args.command](args, parser)
    except (ValueError, RuntimeError, FileNotFoundError, OSError) as e:
        # Report expected failures as a clean one-line error, never a traceback.
        print(f"zombi2: error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
