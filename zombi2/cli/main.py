"""zombi2 command-line entry point — assembles the subcommand parser and dispatches.

Each subcommand lives in its own module (species, genomes, traits, sequences, coevolve, tools,
experimental); this module wires them into one argparse parser via the shared framework and routes
``args.command`` to the module's ``run``. Adding a command is: write a module with an ``_add_*_args``
argument builder and a ``run(args, parser)`` handler, then add one ``_add_subcommand(...)`` call and
one ``_RUN`` entry here.
"""
from __future__ import annotations

import argparse
import sys

from zombi2 import __version__
from zombi2.cli import coevolve, experimental, genomes, sequences, species, tools, traits
from zombi2.cli.framework import (
    _DESCRIPTION, ZombiHelpFormatter, _add_subcommand, _apply_params_file, _banner, _examples,
)

#: command name -> handler; the single source of dispatch (replaces the old if/elif chain)
_RUN = {
    "species": species.run, "genomes": genomes.run, "traits": traits.run,
    "sequences": sequences.run, "coevolve": coevolve.run, "tools": tools.run,
    "experimental": experimental.run,
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="zombi2", description=_banner() + "\n\n" + _DESCRIPTION,
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # 1. a dated species tree (20 extant tips)",
            "  zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o out/",
            "",
            "  # 2. gene families along it",
            "  zombi2 genomes -t out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/",
            "",
            "  # 3. a trait (or sequences, or coupled processes) on the same tree",
            "  zombi2 traits -t out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/",
            "",
            "Run 'zombi2 <command> -h' for a command's options and its own examples.",
        ),
    )
    parser.add_argument("--version", action="version", version=f"ZOMBI2 {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>", required=True)

    _add_subcommand(
        sub, "species", "simulate a dated species tree",
        "Simulate a dated species tree.",
        "zombi2 species -o DIR [--mode MODE] [--diversification PROCESS] [options]",
        species._add_species_args,
        epilog=_examples(
            "  # backward birth-death tree, 50 extant tips at age 5",
            "  zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o out/",
            "",
            "  # forward simulation to a fixed age (complete tree, keeps extinct lineages)",
            "  zombi2 species --mode forward --birth 1 --death 0.4 --age 5 --seed 1 -o out/",
        ))

    _add_subcommand(
        sub, "genomes", "evolve gene families along a species tree",
        "Evolve gene families along a species tree.",
        "zombi2 genomes -t FILE -o DIR [--genome-resolution RESOLUTION] [--rate-per UNIT] "
        "[--write PART ...] [options]",
        genomes._add_rate_args,
        epilog=_examples(
            "  # DTL gene families with a full event log and gene trees",
            "  zombi2 genomes -t out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/",
            "",
            "  # counts-only profiles (scales to very large trees)",
            "  zombi2 genomes -t out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --write profiles --seed 42 -o out/",
        ))

    _add_subcommand(
        sub, "traits", "evolve a phenotypic trait along a given species tree",
        "Evolve a phenotypic trait along a species tree, writing tip and ancestral values.",
        "zombi2 traits -t FILE -o DIR [--model MODEL] [options]", traits._add_trait_args,
        aliases=["trait"],
        epilog=_examples(
            "  # Ornstein-Uhlenbeck continuous trait",
            "  zombi2 traits -t out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/",
            "",
            "  # 3-state discrete Mk trait, 20 replicates",
            "  zombi2 traits -t out/species_tree.nwk --model mk --states 3 --replicates 20 --seed 1 -o out/",
        ))

    _add_subcommand(
        sub, "coevolve", "co-evolve coupled processes (--couple driver:target)",
        "Co-evolve coupled processes over {species, genomes, traits} — pick directed edges with "
        "--couple (e.g. traits:species = SSE, traits:genomes = trait-conditioned gene families).",
        "zombi2 coevolve -o DIR --couple DRIVER:TARGET [-t FILE] [--age T|--tips N] [options]",
        coevolve._add_coevolve_mode_args,
        epilog=_examples(
            "  # trait-conditioned gene families (loss/gain depends on a simulated trait)",
            "  zombi2 coevolve --couple traits:genomes -t out/species_tree.nwk --trait-model mk --states 2 --trait-center --responsive 0.3 --effect-loss 3 --seed 1 -o out/",
            "",
            "  # trait-dependent diversification (BiSSE), grows the tree",
            "  zombi2 coevolve --couple traits:species --sse-model bisse --tips 50 --seed 1 -o out/",
        ))

    _add_subcommand(
        sub, "sequences", "simulate DNA/protein alignments along a genomes run's gene trees",
        "Rescale a 'genomes' run's gene trees from time into substitutions/site under a "
        "gene × lineage clock, then (with --subst-model) simulate a DNA or protein sequence "
        "alignment along each rescaled gene tree.",
        "zombi2 sequences --genomes DIR -o DIR [--subst-model MODEL] "
        "[--clock MODEL [--clock-sigma S]] [options]",
        sequences._add_sequence_args,
        aliases=["sequence"],
        epilog=_examples(
            "  # rescale gene trees into substitutions/site (needs a 'genomes' run done with --write trace)",
            "  zombi2 sequences --genomes out/ --clock autocorrelated-lognormal --clock-sigma 0.4 --family-speed 0.5 --seed 7 -o out/",
            "",
            "  # ...and also simulate DNA alignments under HKY85",
            "  zombi2 sequences --genomes out/ --subst-model hky85 --clock autocorrelated-lognormal --clock-sigma 0.4 --seed 7 -o out/",
        ))

    _add_subcommand(
        sub, "tools", "compute on ZOMBI2 outputs (reconcile, treedist, recon-accuracy, red, parse, export)",
        "Analysis tools that compute on ZOMBI2 outputs — the stable analysis complement to the "
        "simulator (the zombi2.tools layer). Each tool is a sub-subcommand; run "
        "'zombi2 tools <tool> -h' for its options.\n\n"
        "Tools\n"
        "  reconcile            ALE reconciliation likelihood of a gene tree (ALElite)\n"
        "  treedist             tree distances (RF, branch-score, quartet, matching) vs a reference\n"
        "  recon-accuracy       accuracy of an inferred reconciliation vs a known one\n"
        "  red                  Relative Evolutionary Divergence of every node (Parks et al. 2018)\n"
        "  parse                read external ALE / AleRax reconciliation output (reconparser)",
        "zombi2 tools <tool> [options]",
        tools._add_tools_args,
        epilog=_examples(
            "  # ALE reconciliation log-likelihood of a gene tree at given DTL rates",
            "  zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15",
            "",
            "  # Relative Evolutionary Divergence of every node of a tree",
            "  zombi2 tools red -t species_tree.nwk -o out/",
            "",
            "  # summarize an existing ALE / AleRax reconciliation (needs zombi2[reconparser])",
            "  zombi2 tools parse results.ale",
        ))

    _add_subcommand(
        sub, "experimental", "unstable, opt-in models (ils: multispecies-coalescent gene trees)",
        "Experimental, not-yet-validated models (the zombi2.experimental layer) — APIs and outputs "
        "may change. Each is a sub-subcommand; run 'zombi2 experimental <model> -h' for its options.\n\n"
        "Models\n"
        "  ils                  incomplete lineage sorting (multispecies-coalescent gene trees)",
        "zombi2 experimental <model> [options]",
        experimental._add_experimental_args,
        epilog=_examples(
            "  # draw 1000 gene trees under the multispecies coalescent (incomplete lineage sorting)",
            "  zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 -o out/",
        ))

    _apply_params_file(sub, argv)               # --params FILE seeds defaults; CLI flags override
    args = parser.parse_args(argv)              # the banner shows on --help only, not on every run
    # Commands are plural nouns; the singular spellings (trait/sequence) are accepted but deprecated.
    # Normalise so dispatch and the run-manifest filename use the canonical plural.
    _SINGULAR_COMMANDS = {"trait": "traits", "sequence": "sequences"}
    if args.command in _SINGULAR_COMMANDS:
        canonical = _SINGULAR_COMMANDS[args.command]
        print(f"warning: 'zombi2 {args.command}' is deprecated; use 'zombi2 {canonical}'.",
              file=sys.stderr)
        args.command = canonical
    # The sequence-clock --branch-* shortcuts were folded into the --clock interface; warn on use.
    if args.command == "sequences":
        _DEPRECATED_CLOCK_FLAGS = {
            "--branch-speed": "--clock autocorrelated-lognormal --clock-sigma",
            "--branch-bins": "--clock-bins", "--branch-switch-rate": "--clock-switch-rate",
            "--branch-up-bias": "--clock-up-bias",
        }
        for _tok in (argv if argv is not None else sys.argv[1:]):
            _repl = _DEPRECATED_CLOCK_FLAGS.get(_tok.split("=", 1)[0])
            if _repl is not None:
                print(f"warning: {_tok.split('=', 1)[0]} is deprecated; use {_repl}.", file=sys.stderr)
    try:
        return _RUN[args.command](args, parser)
    except (ValueError, RuntimeError, FileNotFoundError, OSError) as e:
        # Report expected failures as a clean one-line error, never a traceback.
        print(f"zombi2: error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
