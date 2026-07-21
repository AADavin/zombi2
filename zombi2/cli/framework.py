"""Shared CLI plumbing: the banner/help formatter, ``--params`` handling, the run-parameters log,
and the subcommand builder every command module leans on."""
from __future__ import annotations

import argparse
import datetime
import os
import shutil
import sys

from zombi2 import __version__


_DESCRIPTION = """\
Simulate each level of evolution on its own. Run 'zombi2 <command> -h' for a command's options.

Levels
  species              simulate a dated species tree
  genomes              evolve gene families along a species tree (unordered or ordered)
  sequences            evolve sequences down each gene tree (a prior genomes run)

Traits and the coupled (conditioned / joint) models land here as those levels are rebuilt on the
clean core; until then they live only in the legacy code.
"""


# ── house style: an IQ-TREE-like grouped, sectioned help ────────────────────────────
_BOLD, _RESET = "\033[1m", "\033[0m"


def _use_color() -> bool:
    """Bold section headers only for an interactive terminal (never when piped/redirected, under
    NO_COLOR, or a dumb terminal) — so redirected help stays plain text."""
    if os.environ.get("NO_COLOR") or os.environ.get("TERM") == "dumb":
        return False
    return sys.stdout.isatty()


def _banner() -> str:
    return (f"ZOMBI2 {__version__} — simulating the evolution of species, genomes, "
            "sequences and traits")


class ZombiHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Grouped help in the IQ-TREE house style: argument-group titles become UPPERCASE section
    headers (bold on a terminal), with a wide, aligned help column. The auto usage line is kept
    short by giving each command a hand-written ``usage=``."""

    def __init__(self, prog: str) -> None:
        width = min(shutil.get_terminal_size((90, 24)).columns - 2, 92)
        super().__init__(prog, max_help_position=32, width=width)

    def start_section(self, heading: str | None) -> None:
        if heading and heading not in ("positional arguments", "options", "optional arguments"):
            heading = heading.upper()
            if _use_color():
                heading = _BOLD + heading + _RESET
        super().start_section(heading)

    def _format_action(self, action: argparse.Action) -> str:
        # Hide the auto subcommand list from the top-level help — the commands are curated, grouped
        # by theme, in the description instead (avoids a duplicate, ungrouped dump).
        if isinstance(action, argparse._SubParsersAction):
            return ""
        return super()._format_action(action)


def _examples(*lines: str) -> str:
    """Build an ``EXAMPLES`` epilog block in the house style: a bold header on a TTY (plain when the
    output is piped), followed by the given lines verbatim. Safe because the parser's formatter is
    ``RawDescription``-based, so these line breaks are kept."""
    header = _BOLD + "EXAMPLES" + _RESET if _use_color() else "EXAMPLES"
    return header + "\n" + "\n".join(lines)


def _add_params_arg(g) -> None:
    """Add ``--params FILE`` (a TOML parameters file) to a subcommand's ``general`` group."""
    g.add_argument("--params", metavar="FILE",
                   help="a TOML parameters file whose keys are this command's long option names "
                        "(hyphens or underscores); applied as defaults, so any flag given on the "
                        "command line overrides it. A '[<command>]' table scopes one file to a "
                        "whole pipeline. Required I/O paths (-o / -t) stay on the command line.")


def _write_params_log(path: str, args: argparse.Namespace, summary: str) -> None:
    """Write the full set of run parameters to ``path`` — always, for reproducibility."""
    lines = ["# ZOMBI2 run parameters",
             f"zombi2_version\t{__version__}",
             f"timestamp\t{datetime.datetime.now().isoformat(timespec='seconds')}",
             f"command_line\t{' '.join(sys.argv)}"]
    for key, value in sorted(vars(args).items()):
        lines.append(f"{key}\t{value}")
    lines.append(f"result\t{summary}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _add_subcommand(sub, name: str, help: str, description: str, usage: str, adder,
                    epilog: str | None = None):
    """Register a subcommand with the house-style formatter and a hand-written compact usage.

    The command list itself is curated (grouped by theme) in the top-level description, so the
    per-command ``help`` is suppressed from argparse's auto listing to avoid a duplicate dump.
    ``epilog`` (built with :func:`_examples`) adds a worked-example block below the options.
    """
    p = sub.add_parser(name, help=help, description=description, usage=usage, epilog=epilog,
                       formatter_class=ZombiHelpFormatter)
    adder(p)
    return p


def _apply_params_file(sub, argv) -> None:
    """If the invocation is ``<command> … --params FILE …`` for a params-aware subcommand, load the
    TOML file and set that subcommand's argument defaults from it — so explicit command-line flags,
    parsed afterwards, still override the file."""
    tokens = argv if argv is not None else sys.argv[1:]
    if not tokens or tokens[0].startswith("-"):
        return
    subp = sub.choices.get(tokens[0])
    if subp is None or not any(a.dest == "params" for a in subp._actions):
        return
    path = None
    for i, tok in enumerate(tokens[1:], 1):
        if tok == "--params" and i + 1 < len(tokens):
            path = tokens[i + 1]
            break
        if tok.startswith("--params="):
            path = tok.split("=", 1)[1]
            break
    if path is None:
        return
    from zombi2.cli._params import load_params_file
    action_by_dest = {a.dest: a for a in subp._actions}
    try:
        overrides = load_params_file(path, set(action_by_dest), tokens[0])
    except (OSError, ValueError) as e:              # missing file, TOML error, or unknown key
        subp.error(str(e))
    # a variable-length option (nargs '+'/'*', e.g. --write) is a list on the command line; accept a
    # bare scalar in the file and wrap it, so `write = "events"` works like `--write events`.
    for dest, val in list(overrides.items()):
        action = action_by_dest.get(dest)
        if action is not None and action.nargs in ("+", "*") and not isinstance(val, list):
            overrides[dest] = [val]
    subp.set_defaults(**overrides)
