"""Shared CLI plumbing: the banner/help formatter, ``--params`` handling, the run-parameters log,
and the subcommand builder every command module leans on."""
from __future__ import annotations

import argparse
import datetime
import os
import shutil
import sys
import textwrap

from zombi2 import __version__


_DESCRIPTION = """\
Simulate each level of evolution on its own. Run 'zombi2 <command> -h' for a command's options.

Levels
  species              simulate a dated species tree
  genomes              evolve gene families along a species tree (unordered or ordered)
  sequences            evolve sequences down each gene tree (a prior genomes run)
  traits               evolve a trait along a species tree (continuous or discrete)

The coupled (conditioned / joint) models land here as those levels are rebuilt on the clean core;
until then they live only in the legacy code.
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


def _rate(text: str):
    """The argparse ``type`` for every rate flag: the written form of a rate (SPEC §5).

    ``--birth 1.0`` and ``--birth "1.0 * OnTime({0: 1.0, 3: 0.3})"`` both come through here, so the
    command line takes exactly the expression the Python API takes. Re-raised as an
    ``ArgumentTypeError`` so argparse prints the parser's own message ("unknown name 'OnDiversity'
    — did you mean …?") instead of burying it under a generic "invalid value".
    """
    from zombi2.rates.parse import parse_rate

    try:
        return parse_rate(text)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e)) from None


#: one gloss and one worked snippet per modifier, for the RATES help block. A modifier with no entry
#: still lists (by name), so the help can never fall behind a level's ``WIRED_MODIFIERS`` declaration.
_MODIFIER_HELP = {
    "OnTime": ("OnTime({0: 1.0, 3: 0.3})", "the rate changes in time — a skyline"),
    "OnTotalDiversity": ("OnTotalDiversity(cap=100)", "the rate slows as the clade fills up"),
    "FromParent": ("FromParent(spread=0.2)", "the rate drifts down the tree (ClaDS)"),
    "ByLineage": ("ByLineage(spread=0.3)", "one draw per lineage — the uncorrelated clock"),
    "DrivenBy": ("DrivenBy('habitat.tsv', {'aquatic': 3.0})", "the rate is driven by another level"),
}


def _wrap_note(note: str, width: int = 86) -> list[str]:
    return textwrap.wrap(note, width=width, initial_indent="  ", subsequent_indent="  ")


def _rates_help(wired, flag: str, *, scopes: str | None = None, note: str | None = None) -> str:
    """The ``RATES`` epilog block for a command, built from that level's ``WIRED_MODIFIERS``.

    Listing what the engine *declares* (rather than a hand-kept list) is what keeps the help honest:
    a modifier the level does not wire is rejected by the engine, so it must not be advertised here
    either — and the worked example is drawn from the same list, so it is always a modifier that runs.
    """
    examples = [_MODIFIER_HELP[m.__name__][0] for m in wired if m.__name__ in _MODIFIER_HELP]
    shown = examples[0] if examples else None
    key = flag.lstrip("-").replace("-", "_")

    header = _BOLD + "RATES" + _RESET if _use_color() else "RATES"
    lines = [header,
             "  Every rate is scope(base) × modifiers (SPEC §5) — a bare number, or the same",
             "  expression you would write in Python, quoted:",
             "",
             f"    {flag} 1.0"]
    if shown:
        lines.append(f'    {flag} "1.0 * {shown}"')
    if scopes:
        lines.append(f'    {flag} "{scopes}"')
    lines += ["", "  Modifiers wired for this level (anything else is an error, never ignored):"]
    for m in wired:
        name = m.__name__
        entry = _MODIFIER_HELP.get(name)
        lines.append(f"    {name:<20}{entry[1]}" if entry else f"    {name}")
    lines.append("")
    if note:
        lines += _wrap_note(note)
        lines.append("")
    if shown:
        lines.append(f'  A --params file takes the same text:  {key} = "1.0 * {shown}"')
    return "\n".join(lines)


def _add_params_arg(g) -> None:
    """Add ``--params FILE`` (a TOML parameters file) to a subcommand's ``general`` group."""
    g.add_argument("--params", metavar="FILE",
                   help="a TOML parameters file whose keys are this command's long option names "
                        "(hyphens or underscores); applied as defaults, so any flag given on the "
                        "command line overrides it. A '[<command>]' table scopes one file to a "
                        "whole pipeline. The run directory stays on the command line.")


def _add_flat_arg(g) -> None:
    """Add ``--flat`` (write everything into one directory) to a subcommand's ``outputs`` group."""
    g.add_argument("--flat", action="store_true",
                   help="write every file straight into the output directory instead of grouping "
                        "them by level. A run of a hundred families writes hundreds of files, so "
                        "the grouped layout is the default; use this when another tool expects one "
                        "flat directory")


def default_outputs(result) -> tuple[str, ...]:
    """What a result writes when ``--write`` is not given — read off its own ``write()`` signature.

    The CLI has to know the default before it calls ``write``, because it routes the one-file-per-family
    outputs into their own directories. Reading it here rather than repeating it keeps the two from
    drifting: a level that changes what it writes by default changes it in one place."""
    import inspect

    return tuple(inspect.signature(result.write).parameters["outputs"].default)


def _add_quiet_arg(g) -> None:
    """Add ``--quiet`` — no progress bar, for a log file or a batch of runs."""
    g.add_argument("--quiet", action="store_true",
                   help="no progress bar. A command shows one while it works, which is noise in a "
                        "log file or a script running hundreds of replicates")


def level_dir(output: str, level: str, flat: bool) -> str:
    """Where one level's files belong: ``<output>/<level>/``, or ``<output>/`` under ``--flat``.

    Grouping is the CLI's business, not the engines': a ``Result.write`` writes whatever it is given
    into the one directory it is handed, and the layout is chosen here by calling it more than once.
    Created on the way out, so a caller can write into it immediately."""
    path = output if flat else os.path.join(output, level)
    os.makedirs(path, exist_ok=True)
    return path


def _add_run_arg(p, what: str) -> None:
    """Add the run directory — the one positional every command takes.

    A run accumulates in one directory: each level reads what the level before it left there and
    writes its own beside it. Naming that directory once, positionally, is the whole invocation's
    shape; ``--from`` is the exception for when the input lives somewhere else."""
    p.add_argument("run", metavar="DIR",
                   help=f"the run directory: {what}. Created if needed, and read from as well as "
                        f"written to, so a pipeline names it once per command")


def _add_from_arg(g, what: str) -> None:
    """Add ``--from`` — where to read the previous level, when it is not the run directory."""
    g.add_argument("--from", dest="source", default=None, metavar="PATH",
                   help=f"read {what} from here instead of from the run directory. Use it for a "
                        f"tree or a run that came from somewhere else, or to write a run separate "
                        f"from the one it reads")


#: What a species tree resolves to inside a run directory, in the order tried: the grouped layout
#: first, then --flat. Both name the *complete* tree — every level runs on it, extinct lineages
#: included.
_TREE_IN_RUN = (os.path.join("species", "species_complete.nwk"), "species_complete.nwk")

#: What a genomes handoff resolves to: the grouped ``genomes/`` first, then a --flat directory.
_GENOMES_IN_RUN = ("genomes", "")


def resolve_tree(path: str) -> str:
    """Give back the species-tree file to open, from either a Newick file or a **run directory**.

    Spelling out ``out/species/species_complete.nwk`` is a detour through a layout the command
    already knows; the run directory says the same thing. A path that is not a directory is returned
    untouched, so any tree from anywhere still works."""
    if not os.path.isdir(path):
        return path
    for candidate in _TREE_IN_RUN:
        full = os.path.join(path, candidate)
        if os.path.exists(full):
            return full
    raise FileNotFoundError(
        f"{path} is a directory but holds no species tree — looked for "
        f"{' and '.join(_TREE_IN_RUN)}. Point it at a 'zombi2 species' run directory, or give a "
        f"Newick file with --from.")


def resolve_genomes(path: str) -> tuple[str, str]:
    """Give back ``(events directory, species-tree file)`` for a genomes run, in either layout.

    A genomes run is identified by its event log; the tree the events index against is the species
    tree of the same run, which `zombi2 genomes` guarantees is there — it writes the canonicalised
    one itself when its tree came from elsewhere."""
    for candidate in _GENOMES_IN_RUN:
        full = os.path.join(path, candidate) if candidate else path
        if os.path.exists(os.path.join(full, "genome_events.tsv")):
            return full, resolve_tree(path)
    raise FileNotFoundError(
        f"{path} holds no genomes run — looked for genome_events.tsv in {path}/genomes/ and in "
        f"{path} itself. Run 'zombi2 genomes' there first, or point --from at a run that has.")


def _log_value(value: object) -> str:
    """Render one parameter for the run log. A rate is recorded in its **written form**, so the log
    line can be pasted straight back into the flag (or a ``--params`` file) rather than being a repr
    the reader has to translate."""
    from zombi2.rates.modifiers import DrivenBy, Modifier
    from zombi2.rates.parse import written_form
    from zombi2.rates.rate import Rate
    from zombi2.rates.scope import Scope

    if isinstance(value, DrivenBy):
        # a bare DrivenBy is how the choice slots are written (--transfer-to), where there is no base
        # number to print; its repr is the same expression the flag takes, and it also round-trips as
        # a rate (a bare modifier is base 1.0), so both readings paste straight back in.
        return repr(value)
    if isinstance(value, (Rate, Scope, Modifier)):
        return written_form(value)
    return str(value)


def _write_params_log(path: str, args: argparse.Namespace, summary: str) -> None:
    """Write the full set of run parameters to ``path`` — always, for reproducibility."""
    lines = ["# ZOMBI2 run parameters",
             f"zombi2_version\t{__version__}",
             f"timestamp\t{datetime.datetime.now().isoformat(timespec='seconds')}",
             f"command_line\t{' '.join(sys.argv)}"]
    for key, value in sorted(vars(args).items()):
        lines.append(f"{key}\t{_log_value(value)}")
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
