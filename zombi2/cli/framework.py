"""Shared CLI plumbing: the banner/help formatter, ``--params`` handling, the run-parameters log,
and the subcommand builder every command module leans on."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import os
import platform
import shlex
import shutil
import sys
import textwrap

import numpy as np

from zombi2 import __version__
from zombi2.rng import draw_seed


_DESCRIPTION = """\
Simulate each level of evolution on its own. Run 'zombi2 <command> -h' for a command's options.

Levels
  species              simulate a dated species tree
  genomes              evolve gene families along a species tree (family, ordered, nucleotide)
  sequences            evolve sequences down each gene tree (a prior genomes run)
  traits               evolve a trait along a species tree (continuous or discrete)

Joint
  joint                grow a species tree and the level driving it, in one pass

Tools
  tools                analyses that read a finished run (homology, markers, recphylo, …)
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
#: still lists (by name), so the help can never fall behind a level's ``IMPLEMENTED_MODIFIERS`` declaration
#: of what it supports.
#: ``DrivenBy`` carries no snippet: its source is a *file* on a conditioned level and a *live level*
#: name on ``joint``, so one snippet would be wrong on one of the two screens (it was — the joint help
#: showed a filename, which ``joint`` rejects). A command that wants one passes ``example=``.
_MODIFIER_HELP = {
    "OnTime": ("OnTime({0: 1.0, 3: 0.3})", "the rate changes in time — a skyline"),
    "OnTotalDiversity": ("OnTotalDiversity(cap=100)", "the rate slows as the clade fills up"),
    "FromParent": ("FromParent(spread=0.2)", "the rate drifts down the tree"),
    # "clock" is reserved for the sequences by-lineage substitution modifier (SPEC §7), and this
    # string now prints on `zombi2 species -h` and `zombi2 genomes -h` too
    "ByLineage": ("ByLineage(spread=0.3)", "one independent draw per lineage — uncorrelated"),
    "DrivenBy": (None, "the rate is driven by another level"),
}


def _wrap_note(note: str, width: int = 86) -> list[str]:
    return textwrap.wrap(note, width=width, initial_indent="  ", subsequent_indent="  ")


def _rates_help(supported, flag: str, *, scopes: str | None = None, note: str | None = None,
                example: str | None = None) -> str:
    """The ``RATES`` epilog block for a command, built from that level's ``IMPLEMENTED_MODIFIERS``.

    Listing what the engine *declares* (rather than a hand-kept list) is what keeps the help honest:
    a modifier the level does not support is rejected by the engine, so it must not be advertised
    here either — and the worked example is drawn from the same list, so it is always a modifier that
    runs. ``example`` overrides that snippet for a level whose form differs (see `_MODIFIER_HELP`).
    """
    snippets = [_MODIFIER_HELP[m.__name__][0] for m in supported if m.__name__ in _MODIFIER_HELP]
    shown = example or next((s for s in snippets if s), None)

    header = _BOLD + "RATES" + _RESET if _use_color() else "RATES"
    intro = "  A rate is a number, or a quoted expression — the same text Python takes"
    lines = [header, intro + (":" if shown or scopes else ".")]
    if shown:
        lines.append(f'    {flag} "1.0 * {shown}"')
    if scopes:
        lines.append(f'    {flag} "{scopes}"')
    lines.append("  Modifiers this level takes (anything else is an error):")
    for m in supported:
        name = m.__name__
        entry = _MODIFIER_HELP.get(name)
        lines.append(f"    {name:<20}{entry[1]}" if entry else f"    {name}")
    if note:
        lines += _wrap_note(note)
    return "\n".join(lines)


def _add_params_arg(g) -> None:
    """Add ``--params FILE`` (a TOML parameters file) to a subcommand's ``general`` group."""
    g.add_argument("--params", metavar="FILE",
                   help="a TOML parameters file; keys are long option names, overridden by any "
                        "flag you also pass")


def _add_flat_arg(g) -> None:
    """Add ``--flat`` (write everything into one directory) to a subcommand's ``outputs`` group."""
    g.add_argument("--flat", action="store_true",
                   help="write every file straight into DIR, not grouped by level")


def default_outputs(result) -> tuple[str, ...]:
    """What a result writes when ``--write`` is not given — read off its own ``write()`` signature.

    The CLI has to know the default before it calls ``write``, because it routes the one-file-per-family
    outputs into their own directories. Reading it here rather than repeating it keeps the two from
    drifting: a level that changes what it writes by default changes it in one place."""
    import inspect

    return tuple(inspect.signature(result.write).parameters["outputs"].default)


def _add_quiet_arg(g) -> None:
    """Add ``--quiet`` — no progress bar, for a log file or a batch of runs."""
    g.add_argument("--quiet", action="store_true", help="no progress bar")


def _add_force_arg(g) -> None:
    """Add ``--force`` — re-run this level even though a later level built from it is already in the run
    directory, removing that now-stale downstream output (see `check_stale_downstream()`)."""
    g.add_argument("--force", action="store_true",
                   help="re-run this level even if a later level in the run was built from its "
                        "output, deleting that now-stale downstream (without it, the command "
                        "refuses)")


def _add_parallel_arg(g) -> None:
    """Add ``--parallel [N]`` — evolve the run's independent units concurrently. Opt-in: omitting it
    runs serially, the default."""
    g.add_argument("--parallel", nargs="?", const="__all__", default=None, metavar="N",
                   help="evolve independent units — a gene family (genomes), a gene tree "
                        "(sequences) — one per worker process: bare = every core, N = N workers "
                        "(default: serial). A separate engine, so the result is the same for any "
                        "worker count but differs from a serial run of the same seed (both valid)")


def parallel_from_args(args, parser):
    """Turn the ``--parallel`` flag into the API's ``parallel`` value: ``False`` (omitted, serial),
    ``True`` (bare, every core), or a positive worker count. A non-positive or non-integer value is a
    usage error, named here rather than left to raise deeper in."""
    v = args.parallel
    if v is None:
        return False
    if v == "__all__":                                  # bare --parallel
        return True
    try:
        n = int(v)
    except (TypeError, ValueError):
        parser.error(f"--parallel takes a worker count (a positive integer) or no value for all "
                     f"cores, got {v!r}")
    if n < 1:
        parser.error(f"--parallel needs a positive worker count, got {n}")
    return n


def level_dir(output: str, level: str, flat: bool) -> str:
    """Where one level's files belong: ``<output>/<level>/``, or ``<output>/`` under ``--flat``.

    Grouping is the CLI's business, not the engines': a ``Result.write`` writes whatever it is given
    into the one directory it is handed, and the layout is chosen here by calling it more than once.
    Created on the way out, so a caller can write into it immediately."""
    path = output if flat else os.path.join(output, level)
    os.makedirs(path, exist_ok=True)
    return path


def guidance(args, *looks: str) -> None:
    """One or more one-line pointers printed after a run — the file(s) worth looking at, for someone
    finding their way around a run. It only says *what was written and where*, never what to do next
    (a run is not one road: after a species tree you might grow genomes, or a trait, or stop). It is
    suppressed by ``--quiet`` (a scripted batch stays quiet); the ``wrote …`` summary still prints."""
    if getattr(args, "quiet", False):
        return
    for look in looks:
        print(f"  → {look}")


def signpost(args, report_path, *level_dirs) -> None:
    """Point the user at what the run just wrote: one line per output directory listing the files it
    holds (a per-family directory named as itself), then the run report. The file-by-file descriptions
    live in ``run.zombi2`` — the names it lists come from the same place (`output_signposts`), so the
    terminal and the report never disagree. Record files (``.log`` / ``_summary.json``) are left out;
    they are named in the report. Suppressed by ``--quiet``; the ``wrote …`` line still prints."""
    if getattr(args, "quiet", False):
        return
    from zombi2._runtime.report import output_signposts
    for d in dict.fromkeys(level_dirs):              # dedup: --flat collapses every level into one dir
        names = [name for name, _gloss in output_signposts(d)]
        if names:
            print(f"  → {d}/  " + " · ".join(names))
    if report_path:
        print(f"  → {report_path}  ← open this first: a report of every file, and how to reproduce")


def defaults_used(args, **fallbacks) -> None:
    """Fill any of ``fallbacks`` the caller left unset on ``args`` (silently).

    A bare ``zombi2 species out/`` runs, so it shows a newcomer the shape of a command instead of a
    list of what they failed to supply — the rates are precisely the part nobody can guess on a first
    run, and they are round, illustrative values calibrated to nothing. They are not announced on
    stderr — a paragraph-long warning on every bare run was noise — but they are not hidden either: the
    ``.log`` and ``run.zombi2`` record every resolved value, and the report's parameters and
    TO REPRODUCE block make the run reconstructable. The number is visible where it is useful (in the
    record) rather than shouted where it is not (on the terminal).

    Not every level runs bare: ``traits`` still requires ``--kind``, because it decides which options
    apply at all, so there is nothing sensible to default it to."""
    for name, value in fallbacks.items():
        if getattr(args, name, None) is None:
            setattr(args, name, value)


def warn_if_fates_were_inferred(tree, args) -> None:
    """Say so when a tree read from ``--from`` had its tip fates *guessed* and some came out extinct.

    A tree whose nodes all carry labels is taken to be a ZOMBI tree, so its tip depths are trusted:
    a tip sitting below the present is read as an extinct lineage. That is right for a complete tree
    ZOMBI wrote, and quietly wrong for one that has been re-estimated, re-dated or edited since,
    where the depths differ by *noise* — every tip then reads as extinct, and the run continues on a
    handful of survivors, or one, without complaint. That is not recoverable from the outputs unless
    you happen to count the genomes.

    So it is announced rather than refused: refusing would reject the complete trees this heuristic
    exists to serve. Passing ``--tip-fates`` states the fates instead of inferring them, and silences
    this — a species run's own ``species_fates.tsv`` is already in that format."""
    if getattr(args, "tip_fates", None) or not getattr(args, "source", None):
        return
    tips = [n for n in tree.nodes.values() if n.children is None]
    dead = [n for n in tips if n.fate != "extant"]
    if dead and len(dead) > len(tips) // 2:
        warn(f"{len(dead)} of {len(tips)} tips in this tree sit below the present, so they were read "
             f"as extinct lineages and only {len(tips) - len(dead)} genome(s) will be simulated. That "
             f"is what a complete ZOMBI tree looks like — but if this tree was re-estimated, re-dated "
             f"or edited, the depths differ by noise rather than by extinction and the fates are "
             f"wrong. Declare them with --tip-fates FILE to be sure.")


def resolve_seed(args) -> None:
    """Draw a seed when the caller did not give one, so the run can be regenerated from its own log.

    Every run is seeded — leaving ``--seed`` out only meant the seed came from the OS and was then
    thrown away, so the ``.log`` recorded ``seed None`` and the dataset could never be reproduced by
    anyone, its author included. Nothing warned them. Drawing it here costs nothing and makes the log
    a complete record: re-running the command it contains reproduces the run exactly.

    A run that wants fresh randomness still gets it — this draws from the OS each time. What changes
    is only that the draw is written down."""
    if getattr(args, "seed", None) is None:
        args.seed = draw_seed()      # a 32-bit value, so the log's number is one a reader can retype
        # …and say so. The seed was always written into the report, but the report is a file nobody
        # has opened yet: on screen an unseeded run looked exactly like a seeded one. A class of
        # students following a worksheet got 17, 19 and 15 families where the sheet said 20, and their
        # only clue was that the answers were wrong. One line on stderr, surviving --quiet for the
        # same reason `warn` does — a batch job's log is exactly where this number needs to end up.
        print(f"zombi2: no --seed given; drew {args.seed} (recorded in the run report)",
              file=sys.stderr)


def warn(message: str) -> None:
    """A diagnostic about the *result* — the run succeeded, but something about what came out is
    likely not what was wanted. Unlike `guidance()` it goes to **stderr** and survives
    ``--quiet``: a scripted batch is precisely the caller who needs to hear that its data is
    degenerate, and stdout stays clean for the ``wrote …`` line. A healthy run prints nothing here,
    so an empty stderr still means "nothing to report"."""
    print(f"zombi2: warning: {message}", file=sys.stderr)


#: The fixed pipeline edges — a level → the levels that read its output *directly*. ``species`` feeds
#: ``genomes`` and ``traits`` (both read the species tree); ``genomes`` feeds ``sequences`` (the gene
#: trees). Re-running a level orphans everything reachable from it here, plus any level that recorded a
#: DrivenBy *conditioning* on it (a dynamic edge — see `_conditioned_on()`).
_STRUCTURAL = {
    "species": ("genomes", "traits"),
    "genomes": ("sequences",),
}

#: Levels whose rates can be conditioned on another level (they take a ``DrivenBy``), so they may carry
#: a ``conditioned_on`` record.
_CONDITIONABLE = ("genomes", "sequences", "traits")

#: Every level, in pipeline order — a stable order for listing them in a message.
_LEVEL_ORDER = ("species", "genomes", "sequences", "traits")

#: The marker a conditioned level writes, naming the levels its rates read via ``DrivenBy`` — the
#: dynamic half of the staleness graph (the fixed half is `_STRUCTURAL`).
_CONDITIONED_ON_FILE = "conditioned_on"


def _level_present(run: str, level: str) -> bool:
    """Whether ``run/<level>/`` holds a level's output — a non-empty grouped sub-directory. (The check
    is grouped-only: ``--flat`` commingles every level in one directory, where they cannot be told
    apart, so a ``--flat`` run is left to the user — see `check_stale_downstream()`.)"""
    d = os.path.join(run, level)
    return os.path.isdir(d) and bool(os.listdir(d))


def _conditioned_on(run: str, level: str) -> set:
    """The levels ``run/<level>/`` recorded a ``DrivenBy`` conditioning on (its rates read their
    output), from the `_CONDITIONED_ON_FILE` marker — empty if it conditioned on nothing."""
    p = os.path.join(run, level, _CONDITIONED_ON_FILE)
    if not os.path.exists(p):
        return set()
    with open(p, encoding="utf-8") as f:
        return {ln.strip() for ln in f if ln.strip()}


def _drivenby_sources(spec):
    """Every ``DrivenBy`` source in a rate spec — whether ``spec`` is a bare ``DrivenBy`` (a choice slot
    like ``transfer_to``) or a rate carrying ``DrivenBy`` modifiers. A plain number yields none."""
    from zombi2.rates.modifiers import DrivenBy
    if isinstance(spec, DrivenBy):
        yield spec.source
    for m in getattr(spec, "modifiers", ()):
        if isinstance(m, DrivenBy):
            yield m.source


def conditioned_levels(run: str, rate_specs) -> set:
    """Which of THIS run's levels the given rate specs are conditioned on — a same-run ``DrivenBy`` file
    maps to the level whose directory holds it (``run/traits/trait_events.tsv`` → ``traits``). A driver
    file outside the run does not count: re-running a level of this run cannot make it stale."""
    run_abs = os.path.abspath(run)
    levels = set()
    for spec in rate_specs:
        for src in _drivenby_sources(spec):
            if isinstance(src, str):
                first = os.path.relpath(os.path.abspath(src), run_abs).split(os.sep)[0]
                if first not in (os.curdir, os.pardir):    # a file under one of this run's level dirs
                    levels.add(first)
    return levels


def record_conditioning(level_out: str, driver_levels) -> None:
    """Write (or clear) the `_CONDITIONED_ON_FILE` marker in a level's output directory: the
    same-run levels its rates read via ``DrivenBy``, so the guard knows re-running one of them orphans
    this level. Removes a stale marker when a re-run conditions on nothing."""
    driver_levels = sorted(set(driver_levels))
    p = os.path.join(level_out, _CONDITIONED_ON_FILE)
    if driver_levels:
        with open(p, "w", encoding="utf-8") as f:
            f.write("\n".join(driver_levels) + "\n")
    elif os.path.exists(p):
        os.remove(p)


def _stale_downstream(args, level: str) -> list:
    """The downstream levels already present that re-running ``level`` would orphan — everything
    reachable from it by a 'reads its output' edge (the fixed pipeline `_STRUCTURAL` plus recorded
    ``DrivenBy`` conditioning), listed in pipeline order. Empty under ``--flat`` (not guarded)."""
    if getattr(args, "flat", False):
        return []
    run = args.run
    edges = {k: set(v) for k, v in _STRUCTURAL.items()}
    for consumer in _CONDITIONABLE:                          # a conditioning edge: driver → consumer
        for driver in _conditioned_on(run, consumer):
            if driver == consumer:
                # A trait driving another trait: both sit under `traits/`, so the edge is
                # traits → traits and this graph, whose nodes are LEVELS, cannot tell the driver
                # from the target. Following it would mark the level stale against itself and
                # `--force` would then delete every trait in the run, including the one just
                # written. The marker is still recorded — it is what the log and the run report
                # show — but a within-level dependency is not invalidated automatically until the
                # graph's nodes are `traits/<name>` rather than `traits`.
                continue
            edges.setdefault(driver, set()).add(consumer)
    seen, stack = set(), list(edges.get(level, ()))
    while stack:
        d = stack.pop()
        if d not in seen:
            seen.add(d)
            stack.extend(edges.get(d, ()))
    return [d for d in _LEVEL_ORDER if d in seen and _level_present(run, d)]


def check_stale_downstream(args, level: str) -> None:
    """Refuse, up front, to re-run ``level`` when the run directory already holds a later level built
    from it: re-running would leave that downstream output silently mismatched with the new one. The
    normal forward pipeline never trips this (each level is run once, downstream not yet there); only a
    re-run does. ``--force`` allows it — the orphaned downstream is removed afterwards by
    `clear_stale_downstream()`, so the run can never end up a stale mix. Raises ``ValueError`` (the
    CLI reports it as ``zombi2: error: …``); does nothing under ``--flat`` (see `_level_present()`)."""
    if getattr(args, "force", False):
        return
    present = _stale_downstream(args, level)
    if present:
        names = ", ".join(present)
        raise ValueError(
            f"{args.run}/ already holds a downstream {names} run built from this {level}; re-running "
            f"{level} would leave {'it' if len(present) == 1 else 'them'} stale. Write to a fresh "
            f"directory, or pass --force to re-run {level} and remove the now-stale {names}.")


def clear_stale_downstream(args, level: str) -> None:
    """With ``--force``, remove the downstream levels re-running ``level`` has now orphaned, so the run
    is left consistent (this level fresh, nothing stale beneath it). Called after the level's own run
    succeeds, so a failed re-run never deletes the old downstream. A no-op without ``--force`` or under
    ``--flat``, and silent unless it actually removed something — but never silent when it did."""
    if not getattr(args, "force", False):
        return
    present = _stale_downstream(args, level)
    for d in present:
        shutil.rmtree(os.path.join(args.run, d))
    if present:
        # NOT gated on --quiet, and on stderr: this is the one place a command deletes work the user
        # already has, and --quiet means "no progress bar", not "say nothing when you delete a run".
        # Every pipeline passes --force (a re-run re-enters a directory that holds its downstream),
        # so this is exactly the caller who must still hear it.
        warn(f"--force removed the now-stale {', '.join(present)} — rebuilt from the new {level}")


def _run_dir(s: str) -> str:
    """Normalise the run-directory argument: drop trailing slashes, so the completion line reads
    ``wrote out/`` and not ``wrote out//`` when the user wrote ``out/`` (as the quickstart does). Path
    semantics are unchanged — ``Path`` and ``os.path.join`` treat ``out`` and ``out/`` alike; this only
    tidies how the directory is echoed back. An all-slashes argument keeps a single slash."""
    return s.rstrip("/") or "/"


def _add_run_arg(p, what: str = "") -> None:
    """Add the run directory — the one positional every command takes.

    ``what`` (a per-command gloss) is accepted but not shown: the help line is the same three words
    on every screen. What a run directory *is* — written to and read from, named once per command —
    is the manual's job (Ch2), not ten repetitions of it in ``--help``."""
    p.add_argument("run", metavar="DIR", type=_run_dir, help="the run directory")


def _add_from_arg(g, what: str, detail: str = "") -> None:
    """Add ``--from`` — where to read the previous level, when it is not the run directory.

    ``what`` is a bare noun phrase and ``detail`` a parenthetical. They are separate because the two
    read badly fused: an appositive inside ``what`` landed mid-sentence ("read the species tree — a
    Newick file, or another run's directory from here instead of...")."""
    g.add_argument("--from", dest="source", default=None, metavar="PATH",
                   help=f"read {what} from PATH" + (f" ({detail})" if detail else ""))


#: What a species tree resolves to inside a run directory, in the order tried: the grouped layout
#: first, then --flat. Both name the *complete* tree — every level runs on it, extinct lineages
#: included.
_TREE_IN_RUN = (os.path.join("species", "species_complete.nwk"), "species_complete.nwk")

#: What a genomes handoff resolves to: the grouped ``genomes/`` first, then a --flat directory.
_GENOMES_IN_RUN = ("genomes", "")


def resolve_tree(path: str, *, is_run_dir: bool = False) -> str:
    """Give back the species-tree file to open, from either a Newick file or a **run directory**.

    Spelling out ``out/species/species_complete.nwk`` is a detour through a layout the command
    already knows; the run directory says the same thing. A path that is not a directory is returned
    untouched, so any tree from anywhere still works.

    ``is_run_dir`` says the path came from the run argument rather than ``--from``, which changes
    what a missing path means: not "your tree file is not there" but "this run has no species tree
    yet". That is the first mistake anyone makes — typing the genomes command before the species
    one — and it used to fall through to a bare ``tree file not found``."""
    if is_run_dir and not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} does not exist yet, so it holds no species tree to evolve along. Run "
            f"'zombi2 species {path} ...' first to grow one, or point --from at a Newick file or "
            f"another run.")
    if not os.path.isdir(path):
        return path
    for candidate in _TREE_IN_RUN:
        full = os.path.join(path, candidate)
        if os.path.exists(full):
            _refuse_partial_species_run(full)
            return full
    raise FileNotFoundError(
        f"{path} is a directory but holds no species tree — looked for "
        f"{' and '.join(_TREE_IN_RUN)}. Point it at a 'zombi2 species' run directory, or give a "
        f"Newick file with --from.")


#: What a *finished* species run leaves beside its tree. A run directory holding the tree and neither
#: of these did not finish writing.
_SPECIES_RUN_MARKERS = ("species_fates.tsv", "species.log")


def _refuse_partial_species_run(tree_path: str) -> None:
    """Refuse a run directory whose species tree is there but whose run plainly did not finish.

    A downstream level reads each tip's fate from ``species_fates.tsv`` and falls back to guessing
    from tip depth when it is absent — which is right for an **external** tree (a published phylogeny
    has no fates file and never will) and wrong for a ZOMBI run that was interrupted, rsynced
    halfway, or tidied. The two were indistinguishable to the reader, and the failure was silent and
    material: under ``--sampling`` every unsampled tip reads back as sampled, so the same seed and
    rates gave 30 extant genomes where the complete run gave 14, with exit 0 and nothing on stderr.

    What tells them apart is that this is a run *directory* with the run's own file names in it. A
    run that finished wrote a log beside the tree; one that did not, did not. So: tree present,
    neither marker present → stop. An external tree still goes in the way it always did, as a file
    path (or ``--from tree.nwk``), which never reaches here."""
    directory = os.path.dirname(tree_path)
    if any(os.path.exists(os.path.join(directory, m)) for m in _SPECIES_RUN_MARKERS):
        return
    raise FileNotFoundError(
        f"{tree_path} is there but {' and '.join(_SPECIES_RUN_MARKERS)} are not, so this looks like "
        f"a species run that did not finish writing (an interrupted run, a partial copy). Reading it "
        f"would treat every tip as sampled and extant, which silently changes the taxon set your "
        f"genomes are built over. Re-run 'zombi2 species' here — or, if this is an external tree, "
        f"pass the .nwk file itself with --from rather than the directory.")


def _read_tip_fates(path: str) -> dict:
    """Parse a ``--tip-fates`` file into ``{tip_name: fate}``: one
    ``tip_name<TAB>extant|extinct|unsampled`` row per tip (whitespace also accepted; blank lines and
    ``#`` comments skipped). This is the same shape ``species_fates.tsv`` is written in, so that output
    feeds straight back in — its ``lineage<TAB>fate`` header row is recognised and skipped. The values
    are checked against the tree by `read_newick()`."""
    fates = {}
    try:
        with open(path, encoding="utf-8") as f:
            for lineno, raw in enumerate(f, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t") if "\t" in line else line.split()
                if parts == ["lineage", "fate"]:
                    continue  # the species_fates.tsv header, so that file is a valid --tip-fates input
                if len(parts) != 2:
                    raise ValueError(f"{path}:{lineno}: expected 'tip_name<TAB>extant|extinct|unsampled', "
                                     f"got {raw.rstrip()!r}")
                fates[parts[0]] = parts[1]
    except FileNotFoundError:
        raise FileNotFoundError(f"tip-fates file not found: {path}") from None
    return fates


def sibling_fates(tree_path: str) -> dict | None:
    """The ``species_fates.tsv`` a species run writes next to its tree, parsed, or ``None`` if absent.
    Lets a downstream level read each tip's fate from the run — telling extinct and unsampled tips
    apart — instead of guessing it from tip depth. A user's explicit ``--tip-fates`` overrides this."""
    cand = os.path.join(os.path.dirname(tree_path), "species_fates.tsv")
    return _read_tip_fates(cand) if os.path.exists(cand) else None


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


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def input_digests(*values) -> list[tuple[str, str]]:
    """``(path, sha256)`` for every input file behind ``values``, deduplicated, in the order given.

    What a run *read* is as much a parameter of it as any flag, and the file's **name** does not
    pin it down: two runs from two different trees log the same ``--from tree.nwk`` and the same
    seed, and nothing in the log tells them apart. The digest does. A ``value`` is either a path or
    a rate, whose ``DrivenBy`` file source — the level a conditioned run was driven by — is as much
    an input as the tree. Anything that is not an existing file is skipped."""
    from zombi2.rates.modifiers import DrivenBy

    paths: list[str] = []
    for value in values:
        if isinstance(value, str):
            paths.append(value)
        elif value is not None:
            mods = getattr(value, "modifiers", (value,))     # a rate's modifiers, or a bare one
            paths.extend(m.source for m in mods
                         if isinstance(m, DrivenBy) and isinstance(m.source, str))
    out, seen = [], set()
    for path in paths:
        if path not in seen and os.path.isfile(path):
            seen.add(path)
            out.append((path, _sha256(path)))
    return out


def _command_line(args) -> str:
    """The command as typed — ``zombi2 <argv>`` — from the argv `main()` captured on ``args._argv``.
    That is the real invocation both from a shell and from a test that calls ``main([...])`` directly
    (where ``sys.argv`` would be pytest's). Lets the report list only the flags the user actually gave.
    Each token is shell-quoted, so a rate expression (``--birth "1.0 * OnTime({0: 1.0})"``) pastes back
    intact rather than globbing on the ``*`` or breaking on the braces."""
    argv = getattr(args, "_argv", None)
    if argv is None:
        return " ".join(sys.argv)
    return "zombi2 " + " ".join(shlex.quote(tok) for tok in argv)


def _write_params_log(path: str, args: argparse.Namespace, summary: str, effective=None,
                      inputs=(), omit=()) -> None:
    """Write the full set of run parameters to ``path`` — always, for reproducibility. ``effective``
    overrides the logged value of the named args with the **resolved** value the run actually used
    (e.g. a model's default ``kappa`` in place of the bare ``None`` that was on the command line), so
    the log reproduces the run without the reader having to know each default. ``inputs`` are the
    files the run read, from `input_digests()`, recorded by content as well as by name. ``omit`` drops
    options this run could not have taken — a command whose flags depend on a mode records the ones
    it *has*, so the log is the run's parameters rather than the parser's."""
    lines = ["# ZOMBI2 run parameters",
             f"zombi2_version\t{__version__}",
             # the software environment, so the record says what produced the run, not only with which
             # flags — a bug report or an irreproducible result usually turns on one of these three.
             f"python_version\t{platform.python_version()}",
             f"numpy_version\t{np.__version__}",
             f"platform\t{platform.system()} {platform.machine()}",
             f"timestamp\t{datetime.datetime.now().isoformat(timespec='seconds')}",
             f"command_line\t{_command_line(args)}"]
    for in_path, digest in inputs:
        lines.append(f"input\t{digest}\t{in_path}")
    for key, value in sorted({**vars(args), **(effective or {})}.items()):
        if key not in omit and not key.startswith("_"):     # `_argv` etc. are bookkeeping, not parameters
            lines.append(f"{key}\t{_log_value(value)}")
    lines.append(f"result\t{summary}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _add_subcommand(sub, name: str, help: str, description: str, usage: str, adder,
                    epilog: str | None = None):
    """Register a subcommand with the house-style formatter and a hand-written compact usage.

    The command list itself is curated (grouped by theme) in the top-level description, so the
    per-command ``help`` is suppressed from argparse's auto listing to avoid a duplicate dump.
    ``epilog`` (built with `_examples()`) adds a worked-example block below the options.
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
    for i, tok in enumerate(tokens[1:], 1):        # last --params wins, matching argparse and the log
        if tok == "--params" and i + 1 < len(tokens):
            path = tokens[i + 1]
        elif tok.startswith("--params="):
            path = tok.split("=", 1)[1]
    if path is None:
        return
    from zombi2.cli._params import load_params_file
    action_by_dest = {a.dest: a for a in subp._actions}
    try:
        overrides = load_params_file(path, set(action_by_dest), tokens[0], set(sub.choices))
    except (OSError, ValueError) as e:              # missing file, TOML error, or unknown key/section
        subp.error(str(e))
    given = {opt for tok in tokens for opt in [tok.split("=", 1)[0]]}  # option strings on the CLI
    for dest, val in list(overrides.items()):
        action = action_by_dest.get(dest)
        if action is None:
            continue
        # a variable-length option (nargs '+'/'*', e.g. --write) is a list on the command line; accept
        # a bare scalar in the file and wrap it, so `write = "events"` works like `--write events`.
        if action.nargs in ("+", "*") and not isinstance(val, list):
            overrides[dest] = val = [val]
        # an 'append' option (e.g. --mass-extinction) appends to its default, so a params default plus
        # a command-line flag would concatenate rather than override; when the flag is on the command
        # line, drop the params default so the command line alone is used.
        if isinstance(action, argparse._AppendAction) and given.intersection(action.option_strings):
            del overrides[dest]
            continue
        # a params value bypasses argparse's `choices=` check (that runs on command-line tokens, not on
        # defaults), so validate it here rather than let a bad value crash deep in the command. A
        # variable-length option holds a list, each element of which must be a valid choice.
        if action.choices is not None:
            bad = [v for v in (val if isinstance(val, list) else [val]) if v not in action.choices]
            if bad:
                subp.error(f"argument {action.option_strings[-1]}: invalid choice {bad[0]!r} in "
                           f"--params (choose from {', '.join(map(str, action.choices))})")
    subp.set_defaults(**overrides)
