"""Read a TOML ``--params`` file into ``{dest: value}`` overrides for a CLI subcommand.

A parameters file collects the *model* parameters of a run (rates, tree size, event knobs) so a
configuration is reproducible and shareable — the modern successor to ZOMBI1's ``*Parameters.tsv``.
It is applied as **defaults**, so any flag given on the command line still wins; the required I/O
paths (``-o``, ``-t``) stay on the command line.

The file is TOML. Keys are the command's long option names (hyphens or underscores both work, since
the option name *is* the API keyword), and values are native TOML scalars/arrays::

    # a flat file for one subcommand
    birth      = 1.0
    death      = 0.3
    total-time = 5.0
    seed       = 42

    # ...or one file for a whole pipeline, with a table per subcommand. A key outside any table is a
    # shared base broadcast to every command (here every level uses seed = 42); the command's own
    # table overrides it.
    seed = 42

    [species]
    birth = 1.0
    death = 0.3

    [genomes]
    duplication = 0.2
    write       = ["events", "profiles"]

``tomllib`` is used on Python 3.11+; the tiny ``tomli`` backport on 3.10.
"""

from __future__ import annotations

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # pragma: no cover - exercised only on 3.10
    import tomli as tomllib  # type: ignore[no-redef]


def _looks_like_an_eaten_path(path: str) -> bool:
    """Whether the file holds a backslash inside an ordinary ``"…"`` TOML string.

    TOML reads a backslash there as an escape, so a Windows path in a value is mangled — or, more
    often, rejected with a message about a hex value that says nothing about paths. The fix is a TOML
    *literal* string, and this is what lets the error say so."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return any("\\" in line.split("=", 1)[1] for line in f if '="' in line.replace(" ", ""))
    except OSError:
        return False


def load_params_file(path: str, valid_dests: set[str], command: str,
                     known_commands: set[str]) -> dict:
    """Parse ``path`` into ``{dest: value}`` for ``command``, validating every key.

    ``valid_dests`` is the set of argument ``dest``s the subcommand accepts; ``known_commands`` is
    every ``zombi2`` subcommand, so a mistyped ``[table]`` can be caught. The file is read as:

    - **``[section]`` tables** each scope keys to one subcommand. A section whose name is not a known
      command is an error (``[speces]`` is a typo, not a silently-dropped block). The ``[command]``
      table's keys are this run's, and each must be a known option of it.
    - **top-level scalar keys** (outside any table) are a **shared base broadcast to every command**;
      the ``[command]`` table overrides them on conflict. In a *pipeline* file (one with tables) a
      top-level key that this command does not accept is skipped — it is meant for another command —
      but in a *flat* file (no tables at all) every top-level key must be this command's, so a typo
      still surfaces instead of being silently dropped.
    """
    try:
        with open(path, "rb") as f:
            data = tomllib.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"params file not found: {path}") from None
    except tomllib.TOMLDecodeError as e:
        hint = ("; a value containing a path needs a TOML literal string ('''...'''), because an "
                "ordinary \"...\" string reads its backslashes as escapes"
                if _looks_like_an_eaten_path(path) else "")
        raise ValueError(f"{path} is not valid TOML: {e}{hint}") from None

    tables = {k: v for k, v in data.items() if isinstance(v, dict)}
    top_scalars = {k: v for k, v in data.items() if not isinstance(v, dict)}

    bad_tables = sorted(t for t in tables if t not in known_commands)
    if bad_tables:
        raise ValueError(
            f"unknown section(s) in {path}: {', '.join('[' + t + ']' for t in bad_tables)} — not a "
            f"zombi2 command; expected one of {', '.join(sorted(known_commands))}")

    out: dict = {}
    unknown: list[str] = []
    # top-level scalars first (the shared base); the command's own table overrides them below
    for key, value in top_scalars.items():
        dest = key.replace("-", "_")
        if dest in valid_dests:
            out[dest] = value
        elif not tables:                           # flat file: a top-level key must be this command's
            unknown.append(key)
        # else: pipeline file — a shared key this command lacks is for another command, so skip it
    for key, value in tables.get(command, {}).items():
        dest = key.replace("-", "_")
        if dest in valid_dests:
            out[dest] = value                      # the [command] table wins over the shared base
        else:
            unknown.append(key)
    if unknown:
        raise ValueError(
            f"unknown parameter(s) in {path} for 'zombi2 {command}': "
            f"{', '.join(sorted(unknown))} — keys must be that command's option names")
    return out
