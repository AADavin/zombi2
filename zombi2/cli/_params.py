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

    # ...or one file for a whole pipeline, with a table per subcommand
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


def load_params_file(path: str, valid_dests: set[str], command: str) -> dict:
    """Parse ``path`` into ``{dest: value}`` for ``command``, validating every key.

    ``valid_dests`` is the set of argument ``dest``s the subcommand accepts. A ``[command]`` table
    scopes the file to that subcommand; otherwise the file's top-level scalar keys are used (any
    other ``[section]`` tables are ignored, so one file can serve a whole pipeline). Every key must
    map to a known argument (a hyphen/underscore-insensitive match) — an unknown key is an error,
    so typos surface instead of being silently dropped.
    """
    with open(path, "rb") as f:
        data = tomllib.load(f)

    if isinstance(data.get(command), dict):        # a [command] table -> use it
        section = data[command]
    else:                                          # flat file -> top-level scalars only
        section = {k: v for k, v in data.items() if not isinstance(v, dict)}

    out: dict = {}
    unknown: list[str] = []
    for key, value in section.items():
        dest = key.replace("-", "_")
        if dest in valid_dests:
            out[dest] = value
        else:
            unknown.append(key)
    if unknown:
        raise ValueError(
            f"unknown parameter(s) in {path} for 'zombi2 {command}': "
            f"{', '.join(sorted(unknown))} — keys must be that command's option names")
    return out
