"""Fail a pull request that changes ZOMBI2's public surface without touching `CHANGELOG.md`.

`CONTRIBUTING.md` asks for a `[Unreleased]` line whenever a user could observe the change, and
nothing enforced it: #319 added `zombi2 tools tree --gamma` and `zombi2.tree.gamma_statistic` and
went green with the changelog untouched. Lint, types, tests and the docs build all pass for a pull
request that adds a flag and says nothing about it.

**What counts as the public surface**, and only this:

- every ``--flag`` handed to ``add_argument()`` under ``zombi2/cli/`` — what a user types;
- every name in an ``__all__`` under ``zombi2/`` — what a user imports.

**It compares the surface, not the diff.** Both revisions are parsed with `ast` and the two sets are
subtracted. Reformatting, a renamed local, a rewritten help string and a reworded docstring all leave
the surface identical and raise nothing — which matters, because the same pull request that added
`--gamma` also rewrote help text across six CLI modules, and a diff-line guard would have fired on
all of it and taught everyone to ignore the check.

Run it locally the way CI does:

    python scripts/check_public_surface.py --base origin/main
"""

from __future__ import annotations

import argparse
import ast
import subprocess
import sys
from collections.abc import Iterable, Mapping

CHANGELOG = "CHANGELOG.md"
PACKAGE = "zombi2/"
CLI = "zombi2/cli/"


def _module(path: str) -> str:
    """``zombi2/genomes/__init__.py`` -> ``zombi2.genomes``; ``zombi2/tree.py`` -> ``zombi2.tree``."""
    stem = path[: -len(".py")] if path.endswith(".py") else path
    if stem.endswith("/__init__"):
        stem = stem[: -len("/__init__")]
    return stem.replace("/", ".")


def _flags(tree: ast.Module) -> set[str]:
    """Every option string passed to an ``add_argument()`` call, anywhere in the module."""
    found: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "add_argument"):
            continue
        for arg in node.args:                      # argparse takes the option strings positionally
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str) \
                    and arg.value.startswith("-"):
                found.add(arg.value)
    return found


def _exports(tree: ast.Module) -> set[str]:
    """The names in a module-level ``__all__``, including any built up with ``+=``.

    A computed ``__all__`` is skipped rather than guessed at: this guard would rather miss a change
    than invent one, since a false alarm is what gets a check disabled.
    """
    found: set[str] = set()
    for node in tree.body:
        targets: Iterable[ast.expr]
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AugAssign):
            targets = [node.target]
        else:
            continue
        if not any(isinstance(t, ast.Name) and t.id == "__all__" for t in targets):
            continue
        try:
            value = ast.literal_eval(node.value)
        except (ValueError, SyntaxError):
            continue
        if isinstance(value, (list, tuple)):
            found |= {str(name) for name in value}
    return found


def public_surface(sources: Mapping[str, str]) -> set[str]:
    """The public surface of a revision, as labels, given ``{path: file text}``.

    Labels rather than bare names so the failure message says what kind of thing moved, and so a flag
    that migrates between CLI modules does not read as a removal and an addition — it is the same
    thing to type either way.
    """
    surface: set[str] = set()
    for path, text in sorted(sources.items()):
        if not path.startswith(PACKAGE) or not path.endswith(".py"):
            continue
        try:
            tree = ast.parse(text)
        except SyntaxError:                        # a revision we cannot parse tells us nothing
            continue
        if path.startswith(CLI):
            surface |= {f"flag {flag}" for flag in _flags(tree)}
        module = _module(path)
        surface |= {f"api {module}.{name}" for name in _exports(tree)}
    return surface


def _git(*args: str) -> str:
    return subprocess.run(("git", *args), check=True, capture_output=True, text=True).stdout


def _sources_at(ref: str | None) -> dict[str, str]:
    """Every `zombi2/**/*.py` at ``ref``, or in the working tree when ``ref`` is None."""
    if ref is None:
        paths = _git("ls-files", f"{PACKAGE}*.py").split()
        out = {}
        for path in paths:
            with open(path, encoding="utf-8") as fh:
                out[path] = fh.read()
        return out
    listing = _git("ls-tree", "-r", "--name-only", ref, PACKAGE).split()
    return {path: _git("show", f"{ref}:{path}")
            for path in listing if path.endswith(".py")}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", required=True,
                        help="the revision to compare against (e.g. origin/main)")
    args = parser.parse_args(argv)

    base = _git("merge-base", args.base, "HEAD").strip()
    before = public_surface(_sources_at(base))
    after = public_surface(_sources_at(None))

    added, removed = sorted(after - before), sorted(before - after)
    if not added and not removed:
        print("public surface unchanged — no CHANGELOG entry required")
        return 0

    changed = _git("diff", "--name-only", base, "--", CHANGELOG).strip()
    for label in added:
        print(f"  + {label}")
    for label in removed:
        print(f"  - {label}")
    if changed:
        print(f"\npublic surface changed and {CHANGELOG} was updated — good")
        return 0

    print(
        f"\nThis branch changes what a user can type or import, and does not touch {CHANGELOG}.\n"
        f"Add a line under `## [Unreleased]`, grouped by `### Added` / `### Changed` / `### Fixed`\n"
        "(CONTRIBUTING.md). If the change really is invisible to users, label the pull request\n"
        "`no changelog` and this check will be skipped.", file=sys.stderr)
    return 1


if __name__ == "__main__":                         # pragma: no cover
    raise SystemExit(main())
