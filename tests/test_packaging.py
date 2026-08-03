"""What the *installed* package promises, as opposed to what this checkout can do.

This file exists because of a bug no other test in this suite could have found. ``--params`` reads
TOML, and on Python 3.10 that means the `tomli` backport, which was not declared as a dependency.
The suite never noticed: `pytest` and `mypy` both require `tomli` below 3.11, so the dev extra
supplied it and CI's 3.10 job passed happily. A plain ``pip install zombi2`` on 3.10 did not supply
it, and the documented ``--params`` workflow died with a raw ``ModuleNotFoundError``.

The lesson generalises past that one package: **a dev environment cannot test its own completeness**.
So these tests read the declared metadata and compare it against what the code actually imports,
rather than against what happens to be importable here.
"""

from __future__ import annotations

import ast
import pathlib
import sys

import pytest

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PKG = _ROOT / "zombi2"

if sys.version_info >= (3, 11):
    import tomllib
else:                                     # the very fallback this file is about
    import tomli as tomllib               # type: ignore[no-redef]


def _pyproject() -> dict:
    path = _ROOT / "pyproject.toml"
    if not path.exists():
        pytest.skip("not running from a source checkout")
    with open(path, "rb") as f:
        return tomllib.load(f)


def _declared_runtime_names() -> set[str]:
    """The distribution names in ``[project] dependencies``, lowercased, markers stripped."""
    names = set()
    for spec in _pyproject()["project"]["dependencies"]:
        head = spec.split(";")[0]
        for sep in ("<=", ">=", "==", "~=", "!=", "<", ">", "["):
            head = head.split(sep)[0]
        names.add(head.strip().lower())
    return names


#: Modules that are standard library on *some* supported Python but not on the one running this test.
#: `sys.stdlib_module_names` answers for this interpreter alone, so on 3.10 `tomllib` — stdlib from
#: 3.11 — reads as a third-party import and the check below fails on the very version it was written
#: to protect. That is the same mistake in miniature as the bug this file exists for: a check that
#: silently asks "what is true here?" when the question is "what is true for the package".
_STDLIB_ON_A_LATER_PYTHON = {"tomllib"}


def _third_party_imports() -> set[str]:
    """Every top-level module the package imports that is neither ours nor in the standard library.

    Read from the AST rather than by importing, so this reports what the *code* needs on a machine
    that does not happen to have it — which is the whole point."""
    ours = {"zombi2"}
    stdlib = set(sys.stdlib_module_names) | _STDLIB_ON_A_LATER_PYTHON
    found: set[str] = set()
    for path in _PKG.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                found.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                found.add(node.module.split(".")[0])
    return {m for m in found if m not in ours and m not in stdlib}


def test_every_third_party_import_is_a_declared_dependency():
    """The guard the tomli bug walked straight past.

    An import that is not declared works fine here — something else in the dev extra dragged it in —
    and fails for the person who installed the package. Comparing the code against the *metadata*
    rather than against this interpreter is what makes the check meaningful."""
    undeclared = _third_party_imports() - _declared_runtime_names()
    # optional extras: imported behind a try/except or inside a function, for a feature that says so
    optional = {"phylustrator", "matplotlib", "scipy", "pytest"}
    assert not (undeclared - optional), (
        f"zombi2 imports {sorted(undeclared - optional)} but does not declare "
        f"{'it' if len(undeclared - optional) == 1 else 'them'} in [project] dependencies. This "
        f"passes here only because something in the dev extra supplies it; a plain "
        f"`pip install zombi2` would not.")


def test_tomli_is_declared_for_the_python_that_needs_it():
    """Named on its own because it is the one that got away, and because the marker matters as much
    as the name: declaring `tomli` unconditionally would install a redundant backport on 3.11+."""
    deps = _pyproject()["project"]["dependencies"]
    tomli = [d for d in deps if d.split(";")[0].strip().lower().startswith("tomli")]
    assert tomli, "tomli must be declared: the CLI's --params reader imports it on Python 3.10"
    assert "python_version" in tomli[0] and "3.11" in tomli[0], (
        f"tomli must carry the marker that keeps it off 3.11+, got {tomli[0]!r}")


def test_the_params_reader_still_needs_the_backport_this_declares():
    """If the floor ever rises to 3.11 the fallback goes, and so should the dependency. This fails
    then, which is the reminder to remove both rather than carry a dependency nobody needs."""
    source = (_PKG / "cli" / "_params.py").read_text(encoding="utf-8")
    floor = _pyproject()["project"]["requires-python"]
    if "3.10" in floor:
        assert "tomli" in source, "the 3.10 fallback vanished but the floor is still 3.10"
    else:
        assert "tomli" not in source, (
            "requires-python no longer includes 3.10, so the tomli fallback and the tomli "
            "dependency can both go")
