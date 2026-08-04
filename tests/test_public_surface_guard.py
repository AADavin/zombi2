"""The guard that asks for a CHANGELOG line when the public surface moves.

`scripts/check_public_surface.py` exists because #319 added `zombi2 tools tree --gamma` and
`zombi2.tree.gamma_statistic` and went green with `[Unreleased]` empty. The two tests that matter are
the two failure modes of a guard like this: missing the thing it was written for, and crying wolf
often enough that people turn it off. Both are pinned below with the real diff that motivated it.
"""

from __future__ import annotations

import importlib.util
import pathlib

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "check_public_surface.py"
_spec = importlib.util.spec_from_file_location("check_public_surface", _SCRIPT)
assert _spec and _spec.loader
check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check)
public_surface = check.public_surface


# The shape of the real thing: a CLI module whose parser adds flags, and a library module with
# __all__. Kept small so the assertions read.
_CLI_BEFORE = '''
def _add_tools_tree_args(p):
    m = p.add_mutually_exclusive_group(required=True)
    m.add_argument("--prune", action="store_true", help="drop dead lineages")
    m.add_argument("--red", action="store_true", help="rescale to RED")
    p.add_argument("-o", "--output", metavar="FILE", help="write here instead of stdout")
'''
_CLI_AFTER = _CLI_BEFORE.replace(
    '    p.add_argument("-o"',
    '    m.add_argument("--gamma", action="store_true", help="Pybus & Harvey\'s gamma")\n'
    '    p.add_argument("-o"')

_TREE_BEFORE = '''
def prune(tree): ...
def red_scaled(tree): ...
__all__ = ["prune", "red_scaled"]
'''
_TREE_AFTER = _TREE_BEFORE.replace(
    '__all__ = ["prune", "red_scaled"]',
    'def gamma_statistic(tree): ...\n__all__ = ["prune", "red_scaled", "gamma_statistic"]')


def _surface(cli: str, tree: str) -> set[str]:
    return public_surface({"zombi2/cli/tools.py": cli, "zombi2/tree.py": tree})


def test_it_catches_the_change_it_was_written_for():
    """#319: a new flag and a new exported function, which shipped with no changelog line."""
    before = _surface(_CLI_BEFORE, _TREE_BEFORE)
    after = _surface(_CLI_AFTER, _TREE_AFTER)
    assert after - before == {"flag --gamma", "api zombi2.tree.gamma_statistic"}
    assert not before - after


def test_rewriting_help_text_and_docstrings_raises_nothing():
    """The other half of #319 rewrote help strings across six CLI modules. A guard that read diff
    lines would have fired on every one of them, which is how a check gets ignored."""
    reworded = (_CLI_BEFORE
                .replace('help="drop dead lineages"', 'help="drop dead and unsampled lineages"')
                .replace('help="rescale to RED"',
                         'help="rescale to RED, Relative Evolutionary Divergence on [0,1]"'))
    documented = '"""Trees, and what you can measure on them."""\n' + _TREE_BEFORE
    assert _surface(reworded, documented) == _surface(_CLI_BEFORE, _TREE_BEFORE)


def test_reformatting_and_reordering_raise_nothing():
    """A surface is a set: reordering `__all__` or splitting a call over lines is not a change."""
    reordered = _TREE_BEFORE.replace('["prune", "red_scaled"]', '[\n    "red_scaled",\n    "prune",\n]')
    rewrapped = _CLI_BEFORE.replace(
        '    m.add_argument("--red", action="store_true", help="rescale to RED")',
        '    m.add_argument(\n        "--red",\n        action="store_true",\n'
        '        help="rescale to RED",\n    )')
    assert _surface(rewrapped, reordered) == _surface(_CLI_BEFORE, _TREE_BEFORE)


def test_a_removed_flag_counts_too():
    """Taking a flag away is at least as observable as adding one."""
    before = _surface(_CLI_BEFORE, _TREE_BEFORE)
    after = _surface(_CLI_BEFORE.replace(
        '    m.add_argument("--red", action="store_true", help="rescale to RED")\n', ""),
        _TREE_BEFORE)
    assert before - after == {"flag --red"}


def test_it_ignores_everything_outside_the_package():
    """The gallery, the analyses and the tests all call argparse and none of them is a promise."""
    assert public_surface({"gallery/build.py": _CLI_AFTER, "tests/test_x.py": _TREE_AFTER}) == set()


def test_flags_are_only_read_from_the_cli_package():
    """`__all__` is public anywhere in zombi2; an `add_argument` outside `zombi2/cli/` is not a
    command-line promise, so only the exports of such a module count."""
    surface = public_surface({"zombi2/tools/homology.py": _CLI_AFTER + _TREE_AFTER})
    assert not any(label.startswith("flag ") for label in surface)
    assert "api zombi2.tools.homology.gamma_statistic" in surface


def test_a_computed_all_is_skipped_rather_than_guessed_at():
    """It would rather miss a change than invent one — a false alarm is what gets a check disabled."""
    assert public_surface({"zombi2/x.py": "__all__ = [n for n in dir() if not n.startswith('_')]"}) \
        == set()


def test_an_unparseable_revision_raises_nothing():
    assert public_surface({"zombi2/x.py": "def broken(:\n"}) == set()


def test_the_real_package_has_a_surface_this_finds():
    """A guard on a pattern the codebase does not actually use would be silently vacuous."""
    import subprocess

    root = pathlib.Path(__file__).resolve().parent.parent
    paths = subprocess.run(["git", "ls-files", "zombi2/*.py"], cwd=root, check=True,
                           capture_output=True, text=True).stdout.split()
    sources = {p: (root / p).read_text(encoding="utf-8") for p in paths}
    surface = public_surface(sources)
    assert "api zombi2.tree.gamma_statistic" in surface, "the exports scan found nothing real"
    assert "flag --birth" in surface, "the flag scan found nothing real"
