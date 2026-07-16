"""``read_newick`` robustness: malformed input is reported, never crashed on (2026-07-16 audit).

Every simulate command loads its species tree through :func:`zombi2.read_newick`, and the CLI's
error handler only catches ``(ValueError, RuntimeError, FileNotFoundError, OSError)`` — so a parse
failure that escapes as anything else reaches the user as a raw traceback. The parser used to index
the string blindly after each clade: a truncated tree ran off the end with ``IndexError``, and a
character that was neither ``,`` nor ``)`` fell through *both* branches and silently mis-parsed.
The empty-string case already raised a clean ``ValueError``; these tests pin the rest to match.
"""

import pytest

from zombi2.tree import read_newick


@pytest.mark.parametrize("bad, why", [
    ("(a,b", "truncated: never closes the clade"),
    ("((a,b)", "truncated: outer clade left open"),
    ("(a,(b,c)", "truncated after a nested clade"),
    ("(a,b(c,d))", "internal label before its clade — used to parse as siblings"),
    ("(a:xx,b:1)", "branch length is not a number"),
    ("", "empty string"),
])
def test_malformed_newick_raises_valueerror(bad, why):
    """A malformed tree is a ``ValueError`` with a message, not an ``IndexError`` traceback."""
    with pytest.raises(ValueError) as exc:
        read_newick(bad)
    assert str(exc.value), f"empty error message for {why}"


def test_malformed_newick_never_raises_indexerror():
    """Specifically: the failure mode the CLI could not catch. ``IndexError`` is not a ``ValueError``."""
    for bad in ("(a,b", "((a,b)", "(a,(b,c)"):
        with pytest.raises(ValueError):      # would be IndexError before the fix
            read_newick(bad)


@pytest.mark.parametrize("good", [
    "((a:1,b:1)i1:1,c:2);",
    "((a:1,b:1)i1:1,c:2)",                    # trailing semicolon is optional
    "  ( (a:1, b:1) i1:1 , c:2 ) ;  ",        # insignificant whitespace
    "('Homo sapiens':1,b:1);",                # quoted label with a space
    "(a,b);",                                 # no branch lengths
])
def test_valid_newick_still_parses(good):
    """The stricter parser must not reject anything legal."""
    tree = read_newick(good)
    assert len(tree.leaves()) >= 2
    assert tree.root is not None


def test_quoted_label_with_structural_char_is_not_mistaken_for_syntax():
    """A ',' or ')' inside a quoted label is data, not structure — the strictness must not break it."""
    tree = read_newick("('a,weird)name':1,b:1);")
    assert {leaf.name for leaf in tree.leaves()} == {"a,weird)name", "b"}
