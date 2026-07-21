"""Tests for zombi2.rates.parse — the written form of a rate (SPEC §5).

The point of this module is that there is exactly *one* way to write a rate: what you type in Python
is what you type on the command line and in a ``--params`` file. So the tests are mostly "this text
produces the object the Python expression produces", plus the guarantees that make it safe to accept
that text from a file: it parses, it never evaluates, and only ``*`` composes.
"""

import pytest

from zombi2.rates import modifiers as mod
from zombi2.rates import scope
from zombi2.rates.rate import Rate
from zombi2.rates.parse import RateSyntaxError, parse_rate, written_form


# --- the text produces the same object as the Python expression -----------

def test_a_bare_number_is_a_rate():
    assert parse_rate("1.0") == 1.0
    assert parse_rate("0") == 0.0
    assert parse_rate("1e-3") == 0.001


def test_an_integer_reads_as_a_float():
    # "1" is the rate 1.0, so the run log and the API see a rate, not a count
    assert isinstance(parse_rate("1"), float)


def test_a_number_passes_through_unparsed():
    # a --params file may hold a TOML float; it needs no special case at the call site
    assert parse_rate(2.5) == 2.5
    assert parse_rate(3) == 3.0


def test_scope_wrapper():
    assert parse_rate("Global(1.0)") == scope.Global(1.0)
    assert parse_rate("PerCopy(0.25)") == scope.PerCopy(0.25)


def test_number_times_modifier_matches_the_python_expression():
    assert parse_rate("1.0 * OnTime({0: 1.0, 3: 0.3})") == 1.0 * mod.OnTime({0: 1.0, 3: 0.3})


def test_keyword_arguments():
    assert parse_rate("1.0 * FromParent(spread=0.2)") == 1.0 * mod.FromParent(spread=0.2)
    assert parse_rate("1.0 * OnTotalDiversity(cap=100)") == 1.0 * mod.OnTotalDiversity(cap=100)


def test_a_string_argument():
    assert parse_rate("1.0 * ByLineage(spread=0.3, dist='gamma')") == \
        1.0 * mod.ByLineage(spread=0.3, dist="gamma")


def test_modifiers_stack():
    r = parse_rate("1.0 * FromParent(spread=0.2) * OnTotalDiversity(cap=100)")
    assert isinstance(r, Rate)
    assert r.modifiers == (mod.FromParent(spread=0.2), mod.OnTotalDiversity(cap=100))


def test_a_scope_and_a_modifier_compose():
    r = parse_rate("Global(1.0) * OnTime({0: 1.0, 3: 0.3})")
    assert r.scope == scope.Global(1.0) and r.modifiers == (mod.OnTime({0: 1.0, 3: 0.3}),)


def test_the_python_qualifiers_are_optional():
    # so a snippet copied out of the manual pastes into a shell unchanged
    assert parse_rate("1.0 * mod.OnTime({0: 1.0})") == parse_rate("1.0 * OnTime({0: 1.0})")
    assert parse_rate("scope.Global(1.0)") == parse_rate("Global(1.0)")


def test_a_driver_reads_as_a_coupling():
    r = parse_rate("0.25 * DrivenBy('habitat.tsv', {'aquatic': 3.0, 'terrestrial': 1.0})")
    assert r == 0.25 * mod.DrivenBy("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})


# --- it parses; it never evaluates ---------------------------------------

@pytest.mark.parametrize("text", [
    "__import__('os').system('echo pwned')",
    "open('/etc/passwd').read()",
    "(1).__class__.__mro__",
    "[x for x in range(3)]",
    "lambda: 1",
])
def test_code_is_not_executed(text):
    # a rate can arrive from a shared --params file, so the grammar must not be a code path
    with pytest.raises(RateSyntaxError):
        parse_rate(text)


def test_only_a_scope_or_modifier_may_be_called():
    with pytest.raises(RateSyntaxError, match="only call a scope or a modifier"):
        parse_rate("os.system('x')")


# --- only '*' composes ----------------------------------------------------

@pytest.mark.parametrize("text,op", [("1.0 + 2.0", r"\+"), ("1.0 - 0.5", "-"),
                                     ("1.0 / 2.0", "/"), ("2.0 ** 2", r"\*\*")])
def test_other_operators_are_rejected(text, op):
    with pytest.raises(RateSyntaxError, match=f"only '\\*' composes a rate, got '{op}'"):
        parse_rate(text)


def test_composing_junk_is_a_readable_error():
    with pytest.raises(RateSyntaxError, match="cannot compose list with OnTime"):
        parse_rate("[1, 2] * OnTime({0: 1.0})")


# --- errors a user will actually hit --------------------------------------

def test_an_unknown_modifier_suggests_the_real_one():
    with pytest.raises(RateSyntaxError, match="did you mean 'OnTotalDiversity'"):
        parse_rate("1.0 * OnDiversity(cap=10)")


def test_an_unknown_name_lists_the_menu():
    with pytest.raises(RateSyntaxError, match="modifiers: .*OnTime"):
        parse_rate("1.0 * Wobble(3)")


def test_a_modifier_used_as_a_value_says_to_call_it():
    with pytest.raises(RateSyntaxError, match=r"write OnTime\(\.\.\.\)"):
        parse_rate("1.0 * OnTime")


def test_a_misspelt_keyword_names_the_modifier():
    with pytest.raises(RateSyntaxError, match="ByLineage:"):
        parse_rate("1.0 * ByLineage(spred=0.3)")


def test_curve_points_at_the_python_api():
    # Curve maps a driver with a callable, which no text grammar can carry
    with pytest.raises(RateSyntaxError, match="use the Python API"):
        parse_rate("1.0 * Curve(lambda x: x)")


def test_empty_and_non_text():
    with pytest.raises(RateSyntaxError, match="cannot be empty"):
        parse_rate("   ")
    with pytest.raises(RateSyntaxError, match="a rate is a number, not text"):
        parse_rate("'hello'")
    with pytest.raises(RateSyntaxError):
        parse_rate(True)


def test_the_rate_classes_still_raise_their_own_domain_errors():
    # the parser does not duplicate validation — a negative base is the scope's error, not a syntax one
    with pytest.raises(ValueError, match="non-negative"):
        parse_rate("Global(-1)")
    with pytest.raises(ValueError, match="non-empty schedule"):
        parse_rate("1.0 * OnTime({})")


def test_a_syntax_error_quotes_the_expression():
    with pytest.raises(RateSyntaxError, match="1.0 \\* OnTime\\(\\{0: 1.0"):
        parse_rate("1.0 * OnTime({0: 1.0")


# --- written_form is the inverse -----------------------------------------

@pytest.mark.parametrize("text", [
    "1.0",
    "Global(1.0)",
    "1.0 * OnTime({0: 1.0, 3: 0.3})",
    "1.0 * FromParent(spread=0.2) * OnTotalDiversity(cap=100)",
    "1.0 * ByLineage(spread=0.3, dist='gamma')",
    "0.25 * DrivenBy('habitat.tsv', {'aquatic': 3.0})",
])
def test_written_form_round_trips(text):
    once = written_form(parse_rate(text))
    assert parse_rate(once) == parse_rate(text)      # the rendering means the same thing
    assert written_form(parse_rate(once)) == once    # and it is a fixed point


def test_written_form_keeps_full_precision():
    # the run log is a reproducibility record, so a base must not be rounded on its way in
    assert written_form(parse_rate("0.123456789")) == "0.123456789"
