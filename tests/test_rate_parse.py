"""Tests for zombi2.params.parse — the written form of a rate (SPEC §5).

The point of this module is that there is exactly *one* way to write a rate: what you type in Python
is what you type on the command line and in a ``--params`` file. So the tests are mostly "this text
produces the object the Python expression produces", plus the guarantees that make it safe to accept
that text from a file: it parses, it never evaluates, and the only attributes it will follow are the
verbs.
"""

import pytest

from zombi2.params import (Between, Drift, Extent, Gamma, LogNormal, PerCopy, PerLineage,
                          Random, Recipients, Scalar, TotalDiversity)
from zombi2.params import driver as drv
from zombi2.params import law as law
from zombi2.params import scope
from zombi2.params.rate import Rate
from zombi2.params.parse import RateSyntaxError, parse_rate, written_form


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


def test_a_scope_reads_as_the_rate_it_builds():
    assert parse_rate("Global(1.0)") == scope.Global(1.0)
    assert parse_rate("PerCopy(0.25)") == scope.PerCopy(0.25)
    assert parse_rate("PerCopy(0.25)").scope is scope.PerCopy


def test_changing_at_matches_the_python_expression():
    assert parse_rate("PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})") \
        == PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})


def test_keyword_arguments():
    assert parse_rate("PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2), bins=8))") \
        == PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2), bins=8))
    assert parse_rate("PerLineage(1.0).scaled_by(TotalDiversity(cap=100))") \
        == PerLineage(1.0).scaled_by(TotalDiversity(cap=100))


def test_a_string_argument():
    assert parse_rate("PerLineage(1.0).varying_among('lineages', LogNormal(0.0, 0.3))") \
        == PerLineage(1.0).varying_among('lineages', LogNormal(0.0, 0.3))


def test_a_distribution_argument():
    """A distribution is built from literals, so it is writable and round-trips — which the
    one-written-form rule needs, now that a law is written as an object rather than a name."""
    assert parse_rate("PerCopy(0.25).varying_among('families', Gamma(shape=4.0, scale=0.25))") \
        == PerCopy(0.25).varying_among('families', Gamma(4.0, 0.25))


def test_verbs_chain():
    r = parse_rate("PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))"
                   ".scaled_by(TotalDiversity(cap=100))")
    assert isinstance(r, Rate)
    assert r.modifiers == (Random('lineages', Drift(LogNormal(0.0, 0.2))),
                           drv.OnTotalDiversity(cap=100))


def test_a_scope_and_a_verb_compose():
    r = parse_rate("Global(1.0).changing_at({0: 1.0, 3: 0.3})")
    assert r.scope is scope.Global and r.modifiers == (drv.OnTime({0: 1.0, 3: 0.3}),)


def test_the_python_qualifiers_are_optional():
    # so a snippet copied out of the manual pastes into a shell unchanged
    assert parse_rate("scope.Global(1.0)") == parse_rate("Global(1.0)")
    assert parse_rate("scope.PerCopy(0.25).changing_at({0: 1.0})") \
        == parse_rate("PerCopy(0.25).changing_at({0: 1.0})")


def test_a_driver_reads_as_a_scaled_rate():
    r = parse_rate("PerCopy(0.25).scaled_by('habitat.tsv', {'aquatic': 3.0, 'terrestrial': 1.0})")
    assert r == PerCopy(0.25).scaled_by("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})


def test_a_replaced_base_reads_from_the_bare_scope():
    assert parse_rate("PerCopy().set_by('habitat.tsv', {'aquatic': 1.0})") \
        == PerCopy().set_by("habitat.tsv", {"aquatic": 1.0})


def test_the_other_two_parameter_kinds_read_too():
    """An extent and a choice are written from their own entry points, and both go through the same
    whitelist — one written form covers all three parameters, not only rates."""
    assert parse_rate("Extent(500).scaled_by('habitat.tsv', {'aquatic': 2.0})") \
        == Extent(500).scaled_by("habitat.tsv", {"aquatic": 2.0})
    assert parse_rate("Recipients().weighted_by('competence.tsv', {'competent': 3.0})") \
        == Recipients().weighted_by("competence.tsv", {"competent": 3.0})


def test_a_between_kernel_reads_as_a_choices_weight():
    r = parse_rate("Recipients().weighted_by('habitat.tsv', "
                   "Between({('marine', 'soil'): 3.0, ('soil', 'marine'): 3.0}))")
    assert r == Recipients().weighted_by(
        "habitat.tsv", Between({("marine", "soil"): 3.0, ("soil", "marine"): 3.0}))


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


def test_only_a_name_from_the_grammar_may_be_called():
    with pytest.raises(RateSyntaxError, match="unknown name 'os'"):
        parse_rate("os.system('x')")


def test_an_attribute_that_is_not_a_verb_never_reaches_getattr():
    """The whitelist that replaced "only ``*`` composes". An expression is a call on an attribute
    now, so checking the names being called is only half the check."""
    with pytest.raises(RateSyntaxError, match="'upper' is not a verb"):
        parse_rate("PerCopy(0.25).upper()")
    with pytest.raises(RateSyntaxError, match=r"'__class__' is not a verb"):
        parse_rate("PerCopy(0.25).__class__()")


def test_a_misspelt_verb_suggests_the_real_one():
    with pytest.raises(RateSyntaxError, match="did you mean 'scaled_by'"):
        parse_rate("PerCopy(0.25).scald_by('h', {'a': 1.0})")


def test_a_verb_only_goes_on_a_parameter():
    with pytest.raises(RateSyntaxError, match="a verb goes on a parameter"):
        parse_rate("[1, 2].scaled_by('h', {'a': 1.0})")


def test_a_parameter_that_refuses_a_verb_says_which_one_fits():
    with pytest.raises(RateSyntaxError, match="the verb is scaled_by"):
        parse_rate("PerCopy(0.25).weighted_by('h', {'a': 1.0})")


# --- only a verb composes -------------------------------------------------

@pytest.mark.parametrize("text,op", [("1.0 + 2.0", r"\+"), ("1.0 - 0.5", "-"),
                                     ("1.0 / 2.0", "/"), ("2.0 ** 2", r"\*\*")])
def test_other_operators_are_rejected(text, op):
    with pytest.raises(RateSyntaxError, match=f"only a verb composes a rate, got '{op}'"):
        parse_rate(text)


def test_star_is_still_read_far_enough_to_name_the_replacement():
    """An old ``--params`` file starts by multiplying. If ``*`` stopped parsing entirely, the file
    would fail with a syntax error instead of the retired-name sentence that says what to write."""
    with pytest.raises(RateSyntaxError, match="a random value has its own verb"):
        parse_rate("0.25 * Drawn(per='family', dist=LogNormal(0.0, 0.5))")
    with pytest.raises(RateSyntaxError, match="the run's clock has its own verb"):
        parse_rate("1.0 * OnTime({0: 1.0, 3: 0.3})")
    # only when both sides are still readable does '*' itself get the blame
    with pytest.raises(RateSyntaxError, match=r"'\*' no longer composes a rate"):
        parse_rate("PerCopy(1.0) * PerCopy(2.0)")


@pytest.mark.parametrize("text, names", [
    ("ScaledBy('h', {'a': 1.0})", "the verbs are methods on the parameter now"),
    ("SetBy('h', {'a': 1.0})", "written from the bare scope"),
    ("Weights('h', {'a': 1.0})", "Recipients\\(\\).weighted_by"),
    ("OnTotalDiversity(cap=100)", r"scaled_by\(TotalDiversity\(cap=100\)\)"),
    ("Inherited(per='lineage', dist=LogNormal(0.0, 0.2))", r"Drift\(LogNormal"),
    ("ByFamily(0.5)", r"varying_among\('families'"),
    ("ByLineage(0.5)", r"varying_among\('lineages'"),
    ("FromParent(0.5)", r"Drift\(LogNormal"),
    ("DrivenBy('h', {'a': 1.0})", "write the verb that says what the number does"),
])
def test_a_retired_name_says_what_replaced_it(text, names):
    """A difflib guess is no help here: the replacement is a **verb** on the parameter rather than
    another name, so there is nothing close enough to guess at."""
    with pytest.raises(RateSyntaxError, match=names):
        parse_rate(text)


@pytest.mark.parametrize("text, names", [
    ("PerCopy(0.25).varying_among(per='families')", "units are plural"),
    ("PerCopy(0.25).varying_among('families', spread=0.5)", "write the law out"),
    ("Random(per='family', dist=LogNormal(0.0, 0.5))", "units are plural"),
])
def test_a_retired_keyword_says_what_replaced_it(text, names):
    with pytest.raises(RateSyntaxError, match=names):
        parse_rate(text)


# --- errors a user will actually hit --------------------------------------

def test_an_unknown_driver_suggests_the_real_one():
    with pytest.raises(RateSyntaxError, match="did you mean 'TotalDiversity'"):
        parse_rate("PerLineage(1.0).scaled_by(TotalDivrsity(cap=10))")


def test_an_unknown_name_lists_the_menu():
    with pytest.raises(RateSyntaxError, match="names:  .*Random"):
        parse_rate("PerLineage(1.0).scaled_by(Wobble(3))")
    with pytest.raises(RateSyntaxError, match="verbs:  .*varying_among"):
        parse_rate("PerLineage(1.0).scaled_by(Wobble(3))")


def test_a_name_used_as_a_value_says_to_call_it():
    with pytest.raises(RateSyntaxError, match=r"write Clade\(\.\.\.\)"):
        parse_rate("PerLineage(1.0).scaled_by(Clade, {'a': 1.0})")


def test_a_misspelt_keyword_names_the_verb():
    with pytest.raises(RateSyntaxError, match="varying_among:"):
        parse_rate("PerLineage(1.0).varying_among('lineages', LogNormal(0.0, 0.3), spred=0.3)")


def test_curve_points_at_the_python_api():
    # Curve maps a driver with a callable, which no text grammar can carry
    with pytest.raises(RateSyntaxError, match="use the Python API"):
        parse_rate("PerLineage(1.0).scaled_by('h', Curve(lambda x: x))")


def test_a_value_written_on_its_own_is_not_a_rate():
    """A `Random` says how something varies without saying what varies or per what, so it needs the
    scope and the verb that read it — and the message writes both."""
    with pytest.raises(RateSyntaxError, match=r"PerCopy\(0.25\).varying_among\('families'"):
        parse_rate("Random('families', LogNormal(0.0, 0.5))")


def test_empty_and_non_text():
    with pytest.raises(RateSyntaxError, match="cannot be empty"):
        parse_rate("   ")
    with pytest.raises(RateSyntaxError, match="a rate is a number, not text"):
        parse_rate("'hello'")
    with pytest.raises(RateSyntaxError):
        parse_rate(True)


def test_the_rate_classes_still_raise_their_own_domain_errors():
    # the parser does not duplicate validation — a negative base is the rate's error, not a syntax one
    with pytest.raises(ValueError, match="non-negative"):
        parse_rate("Global(-1)")
    with pytest.raises(ValueError, match="schedule cannot be empty"):
        parse_rate("PerLineage(1.0).changing_at({})")
    with pytest.raises(ValueError, match="its own verb"):
        parse_rate("PerLineage(1.0).scaled_by(Time(), {0: 1.0})")


def test_a_syntax_error_quotes_the_expression():
    with pytest.raises(RateSyntaxError, match=r"PerLineage\(1.0\).changing_at\(\{0: 1.0"):
        parse_rate("PerLineage(1.0).changing_at({0: 1.0")


# --- written_form is the inverse -----------------------------------------

@pytest.mark.parametrize("text", [
    "1.0",
    "Global(1.0)",
    "PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})",
    "PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})",
    "PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))"
    ".scaled_by(TotalDiversity(cap=100))",
    "PerLineage(1.0).varying_among('lineages', Gamma(shape=11.11, scale=0.09))",
    "PerCopy(0.25).scaled_by('habitat.tsv', {'aquatic': 3.0})",
    "PerCopy().set_by('habitat.tsv', {'aquatic': 1.0})",
    "Extent(500).changing_at({0: 1.0, 3: 0.3})",
    "Recipients().weighted_by('competence.tsv', {'competent': 3.0})",
])
def test_written_form_round_trips(text):
    once = written_form(parse_rate(text))
    assert parse_rate(once) == parse_rate(text)      # the rendering means the same thing
    assert written_form(parse_rate(once)) == once    # and it is a fixed point


def test_written_form_keeps_full_precision():
    # the run log is a reproducibility record, so a base must not be rounded on its way in
    assert written_form(parse_rate("0.123456789")) == "0.123456789"


def test_a_windows_path_in_a_rate_is_taken_as_written():
    # the strings in a rate are paths and state labels, never escape sequences — but the expression
    # is read by Python's own parser, which sees C:\Users and reports a truncated \UXXXXXXXX escape.
    # A pasted path is the normal way to write one, so it has to mean itself.
    from zombi2.params.connection import Driven

    rate = parse_rate(r"PerCopy(0.1).scaled_by('C:\Users\me\trait_events.tsv', {'a': 2.0})")
    driver = next(m for m in rate.modifiers if isinstance(m, Driven))
    assert driver.driver == r"C:\Users\me\trait_events.tsv"

    unc = parse_rate(r"PerCopy(0.1).scaled_by('\\server\share\trait.tsv', {'a': 2.0})")
    assert next(m for m in unc.modifiers if isinstance(m, Driven)).driver == r"\\server\share\trait.tsv"

    posix = parse_rate("PerCopy(0.1).scaled_by('/home/me/trait.tsv', {'a': 2.0})")
    assert next(m for m in posix.modifiers if isinstance(m, Driven)).driver == "/home/me/trait.tsv"


def test_an_already_escaped_path_is_left_as_written():
    # repr() of a path is the natural way to build an expression in Python, and it escapes the
    # backslashes properly — reading those literally as well would double them.
    from zombi2.params.connection import Driven

    path = r"C:\Users\me\trait_events.tsv"
    rate = parse_rate(f"PerCopy(0.1).scaled_by({path!r}, {{'a': 2.0}})")
    assert next(m for m in rate.modifiers if isinstance(m, Driven)).driver == path


def test_a_path_whose_every_backslash_is_a_valid_escape_still_means_itself():
    # the dangerous one: \t \n \f are all real escapes, so this PARSES and silently becomes control
    # characters. A path never contains one, which is how it is caught.
    from zombi2.params.connection import Driven

    rate = parse_rate(r"PerCopy(0.1).scaled_by('C:\temp\new\file.tsv', {'a': 2.0})")
    assert next(m for m in rate.modifiers if isinstance(m, Driven)).driver == r"C:\temp\new\file.tsv"


def test_a_rate_is_written_to_full_precision():
    """The written form is a reproducibility record: a run's log holds it, and a reader pastes it
    back into a flag. Every mapping printed its numbers with `:g`, which rounds to six significant
    figures, so a factor of 0.0123456789012 was logged as 0.0123457 — the same run recorded as a
    different model, with nothing to say so. Only the base was exact.
    """
    exact = 0.0123456789012
    for rate in (PerLineage(1.0).changing_at({0: 1.0, 3: exact}),
                 PerLineage(1.0).scaled_by("h.tsv", {"cave": exact}),
                 PerLineage(1.0).scaled_by("h.tsv", {"cave": 1.0}),
                 PerLineage(1.0).scaled_by("h.tsv", Scalar(exact)),
                 PerLineage(1.0).scaled_by("h.tsv", Between({("a", "b"): exact}, default=exact))):
        text = written_form(rate)
        assert parse_rate(text) == rate, text

    # and the digits really are all there, not merely equal after a lucky round
    assert repr(exact) in written_form(PerLineage(1.0).changing_at({0: 1.0, 3: exact}))


def test_the_written_form_keeps_the_grammars_own_message_for_a_misplaced_set_by():
    """`Rate.set_by` raises with a sentence written for exactly this mistake. The parser caught
    every `TypeError` from a call and replaced it with the class's name, so
    ``--loss "PerCopy(0.25).set_by(...)"`` lost the sentence that says what to write instead. Ours
    are kept; CPython's own errors, which are about types rather than the grammar, are not."""
    # matched on the GUIDANCE, not on a class name: the generic message names the call too, so
    # "set_by" alone would pass whether the sentence survived or not
    for text in ("PerCopy(0.25).set_by('h', {'c': 1.0})",
                 "PerCopy(0.25).scaled_by('s', {'b': 1.0}).set_by('h', {'c': 1.0})",
                 "PerCopy().set_by('h', {'c': 1.0}).set_by('h', {'c': 2.0})"):
        with pytest.raises(RateSyntaxError, match="silently discard"):
            parse_rate(text)


def test_a_call_cpython_refuses_gets_the_parsers_own_wording():
    """The other half of the rule above. A whitelisted name can also be called wrongly — a missing
    argument, an unknown keyword — and CPython's message is about the signature rather than about
    the rate, so the parser answers it with the name of what was being called."""
    for text, name in (("PerLineage(1.0).changing_at()", "changing_at:"),
                       ("Clade()", "Clade:"),
                       ("LogNormal(0.0, 0.3, 9.9)", "LogNormal:")):
        with pytest.raises(RateSyntaxError, match=name):
            parse_rate(text)
