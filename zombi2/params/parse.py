"""The **written form** of a parameter (SPEC §5): text in, a parameter spec out.

A rate is written the same way everywhere — in Python, on the command line, and in a ``--params``
file. This module is what makes the last two true: it reads the expression a user would type in
Python and returns the same object, so a snippet pastes between the three unchanged::

    parse_rate("1.0")                                     -> 1.0
    parse_rate("Global(1.0)")                             -> a Rate scoped to the whole run
    parse_rate("PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})")   -> a Rate carrying that schedule

The ``scope.`` qualifier Python needs is optional here, so ``scope.PerCopy(0.25)`` and
``PerCopy(0.25)`` both read.

**It parses, it does not evaluate.** The text is parsed to a syntax tree and walked against two
whitelists — the names it may **call** (the scopes, the drivers, the laws, the distributions, the
mappings) and the verbs it may follow as an **attribute**. Numbers, strings, dicts, lists and
keyword arguments come through; nothing else does. There is no ``eval``, no builtins, and no
attribute access beyond those verbs and the optional qualifier, so a parameters file from a
colleague cannot run code.

The verb whitelist is the shape this grammar needs and the older one did not: an expression is now a
call **on an attribute** (``PerCopy(0.25).scaled_by(...)``) rather than a product, so checking the
names being called is only half the check. Whether a given verb is *supported* by a level is not
this module's business — each level declares what it takes and rejects the rest, with a message
naming the alternatives.
"""

from __future__ import annotations

import ast
import warnings
import difflib
from typing import Any, cast

from . import choice as _choice
from . import extent as _extent
from . import mapping as _mapping
from . import distributions as _distributions
from . import modifiers as _modifiers
from . import driver as _driver
from . import verbs as _verbs
from . import scope as _scope
from .rate import Rate, RateCompositionError
from .retired import RETIRED, RETIRED_KEYWORDS, keyword_message, name_message

#: the names an expression may **call** — the scopes, the drivers and laws, the distributions, the
#: mappings, and the two entry points the other parameter kinds are written from. The abstract bases
#: (``Scope``, ``Modifier``, ``Mapping``) are deliberately absent: they are not things a user writes.
#: ``Curve`` needs a callable, which this grammar cannot express, so it is excluded here and
#: reported with a pointer to the Python API.
_NAMES: dict[str, Any] = {
    **{n: getattr(_scope, n) for n in _scope.__all__ if n != "Scope"},
    **{n: getattr(_modifiers, n) for n in _modifiers.WRITABLE},
    **{n: getattr(_driver, n) for n in _driver.WRITABLE},
    **{n: getattr(_distributions, n) for n in _distributions.WRITABLE},
    **{n: getattr(_choice, n) for n in _choice.WRITABLE},
    "Extent": _extent.Extent,
    "Table": _mapping.Table,
    "Scalar": _mapping.Scalar,
    "Between": _mapping.Between,  # the choice's kernel: weighted_by(driver, Between({(a, b): w}))
}

#: what a verb may be written on, and what to call each in a message
_PARAMETERS: dict[type, str] = {
    Rate: "a rate", _extent.Extent: "an extent", _choice.Choice: "a choice",
}

#: the optional Python qualifiers — ``scope.Global(...)`` reads as itself
_QUALIFIERS = frozenset({"mod", "modifiers", "scope", "scopes"})

_OP_NAMES = {ast.Add: "+", ast.Sub: "-", ast.Div: "/", ast.FloorDiv: "//", ast.Mod: "%",
             ast.Pow: "**", ast.MatMult: "@"}


class RateSyntaxError(ValueError):
    """A rate expression that could not be read. A ``ValueError``, so every caller that already
    reports the rate classes' own domain errors reports this the same way."""


def _paths_are_literal(text: str) -> str:
    """Double every backslash inside a quoted string, so a path in a rate expression means itself.

    Only the inside of a literal is touched, so the expression's own syntax is untouched. See
    `parse_rate()` for when this is applied — never unconditionally, because an expression whose
    backslashes are *already* escaped (what ``repr()`` of a path gives) must be left as written."""
    out, quote = [], None
    for ch in text:
        if quote is None:
            quote = ch if ch in "\"'" else None
        elif ch == "\\":
            out.append("\\")                 # keep it as the character it is, not as an escape
        elif ch == quote:
            quote = None
        out.append(ch)
    return "".join(out)


def _fail(message: str, text: str) -> RateSyntaxError:
    return RateSyntaxError(f"{message}\n  in the rate {text!r}")


def _unknown_name(name: str, text: str) -> RateSyntaxError:
    if name in RETIRED:
        return _fail(name_message(name), text)
    close = difflib.get_close_matches(name, _NAMES, n=1, cutoff=0.6)
    hint = f" — did you mean {close[0]!r}?" if close else ""
    scopes = ", ".join(n for n in _NAMES if n in _scope.__all__)
    others = ", ".join(n for n in _NAMES if n not in _scope.__all__)
    return _fail(
        f"unknown name {name!r}{hint}\n"
        f"  scopes: {scopes}\n"
        f"  names:  {others}\n"
        f"  verbs:  {', '.join(_verbs.VERBS)}", text)


def _node(node: ast.AST, text: str):
    """Evaluate one whitelisted node. Anything outside the grammar raises."""
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or node.value is None:
            raise _fail(f"{node.value!r} is not a rate value", text)
        if isinstance(node.value, (int, float, str)):
            return node.value
        raise _fail(f"{node.value!r} is not allowed in a rate", text)

    if isinstance(node, ast.BinOp):
        if not isinstance(node.op, ast.Mult):
            op = _OP_NAMES.get(type(node.op), type(node.op).__name__)
            raise _fail(
                f"only a verb composes a rate, got {op!r} — a rate is a scope with verbs chained "
                f"onto it, and each verb reads one thing (SPEC §4)", text)
        # `*` composed a rate until the verbs replaced it, and it is still recognised this far so
        # that an old expression reaches the name inside it: `0.25 * Drawn(per='family', ...)` is
        # answered by the retired-name entry for `Drawn`, which says what to write, rather than by a
        # sentence about `*`. Only when both sides are still readable does `*` itself get the blame.
        _node(node.left, text)
        _node(node.right, text)
        raise _fail(
            "'*' no longer composes a rate — the verbs do, and each returns a new rate so they "
            "chain: PerCopy(0.25).scaled_by(driver, mapping), "
            ".varying_among('families', LogNormal(0.0, 0.5)), .changing_at({0: 1.0, 3: 0.3})", text)

    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        # parsed so that a negative number is rejected by the rate itself ("must be non-negative"),
        # which says far more than a syntax error would
        value = _node(node.operand, text)
        if not isinstance(value, (int, float)):
            raise _fail("a sign may only be applied to a number", text)
        return -value if isinstance(node.op, ast.USub) else value

    if isinstance(node, ast.Call):
        # the callable first, so that a retired *name* is answered by its own entry rather than by
        # the retired keyword it was usually written with: `Inherited(per=…)` is one mistake, and
        # the sentence that fixes it names the whole replacement, not just the argument
        fn = (_verb(cast(ast.Attribute, node.func), text) if _is_verb(node.func)
              else _callable(node.func, text))
        args = [_node(a, text) for a in node.args]
        kwargs = {}
        for kw in node.keywords:
            if kw.arg is None:
                raise _fail("'**' unpacking is not allowed in a rate", text)
            if kw.arg in RETIRED_KEYWORDS:
                raise _fail(keyword_message(kw.arg), text)
            kwargs[kw.arg] = _node(kw.value, text)
        return _call(fn, args, kwargs, text)

    if isinstance(node, ast.Dict):
        keys = node.keys
        if any(k is None for k in keys):          # {**other}
            raise _fail("'**' unpacking is not allowed in a rate", text)
        return {_node(k, text): _node(v, text)
                for k, v in zip(cast("list[ast.expr]", keys), node.values)}

    if isinstance(node, (ast.List, ast.Tuple)):
        values = [_node(e, text) for e in node.elts]
        return values if isinstance(node, ast.List) else tuple(values)

    if isinstance(node, ast.Name):
        if node.id in _NAMES:
            raise _fail(f"{node.id!r} is used as a value, but it is built by calling it — write "
                        f"{node.id}(...)", text)
        raise _unknown_name(node.id, text)

    raise _fail(f"{type(node).__name__} is not allowed in a rate expression", text)


def _call(fn, args: list, kwargs: dict, text: str):
    """Apply one whitelisted callable, reporting a wrong arity or an unknown keyword as a rate
    error rather than as a traceback from inside the class."""
    try:
        return fn(*args, **kwargs)
    except TypeError as e:
        if isinstance(e, RateCompositionError):
            raise _fail(str(e), text) from None
        name = getattr(fn, "__name__", type(fn).__name__)
        raise _fail(f"{name}: {e}", text) from None


def _is_verb(func: ast.AST) -> bool:
    """Whether this call is a verb on a parameter (``PerCopy(0.25).scaled_by(...)``) rather than a
    name being called, optionally qualified (``scope.PerCopy(...)``)."""
    return (isinstance(func, ast.Attribute)
            and not (isinstance(func.value, ast.Name) and func.value.id in _QUALIFIERS))


def _verb(func: ast.Attribute, text: str):
    """Resolve a verb call: the parameter to its left, then the verb name against `verbs.VERBS`.

    This is the whitelist that replaced "only ``*`` composes". A name outside the verb list never
    reaches ``getattr``, and a verb is looked up only on the three parameter classes, so no
    expression can walk out of the grammar through an attribute.
    """
    target = _node(func.value, text)
    name = func.attr
    if name not in _verbs.VERBS:
        close = difflib.get_close_matches(name, _verbs.VERBS, n=1, cutoff=0.6)
        hint = f" — did you mean {close[0]!r}?" if close else ""
        raise _fail(f"{name!r} is not a verb{hint}\n  verbs: {', '.join(_verbs.VERBS)}", text)
    kind = _PARAMETERS.get(type(target))
    if kind is None:
        raise _fail(
            f"a verb goes on a parameter — a scope, Extent(...) or Recipients() — and "
            f"{type(target).__name__} is not one", text)
    verb = getattr(target, name, None)
    if verb is None:
        # not enumerated from `target`, because the verbs a parameter refuses are still methods on
        # it — they exist to carry the sentence that says which verb to write instead
        raise _fail(
            f"{kind} takes no {name!r}: a rate is scaled, replaced or varied, an extent is scaled "
            f"or varied, and a choice is weighted (SPEC §3)", text)
    return verb


def _callable(func: ast.AST, text: str):
    """Resolve the callable of a ``Call``: a whitelisted name, optionally qualified ``scope.X``."""
    if isinstance(func, ast.Attribute):
        if not isinstance(func.value, ast.Name) or func.value.id not in _QUALIFIERS:
            raise _fail(
                "a rate may only call a name from the grammar (optionally qualified 'scope.')", text)
        name = func.attr
    elif isinstance(func, ast.Name):
        name = func.id
    else:
        raise _fail("a rate may only call a name from the grammar", text)

    if name == "Curve":
        raise _fail(
            "Curve maps a driver with a function, which cannot be written on the command line — "
            "use the Python API for a continuous mapping, or a Table for a discrete one", text)
    if name not in _NAMES:
        raise _unknown_name(name, text)
    return _NAMES[name]


def parse_rate(text: object):
    """Read a parameter in its written form and return the spec the ``simulate_*`` functions take.

    ``text`` is the expression (``"1.0"``, ``"Global(1.0)"``,
    ``"PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})"``); a number passes through, so a ``--params``
    value that is already a TOML float needs no special case. The result is a number, a ``Rate``, an
    ``Extent``, a choice, or a distribution — whichever the expression built — and each level takes
    the one it asked for.

    Raises `RateSyntaxError` (a ``ValueError``) for anything outside the grammar, and lets the
    grammar's own classes raise their domain errors (a negative base, an empty schedule, …).
    """
    if isinstance(text, bool):
        raise RateSyntaxError(f"a rate must be a number or an expression, got {text!r}")
    if isinstance(text, (int, float)):
        return float(text)
    if not isinstance(text, str):
        raise RateSyntaxError(f"a rate must be a number or an expression, got {text!r}")
    if not text.strip():
        raise RateSyntaxError("a rate cannot be empty")
    # A backslash in a path is not an escape. `scaled_by('C:\\Users\\me\\t.tsv', …)` is well formed to
    # the person who pasted it and a truncated \\UXXXXXXXX escape to Python's parser, and `C:\\temp`
    # is worse — it parses, silently, as a tab. But an expression may equally have its backslashes
    # already escaped, which is what repr() of a path gives, and that must be left alone. There is no
    # telling the two apart from the text, so: read it as written first, and fall back to reading the
    # backslashes literally when that fails, or when it succeeded and produced a control character —
    # which a file path and a state label never contain, and an eaten escape always does.
    with warnings.catch_warnings():
        # an unknown escape (\\s of a UNC \\\\server\\share) is a SyntaxWarning on the first reading and
        # nothing to tell the user about, since the literal reading below is the one that is kept
        warnings.simplefilter("ignore", SyntaxWarning)
        try:
            tree = ast.parse(text, mode="eval")
        except SyntaxError as first:
            try:
                tree = ast.parse(_paths_are_literal(text), mode="eval")
            except SyntaxError:
                raise _fail(f"could not read the expression ({first.msg})", text) from None
        else:
            if any(isinstance(n, ast.Constant) and isinstance(n.value, str)
                   and any(c < " " for c in n.value) for n in ast.walk(tree)):
                try:
                    tree = ast.parse(_paths_are_literal(text), mode="eval")
                except SyntaxError:
                    pass                        # keep the first reading; _node will report on it
    value = _node(tree.body, text)
    if isinstance(value, str):
        raise _fail("a rate is a number, not text", text)
    if isinstance(value, _modifiers.Modifier):
        # a `Random(...)` written on its own. It is a value, not a parameter: it says how something
        # varies without saying what varies or per what.
        raise _fail(
            f"{value!r} is a value, not a rate — write it on a scope, with the verb that reads it: "
            f"PerCopy(0.25).{value.written_call()}", text)
    if isinstance(value, int):        # "1" is the rate 1.0, not the integer 1
        return float(value)
    return value


def written_form(spec: object) -> str:
    """Render a parameter spec back as the expression that produced it — the inverse of
    `parse_rate()`.

    Used where a run records what it was given (the ``*.log`` every command writes), so the record is
    something you can paste straight back into a flag or a ``--params`` file rather than a repr you
    would have to translate.

    There is barely anything left to do here, and that is the point: each parameter's ``__repr__``
    **is** its written form, so one renderer serves the log, the error messages and the debugger,
    and there is no second one to drift out of step with it. What remains is the bare number, which
    is recorded as a float so that a rate given as ``1`` reads back as the rate ``1.0``.
    """
    # bases are rendered with repr(float), not a fixed precision: this is a reproducibility record,
    # so 0.123456789 must come back as itself rather than rounded to six significant digits
    if isinstance(spec, bool):
        return repr(spec)
    if isinstance(spec, (int, float)):
        return repr(float(spec))
    return repr(spec)


def written_choice(spec: object) -> str:
    """A ``transfer_to`` rule as the text that flag takes back.

    A choice is not a rate (SPEC §5), and a **named** rule is written bare — ``uniform``, not
    ``'uniform'`` — because that is what ``--transfer-to`` accepts and a quoted string is not. Every
    other rule renders itself: `Choice` writes ``Recipients().weighted_by(...)`` with no base in
    front, because a choice has none and the flag refuses one there.

    `Distance` and `Clades` render as the constructor calls the API takes.
    """
    if isinstance(spec, str):
        return spec                                   # 'uniform' / 'distance', bare
    return repr(spec)


__all__ = ["parse_rate", "written_form", "written_choice", "RateSyntaxError"]
