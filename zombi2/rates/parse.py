"""The **written form** of a rate (SPEC §5): text in, a rate spec out.

``scope(base) × modifiers`` is written the same way everywhere — in Python, on the command line, and
in a ``--params`` file. This module is what makes the last two true: it reads the expression a user
would type in Python and returns the same object, so a snippet pastes between the three unchanged::

    parse_rate("1.0")                              -> 1.0
    parse_rate("Global(1.0)")                      -> scope.Global(1.0)
    parse_rate("1.0 * OnTime({0: 1.0, 3: 0.3})")   -> a Rate carrying that modifier

The ``mod.`` / ``scope.`` qualifiers Python needs are optional here, so
``1.0 * mod.FromParent(spread=0.2)`` and ``1.0 * FromParent(spread=0.2)`` both read.

**It parses, it does not evaluate.** The text is parsed to a syntax tree and walked against a
whitelist — the scope wrappers, the modifiers, numbers, strings, dicts/lists, keyword arguments, and
``*``. There is no ``eval``, no builtins, no attribute access beyond the two optional qualifiers, so
a parameters file from a colleague cannot run code.

Only ``*`` composes, because that is the only operation the grammar defines: a rate is ``time⁻¹`` and
a modifier dimensionless, so ``+`` and ``/`` between them mean nothing (SPEC §5). Whether a given
modifier is *supported* is not this module's business — each level declares what it wires and rejects
the rest, with a message naming the alternatives.
"""

from __future__ import annotations

import ast
import difflib

from . import mapping as _mapping
from . import modifiers as _modifiers
from . import scope as _scope
from .rate import Rate

#: the names an expression may call — the scope wrappers, the modifiers, and the driver responses.
#: The abstract bases (``Scope``, ``Modifier``, ``Mapping``) are deliberately absent: they are not
#: things a user writes. ``Curve`` needs a callable, which this grammar cannot express, so it is
#: excluded here and reported with a pointer to the Python API.
_NAMES: dict[str, type] = {
    **{n: getattr(_scope, n) for n in _scope.__all__ if n != "Scope"},
    **{n: getattr(_modifiers, n) for n in _modifiers.__all__ if n != "Modifier"},
    "Table": _mapping.Table,
    "Scalar": _mapping.Scalar,
}

#: the optional Python qualifiers — ``mod.OnTime(...)`` / ``scope.Global(...)`` read as themselves
_QUALIFIERS = frozenset({"mod", "modifiers", "scope", "scopes"})

_OP_NAMES = {ast.Add: "+", ast.Sub: "-", ast.Div: "/", ast.FloorDiv: "//", ast.Mod: "%",
             ast.Pow: "**", ast.MatMult: "@"}


class RateSyntaxError(ValueError):
    """A rate expression that could not be read. A ``ValueError``, so every caller that already
    reports the rate classes' own domain errors reports this the same way."""


def _fail(message: str, text: str) -> RateSyntaxError:
    return RateSyntaxError(f"{message}\n  in the rate {text!r}")


def _unknown_name(name: str, text: str) -> RateSyntaxError:
    close = difflib.get_close_matches(name, _NAMES, n=1, cutoff=0.6)
    hint = f" — did you mean {close[0]!r}?" if close else ""
    scopes = ", ".join(n for n in _NAMES if n in _scope.__all__)
    mods = ", ".join(n for n in _NAMES if n in _modifiers.__all__)
    return _fail(
        f"unknown name {name!r}{hint}\n"
        f"  scopes:    {scopes}\n"
        f"  modifiers: {mods}", text)


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
                f"only '*' composes a rate, got {op!r} — a rate is scope(base) × modifiers, and a "
                f"modifier is a dimensionless multiplier (SPEC §5)", text)
        left, right = _node(node.left, text), _node(node.right, text)
        try:
            return left * right
        except TypeError:                         # e.g. a list or a string on one side of the '*'
            raise _fail(
                f"cannot compose {type(left).__name__} with {type(right).__name__} — '*' puts a "
                f"modifier on a base or a scope", text) from None

    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        # parsed so that a negative number is rejected by the rate itself ("must be non-negative"),
        # which says far more than a syntax error would
        value = _node(node.operand, text)
        if not isinstance(value, (int, float)):
            raise _fail("a sign may only be applied to a number", text)
        return -value if isinstance(node.op, ast.USub) else value

    if isinstance(node, ast.Call):
        cls = _callable(node.func, text)
        args = [_node(a, text) for a in node.args]
        kwargs = {}
        for kw in node.keywords:
            if kw.arg is None:
                raise _fail("'**' unpacking is not allowed in a rate", text)
            kwargs[kw.arg] = _node(kw.value, text)
        try:
            return cls(*args, **kwargs)
        except TypeError as e:                    # wrong arity / unknown keyword
            raise _fail(f"{cls.__name__}: {e}", text) from None

    if isinstance(node, ast.Dict):
        if any(k is None for k in node.keys):     # {**other}
            raise _fail("'**' unpacking is not allowed in a rate", text)
        return {_node(k, text): _node(v, text) for k, v in zip(node.keys, node.values)}

    if isinstance(node, (ast.List, ast.Tuple)):
        values = [_node(e, text) for e in node.elts]
        return values if isinstance(node, ast.List) else tuple(values)

    if isinstance(node, ast.Name):
        if node.id in _NAMES:
            raise _fail(f"{node.id!r} is used as a value, but a scope or modifier is built by "
                        f"calling it — write {node.id}(...)", text)
        raise _unknown_name(node.id, text)

    raise _fail(f"{type(node).__name__} is not allowed in a rate expression", text)


def _callable(func: ast.AST, text: str) -> type:
    """Resolve the callable of a ``Call``: a whitelisted name, optionally qualified ``mod.X``."""
    if isinstance(func, ast.Attribute):
        if not isinstance(func.value, ast.Name) or func.value.id not in _QUALIFIERS:
            raise _fail(
                "a rate may only call a scope or a modifier by name (optionally qualified "
                "'mod.' / 'scope.')", text)
        name = func.attr
    elif isinstance(func, ast.Name):
        name = func.id
    else:
        raise _fail("a rate may only call a scope or a modifier by name", text)

    if name == "Curve":
        raise _fail(
            "Curve maps a driver with a function, which cannot be written on the command line — "
            "use the Python API for a continuous response, or a Table for a discrete one", text)
    if name not in _NAMES:
        raise _unknown_name(name, text)
    return _NAMES[name]


def parse_rate(text: object):
    """Read a rate in its written form and return the spec the ``simulate_*`` functions take.

    ``text`` is the expression (``"1.0"``, ``"Global(1.0)"``,
    ``"1.0 * OnTime({0: 1.0, 3: 0.3})"``); a number passes through, so a ``--params`` value that is
    already a TOML float needs no special case. The result is a number, a scope wrapper, a modifier,
    or a ``Rate`` — all four are what ``as_rate`` accepts, so every level takes it as it is.

    Raises :class:`RateSyntaxError` (a ``ValueError``) for anything outside the grammar, and lets the
    scope/modifier classes raise their own domain errors (a negative base, an empty schedule, …).
    """
    if isinstance(text, bool):
        raise RateSyntaxError(f"a rate must be a number or an expression, got {text!r}")
    if isinstance(text, (int, float)):
        return float(text)
    if not isinstance(text, str):
        raise RateSyntaxError(f"a rate must be a number or an expression, got {text!r}")
    if not text.strip():
        raise RateSyntaxError("a rate cannot be empty")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as e:
        raise _fail(f"could not read the expression ({e.msg})", text) from None
    value = _node(tree.body, text)
    if isinstance(value, str):
        raise _fail("a rate is a number, not text", text)
    if isinstance(value, int):        # "1" is the rate 1.0, not the integer 1
        return float(value)
    return value


def written_form(spec: object) -> str:
    """Render a rate spec back as the expression that produced it — the inverse of :func:`parse_rate`.

    Used where a run records what it was given (the ``*.log`` every command writes), so the record is
    something you can paste straight back into a flag or a ``--params`` file rather than a repr you
    would have to translate. Anything that is not a rate spec is returned as its ``repr``.
    """
    # bases are rendered with repr(float), not a fixed precision: this is a reproducibility record,
    # so 0.123456789 must come back as itself rather than rounded to six significant digits
    if isinstance(spec, bool):
        return repr(spec)
    if isinstance(spec, (int, float)):
        return repr(float(spec))
    if isinstance(spec, _scope.Scope):
        return f"{type(spec).__name__}({float(spec.base)!r})"
    if isinstance(spec, _modifiers.Modifier):
        return f"1.0 * {spec!r}"
    if isinstance(spec, Rate):
        head = (f"{type(spec.scope).__name__}({float(spec.base)!r})" if spec.scope is not None
                else repr(float(spec.base)))
        return " * ".join([head, *(repr(m) for m in spec.modifiers)])
    return repr(spec)


__all__ = ["parse_rate", "written_form", "RateSyntaxError"]
