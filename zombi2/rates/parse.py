"""The **written form** of a rate (SPEC §5): text in, a rate spec out.

``scope(base) × modifiers`` is written the same way everywhere — in Python, on the command line, and
in a ``--params`` file. This module is what makes the last two true: it reads the expression a user
would type in Python and returns the same object, so a snippet pastes between the three unchanged::

    parse_rate("1.0")                              -> 1.0
    parse_rate("Global(1.0)")                      -> scope.Global(1.0)
    parse_rate("1.0 * OnTime({0: 1.0, 3: 0.3})")   -> a Rate carrying that modifier

The ``mod.`` / ``scope.`` qualifiers Python needs are optional here, so
``1.0 * mod.Inherited(per="lineage", spread=0.2)`` and ``1.0 * Inherited(per="lineage", spread=0.2)`` both read.

**It parses, it does not evaluate.** The text is parsed to a syntax tree and walked against a
whitelist — the scope wrappers, the modifiers, numbers, strings, dicts/lists, keyword arguments, and
``*``. There is no ``eval``, no builtins, no attribute access beyond the two optional qualifiers, so
a parameters file from a colleague cannot run code.

Only ``*`` composes, because that is the only operation the grammar defines: a rate is ``time⁻¹`` and
a modifier dimensionless, so ``+`` and ``/`` between them mean nothing (SPEC §5). Whether a given
modifier is *supported* is not this module's business — each level declares what it takes and rejects
the rest, with a message naming the alternatives.
"""

from __future__ import annotations

import ast
import warnings
import difflib
from typing import cast

from . import mapping as _mapping
from . import distributions as _distributions
from . import modifiers as _modifiers
from . import values as _values
from . import verbs as _verbs
from . import scope as _scope
from .rate import Rate, RateCompositionError

#: the names an expression may call — the scope wrappers, the modifiers, and the mappings.
#: The abstract bases (``Scope``, ``Modifier``, ``Mapping``) are deliberately absent: they are not
#: things a user writes. ``Curve`` needs a callable, which this grammar cannot express, so it is
#: excluded here and reported with a pointer to the Python API.
_NAMES: dict[str, type] = {
    **{n: getattr(_scope, n) for n in _scope.__all__ if n != "Scope"},
    **{n: getattr(_modifiers, n) for n in _modifiers.WRITABLE},
    **{n: getattr(_values, n) for n in _values.WRITABLE},
    **{n: getattr(_distributions, n) for n in _distributions.WRITABLE},
    **{n: getattr(_verbs, n) for n in _verbs.WRITABLE},
    "Table": _mapping.Table,
    "Scalar": _mapping.Scalar,
    "Between": _mapping.Between,  # the choice's kernel: Driven(driver, Between({(a, b): w}))
}

#: the optional Python qualifiers — ``mod.OnTime(...)`` / ``scope.Global(...)`` read as themselves
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


#: Names that were the written form and are not any more, each with the sentence a reader needs.
#: A difflib guess would offer one of the three verbs at random; which one is right depends on what
#: the rate is attached to, so the message says all three and what each is for.
_RETIRED = {
    "DrivenBy": ("write the verb that says what the number does: ScaledBy(driver, mapping) on a "
                 "rate or an extent, Weights(driver, mapping) on transfer_to, SetBy(driver, "
                 "mapping) to replace the base rather than scale it"),
    "ByFamily": ("write Drawn(per='family', spread=...)"),
    "ByLineage": ("write Drawn(per='lineage', spread=...)"),
    "FromParent": ("write Inherited(per='lineage', spread=...)"),
}


def _unknown_name(name: str, text: str) -> RateSyntaxError:
    if name in _RETIRED:
        return _fail(f"{name} is no longer a rate name — {_RETIRED[name]}", text)
    close = difflib.get_close_matches(name, _NAMES, n=1, cutoff=0.6)
    hint = f" — did you mean {close[0]!r}?" if close else ""
    scopes = ", ".join(n for n in _NAMES if n in _scope.__all__)
    mods = ", ".join(n for n in _NAMES
                     if n in _modifiers.WRITABLE + _values.WRITABLE + _verbs.WRITABLE
                     + _distributions.WRITABLE)
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
        except TypeError as e:
            # `SetBy.__rmul__` and `Rate.__mul__` raise TypeError with a message written for exactly
            # this mistake — "SetBy replaces the base, so there is no base to write in front of it".
            # Replacing it with the generic one below threw away the sentence that says what to do.
            # Only ours, told apart by its class rather than by the operand types: two grammar
            # objects can also fail with CPython's own "unsupported operand type(s)" — `Rate * Rate`,
            # `PerCopy * PerCopy` — and that message says nothing a reader of a rate can act on.
            if isinstance(e, RateCompositionError):
                raise _fail(str(e), text) from None
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
            "use the Python API for a continuous mapping, or a Table for a discrete one", text)
    if name not in _NAMES:
        raise _unknown_name(name, text)
    return _NAMES[name]


def parse_rate(text: object):
    """Read a rate in its written form and return the spec the ``simulate_*`` functions take.

    ``text`` is the expression (``"1.0"``, ``"Global(1.0)"``,
    ``"1.0 * OnTime({0: 1.0, 3: 0.3})"``); a number passes through, so a ``--params`` value that is
    already a TOML float needs no special case. The result is a number, a scope wrapper, a modifier,
    or a ``Rate`` — all four are what ``as_rate`` accepts, so every level takes it as it is.

    Raises `RateSyntaxError` (a ``ValueError``) for anything outside the grammar, and lets the
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
    # A backslash in a path is not an escape. `Driven('C:\\Users\\me\\t.tsv', …)` is well formed to
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
    if isinstance(value, int):        # "1" is the rate 1.0, not the integer 1
        return float(value)
    return value


def written_form(spec: object) -> str:
    """Render a rate spec back as the expression that produced it — the inverse of `parse_rate()`.

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
        # A `SetBy` replaces the base, so writing one in front of it produces the very expression
        # `SetBy` refuses — and that expression is what a run's log records, so the record of the
        # model was one that could not be run again.
        if isinstance(spec, _modifiers.SetBy):
            return repr(spec)
        return f"1.0 * {spec!r}"
    if isinstance(spec, Rate):
        # SPEC §5: a replaced base is written first, because anything to its left is a base it would
        # discard. `sorted` is stable, so the factors keep the order they were written in.
        mods = sorted(spec.modifiers, key=lambda m: not isinstance(m, _modifiers.SetBy))
        if mods and isinstance(mods[0], _modifiers.SetBy):
            # no head: the driver supplies the number. The scope goes with it and loses nothing —
            # a scope cannot wrap a `SetBy` (there is no base for it to wrap), so the scope here is
            # always the level's default, which reading the text back at that level restores.
            return " * ".join(repr(m) for m in mods)
        head = (f"{type(spec.scope).__name__}({float(spec.base)!r})" if spec.scope is not None
                else repr(float(spec.base)))
        return " * ".join([head, *(repr(m) for m in mods)])
    return repr(spec)


__all__ = ["parse_rate", "written_form", "RateSyntaxError"]
