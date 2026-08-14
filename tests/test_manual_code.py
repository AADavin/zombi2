"""Every code block in the manual is executed, and every ``zombi2`` command line in it is parsed.

The manual drifted away from the code and nothing noticed. Chapter 4 documented a `max_family_size`
API that had been removed three releases earlier — the default it named raised, and the spelling it
called "refused, and deliberately so" was the one that worked — and that text is included verbatim
into the published website. The modifier tables were already checked against the engines
(`test_validation.py`); the prose around them was not.

So: **a code block in the manual is a test**. Python blocks run, in order, in one namespace per
chapter, so a block may use what an earlier one defined. ``zombi2`` command lines in bash blocks are
handed to the real parser, which rejects an unknown flag or an invalid choice — the cheap half of
this, and the half that catches a renamed subcommand or a `--format` value that no longer exists.

**A block that only assigns is not a test.** Chapter 4's stale block would have passed a naive
"execute it" check, because ``scope.PerLineage(10)`` is still a perfectly good object — it just is
not a cap any more. A block earns its keep by *calling* the thing it documents, which is also what a
reader wants to copy. Where a block here is illustrative rather than runnable, it says so with

    <!-- doc-test: skip -->

on the line before the fence, and the reason belongs in the same comment. Reach for that sparingly:
every skip is a paragraph nothing is checking.
"""

from __future__ import annotations

import ast
import contextlib
import io
import os
import pathlib
import re
import shlex

import pytest

MANUAL = pathlib.Path(__file__).resolve().parent.parent / "manual" / "book"

#: The vocabulary a rate expression may use: what `zombi2.params` exports, and what a rate answers to.
_KNOWN_CTORS = {n for n in dir(__import__("zombi2.params", fromlist=["x"])) if n[:1].isupper()}
_KNOWN_VERBS = {a for a in dir(__import__("zombi2.params", fromlist=["x"]).PerLineage(1.0))
                if not a.startswith("_")}

#: blocks whose *language* is not something we can run
#: Fenced languages the tests read. ``toml`` is here because a page that documents ``--params`` shows
#: the file as a toml block, and the command test writes it to disk so the flag resolves against the
#: real file; each test filters to the language it wants.
_FENCED = ("python", "bash", "toml")

_FENCE = re.compile(r"^```(\w+)?\s*$")
_SKIP = re.compile(r"<!--\s*doc-test:\s*skip")



#: A backticked cell of a manual table. The reference tables are where most of the rate grammar is
#: written down — chapter 7's clocks, and the "From the literature" tables of chapters 3, 4, 7 and 8 —
#: and a table is not a python block, so nothing checked a word of them. `ByFamily` sat in one of
#: those tables for two releases after the class was renamed, and `*` stayed in another after it
#: stopped composing a rate.
_CELL = re.compile(r"`([^`]+)`")

#: The names a rate expression starts with. A cell headed by one of these is meant to be code; a cell
#: headed by anything else is prose in backticks — a flag, a file name — and is not this test's
#: business. `simulate_*(...)` cells are signatures rather than values, so they are not here either.
_GRAMMAR = ("PerCopy", "PerLineage", "PerSite", "PerChromosome", "Extent", "Clade", "Clades",
            "Between", "Drift", "LogNormal", "Gamma", "Random", "TotalDiversity", "Curve", "Table")


def _grammar_cells(text: str) -> list[tuple[int, str]]:
    """``(line, expression)`` for every rate-grammar expression written in a table cell."""
    out = []
    for n, line in enumerate(text.splitlines(), start=1):
        if not line.startswith("|"):
            continue
        for cell in _CELL.findall(line):
            expr = cell.strip()
            # `substitution = PerSite(…)` is an assignment and the value is what matters; the `=` in
            # `simulate_continuous(rate=…)` is a keyword argument and splitting on it leaves an
            # unbalanced tail, so only a bare name on the left counts as an assignment
            head, _, tail = expr.partition("=")
            if tail and head.strip().isidentifier():
                expr = tail.strip()
            if expr.startswith(_GRAMMAR) and "(" in expr:
                out.append((n, expr))
    return out


def _names_used(expr: str) -> tuple[set, set]:
    """The constructors and the chained verbs an expression names, without running it.

    The cells are **schemas**, not expressions — the ellipsis in ``changing_at({…})`` stands for a
    dict and the one in ``Gamma(…)`` for two arguments — so substituting a value and evaluating
    invents failures. What can be checked without guessing is the vocabulary: a renamed class or a
    removed verb is a name that no longer resolves, which is exactly what rotted before."""
    tree = ast.parse(expr.replace("…", "None").replace("...", "None"), mode="eval")
    ctors = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and n.id[:1].isupper()}
    verbs = {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    return ctors, verbs


def _blocks(path: pathlib.Path) -> list[tuple[int, str, str]]:
    """``(line number, language, source)`` for each fenced block, skips excluded."""
    lines = path.read_text(encoding="utf-8").splitlines()
    out, i = [], 0
    while i < len(lines):
        m = _FENCE.match(lines[i])
        if not m or m.group(1) not in _FENCED:
            i += 1
            continue
        lang, start = m.group(1), i
        body, i = [], i + 1
        while i < len(lines) and not _FENCE.match(lines[i]):
            body.append(lines[i])
            i += 1
        i += 1
        skipped = any(_SKIP.search(ln) for ln in lines[max(0, start - 2):start])
        if not skipped:
            out.append((start + 1, lang, "\n".join(body)))
    return out


#: docs/ pages are mostly one-line `--8<--` includes of a manual chapter, and testing those twice
#: would be pointless. A page with real prose in it is its own source, though, and gets checked here
#: like a chapter — otherwise "move it out of the book" quietly means "stop testing it".
DOCS = pathlib.Path(__file__).resolve().parent.parent / "docs"


def _standalone_docs() -> list[pathlib.Path]:
    """Pages the **site publishes** that are their own source, rather than an include of a chapter.

    ``docs/design/`` is out of scope, as it is for mkdocs: SPEC is a design document, and it names dead
    flags on purpose — it keeps a table of retired vocabulary, so a check for "every flag named here
    exists" is precisely wrong about it."""
    return sorted(p for p in DOCS.rglob("*.md")
                  if "design" not in p.parts
                  and "--8<--" not in p.read_text(encoding="utf-8"))


#: the root README: the first code anyone runs, and the one page here that is neither a chapter nor a
#: docs page. Nothing executed it, so its quickstart could break while every chapter stayed green.
README = pathlib.Path(__file__).resolve().parent.parent / "README.md"


def _chapters() -> list[pathlib.Path]:
    """Everything with prose to check: the front page, the book, plus any docs page that is its own
    source."""
    return ([README] + sorted(p for p in MANUAL.glob("*.md") if p.name != "README.md")
            + _standalone_docs())


def _zombi2_lines(source: str) -> list[list[str]]:
    """The ``zombi2 …`` invocations in a bash block, as token lists, with line continuations joined.

    Other lines — ``pip install``, ``cd``, a pipe into ``awk`` — are not ours to validate. A pipeline
    that *starts* with zombi2 has its own command taken and the rest dropped, so the flags still get
    checked."""
    joined = re.sub(r"\\\s*\n\s*", " ", source)
    out = []
    for raw in joined.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line.startswith("zombi2 "):
            continue
        head = re.split(r"[|>]", line)[0]
        try:
            tokens = shlex.split(head)[1:]
        except ValueError:                    # an unbalanced quote is a doc bug of its own
            pytest.fail(f"could not tokenise: {head!r}")
        if tokens:
            out.append(tokens)
    return out


@pytest.mark.parametrize("chapter", _chapters(), ids=lambda p: p.name)
def test_the_python_in_the_manual_runs(chapter, tmp_path, monkeypatch):
    # one namespace per chapter, blocks in order: the manual builds up an example across a section,
    # and a block that needs `sp` from three paragraphs earlier is a normal thing to write.
    blocks = [(n, src) for n, lang, src in _blocks(chapter) if lang == "python"]
    if not blocks:
        pytest.skip(f"no python blocks in {chapter.name}")
    monkeypatch.chdir(tmp_path)               # blocks that write("out/…") land in a temp directory
    # `mytree.nwk` is the manual's stock name for "a tree you already have", so provide one: that
    # turns every block that reads it from an untestable illustration into a live check.
    (tmp_path / "mytree.nwk").write_text(
        "((n2:1.0,n3:1.0)n1:0.5,(n5:1.0,n6:1.0)n4:0.5)n0:0.2;\n", encoding="utf-8")
    ns: dict = {"__name__": "__main__"}
    for line, src in blocks:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(src, f"{chapter.name}:{line}", "exec"), ns)   # noqa: S102
        except Exception as e:
            pytest.fail(f"{chapter.name}:{line} — the manual's own code raised "
                        f"{type(e).__name__}: {e}\n\n{src}")


@pytest.mark.parametrize("chapter", _chapters(), ids=lambda p: p.name)
def test_the_zombi2_commands_in_the_manual_parse(chapter, tmp_path, monkeypatch):
    # the real parser, with the handlers stubbed out: this validates every flag name, every choice
    # and every type without running a simulation, so a renamed subcommand or a --format value that
    # no longer exists fails here rather than in someone's terminal.
    from zombi2.cli import main as cli_main

    cmds = [(n, c) for n, lang, src in _blocks(chapter) if lang == "bash"
            for c in _zombi2_lines(src)]
    if not cmds:
        pytest.skip(f"no zombi2 command lines in {chapter.name}")
    monkeypatch.setattr(cli_main, "_RUN", {k: (lambda args, parser: 0) for k in cli_main._RUN})
    monkeypatch.chdir(tmp_path)
    # A page that shows a --params file writes it as a ```toml block whose first line names it. Put
    # those on disk so the commands below read the very file the page shows: --params is resolved by
    # the parser, so this checks the file's keys against the real options rather than eyeballing them.
    for _, lang, src in _blocks(chapter):
        if lang != "toml":
            continue
        first = src.splitlines()[0].strip() if src.strip() else ""
        if first.startswith("#") and first.endswith(".toml"):
            (tmp_path / first.lstrip("#").strip()).write_text(src, encoding="utf-8")
    for line, tokens in cmds:
        try:
            rc = cli_main.main(tokens)
        except SystemExit as e:
            pytest.fail(f"{chapter.name}:{line} — the manual's own command line is rejected by the "
                        f"parser (exit {e.code}): zombi2 {' '.join(tokens)}")
        assert rc == 0, f"{chapter.name}:{line} — zombi2 {' '.join(tokens)} returned {rc}"


#: the public entry points, and a setup small enough to probe one keyword against in milliseconds
def _entry_points():
    from zombi2 import genomes, joint, sequences, species, traits
    return {"joint": joint.simulate_joint,
            "species": species.simulate_species_tree,
            "genomes": genomes.simulate_genomes_family,
            "ordered": genomes.simulate_genomes_ordered,
            "nucleotide": genomes.simulate_genomes_nucleotide,
            "sequences": sequences.simulate_sequences,
            "continuous": traits.simulate_continuous,
            "discrete": traits.simulate_discrete}


def _probe(name: str, value) -> str | None:
    """Try ``value`` as the keyword ``name`` on whichever entry point takes it. Returns the rejection
    message, or ``None`` if some entry point accepted it.

    This is the check chapter 4 needed. Executing ``max_family_size = scope.PerLineage(10)`` succeeded
    for three releases after the API stopped accepting it, because binding a name cannot fail — only
    *passing* it can. So the harness passes it."""
    import inspect

    from zombi2 import species
    tree = species.simulate_species_tree(birth=1.0, n_extant=4, seed=1)
    from zombi2.traits import DiscreteTrait
    base = {"species": dict(birth=1.0, n_extant=4, seed=1),
            # a joint run takes a rate driven by a LIVE level, which no other engine accepts
            "joint": dict(birth=1.0, n_extant=4, seed=1,
                          trait=DiscreteTrait(states=("small", "large"), switch=0.3)),
            "genomes": dict(tree=tree, duplication=0.1, initial_families=2, seed=1),
            "ordered": dict(tree=tree, duplication=0.1, initial_families=2, seed=1),
            "nucleotide": dict(tree=tree, root_length=200, seed=1),
            "continuous": dict(tree=tree, rate=1.0, seed=1),
            "discrete": dict(tree=tree, states=("a", "b"), switch=0.1, seed=1)}
    last = None
    for key, fn in _entry_points().items():
        if key not in base or name not in inspect.signature(fn).parameters:
            continue
        try:
            fn(**{**base[key], name: value})
            return None                        # some entry point takes it: the manual is current
        except FileNotFoundError:
            return None      # a conditioned rate went looking for its driver: the value was taken
        except (TypeError, ValueError) as e:
            last = f"{fn.__name__}({name}=…) -> {type(e).__name__}: {e}"
    return last


def test_a_value_the_manual_assigns_to_a_parameter_is_still_accepted(tmp_path, monkeypatch):
    """Every ``&lt;parameter&gt; = &lt;value&gt;`` in the manual is passed to the API that names it.

    Chapter 4 is why. It showed three ``max_family_size`` values, every one of them a valid object and
    none of them a value the API would take any more; executing the block succeeded, because an
    assignment cannot fail. Binding is not using. So the harness binds *and* uses."""
    import inspect

    params = set()
    for fn in _entry_points().values():
        with contextlib.suppress(TypeError, ValueError):
            params |= set(inspect.signature(fn).parameters)
    params -= {"seed", "progress", "tree", "genomes", "model", "record", "parallel", "flat",
               "outputs", "stream_to"}                    # plumbing, not model description

    monkeypatch.chdir(tmp_path)
    failures, checked = [], 0
    for chapter in _chapters():
        ns: dict = {"__name__": "__main__"}
        for line, lang, src in _blocks(chapter):
            if lang != "python":
                continue
            with contextlib.suppress(Exception), contextlib.redirect_stdout(io.StringIO()):
                exec(compile(src, f"{chapter.name}:{line}", "exec"), ns)   # noqa: S102
            for node in ast.walk(ast.parse(src)):
                if not isinstance(node, ast.Assign):
                    continue
                for target in node.targets:
                    if not isinstance(target, ast.Name) or target.id not in params:
                        continue
                    if target.id not in ns:              # the block did not run: covered elsewhere
                        continue
                    checked += 1
                    if (why := _probe(target.id, ns[target.id])) is not None:
                        failures.append(f"{chapter.name}:{line} — {why}")
        # and the same spellings written in the chapter's *tables*, which no block executes
        # and the rate grammar written in the chapter's *tables*, which no block executes
        for line, expr in _grammar_cells(chapter.read_text(encoding="utf-8")):
            try:
                ctors, verbs = _names_used(expr)
            except SyntaxError:
                continue        # a schema an ellipsis cannot be filled into; nothing to resolve
            checked += 1
            for missing in sorted(ctors - _KNOWN_CTORS):
                failures.append(f"{chapter.name}:{line} (table) — `{expr}` names {missing}, "
                                f"which zombi2.params no longer has")
            for missing in sorted(verbs - _KNOWN_VERBS):
                failures.append(f"{chapter.name}:{line} (table) — `{expr}` chains .{missing}(), "
                                f"which a rate no longer has")
    # A canary on the walk, not a coverage target: it fires when the walk stops finding the manual's
    # values at all. It counts the tables as well as the blocks — chapter 7's clock section traded
    # its snippets for a reference table, and a table is where most of the rate grammar is written
    # down now. Raise this when a chapter adds spellings; do not lower it to accommodate a walk that
    # has quietly broken.
    assert checked >= 15, f"only probed {checked} values; the walk is not finding the manual's values"
    assert not failures, ("the manual assigns values the API no longer accepts:\n  "
                          + "\n  ".join(failures))


def test_the_harness_sees_the_manual():
    """A guard on the guard: if the glob or the fence regex breaks, every test above turns green by
    finding nothing. This is the tripwire for that."""
    assert len(_chapters()) >= 10
    assert README in _chapters(), "the front page has dropped out of the sweep"
    assert _standalone_docs(), "docs/ pages are all includes — or the include check has broken"
    total = sum(len(_blocks(p)) for p in _chapters())
    assert total >= 60, f"expected the manual's ~80 code blocks, found {total}"
    assert os.environ.get("ZOMBI2_SKIP_MANUAL") is None


def test_every_flag_the_manual_names_in_prose_exists():
    """Flags mentioned in **prose** are checked too, not only ones inside a bash block.

    The migration appendix is mostly tables mapping a ZOMBI v1 key to a v2 flag, and a table is prose:
    the command-line check above never sees it. A mapping document whose right-hand column has rotted
    is worse than no mapping document, because a reader trusts it and then goes hunting for a flag that
    was renamed two releases ago — which is exactly how chapter 4 came to document an API that had been
    removed."""
    import contextlib as _c
    import io as _io

    from zombi2.cli import main as cli_main

    declared: set[str] = set()
    for argv in (["species"], ["genomes"], ["sequences"], ["traits"], ["joint"],
                 ["tools"], ["tools", "tree"], ["tools", "treedist"], ["tools", "format"]):
        buf = _io.StringIO()
        with _c.suppress(SystemExit), _c.redirect_stdout(buf), _c.redirect_stderr(buf):
            cli_main.main([*argv, "--help"])
        declared |= set(re.findall(r"(--[a-z][a-z0-9-]+)", buf.getvalue()))
    assert "--n-extant" in declared, "the help scrape found nothing; this test would pass vacuously"

    missing = []
    for chapter in _chapters():
        text = chapter.read_text(encoding="utf-8")
        for flag in sorted(set(re.findall(r"`(--[a-z][a-z0-9-]+)", text))):
            if flag not in declared:
                missing.append(f"{chapter.name}: {flag}")
    assert not missing, ("the manual names flags no command declares — renamed, or never existed:\n  "
                         + "\n  ".join(missing))

def test_every_figure_is_well_formed_xml():
    """An SVG that is not well-formed XML does not render — the browser refuses the whole file.

    `conditioning.svg` and `joining.svg` both shipped broken for exactly this reason: a comment
    explaining the colour convention mentioned ``var(-``\u2011``ink)``, and a double hyphen is illegal
    inside an XML comment. Nothing complained; the images were simply blank in the README, which is
    the sort of thing only a reader notices."""
    import xml.etree.ElementTree as ET

    figures = sorted((MANUAL / "figures").glob("*.svg")) + sorted(
        (pathlib.Path(__file__).resolve().parent.parent / "assets").glob("*.svg"))
    assert figures, "no figures found — has the layout moved?"
    broken = []
    for svg in figures:
        try:
            ET.parse(svg)
        except ET.ParseError as e:
            broken.append(f"{svg.name}: {e}")
    assert not broken, "these figures will not render:\n  " + "\n  ".join(broken)


def test_every_figure_the_book_names_can_be_built_and_served():
    """A chapter names a figure as the artefact the *book* builds — `figures/NAME_print.png` for a
    hand-authored chapter diagram, `figures/NAME.pdf` for a generated one. The website reaches the
    same figure by a different route: a docs hook rewrites the name to `../img/NAME.svg`, which has
    to be a symlink under `docs/`, because MkDocs serves only what is there.

    Four figures had no symlink, so a reader of the Conditioning chapter on the site saw one of its
    four figures and a reader of Joining one of two, while the PDF had them all — for weeks.
    `mkdocs --strict` does not catch it: `validation.links` reports a missing target at `info`.
    """
    named, missing = set(), []
    for chapter in sorted(MANUAL.glob("*.md")):
        for m in re.finditer(r"\]\(figures/([A-Za-z0-9_-]+?)(_print\.png|\.pdf)\)", chapter.read_text(encoding="utf-8")):
            named.add((m.group(1), m.group(2), chapter.name))

    root = MANUAL.parent.parent
    for name, kind, chapter in sorted(named):
        source = (MANUAL / "figures" / f"{name}.svg" if kind == "_print.png"
                  else root / "figures" / "svg" / f"{name}.svg")
        if not source.exists():
            missing.append(f"{chapter} names {name}{kind}, but {source.relative_to(root)} is not there")
        link = root / "docs" / "img" / f"{name}.svg"
        if not link.exists():
            missing.append(f"{chapter} names {name}{kind}, but docs/img/{name}.svg is not there, "
                           f"so the figure is missing on the website")

    assert named, "no figures found; this test would pass vacuously"
    assert not missing, "\n  ".join(["", *missing])
