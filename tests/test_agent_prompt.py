"""The agent-setup prompt, ``docs/agent-setup/prompt.txt``, is a test: an agent fetches it and follows
it, so a flag it names that the parser no longer takes, a Python call that no longer exists, a rate
written in a form the grammar no longer parses or a URL that 404s would send that agent wrong in
someone else's project. Every command line in it goes through the real parser, every indented
Python block runs, every rate expression parses, every output file it names is in Appendix B, and
every URL it gives is one the site serves (checked by path against the docs tree, not the network).
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import re
import shlex

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
PROMPT = ROOT / "docs" / "agent-setup" / "prompt.txt"
APPENDIX_B = ROOT / "manual" / "book" / "appendix-b.md"
SITE = "https://aadavin.github.io/zombi2/"


def _indented_blocks(text: str) -> list[tuple[int, str]]:
    """The four-space-indented blocks of the prompt, each with its first line number."""
    blocks, current, start = [], [], None
    for n, line in enumerate(text.splitlines(), 1):
        if line.startswith("    "):
            if start is None:
                start = n
            current.append(line[4:])
        elif current:
            blocks.append((start, "\n".join(current)))
            current, start = [], None
    if current:
        blocks.append((start, "\n".join(current)))
    return blocks


def _commands() -> list[tuple[int, list[str]]]:
    out = []
    for start, block in _indented_blocks(PROMPT.read_text(encoding="utf-8")):
        for i, raw in enumerate(block.splitlines()):
            line = raw.split("#", 1)[0].strip()
            if line.startswith("zombi2 "):
                out.append((start + i, shlex.split(line)[1:]))
    return out


def test_the_prompt_exists_and_names_the_version_that_is_installed():
    import zombi2
    text = PROMPT.read_text(encoding="utf-8")
    assert "pip install zombi2" in text
    stated = re.search(r"Version (\d+\.\d+\.\d+) or later", text)
    assert stated, "the prompt says which version it describes"
    assert tuple(map(int, stated.group(1).split("."))) <= tuple(map(int, zombi2.__version__.split(".")[:3]))


def test_every_command_line_in_the_prompt_parses(tmp_path, monkeypatch):
    from zombi2.cli import main as cli_main
    cmds = _commands()
    assert len(cmds) >= 12
    monkeypatch.setattr(cli_main, "_RUN", {k: (lambda args, parser: 0) for k in cli_main._RUN})
    monkeypatch.chdir(tmp_path)
    for line, tokens in cmds:
        if tokens[:1] == ["--version"]:
            continue
        try:
            rc = cli_main.main(tokens)
        except SystemExit as e:
            pytest.fail(f"prompt.txt:{line} — rejected by the parser (exit {e.code}): "
                        f"zombi2 {' '.join(tokens)}")
        assert rc == 0, f"prompt.txt:{line} — zombi2 {' '.join(tokens)} returned {rc}"


def test_every_rate_expression_in_the_prompt_parses():
    from zombi2.params.parse import parse_rate
    text = PROMPT.read_text(encoding="utf-8")
    exprs = re.findall(r'"((?:Per\w+|Global)\([^"]*\))"', text)
    assert len(exprs) >= 4
    for expr in exprs:
        parse_rate(expr)


def test_the_python_in_the_prompt_runs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    blocks = [b for _, b in _indented_blocks(PROMPT.read_text(encoding="utf-8"))
              if b.lstrip().startswith(("from zombi2", "import zombi2"))]
    assert blocks
    ns: dict = {"__name__": "__main__"}
    for src in blocks:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(src, "prompt.txt", "exec"), ns)      # noqa: S102
    assert ns["observed"] and ns["trees"]


def test_every_output_file_the_prompt_names_is_in_appendix_b():
    appendix = APPENDIX_B.read_text(encoding="utf-8")
    text = PROMPT.read_text(encoding="utf-8")
    section = text.split("5. What a run writes")[1].split("From Python")[0]
    names = set(re.findall(r"[A-Za-z_<>]+\.(?:tsv|nwk|json|fasta|zombi2)", section))
    names |= {"gene_tree_fam<f>_complete.nwk"}
    missing = [n for n in names if n.replace("<f>", "<f>") not in appendix]
    assert not missing, missing


def test_every_url_in_the_prompt_is_a_page_of_this_site():
    text = PROMPT.read_text(encoding="utf-8")
    urls = sorted(set(re.findall(r"https://aadavin\.github\.io/zombi2/\S+", text)))
    assert urls
    for url in urls:
        path = url[len(SITE):].rstrip(".,)").rstrip("/")     # a sentence's end, then the slash
        if path == "gallery.html":
            assert (ROOT / "web" / "gallery.html").exists()
        elif path.startswith("docs/"):
            rel = path[len("docs/"):]
            candidates = [ROOT / "docs" / rel, ROOT / "docs" / (rel + ".md"),
                          ROOT / "docs" / rel / "index.md"]
            assert any(c.exists() for c in candidates), url
        else:
            pytest.fail(f"unexpected URL shape: {url}")
