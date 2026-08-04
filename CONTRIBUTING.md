# Contributing

Thanks for looking. Bug reports, corrections to the documentation and pull requests are all
welcome, and so are suggestions for models ZOMBI2 does not have yet.

Everything merged into ZOMBI2 is approved by Adrián A. Davín.

## Reporting a bug

Bug reports are always welcome. Open an [issue](https://github.com/AADavin/zombi2/issues) with the
command or the few lines of Python that reproduce it, the **seed**, and `zombi2 --version`. A run
with a seed is reproducible, so that is usually enough to see the same thing you saw.

If you ran it from the command line, the quickest thing to paste is the `TO REPRODUCE` block out of
`run.zombi2` in your run directory: it already holds the version, every resolved argument and the
commands that regenerate the run.

## Telling us the documentation is unclear

**Please do this.** If something feels under-documented, or a page uses a word without ever saying
what it means, that is worth an issue on its own.

ZOMBI2 has been developed with heavy use of large language models. Everything is checked, but plenty
will have been missed: LLMs tend to overcomplicate an explanation and to be loose with terminology,
introducing a word in one chapter and quietly using a different one in the next. If you hit a
sentence you had to read three times, or a term you cannot find defined anywhere, flag it. A
one-line issue saying "I don't know what X means here" is a real contribution and an easy one to
act on.

## Suggesting a model

Some things are missing from ZOMBI2 only because nobody has added them yet. As of 2026-08-04 there
is no model of gene conversion, for instance, and no reason there could not be. If you want one,
open an issue using this template:

```
**The model.** One sentence: what evolves, and how.

**Where it belongs.** Which level (Species, Genomes, Sequences, Traits), and at which genome
resolution if it is a genome model (family, ordered, nucleotide).

**Parameters.** The rates and extents it adds, written the way ZOMBI2 writes rates
(SPEC §5) — e.g. `conversion = 0.1 * ByFamily(spread=0.5)`.

**Output.** Anything new it would have to write, or "nothing new".

**Reference.** One paper that defines or uses the model.

**Independence.** Does it need the units of its level to affect each other — families that
change each other's rates, sites that evolve in each other's context? If you are not sure,
say so.
```

## Setting up

```bash
git clone https://github.com/AADavin/zombi2 && cd zombi2
pip install -e ".[dev]"
```

## Before opening a pull request

Run what CI runs:

```bash
pytest -q && ruff check . && mypy zombi2
```

and, if you touched the documentation, `pip install -e ".[docs]" && mkdocs build --strict`.

A pull request should also:

- **add a test** for a fix or a new behaviour — one that fails without the change;
- **add a line to `CHANGELOG.md`** under `## [Unreleased]` if a user could observe the change
  (new behaviour, a flag, an output file, a public function, a bug fix), grouped under
  `### Added` / `### Changed` / `### Fixed`. CI checks this for the two cases it can see for
  itself — a `--flag` added or removed under `zombi2/cli/`, and a name added to or removed from an
  `__all__` — and you can run the same check with
  `python scripts/check_public_surface.py --base origin/main`. If the change really is invisible to
  users, label the pull request `no changelog`;
- **update [`manual/book/appendix-b.md`](manual/book/appendix-b.md)** if it adds or changes a file
  ZOMBI2 writes — that table is the catalogue of every output, and the docs site single-sources it.

New models belong to a level and are added as described in
[`docs/design/SPEC.md`](docs/design/SPEC.md) §9.

## Code written with an AI agent

Contributions written with an agent are **very** welcome. Not using AI in these times is a red
flag!!!

There is one thing to be careful about, and it is not the code being wrong. Left to grow without
control, agent-written code makes a library bloated: three functions that do the same thing, a new
word for a concept that already had one, an explanation four paragraphs long where one line would
do. Each piece looks reasonable on its own, which is what makes it hard to catch later. So:

- **Read what you are submitting**, and be ready to explain why each part of it is there.
- **Reuse the vocabulary that exists.** [`docs/design/SPEC.md`](docs/design/SPEC.md) fixes the words
  ZOMBI2 uses, and a pull request that invents a synonym will be asked to drop it.
- **Prefer the smaller change.** If a diff is large, say in the description what is essential and
  what is incidental.

We are committed to clarity and to good documentation. That is the bar a contribution is held to,
whoever or whatever wrote it.

## Experimental features

A feature that is elaborate, or whose interface we are not sure of yet, does not have to be settled
to be merged. It can go in marked **experimental**: an `experimental` marker in its docstring. It
lives in its level's package like anything else and is reachable from Python but not from the
command line, because a flag is a promise to keep something.

Promotion out of experimental is a **human decision**, and it is mine. What it takes is the feature
being properly documented: a manual section, its outputs in Appendix B, and ideally an example in
the [gallery](https://aadavin.github.io/zombi2/gallery.html). Graduating then moves no files — the
marker comes off and the flag goes on, and the code was in the right place all along.

## What ZOMBI2 promises

The same seed and the same version give the same output, on any machine and any supported Python.
Anything that changes what a seed produces is a breaking change; see
[Reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/).
