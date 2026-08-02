# Contributing

Thanks for looking. Bug reports and pull requests are welcome.

## Reporting a bug

Open an [issue](https://github.com/AADavin/zombi2/issues) with the command or the few lines of
Python that reproduce it, the **seed**, and `zombi2 --version`. A run with a seed is reproducible,
so that is usually enough to see the same thing you saw.

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
  `### Added` / `### Changed` / `### Fixed`;
- **update [`manual/book/appendix-b.md`](manual/book/appendix-b.md)** if it adds or changes a file
  ZOMBI2 writes — that table is the catalogue of every output, and the docs site single-sources it.

New models belong to a level and are added as described in
[`docs/design/SPEC.md`](docs/design/SPEC.md) §8. If you are proposing something larger, open an
issue first — it is cheaper to agree on the shape before the code exists.

## What ZOMBI2 promises

The same seed and the same version give the same output, on any machine and any supported Python.
Anything that changes what a seed produces is a breaking change; see
[Reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/).
