# Publishing a release

ZOMBI2 ships as **two** PyPI packages, released together in lockstep:

| Package       | What it is                              | Wheels |
| ------------- | --------------------------------------- | ------ |
| `zombi2`      | the pure-Python library and CLI         | one universal wheel + sdist |
| `zombi2_core` | the compiled Rust engine (`rust/`)      | one **abi3** wheel per platform (CPython 3.10+) + sdist |

`zombi2` pins `zombi2_core` to the exact same version, so both must be published
at the same version for `pip install zombi2` to resolve.

Publishing is driven by [`.github/workflows/release.yml`](https://github.com/AADavin/zombi2/blob/main/.github/workflows/release.yml)
and uses **PyPI Trusted Publishing** (OIDC) — there are no API tokens to store or
rotate. The build jobs run on every trigger; the *publish* jobs are gated:

- **PyPI** — only when a **GitHub Release is published**.
- **TestPyPI** — only on a **manual run** (`workflow_dispatch`) with the *testpypi* box ticked.

Until the one-time setup below is done, the publish jobs fail safely at the OIDC
step and nothing is uploaded.

## One-time setup (maintainer)

Do this once per index (PyPI and, if you want a dry run, TestPyPI).

1. **Create the GitHub environments.** In the repo, *Settings → Environments*,
   add all four: `pypi`, `pypi-core`, `testpypi`, and `testpypi-core` (one per
   package per index — see step 2 for why). Optionally add required reviewers to
   `pypi` / `pypi-core` for an extra manual gate before anything is uploaded.

2. **Register a trusted publisher for each package.** On
   [PyPI](https://pypi.org/manage/account/publishing/) and
   [TestPyPI](https://test.pypi.org/manage/account/publishing/), add a *pending*
   GitHub publisher for **both** `zombi2` and `zombi2_core`. They share the same
   owner `AADavin`, repository `zombi2`, and workflow filename `release.yml`, but
   each package uses its **own environment**. This is required: PyPI refuses two
   *pending* publishers that match on owner/repo/workflow/environment, so the
   environment is what tells them apart.

   | PyPI project  | Index         | Environment     |
   | ------------- | ------------- | --------------- |
   | `zombi2`      | pypi.org      | `pypi`          |
   | `zombi2_core` | pypi.org      | `pypi-core`     |
   | `zombi2`      | test.pypi.org | `testpypi`      |
   | `zombi2_core` | test.pypi.org | `testpypi-core` |

   That's four registrations total (two packages × two indexes). PyPI turns each
   pending publisher into a real project the first time it receives an upload.

## Cutting a release

1. **Bump the version.** Each package's version is single-sourced, so there are
   just two files to edit, plus the lockstep pin:
   - `zombi2/__init__.py` → `__version__` — the zombi2 version (`pyproject.toml`
     declares `dynamic = ["version"]` and hatchling reads it from here).
   - `rust/Cargo.toml` → `[package].version` — the engine version (`rust/pyproject.toml`
     declares `dynamic = ["version"]` and maturin reads it from here). Plain
     SemVer, e.g. `0.2.0`.
   - `pyproject.toml` → the `zombi2_core==<version>` pin, bumped to match
     `rust/Cargo.toml` so the two stay in lockstep.

   All three must land on the same version number.

2. **Dry run on TestPyPI (recommended).** Push the bump, then run the *Release*
   workflow manually (*Actions → Release → Run workflow*) with **testpypi**
   ticked. Confirm both packages appear on TestPyPI and install cleanly:

   ```bash
   pip install --index-url https://test.pypi.org/simple/ \
       --extra-index-url https://pypi.org/simple/ zombi2
   ```

3. **Publish.** Draft a **GitHub Release** for the tag and publish it. That
   triggers the wheels build across all platforms and uploads both packages to
   PyPI.

4. **Verify.**

   ```bash
   pip install zombi2 && python -c "import zombi2; print(zombi2.__version__, zombi2.rust_available())"
   ```

## Notes

- The engine builds a **single abi3 wheel per platform** (`cp310-abi3`), usable on
  CPython 3.10 and every later version without a rebuild — so new Python releases
  do not need a new wheel matrix.
- Platforms built: Linux `x86_64` + `aarch64` (manylinux), macOS `x86_64` +
  `arm64`, Windows `x64`. Users on other platforms (e.g. musllinux) fall back to
  building `zombi2_core` from the sdist, which needs a Rust toolchain.
