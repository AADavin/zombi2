# Reproducibility

**A seed names a run.** The same ZOMBI2 version, the same arguments and the same `seed` give the same
dataset — the same tree, events, sequences and trait values — on any machine, any supported Python,
Linux, macOS or Windows.

```bash
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 42
```

That command and the version it ran under are the whole record.

## What is guaranteed

- **Same version + same arguments + same seed ⇒ identical output**, across operating systems, CPU
  architectures and Python 3.10 – 3.13. This is checked rather than assumed:
  `tests/test_reproducibility.py` hashes a run of every level against a recorded digest, on six CI
  configurations.
- **Every run records how to repeat itself.** Each command writes a `run.zombi2` report with the
  arguments, the version and a `TO REPRODUCE` block. A conditioned run lists its driver's command
  first, so the block runs top to bottom.
- **A seed you did not give is still written down.** Leave `--seed` off and one is drawn from the
  operating system and printed into the report; from Python, `seed=None` puts it on `result.seed`.
  There is no such thing as an unrepeatable ZOMBI2 run.
- **Each level is seeded separately.** `zombi2 genomes … --seed 7` on a tree grown with `--seed 42`
  gives the same genomes whatever else happened in between.
- **One seed on two levels is not the same numbers twice.** Each level draws from its own stream,
  spawned from your seed under a per-level key, so `--seed 42` on both commands gives two
  *independent* runs. Before 0.28.0 both opened the same generator from the same integer, so a shared
  seed quietly coupled them.
- **A run written to disk is the same run in memory.** Branch lengths are written at full precision,
  so two commands and two Python calls give the same history.

## What is not

- **Not across versions.** A change to an engine moves the random stream, and every seed then means
  something different. *"ZOMBI2 0.23.0, seed 42"* is a reproducible reference; *"seed 42"* is not.
  When a release changes what a seed produces it is a **breaking** change, called out in
  [`CHANGELOG.md`](https://github.com/AADavin/zombi2/blob/main/CHANGELOG.md) with a minor version
  bump. Pin the version for anything you intend to rerun in a year.

- **Not between a serial run and a `--parallel` one.** `--parallel` is a separate engine: it gives
  each unit of work its own stream spawned from your seed. So the same seed gives the same output at
  **any** worker count, but a different — equally valid — realisation from the serial default. Choose
  once at the start of a study and record which you used.

    From Python a parallel run needs its call under a `__main__` guard, the standard requirement for
    anything that starts worker processes. Leave it out and the run stops with a message naming the
    guard. The `zombi2` command already has one, and a notebook has no script to re-import.

    ```python
    if __name__ == "__main__":
        g = genomes.simulate_genomes_family(sp, duplication=0.2, seed=1, parallel=8)
    ```

- **Not, strictly, across numpy major versions.** numpy does not promise a `Generator` method keeps
  its stream forever. In practice the digests are unchanged from numpy 1.26 to 2.5, and ZOMBI2 pins
  `numpy<3`. If you are reproducing an old run and the numbers are close rather than equal, check
  this after the ZOMBI2 version; `run.zombi2` records the numpy version for that reason.

- **Not bit-identical floating point across platforms.** `exp`, `log` and the matrix routines behind
  a substitution model differ in their last bit or two between platforms, so branch lengths and trait
  values can differ around the sixteenth significant digit. Nothing discrete depends on that, which
  is why the digests are taken at nine significant digits.

## Reproducing someone else's run

```bash
cat out/run.zombi2
```

The report holds the version, every argument as it was resolved, and the `TO REPRODUCE` block.
Install that version and run the block. If the version in the report is not the one you have, say so
when you report a difference: it explains most of them.
