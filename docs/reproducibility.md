# Reproducibility

A simulator's output is only useful if someone else can get it back. ZOMBI2's contract is short:

**A seed names a run.** The same ZOMBI2 version, given the same arguments and the same `seed`,
produces the same dataset — the same tree, the same events, the same sequences, the same trait
values — on any machine, on any supported Python, on Linux, macOS or Windows.

```bash
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 42
```

That command is the whole record. Anyone who has it, and the version it ran under, has the dataset.

## What is guaranteed

- **Same version + same arguments + same seed ⇒ identical output.** Across operating systems, CPU
  architectures and Python 3.10 – 3.13. This is checked, not assumed: `tests/test_reproducibility.py`
  hashes a run of every level and compares it against a recorded digest, and CI runs that test on six
  configurations. Six machines, six independently randomised Python hash seeds, one number.
- **Every run records how to repeat itself.** Each command writes a `run.zombi2` report holding the
  arguments, the version, and a `TO REPRODUCE` block of commands that regenerate it. A run that
  conditions on another level lists the driver's command first, so the block runs top to bottom.
- **A seed you did not give is still written down.** Leave `--seed` off and one is drawn from the
  operating system — and printed into the report, so the run is reproducible afterwards even though
  it was not planned to be. From Python, `seed=None` does the same thing: the drawn seed is on the
  result as `result.seed`, so an interesting realisation found while exploring is never lost. There
  is no such thing as an unrepeatable ZOMBI2 run.
- **Each level is seeded separately.** `zombi2 genomes … --seed 7` on a tree grown with `--seed 42`
  gives the same genomes whatever else has happened in between, because the genome stream depends on
  the tree and its own seed, not on the order you ran things in.
- **One seed on two levels does not mean the same numbers twice.** Each level draws from its own
  stream, spawned from your seed under a per-level key, so `--seed 42` on both the species and the
  genomes command gives you two *independent* runs — which is what SPEC §2 means by two levels being
  independent. Before 0.28.0 both levels opened the same generator from the same integer, so a shared
  seed quietly coupled them; if you have a study built that way, its levels are correlated.
- **A run written to disk is the same run in memory.** Branch lengths are written at full precision,
  so the tree the CLI hands to the next level is exactly the tree the previous one produced. The same
  tree and the same seed give the same history whether you run two commands or two Python calls.

## What is not

- **Not across versions.** A change to an engine — a bug fix, a new event, a different way of drawing
  a waiting time — moves the random stream, and every seed then means something different. This is
  why a version number is part of the record and not a footnote: *"ZOMBI2 0.23.0, seed 42"* is a
  reproducible reference and *"seed 42"* is not.

    When a release does change what a seed produces, it is a **breaking** change: it is called out in
    [`CHANGELOG.md`](https://github.com/AADavin/zombi2/blob/main/CHANGELOG.md) and carries a minor
    version bump. Pin the version (`pip install zombi2==0.23.0`) for anything you intend to be able to
    rerun in a year.

- **Not between a serial run and a `--parallel` one.** `--parallel` is a *separate engine*: it gives
  each unit of work — a gene family, a gene tree — its own stream spawned from your seed, rather than
  drawing them all from one. So the same seed gives the same output at **any** worker count, one core
  or thirty-two, but a different (equally valid) realisation from the serial default. Choose once at
  the start of a study rather than switching partway through, and record which you used.

    From Python, a parallel run needs its call under a `__main__` guard, which is the standard
    requirement for anything that starts worker processes — they re-import your script, and without
    the guard they run it again from the top:

    ```python
    if __name__ == "__main__":
        g = genomes.simulate_genomes_family(sp, duplication=0.2, seed=1, parallel=8)
    ```

    Leave it out and the run stops with a message naming the guard. The `zombi2` command already has
    one, and a notebook or `python -c` has no script to re-import, so both are unaffected — the
    library notices and runs single-process there instead.

- **Not, strictly, across numpy major versions.** Every random number ZOMBI2 draws comes from numpy's
  `Generator`, and numpy does not promise that a `Generator` method keeps producing the same stream
  forever — a better algorithm for a distribution is allowed to land in a major release. In practice
  the digests above are unchanged from numpy 1.26 to 2.5, and ZOMBI2 pins `numpy<3` so a major bump
  cannot arrive without a release of ours. But if you are reproducing a years-old run and the numbers
  are close rather than equal, this is the second thing to check after the ZOMBI2 version. The
  `run.zombi2` report records the numpy version for exactly that reason.

- **Not bit-identical floating point where the platform's own maths differs.** `exp`, `log` and the
  matrix routines behind a substitution model are allowed a rounding error in their last bit or two,
  and different platforms use different implementations. Branch lengths and trait values can differ
  around the sixteenth significant digit between a Mac and a Linux box. Nothing discrete — no
  topology, no event, no sequence — depends on that, which is why the digests above are taken at nine
  significant digits: far more precision than any simulated quantity is claimed to, and far away from
  the noise.

## Reproducing someone else's run

Everything you need is in the run directory:

```bash
cat out/run.zombi2
```

The report holds the ZOMBI2 version, every argument as it was resolved (including defaults and the
drawn seed), and the `TO REPRODUCE` block. Install that version and run the block.

If the version in the report is not the one you have, say so when you report a difference — that is
the first thing to check, and it explains most of them.
