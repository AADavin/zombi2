# ZOMBI2 · performance analysis

A clean, reproducible workspace for benchmarking ZOMBI2 and producing
publication-quality figures. Built to be **iterated on**: benchmarks, rate
regimes, and figures will change as the library evolves, so the moving parts are
kept small and decoupled.

## The one idea

**Measurement and plotting never touch each other.** They meet only through
`results/*.json`.

```
benchmarks.py ──run.py──▶ results/*.json ──plot.py──▶ figures/*.{svg,png}
   (defines)    (measures)   (raw timings)   (renders)    (deliverables)
```

* Re-running one benchmark rewrites only its own JSON file.
* Restyling a figure never re-runs a simulation — `plot.py` reads the JSON.
* Every result file embeds a provenance header (git commit, Rust on/off,
  interpreter, CPU, timestamp) so a curve is never ambiguous.

## Quick start

```bash
cd analyses/performance

python run.py --quick        # fast smoke test (small grids, few repeats)
python run.py                # standard publishable run (up to 3M tips)
python plot.py               # render every figure it has data for
```

Then open `figures/overview.png` (or any standalone figure).

```bash
python run.py --list                 # what benchmarks exist
python run.py species_tree           # run just one
python run.py --full                 # push to 10M tips (heavier; see "Scale")
python run.py write_output           # opt-in extra benchmark
ZOMBI1_DIR=/path/to/ZOMBI python run.py vs_zombi1              # opt-in: vs legacy ZOMBI 1 (scaling curve)
ZOMBI1_DIR=/path/to/ZOMBI python run.py vs_zombi1_fixedtree    # opt-in: vs ZOMBI 1 on one shared tree (box-plot)
```

## Report

A written write-up of the findings (incorporating the figures) lives in
[`report/report.tex`](report/report.tex) → `report/report.pdf`. Rebuild with:

```bash
cd report && latexmk -pdf report.tex   # needs a LaTeX toolchain (pdflatex)
```

## Layout

```
analyses/performance/
├── config.py           the ONE simulation regime (model + rates), shared by all
├── perfkit/            reusable core — stable, rarely edited
│   ├── timing.py         measure(fn) → per-repeat times; the Point record
│   ├── environment.py    describe() → git / interpreter / Rust snapshot
│   ├── io.py             Result → save/load self-describing JSON
│   ├── memory.py         isolated-subprocess peak-RSS measurement
│   └── style.py          the monochrome (B&W) matplotlib house style
├── benchmarks.py       the benchmark registry — edit this to add/change tasks
├── run.py              measure → results/     (measurement only)
├── plot.py             results → figures/     (plotting only)
├── results/            raw timings (JSON, git-committable, reproducible)
└── figures/            publication figures (SVG + PNG @ 300 dpi)
```

## What is measured

ZOMBI2 runs a **single engine** by default (the built-in model is Rust; there is
no Python-vs-Rust comparison anymore), so these measure the library as a user
actually runs it, pushed to where each task stops being cheap.

| Benchmark          | Task                                          | The story |
|--------------------|-----------------------------------------------|-----------|
| `species_tree`     | Simulate a species tree vs #tips              | O(N) → tens of millions |
| `gene_families`    | Evolve D/T/L/O families vs #tips              | full log (~10⁵) vs **event trace** vs sparse counts (both → millions) |
| `memory_scaling`   | Peak RSS vs #tips (isolated subprocess)       | what fits in RAM |
| `parallel_scaling` | N independent sims vs #worker processes       | compute-bound speed-up |
| `write_output`     | Serialise a full simulation to disk (opt-in)  | reconstruct + write cost |
| `vs_zombi1`        | Same task in ZOMBI2 vs legacy ZOMBI 1 (opt-in)| scaling curve; >1000× faster; ZOMBI 1 ceiling ~1200 tips |
| `vs_zombi1_fixedtree` | Both engines' genome step on **one shared** 1000-tip tree (opt-in) | box-plot of per-run times; ≈580× faster on the identical tree (overview panel c) |

All benchmarks share one regime (`config.py`): `BirthDeath(λ=1.0, μ=0.3)`, tree
age 2.0, `D=0.2 T=0.1 L=0.25 O=0.5`, `initial_families=20`. Trees are built once per
size and only the timed call sits inside the timer. (`vs_zombi1` uses a matched
pure-birth Yule tree and rate regime so the two tools do comparable work — see
[`vs_zombi1/NOTES.md`](vs_zombi1/NOTES.md); it is opt-in and needs the ZOMBI 1
checkout via `$ZOMBI1_DIR`.) `vs_zombi1_fixedtree` goes further: it builds **one**
1000-tip tree and runs *both* engines' genome step on that identical tree (ZOMBI2 reads
it via `read_newick`), so the only variable is the engine — this removes the tree-depth
confound and is the box-plot in overview panel c.

Methodology: a monotonic clock (`time.perf_counter`), one untimed warm-up (none
for the multi-second giants), several repeats **with the cyclic GC disabled inside
each timed region** (a full `gc.collect()` between repeats) — so a collection can't
land in the timed call and inflate the object-heavy paths, exactly as `timeit` does.
Figures plot the **median** with a shaded **interquartile (25–75%)** band. Raw
per-repeat times stay in the JSON, so any estimator recomputes without re-running.

## Style

Monochrome, print-first (`perfkit/style.py`): every series is black or grey and
told apart by **line style** + **open marker shape**, never colour — reproduces
perfectly in greyscale and single-ink print.

## Scale — how big can we go?

Measured on this machine (Apple Silicon, 4 P + 6 E cores, 34 GB):

* **Species tree** is O(N) and featherweight: **3M tips in ~18 s using ~2 GB**
  (1M in ~6 s / 0.7 GB). This scales to tens of millions.
* **Gene-family counts (profiles)** are **sparse (COO) output → O(N)**:
  **1M tips in ~18 s using ~3.2 GB**, where the old *dense* N×N matrix hit a wall
  at ~100k (it would have needed ~8 TB of address space at 1M). See the sparse
  refactor of `zombi2/genomes/profiles.py`.
* **Gene-family event trace** (`output="trace"`) drops the redundant speciation rows
  (recovered from the species tree by replay) and keeps only per-family leaf counts, so
  gene trees stay reconstructable yet it **reaches 1M tips (~19 s, ~5 GB)** — barely above
  the counts-only floor, the intermediate between counts and the full log.
* **Gene-family full event log** grows with the *event* count (one Python object
  per event → ~3.5 GB at 100k) → practical to ~10⁵ tips; the heaviest path, and why
  the trace exists.

**Running 10M+ :** `python run.py --full`. Species trees, event traces and sparse
profiles all reach the millions here; the full event log is the remaining heavy path.
The workspace is cluster-ready (plain-JSON results, a `--full` profile, the regime in
one `config.py`) for pushing further on Euler (SLURM).

## Adding a benchmark

Write a function returning a `Result`, decorate it, return measured `Point`s:

```python
@benchmark("my_task")
def my_task(profile: str) -> Result:
    points = []
    for n in _grid(profile, quick=[...], standard=[...], full=[...]):
        times, result = measure(lambda: do_work(n), repeat=PROFILES[profile]["reps"])
        points.append(Point(series="my series", x=n, times=times, work={...}))
    return Result(name="my_task", title="...", x_label="...", points=points)
```

`run.py my_task` picks it up automatically; add a drawer in `plot.py` for a figure.

## Notes

* Requires the compiled Rust extension (`zombi2.rust_available()` must be `True`).
* `results/` is meant to be committed — the raw data *is* the record; `figures/`
  regenerates from it at any time.
* Parallel scaling on macOS is bounded by *spawn* (each worker re-imports zombi2)
  and the P/E-core split; Linux/Euler use *fork* and scale closer to ideal.
```
