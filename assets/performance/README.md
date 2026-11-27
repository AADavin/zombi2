# Performance results

The current README and documentation figure uses the **27 November 2026 ZOMBI2 0.64.2 study**
in [`v0.64.2/`](v0.64.2/). It includes species-tree scaling, genome scaling at three resolutions,
and the shared complete-gene-tree comparison with AsymmeTree 2.3.1. Each dataset includes raw
repeats and provenance. Timing protocols and qualifications are documented on the
[performance page](https://aadavin.github.io/ZOMBI2/docs/performance/).

Regenerate the current figure with `python figures/scripts/fig_performance.py`.

## Historical August 2026 records

The files below are retained from the earlier study. They are not the source of the current
performance figure or its AsymmeTree comparison.

Measured on 2026-08-16, on ZOMBI2 main at commit `fcfe57f` (CPython 3.12.2, NumPy 1.26.4), on
macOS / Apple Silicon (10 cores), each benchmark run serially on an otherwise idle machine.
Absolute times are machine-specific; the scaling is not.

| File | What it holds |
|---|---|
| `species_tree.json` | Wall-clock to grow a birth–death species tree vs number of extant tips |
| `gene_families.json` | Genome (D/T/L/O) simulation vs tip count, at the family, ordered and nucleotide resolutions |
| `memory_scaling.json` | Peak resident memory vs run size, one isolated subprocess per point |
| `parallel_scaling.json` | Wall-clock of a batch of independent runs vs worker-process count |
| `write_output.json` | Time to reconstruct every gene tree and write a run's output folder |
| `vs_zombi1.json` | ZOMBI2 vs the legacy ZOMBI v1 on the same gene-family task, over tip counts |
| `vs_zombi1_fixedtree.json` | The same head-to-head at a fixed size, on one shared 1,000-tip tree |
