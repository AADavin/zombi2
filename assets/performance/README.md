# Performance results

The raw timing records behind the README's performance figure and the docs site's
[Performance page](https://aadavin.github.io/zombi2/docs/performance/): one JSON file per benchmark,
holding the raw per-repeat times point by point, with a provenance header (library version,
interpreter, NumPy, platform, timestamp) embedded in each file.

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
