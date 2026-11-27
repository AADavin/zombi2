# Performance

Measurements for ZOMBI2 **0.64.2**, commit `047ce05fd96362439a9b9f633b13ce0cfd45f94e`,
on 27 November 2026. These are timings for the specified workloads on one laptop.

<figure markdown="span">
  ![Performance of ZOMBI2 0.64.2: species-tree scaling, genome scaling at three resolutions, and complete-gene-tree runtimes compared with AsymmeTree 2.3.1](img/performance-overview.svg)
</figure>

## Species trees

A birth–death process with birth rate 1 and death rate 0.3 was grown to each requested number
of extant tips, from 1,000 to 1,000,000. Each size used three fresh processes with seeds
1001–1003. Timing stopped when simulation returned; conversion and file writing were excluded.
Cyclic garbage collection was disabled during the timed call. One million tips took a median
of **5.90 s**; 100,000 tips took **0.524 s**.

## Genome resolutions

Each size used a shared dated tree, generated with seed `100 + N`, and three genome seeds
7001–7003. Family and gene-order genomes started with 20 families. Nucleotide genomes started
with one synthetic circular chromosome of 4,411,532 bases and 4,000 genes of 1,000 bases each.
The nucleotide input was synthetic DNA, not a biological reference genome.

Duplication, transfer, loss and origination had numerical rates 0.2, 0.1, 0.25 and 0.5.
Duplication, transfer and loss were per copy for family and gene-order runs, but per lineage
for nucleotide runs. Origination was per lineage. Gene-order and nucleotide runs also used
inversion rate 0.1 per lineage. All rates are per unit time. Species-tree generation, gene-tree
construction and writing were excluded from genome timings.

At 1,000 tips, median runtimes were **0.303 s** for family, **0.271 s** for gene order and
**11.81 s** for nucleotide genomes. At 3,000 tips, they were **2.41 s**, **1.20 s** and **34.79 s**.
The family run at 10,000 tips took **22.44 s**. Error bars show the range across three repeats.

The three curves compare different representations, initial genome sizes and event-rate scopes.
They do not compare equal event counts or equal gene turnover. In particular, the nucleotide
workload retains 50-base event extents against 1,000-base declared genes. Structural-event
boundaries cannot split a declared gene, so most short events do not alter complete genes.
ZOMBI2 warns about this configuration. These settings were retained to repeat the earlier workload.

## Shared task with AsymmeTree

Both ZOMBI2 and AsymmeTree **2.3.1** received the same dated species tree for each workload,
with root age 1 and a recorded tree checksum. Duplication, transfer and loss rates were 0.5,
0.3 and 0.6 per copy; origination was disabled. Seven workloads varied tips and initial families.
Each workload used seven paired runs with seeds 1000–1006 and alternating tool order, after a
small warm-up. Both timers included simulation and materialisation of **all complete gene trees**.
Tree generation, output traversal/counting and writing were excluded.

The ratio of AsymmeTree's median runtime to ZOMBI2's ranged from **1.14 to 1.47**. Median output
leaf counts differed by at most **1.63%**, a workload check rather than evidence that the stochastic
output distributions are equivalent. Peak memory was measured in three fresh processes per tool
and size, with 200 initial families. Median peaks were **306.78 versus 197.20 MiB** at 400 tips,
and **540.27 versus 331.59 MiB** at 800 tips, for ZOMBI2 and AsymmeTree respectively. These values
include the interpreter and imported libraries.

## Environment and reproduction

Runs used a 10-core arm64 laptop with 32 GiB RAM, macOS 26.7.1 and CPython 3.12.2. Species-tree
and genome scaling used NumPy 1.26.4. Both tools in the comparison used the same environment
with NumPy 2.5.1. Raw records and dependency freezes are under
[`assets/performance/v0.64.2/`](https://github.com/AADavin/ZOMBI2/tree/main/assets/performance/v0.64.2).
The scripts in that directory take `--source-root` and `--out` arguments. Select a checkout of
commit `047ce05fd96362439a9b9f633b13ce0cfd45f94e` explicitly; they assert version 0.64.2.

The laptop was on AC power with low-power mode disabled. Startup indexing/compilation was
allowed to settle. Each timing component began after three observations with CPU idle at least
85% and one-minute load at most 5. Components ran serially, without concurrent validation or
another analysis job from this study. System observations and power/thermal records accompany
the measurements. Absolute times depend on the computer, operating system and session.

Regenerate the figure from the raw records with `python figures/scripts/fig_performance.py`
(requires Matplotlib and NumPy). The editable SVG and PNG/PDF exports are in `assets/`.

Earlier August 2026 benchmark records, including comparisons with the legacy ZOMBI1, remain
in the parent `assets/performance/` directory. They are historical measurements and are not the
source of the current figure or its AsymmeTree comparison.
