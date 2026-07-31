# Performance

Timings for the two operations that dominate a run: growing a species tree, and evolving genomes
along it.

<figure markdown="span">
  ![ZOMBI2 performance overview: (a) species-tree simulation scaling to millions of tips; (b) genome simulation at the family, ordered and nucleotide resolutions; (c) ZOMBI2 about 183 times faster than the legacy ZOMBI v1 on one shared 1,000-tip tree](img/performance-overview.svg)
</figure>

## What was measured

Every benchmark uses one regime: a birth–death species tree (λ = 1.0, μ = 0.3) grown to *N* extant
tips, then genomes under D = 0.2, T = 0.1, L = 0.25, O = 0.5 from 20 initial families.

**(a) Species trees.** Linear in the number of tips: 0.96 s at 100,000 tips, 9.0 s at 1,000,000
(peak memory 1.9 GB).

**(b) Genomes, at the three resolutions.** `family` takes 0.33 s at 1,000 tips and 30 s at 10,000.
`ordered` — the same rates plus inversion = 0.1 — takes 1.2 s at 1,000 tips and 10 s at 3,000.
`nucleotide`, run on an *M. tuberculosis* H37Rv genome (4.41 Mbp), takes 11 s at 1,000 tips and 31 s
at 3,000. The three share rate *values*, not event counts: family and ordered rates are per gene
copy, nucleotide rates are per lineage. The panel therefore compares cost at matched rates, not at
matched work.

Family and ordered are super-linear in tip count. A tree grown to *N* extant tips deepens as ln *N*,
so genome work grows faster than *N*.

**(c) ZOMBI2 against the legacy ZOMBI v1.** One shared 1,000-tip species tree, both engines running
their genome step on it, 10 runs each. The ZOMBI2 median is 0.25 s against ZOMBI v1's 44.9 s —
**183× faster**. Both are pure Python, so this measures the rewrite, not a change of language.

## Provenance

Measured on 2026-07-30 with ZOMBI2 0.18.0, CPython 3.12.2 and NumPy 1.26.4, on macOS / Apple Silicon
(10 cores). Each point is the median of five timed repeats — fewer at the largest sizes — with the
cyclic garbage collector disabled inside each timed call. Memory is peak resident set size from an
isolated subprocess, one measurement per size. Absolute times are machine-specific; the scaling is
not.
