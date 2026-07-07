# Running many replicates in parallel

A single ZOMBI2 simulation is a sequential process, but **replicate** simulations
(independent seeds) are embarrassingly parallel — and Monte-Carlo studies usually want
many. `run_replicates` runs them across CPU cores: each worker runs one full simulation and
**writes its own output to disk**, so the millions of log records never travel back between
processes.

```python
from zombi2.species import BirthDeath
from zombi2.genomes import run_replicates

summaries = run_replicates(
    100, "runs/",                         # 100 replicates -> runs/replicate_0000/, ...
    BirthDeath(birth=1.0, death=0.3),     # species-tree model
    n_tips=50, age=8.0,
    duplication=0.1, transfer=0.05, loss=0.15, origination=0.5,  # or rates=SharedRates(...)
    initial_families=100, max_family_size=0.5,
    seed=42,                              # base seed -> reproducible, independent per replicate
    processes=None,                       # default: all cores; processes=1 = serial
)

for s in summaries:                       # small per-replicate summaries come back
    print(s["replicate"], s["seed"], s["n_families"], s["n_events"], s["path"])
```

Each replicate `i` is written to `outdir/replicate_<i>/` (the full ZOMBI1-style output —
trees, event tables, transfers, profiles) with an independent seed derived from `seed`, so
the batch is **reproducible and independent of the number of processes** (serial and
parallel give identical results).

The species-tree and gene-family parameters mirror
[`simulate_species_tree`](species-trees.md) and [`simulate_genomes`](gene-families.md)
(`rates`, `transfers`, `max_family_size`, `genome_factory`, …).

!!! warning "Two multiprocessing rules"
    The default start method is *spawn* (macOS/Windows), so:

    1. **Call `run_replicates` under an `if __name__ == "__main__":` guard** in your script.
    2. **Every argument must be picklable** — use the built-in models/distributions, and
       `functools.partial(OrderedGenome, extension=0.5)` rather than a `lambda` for a
       `genome_factory`.

## When it helps

Speedup approaches the core count once each replicate does real work; for very short
simulations the pool startup and per-replicate disk writes dominate. Rough example (8
replicates of 500-tip genomes on a 10-core machine): ~13.9 s serial → ~3.6 s parallel.
