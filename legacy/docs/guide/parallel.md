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
    duplication=0.1, transfer=0.05, loss=0.15, origination=0.5,  # or rates=Rates(...)
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
[`simulate_species_tree`](species-trees.md) and [`simulate_genomes`](genomes.md)
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

## Running on a cluster (SLURM)

`run_replicates` parallelises across the cores of **one machine**. To scale a sweep across a
cluster, drive ZOMBI2 from a workflow manager and let the scheduler place the jobs — the shipped
`examples/red_benchmark/` **Snakemake** workflow is a copy-paste template that does exactly this.
Two patterns cover most needs:

- **Snakemake + the SLURM executor** (best for sweeps). Each rule instance becomes its own
  `sbatch` job. Install the executor with `pip install "zombi2[bench]"` and run with a SLURM
  profile:

  ```bash
  snakemake --profile workflow/profiles/slurm --config cfg=config/sweep.yaml
  ```

  The profile (`examples/red_benchmark/workflow/profiles/slurm/config.yaml`) sets `executor:
  slurm`, the max concurrent `jobs`, and per-rule `runtime` / `mem_mb` / `cpus_per_task`. Fill in
  your `slurm_account` (the `CHANGE_ME` placeholder) and, only if your cluster requires one,
  `slurm_partition`.

- **An `sbatch` array over seeds** (simplest for one config, many replicates). Submit an array job
  and map the array index to a seed — e.g. `zombi2 species … --seed $SLURM_ARRAY_TASK_ID`. Because
  every ZOMBI2 run is fully determined by its seed, replicate *i* is reproducible and independent of
  where or when it runs. This composes with `run_replicates`' per-replicate seeding: give each array
  task a **disjoint base seed** so the derived seeds never collide.

Keep the per-job work meaningful — many sub-second jobs are dominated by scheduler latency and
disk writes, so group tiny replicates into one job rather than submitting one job per replicate.
