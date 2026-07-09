# ZOMBI2 × Snakemake — the RED node-age benchmark (parameter sweep)

A reproducible, parallel, cluster-ready **Snakemake** workflow built on ZOMBI2. Its flagship
analysis is the **RED node-age benchmark**: does *Relative Evolutionary Divergence* (Parks et al.
2018, the measure GTDB uses to normalize taxonomic ranks) recover a tree's true node ages after
its branch rates are perturbed across lineages? ZOMBI2 knows the true ages, so the ground truth is
built in. It reproduces the single-tree result of Rinke et al. (2021, GTDB) and generalizes it into
a full sweep over tree size, extinction, clock model, and perturbation strength.

It is a shipped example — copy it as the template for your own ZOMBI2 parameter sweeps. See the
[documentation page](../../docs/examples/red-benchmark.md) for the narrated version.

## Pipeline (the Snakemake DAG)

1. **`simulate_timetree`** — a forward Yule/birth-death tree, pruned to the extant (ultrametric)
   tree, with known node ages. Simulated once per (tree config, seed), reused downstream.
2. **`perturb_and_red`** — apply a relaxed clock (time → substitution branch lengths), compute RED,
   and compare RED-recovered ages to the truth → one-row `metrics.tsv` (Pearson/Spearman r, nRMSE,
   rate spread) + `points.csv` (per-node true vs recovered).
3. **`summarize`** — aggregate every run into `results/summary.tsv` and two figures.

## Run

```bash
pip install "zombi2[bench]"                # snakemake + matplotlib (+ SLURM executor plugin)

snakemake --cores 8                        # full sweep  (config/sweep.yaml)
snakemake --cores 4 --config cfg=config/test.yaml   # tiny smoke grid (CI)
snakemake -n                               # dry run — inspect the DAG
snakemake --profile workflow/profiles/slurm --config cfg=config/sweep.yaml   # Euler/SLURM

pytest tests/                              # unit tests for the RED core (no Snakemake needed)
```

The grid is one config file, chosen with `--config cfg=<path>` (loaded directly, so a test grid
*replaces* the default rather than merging with it). Edit `config/sweep.yaml` freely — `trees`
(size / birth-death / seeds), `perturbations` (clock + params), `replicates`, and the `showcase`
tree for the scatter figure.

## Result (config/sweep.yaml, 245 runs, ~10 s on 10 cores)

| clock | across-lineage spread (95/5) | Pearson r | nRMSE |
|---|---|---|---|
| strict (anchor) | 1× | **1.0000** | 0.00% |
| ratevar (GTDB bins) | 1.5–4.7× | 0.999 → 0.991 | ≤ 3% |
| cir | ~3.8× | ~0.986 | ~3.5% |
| aln (heavy-tailed) | 19–3500× | 0.98 → 0.90 | up to ~16% |

**RED robustly recovers node ages under bounded / moderate across-lineage rate variation** (r ≥
0.99, the Rinke/GTDB result) and **degrades only under extreme heavy-tailed variation** — the
generalization the single-tree script could not show. Two invariants are asserted every run: RED on
the time tree equals `node.time / total_age` to machine precision, and the strict clock recovers
ages exactly (r = 1.000).

`results/figures/red_scatter.png` (true vs recovered, showcase tree) and
`results/figures/red_curve.png` (r & nRMSE vs spread, per tree × clock).

## Layout

```
config/            sweep.yaml (default grid), test.yaml (CI smoke)
workflow/
  Snakefile        config load, grid expansion, `rule all`
  rules/           simulate.smk, analyze.smk, summarize.smk
  scripts/         red.py (clock factory, time-tree sim, RED-recovery metrics),
                   simulate_timetree.py, perturb_and_red.py, summarize.py
  profiles/slurm/  Euler/SLURM profile (fill in the account)
tests/             test_red.py — RED invariants, clock factory, determinism
results/           outputs (git-ignored)
```

RED itself is the shipped tool: `red.py` calls `zombi2.tools.relative_evolutionary_divergence`
(the same code behind `zombi2 tools red`) — the example only orchestrates and measures, and treats
ZOMBI2 as read-only.
