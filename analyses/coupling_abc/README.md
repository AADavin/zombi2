# Coupling ABC analysis

Simulation-based validation of **ABC inference for the coupled (Potts) gene-family model**
(`zombi2.match_coupled` + `zombi2.cooccurrence_summary`), and the report explaining it.

## What it shows

Can the strength of gene-family coupling (`J`) be recovered from a presence/absence table by
Approximate Bayesian Computation? Uses the full **gain+loss** coupling model (`gain_coupling > 0`:
coupling lowers loss *and* biases HGT establishment). On simulated data: **yes** (error ≈ 0.056).
The analysis also teases apart *how* it works — coupling makes genes both more common and more
co-occurring, so a summary of gene *commonness* alone recovers `J` in ordinary regimes; when genes
are tuned to stay equally common regardless of coupling, only the **co-occurrence** signal recovers
it (halving the commonness-only error).

It also includes an **empirical test** on a real 43-genome eggNOG dataset (`ZOMBI2_DATA`): genes of
the same COG functional category co-occur far more than genes of different functions (p≈0.002 vs a
label-permutation null that keeps each gene's real presence pattern), led overwhelmingly by **cell
motility** (12/18 functions individually significant). The non-functional "function unknown" /
"general prediction" buckets are set aside from the start. This supports **using gene function as the
coupling blueprint** — take same-function genes as coupled groups and fit their strength.

## Regenerate

```bash
python analyses/coupling_abc/run_analysis.py          # simulations: figures 1-4 + results.json (~2-3 min)
python analyses/coupling_abc/empirical_analysis.py    # eggNOG: figures 5-6 + empirical_results.json
cd analyses/coupling_abc/report && latexmk -pdf report.tex   # -> report.pdf
```

Everything is deterministic (fixed seeds). The empirical step reads `ZOMBI2_DATA/` (the tree,
the genomes table, and the eggNOG annotations).

## Layout

- `run_analysis.py` — the simulation study; writes `results.json` and figures 1-4.
- `empirical_analysis.py` — the eggNOG / COG-category analysis; writes `empirical_results.json`
  and figures 5-6.
- `eggnog.py` — reusable parser: eggNOG-mapper annotations → a `zombi2.ProfileMatrix`
  (orthologous groups × genomes) + a `{OG → COG category}` map.
- `figures/` — `fig1` co-occurrence heatmaps, `fig2` recover-J (realistic), `fig3` summary
  comparison, `fig4` prevalence-neutral payoff, `fig5` empirical pangenome + COG test,
  `fig6` per-category cohesion.
- `report/report.tex` → `report/report.pdf` — the write-up.
