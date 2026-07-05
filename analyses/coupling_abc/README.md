# Coupling ABC analysis

Simulation-based validation of **ABC inference for the coupled (Potts) gene-family model**
(`zombi2.match_coupled` + `zombi2.cooccurrence_summary`), and the report explaining it.

## What it shows

Can the strength of gene-family coupling (`J`) be recovered from phylogenetic profiles by
Approximate Bayesian Computation? On simulated data: **yes** (RMSE ≈ 0.056). The analysis also
teases apart *how* it works — coupling is entangled with prevalence, so the marginal frequency
spectrum recovers `J` in ordinary regimes; when prevalence is held fixed the marginal goes blind
and **gene-family co-occurrence (module) structure** becomes essential, roughly halving the error.

It also includes an **empirical test** on a real 43-genome eggNOG dataset (`ZOMBI2_DATA`):
**17 of 18 functional COG categories are significant co-occurrence modules** (cell motility most of
all), and same-category groups co-occur more than different-category ones (p≈0.002 vs a
phylogeny-controlled null) once the "function unknown" catch-all (S, ~60% of variable groups, not a
module) is excluded. This supports **scaffolding the coupling structure J on *functional* COG
categories** — turning "infer which families couple" into a principled prior.

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
