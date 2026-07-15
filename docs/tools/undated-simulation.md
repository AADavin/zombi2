# Undated reconciliation simulation

`zombi2 tools simulate` samples gene families **directly from the ALE undated model** — it is the
**generative twin** of the [reconciliation likelihood](reconciliation-likelihood.md). Where
`tools reconcile` computes $P(\text{gene tree} \mid \text{species tree}, \text{DTL})$ under the
undated / reldated model, `tools simulate` draws gene-family histories *from* that same model. The
two are a matched pair — score one direction, generate the other — and the odds round-trip:
simulate under some odds, fit the likelihood, and its maximum lands back on the generating odds.

It is admitted as a tool (rather than a mode of `zombi2 genomes`) for the same reason its
likelihood twin is: it is a **bounded, validated computation** under one fixed model, and it is not
the dated, contemporaneous-transfer process the core simulator runs. For a forward simulation with
real time and time-consistent transfers, use [`zombi2 genomes`](../guide/genomes.md); reach
for `tools simulate` when you specifically want the *undated* model people run on real data.

## The parameters are per-branch odds, not rates

This is the whole point. The undated model has **no time**: the species tree carries no meaningful
dates, and every branch has the same event *odds*. Given duplication / transfer / loss values
$d, t, l$, one traversal step of a branch consumes exactly one of four slots,

$$\text{denom} = 1 + d + t + l,\quad p_D = \tfrac{d}{\text{denom}},\; p_T = \tfrac{t}{\text{denom}},\; p_L = \tfrac{l}{\text{denom}},\; p_S = \tfrac{1}{\text{denom}},$$

with $p_S$ the "speciate/sample" reference the others are odds against. So the `--dup/--trans/--loss`
you pass here are **dimensionless per-branch odds** — the very numbers ALEml_undated / GeneRax
report — *not* the per-unit-time δ/τ/λ of a dated simulation. There is nothing to convert and no
time axis to invent: **a cladogram with no branch lengths is a valid input** (unit branches are
assumed). A duplication or transfer puts a copy back on the *same* branch, so events stack
geometrically along a branch — exactly the structure the likelihood's dynamic program resums.

| Model | Transfers |
| --- | --- |
| `undated` (default) | a transfer may land on **any** branch (plain ALEml_undated) |
| `reldated` | a transfer may land only on a branch that **overlaps the donor in time** (needs a dated tree) — the same time-overlap rule the `reldated` likelihood uses |

## When to use

- **Benchmark reconciliation inference.** Simulate ground-truth families under known undated odds,
  run ALE / AleRax (or any reconciler) on them, and measure recovery — the model that generated the
  data is exactly the model the inference assumes.
- **Validate ALElite.** The simulated survivors, scored by [`tools reconcile`](reconciliation-likelihood.md)
  under the same odds, give a likelihood whose maximum recovers the generating odds (see
  [Validation](#validation)).
- **Produce recon-accuracy ground truth.** The reconciliations are written in ZOMBI2's native
  annotated-Newick format, so they drop straight into [`tools recon-accuracy`](recon-accuracy.md)
  as the `--truth` argument.

## Parameters

`simulate_undated(tree, model, *, n_families=100, origination="root", transfers="global",
seed=None, max_events=1_000_000)`

| Parameter | Meaning |
| --- | --- |
| `tree` | A `zombi2.tree.Tree`. A dateless cladogram is fine for `undated` (unit branches assumed); `reldated` needs real branch lengths. |
| `model` | A `UndatedDTL(dup, transfer, loss)` — the per-branch odds. |
| `n_families` | How many families to originate and evolve. |
| `origination` | `"root"` (default; every family enters on the root branch) or `"uniform"` (each family enters on a uniformly chosen branch). |
| `transfers` | `"global"` (plain undated — any branch) or `"dated"` (reldated — time-overlapping branches only). The CLI exposes these as `--model undated` / `--model reldated`. |
| `seed` | Seed for a reproducible draw. |
| `max_events` | Per-family event cap; guards against runaway families at supercritical odds (`d + t` large relative to `l`), which the undated model genuinely produces. |

## Python usage

```python
from zombi2 import read_newick
from zombi2.tools.reconciliation import UndatedDTL, simulate_undated, undated_joint_loglik

tree = read_newick("((A,B)X,(C,D)Y)root;")     # a cladogram — undated needs no dates
res = simulate_undated(tree, UndatedDTL(dup=0.2, transfer=0.1, loss=0.3),
                       n_families=200, origination="root", seed=1)

print(res.n_surviving, "survived,", res.n_extinct, "went extinct")
for recon in res.reconciliations[:3]:
    if recon.extant is not None:               # a Reconciliation: complete / extant / events
        print(recon.extant)                    # ground-truth reconciled Newick

# round-trip: the undated likelihood of the survivors is maximised at the generating odds
ll = undated_joint_loglik(res.gene_trees(), res.species_tree,
                          UndatedDTL(0.2, 0.1, 0.3), n_extinct=res.n_extinct)
```

`simulate_undated` returns an `UndatedSimResult`: `reconciliations` (one ZOMBI2
[`Reconciliation`](reconciliation-likelihood.md) per family), `gene_trees()` (the extant trees,
ready for the likelihood), `profile_rows()` (per-species copy number for each surviving family),
`event_counts`, `n_surviving` / `n_extinct`, and the raw per-family `records`.

## Command line

```bash
# sample 200 families under undated odds on a (dateless) cladogram; write the ground truth
zombi2 tools simulate -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.3 -n 200 --seed 1 -o truth/

# a quick round-trip check, no files written: prints the joint undated log-likelihood
zombi2 tools simulate -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.3 -n 200 --score

# feed the ground truth to recon-accuracy to score an inferred reconciliation against it
zombi2 tools recon-accuracy -t truth/reconciled_extant.nwk -i inferred.nwk
```

`--model` selects `undated` (default) or `reldated`; `--origination` is `root` (default) or
`uniform`; `-n/--families`, `--seed` and `--max-events` are as in the API. Without `-o` a one-line
summary is printed; `--score` adds the joint undated log-likelihood of the survivors (with the true
extinct count) under the generating odds.

## Output

With `-o DIR`, four files are written:

| File | Contents |
| --- | --- |
| `reconciled_extant.nwk` | the **survivors-only** reconciled gene tree of each surviving family, one bare Newick per line — the format [`tools recon-accuracy`](recon-accuracy.md) reads. Internal labels `branch\|EVENT` (`donor\|T>recipient` for transfers), tips `species\|gid`. |
| `reconciled_complete.nwk` | the **complete** history including `LOSS\|branch` tips — the full ground truth. |
| `reconciliation_events.tsv` | a flat event table, one row per S/D/T/L event: `family`, `event`, `species`, `recipient`, `time`, `gene`. |
| `gene_family_profiles.tsv` | the **phyletic pattern**: a `family × species` copy-number matrix over the surviving families (the classic undated observable). |

## Validation

`simulate_undated` is the exact generative dual of the undated likelihood, and it is pinned three
ways (`tests/test_undated_sim.py`):

- **The likelihood peaks at the truth.** Families simulated under known odds give a joint undated
  log-likelihood that is higher at the generating $d, t, l$ than at odds off by 2× either way, on
  every axis — `::test_likelihood_peaks_at_truth`. (A finer grid recovers the generating odds as
  the argmax on all three axes.)
- **The `pS^{(2k-1)}` oracle, as a frequency.** With $d = t = 0$ on `(A,B)`, a family survives in
  both tips with probability $p_S^3$ (root speciates, each tip is sampled) — the same closed form
  the likelihood is checked against — and the sampler reproduces it as an empirical frequency
  (`::test_matches_undated_oracle_frequency`).
- **Ground truth scores 1.0 against itself.** Every emitted reconciliation, scored through
  `recon-accuracy` against itself, has event and species-mapping accuracy 1.0 and full
  transfer-donor/recipient recovery — proving the annotation is well-formed
  (`::test_self_recon_accuracy_is_perfect`).

The draw is deterministic under `seed`, and the `max_events` guard fires on supercritical odds
(`::test_deterministic`, `::test_max_events_guard`).

## References

- Szöllősi, G. J., Tannier, E., Lartillot, N. & Daubin, V. (2013). Lateral gene transfer from the
  dead. *Systematic Biology* 62(3): 386–397. (ALE / amalgamated likelihood.)
- Morel, B., Kozlov, A. M., Stamatakis, A. & Szöllősi, G. J. (2020). GeneRax: a tool for
  species-tree-aware maximum likelihood-based gene family tree inference under gene duplication,
  transfer, and loss. *Molecular Biology and Evolution* 37(9): 2763–2774. (The `UndatedDTL` model
  simulated here.)
- Morel, B., Schade, P., Lutteropp, S., Williams, T. A., Szöllősi, G. J. & Stamatakis, A. (2022).
  AleRax: species-tree-aware maximum likelihood gene tree inference. (The undated reconciliation
  pipeline these simulations benchmark.)
- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a
  phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.
  *Bioinformatics* 36(4): 1286–1288.
