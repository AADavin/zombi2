# Reconciliation likelihood (ALElite)

**ALElite** computes the **marginal reconciliation likelihood** of a simulated gene family:

$$P(\text{gene tree} \mid \text{species tree}, \text{DTL rates})$$

— the quantity [ALE](https://github.com/ssolo/ALE) (Amalgamated Likelihood Estimation) reports.
In ZOMBI2 the gene tree and species tree are *exact* (no gene-tree uncertainty / no CCPs), so
the "amalgamation" over a sample of gene trees collapses to one tree at probability 1, and the
likelihood is a plain dynamic-programming sum over **every reconciliation** of that one gene
tree against the species tree, marginalising losses analytically.

## What it computes

Given the **extant** (survivors-only) reconciled gene tree of a family and duplication /
transfer / loss rates, ALElite returns the log of `P(gene tree | species tree)` under one of
three ALE models:

| Model | Transfers | Rates | Notes |
| --- | --- | --- | --- |
| `dated` (default) | only to lineages alive at the **same instant** (time-sliced) | per-unit-time δ/τ/λ | The Szöllősi-2013 dated model — **faithful to ZOMBI2's simulator** (`_choose_recipient` draws a contemporaneous recipient). Rates are ZOMBI2's native units, directly comparable to what a run was generated under. |
| `undated` | to **any** branch | per-branch odds | The ALEml_undated / GeneRax `UndatedDTL` model — what people run on real data. Fast, but not time-consistent. |
| `reldated` | only to branches that **overlap the donor in time** | per-branch odds | The middle ground: undated per-branch odds, but transfers restricted to contemporaneous branches using the species-tree dates (no full slicing). |

The likelihood is a **probability**, not a density: the returned value is always `≤ 0` (log of
a number in `(0, 1]`), and `-inf` for a reconciliation impossible under the given rates (e.g. a
gene tree discordant with the species tree when `transfer = 0`).

Losses are marginalised via the extinction probability, so ALElite consumes the **extant** tree
— never the complete tree with its `LOSS…` tips (feeding the complete tree is refused with a
clear error).

## When to use

- **Score a simulated scenario under its generating model** — the log-likelihood ALE inference
  would assign in *perfect* conditions (perfect species tree, perfect gene trees). Use `dated`
  with the simulation's own δ/τ/λ.
- **Compare models or rate settings** on a fixed simulated gene tree (e.g. how much a transfer
  rate raises the likelihood of a discordant family).
- **Build a test bed for reconciliation inference** — a ground-truth likelihood to validate an
  external tool against.

It is **not** an inference method: it evaluates the likelihood at rates you supply; it does not
*fit* rates. Simulation-based inference (ABC/MCMC) is [deferred to a future release](index.md).

## Parameters

`reconciliation_likelihood(reconciliation=None, species_tree=None, *, gene_tree=None,
duplication=0.0, transfer=0.0, loss=0.0, model="dated", origination="root", n_steps=100)`

| Parameter | Meaning |
| --- | --- |
| `reconciliation` | A ZOMBI2 `Reconciliation` (one value of `Genomes.reconciliations()`); its **extant** tree is scored. Pass this *or* `gene_tree`. |
| `gene_tree` | An explicit gene tree: a `GeneTree`, or a reconciled-extant Newick string (leaf labels `"<species>\|<gid>"`). |
| `species_tree` | The extant/reconstructed species tree — a `zombi2.tree.Tree` or a prebuilt `SpeciesTree`. |
| `duplication`, `transfer`, `loss` | DTL rates. For `dated`: per-unit-time δ/τ/λ. For `undated`/`reldated`: per-branch odds. |
| `model` | `"dated"` (default), `"undated"`, or `"reldated"`. |
| `origination` | Where the family enters the tree: `"root"` (default; exact for ZOMBI2's root-seeded families) or `"uniform"` (root gene node averaged over all branches). |
| `n_steps` | `dated` only: per-slice time-grid resolution (raise to check convergence; default 100). |

## Python usage

The one-call entry point on `zombi2.tools`:

```python
import zombi2 as z
from zombi2.tools import reconciliation_likelihood

# a small simulated scenario
tree = z.simulate_species_tree(z.Yule(1.0), n_tips=6, age=2.0, seed=7)
genomes = z.simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.15,
                             origination=0.0, initial_families=8, seed=7)

# score every extant family under the (faithful) dated model at the generating rates
for family, recon in genomes.reconciliations().items():
    if recon.extant is None:        # fully-extinct family: no observable gene tree
        continue
    ll = reconciliation_likelihood(recon, tree, model="dated",
                                   duplication=0.1, transfer=0.05, loss=0.15)
    print(family, round(ll, 4))     # a finite log-likelihood <= 0
```

To score many families against the same background efficiently, bind the settings once with the
**`ReconciliationLikelihood`** class (it builds the species-tree index and extinction background
once):

```python
from zombi2.tools import ReconciliationLikelihood

scorer = ReconciliationLikelihood(tree, duplication=0.1, transfer=0.05, loss=0.15,
                                  model="dated")
rows = scorer.score_all(genomes.reconciliations())   # list[FamilyScore]
```

You can also drive the underlying engines directly (`undated_loglik`, `dated_loglik`,
`reldated_loglik`, and the batched `*_joint_loglik` / `score_reconciliations`) — see
`zombi2.tools.reconciliation`.

### Command line

Two entry points, for the two situations.

#### Score gene trees you already have — `zombi2 tools reconcile`

Given a gene tree and a species tree (from a ZOMBI2 run, or anywhere), evaluate the
reconciliation log-likelihood at a **given** set of DTL rates — the direct CLI counterpart of
`reconciliation_likelihood()`. It *evaluates* the likelihood at the rates you pass; it does
**not** fit or optimise them.

```bash
# one log-likelihood, printed to stdout (dated model)
zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15

# all three models, written into out/Reconciliation_likelihoods.tsv
zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk \
    --dup 0.1 --trans 0.05 --loss 0.15 --model dated undated reldated -o out/
```

`-g` is a Newick file of one or more reconciled gene trees (one per line; tip labels
`<species>|<gid>`, the species matching a leaf of the dated `-t` species tree). With a single
tree and a single model the bare log-likelihood is printed — handy for scripting,
`ll=$(zombi2 tools reconcile …)`; otherwise a `family / extant_copies / <model>_loglik…` table
is printed, and `-o DIR` writes that table as `DIR/Reconciliation_likelihoods.tsv` (the same
file `genomes --score-likelihoods` produces). Flags mirror the API: `--model` (one or more of
`dated undated reldated`, default `dated`), `--n-steps` (dated grid resolution),
`--origination` (`root`|`uniform`).

#### Score a whole simulation in one pass — `genomes --score-likelihoods`

When you are *simulating* the families anyway, `genomes` can score every extant family it
generates in the same run:

```bash
zombi2 genomes --tree species_tree.nwk \
    --dup 0.1 --trans 0.05 --loss 0.15 --initial-families 8 --seed 7 \
    --write trees --score-likelihoods --score-model dated undated -o out/
```

This writes `out/Reconciliation_likelihoods.tsv`. Flags: `--score-model` (one or more of
`dated undated reldated`, default `dated undated`), `--score-nsteps` (dated grid resolution),
`--score-origination` (`root`|`uniform`). Scoring forces the full gene-family path (the
counts-only / trace fast paths do not reconstruct gene trees).

## Output

- **Python** — `reconciliation_likelihood(...)` returns a `float` (log-likelihood). The
  per-family helpers return `FamilyScore(family, extant_tips, logliks)` rows, where `logliks`
  maps each model name to its log-likelihood.
- **CLI / `write_scores_tsv`** — `Reconciliation_likelihoods.tsv`, one row per extant family:
  `family`, `extant_copies`, and a `<model>_loglik` column for each scored model.
- **CLI `tools reconcile`** — prints the bare log-likelihood to stdout for a single tree and a
  single model, else the same table; `-o DIR` writes it as `DIR/Reconciliation_likelihoods.tsv`.

## Validation

The undated model has exact closed-form limits, so the DP is checked against **hand-derived
probabilities** (an oracle), not merely against a reference binary:

- **Matching tree, no D/T.** A gene tree perfectly matching a `k`-tip species subtree with
  `d = t = 0` has probability `pS^(2k−1)` (one slot per speciation and per tip sampling, nothing
  to lose) — `test_tools_alelite.py::test_matching_tree_is_pS_to_the_2k_minus_1`.
- **Gene in one of two sisters.** With `d = t = 0` and the gene present in only one of two sister
  species, `P = l/(1+l)^3` (root speciates, present copy sampled, absent copy lost) —
  `test_tools_alelite.py::test_gene_present_in_one_of_two_species`.

The **dated** engine is pinned against **birth–death closed forms** (the `τ = 0` limit) and, most
decisively, by **inject-recover**:

- **Extinction → birth–death.** With `τ = 0` the coupled extinction ODE reduces to the
  birth–death extinction probability —
  `test_tools_alelite.py::test_dated_extinction_matches_birth_death`.
- **Single gene / matching pair → birth–death.** `P = p₁(s)·E_sister` and `P = p₁(s)²`
  respectively — `test_tools_alelite.py::test_dated_single_gene_matches_birth_death_closed_form`
  and `::test_dated_matching_pair_matches_birth_death_closed_form`.
- **Inject-recover (the decisive test).** Families simulated under known δ/τ/λ on a Yule tree
  (no species extinction ⇒ ZOMBI2's contemporaneous transfers match the dated model exactly)
  give a joint dated log-likelihood that is **higher at the true rates** than at rates off by
  ~2.5× either way — `test_tools_alelite.py::test_dated_inject_recover_prefers_true_rates`.

The optional Rust fast path is checked **bit-for-bit** against the pure-Python reference for all
three models (`::test_dated_rust_matches_python`, `::test_undated_reldated_rust_matches_python`).
Every extant ZOMBI2 family is confirmed to score to a finite log-likelihood `≤ 0`
(`::test_zombi_reconciled_gene_tree_scores_finite`).

## References

- Szöllősi, G. J., Rosikiewicz, W., Boussau, B., Tannier, E. & Daubin, V. (2013). Efficient
  exploration of the space of reconciled gene trees. *Systematic Biology* 62(6): 901–912.
  (The dated, time-sliced DTL reconciliation model.)
- Szöllősi, G. J., Tannier, E., Lartillot, N. & Daubin, V. (2013). Lateral gene transfer from the
  dead. *Systematic Biology* 62(3): 386–397. (ALE / amalgamated likelihood.)
- Morel, B., Kozlov, A. M., Stamatakis, A. & Szöllősi, G. J. (2020). GeneRax: a tool for
  species-tree-aware maximum likelihood-based gene family tree inference under gene duplication,
  transfer, and loss. *Molecular Biology and Evolution* 37(9): 2763–2774. (The `UndatedDTL`
  model implemented here as `undated`.)
- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a
  phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.
  *Bioinformatics* 36(4): 1286–1288.
