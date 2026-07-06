# Birth–death (species tree)

**What it is.** The constant-rate birth–death process: lineages speciate at rate `birth` (λ) and
go extinct at rate `death` (μ). ZOMBI2 draws it **backward** as a reconstructed process
conditioned on the number of extant tips — internal-node ages are sampled i.i.d. from the
reconstructed-process CDF and assembled by uniform coalescence (Hartmann, Wong & Stadler 2010) —
or grows it **forward** in time, keeping extinct lineages. `Yule` is the pure-birth special case
(`death = 0`).

**When to use it.** The default way to make a dated species tree to simulate genomes, traits, or
sequences along. Reach for a sibling when you need heterogeneity the constant-rate process lacks:
[episodic](../guide/species-trees.md#episodic-skyline-birthdeath) rates through time,
[ClaDS](../guide/species-trees.md#per-lineage-rates-clads) per-lineage rates, or
[diversity-dependent](../guide/species-trees.md#diversity-dependent-diversification) slow-down.

## Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `birth` (λ) | speciation rate | `1.0` |
| `death` (μ) | extinction rate | `0.3` |
| `n_tips` | condition on this many extant leaves (≥ 2) | — |
| `age` / `age_type` | tree age; `crown` (root) or `stem` (origin) | — / `crown` |

Backward mode needs an explicit `age`. Forward mode (`direction="forward"` / `--mode forward`)
retains extinct `e*` leaves and accepts `mass_extinctions` and `sampling_fraction` overlays.

## Command line

```bash
zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o run/
```

Writes `species_tree.nwk` (and `species_tree_extant.nwk` when extinct/unsampled tips exist),
`species_nodes.tsv`, and the `species_tree.log` manifest.

## Python

```python
import zombi2 as z

model = z.BirthDeath(birth=1.0, death=0.3)          # z.Yule(1.0) == BirthDeath(1.0, 0.0)
tree = z.simulate_species_tree(model, n_tips=20, age=5.0, seed=1)
tree.to_newick()
```

## Output

A timed Newick tree with the [standard node names](../contributing/conventions.md#naming)
(`n*` extant, `i*` internal, `e*` extinct, `u*` unsampled) and its `species_nodes.tsv`
(`name`, `time`, `is_leaf`, `is_extant`).

## Validation

- **Oracle.** The backward age sampler matches the analytic reconstructed-process CDF by a
  Kolmogorov–Smirnov test across the Yule, critical, and subcritical regimes —
  `test_species_bd.py::test_sample_age_matches_cdf`.
- **Statistical.** Forward Yule growth matches `E[extant] = 2·e^{λ·age}` over 400 replicates —
  `test_species_forward.py::test_yule_has_no_extinction_and_matches_theory`.
- **Determinism.** A fixed seed reproduces the Newick exactly —
  `test_species_forward.py::test_reproducible`.

## Reference

Hartmann, K., Wong, D. & Stadler, T. (2010). Sampling trees from evolutionary models.
*Systematic Biology* 59(4): 465–476.
