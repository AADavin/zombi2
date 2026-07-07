# Advanced diversification (heterogeneous rates)

These models extend the [basic birth–death](birth-death.md) process by letting the rates vary —
through time, across lineages, or at a carrying capacity — plus a mass-extinction overlay that lays
survival pulses on top of any of them. Unlike constant-rate birth–death and Yule, none of these have
a closed-form reconstructed process, so they are **forward only**: ZOMBI2 grows the tree in time and
keeps extinct lineages. Speciation still happens at rate `birth` (λ) and extinction at rate `death`
(μ), but now those rates are allowed to change.

| Model | Rates | Reach for it when |
| --- | --- | --- |
| **Episodic** (skyline) | piecewise-constant through epochs | rates changed at known times |
| **ClaDS** | per-lineage, drifting at each split | rate heterogeneity scattered across the tree |
| **Diversity-dependent** | λ declines toward a carrying capacity `K` | a diversity-dependent slowdown |
| **Clade-shift** | a scheduled shift on part of the tree | one subclade radiates or slows |
| **Mass extinctions** | a survival pulse overlaid on any of the above | a pulse of extinction at set times |

## The models

### Episodic (skyline) birth–death

Piecewise-constant rates: pass a list of `birth`/`death` values and the `shifts` ages that separate
the epochs (present → past). Use it when rates are known to have changed at particular times.

### ClaDS

Each daughter lineage inherits its parent's rate times a lognormal factor (`alpha`, `sigma`), so
speciation rates drift across the tree; `turnover` sets μ/λ. Use it for rate heterogeneity with no
single shift point (Maliet, Hartig & Morlon 2019).

### Diversity-dependent

Speciation declines as standing diversity approaches the carrying capacity `carrying_capacity` (`K`),
with constant `death`. Use it for a diversity-dependent slowdown toward a diversity limit.

### Clade-shift birth–death

Constant background rates with one or more scheduled `clade_shifts` — each an `(age, birth, death)`
that re-sets the rates on a clade from that age. Use it when a single subclade radiates or slows.

### Mass extinctions (overlay)

Not a standalone model: a list of `(age, fraction)` survival pulses layered on any forward model —
at each `age`, each lineage survives with probability `fraction`. Use it for punctual extinction
events on top of any diversification process.

## Command line

`--diversification` selects the rate process. All of these run **forward only** (`--mode forward`).

```bash
# ClaDS, forward
zombi2 species --mode forward --diversification clads \
    --birth 1 --clads-alpha 0.9 --clads-sigma 0.2 --turnover 0.1 --age 5 --seed 1 -o clads/

# diversity-dependent (carrying capacity K), forward
zombi2 species --mode forward --diversification diversity-dependent \
    --birth 2 --death 0.2 -K 50 --age 15 --seed 1 -o dd/

# episodic skyline (two epochs, rates shift at age 3), forward
zombi2 species --mode forward --birth 2 1 --death 0.2 0.4 --shifts 3 --age 6 --seed 1 -o epi/

# a clade shift and a mass-extinction pulse, forward
zombi2 species --mode forward --birth 1 --death 0.3 --clade-shift 3 2.5 0.3 \
    --mass-extinction 2 0.5 --age 6 --seed 1 -o shift/
```

## Python

Models live in `zombi2.species` (and re-export at the top level, so `zombi2.ClaDS` also works). All
require `direction="forward"`:

```python
from zombi2.species import (
    BirthDeath, EpisodicBirthDeath, ClaDS, DiversityDependent, CladeShiftBirthDeath,
    simulate_species_tree,
)

# episodic (skyline): two epochs, rates shift at age 3 (present -> past)
epi = simulate_species_tree(
    EpisodicBirthDeath([2.0, 1.0], [0.2, 0.4], shifts=[3.0]),
    age=6.0, direction="forward", seed=1)

# ClaDS: per-lineage rates drift at every split
clads = simulate_species_tree(
    ClaDS(1.0, alpha=0.9, sigma=0.2, turnover=0.1),
    age=5.0, direction="forward", seed=1)

# diversity-dependent: speciation declines toward carrying capacity K
dd = simulate_species_tree(
    DiversityDependent(2.0, death=0.2, carrying_capacity=50),
    age=15.0, direction="forward", seed=1)

# clade-shift: one subclade switches rates at age 3
shift = simulate_species_tree(
    CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(3.0, 2.5, 0.3)]),
    age=6.0, direction="forward", seed=1)

# a mass-extinction overlay: (age, survival-fraction) pulses on any forward model
me = simulate_species_tree(
    BirthDeath(1.0, 0.3, mass_extinctions=[(2.0, 0.5)]),
    age=6.0, direction="forward", seed=1)
```

## Output

`species_tree.nwk` (the complete tree, extinct/ghost tips kept), and — when dead tips are present —
`species_tree_extant.nwk` pruned to the sampled-extant leaves. `species_nodes.tsv` gives per-node
metadata (`name`, `time`, `is_leaf`, `is_extant`), and `species_tree.log` is the run manifest.
Leaves use the [standard node names](../contributing/conventions.md#naming): `n*` extant, `e*`
extinct, `u*` unsampled, `i*` internal.

## Validation

- **Oracle.** A single epoch of the skyline model matches the closed-form CDF
  (`test_episodic.py::test_episodic_single_epoch_matches_analytic_cdf`); diversity-dependent growth
  saturates at exactly `K` when `death 0`
  (`test_species_forward.py::test_diversity_dependent_saturates_at_K`).
- **Statistical / reduction.** A single skyline epoch reduces to constant-rate birth–death
  (`test_species_forward.py::test_episodic_single_epoch_matches_constant`); a clade shifted to a
  faster regime raises the tip count
  (`test_species_forward.py::test_clade_shift_to_fast_regime_raises_tip_count`); a severe pulse
  lowers the extant count (`test_species_forward.py::test_mass_extinction_reduces_extant_tips`).
- **Invariant.** ClaDS trees are binary and reach the requested age with extant tips at the present,
  and `turnover 0` produces no extinction
  (`test_species_forward.py::test_clads_age_mode_complete_tree`,
  `::test_clads_turnover_zero_has_no_extinction`); mass-extinction victims die exactly at the pulse
  age (`test_species_forward.py::test_mass_extinction_kills_at_the_pulse_time`).
- **Determinism.** A fixed `seed` reproduces the Newick exactly
  (`test_species_forward.py::test_reproducible`).

## Not yet implemented

A short menu of sampling extensions that fit the existing forward/episodic machinery but are not
built yet:

- **Time-varying sampling (BDST).** Piecewise sampling `ψ(t)` conditioned on *sampled* tips —
  important for phylodynamics (epidemic trees) and variable fossil preservation. Fits the existing
  episodic machinery (a numerically inverted piecewise CDF) (Stadler et al. 2013).
- **Skyline (episodic) FBD** and the **occurrence birth–death process (OBDP)** — combined-evidence
  extensions of the fossilized birth–death already shipped (Andréoletti et al. 2022).

## References

- Maliet, O., Hartig, F. & Morlon, H. (2019). A model with many small shifts for estimating
  species-specific diversification rates. *Nature Ecology & Evolution* 3: 1086–1092.
- Etienne, R. S. et al. (2012). Diversity-dependence brings molecular phylogenies closer to
  agreement with the fossil record. *Proceedings of the Royal Society B* 279: 1300–1309.
