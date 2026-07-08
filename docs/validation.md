# Validation

A simulator is only as useful as the ground truth it produces. ZOMBI2 treats validation as a
**first-class, visible guarantee**: models are checked against theory, and that check is part of
the code, run on every change.

!!! quote "The rule"
    **No model enters the core without an oracle or a statistical test.** A test that only
    asserts "it runs without error" is not validation. This is enforced in review — see
    [Adding a model](contributing/adding-a-model.md).

## How models are checked

Validation comes in four kinds, strongest first. Most models carry several.

### Oracle — against a closed form

Where the model has a known analytic expectation, many seeded replicates are compared to it.

| Model | Checked against | Test |
|---|---|---|
| Brownian motion | tip mean `x₀ + trend·depth` and covariance `σ²·C` | `test_traits.py::test_bm_tip_moments_match_theory` |
| Ornstein–Uhlenbeck | reversion mean and variance `σ²/2α·(1−e^{−2αt})` | `test_traits.py::test_ou_transition_moments_match_theory` |
| Mk (equal rates) | transition matrix `P(t)` closed form, to `1e-9` | `test_traits.py::test_mk_equal_rates_transition_closed_form` |
| Birth–death | reconstructed internal-age CDF, via a KS test | `test_species_bd.py::test_sample_age_matches_cdf` |
| JC69 sequences | Jukes–Cantor distance recovered from 150k sites | `test_sequence_sim.py::test_jc_distance_recovered` |
| HKY85 sequences | stationary base frequencies recovered to ±0.01 | `test_sequence_sim.py::test_stationary_frequencies_recovered` |
| LG/WAG/JTT/Dayhoff | published amino-acid frequencies, to `1e-4` | `test_sequence_sim.py::test_empirical_aa_frequencies_match_published` |
| Gene conversion | mean within-family coalescence depth `1/(2c)` on a stable two-copy family | `test_gene_conversion.py::test_conversion_homogenizes_coalescence_depth_matches_theory` |

### Statistical — reduction to a known case

Where there is no closed form, the model is shown to reduce to one that has, within a
Monte-Carlo tolerance:

- A single epoch of `EpisodicBirthDeath` matches constant-rate `BirthDeath`
  (`test_species_forward.py::test_episodic_single_epoch_matches_constant`).
- A state-independent SSE recovers the Yule mean `2·e^{λt}`
  (`test_sse.py::test_sse_reduces_to_yule_mean_count`), and a faster-speciating state biases
  the tips toward it (`test_sse.py::test_sse_faster_speciation_biases_tips`).
- The **Rust and pure-Python** gene-family engines agree on mean family count within 15 %
  (`test_rust.py::test_rust_matches_python_engine`) — the compiled fast path is held to the
  reference implementation.
- **Transfer receptivity** lands about the same share of transfers on a boosted branch under both
  the Rust and Python engines, and a run with no receptivity is byte-identical to a plain run
  (`test_rate_controls.py::test_rust_and_python_agree_on_receptivity`,
  `::test_receptivity_off_is_byte_identical`) — recipient choice is also a frequency oracle:
  candidates are chosen in proportion to their weight (`::test_receptivity_makes_selection_proportional_to_weight`).

### Invariants — structural laws that must hold

- **Reconciliation**: every family's extant gene-tree leaf count equals its extant copy number
  (`test_rust.py::test_full_log_gene_tree_invariant`).
- **Stochastic maps tile their branch**: segment durations sum to the branch length and states
  stay continuous across nodes (`test_traits.py::test_mk_history_tiles_branches_and_matches_values`).
- **Genome assembly**: concatenating a leaf's blocks by its mosaic reconstructs the leaf genome
  exactly, and a zero-divergence run returns the input genome
  (`test_sequence_sim.py`, `test_nucleotide_rust.py`).

### Determinism

Every model has a test that a fixed `seed` reproduces its output exactly — the guarantee behind
every published dataset.

## Tolerances

Thresholds are chosen for the kind of claim, in the spirit of the whole suite:

| Claim | Tolerance |
|---|---|
| Exact identity (transition matrix, age match, assembly) | `1e-9` |
| Recovered frequency / distance | ~1 % |
| Monte-Carlo moment (mean, variance, cross-engine) | 5–15 % |

## Continuous integration

The full suite (**over 1,200 tests**) runs on every push and pull request, on Python 3.10, 3.11,
and 3.12. The compiled engine is a required part of that run: CI asserts `zombi2.rust_available()`
and **fails rather than skips** if it is missing, so the Rust path can never pass untested. A
strict documentation build runs alongside.
