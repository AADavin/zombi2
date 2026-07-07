# Continuous-trait models

A quantitative trait — a body size, an expression level, a latent liability — evolves down the
tree as a **diffusion**: along each branch it performs a random walk, and ZOMBI2 integrates that
walk **exactly**, node by node in preorder from a root value, so the internal nodes are true
ancestral states rather than an inference. The models below are variants of that diffusion: a
free random walk (BM), a pull toward an optimum (OU), a rate that changes through time
(EarlyBurst), a per-clade optimum (MultiOptimumOU), and vector-valued versions that let several
traits evolve in a correlated way (MultivariateBrownian, MultivariateOU). All are exact —
every one reproduces the closed-form tip law, so tips come out multivariate-normal.

| Model | Diffusion | Reach for it when |
| --- | --- | --- |
| **BM** | free random walk (optional trend) | the default null model for a neutral quantitative trait |
| **OU** | pulled toward one optimum `theta` | stabilizing selection / adaptation to a fixed optimum |
| **EarlyBurst** | rate changes exponentially through time | an adaptive radiation (early burst) or an accelerating rate |
| **MultiOptimumOU** | a different optimum on each painted regime | different clades adapt toward different optima |
| **MultivariateBrownian** | vector trait, rate matrix `R` | several traits evolve together, correlated |
| **MultivariateOU** | vector trait pulled toward `theta` | correlated traits under multivariate stabilizing selection |

## The models

### Brownian motion (BM)

The classic model (Felsenstein 1985): along a branch `dX = trend·dt + σ·dW`, so the endpoint
from `x` is normal with mean `x + trend·t` and variance `sigma2·t`. Simulated in preorder this
reproduces the exact tip law — tips are multivariate-normal with mean `x0 + trend·depth` and
covariance `sigma2 · C`, where `C` is the tree's shared-path-length matrix. Parameters: `sigma2`
(diffusion rate, ≥ 0), `x0` (root value, default `0.0`), `trend` (directional drift, default
`0.0`). With `trend = 0` it is the driftless random walk; a non-zero `trend` biases it. The
default null for a trait wandering neutrally down a tree.

### Ornstein–Uhlenbeck (OU)

A trait under stabilizing selection: it diffuses but is pulled toward an optimum `theta` (θ) with
strength `alpha` (α), `dX = alpha·(theta − X)·dt + σ·dW` (Hansen 1997; Butler & King 2004). Over a
branch of duration `t` the endpoint from `x` is normal with mean `theta + (x − theta)·e^{−alpha·t}`
and variance `sigma2/(2·alpha)·(1 − e^{−2·alpha·t})`, so tips settle into the stationary law
`N(theta, sigma2/(2·alpha))`. Parameters: `sigma2` (default `1.0`), `alpha` (> 0), `theta`, and
`x0` (defaults to `theta` — start at the optimum). `alpha` must be strictly positive — use BM for
`alpha = 0`. Reach for it when a trait is constrained rather than free to wander.

### Early burst (EarlyBurst)

Brownian motion whose rate changes exponentially through time (Blomberg et al. 2003; Harmon et al.
2010): the diffusion rate at absolute time `t` (root at 0) is `σ²(t) = sigma2 · e^{rate·t}`. With
`rate < 0` the rate **decays** — most divergence happens early, the signature of an adaptive
radiation (an *early burst*); with `rate > 0` it **accelerates** (the AC of ACDC); `rate = 0` is
plain BM. The variance over a branch spanning `[t1, t2]` is the exact integral
`sigma2·(e^{rate·t2} − e^{rate·t1})/rate`, so tips stay multivariate-normal. Parameters: `sigma2`
(rate **at the root**, ≥ 0), `rate`, `x0` (default `0.0`), `trend` (default `0.0`).

### Multi-optimum OU (MultiOptimumOU)

OU with a different optimum on each painted regime of the tree (`OUwie`, `ouch`): each branch
belongs to a regime, and the trait follows an OU pulled toward that regime's optimum. The regimes
come from a **discrete stochastic map** — typically an `Mk` trait simulated on the *same* tree —
so a regime may switch partway along a branch, and the OU is integrated exactly piece by piece.
Parameters: `regimes` (a discrete `TraitResult` carrying the map), `theta` (one optimum per
regime), `alpha` (> 0, scalar or one per regime), `sigma2` (≥ 0, scalar or one per regime), `x0`
(defaults to the optimum of the root's regime). Optionally `alpha`/`sigma2` also vary by regime
(the `OUMV` / `OUMA` / `OUMVA` variants); by default only `theta` differs (`OUM`).

### Multivariate Brownian motion (MultivariateBrownian)

Brownian motion of a **vector-valued** trait with a rate (covariance) matrix `R`: a length-`k`
trait diffuses so the increment over a branch is `MVN(trend·t, R·t)`. The off-diagonal `R[a, b]`
couples dimensions `a` and `b`, so this is the model of **correlated** continuous-trait evolution
(`mvMORPH`, `Rphylopars`): tips are jointly multivariate-normal with covariance `R ⊗ C`. Each
node's value is a length-`k` array. Parameters: `R` (k×k symmetric positive-semidefinite rate
matrix), `x0` (root vector, default zeros), `trend` (per-dimension drift, default zeros).

### Multivariate OU (MultivariateOU)

Multivariate stabilizing selection (`mvMORPH`): `dX = A·(theta − X)·dt + Σ^{1/2}·dW`, with
mean-reversion matrix `A` (`alpha`), optimum vector `theta`, and diffusion covariance `R` (`Σ`).
The exact branch transition has mean `theta + e^{−A·t}·(x − theta)` and covariance
`V − e^{−A·t}·V·e^{−Aᵀ·t}`, where the stationary covariance `V` solves the Lyapunov equation
`A·V + V·Aᵀ = R`. Parameters: `R` (k×k PSD diffusion covariance), `alpha` (mean reversion as a
scalar `alpha·I`, a length-`k` diagonal, or a k×k matrix with eigenvalues of positive real part),
`theta` (optimum vector), `x0` (default `theta`). For a scalar `alpha` this reduces to per-dimension
OU with correlated diffusion.

## Command line

`zombi2 trait` needs a tree — make one first with `zombi2 species` (writing `species_tree.nwk`).
The CLI covers the scalar continuous models `--model bm | ou | eb`; the vector-valued models
(MultivariateBrownian, MultivariateOU) and MultiOptimumOU are Python-only. Shared flags: `--sigma2`
(default `1.0`), `--x0` (root value; OU defaults it to `--theta`), `--trend` (bm/eb). OU adds
`--alpha` (default `1.0`) and `--theta` (default `0.0`); EB adds `--rate` (negative = early burst,
default `1.0`).

```bash
# Brownian motion
zombi2 trait -t run/species_tree.nwk --model bm --sigma2 0.5 --seed 1 -o run/

# Ornstein–Uhlenbeck toward an optimum
zombi2 trait -t run/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o run/

# early burst (rate decays through time)
zombi2 trait -t run/species_tree.nwk --model eb --sigma2 1 --rate -1.5 --seed 1 -o run/
```

Add `--replicates N` to write `traits.tsv` with one column per replicate.

## Python

Models live in `zombi2.traits` (and re-export at the top level, so `zombi2.BrownianMotion` also
works):

```python
import numpy as np
from zombi2.traits import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst,
    MultivariateBrownian, MultivariateOU, MultiOptimumOU, Mk,
    simulate_traits, replicate_traits,
)
from zombi2.species import BirthDeath, simulate_species_tree

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

# scalar continuous models
res = simulate_traits(tree, BrownianMotion(sigma2=0.5, x0=0.0, trend=0.0), seed=1)
res.values                 # {extant leaf: value} — the observable tip data
res.ancestral_states()     # {internal node: value} — exact, not inferred

simulate_traits(tree, OrnsteinUhlenbeck(sigma2=0.4, alpha=2.0, theta=10.0), seed=1)
simulate_traits(tree, EarlyBurst(sigma2=1.0, rate=-1.5), seed=1)
replicate_traits(tree, BrownianMotion(0.5), 100, seed=7)   # 100 independent draws

# correlated vector-valued traits (each node value is a length-k array)
R = [[1.0, 0.6], [0.6, 0.8]]
simulate_traits(tree, MultivariateBrownian(R), seed=1)
simulate_traits(tree, MultivariateOU(R, alpha=2.0, theta=[0.0, 0.0]), seed=1)

# per-regime optima: paint regimes with an Mk map on the SAME tree, then run OU on it
regimes = simulate_traits(tree, Mk.equal_rates(2, 0.4), seed=1)
simulate_traits(tree, MultiOptimumOU(regimes, theta=[-5.0, 5.0], alpha=4.0, sigma2=0.4), seed=2)
```

## Output

A `TraitResult`: `values` at the extant tips (the observable comparative data),
`ancestral_states()` at every internal node (exact, not inferred), and `node_values` for all
nodes. Being continuous, it carries no per-branch history (`history is None`). `to_tsv()` gives a
`node`/`trait` table (pass `nodes="all"` for ancestral rows too) and `to_newick()` a Newick with
`[&trait=…]` on every node. The CLI writes `traits.tsv` (tip **and** ancestral values),
`trait_tree.nwk` (values annotated on every node), and the `trait.log` manifest; `--replicates N`
writes one `traits.tsv` column per replicate instead.

## Validation

- **BM.** Over many replicates the empirical tip mean matches `x0 + trend·depth` and the tip
  covariance matches `sigma2·C` (`C` = the shared-path-length / MRCA-time matrix) —
  `test_traits.py::test_bm_tip_moments_match_theory`.
- **OU.** Over many single-branch replicates the endpoint mean reverts as
  `theta + (x0 − theta)·e^{−alpha·t}` and the variance matches
  `sigma2/(2·alpha)·(1 − e^{−2·alpha·t})` —
  `test_traits.py::test_ou_transition_moments_match_theory`.
- **EarlyBurst.** The single-branch tip variance matches the exact integral
  `sigma2·(e^{rate·t} − 1)/rate` for both a decaying and an accelerating rate —
  `test_traits.py::test_early_burst_variance_matches_integral`.
- **MultiOptimumOU.** With strong `alpha`, tips in each painted regime concentrate near that
  regime's own optimum — `test_traits.py::test_multi_optimum_ou_tracks_local_optima`.
- **MultivariateBrownian.** The per-tip covariance matches `R·depth` and the cross-tip same-dimension
  covariance matches `R[a,a]·MRCA-time` — `test_traits.py::test_mvbm_tip_covariance_matches_R_times_C`.
- **MultivariateOU.** Over many single-branch replicates the endpoint mean matches
  `theta + e^{−A·t}·(x0 − theta)` and the covariance matches `V − e^{−A·t}·V·e^{−Aᵀ·t}` —
  `test_traits.py::test_mvou_single_branch_moments`.

## References

- Felsenstein, J. (1985). Phylogenies and the comparative method. *The American Naturalist*
  125(1): 1–15.
- Hansen, T. F. (1997). Stabilizing selection and the comparative analysis of adaptation.
  *Evolution* 51(5): 1341–1351.
- Butler, M. A. & King, A. A. (2004). Phylogenetic comparative analysis: a modeling approach for
  adaptive evolution. *The American Naturalist* 164(6): 683–695.
- Blomberg, S. P., Garland, T. & Ives, A. R. (2003). Testing for phylogenetic signal in
  comparative data. *Evolution* 57(4): 717–745.
- Clavel, J., Escarguel, G. & Merceron, G. (2015). mvMORPH: an R package for fitting multivariate
  evolutionary models to morphometric data. *Methods in Ecology and Evolution* 6(11): 1311–1319.
