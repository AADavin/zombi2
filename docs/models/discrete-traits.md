# Discrete-trait models

A **discrete trait** is a character that takes one of finitely many states — a habitat, the
presence or absence of a structure, a diet class. ZOMBI2 evolves it as a **continuous-time Markov
chain** along the tree: given a rate matrix `Q` (off-diagonal `Q[i,j] ≥ 0` is the instantaneous
rate `i → j`), every branch's jumps are simulated **exactly** by the Gillespie algorithm, so the
realized `(state, duration)` history — a *stochastic character map* — comes for free (`result.history`,
`result.changes()`). All the models below are one shared `Mk` engine wearing different `Q` matrices:
they differ only in what transition structure they impose. One function, `simulate_traits`, runs
them all.

| Model | Structure | Reach for it when |
| --- | --- | --- |
| **Mk** | one `k×k` rate matrix `Q` | any single discrete character (ER / SYM / ordered / ARD) |
| **CorrelatedBinary** | two binary traits, one flips at a time | testing whether two binary characters coevolve (Pagel) |
| **CorrelatedBinaryK** | `k` binary traits, one flips at a time | the `k`-trait generalization of the above |
| **HiddenStateMk** | observed states × hidden rate classes | rate heterogeneity a plain `Mk` cannot capture (corHMM) |
| **ThresholdModel** | discrete state cut from a latent BM liability | a discrete character with an underlying continuous cause |

## The models

### Mk

A `k`-state chain with an arbitrary rate matrix `Q` (Lewis 2001). Convenience constructors cover the
standard sub-models: `Mk.equal_rates(k, rate)` is all-to-all at one shared rate (**ER**),
`Mk.symmetric(R)` makes `i→j` and `j→i` equal (**SYM**), `Mk.ordered(k, rate)` is the tridiagonal
nearest-neighbour chain (adjacent steps only, `i ↔ i±1` — a meristic character), and the raw
constructor `Mk(Q)` takes any matrix, including all-rates-different (**ARD**). The root state is
`"uniform"` by default (`"stationary"`, an index, or a probability vector also work), and every `Mk`
exposes its analytic `transition_matrix(t)` = `exp(Q·t)` and `stationary_distribution()`. This is the
default discrete model and the base class of every other model on this page.

### CorrelatedBinary

Two binary characters **X** and **Y** evolving jointly over the four states `(X, Y)`, with **one
trait changing at a time** (simultaneous double flips have rate 0). Each trait's gain/loss rate may
depend on the *other* trait's current state — that dependence *is* correlated evolution (Pagel 1994).
Pass the eight directional rates (`x_gain_y0`, `x_gain_y1`, … named by the changing trait, its
direction, and the other trait's state); `CorrelatedBinary.independent(x_gain, x_loss, y_gain,
y_loss)` builds the null model where the two evolve independently, against which the dependent fit is
tested. Each node's value is the `(X, Y)` pair. Reach for it to ask whether two binary characters
coevolve.

### CorrelatedBinaryK

The `k`-trait sibling of `CorrelatedBinary` (Pagel 1994, generalized to `k ≥ 2`). `k` binary traits
evolve over the `2^k` configurations, still flipping **exactly one bit at a time**, and each flip may
depend on the other traits' states. The compact entry points are `CorrelatedBinaryK.independent(gains,
losses)` (the null model — `Q` is the Kronecker sum of `k` independent 2-state chains),
`.equal_rates(k, gain, loss)` (one shared gain/loss), and `.partner_coupling(gains, losses, partners,
boost_gain, boost_loss)` (each trait multiplicatively boosted by one designated partner — `O(k)`
parameters that still induce genuine pairwise correlation). `.from_table` / the raw `rate_fn`
constructor give full generality. Each node's value is the `k`-tuple. Reach for it when you have more
than two coevolving binary characters.

### HiddenStateMk

An observed character whose transition rates depend on an unobserved **hidden rate class** that itself
switches along the tree (corHMM; Beaulieu et al. 2013). Pass one `O×O` observed-rate matrix **per
hidden class** (`observed_rates`, e.g. a slow class and a fast class) and the `hidden_rate` between
classes (a scalar for a symmetric all-to-all rate, or an `H×H` matrix). The full state is the
`(observed, hidden)` pair; `result.labeled_values()` reports only the observed part, while
`result.full_label(v)` and `result.changes()` expose the hidden dimension. Reach for it when a plain
`Mk` cannot absorb the rate heterogeneity in a character.

### ThresholdModel

A discrete state derived from a latent continuous **liability** that evolves by Brownian motion; the
observed state is whichever interval the liability currently falls in, cut by an ordered set of
`thresholds` (`k−1` cuts give `k` states; `[0.0]` is a binary trait) (Felsenstein 2012). Only the
ratio of thresholds to the diffusion scale is identifiable, so `sigma2` is fixed to `1.0` by default;
`x0` sets the root liability. The evolving value at each node is the liability (`result.values`,
continuous); the observed discrete state comes from `result.labeled_values()`. Reach for it to model a
discrete character with an underlying continuous cause (and the correlated / polymorphic behaviour that
implies). Its `kind` is `"continuous"` because the liability, not the state, is what diffuses.

## Command line

The `trait` command covers the two discrete models that need no auxiliary structure — **Mk** and
**ThresholdModel**. (`CorrelatedBinary`, `CorrelatedBinaryK`, and `HiddenStateMk` are Python-only, as
their multi-trait / hidden-class parameterizations have no flat CLI form.) `--model mk` is
equal-rates by default; `--ordered` gives the adjacent-only chain and `--q-matrix FILE` reads an
arbitrary `Q` (overriding `--states`/`--rate`/`--ordered`).

```bash
# a tree to evolve the trait along
zombi2 species --birth 1 --death 0.3 --tips 30 --age 5 --seed 1 -o run/

# Mk: 3-state equal-rates character
zombi2 trait -t run/species_tree.nwk --model mk --states 3 --rate 0.4 --seed 1 -o mk/

# Mk: ordered (meristic) 4-state character, 20 replicates
zombi2 trait -t run/species_tree.nwk --model mk --states 4 --ordered --replicates 20 --seed 1 -o mko/

# threshold: binary (default) and a 3-state ordered character
zombi2 trait -t run/species_tree.nwk --model threshold --seed 1 -o th/
zombi2 trait -t run/species_tree.nwk --model threshold --thresholds=-1,1 --seed 1 -o th3/
```

(Invoke the CLI as `python -m zombi2 trait …`.) `--replicates N` writes `traits.tsv` with one column
per replicate.

## Python

Models live in `zombi2.traits` (and re-export at the top level, so `zombi2.Mk` also works):

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import (Mk, CorrelatedBinary, CorrelatedBinaryK,
                           HiddenStateMk, ThresholdModel, simulate_traits)

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=5.0, seed=1)

# Mk: equal-rates over labeled states; the stochastic map comes for free
mk = simulate_traits(tree, Mk.equal_rates(3, 0.4, states=["marine", "brackish", "fresh"]), seed=2)
mk.labeled_values()                    # {extant leaf: "marine" | "brackish" | "fresh"}
mk.history[node]                       # [(state, duration), ...] — the stochastic map
Mk.symmetric([[0, 2, 1], [2, 0, 3], [1, 3, 0]])     # SYM
Mk.ordered(4, 0.5)                                   # ordered / meristic
Mk([[0, 1, 2], [3, 0, 1], [1, 1, 0]])               # ARD (any Q)

# CorrelatedBinary: Y tracks X (gained fast when X=1, lost fast when X=0)
cb = CorrelatedBinary(x_gain_y0=0.5, x_gain_y1=0.5, x_loss_y0=0.5, x_loss_y1=0.5,
                      y_gain_x0=0.05, y_gain_x1=2.0, y_loss_x0=2.0, y_loss_x1=0.05)
simulate_traits(tree, cb, seed=1).labeled_values()          # {leaf: (X, Y)}
CorrelatedBinary.independent(x_gain=0.5, x_loss=0.5, y_gain=0.5, y_loss=0.5)   # null model

# CorrelatedBinaryK: three binary traits (independent null, or partner-coupled)
CorrelatedBinaryK.independent(gains=[0.4, 0.4, 0.4], losses=[0.6, 0.6, 0.6])
cbk = CorrelatedBinaryK.partner_coupling(gains=0.3, losses=0.3,
                                         partners=[1, None, None], boost_gain=5.0)  # trait 0 tracks trait 1
simulate_traits(tree, cbk, seed=1).labeled_values()         # {leaf: (t0, t1, t2)}

# HiddenStateMk: a slow and a fast hidden class over an observed binary character
hmm = HiddenStateMk(observed_rates=[[[0, 0.1], [0.1, 0]], [[0, 3.0], [3.0, 0]]],
                    hidden_rate=0.5, observed_states=[0, 1], hidden_states=["slow", "fast"])
res = simulate_traits(tree, hmm, seed=1)
res.labeled_values()                                        # observed 0/1 (hidden collapsed)
res.full_label(res.node_values[tree.extant_leaves()[0]])    # (observed, hidden), e.g. (1, 'fast')

# ThresholdModel: binary state from a latent BM liability
th = simulate_traits(tree, ThresholdModel(thresholds=[0.0]), seed=1)
th.values                              # liabilities (continuous, latent)
th.labeled_values()                    # observed 0/1 states
```

## Output

`simulate_traits` returns a `TraitResult`: `values` are the extant-tip states (the observable data),
`labeled_values()` decodes state indices to their labels (and returns the `(X, Y)` / `k`-tuple pair
for the correlated models), `ancestral_states()` gives the exact internal-node states, and
`history` / `changes()` give the per-branch stochastic map and the realized transition events.
`to_tsv()` and `to_newick()` write the tip table and a `[&trait=…]`-annotated Newick.

From the CLI, `trait` writes `traits.tsv` (`node <TAB> trait`, or one `rep_*` column per replicate),
`trait_tree.nwk` (the Newick with `[&trait=…]` on every node), and `trait.log`. For `ThresholdModel`,
the reported trait is the discrete **state**, not the latent liability.

## Validation

- **Mk.** `P(t)` for the equal-rates model matches the closed form
  `P_ii = 1/k + (1−1/k)e^{−kqt}`, `P_ij = (1/k)(1−e^{−kqt})`
  (`test_traits.py::test_mk_equal_rates_transition_closed_form`), and the empirical end-state
  distribution over a single branch matches `P(t)` from the start state across 20 000 replicates
  (`test_traits.py::test_mk_simulation_matches_transition_matrix`).
- **CorrelatedBinary.** The assembled `Q` zeroes both double-transitions (no simultaneous flips),
  sums to zero per row, and places each named rate in the right off-diagonal cell
  (`test_traits.py::test_correlated_binary_Q_structure`).
- **CorrelatedBinaryK.** The `independent` constructor produces exactly the Kronecker sum of the `k`
  per-trait 2-state chains, for `k = 2..5`
  (`test_correlated_binary_k.py::test_independent_is_kronecker_sum`).
- **HiddenStateMk.** With the same observed-rate matrix in every hidden class, the hidden dimension is
  irrelevant and the observed binary character stays symmetric (tip mean ≈ 0.5)
  (`test_traits.py::test_hidden_state_mk_same_rates_are_observed_symmetric`).
- **ThresholdModel.** A symmetric binary threshold (`thresholds=[0.0]`, `x0=0`) yields balanced tip
  states (mean ≈ 0.5) (`test_traits.py::test_threshold_binary_symmetric_is_balanced`).

## References

- Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from discrete morphological
  character data. *Systematic Biology* 50(6): 913–925.
- Pagel, M. (1994). Detecting correlated evolution on phylogenies: a general method for the
  comparative analysis of discrete characters. *Proceedings of the Royal Society B* 255: 37–45.
- Beaulieu, J. M., O'Meara, B. C. & Donoghue, M. J. (2013). Identifying hidden rate changes in the
  evolution of a binary morphological character. *Systematic Biology* 62(5): 725–737.
- Felsenstein, J. (2012). A comparative method for both discrete and continuous characters using the
  threshold model. *The American Naturalist* 179(2): 145–156.
