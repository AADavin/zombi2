# Trait API — design target

**Status: the design to build.** The target for rewriting `zombi2/traits`, the detailed consequence of
`SPEC.md` for the trait level. Designed with Adrián on 2026-07-18. **Not built yet** — today's code ships
a thirteen-class model zoo (§ *What to delete*). Parallels `species-api.md`, `genome-api.md`, and
`sequence-api.md`; read those first.

---

## The trait level is a different kind of object (and that is fine)

The other three levels are **genealogies**: lineages, gene copies, and sites are *things born and lost*
along a tree, counted as events. A trait is not born or lost — it is a **value that rides the tree** (a
body size, a habitat, a presence/absence), and you observe the value itself, not an event count. So the
trait level has no "rate of events" the way the others do. That is a real seam, named rather than papered
over. What makes it still belong: the *ways* a value evolves reuse the exact modifier vocabulary the other
levels already have.

## The problem it fixes

`zombi2/traits` ships **thirteen** model classes — `BrownianMotion`, `OrnsteinUhlenbeck`, `EarlyBurst`,
`Mk`, `HiddenStateMk`, `CorrelatedBinary`, `CorrelatedBinaryK`, `MultivariateBrownian`, `MultivariateOU`,
`MultiOptimumOU`, `ThresholdModel`, `DEC`, `Cladogenesis` — a zoo exactly like the seven species processes
and the eight-class clock hierarchy. Almost all of it is the same handful of ideas wearing literature
names.

## Two entry points, split by state space

The one thing that genuinely differs is the state space — a real value vs a discrete state — and it gives
genuinely different argument sets. So the *kind is the function*, not an argument (same principle as the
genome level's three functions):

```python
traits.simulate_continuous(tree, …)   # real-valued: BM / OU / EB
traits.simulate_discrete(tree, …)     # finite states: Mk / threshold
```

A trait rides any tree, so `tree` is a species tree (or any tree). One value per branch, full stop — no
per-family multiplicity, which is why traits are *cleaner* than sequences (no "which tree does it ride"
question).

## Continuous traits: Brownian motion is the native process

A continuous trait does Brownian motion natively — that is what it *is*, not a modifier bolted on. The
variants are the same knobs the rate-modifiers expose one level over:

```python
# BM — a body size diffusing at variance-rate σ² = 1.0
traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)

# OU — same diffusion, pulled toward an optimum
#   reverts_to + pull are the SAME knobs that turn the autocorrelated clock into CIR
traits.simulate_continuous(tree, start=0.0, rate=1.0, reverts_to=2.0, pull=0.5, seed=1)

# EB (early burst) — the diffusion rate decays through time
#   ...the SAME Time modifier that gives species its skyline
traits.simulate_continuous(tree, start=0.0, rate=1.0 * mod.Time({0: 1.0, 5: 0.2}), seed=1)
```

The unification is at the level of *knobs*, not a shared wrapper: `reverts_to`/`pull` (→ CIR clock) and
`Time` (→ species skyline) are literally the same knobs. `rate` is the BM variance-rate σ², and it takes
modifiers like any other rate.

## Discrete traits: a state switching along the tree (Mk)

```python
# Mk — habitat flips between two states at rate 0.1
#   ...the trait twin of the discrete clock, mod.Markov(states=…, switch=…)
traits.simulate_discrete(tree, states=["marine", "terrestrial"], switch=0.1,
                         start="marine", seed=1)

# asymmetric — a rate matrix instead of one symmetric rate
traits.simulate_discrete(tree, states=["absent", "present"],
                         switch={"absent->present": 0.2, "present->absent": 0.05}, seed=1)

# threshold — a discrete state driven by an underlying continuous liability (BM)
traits.simulate_discrete(tree, states=["absent", "present"],
                         liability=1.0, threshold=0.0, seed=1)
```

`simulate_discrete(states=…, switch=…)` is the literal twin of `mod.Markov(rates=…, switch=…)`: a set of
states plus a switching rate.

## Correlated traits: the joint rule, inside a level

Independent traits are separate calls. Two traits that drift *together* cannot be separated — that is
`P(size, limb)`, **joint** — so they go in **one call**. It is the Ch2 joint rule ("joint → one command"),
applied within the trait level rather than across levels; nothing new to teach.

**How correlation is specified (decided: per-trait rates + a correlation overlay, not a full Σ matrix).**
Keeping the rates per-trait preserves the whole grammar (each trait keeps its own modifiers); correlation
is a clean, dimensionless, bounded overlay (∈ [−1, 1]), easier to validate than a positive-semi-definite
covariance matrix, and reads the way people think:

```python
traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8 * mod.Time({0: 1, 5: 0.3})},   # each keeps its own modifiers
    correlation={("size", "limb"): 0.6},                           # overlay, ∈ [−1, 1]
    seed=1)
```

The **same `correlation=` overlay handles discrete correlation for free**, via the threshold model —
correlated liabilities underneath, thresholds on top (this is the Felsenstein threshold route, and it
replaces the `CorrelatedBinary` CTMC machinery):

```python
traits.simulate_discrete(tree, states=["absent", "present"],
    liability={"wings": 1.0, "flight": 1.0},
    correlation={("wings", "flight"): 0.7}, threshold=0.0, seed=1)
```

One overlay, working on continuous traits directly and on discrete traits through their liabilities. A
full Σ matrix may still be *accepted* as an alternative input for the comparative-methods crowd, but the
per-trait + `correlation=` form is the surface.

## The literature → command bridge (goes in the chapter)

| Literature name | What it does | ZOMBI2 |
|---|---|---|
| Brownian motion (BM) | a value diffusing | `simulate_continuous(rate=…)` |
| Ornstein–Uhlenbeck (OU) | diffusion pulled to an optimum | `simulate_continuous(rate=…, reverts_to=…, pull=…)` |
| Early burst (EB / ACDC) | diffusion rate decays through time | `simulate_continuous(rate=1.0 * mod.Time({…}))` |
| Multivariate BM / OU | traits evolving together | one `simulate_continuous(rate={…}, correlation={…})` call |
| Mk (k-state Markov) | a discrete state switching | `simulate_discrete(states=…, switch=…)` |
| Threshold / liability (Wright–Felsenstein) | discrete driven by continuous liability | `simulate_discrete(liability=…, threshold=…)` |
| Correlated binary / Pagel | discrete traits evolving together | `simulate_discrete(liability={…}, correlation={…})` |
| DEC biogeography | range = set of areas | → **experimental** (purged for now) |
| BiSSE / MuSSE / QuaSSE | trait drives speciation | **not a trait model** — trait↔species *joint*, Part III |

## Still to design

- **Decided: unify.** `Inherited(spread=, reverts_to=, pull=)` everywhere — plain `spread` = pure drift
  (BM / ClaDS / autocorrelated clock); add `reverts_to` (target) + `pull` (strength) = mean-reverting (OU
  trait / CIR clock). The CIR clock grows a `pull`; OU and CIR share the same two knobs.
- **`MultiOptimumOU`** — the optimum shifts on certain branches (regime painting). An advanced case;
  probably a `regimes=` argument on `simulate_continuous`. Deferred, named honestly.
- **`Cladogenesis` — spelling only; placement is decided.** A trait that jumps *at speciation nodes*
  rather than along branches. **Within-level, not Part III (SPEC §4):** a jump at speciation *reads* the
  tree it already lives on — it does not change *which* tree exists — so it is an option of the trait's own
  model (an `at_speciation=` jump rule on `simulate_discrete`/`simulate_continuous`), not a coupling. Only
  the argument name is open. It becomes *joint* only when the same trait *also drives* speciation (SSE with
  cladogenetic change), which is Part III.
- **`HiddenStateMk`** — hidden rate categories under an Mk trait (the trait twin of the `Markov` clock's
  hidden classes). Likely a hidden-state option on `simulate_discrete`; deferred.
- **Decided:** `switch=` accepts **both** a string-keyed `{"a->b": rate}` dict (readable, few states) and
  a numeric matrix with `states=[...]` (many states); `switch=0.1` stays the symmetric shortcut.
  (`pull`/`reverts_to` names stand.)

## What to delete / change in `zombi2/traits`

- Delete the model-class zoo (`BrownianMotion`, `OrnsteinUhlenbeck`, `EarlyBurst`, `Mk`, `HiddenStateMk`,
  `CorrelatedBinary`, `CorrelatedBinaryK`, `MultivariateBrownian`, `MultivariateOU`, `MultiOptimumOU`,
  `ThresholdModel`, `Cladogenesis`). They become two functions + shared knobs + the `correlation=` overlay.
- Move `DEC` / `biogeography.py` to `zombi2.experimental` — purged from the trait level for now,
  recoverable later.
- SSE (BiSSE/MuSSE/QuaSSE) leaves the trait level entirely — it is trait↔species *joint*, Part III.
- Rates follow the cross-level `count(base) × modifiers` grammar; the trait's `rate` is the BM variance-rate.
- Keep the `TraitResult` output object.
