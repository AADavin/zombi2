# Sequence API — design target

**Status: the design to build.** The target for rewriting `zombi2/sequences`, the detailed consequence
of `SPEC.md` for the sequence level. Designed with Adrián on 2026-07-18. **Not built yet** — today's code
ships an eight-class `Clock` hierarchy (§ *What to delete*). Parallels `species-api.md` and
`genome-api.md`; read those first.

---

## The problem it fixes

`zombi2/sequences/clocks.py` ships **eight** clock classes — `Clock`, `StrictClock`,
`UncorrelatedLogNormalClock`, `UncorrelatedGammaClock`, `WhiteNoiseClock`,
`AutocorrelatedLogNormalClock`, `CIRClock`, `RateVariation` — a zoo exactly like the seven species
processes and the `RateModel` hierarchy. But a clock is not a special kind of object: it is a **modifier
on the substitution rate**, the same `count(base) × modifiers` grammar as every other level (`SPEC §5`).
"That is the one place the word *clock* belongs" (Ch2). Once clocks are modifiers, the zoo collapses to
the same small shared vocabulary the rest of ZOMBI2 uses.

## The entry point

One function, taking the gene trees a genome run produced:

```python
sequences.simulate_sequences(gene_trees, model=hky85(kappa=2.0), length=1000, seed=1)
```

(Today the surface is a `SequenceEvolution` class plus `evolve_on_tree`; the function is the target.)

## Two things vary, and they are different axes

A substitution rate can vary **across branches** (the clock) and **across sites** (rate heterogeneity).
They are orthogonal and must not be conflated:

- **Across branches — the clock.** A modifier on the substitution rate, exactly like a modifier at any
  other level. This is where the zoo lived.
- **Across sites — +Γ.** Some sites evolve faster than others, drawn from a Gamma. This is the classic
  `+Γ` of phylogenetics and stays its own argument (`gamma=0.5`, the shape α), not a clock.

```python
sequences.simulate_sequences(gene_trees,
    model=gtr(...),                       # the substitution model (a menu, see below)
    substitution=1.0 * mod.ByBranch(spread=0.3),   # the clock: across-branch variation
    gamma=0.5,                            # +Γ: across-site variation (shape α)
    length=1000, seed=1)
```

## The substitution model is a MENU, not a zoo

Not everything collapses, and that is correct. JC69 / HKY / GTR / LG are genuinely different chemistry —
different rate matrices — so they stay a **menu of constructors**, each taking its own physical
parameters:

```python
model=jc69()                       # no free parameters
model=hky85(kappa=2.0)             # transition/transversion bias
model=gtr(rates=..., freqs=...)    # six exchangeabilities + base frequencies
model=lg()                         # empirical amino-acid matrix
model=gy94(omega=0.2)              # codon model, dN/dS
```

Faking a grammar over the matrices would be worse than a menu. The menu is the honest shape here.

## The clock collapses to three modifiers + the strict default

Every clock in the literature is one of three things happening to the rate along the tree, expressed with
the **same modifiers the other levels use**:

```python
# strict clock — no across-branch variation (the default; write nothing)
substitution = 1.0

# uncorrelated / relaxed — each branch draws its rate independently (i.i.d.)
substitution = 1.0 * mod.ByBranch(spread=0.3)                 # lognormal (default)
substitution = 1.0 * mod.ByBranch(spread=0.3, dist="gamma")   # gamma; white-noise is another dist

# autocorrelated — the rate drifts continuously along the tree (geometric Brownian)
substitution = 1.0 * mod.Inherited(spread=0.3)

# CIR — the same, but mean-reverting (Ornstein–Uhlenbeck)
substitution = 1.0 * mod.Inherited(spread=0.3, reverts_to=1.0)

# the Markov clock — the rate hops between discrete categories along the branches
substitution = 1.0 * mod.Markov(rates=[0.5, 1.0, 2.0], switch=0.1)
```

Three modifiers, grouped by what memory the rate has:

- **`ByBranch`** — *no memory*: each branch independent. (The uncorrelated / relaxed family.)
- **`Inherited`** — *continuous memory*: the rate drifts, parent to child. (Autocorrelated; CIR is this
  with `reverts_to`.)
- **`Markov`** — *discrete memory*: the rate switches between a few states via a CTMC on rate categories.

Two of the three are **shared across levels**, which is the whole point of the grammar:

- `ByBranch` is the branch-twin of the genome level's `ByFamily` — the same i.i.d.-heterogeneity idea.
- `Inherited` is **literally the species `Inherited`** (ClaDS): a rate that drifts along the tree. The
  autocorrelated molecular clock and ClaDS diversification are the same modifier at two levels.
- `Markov` is new to sequences, but even it echoes species: a clade shift is one discrete rate jump;
  `Markov` is that, happening repeatedly at a rate.

## The literature → command bridge (goes in the chapter)

The deprecated model names survive **only** in this table — a reader who knows "I want a CIR clock" finds
the command, and no one has to memorise the acronyms anywhere else. Every chapter carries one of these;
this is the sequence chapter's.

| Literature name | What it does | ZOMBI2 |
|---|---|---|
| Strict / global clock | one rate everywhere | `substitution = 1.0` (default) |
| Uncorrelated lognormal (UCLN) | each branch i.i.d. lognormal | `1.0 * mod.ByBranch(spread=…)` |
| Uncorrelated gamma (UGAM) | each branch i.i.d. gamma | `1.0 * mod.ByBranch(spread=…, dist="gamma")` |
| White-noise clock | each branch i.i.d., short memory | `1.0 * mod.ByBranch(spread=…, dist=…)` |
| Autocorrelated lognormal (Thorne–Kishino) | rate drifts along the tree | `1.0 * mod.Inherited(spread=…)` |
| CIR clock | drift with mean-reversion | `1.0 * mod.Inherited(spread=…, reverts_to=…)` |
| Discrete-category / random local clock | rate hops between categories | `1.0 * mod.Markov(rates=[…], switch=…)` |
| +Γ rate heterogeneity | variation across sites | `gamma=α` (not a clock) |

## Still to design

- Whether `ByBranch` exposes the distribution as `dist=` (lognormal / gamma / white-noise) or as separate
  helpers; the default is lognormal.
- Whether CIR's mean-reversion is `reverts_to=` on `Inherited` or its own modifier (leaning: a parameter
  on `Inherited`, since it is the same continuous-drift process with a pull term).
- `Markov`'s exact signature (`rates=` as multipliers vs absolute; `switch=` one global rate vs a matrix).
- Where a **trait- or driver-conditioned** clock lives — that is `DriverClock` today, and it belongs to
  Part III (a `traits → sequences` conditioning), not to this menu.
- Codon-model surface (`gy94`/`mg94`, the M-series site models) — already shipped; confirm they read as
  menu constructors here.

## What to delete / change in `zombi2/sequences`

- Delete the `Clock` hierarchy (`StrictClock`, `UncorrelatedLogNormalClock`, `UncorrelatedGammaClock`,
  `WhiteNoiseClock`, `AutocorrelatedLogNormalClock`, `CIRClock`, `RateVariation`). Clocks become the
  shared `mod.ByBranch` / `mod.Inherited` / `mod.Markov` modifiers on the `substitution` rate.
- Keep the substitution models as a **menu** of constructors (`jc69`, `k80`, `hky85`, `gtr`, `lg`, …);
  they are genuinely different matrices, not a zoo.
- Keep `+Γ` as its own `gamma=` argument (across-site variation ≠ across-branch clock).
- Add the `sequences.simulate_sequences(gene_trees, …)` entry point over the existing evolution core.
- `DriverClock` → a Part III conditioning (`traits`/level → sequences), not a menu member.
