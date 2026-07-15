# Design: the opportunity axis ("how many clocks?")

**Status:** proposed (2026-07-15). A clarity upgrade: make the *opportunity* of a rate — the "per
what?" — a first-class idea in the docs and, eventually, in the code. Extends the two-term rate rule
([naming-consolidation](naming-consolidation.md) C3): a **rate** is `time⁻¹`, a **modifier** is
dimensionless — and now the second axis, **opportunity**, gets the same treatment.

## The idea, in one metaphor

Think of every pending event as a **clock** that ticks at some rate; a tick fires the event. A whole
model is two answers:

1. **How many clocks?** — one per lineage / one per gene copy / one *shared* by the whole tree. This
   is the **opportunity**.
2. **How fast does each tick?** — the base rate, possibly varying with time, crowding, or lineage.

`propensity = base × opportunities × modifiers`. The base carries the units (`time⁻¹`); opportunities
(a count) and modifiers (a multiplier) are dimensionless. **"How many clocks" is the opportunity, and
it decides exponential vs linear** — because a per-lineage count *tracks the growing quantity* and a
shared count does not.

## Where ZOMBI2 is today

- **Genomes** already expose the opportunity: `PerCopyRates` (one clock per gene copy → exponential
  families) vs `PerLineageRates` (one per genome → linear). Good.
- **Species** models *all* answer "one clock per lineage" and only vary the *speed*:

  | Model | clocks | speed | tree |
  |---|---|---|---|
  | `BirthDeath` / `Yule` | per lineage | fixed | exponential |
  | `DiversityDependent` | per lineage | slows as tree fills (`λ₀(1−n/K)`) | saturating |
  | `ClaDS` | per lineage | each its own, inherited + jittered | uneven |
  | `EpisodicBirthDeath` | per lineage | fixed within an epoch | piecewise |
  | `CladeShiftBirthDeath` | per lineage | a clade's jump at a set time | clade burst |
  | **(missing)** | **shared** | fixed | **linear** |

  The empty row is the gap: there is no "one shared clock" species model, even though it's a natural,
  different process (a fixed diversification *budget* the lineages compete for).

## The plan

**Part 1 — docs (byte-identical).** Make "how many clocks / how fast" the spine:
- `rates.md` — already states the two axes and the per-lineage↔shared birth–death example; add the
  clock metaphor as the lead and the "place every model" table.
- `species-trees.md` — reframe the diversification models from *a list* into *one axis* (clocks ×
  speed), so the reader sees the missing shared row for themselves.

**Part 2 — code (additive, forward-only).** Add the missing model so the example is runnable:
- `SharedBirthDeath(birth, death)` — one shared clock. In the existing Gillespie view pattern it is
  simply `lineage_rates(state, n) = (birth/n, death/n)`, so the totals are constant (`birth`, `death`)
  regardless of `n` → linear diversification (equivalently, per-lineage `λ(n)=Λ/n` — diversity-
  dependence in disguise). Registered like `DiversityDependent`: `SpeciesCaps(GILLESPIE,
  supports_n_tips=True)`, a `_SharedView`, `--diversification shared`. Forward-only (no closed-form
  backward sampler, exactly like `DiversityDependent`/`ClaDS`). Touches nothing existing → the current
  suite stays byte-identical.

**Part 3 — code, later (the elegant endpoint).** Today "how many clocks" is *implicit* (whether a
view divides by `n`). Make **opportunity a named, per-event knob** (`per-lineage | shared`, per birth
/ death / D / T / L …), so the code reads the way the docs teach and the whole matrix — including
*shared birth + per-lineage death* (a self-limiting tree: births capped, deaths grow with N) — is
expressible on purpose, not by a `/n` trick. A real interface refactor; deferred until Parts 1–2 land.

## Naming (one decision)

Proposed: **`SharedBirthDeath`** / **`--diversification shared`** — names the *mechanism* (a shared
diversification budget), parallels `PerLineageRates` (mechanism, not outcome), and matches the
"shared clock / shared budget" language in the docs. Alternatives: `global` (overloaded elsewhere),
`linear` (names the *outcome*, not the cause). Flagged for Adrián to veto before it ships.

## The rule, in one line

> **A rate is `time⁻¹`; a modifier is dimensionless; the opportunity is *how many clocks* — per
> lineage, per copy, or shared — and that count, not the units, is what makes diversification
> exponential or linear.**
