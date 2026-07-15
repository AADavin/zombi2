# Design: opportunity as a first-class knob (opportunity, Part 3)

**Status:** proposed (2026-07-15). The third and final part of the opportunity line
([opportunity](opportunity.md) — Parts 1 (docs) and 2 (`SharedBirthDeath`) shipped in PR #146). Where
Part 1 *taught* the opportunity axis and Part 2 *added the one missing model*, Part 3 makes
**opportunity a thing you select**, not a thing you pick a class for. The goal, in Adrián's words:

> when you pick your model you just have an option to select the unit the clock is attached to.

Three words drive it: **clarity** (the code says what the docs teach), **elegance** (one knob replaces a
family of classes), **flexibility** (mix the unit per event → models you can't build today).

## The idea in one line

> Every rate is a **clock**. The **opportunity** is the *unit the clock is attached to* — a site, a
> gene copy, a lineage, or the whole thing at once (`shared`). Today that choice is hidden inside a
> class name; Part 3 turns it into a knob you set.

## Governing principles

- **P1 — Opportunity is one axis, orthogonal to the rate.** A model is `base × opportunity × modifiers`
  (the rule from [rate-vocabulary](rate-vocabulary.md) / [rate-modifiers](rate-modifiers.md)). The base
  carries the units (`time⁻¹`); the opportunity is a **count**; modifiers are multipliers. Part 3 only
  touches the middle term.
- **P2 — One vocabulary, every level.** The same word means the same idea for species, genomes,
  sequences. `lineage` is a lineage whether it is speciating or duplicating; `shared` is one clock for
  the whole process wherever it appears.
- **P3 — Three per-unit rungs nest; `shared` sits apart.** The per-unit opportunities form one ladder,
  finest to coarsest: **`site ⊂ copy ⊂ lineage`**. **`shared`** is the odd one out — *not* per-unit but
  *one clock for the whole thing*. (We deliberately do **not** call the top rung `clade`: that word
  implies you must define or select a group, and you never do — `shared` just means "one, for
  everything.") A level advertises which rungs it allows — **scoped by representation** — and nonsense
  (per-copy speciation, per-site duplication) is rejected, not silently ignored.
- **P4 — Zero behaviour change; classes become presets.** `SharedBirthDeath`, `PerCopyRates`,
  `PerLineageRates` do not vanish — they become named shorthands for a `per=…` setting, via the
  deprecation machinery already built ([naming-consolidation](naming-consolidation.md)). Same seed →
  byte-identical output.

## The unified ladder

The opportunity is *what the clock's count tracks*. That single choice sets exponential-vs-linear,
because growth is exponential exactly when the count **tracks the growing quantity**.

| opportunity | the count tracks… | rides events like | lives at |
|---|---|---|---|
| **`site`** (per nt) | number of sites / nucleotides | substitution, insertion, deletion | sequences, nucleotide genome |
| **`copy`** (per gene) | number of gene copies | duplication, loss, transfer, conversion | genome content |
| **`lineage`** (per genome) | number of lineages | speciation, extinction; genome-level D/T/L | species, genome content |
| **`shared`** (one for all) | nothing — fixed at 1 | shared birth–death; a family's shared budget | species, genome content |

Reading the ladder as growth laws for one focal quantity (say, a gene family's copy number):
`copy` → **exponential** (count tracks copies), `lineage` → **linear** (one clock per genome,
regardless of copies), `shared` → **linear and pooled** (one clock for the whole family).
The finest rung, `site`, is the substitution/indel world Adrián flagged: a per-nucleotide clock, the
same axis one level down.

## Opportunity is not heterogeneity

Opportunity is easy to confuse with *per-family rates*, but they are different slots of
`base × opportunity × modifiers`:

- **opportunity** (`per=`) — *how the rate scales with counts* (per site / copy / lineage / shared).
  The middle term.
- **heterogeneity** — *which families (or branches, or pairs) get which values*. It lives in the
  **base** (`FamilySampledRates` — each family its own rate) or a **modifier** (`FamilyModifier` — a
  shared base times a per-family factor), never in the opportunity.

A family can be `per="copy"` **and** carry its own rate; the slots are orthogonal. "Per family" is
therefore never an opportunity value — there is no `per="family"`. (Worked table in
[the rates primer](../guide/rates.md).)

One wrinkle this design should close: `FamilySampledRates` currently hard-codes opportunity = `copy`,
so *per-family rates with a per-lineage opportunity* is reachable only through the modifier route
(`FamilyModifier` on a per-lineage base). Giving the base models a `per=` of their own — part of
phase B — makes all three slots independently selectable, so the orthogonality the docs teach is the
orthogonality the code offers.

## Where ZOMBI2 is today

The knob is **half-born**, and asymmetrically:

- **Genomes** already expose two rungs on the CLI — `--rate-per {copy, lineage}` (i.e.
  `PerCopyRates` vs `PerLineageRates`). But there is no `shared`, and the choice is a *class*, not a
  per-event setting.
- **Species** hide the rung entirely inside class names: `BirthDeath` *is* per-lineage,
  `SharedBirthDeath` *is* shared — but you'd never know the word "opportunity" was involved.
- **Sequences / nucleotide** bake `site` in as the only option (substitution and indel rates are
  per-nucleotide by construction); the axis is there but not named.

So Part 3 is less "invent" than "**finish and unify**": lift the existing `--rate-per` idea to a single
knob, give it the missing `shared` rung, and make it the same concept at every level.

## The design

**The knob.** A per-model default, with an optional per-event override.

```python
# --- species: the rung comes out of the class name ---
z.BirthDeath(1.0, 0.2)                    # per="lineage"  (default → exponential)
z.BirthDeath(1.0, 0.2, per="shared")      # one shared clock → linear   (== SharedBirthDeath)

# --- genomes: one rate model, opportunity as a knob ---
z.Rates(duplication=0.5, loss=0.3)               # per="copy"     (default → exponential families)
z.Rates(duplication=0.5, loss=0.3, per="lineage")# one clock per genome → linear
z.Rates(duplication=0.5, per="shared")           # one budget for the whole family → pooled/linear

# --- the flexibility peak: mix the unit per event ---
z.Rates(duplication=z.Per("shared", 0.5),        # births capped family-wide …
        loss=z.Per("copy", 0.3))                 # … deaths ride every copy  → a self-limiting family
```

```
zombi2 species  --birth 1 --death 0.2 --per shared
zombi2 genomes  --duplication 0.5 --per copy|lineage|shared
```

**Presets stay.** `SharedBirthDeath`, `PerCopyRates`, `PerLineageRates` remain as friendly names that
set `per=` under the hood (deprecated only if we decide to; they read well, so likely kept as
sugar). Beginners never meet the knob; power users turn it.

**Validation is part of the clarity — and it is scoped by representation.** Each level *and genome
representation* publishes the rungs it allows; an illegal combination (per-copy speciation, per-site
duplication) is a clear error, not a silent no-op. Crucially, the `shared` **unit follows the
representation**: in the *unordered* genome a family is the unit of duplication, so `shared` pools
**per family** (the clean `SharedBirthDeath` analogue); in the *ordered / nucleotide* genomes
duplication is *positional* (a segment, possibly spanning families), so `shared` pools **per genome**
(and the finest rung, `site`, unlocks).

| level / representation | legal opportunities | `shared` pools over |
|---|---|---|
| species | `lineage`, `shared` | the whole tree |
| genome — unordered | `copy`, `lineage`, `shared` | one **family** |
| genome — ordered | `copy`, `lineage`, `shared` | the **genome** |
| genome — nucleotide | `site`, `copy`, `lineage`, `shared` | the **genome** |
| sequences (substitution) | `site` | — |

## What each opportunity *means* (the exact totals)

For an event with base rate `r`, the total firing rate over the relevant unit is:

- **`site`** — `r × (#sites)` — one clock per nucleotide.
- **`copy`** — `r × (#copies)` — one clock per gene copy (in this genome).
- **`lineage`** — `r` per lineage — one clock per genome; across `L` live lineages the process fires at
  `r × L`.
- **`shared`** — `r` for the **whole process**, once — one clock total; as the tree/family grows, each
  unit's share shrinks as `r / count`. This is the diversity-dependence-in-disguise case
  (`SharedBirthDeath` is exactly this for speciation).

The single sentence that unifies it: **`lineage` puts one clock on each unit; `shared` puts one clock
on all units together.** Everything else is choosing *which* unit.

## Implementation reality (honest about the work)

Not all rungs cost the same. In rough order of effort:

1. **Species `lineage`/`shared`** — *done.* The species Gillespie loop already tracks the global
   lineage count `n`; `BirthDeath` and `SharedBirthDeath` are the two views. Part 3 here is only a
   surface change: expose `per=` and fold `SharedBirthDeath` into it as a preset.
2. **Genome `copy`/`lineage`** — *done.* `PerCopyRates` / `PerLineageRates`. Part 3 folds them under
   one `Rates(per=…)` and the `--rate-per` flag (already exists).
3. **Genome `shared`** — **the real engine work.** The gene simulator is a joint Gillespie over all
   branches (a Fenwick tree picks the next event across every lineage), but its rate interface,
   `event_weights(genome, branch, time)`, is strictly **per-branch** — a model only sees one genome. A
   `shared` family clock needs the per-family weights to **sum to `r` across branches**, which the
   per-branch view can't express. Two viable mechanisms:
   - **(a) A per-family global pool.** The sim maintains, per family, the branches holding it (and
     their copy counts). A `shared` event contributes *one* weight `r` to a global pool; when it fires,
     the sim localises it (pick a branch ∝ copies, then a copy). Clean, O(1) per fire; the cost is a
     second event pool beside the Fenwick tree.
   - **(b) Normalised per-branch weights.** Each branch `b` holding family `f` gets weight
     `r × copies_b / total_copies_f`, so the branch weights sum to `r`. Keeps one pool, but every
     duplication/loss shifts *all* of `f`'s branch weights (re-normalisation), which is O(branches
     holding `f`) per event.
   Recommendation: **(a)**, the global pool — it matches the mental model ("one clock for the family")
   and keeps the hot loop cheap.
4. **The declarative interface** — the enabling refactor. Instead of `event_weights` returning a
   weight with the opportunity *baked in* (`r × n` for per-copy today), a model declares
   `(event, base_rate, opportunity)` and the **simulator owns the count**. This is what makes the knob
   real: the opportunity becomes data the engine reads, not arithmetic hidden in each class. Species
   and the existing genome rungs are re-expressed through it; `shared` and per-event mixing then fall
   out.
5. **Sequences `site`** — mostly *naming*: substitution/indel rates are already per-site; Part 3 just
   lets the docs and (optionally) the API say `per="site"` so the axis is visible at every level.

## Phasing (each independently shippable)

- **A — Species knob.** ✅ `BirthDeath(per="lineage"|"shared")`; `--per` on `zombi2 species`;
  `SharedBirthDeath` → deprecated preset. Pure surface; byte-identical. *(shipped)*
- **B — Genome knob (existing rungs).** ✅ `Rates(per="copy"|"lineage")` unifies `PerCopyRates` /
  `PerLineageRates` under one model + `--rate-per`/`--per`; `FamilySampledRates` gains `per=` so
  per-family heterogeneity and opportunity are independently selectable (closes the hard-coded-`copy`
  wrinkle above); Rust gate generalised to `isinstance(rates, Rates) and per=="copy"`. Byte-identical.
  *(shipped; the declarative interface — step 4 — is deferred to phase C, where `shared` needs it.)*
- **C — Genome `shared`.** ✅ `Rates(per="shared")` / `zombi2 genomes --per shared`: one tree-wide
  clock per family for duplication and loss (constant total rate; a fire localises to a copy chosen
  uniformly across the family). Implemented as a "shared pool" beside the per-branch Gillespie —
  inert for every other model, so those stay byte-identical (mechanism (a), the global pool). v1 scope:
  unordered genomes, dup/loss (transfer/rearrangements rejected), origination stays per-lineage;
  `FamilySampledRates(per="shared")` and the full declarative interface (step 4) deferred. *(shipped)*
- **D — Per-event mixing.** `Per(unit, rate)` overrides → self-limiting and mixed models.
- **E — Sequences `site`.** Name the axis at the nucleotide level.

## Decisions (resolved with Adrián, 2026-07-15)

1. **Coarsest rung = `shared`; drop `clade`.** "Clade" wrongly implies you must define or select a
   group — you never do. `shared` simply means "one clock for the whole thing" (tree or family). The
   per-unit rungs nest `site ⊂ copy ⊂ lineage`; `shared` sits apart as the "one-for-all" option.
2. **Knob name = `--rate-per`** (the shipped genomes flag), with **`--per`** as an alias. On the API,
   the keyword is `per=`.
3. **Deprecate the preset classes** (`SharedBirthDeath`, `PerCopyRates`, `PerLineageRates`) — one axis,
   one knob, for clarity. Migration is smooth because the *defaults* preserve today's behaviour
   (`BirthDeath(…)` = `per="lineage"`, `Rates(…)` = `per="copy"`); only *explicit* use of a retired
   class name warns, removed in 0.4.0. `SharedBirthDeath` ships in the open PR #146, so it is simply
   born into the knob and deprecated there — since 0.3.0 is unreleased, no released version ever had it
   as the only way.
4. **Genome `shared` follows the representation** (see the validity matrix): **per family** in the
   *unordered* genome (a family is the unit of duplication — the clean `SharedBirthDeath` analogue),
   **per genome** in the *ordered / nucleotide* genomes (duplication is positional). The deeper point:
   *opportunity is scoped by representation.*

## The rule, in one line

> **A rate is `time⁻¹`; a modifier is dimensionless; the opportunity is the unit the clock rides —
> `site ⊂ copy ⊂ lineage`, or `shared` — and it is a knob you turn, the same knob at every level.**
