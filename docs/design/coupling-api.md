# Coupling API — design target (Part III)

**Status: building.** The target for Part III of the manual and for rewriting the `coevolve`
machinery. Designed with Adrián on 2026-07-18. **Built (2026-07-20):** the `DrivenBy` modifier
(`zombi2/rates/modifiers.py`) + `Table`/`Curve`/`Scalar` mappings (`zombi2/rates/mapping.py`);
**conditioned** trait→gene loss (`zombi2/rates/driver.py` + the per-lineage genome retrofit); and
**joint** discrete-trait→speciation (BiSSE/MuSSE) via `zombi2.joint.simulate_joint` (there is no
`coupling` package — conditioning folds into the target level, so `joint` is the only engine). Next:
joint gene-content→speciation. Parallels the four level design docs (`species-api.md`,
`genome-api.md`, `sequence-api.md`, `trait-api.md`); read those and `SPEC.md` §2–§4 first.

---

## The one idea: a rate driven by another level

A **coupling** is Ch2's definition made literal: *a parameter that reads its value from another level
instead of a number you type.* There is exactly one mechanism, a modifier:

```python
loss = 0.25 * mod.DrivenBy(source, mapping)
```

`DrivenBy` reads the driver's value on each lineage and multiplies the base rate by the mapped factor. That
is the whole of Part III. Everything below is where the driver comes from.

**`DrivenBy` targets a rate (a "how often") and multiplies — value-driving is deferred (v1).** Driving a
*value* (a "what") — the one case being an OU optimum, SPEC §4's "gene content drives a trait's optimum" —
is **deferred to experimental for v1.** A destination lives on a real line and is *set or shifted*, not
scaled by a positive factor, so it does not fit `DrivenBy`'s multiply. When it lands it gets its **own verb**
(a value-reader, e.g. `mod.Reads(source, mapping)` — "Option B"), never an overload of `DrivenBy`. For v1,
`DrivenBy` is rate-only: one clean mechanism, no exception.

## Conditioned vs joint is one distinction: can the driver be grown first?

**This is the spine of the chapter — not "does it change the tree".** The organizing question is whether
the driver can be simulated on its own and handed over:

- **Conditioned** — the driver **can** be grown first, so `source` is a **file**. Two commands, ordered.
  ```python
  habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.1, seed=1)
  habitat.write("out", outputs=("events",))                      # → out/trait_events.tsv
  genomes.simulate_genomes_unordered(tree,
      loss = 0.25 * mod.DrivenBy("out/trait_events.tsv",
                                 {"aquatic": 3.0, "terrestrial": 1.0}), seed=2)
  ```
  (The driver is written through the standard `write(dir, outputs=(...))` spine, like every level's
  outputs; the file it produces is `trait_events.tsv`, the trait event log (a `root` row then the switches), which a driven run replays against the shared tree — a per-branch `node · start · end · state` segment
  table — exact even when a discrete driver switches mid-branch.)

  **In-memory shortcut.** In a single Python session you can skip the file and pass the grown result
  object straight to `DrivenBy` — `DrivenBy(habitat, {…})` — which is *the same conditioning* (the
  driver grown first, held fixed), just handed over in memory. The file round-trip is lossless, so the
  two give an identical run; the file matters for reproducibility, the CLI, and cross-process pipelines.
- **Joint** — the driver **cannot** be grown first, because it is entangled with what it drives, so
  `source` is a **live level name** and both are grown in one call.
  ```python
  joint.simulate_joint(
      birth = 1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.0}),  # trait drives speciation
      death = 0.2,                                                       # death can be DrivenBy too → BiSSE
      trait = traits.discrete(states=["small", "large"], switch=0.1),    # grown WITH the tree
      n_extant = 100, seed = 1)
  ```

Same modifier, `mod.DrivenBy`; the only difference is `source` = filename (conditioned) vs level-name
(joint). So **conditioning and joint models are one chapter, not two** — one mechanism split by
file-vs-live driver. (Manual: **all of Part III — 9 *Conditioning* + 10 *Joint* + 11 *Nulls* — collapses
to a single "Coupling levels" chapter**, with nulls as its closing section. See below.)

### Why "changes the tree" is NOT the framing

For v1 the *only* live-driver models we ship are tree-changing (a trait or gene content drives
speciation), so "joint = changes the tree" is *accidentally* exact. It must not become the chapter's
organizing idea: the moment a **tree-fixed** joint model lands (a trait and gene content driving each
other, SPEC §4), that framing breaks. The future-proof line is *can the driver be grown first?* — which
already covers the tree-fixed case with no reframing. v1 shipping only tree-changing joints is a **scope
note**, mentioned once, never the spine.

## The mapping is the "response" (from the coevolve grammar)

`mapping` says how the driver's value becomes a factor. It reuses the four response shapes of the coevolve
grammar (`coevolve-grammar.md`):

- **Table** — a discrete driver → a dict: `{"aquatic": 3.0, "terrestrial": 1.0}`.
- **Curve** — a continuous driver → a function/curve: `lambda x: exp(0.5 * x)` (e.g. QuaSSE, `dN/dS` vs a
  continuous trait).
- **Scalar** — a single multiplier when the driver is already 0/1 or binary.
- **Jump** — the response fires *at an event* (a burst of gene change at each split), not continuously.

## What `joint.simulate_joint` grows

`joint.simulate_joint(...)` runs a single Gillespie over both levels' events at once: speciation uses the
trait-dependent `birth`/`death`, and trait-change events evolve the trait on the growing tree. It produces
**both** levels (the grown tree + the trait history). The v1 live pairs:

- a trait drives speciation → `P(Species, Traits)` (BiSSE / MuSSE / QuaSSE / HiSSE)
- gene content drives speciation → `P(Species, Genomes)`

(Sequence-driven diversification and tree-fixed mutual driving are deferred to experimental, SPEC §10.)

## Nulls — a recipe, not an API (the Coupling chapter's closing section)

**Decided:** there is **no null function and no null subsystem.** A null is just the simulator run under a
model where the coupling isn't real, built from primitives that already exist. This **retires the separate
Nulls chapter** — it becomes the closing section of the Coupling chapter.

The one idea worth teaching: *a pattern on a tree is not evidence of a coupling until you know what no
coupling produces on that same tree* — the tree manufactures associations through shared ancestry and
timing, so the baseline must be simulated on the same tree. Three recipes, each a one-line change to the
rate:

```python
loss = 0.25 * mod.DrivenBy("habitat.tsv", {"cave": 4.0, "surface": 1.0})   # the coupling under test

loss = 0.25                                    # independent null — drop the coupling
loss = 0.25 * mod.ByLineage(spread=0.5)        # CID null — background heterogeneity, NOT the trait
loss = 0.25 * mod.DrivenBy(shuffle("habitat.tsv"), {...})   # shuffle null — permute the pairing
```

> **Build note (2026-07-20).** Of the three, *independent* and *shuffle* run today; the **CID null does
> not yet** — `mod.ByLineage` is not wired on genome rates (no clean-core engine threads a per-lineage
> `bylineage` draw yet). Wiring `ByLineage` on the genome events is a small follow-up (the per-lineage
> machinery now exists from `DrivenBy`). `shuffle()` stays a recipe (permute the driver file's per-node
> state assignment), not an API — a tiny `tools.shuffle` helper is the most it might warrant.

Then a plain `for seed in range(100)` gives the distribution. CID is *literally* `ByLineage` (rate varies
across the tree, not by the trait), reusing the clock-collapse modifier. **ZOMBI2 does not own the test
statistic or p-value** — choosing the association measure and the test is inference, the user's job; the
simulator only generates the baseline. The only thing that might warrant a `tools` utility is a tiny
`shuffle()` helper for permuting a tip-value file.

## Decided (2026-07-18)

- **Gene-content driver source.** When gene content drives speciation, `source` names a summary of gene
  content: **presence of a named family** — `DrivenBy("genes:toxin", {"present": 2.0, "absent": 1.0})`, a
  Table — or **total gene count** — `DrivenBy("genes:count", curve)`, a Curve. Richer multi-family profiles
  are deferred.
- **Driving both birth *and* death.** Trivially supported: both are rates, so `birth = 1.0 * DrivenBy(...)`
  and `death = 0.2 * DrivenBy(...)` each work independently — full state-dependent diversification (BiSSE's
  λ *and* μ).
- **Process spec vs runner.** `traits.discrete(...)` / `traits.continuous(...)` are thin **process specs**
  (parameters bundled, unexecuted); `traits.simulate_discrete(tree, ...)` is the runner. `joint` takes the
  spec and runs it as the tree grows.

## Still to design (naming only)

- **Decided (2026-07-18): the modifier is `DrivenBy`** — `loss = 0.25 * mod.DrivenBy("habitat", {…})`
  reads "loss, driven by habitat" (the *source* is the driver; the *rate* is what's driven).
- Open: the level-name referencing convention for live drivers (`"trait"`, `"genes:toxin"`, `"genes:count"`).

## What to delete / change

- The `coevolve` command and its `--couple driver:target` grammar → the `DrivenBy` modifier + `joint`
  command. Conditioned couplings fold into the **target level's** command (`--loss-driven-by file`);
  only genuinely-joint (live-driver) models keep the dedicated `joint` command.
- **No null layer.** Nulls are recipes built from existing modifiers (drop / swap for `ByLineage` /
  shuffle the driver file), not a `.null()` API. ZOMBI2 never computes the test statistic.
