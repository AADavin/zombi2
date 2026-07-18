# Species-tree API — design target

**Status: the design to build.** This is the agreed target for rewriting `zombi2/species`. It is the
detailed, level-specific consequence of `SPEC.md` (one process per level; collapse the model zoo;
`scope(base) × modifiers`). Designed with Adrián on 2026-07-18 while writing the Species Trees
chapter. It is **not built yet** — today's code still ships the seven-class zoo (§ *What to delete*).

Principle: **concepts → code → chapter.** Build this, then the chapter documents it.

---

## The problem it fixes

Today `zombi2/species` ships **seven** process classes — `BirthDeath`, `Yule`, `ClaDS`,
`EpisodicBirthDeath`, `DiversityDependent`, `SharedBirthDeath`, `CladeShiftBirthDeath` — and
`simulate_species_tree(model, ...)` **requires** you to construct and pass one. That reifies a zoo the
ontology says does not exist: there is **one** species process (lineages are born and die), and
everything else is that process with the rate depending on something.

---

## One entry point

```python
simulate_species_tree(
    birth,                 # a Rate: an optional scope wrapper around a base, optionally × modifiers
    death=0.0,             # a Rate; default 0.0 = pure birth (Yule)
    *,
    n_tips=None,           # stop at a number of tips ...
    age=None,              # ... or at an age
    age_type="crown",      # "crown" (default) or "stem"
    mass_extinctions=None, # [(time, fraction), ...]  — interventions
    sampling=1.0,          # fraction of extant species observed
    fossils=0.0,           # fossil recovery rate along branches
    seed=None,
) -> Tree
```

No model object. `Yule` is `death=0`. The seven classes are gone.

**Module structure.** Each level is its own module — `zombi2.species`, `zombi2.genomes`,
`zombi2.sequences`, `zombi2.traits` — plus `zombi2.modifiers` for the rate modifiers. The package
mirrors the four-level spine. Call it as `zombi2.species.simulate_species_tree(...)`. **No top-level
re-exports** (`z.simulate_species_tree` goes away): one canonical path per name, so nothing drifts.
Open sub-question: whether the function is `simulate_species_tree` (explicit) or just
`species.simulate(...)` (the module already names the level).

## Rates: `scope(base) × modifiers`

`birth` and `death` each hold a **Rate**, built the cross-level way (`SPEC §5`): an optional **scope
wrapper** around a base number, optionally times dimensionless **modifiers**.

- **Scope** — *per what?* Species has one unit, the lineage, so the bare number is per lineage (the
  default) and there is **no `per=` argument**. To make a rate a single shared clock, wrap it:
  `birth = Global(1.0)` (constant total; `Global` capitalised — `global` is a Python keyword). This
  replaces the earlier `per="global"` argument, and it is better: `Global` tags exactly *which* rate is
  global, so `birth=Global(1.0), death=0.3` is unambiguous.
- **Modifiers** multiply the base (`*`) and are dimensionless: `1.0 * Diversity(cap=100)`. `*` is the
  only composition operator (a rate is `time⁻¹`, a modifier is dimensionless, so multiplying two rates is
  impossible by construction).
- **Naming rule:** "per" is reserved for scopes. A scope is `PerCopy` / `PerLineage` / `PerSite`; a
  modifier never starts with "per".

Modifiers live in **`zombi2.modifiers`** (never top-level `z`). Named by *what the rate depends on*:

| Modifier | The rate depends on | Example |
|---|---|---|
| `Time({t: factor, …})` | time (skyline / episodic) | `1.0 * Time({0: 1.0, 3: 0.3})` |
| `Diversity(cap=K)` | standing diversity (slows toward `K`) | `1.0 * Diversity(cap=100)` |
| `Inherited(spread=σ)` | ancestry (drifts at each split, descendants inherit; = ClaDS) | `1.0 * Inherited(spread=0.2)` |

Modifiers are **relative factors**, so `Time({0: 1.0, 3: 0.5})` on base `1.0` reads as absolute
(`1.0`, then `0.5`) and on base `2.0` scales (`2.0`, then `1.0`). Stack with `*`; `birth` and `death`
are bent independently. `Global` gives one budget for the whole tree instead of per lineage
(Ch2's scope).

## Interventions vs observation (not rates)

- **`mass_extinctions=[(3.0, 0.75)]`** — at time 3.0, 75% of living lineages die. A point-in-time
  intervention on the *process*, not a rate. (Clade shift, deferred below, is the same shape.)
- **`sampling=0.5`** — observe only half the surviving species (classic incomplete sampling).
- **`fossils=0.1`** — recover fossils along the branches, a bare rate. **v1: fossils are a side output**
  (the sampled lineages + their ages); the lineage is *not* removed, and fossils do *not* appear in the
  reconstructed tree. Full fossilised birth–death (fossils as dated sampled-ancestors in the tree) is a
  future `tools` command.

## The direction is inferred, never chosen

The user does not pick forward vs backward. It follows from what was asked:

- constant or `Time` rates, with at most extant `sampling` → the reconstructed tree can be **sampled
  backward** (fast);
- `Diversity`, `Inherited`, `mass_extinctions`, or `fossils` need the process to actually play out →
  **forward**.

This is why the chapter does not open on forward-vs-backward: it is a consequence, not a knob. The
complete-vs-reconstructed tree, backward sampling, and ghost lineages are a later section (§4 of the
chapter), still to be designed.

## The bridge table (literature → this API)

The single place the acronyms appear (SPEC §4). Section headings say what each does.

| Literature | What it does | This API |
|---|---|---|
| constant-rate birth–death | fixed rates | `birth=1.0, death=0.3` |
| Yule | pure birth | `death=0.0` |
| skyline / episodic | rates change in time | `birth=1.0 * Time({…})` |
| diversity-dependent | slows with diversity | `birth=1.0 * Diversity(cap=…)` |
| ClaDS | rates drift, inherited | `birth=1.0 * Inherited(spread=…)` |
| mass extinction | a cull at a time | `mass_extinctions=[(t, f)]` |
| incomplete sampling | see a fraction | `sampling=ρ` |
| fossilised birth–death | fossils of the dead | `fossils=ψ` |

## Outputs

- **Both trees, by default:** the `_extant` tree (survivors) **and** the `_complete` tree (with the
  dead), in Newick, using the `_complete` / `_extant` naming (SPEC §D6). Tips are labelled
  **extant / extinct / unsampled** so the three are told apart.
- **The event log, always written:** every speciation and extinction with its time and lineages. This
  is the ground truth the simulator exists to record, so it is not opt-in.
- **Fossils, optional:** written only when `fossils` is set — the sampled fossil lineages and their ages.

The full cross-level list of output files lives in Appendix B of the manual.

## Still to design

- **Clade shift** — a lineage and its descendants switch rates partway through. It must name a clade
  that does not exist yet, so it is referenced by `(time, a rule for picking the lineage)`, not by a
  `birth=` value. Adrián's observation: it is the same *point-in-time intervention* shape as a mass
  extinction, so the two may share one form. **Decided 2026-07-18: set aside for v1** — not in the first
  release; when it returns it uses the `(time, rule)` shape, likely sharing the mass-extinction form.
- **Ghost lineages — set aside for v1 (Adrián, 2026-07-18).** They run on a different paradigm; we will
  do them, just not now.
- **Sampling / reconstructed tree / fossils — decided (2026-07-18, v1).** Incomplete extant sampling (ρ)
  prunes the reconstructed tree; **fossils are a side output** (sampled lineages + ages), the lineage is
  **not** removed, and fossils do **not** appear in the reconstructed tree. Full fossilised birth–death
  (fossils as dated sampled-ancestors *in* the tree) is deferred to a future `tools` command.

## What to delete / change in `zombi2/species`

- Delete the seven process classes; there is no `model` object.
- `simulate_species_tree` takes `birth`/`death` (numbers or `number × modifier`) plus the arguments
  above; it infers the direction.
- Add `zombi2.modifiers` with `Time`, `Diversity`, `Inherited` (dimensionless, `*`-composable).
- Replace the `--rate-model` selection mess with the scope wrappers `Global` / `PerLineage` (retire
  `"shared"`, SPEC §12); there is no `per=` argument.
