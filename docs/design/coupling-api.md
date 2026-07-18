# Coupling API — design target (Part III)

**Status: the design to build.** The target for Part III of the manual and for rewriting the `coevolve`
machinery. Designed with Adrián on 2026-07-18. **Not built yet.** Parallels the four level design docs
(`species-api.md`, `genome-api.md`, `sequence-api.md`, `trait-api.md`); read those and `SPEC.md` §2–§4 first.

---

## The one idea: a rate driven by another level

A **coupling** is Ch2's definition made literal: *a parameter that reads its value from another level
instead of a number you type.* There is exactly one mechanism, a modifier:

```python
loss = 0.25 * mod.Driven(source, mapping)
```

`Driven` reads the driver's value on each lineage and multiplies the base rate by the mapped factor. That
is the whole of Part III. Everything below is where the driver comes from.

## Conditioned vs joint is one distinction: can the driver be grown first?

**This is the spine of the chapter — not "does it change the tree".** The organizing question is whether
the driver can be simulated on its own and handed over:

- **Conditioned** — the driver **can** be grown first, so `source` is a **file**. Two commands, ordered.
  ```python
  habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.1, seed=1)
  habitat.write("habitat.tsv")
  genomes.simulate_unordered(tree,
      loss = 0.25 * mod.Driven("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0}), seed=2)
  ```
- **Joint** — the driver **cannot** be grown first, because it is entangled with what it drives, so
  `source` is a **live level name** and both are grown in one call.
  ```python
  joint.simulate(
      birth = 1.0 * mod.Driven("trait", {"small": 1.0, "large": 2.0}),  # trait drives speciation
      death = 0.2,                                                       # death can be Driven too → BiSSE
      trait = traits.discrete(states=["small", "large"], switch=0.1),    # grown WITH the tree
      n_tips = 100, seed = 1)
  ```

Same modifier, `mod.Driven`; the only difference is `source` = filename (conditioned) vs level-name
(joint). So **conditioning and joint models are one chapter, not two** — one mechanism split by
file-vs-live driver. (Manual: Part III 9 *Conditioning* + 10 *Joint* collapse to a single **Coupling
levels** chapter; **11 Nulls** stays separate — see below.)

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

## What `joint.simulate` grows

`joint.simulate(...)` runs a single Gillespie over both levels' events at once: speciation uses the
trait-dependent `birth`/`death`, and trait-change events evolve the trait on the growing tree. It produces
**both** levels (the grown tree + the trait history). The v1 live pairs:

- a trait drives speciation → `P(Species, Traits)` (BiSSE / MuSSE / QuaSSE / HiSSE)
- gene content drives speciation → `P(Species, Genomes)`

(Sequence-driven diversification and tree-fixed mutual driving are deferred to experimental, SPEC §10.)

## Nulls (in design this session — provisional)

**The different question 11 asks:** a pattern on a tree is not evidence of a coupling until you know what
*no* coupling looks like *on that same tree* — the tree manufactures apparent associations through shared
ancestry and shared timing. A null generates the "no coupling" world to compare against. Three kinds:

- **`independent`** (neutral) — cut the link entirely; both levels evolve unlinked on the same tree. The
  plain baseline: is there *any* coupling?
- **`cid`** (character-independent) — the SSE null: let diversification vary across the tree, but *not*
  tied to the focal trait. Does *this* trait explain more than a nameless background rate shift? (Guards
  the BiSSE false-positive: neutral traits look significant when rate heterogeneity is unmodeled.)
- **`shuffle`** (timing / marginal) — keep each level exactly as observed (same values, same counts, same
  tree) but scramble the pairing. Is the *association* real, or just marginals + tree?

Proposed surface — the result of a coupled run can produce its own matching null distribution:

```python
real = joint.simulate(birth=1.0*mod.Driven("trait", {"small":1, "large":2}),
                      trait=traits.discrete(states=["small","large"], switch=0.1), n_tips=100, seed=1)
null = real.null("cid", replicates=100)   # same processes, coupling replaced by the null → a distribution
```

CLI: `--null {independent,cid,shuffle} --replicates N`.

## Still to design

- Nulls surface (`real.null(...)` vs a top-level `joint.null(...)`); which statistic is reported; whether
  the conditioned (fixed-tree) case uses the same `.null()` path.
- `Driven` naming (`Driven` vs `DependsOn`); the level-name referencing convention for live drivers.
- The `traits.discrete(...)` **process spec** (a description passed to `joint`) vs `simulate_discrete(...)`
  (a runner) — confirm the spec/run split.
- Driving *both* birth and death cleanly (state-dependent extinction).

## What to delete / change

- The `coevolve` command and its `--couple driver:target` grammar → the `Driven` modifier + `joint`
  command. Conditioned couplings fold into the **target level's** command (`--loss-driven-by file`);
  only genuinely-joint (live-driver) models keep the dedicated `joint` command.
- The null layer (`.null(kind=)` / `--null`) attaches to coupled runs, not a separate subsystem.
