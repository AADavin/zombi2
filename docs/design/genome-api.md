# Genome API — design target

**Status: the design to build.** The target for rewriting `zombi2/genomes`, the detailed consequence
of `SPEC.md` for the genome level. Designed with Adrián on 2026-07-18. **Not built yet** — today's code
ships the `RateModel` zoo (§ *What to delete*). Parallels `species-api.md`; read that first.

---

## The problem it fixes

`zombi2/genomes` ships a whole rate hierarchy — `RateModel`, `Rates`, `PerCopyRates`, `PerLineageRates`,
`FamilySampledRates`, `ModifiedRates`, `LineageRates`, plus a `Modifier` family — and `simulate_genomes`
takes *both* a rate object *and* keyword rates, and picks the resolution by passing a **class**
(`genome_factory=UnorderedGenome`). Two axes are tangled in that hierarchy (C7.5): per-copy-vs-per-lineage
(the count) and shared-vs-per-family (heterogeneity).

## Three functions, layered over one core

Resolution is not an argument; it is the choice of function, because each resolution has a genuinely
different (growing) argument set:

```
genomes.simulate_unordered(tree, …)     # D, T, L, O — a multiset of gene families (Chapter 5)
genomes.simulate_ordered(tree, …)       # + position: inversions, transpositions, chromosomes (Chapter 6)
genomes.simulate_nucleotide(tree, …)    # + DNA: genes, intergenes, indels (Chapter 6)
```

**unordered ⊂ ordered ⊂ nucleotide.** They are LAYERS over one shared core, not three copies:
`simulate_unordered` *is* the core; the upper two delegate to it and add their layer. The behaviour
therefore cannot drift — only the shared arguments repeat across the three signatures, pinned by a test
that the shared params match. Each returns the matching object (unordered / ordered / nucleotide result).

## Rates: keyword arguments, `count(base) × modifiers`

The events are keyword rates — no `RateModel` object:

```python
genomes.simulate_unordered(tree, duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
                           initial_families=20, seed=1)
```

Rates follow the **same cross-level grammar** as species (`SPEC §5`): a rate is an optional **count
wrapper** around a base, optionally times **modifiers**.

- **Count** answers *per what*, and the default is right per event: D/T/L are **per copy**, origination is
  **per lineage**. Override with a wrapper: `origination = PerLineage(0.5)` is the default;
  `origination = Global(0.5)` makes it constant over time. `loss = PerLineage(0.25)` counts loss per
  lineage instead of per copy. (`Global`, capitalised — `global` is a Python keyword.)
- **Modifiers** multiply the base (`*`) and are dimensionless: `loss = 0.25 * Time({…})` etc.
- **Naming rule:** *"per" is reserved for counts.* A count is `PerCopy` / `PerLineage` / `PerSite`; a
  modifier never starts with "per" (so the per-family modifier below is `ByFamily`, not `PerFamily`).

## Per-family heterogeneity (the genome-specific piece; resolves C7.5)

This is a **modifier**, cleanly separate from the count. Two layers, which compose:

- **Base rates, per event, independently (default).** A bare number is shared by every family; a
  `ByFamily` modifier makes each family draw its own, independently per event (a family that loses fast is
  not automatically duplicating fast). At the extreme, set them explicitly:
  ```python
  loss = 0.25 * mod.ByFamily(spread=0.5)                 # each family its own loss rate
  families = [dict(duplication=0.5, transfer=0.8, loss=0.3), …]   # explicit specific families
  ```
- **A per-family speed (correlated), with its own name.** One factor per family scaling *all* its rates
  together — a "fast" family is fast at everything. **`family_speed=Speed(spread=0.5)`** (decided
  2026-07-18: `Speed`; NOT "PerFamily"). Effective loss = (base loss) × (family speed).

## Transfers: rate + mechanics (resolves C7.6, C7.8)

The transfer **rate** is a rate like the others. What is special is what happens when it fires — the
**mechanics**, expressed as arguments. **`TransferModel` dissolves into these** (C7.8: mechanics become
arguments, not a passed object):

```python
transfer=0.1,
transfer_to="uniform",   # or "distance": closer relatives more likely (with a strength knob)
replacement=False,       # additive (default) vs overwrite an existing copy
self_transfer=False,
```

The per-lineage biases (C7.6) split by the ontology: a **donor weight IS a modifier** on the transfer
rate (a lineage that donates more just has `transfer = 0.1 * mod.PerLineage(...)`), so "emission" needs
no special word. A **recipient weight is NOT a modifier** — it biases who receives once a transfer fired,
so it is a weight in the `transfer_to` rule (the mechanic). "Receptivity" → **recipient weight**.

## Decided (2026-07-18)

- **Length distributions are configurable, per event type.** Indel and rearrangement/segment lengths
  choose a **distribution**, per event type — `inversion_length=Geometric(mean=5)` (a `<event>_length=`
  argument taking a distribution constructor: `Geometric`, `Poisson`, …), replacing the single global
  geometric knob.

## Still to design

- **Conversion** mechanics (which copy overwrites which; directionality).
- The **resolution-specific arguments**: rearrangements and chromosomes (ordered — number, topology,
  fission/fusion/translocation, C6.10/C6.11); genes/intergenes and indels (nucleotide, C10.x).
- **Decided (2026-07-18): the count-wrapper namespace is `scope`** — `scope.Global`, `scope.PerCopy`,
  `scope.PerLineage`, `scope.PerSite` (reads "global scope", "per-copy scope"). NOT `mod` (wrappers wrap;
  they do not multiply). (`Speed` also decided.)
- Per-event count **override** exact syntax (rides on the rate via the wrapper).

## What to delete / change in `zombi2/genomes`

- Delete the `RateModel` hierarchy and `genome_factory`; there is no rate object and no class-passing.
- Three entry points (`simulate_unordered` / `simulate_ordered` / `simulate_nucleotide`) over one core.
- `TransferModel` → arguments. Per-family heterogeneity → the `ByFamily` modifier + a named speed.
- Rates use the cross-level `count(base) × modifiers` grammar; "per" reserved for counts.
