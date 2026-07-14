# Design: rate vocabulary — opportunity granularity and the word "lineage"

**Status:** accepted (2026-07-14). Finishes the rate-clarity line of work: the primer
([Rates: a primer](../guide/rates.md)) and the modifier layer
([rate-modifiers](rate-modifiers.md)) established *rate = base × opportunities × modifiers, plus
extent*; this note pins the **vocabulary** so no two terms can be conflated, and generalises the
opportunity into a **per-event granularity**.

## Two orthogonal axes (not one)

A rate answers two independent questions, and they must never share a word:

- **Opportunity** — *how many independent chances does the event have right now?* A **count** of
  substrate units.
- **Modifier** — *how much faster/slower is it in this context?* A dimensionless **multiplier**,
  keyed on a family, a lineage, a lineage-pair, a site.

They are genuinely orthogonal, and one event proves it: a **substitution** is counted **per
nucleotide** (opportunity) *and* scaled by a relaxed **clock that varies per lineage** (modifier),
at the same time. So "per nucleotide" and the clock are not two names for one thing — they are
values on two different axes. Keep the axes in separate vocabularies and the whole model reads
cleanly.

## The opportunity axis is a granularity ladder

Opportunity is the **grain of the substrate** the event is counted on. The grains nest — coarse to
fine — and each is a strict container of the next:

```
per lineage  ⊃  per chromosome  ⊃  per copy  ⊃  per nucleotide        (per site, for alignments)
 (the whole        (a karyotype     (a gene       (a base)
  organism)         unit)            copy)
```

`rate = base × (count at the event's grain) × modifiers`. The two ends:

- **per lineage** — count = 1, *independent of contents*. The lineage/organism experiences the event
  as a unit.
- **per copy / per nucleotide** — count scales with the genome's contents (**extensive**).

Because a lineage carries exactly **one** genome, **"per genome" = "per lineage"** — they were always
the same measure. And it is the *same* measure that diversification already uses: speciation and
extinction are per lineage. So the species axis and the genome axis are one axis.

## The word: "lineage" (retire "genome" and "branch")

For maximal clarity we standardise on **one word for the tree-entity — `lineage` — and use it for
the coarsest opportunity grain across every level.**

- **Retire "genome"** from the rate vocabulary: it is conflatable (*per gene in the genome? the
  genome's total?*) and it is just per-lineage at the genome level.
- **Retire "branch"** from the rate vocabulary too, *and reserve it for nothing rate-related*: a
  "branch" is a *finished* tree edge, whereas a rate acts on the *ongoing* lineage and integrates
  over its life. "Per lineage per unit time" is exact; "per branch" faintly implies "one per edge".
  The heterogeneity modifier that used to be called *branch* is renamed **lineage** as well, so the
  single word covers **both** axes and nothing is described two ways.

Why "lineage" over "branch" as the single word:

1. **The ladder reads as containment.** *lineage ⊃ copy ⊃ nucleotide* works because an organism
   contains its genes and bases; *branch ⊃ copy ⊃ nucleotide* is a category slip (a tree edge does
   not contain nucleotides).
2. **Convention.** Birth–death rates are quoted *per lineage per unit time*; nobody says
   "speciation per branch".
3. **Semantics.** A lineage is *ongoing* — the right image for a rate that accrues over time.

## Rename table (all old names kept as deprecated aliases)

| Old | New | Kind |
|---|---|---|
| `--rate-per {copy, genome}` | `--rate-per {copy, lineage}` | CLI value (`genome` aliased) |
| `PerGenomeRates` | `PerLineageRates` | class (alias kept, like `SharedRates`) |
| `BranchRates` | `LineageRates` | class (alias kept) |
| `BranchModifier` | `LineageModifier` | class (alias kept) |
| `--branch-rates FILE` | `--lineage-rates FILE` | CLI flag (`--branch-rates` aliased) |
| "per-branch receptivity", "per branch" (docs) | "per-lineage …" | prose |
| `--orig … (per branch)` help | `… (per lineage)` | prose |

No behaviour changes and no RNG-stream changes: every rename is a pure alias, exactly the trick used
for `SharedRates → PerCopyRates`. `--rate-model` remains a (doubly-)deprecated spelling that maps
through to `--rate-per lineage`.

## Per-event granularity (the general mechanism)

The opportunity grain is a **per-event property**, not a per-model one. The model already mixes
grains — origination is per lineage, D/T/L are per copy, substitution is per nucleotide — so this is
naming an existing fact, then letting it be **chosen** where more than one grain is meaningful.

Consequences:

- `PerCopyRates` and `PerLineageRates` stop being rival classes and become **presets**: "every event
  at copy grain" vs "every event at lineage grain".
- A specific event can be moved off its default grain — e.g. an **inversion counted per lineage**
  (the genome reorganises at a fixed pace) instead of its default per copy (breakpoints ∝ length).
  This generalises to **any** event (duplication, loss, transfer, transposition, translocation, …);
  inversion is only the motivating example.

**Defaults never move.** Every event keeps its canonical grain (rearrangements stay per copy), so
existing runs and the synteny cookbook's per-copy calibration are untouched. Per-event grain is an
**opt-in override**.

### What ships now vs. later

- **Now:** the vocabulary (this note + the rename with aliases + docs + manual) — pure clarity, no
  behaviour change.
- **Deferred until a real use case:** the per-event grain **override API** (letting a chosen event
  switch grain). It is ~15 lines at the weight-emission site — emit a constant (per-lineage) weight
  instead of `× n` (per-copy) for the chosen event — with no engine or Rust impact (the affected
  events are ordered/nucleotide-only, already Python). We build it when a variable-size-genome
  rearrangement study actually needs size-independent rates, so the API is shaped by a real need.

## The rule, in one line

> **Opportunity says *how many* (per lineage → per copy → per nucleotide); modifiers say *how fast,
> here*. One word — *lineage* — for the tree-entity, on both axes; never "genome", never "branch".**
