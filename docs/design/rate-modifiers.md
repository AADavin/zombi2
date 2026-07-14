# Design: a generic rate-modifier layer

**Status:** proposed (2026-07-13). Phase 3 of the rate-clarity work — after the Phase 1 CLI/API
cleanup (`PerCopyRates`, `--rate-per`) and the Phase 2 docs rewrite ([Rates: a primer](../guide/rates.md)).

## Goal

One `Modifier` abstraction — **a context-keyed multiplier on an event's rate** — that every level
composes, replacing the bespoke `LineageRates` wrapper and the special-cased family / receptivity
paths with a single stackable mechanism. This makes the code say exactly what the primer says:
*modifiers multiply, and they are the same idea everywhere*. It also unlocks the two flexibility
asks (a family's transfer biased between two clades; family × branch on sequence speed) as ordinary
compositions rather than new one-off classes.

## The protocol

```python
class Modifier(ABC):
    time_dependent = False
    def bind(self, rng, tree=None) -> None: ...          # optional setup (autocorr needs the tree)
    @abstractmethod
    def factor(self, event, family, branch, time) -> float: ...   # a multiplier; default 1.0
    def refresh_times(self, t0, t1) -> list: return []   # when the factor changes on its own
```

A modifier reads whatever context keys it cares about (`family`, `branch`, a `(donor, recipient)`
pair, a site) and ignores the rest. Composition is `∏ factor`, order-independent, default 1.

## The one subtlety: two application points

A rate in ZOMBI2 is realized in **two stages**, and every modifier attaches to exactly one:

1. **Emission** — `RateModel.event_weights(genome, branch, time)` produces the per-branch weight for
   each event on a *donor* lineage. Branch, family, and size/carrying-capacity modifiers act here:
   they scale *how often* an event fires.
2. **Recipient choice** — for a **transfer**, the recipient lineage is chosen *after* emission,
   weighted by `TransferModel.receptivity`. A **lineage-pair** transfer modifier (donor → recipient)
   *must* act here, because the recipient does not exist yet at emission time.

So the generic layer has two seams that share one `factor(context)` shape but differ in where they
apply:

| Seam | Applied by | Modifiers that live here |
|---|---|---|
| **Emission** | `ModifiedRates` (a `RateModel`) | `LineageModifier`, `FamilyModifier`, size/carrying-capacity |
| **Recipient choice** | `TransferModel` | `PairModifier` (donor→recipient), per-recipient receptivity |

This is the honest version of "one Modifier idea": **one concept, one `factor` signature, two
attachment points** — not one code path pretending transfer recipients are known at emission.

## The composer

`ModifiedRates(base, modifiers)` is a `RateModel`:

- `event_weights` = each base weight × `∏` emission-modifier factors;
- `target_params` / `establishment_probability` delegate to `base`;
- `refresh_times` = union of base + modifiers (so autocorrelated / trait-drifting factors stay exact);
- `bind` binds base and every modifier.

## Mapping the existing classes (back-compat + byte-identity)

- **`LineageRates(base, …)`** — keep the public class; reimplement it as
  `ModifiedRates(base, [LineageModifier(…)])`. **Byte-identity constraint:** `LineageModifier` must
  draw from the RNG in the *exact* order `LineageRates` does today (autocorrelated: tree-preorder
  normals at `bind`; per-branch: lazy i.i.d. per branch). The existing byte-identity tests are the
  guard.
- **`FamilySampledRates`** — keep the public class (the ZOMBI1-style *base* form). Add an equivalent
  `FamilyModifier` (the *overlay* form: base × per-family multiplier, drawn once per family and
  cached). This is what makes **per-genome × per-family** expressible — the corner the old enum
  couldn't name.
- **`PerCopyRates` / `PerLineageRates`** — unchanged; they are the two *base* rate laws (the
  opportunity axis). Modifiers stack on either.
- **`TransferModel.receptivity`** — becomes the first `RecipientModifier`; generalise it to accept a
  per-`(donor, recipient)` **pair** matrix in addition to per-recipient weights.

## The flexibility asks fall out

- **A family's transfer biased between two clades** → `PairModifier({(donor, recipient): factor})` on
  `TransferModel`, optionally family-scoped (recipient-choice seam).
- **Transfer highways** (preferential HGT routes, à la Beiko et al. 2005, "highways of gene sharing")
  → the same `PairModifier`: a highway is a high-weight cell — or a clade→clade *block* — in the
  donor×recipient matrix, and several highways are several blocks. Because the pair factor is a
  `factor(event, family, branch, time)` like any modifier, a highway can later be made
  **epoch-specific** (open only in a time window) via the `refresh_times` hook, and **family-scoped**
  (a highway that carries only some families) by keying on `family`. This is the motivating use case
  for the recipient seam.
- **Family × branch on sequence speed** → stack `FamilyModifier` + `LineageModifier` on the
  substitution rate (emission seam). Composes for free once both exist.

## Rust fast path — untouched

`_rust.eligible` still fires *only* for a plain `PerCopyRates` (no modifiers, no carrying capacity,
no rearrangements). Any `ModifiedRates`, or a receptivity **pair** matrix, routes to the Python
engine — exactly as family / emission rates do today. Per-recipient receptivity keeps its existing
Rust parity; the pair matrix is Python-only to start.

## Decision taken

**Hold byte-identity when porting `LineageRates`** (replicate the exact RNG draw order) rather than
accept a re-seeded stream. Reproducibility across a seed is load-bearing, and the existing
byte-identity tests should stay green with zero churn. Slightly fiddlier, worth it.

## Implementation order (one pass, each step keeps the suite green)

1. `Modifier` protocol + `ModifiedRates` composer (emission seam); no behaviour change yet.
2. `LineageModifier`; re-express `LineageRates` on top of it (byte-identical); tests green.
3. `FamilyModifier` (overlay form); wire `--family-rates` / composition; per-genome × per-family.
4. `PairModifier` + `TransferModel` pair matrix (recipient seam); the transfer-between-clades ask.
5. Sequence side: expose `FamilyModifier` × `LineageModifier` on the substitution rate.
