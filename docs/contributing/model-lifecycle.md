# Model lifecycle

ZOMBI2 keeps a **lean, coherent core**: every model in it is reviewed, validated against an
oracle or statistical test, and carries a stable API. To add new models *without* eroding that
guarantee — and *without* making contribution hard — a model moves through three tiers.

## The three tiers

| Tier | Where | Guarantee |
|---|---|---|
| **Core** | `zombi2.*`, in the model catalog | reviewed · validated · stable API |
| **Experimental** | `zombi2.experimental` | usable and iterated on, but **not yet validated or reviewed**; API may change; not in the CLI or catalog |
| **Archive** | `ZOMBI2_FUTURE/` (outside the repo) | removed or deferred — a different paradigm, or not on the near path in |

`zombi2.experimental` is the **inbound door** — a model on its way *into* the core. The archive is
the **parking lot** — a model removed or deferred. A model is in exactly one place at a time.

## The pipeline

```text
idea ──▶ zombi2.experimental ──▶ core
                  │
                  └──────────────▶ ZOMBI2_FUTURE   (if it doesn't make it)
```

### Landing in `experimental`

The bar to *land* is deliberately low — this is a staging area, not the core:

- it runs and has at least a **smoke test**;
- it lives as a module under `zombi2/experimental/` and is re-exported from `zombi2.experimental`
  (never from top-level `zombi2`);
- its constructor calls `zombi2.experimental.warn_experimental("<Name>")`, so a user is told —
  the first time they build one — that it is experimental.

It is **not** wired into the `zombi2` command line or the [model catalog](adding-a-model.md). That
is exactly what keeps the core surface clean, while the model is still fully usable:

```python
from zombi2.experimental import MyNewModel   # explicit — nothing leaks into zombi2.*
```

### Engine-integrated features

Most experimental models are self-contained classes that plug into an existing seam (a rate model, a
trait model) and never touch the engine. A few are different: a new **event kind** — where the engine
itself must learn to *apply* and *reconstruct* a new kind of event — cannot live entirely under
`zombi2/experimental/`, because event dispatch is in the core simulator loop.

For those, **split it**: the **dormant engine capability** goes in the core (the `EventType` member,
the apply method, the dispatch branch, the reconciliation), inert until something fires it; only the
**user-facing surface** — the model that emits the event — stages in `zombi2.experimental` with its
`warn_experimental`. The core stays byte-identical when the feature is off, and the surface still
gates behind the experimental door.

The worked example is **intra-genome gene conversion**, which travelled this whole path: the
`CONVERSION` event, its application, and its reconciliation went into the core engine as a dormant
capability, while the rate model that turns them on first staged in `zombi2.experimental`. Once it
cleared the bar it was **promoted**: the rate folded into `SharedRates(conversion=...)`, the mechanics
became `zombi2.ConversionModel`, it gained a `--conversion` CLI flag, dropped `warn_experimental`, and
the experimental module was removed — the engine capability was already in place.

### Promotion to the core

A model graduates out of `experimental` once it clears the full bar:

- it passes the **[validation hard rule](../validation.md)** — an oracle or statistical test, not
  just a smoke test;
- it has a **[catalog page](adding-a-model.md)** and follows the
  **[conventions](conventions.md)**;
- it meets the **core-vs-extension bar** (it is a genome-evolution model, and general enough to
  belong in the core);
- it has been **reviewed**.

On promotion it moves out of `zombi2.experimental` into a core module, is re-exported from
`zombi2`, drops the `warn_experimental` call, and gains its CLI surface and catalog page.

## The one rule: a waiting room, not a parking lot

`experimental` exists to move models *through*, not to store them. Every model in it should be on
a path to promotion or removal — if it stalls, take it out (to `ZOMBI2_FUTURE` if it is worth
keeping, otherwise delete it). Keeping the namespace small and moving is what stops it from
becoming the graveyard that experimental namespaces tend to become.
