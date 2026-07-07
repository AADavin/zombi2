# Tools

`zombi2.tools` is the **analysis / interop complement** to the simulator. `zombi2` proper
*generates* data — species trees, gene families, traits, sequences. `zombi2.tools` provides
**bounded, validated computations on those outputs**: likelihoods, scores, statistics,
distances, and format conversions.

Everything here is a **distinct, labelled surface**. Nothing in `zombi2.tools` is re-exported
into the top-level `zombi2` namespace — you import tools explicitly:

```python
from zombi2.tools import reconciliation_likelihood
```

That separation is deliberate: the simulation core stays small and stable, and the tools layer
grows alongside it without enlarging the core's surface.

Tools with a command-line surface live under the `zombi2 tools` command group, mirroring this
namespace — e.g. `zombi2 tools reconcile` evaluates the reconciliation likelihood of a gene tree
you already have. Run `zombi2 tools -h` to list them.

## The bar

A tool is admitted only if it is a **computation with a right answer** that is:

- **bounded** — a well-posed calculation on a ZOMBI2 output, not an open-ended search;
- **built** — implemented, not aspirational;
- **validated** — checked against a closed-form oracle or a reference implementation.

This bar is what keeps `tools` trustworthy. It also draws a firm line:

!!! note "What is *not* a tool"
    **Open-ended inference frameworks — ABC, MCMC, neural density estimation — are out of
    scope for `zombi2.tools`.** Fitting a model's parameters to empirical data by simulation-
    based inference is a different kind of activity (a methodology, not a bounded computation),
    and it remains **deferred to a future inference release**. A tool computes a number you can
    check; an inference framework explores a space you cannot enumerate.

Contrast with [`zombi2.experimental`](../contributing/model-lifecycle.md): experimental *models*
are shipped so you can iterate on them but have **not** cleared validation; tools are the
opposite — **stable and validated**, just adjacent to (rather than part of) the simulation core.

## Available tools

| Tool | What it computes | Page |
| --- | --- | --- |
| **Reconciliation likelihood** (ALElite) | the ALE marginal likelihood `P(gene tree \| species tree, DTL rates)` of a simulated gene family | [Reconciliation likelihood](reconciliation-likelihood.md) |

## Roadmap

The following are **candidate tools — not yet built**. They are listed to show the intended
shape of the layer (bounded, validated computations on ZOMBI2 outputs), not to promise a
schedule. Each would enter only after clearing the bar above.

- **Gene-family gain/loss likelihood (Count-lite).** A per-family two-state (presence/absence)
  gain–loss likelihood on the phylogenetic-profile matrix — the fast, likelihood-based
  complement to the full reconciliation likelihood, for scoring *origins* without gene trees.
- **Phylogenetic-profile / co-occurrence statistics.** Bounded summaries of a `Profiles.tsv`:
  pairwise presence/absence co-occurrence scores, phylogenetic-profile distances, and clustering
  of co-occurring families.
- **Interop exporters.** Convert ZOMBI2 reconciliations and gene trees into the formats other
  tools read — **recPhyloXML**, and **NHX** for ALE / GeneRax / Count — so a simulated scenario
  can be fed straight into an external reconciliation pipeline.
- **Benchmarking tree / event distances.** Compare a simulated (true) tree or reconciliation to
  an inferred one: Robinson–Foulds and related topology distances, plus event-level
  (D/T/L/S mapping) agreement — the metrics a benchmark of an inference method needs.
