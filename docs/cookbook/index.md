# Cookbook

A simulator is only as trustworthy as its parameters. Rather than guess a rate out of thin air,
every recipe here does the same thing: take a **pattern in real data**, find the ZOMBI2 knob that
controls it, and tune the knob until simulations reproduce the pattern — learning, along the way,
**what the data can and cannot constrain**. Each recipe is worked end-to-end, so you can copy the
approach for your own system.

| Recipe | Pattern in real data | What you recover | Take-away |
|---|---|---|---|
| [Inversion rate from synteny](synteny.md) | gene-order conservation vs divergence time | genome **inversion rate** | the rate is identifiable; event *size* is not |
| [RED node-age benchmark](red-benchmark.md) | *(ground-truth stress test)* | how well **RED** recovers node ages | a method's limits, not a parameter |

!!! abstract "Rates pinned so far"
    **Budding-yeast inversion rate — ≈ 3–4 × 10⁻⁴ per gene · Myr** ([derived here](synteny.md);
    consistent with Fischer et al. 2006). Drop it straight into the matching rate flag. This list
    grows as recipes are added.
