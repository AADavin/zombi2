# Cookbook

A simulator is only as trustworthy as its parameters. Rather than guess a rate out of thin air,
every recipe here does the same thing: take a **pattern in real data**, find the ZOMBI2 knob that
controls it, and tune the knob until simulations reproduce the pattern — learning, along the way,
**what the data can and cannot constrain**. Each recipe is worked end-to-end, so you can copy the
approach for your own system.

| Recipe | Pattern in real data | What you recover | Take-away |
|---|---|---|---|
| [Inversion rate from synteny](synteny.md) | gene-order conservation vs divergence time | genome **inversion rate** | the rate is identifiable; event *size* is not |
| [Can you trust RED?](red-validation.md) | root-to-tip substitution spread in the GTDB tree | whether **RED** recovers node ages (and the realistic clock σ) | RED holds at real heterogeneity — an ordinal proxy, few-% error |

!!! abstract "Parameters pinned so far"
    - **Budding-yeast inversion rate — ≈ 3–4 × 10⁻⁴ per gene · Myr** ([synteny](synteny.md);
      consistent with Fischer et al. 2006).
    - **Archaeal clock heterogeneity — σ ≈ 0.56** under the uncorrelated-lognormal relaxed clock
      ([can you trust RED?](red-validation.md)).

    Drop each into the matching ZOMBI2 knob. This list grows as recipes are added.
