# Adding a model

ZOMBI2 grows by **adding models, not editing the engine**. The simulator programs against a
small set of interfaces, so new science drops in as a subclass. This page is the contract: what
you implement, and what must be true before it lands.

The value of ZOMBI2 is coherence — a broad library that behaves like one tool. So a model is
not "done" when it runs. It is done when it plugs into the right seam, follows the
[conventions](conventions.md), ships a validation test, and has a catalog page. All four.

## 1. Implement the interface

Pick the level you are adding to and implement its interface. The deep architecture of the
gene-family seams (`Genome` / `RateModel` / `EventSampler`) is in
[Extending ZOMBI2](../guide/extending.md); this is the map.

| Level | Implement | Key method(s) | Passed to | Lives in |
|---|---|---|---|---|
| Species tree | a model class | `validate()`, and `sample_internal_age(u, A)` (backward) *or* forward-simulable rates | `simulate_species_tree` | `species_model.py` |
| Gene-family rates | subclass `RateModel` | `event_weights(genome, branch, time)` | `simulate_genomes` | `rates.py` |
| Genome representation | the `Genome` protocol | `draw_target`, `apply`, `supported_events`, … | `simulate_genomes(genome_factory=…)` | `genome.py` |
| Trait | a duck-typed model | `kind`, `root_value(rng)`, `evolve(state, dt, t0, rng)` | `simulate_traits` | `_traits_impl.py` |
| Substitution | a `SubstitutionModel` | frozen dataclass `name, Q, stationary, alphabet` | `evolve_on_tree` / `SequenceEvolution` | `sequence_sim.py` |
| Molecular clock | subclass `Clock` | `_branch_rate(node, parent_rate, rng)` | `SequenceEvolution` / `clock.scale` | `rate_variation.py` |
| Event sampler | subclass `EventSampler` | `next_waiting_time`, `choose_index` | `simulate_genomes(sampler=…)` | `_sampling.py` |

Most levels are **duck-typed protocols** (implement the methods, no base class needed);
`RateModel`, `EventSampler` are ABCs; `SubstitutionModel` is a frozen dataclass. Add a
`validate()` that raises a clear `ValueError` on bad parameters, and call it on entry.

## 2. Follow the conventions

Every rule in [Conventions](conventions.md) applies — names, the `seed=`/`rng=` signature, the
per-unit-time rate semantics, and (if the model has a CLI surface) the argument-group grammar.
The fastest way to get this right is to copy the nearest existing model and change the science.

## 3. Export it

A public model is reachable **two ways, as the same object**: from the top level
(`zombi2.MyModel`) and from its level's submodule (`zombi2.species.MyModel`). Add it to both
re-export lists and to `__all__` in `zombi2/__init__.py`.

## 4. Add the CLI surface

If the model is user-facing, wire it into its subcommand: a value on the model-selector flag
(`--model`, `--diversification`, …) and a dedicated argument group whose description names that
value (see the [CLI grammar](conventions.md#cli-grammar)). Reuse `--seed`, `-o/--out`,
`--write`; add a runnable line to the command's `EXAMPLES` epilog. The run manifest records new
flags automatically.

## 5. Validate it — the hard rule

> **No model enters the core without an oracle or a statistical test.** A test that only
> asserts "it runs without error" is not validation.

This is the trust guarantee, and it is enforced in review. Pick the strongest that applies (see
[Validation](../validation.md) for worked examples and tolerance conventions):

1. **Oracle** — compare many seeded replicates to a closed-form expectation (e.g. Brownian
   tip variance `σ²·t`; the Mk transition matrix; a reconstructed-age CDF via a KS test).
2. **Statistical reduction** — if there is no closed form, show the model reduces to a known
   case (e.g. one epoch of an episodic process equals constant-rate; a state-independent SSE
   recovers the Yule mean) within a stated Monte-Carlo tolerance.
3. **Invariants + determinism** — always also assert the structural invariants (e.g. a gene
   tree's extant leaves equal the family's extant copy number; a stochastic map tiles its
   branch) and that a fixed seed reproduces the output exactly.

Always include a determinism test. Choose tolerances in the spirit of the existing suite:
`1e-9` for exact identities, ~1 % for recovered frequencies, 5–15 % for Monte-Carlo moments.

## 6. Document it — the catalog page

Every model gets one page in the **model catalog**, all in the same shape, so breadth stays
legible. Copy this template (see [birth–death](../models/birth-death.md) for a filled example):

```markdown
# <Model name>

**What it is.** One paragraph: the process, and one citation.

**When to use it.** The scientific question it answers, and how it differs from its siblings.

## Parameters
| Parameter | Meaning | Default |
| --- | --- | --- |

## Command line
`zombi2 <cmd> --model <x> ...`  + one runnable example.

## Python
A short `simulate_*` snippet.

## Output
Which files/columns it produces (link the conventions).

## Validation
How we know it is correct — name the test and what it checks.

## Reference
The citation.
```

## What belongs in the core

Not every model belongs in the core suite. A model enters the **core** when all of these hold:

1. **It is a genome-evolution model** — one of the levels ZOMBI2 simulates (species trees,
   gene families, genome structure, traits, sequences, clocks) or a *coupling* among them.
2. **It ships validation** — clears the hard rule in step 5.
3. **It follows the conventions** — clears steps 2–4.
4. **It is general** — a reusable model, not a one-off analysis script.

Anything else is an **Extension**: it lives in a separate release with its own manual, sharing
these conventions but not the core's release cadence. Inference (ABC/NDE) and population
genetics are Extensions by design; gene-family coupling was moved to one when it could not yet
clear the bar. Per-lineage diversification (ClaDS) cleared it and is in the core.

This is a **curation call**, made by the maintainers — the point of the bar is that breadth
never comes at the cost of coherence. When unsure, open an issue before writing the model.
