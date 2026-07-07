# Adding a model

ZOMBI2 grows by **adding models, not editing the engine**. The simulator programs against a
small set of interfaces, so new science drops in as a subclass. This page is the contract: what
you implement, and what must be true before it lands.

The value of ZOMBI2 is coherence — a broad library that behaves like one tool. So a model is
not "done" when it runs. It is done when it plugs into the right seam, follows the
[conventions](conventions.md), ships a validation test, and has a catalog page. All four.

## 1. Implement the interface

ZOMBI2 has **four levels** of evolution, each a subpackage under `zombi2/`. Pick the level you
are adding to and implement its interface — a level can offer more than one kind of seam. The
deep architecture of the gene-family seams (`Genome` / `RateModel` / `EventSampler`) is in
[Extending ZOMBI2](../guide/extending.md); this is the map.

| Level | To add… | Implement (file) | Key method(s) | Passed to |
|---|---|---|---|---|
| **species** <br>`zombi2/species/` | a species-tree model | a model class (`model.py`) | `validate()`, and `sample_internal_age(u, A)` (backward) *or* forward-simulable rates | `simulate_species_tree` |
| **genomes** <br>`zombi2/genomes/` | a rate model | subclass `RateModel` (`rates.py`) | `event_weights(genome, branch, time)` | `simulate_genomes` |
| | a genome representation | the `Genome` protocol (`genome.py`) | `draw_target`, `apply`, `supported_events`, … | `simulate_genomes(genome_factory=…)` |
| **traits** <br>`zombi2/traits/` | a trait model | a duck-typed model (`models.py`) | `kind`, `root_value(rng)`, `evolve(state, dt, t0, rng)` | `simulate_traits` |
| **sequences** <br>`zombi2/sequences/` | a substitution model | a `SubstitutionModel` (`models.py`) | frozen dataclass `name, Q, stationary, alphabet` | `evolve_on_tree` / `SequenceEvolution` |
| | a molecular clock | subclass `Clock` (`clocks.py`) | `_branch_rate(node, parent_rate, rng)` | `SequenceEvolution` / `clock.scale` |

Two seams are cross-cutting rather than levels: the **event sampler** (`EventSampler` in
`_sampling.py` — the Gillespie numeric core, shared by every level) and **coevolution**
(`zombi2/coevolve/`), which *couples* two levels rather than adding one — see
[coevolution models](../coevolution_models.md).

Most levels are **duck-typed protocols** (implement the methods, no base class needed);
`RateModel` and `EventSampler` are ABCs, `SubstitutionModel` is a frozen dataclass. Add a
`validate()` that raises a clear `ValueError` on bad parameters, and call it on entry.

## 2. Follow the conventions

Every rule in [Conventions](conventions.md) applies — names, the `seed=`/`rng=` signature, the
per-unit-time rate semantics, and (if the model has a CLI surface) the argument-group grammar.
The fastest way to get this right is to copy the nearest existing model and change the science.

## 3. Export it

A public model is reachable **two ways, as the same object**: from the top level
(`zombi2.MyModel`) and from its level's submodule (`zombi2.species.MyModel`). Add it to both
re-export lists and to `__all__` in `zombi2/__init__.py`. The docs and examples use the
**submodule form** (`from zombi2.species import MyModel`) as canonical — it mirrors the four levels
and keeps `zombi2.<TAB>` from becoming an undifferentiated dump as the catalog grows; the top-level
alias stays for quick interactive use.

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

Models are documented in **family pages**: one page per family of related models (e.g. *Diversification
models*, *Continuous-trait models*, *DNA substitution models*), grouped by level in the catalog nav.
Related variants share a page so breadth stays legible — the way IQ-TREE lists its substitution
models. The page **title is the family name with no level suffix** (the nav section already carries the
level); name each model by its common acronym, expanded on first use — never "(trait)"/"(species tree)".

Copy this shape (see [basic diversification](../models/birth-death.md) for a filled example):

```markdown
# <Family name>

<One paragraph: the shared process, and how the variants extend it.>

| Model | What varies | Reach for it when |    (a scannable overview, one row per model)
| --- | --- | --- |

## The models
### <Model A (ACRONYM)>    (one compact paragraph each: what it is, key params, when to use it)
### <Model B>

## Command line
`zombi2 <cmd> --<selector> <x> ...` + a runnable example per model.

## Python
A `from zombi2.<level> import ...` snippet — the canonical namespaced form.

## Output
Which files/columns it produces (link the conventions).

## Validation
One bullet per model, naming its test and what it checks — the hard rule, kept visible.

## References
The citations.
```

A single-model family (e.g. biogeography's DEC) still gets its own page in this shape.

## What belongs in the core

Not every model belongs in the core suite. A model enters the **core** when all of these hold:

1. **It is a genome-evolution model** — one of the four levels ZOMBI2 simulates (species trees,
   genomes, traits, sequences) or a *coupling* among them.
2. **It ships validation** — clears the hard rule in step 5.
3. **It follows the conventions** — clears steps 2–4.
4. **It is general** — a reusable model, not a one-off analysis script.

Anything else is an **Extension**: it lives in a separate release with its own manual, sharing
these conventions but not the core's release cadence. **Inference** — ABC profile-matching and
ALElite reconciliation likelihood — and **population genetics** are Extensions by design; the
inference code was moved out to `ZOMBI2_FUTURE/` because likelihood is not simulation. By
contrast, per-lineage diversification (ClaDS) is a genome-evolution model that clears the bar,
so it is in the core.

This is a **curation call**, made by the maintainers — the point of the bar is that breadth
never comes at the cost of coherence. When unsure, open an issue before writing the model.
