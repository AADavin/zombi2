# Adding a model

ZOMBI2 grows by **adding models, not editing the engine**. The simulator programs against a
small set of interfaces, so new science drops in as a subclass. This page is the contract: what
you implement, and what must be true before it lands.

The value of ZOMBI2 is coherence — a broad library that behaves like one tool. So a model is
not "done" when it runs. It is done when it plugs into the right seam, follows the
[conventions](conventions.md), ships a validation test, and has a catalog page. All four.

New models usually **stage in [`zombi2.experimental`](model-lifecycle.md) first** — a low-bar
inbound door — and graduate into the core once they clear everything on this page; a model that
already clears the full bar can also land in the core directly. Either way, **this page is the bar**.
Most models "drop in as a subclass" and never touch the engine, but a few — a new **event kind**, say
— also need a small dormant capability in the core engine; those stage only their *user-facing
surface* in `zombi2.experimental` (see [engine-integrated features](model-lifecycle.md#engine-integrated-features)).

## 1. Implement the interface

ZOMBI2 has **four levels** of evolution, each a subpackage under `zombi2/`. Pick the level you
are adding to and implement its interface — a level can offer more than one kind of seam. The
map:

| Level | To add… | Implement (file) | Key method(s) | Passed to |
|---|---|---|---|---|
| **species** <br>`zombi2/species/` | a species-tree model | a model class (`model.py`) | `validate()`, and `sample_internal_age(u, A)` (backward) *or* forward-simulable rates | `simulate_species_tree` |
| **genomes** <br>`zombi2/genomes/` | a rate model | subclass `RateModel` (`rates.py`) | `event_weights(genome, branch, time)` | `simulate_genomes` |
| | a genome representation | the `Genome` protocol (`genome.py`) | `draw_target`, `apply`, `supported_events`, … | `simulate_genomes(genome_factory=…)` |
| **traits** <br>`zombi2/traits/` | a trait model | a duck-typed model (`models.py`) | `kind`, `root_value(rng)`, `evolve(state, dt, t0, rng)` | `simulate_traits` |
| **sequences** <br>`zombi2/sequences/` | a substitution model | a `SubstitutionModel` (`models.py`) | frozen dataclass `name, Q, stationary, alphabet` | `evolve_on_tree` / `SequenceEvolution` |
| | a molecular clock | subclass `Clock` (`clocks.py`) | `_branch_rate(node, parent_rate, rng)` | `SequenceEvolution` / `clock.scale` |

A trait model normally implements `evolve(state, dt, t0, rng)` — evolve the state forward by a
time step `dt`. If it instead needs the **whole branch** (a per-branch regime history, say), it
may implement the optional hook `evolve_branch(node, x, rng) -> (new_value, map_or_None)` and the
engine calls that instead of `evolve` when it is present (`getattr(model, "evolve_branch", None)`).
`MultiOptimumOU` in `zombi2/traits/models.py` is the worked example — it integrates an OU process
across the regime segments painted on each branch.

Two seams are cross-cutting rather than levels: the **event sampler** (`EventSampler` in
`_sampling.py` — the Gillespie numeric core, shared by every level) and **coevolution**
(`zombi2/coevolve/`), which *couples* two levels rather than adding one — see
[coevolution models](../guide/coevolution.md).

Most levels are **duck-typed protocols** (implement the methods, no base class needed);
`RateModel` and `EventSampler` are ABCs, `SubstitutionModel` is a frozen dataclass. Add a
`validate()` that raises a clear `ValueError` on bad parameters, and call it on entry.

A **species** model additionally declares a `_caps` (a `SpeciesCaps`) class attribute — its
growth engine (`thinning`/`gillespie`), and which modes it supports (backward reconstruction,
ghost grafting, `n_tips` stop mode, incomplete sampling, and any forward-only features). The
simulator dispatches on this declared capability instead of `isinstance`, so a model that omits
it raises a clear error rather than silently routing into the wrong growth loop. Copy the nearest
model (e.g. `ClaDS`) for the shape.

### The gene-family seams in depth

The genomes level is the most developed set of seams, and it is worth seeing in detail. ZOMBI2
is deliberately **interface-first**: one Gillespie simulator programs only against three
protocols, so new science drops in as subclasses rather than edits to the engine.

| Protocol | Responsibility |
|---|---|
| `Genome` | how a genome is represented and mutated (`UnorderedGenome`, `OrderedGenome`, …) |
| `RateModel` | how fast events happen — turns a genome into weighted candidate events |
| `EventSampler` | the numeric hot loop (waiting time + weighted choice); the Rust-swap point |

The rule of thumb: **the simulator loop, the sampler, the profile matrix and the output
never change** when you add a representation, a rate model or an event type.

#### Add a rate model

A rate model implements `event_weights(genome, branch, time)`, returning a list of
`(event, family_or_None, rate)` candidates. `family=None` means "act on a uniformly chosen
copy"; a specific family means "weight the target by this family's own rate".

```python
from zombi2 import EventType
from zombi2.genomes import RateModel, EventWeight, simulate_genomes

class PerLineageRates(RateModel):
    """D/T/L totals independent of genome size."""
    def __init__(self, dup, trans, loss, orig):
        self.d, self.t, self.l, self.o = dup, trans, loss, orig
    def event_weights(self, genome, branch, time):
        out = []
        if genome.size() > 0:
            out += [EventWeight(EventType.DUPLICATION, None, self.d),
                    EventWeight(EventType.TRANSFER,    None, self.t),
                    EventWeight(EventType.LOSS,        None, self.l)]
        out.append(EventWeight(EventType.ORIGINATION, None, self.o))
        return out

genomes = simulate_genomes(tree, PerLineageRates(0.5, 0.2, 0.5, 0.4), seed=1)
```

`RateModel.bind(rng, max_family_size)` is called once per run — override it for stateful
models (e.g. per-family sampling, or a future non-independence model that reads
`genome.presence_vector(order)` to couple families).

#### Add a genome representation

Implement the `Genome` interface (queries, `draw_target`, `apply`, the transfer handoff,
`clone_reminting`, `supported_events`) and pass a factory:

```python
from zombi2.genomes import simulate_genomes

genomes = simulate_genomes(tree, rates,
                           genome_factory=lambda ids: MyGenome(ids))
```

`OrderedGenome` is a worked example: it added gene order plus inversion/transposition with
no change to the engine. Declaring the new events in `supported_events()` is what lets the
simulator fire them.

#### Add an event type

Append a member to `EventType`, emit it from a rate model, and handle it in a genome's
`apply` + `supported_events`. The loop keeps only events a genome supports, so other
representations ignore it automatically.

#### Built through these seams

The interface-first design has let these models drop in as subclasses, without touching the
engine:

- **Non-independence** — a future `RateModel` whose gain/loss rates would read the
  genome's presence vector, so functionally coupled families gain/lose together.
- **Ghost lineages** and **fossilized birth–death** — a forward species-tree simulator that
  retains extinct lineages (see [species trees](../guide/species-trees.md)).
- **Gene length / intergenes** and **genome-wise rates** — further `Genome` / `RateModel`
  subclasses.

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

Copy this shape (see [species trees](../guide/species-trees.md) for a filled example):

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

## PR checklist

Before opening a pull request, confirm the four gates — a reviewer will:

- [ ] **Seam (§1)** — subclasses the right interface; the simulator loop, sampler, profile matrix
  and output are untouched (or, for a new event kind, only a small dormant core capability is added).
- [ ] **Conventions (§2–4)** — names, the `seed=` / `rng=` signature, per-unit-time rate semantics,
  and (if user-facing) the CLI argument-group grammar; exported both ways and added to `__all__` (§3).
- [ ] **Validation (§5)** — an oracle *or* statistical-reduction test **plus** a determinism test.
  "It runs without error" is not validation.
- [ ] **Catalog page (§6)** — a family page in the documented shape, wired into the docs nav.

New or not-yet-fully-validated models can stage in [`zombi2.experimental`](model-lifecycle.md)
first and graduate once every box is ticked.

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