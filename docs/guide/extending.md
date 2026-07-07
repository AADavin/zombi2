# Extending ZOMBI2

ZOMBI2 is deliberately **interface-first**. One Gillespie simulator programs only against
three protocols, so new science drops in as subclasses rather than edits to the engine.

| Protocol | Responsibility |
|---|---|
| `Genome` | how a genome is represented and mutated (`UnorderedGenome`, `OrderedGenome`, …) |
| `RateModel` | how fast events happen — turns a genome into weighted candidate events |
| `EventSampler` | the numeric hot loop (waiting time + weighted choice); the Rust-swap point |

The rule of thumb: **the simulator loop, the sampler, the profile matrix and the output
never change** when you add a representation, a rate model or an event type.

## Add a rate model

A rate model implements `event_weights(genome, branch, time)`, returning a list of
`(event, family_or_None, rate)` candidates. `family=None` means "act on a uniformly chosen
copy"; a specific family means "weight the target by this family's own rate".

```python
from zombi2 import EventType
from zombi2.genomes import RateModel, EventWeight, simulate_genomes

class PerGenomeRates(RateModel):
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

genomes = simulate_genomes(tree, PerGenomeRates(0.5, 0.2, 0.5, 0.4), seed=1)
```

`RateModel.bind(rng, max_family_size)` is called once per run — override it for stateful
models (e.g. per-family sampling, or the future Potts model, which will read
`genome.presence_vector(order)` to couple families).

## Add a genome representation

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

## Add an event type

Append a member to `EventType`, emit it from a rate model, and handle it in a genome's
`apply` + `supported_events`. The loop keeps only events a genome supports, so other
representations ignore it automatically.

## Built through these seams

The interface-first design has let these models drop in as subclasses, without touching the
engine:

- **Non-independence (Potts model)** — a `RateModel` whose gain/loss rates read the genome's
  presence vector, so functionally coupled families gain/lose together.
- **Ghost lineages** and **fossilized birth–death** — a forward species-tree simulator that
  retains extinct lineages (see [species-tree models](../species_tree_models.md)).
- **Gene length / intergenes** and **genome-wise rates** — further `Genome` / `RateModel`
  subclasses.
