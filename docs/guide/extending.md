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
import zombi2 as z
from zombi2.events import EventType

class GenomeWiseRates(z.RateModel):
    """D/T/L totals independent of genome size."""
    def __init__(self, dup, trans, loss, orig):
        self.d, self.t, self.l, self.o = dup, trans, loss, orig
    def event_weights(self, genome, branch, time):
        out = []
        if genome.size() > 0:
            out += [z.EventWeight(EventType.DUPLICATION, None, self.d),
                    z.EventWeight(EventType.TRANSFER,    None, self.t),
                    z.EventWeight(EventType.LOSS,        None, self.l)]
        out.append(z.EventWeight(EventType.ORIGINATION, None, self.o))
        return out

genomes = z.simulate_genomes(tree, GenomeWiseRates(0.5, 0.2, 0.5, 0.4), seed=1)
```

`RateModel.bind(rng, max_family_size)` is called once per run — override it for stateful
models (e.g. per-family sampling, or the future Potts model, which will read
`genome.presence_vector(order)` to couple families).

## Add a genome representation

Implement the `Genome` interface (queries, `draw_target`, `apply`, the transfer handoff,
`clone_reminting`, `supported_events`) and pass a factory:

```python
genomes = z.simulate_genomes(tree, rates,
                             genome_factory=lambda ids: MyGenome(ids))
```

`OrderedGenome` is a worked example: it added gene order plus inversion/transposition with
no change to the engine. Declaring the new events in `supported_events()` is what lets the
simulator fire them.

## Add an event type

Append a member to `EventType`, emit it from a rate model, and handle it in a genome's
`apply` + `supported_events`. The loop keeps only events a genome supports, so other
representations ignore it automatically.

## What's coming through these seams

- **Non-independence (Potts model)** — a `RateModel` whose gain/loss rates read the genome's
  presence vector, so functionally coupled families gain/lose together (in design).
- **Ghost lineages** — a forward species-tree simulator retaining extinct lineages, which
  also unlocks fossilized birth–death (see the
  [species-tree roadmap](../species_tree_models.md)).
- **Gene length / intergenes**, **genome-wise rates** — further `Genome` / `RateModel`
  subclasses.
