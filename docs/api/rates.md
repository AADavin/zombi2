# zombi2.rates

Every event fires at a **rate**, and every rate is written the same way at every level:

```
effective rate  =  scope(base) × modifiers
```

The **base** is the speed of one event, in inverse time. The **scope** answers "per what?" — how
many copies, lineages or sites the event applies to right now — and wraps the base. The
**modifiers** are dimensionless context multipliers, and they multiply. There is no `per=`
argument: the scope lives on each rate.

That expression is not Python syntax the CLI translates; it is *the* way a rate is written, and
the command line and a `--params` file take it verbatim. The same text, three places:

```python
from zombi2 import species
from zombi2.rates import modifiers as mod

birth = 1.0 * mod.OnTime({0: 1.0, 3: 0.3})
species.simulate_species_tree(birth=birth, n_extant=10, seed=1)
```

```toml
# params.toml
birth = "1.0 * OnTime({0: 1.0, 3: 0.3})"
n-extant = 10
```

```bash
zombi2 species out/ --birth "1.0 * OnTime({0: 1.0, 3: 0.3})" --n-extant 10 --seed 1
zombi2 species out/ --params params.toml --seed 1
```

The `mod.` qualifier Python needs is optional in the other two, so a snippet pastes across
unchanged. A bare number stays a bare number everywhere (`birth = 1.0`, `--birth 1.0`).

A level **rejects** the modifiers it does not support rather than ignoring them, so a run is never
quietly not the model you asked for.

## Scopes

::: zombi2.rates.scope

## Modifiers

A modifier's name begins with the preposition that fixes its family: `On` is a covariate (a
deterministic function of a measured quantity), `By` is an independent i.i.d. draw per unit,
`From` is inherited along a genealogical edge. `DrivenBy` sits outside that scheme deliberately —
it says the factor comes from another simulated level entirely.

### Writing your own

Every engine takes a fixed set of modifiers and refuses the rest, because one it never reads would
return its default factor of 1.0 and give a run that is quietly not the model you asked for. A
modifier of your own opens that gate by naming the engines it is wired for:

```python
from zombi2.rates.modifiers import Modifier

class OnLogTime(Modifier):
    wired_for = ("species",)
    def factor(self, *, time: float = 0.0, **_) -> float:
        return 1.0 / (1.0 + time)

species.simulate_species_tree(birth=2.0 * OnLogTime(), n_extant=20, seed=1)
```

The engine names are `species`, `genomes.family`, `genomes.ordered`, `genomes.nucleotide`,
`sequences`, `traits.continuous`, `traits.discrete` and `joint`. Take `**_` and default every keyword
you read: the engine passes whatever context it happens to have, and the set differs between them.
Naming an engine is a claim you are making, not a check the library can do for you — everything you
have not named still refuses your modifier, by name.

The rate *text* grammar (a `--birth` flag, a `--params` file) knows only the built-in names, so a
modifier of your own is Python-only, as an object you construct has to be.

::: zombi2.rates.modifiers
