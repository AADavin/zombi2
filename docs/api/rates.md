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
it says the number comes from an evolved value read on the lineage — a level grown before this run, a
level growing beside it, or another object at the same level (a trait can drive a second trait). On a
rate or an extent that number multiplies; on `transfer_to` it is a weight normalised across the
candidates, which is why that one takes the modifier on its own — `transfer_to = 1.0 * DrivenBy(...)`
is an error.

### Writing your own

Every engine takes a fixed set of modifiers and refuses the rest, because one it never reads would
return its default factor of 1.0 and give a run that is quietly not the model you asked for. A
modifier of your own opens that gate by naming the engines you implemented it for:

```python
from zombi2.rates.modifiers import Modifier

class OnLogTime(Modifier):
    implemented_for = ("species",)
    def factor(self, *, time: float = 0.0, **_) -> float:
        return 1.0 / (1.0 + time)

species.simulate_species_tree(birth=2.0 * OnLogTime(), n_extant=20, seed=1)
```

Each engine supplies a different context, and `sequences` does not take a modifier of your own at all
— it reads its modifiers itself rather than through the rate, so one it did not ship could never be
called, and it refuses rather than ignoring you:

| Engine | Context passed to `factor` |
|---|---|
| `species` | `time`, `lineages`, `diversity` |
| `genomes.family` | `time`, `lineages`, `copies`, `drivers` |
| `genomes.ordered` | `time`, `lineages`, `copies`, `chromosomes`, `drivers` |
| `genomes.nucleotide` | `time`, `lineages`, `copies`, `chromosomes`, `drivers` |
| `traits.continuous` | `time`, `lineages`, `diversity`, `inherited`, `drivers` |
| `traits.discrete` | `time`, `lineages`, `drivers` |
| `joint` | `time`, `lineages`, `diversity`, `drivers` |

Take `**_` and default every keyword you read. Naming an engine is a claim you are making, not a check
the library can do for you — everything you have not named still refuses your modifier, by name.

If your factor varies continuously with `time`, override `next_change` to return the next point at
which it should be re-evaluated. The engine holds a rate constant between events, so without it the
curve is frozen at whatever it was when the last event fired.

The rate *text* grammar (a `--birth` flag, a `--params` file) knows only the built-in names, so a
modifier of your own is Python-only, as an object you construct has to be.

**Three worked examples** — a rate following a measured curve, density dependence in the gene pool,
and rearrangement scaling with the karyotype — are in
[Appendix A, "Writing your own"](../rates.md#writing-your-own).

::: zombi2.rates.modifiers

## Mappings

What a `DrivenBy` carries — the shape that turns the driver's value into a number.

::: zombi2.rates.mapping

## Values

::: zombi2.rates.values

## Verbs

::: zombi2.rates.verbs

## Clades

::: zombi2.rates.clade
