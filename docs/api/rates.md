# zombi2.params

Every event fires at a **rate**, and every rate is written the same way at every level:

```
effective rate  =  scope(base) × modifiers
```

The **base** is the speed of one event, in inverse time. The **scope** answers "per what?" — how
many copies, lineages or sites the event applies to right now — and is the entry point the rate is
written from. The **modifiers** are dimensionless context multipliers, and they multiply. There is
no `per=` argument: the scope lives on each rate.

A modifier is not a class you call: it is what a **verb** records. `scaled_by` multiplies the base,
`set_by` replaces it, and `weighted_by` compares the candidates of a choice. Two drivers are written
so often that each has a verb of its own, which is the only spelling for it: `varying_among` for a
value drawn for each unit of one kind, `changing_at` for the run's clock. Verbs chain, and their
factors multiply.

That expression is not Python syntax the CLI translates; it is *the* way a rate is written, and
the command line and a `--params` file take it verbatim. The same text, three places:

```python
from zombi2 import species
from zombi2.params import PerLineage

birth = PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})
species.simulate_species_tree(birth=birth, n_extant=10, seed=1)
```

```toml
# params.toml
birth = "PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})"
n-extant = 10
```

```bash
zombi2 species out/ --birth "PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})" --n-extant 10 --seed 1
zombi2 species out/ --params params.toml --seed 1
```

Every name the written form may call is importable from `zombi2.params`, so a snippet pastes across
unchanged. Four qualifiers are tolerated where Python needs one — `mod.`, `modifiers.`, `scope.` and
`scopes.` — and nothing else, so `rates.PerCopy(...)` or a dotted `zombi2.params.PerCopy(...)` is
refused: the parser reads a whitelist of names, not attribute paths. A bare number stays a bare
number everywhere (`birth = 1.0`, `--birth 1.0`).

A level **rejects** the modifiers it does not support rather than ignoring them, so a run is never
quietly not the model you asked for.

## Scopes

::: zombi2.params.scope

## Rates, extents and choices

The three chainable parameters, and the verbs that live on them. A rate is written from its scope, an
extent from `Extent(...)`, and a choice from `Recipients()` — the last two only when a verb is
chained onto them, since a bare distribution is already an extent and `"uniform"` is already a
choice.

::: zombi2.params.rate

::: zombi2.params.extent

::: zombi2.params.choice

## Modifiers

A modifier's **kind** says who produces its number, and there are four (SPEC §5):

| Kind | The factor is… | Written |
|---|---|---|
| covariate | a deterministic function of a measured quantity | `changing_at({…})`, `scaled_by(TotalDiversity(cap=…))` |
| drawn | an i.i.d. draw, one per unit — no memory | `varying_among(unit, dist)` |
| inherited | the parent's, perturbed — continuous memory | `varying_among(unit, Drift(dist))` |
| driven | the state of another simulated thing, read as the run walks the tree | `scaled_by`, `set_by`, `weighted_by` |

A **driven** value comes from a level grown before this run, a level growing beside it, another
object at the same level (a trait can drive a second trait), or the tree itself (`Clade`). Which of
the three verbs you write follows from what you attach it to: on a rate the number multiplies the base
(`scaled_by`) or replaces it (`set_by`); on an extent it multiplies only, an extent being an absolute
size with no base to replace. On `transfer_to` it is a weight normalised across
the candidates (`weighted_by`), which is why that one is written from `Recipients()` with no base —
`transfer_to = PerLineage(1.0).weighted_by(...)` is an error.

The **unit** a drawn or inherited value is attached to is an argument, not a class, so a draw among
families and a draw among lineages are one class and two cells of a grid. The unit is **plural**,
because a value varies *among* families rather than being counted *per* one.

### Writing your own

Every engine takes a fixed set of modifiers and refuses the rest, because one it never reads would
return its default factor of 1.0 and give a run that is quietly not the model you asked for. A
modifier of your own opens that gate by naming the engines you implemented it for:

```python
from zombi2 import species
from zombi2.params import PerLineage
from zombi2.params.modifiers import Modifier
from zombi2.params.rate import Rate

class OnLogTime(Modifier):
    implemented_for = ("species",)
    def factor(self, *, time: float = 0.0, **_) -> float:
        return 1.0 / (1.0 + time)

birth = Rate(2.0, PerLineage, (OnLogTime(),))
species.simulate_species_tree(birth=birth, n_extant=20, seed=1)
```

The verbs build the built-in modifiers, so there is no verb that attaches one of yours; a rate
carrying it is built from `Rate` directly, with the scope class and the modifiers in the order they
should be drawn in.

Each engine supplies a different context, and `sequences` does not take a modifier of your own at all
— it reads its modifiers itself rather than through the rate, so one it did not ship could never be
called, and it refuses rather than ignoring you:

| Engine | Context passed to `factor` |
|---|---|
| `species` | `time`, `lineages`, `diversity` |
| `genomes.family` | `time`, `lineages`, `copies`, `drivers` |
| `genomes.ordered` | `time`, `lineages`, `copies`, `chromosomes`, `drivers` |
| `genomes.nucleotide` | `time`, `lineages`, `copies`, `chromosomes`, `drivers` |
| `traits.continuous` | `time`, `lineages`, `diversity`, `drivers` |
| `traits.discrete` | `time`, `lineages`, `drivers` |
| `joint` | `time`, `lineages`, `diversity`, `drivers` |

Take `**_` and default every keyword you read. Naming an engine is a claim you are making, not a check
the library can do for you — everything you have not named still refuses your modifier, by name.

If your factor varies continuously with `time`, override `next_change` to return the next point at
which it should be re-evaluated. The engine holds a rate constant between events, so without it the
curve is frozen at whatever it was when the last event fired.

The rate *text* grammar (a `--birth` flag, a `--params` file) knows only the built-in names, so a
modifier of your own is Python-only, as an object you construct has to be.

**A worked example** — `OnCrowding`, a death rate that rises as the tree fills — is in
[Appendix A, "Writing your own"](../rates.md#writing-your-own), with the two things a modifier of
your own has to provide and the one it may.

::: zombi2.params.modifiers

## Mappings

What a driven parameter carries — the shape that turns the driver's value into a number.

::: zombi2.params.mapping

## Drivers

::: zombi2.params.driver

## Verbs

::: zombi2.params.verbs
