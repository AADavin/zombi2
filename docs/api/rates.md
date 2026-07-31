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
the command line and a `--params` file take it verbatim:

```python
birth = 1.0 * OnTime({0: 1.0, 3: 0.3})          # Python
```

```bash
zombi2 species --birth "1.0 * OnTime({0: 1.0, 3: 0.3})"
```

A level **rejects** the modifiers it does not support rather than ignoring them, so a run is never
quietly not the model you asked for.

## Scopes

::: zombi2.rates.scope

## Modifiers

A modifier's name begins with the preposition that fixes its family: `On` is a covariate (a
deterministic function of a measured quantity), `By` is an independent i.i.d. draw per unit,
`From` is inherited along a genealogical edge. `DrivenBy` sits outside that scheme deliberately —
it says the factor comes from another simulated level entirely.

::: zombi2.rates.modifiers
