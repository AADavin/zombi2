# Bounding gene-family growth

A family's copy number is a birth–death process in disguise: with duplication rate `d`
above loss rate `l` its expected size grows like `e^{(d−l)t}` without bound. **Both**
duplication and transfer create copies, so both must be reined in.

## Hard cap: `max_family_size` (recommended)

A single ceiling on family size, enforced across **all** copy-creating events:

```python
z.simulate_genomes(tree, duplication=0.5, transfer=0.2, loss=0.1, origination=0.3,
                   max_family_size=0.5)     # cap = round(0.5 · N_species)
```

- an **integer** is an absolute cap (`max_family_size=20`);
- a **float** is a fraction/multiple of the number of species
  (`max_family_size=0.5` → half the tips; `2.0` → twice the tips).

Duplication is rate-suppressed once a family reaches the cap; an additive transfer that
would overflow is turned into a **replacement** (net zero). The result is a family size that
never exceeds the cap.

## Soft cap: `carrying_capacity`

A logistic, duplication-only alternative living on the rate model: the per-copy duplication
rate is scaled by `max(0, 1 − n/K)`, so family size settles *around* `K` with a proper
stationary distribution.

```python
z.UniformRates(duplication=0.5, loss=0.1, origination=0.3, carrying_capacity=20)
```

!!! note "Which to use?"
    `carrying_capacity` shapes duplication smoothly but does **not** bound transfers; use
    `max_family_size` when you need a firm ceiling that also accounts for transfer-driven
    growth. They compose — you can set both.

## The safety guard

If a run still diverges (e.g. `allow_self=True` with no cap), the simulator raises a clear
error rather than hanging, pointing you to `max_family_size`.
