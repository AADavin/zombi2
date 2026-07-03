# Ghost lineages

A reconstructed species tree contains only the **sampled, extant** lineages — that is the
whole point of the backward, conditioned birth–death process behind
[`simulate_species_tree`](species-trees.md). But the real diversification process also
produced lineages that **went extinct before the present** (or were simply not sampled).
Those lineages are invisible in the reconstructed tree, yet they were part of the history —
for example, an extinct lineage could have been the donor of a horizontal transfer.

`z.add_ghost_lineages` **un-prunes** the tree: it grafts those dead ("ghost") lineages back
on, so the tree again reflects the full process rather than just the survivors.

!!! note "Un-pruning vs. simulating forward"
    Un-pruning is for when you *already have* (or want) a reconstructed tree conditioned on
    `N` extant tips and need the dead lineages added back. If you instead want the complete
    tree from the start, `simulate_species_tree(..., direction="forward")` grows it forward
    with extinct lineages included natively — no un-pruning needed. See
    [species trees](species-trees.md). The two are distributionally equivalent for the
    reconstructed part; pick by which you have in hand.

This page is the how-to. For the exact conditional law behind it — the `λ(t)·E(t)`
attachment intensity, the h-transform sampler, and why un-pruning is exact — see
[Ghost lineages — how it works](../ghost_lineages.md).

## Basic use

Pass the **same model** you used to build the tree:

```python
model = z.BirthDeath(birth=1.0, death=0.5)
tree = z.simulate_species_tree(model, n_tips=50, age=5.0, seed=1)

z.add_ghost_lineages(tree, model, seed=7)   # grafts ghosts in place, and returns the tree
```

The function mutates `tree` **in place** (and also returns it, so you can chain).

## What you get

Ghost lineages attach along each edge and each roots a birth–death subtree conditioned on
leaving no sampled descendant. Every new node is marked `is_extant=False`; the grafted dead
leaves are named `e*` (extinct) and the new internal nodes `i*` — ZOMBI2's standard
convention (extant leaves stay `n*`). The **sampled leaves are left untouched**, so pruning
back to the extant tips recovers the original reconstructed tree exactly.

```python
extant = [n for n in tree.leaves() if n.is_extant]        # the original sampled tips (n*)
ghosts = [n for n in tree.leaves() if not n.is_extant]    # the grafted dead tips (e*)
print(len(extant), "extant +", len(ghosts), "ghost leaves")

print(tree.to_newick())   # the Newick now includes the e* / i* ghost nodes
```

## When do ghosts appear?

Only where lineages could actually have been lost — i.e. the process must allow extinction
or incomplete sampling:

| Model | Ghosts? |
| --- | --- |
| `Yule` / `BirthDeath(death=0)` with full sampling | **none** — nothing goes extinct |
| `BirthDeath(death>0)` | yes — extinct lineages are grafted back |
| `EpisodicBirthDeath(..., sampling_fraction<1)` | yes — extinct **and** unsampled lineages |

On a tree with no possible extinction, `add_ghost_lineages` returns it unchanged.

## Supported models

`BirthDeath`/`Yule` and `EpisodicBirthDeath` (time-varying rates **and** incomplete
sampling with `sampling_fraction < 1`). Pass the model instance you built the tree with, so
the ghost process uses the matching rates.

## Choosing a sampler (`method=`)

Each ghost subtree is grown conditioned on leaving no sampled descendant. Two equivalent
samplers are available:

- **`"rejection"`** (default) — simple and exact: grow a birth–death subtree and reject any
  that leaves a survivor. When extinction or incomplete sampling is heavy, it may retry
  often. Two guards apply only to this method: `max_subtree_size` (reject a runaway subtree)
  and `max_attempts` (cap on retries per attachment).
- **`"htransform"`** — rejection-free, via Doob's h-transform. Faster in the heavy-retry
  regime, and statistically equivalent.

```python
z.add_ghost_lineages(tree, model, method="htransform", seed=7)
```

## Why un-prune?

Ghost lineages let you study what the reconstructed tree hides: extinct or unsampled
lineages as sources/sinks of horizontal transfer, the effect of sampling on downstream
inference, and the difference between the complete and reconstructed histories. The result
is an ordinary `Tree` — the extra tips are just leaves with `is_extant=False`.

See also the [cookbook](../cookbook.md#add-ghost-extinct-lineages) for the short recipe and
[species trees](species-trees.md) for the underlying model.
</content>
