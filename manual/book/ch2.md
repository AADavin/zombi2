# Species trees

The species tree is the backbone every other level runs on, so it is where almost every workflow begins.

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give a **speciation rate** (`birth`) and an **extinction rate** (`death`), and every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest. ZOMBI2 does not support multifurcations.

The two rates determine the overall shape of the tree and how quickly diversity increases. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process. With extinction above zero, the lineages that died are still kept: they are in the complete tree, which is the tree the next level runs along, and that is what lets a gene be transferred out of a lineage that dies at a later time. The extant tree (also reconstructed tree) is the tree that can be recovered in theory from the lineages that have survived until the present time.

![A species tree created by a birth–death process. Every lineage alive at a given moment has the same chance per unit time of speciating or of dying. The lineages that died are drawn dashed and stop where they died; the survivors reach the present. Both are in the complete tree, and only the solid ones are in the extant tree.](figures/species_tree.pdf){width=100%}

There are two main ways to determine the stopping time of a birth–death process. Either you stop it at a given time (`total_time`), or you stop it when the number of lineages alive reaches a certain number (`n_extant`).[^stopping]

`n_extant` limits the number of leaves by construction; `total_time` does not control directly the number of leaves, which is dependent on the difference between the speciation and extinction rates. A run that passes **100,000 standing lineages** stops with an error by default. If you really want to simulate very large trees, then just change the `max_lineages` parameter, or set it to `None` to put no limit.

`total_time` is not conditioned on survival. A run in which no lineages are alive at the stopping time raises an error.

[^stopping]: The two rules also give different tree shapes, which matters if you are going to estimate rates from the trees. `n_extant` stops the first moment that many lineages are alive together, then draws one more waiting time and puts the present where the next event would have fired, so the two newest tips have a real branch length instead of a zero-length one. With `death=0` that is exactly the standard way to sample a tree of a given size [@hartmann2010sampling]. With extinction it is not: the run stops the *first* time it touches the target, so the trees come out shallower than a birth–death process conditioned on that many tips. The bias in tree height: nothing measurable at `death=0`, about a tenth at 10 tips with `death` at 0.4 of `birth`, a third to a half at 10 tips with `death` at 0.8 of `birth`, and back to nothing by 50 tips at moderate turnover. The same rule also touches measured rates: the run ends where the next event would have fired, without applying it, so the last stretch of every branch contains no event by construction, and a rate measured back from an `n_extant` run comes out slightly low, about one percent at 150 tips and less on larger trees. Use `total_time` if you need the conditioned distribution exactly or you are checking realised rates against their settings; otherwise say in your methods which rule made the trees.

## Non-constant rates

A birth or death rate can change along the tree. It can depend on **time**, on **the number of lineages that are alive at any given time**, or on a lineage's **ancestry**, inherited from its parent or not. Each of these possibilities is an ordinary **modifier** that affects the rates (see Appendix A or the gallery for code examples). The table in 2.4 names each one as the field does, beside the example that shows it.

- **On time.** The rates change at set points in time: the skyline, or episodic, tree. It is written as a schedule of factors on the base: `birth=PerLineage(1.0).changing_at({0: 1.0, 3: 0.33})` runs at the full rate until time 3 and a third of it after. Each entry holds from its own time up to the next, and the earliest entry also applies *backwards* to the origin, so a schedule that starts at time 3 runs at that factor for the whole tree, not only after time 3. Start it at 0 whenever you mean "full rate until".
- **On total diversity.** The rate slows as the tree fills up, so diversity levels off at a carrying capacity instead of growing without bound: `birth=PerLineage(1.0).scaled_by(TotalDiversity(cap=200))`.
- **On the parent's rate.** Each lineage inherits its parent's rate, nudged at every split, so rates wander across the tree and close relatives resemble each other: `birth=PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.3)))`. The distribution inside `Drift` is the step taken at each split; the steps accumulate, so lineages deeper in the tree spread further apart.
- **By lineage.** Each lineage draws its own rate instead, with no memory of its parent, the same spelling without `Drift`: `birth=PerLineage(1.0).varying_among('lineages', LogNormal(0.0, 0.5))`. Bare, the distribution is the spread of rates across lineages, and nothing accumulates.

## Other models

One model does not fit the modifier framework: a **mass extinction**, where at one instant each standing lineage dies with a probability you set. `mass_extinctions=[(3, 0.75)]` kills three quarters of the living lineages at time 3; diversity crashes at the pulse and recovers after it. It needs a run with a fixed end, so it goes with `total_time` and is refused beside `n_extant`.

## Literature

| What it does | From the literature | Gallery |
|-------------------|--------------------------------|---|
| rates change at set times | skyline / episodic birth–death [@stadler2011mammalian] | [Sp3](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:rateshift--> |
| rate slows as the tree fills | diversity-dependent diversification [@rabosky2008densitydependent; @etienne2012diversitydependence] | [Sp4](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:diversity--> |
| rates drift, inherited at each split | ClaDS [@maliet2019clads] | [Sp5](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:inherited--> |
| rates vary, drawn afresh per lineage | uncorrelated ("relaxed") rates | [Sp6](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:perlineage--> |
| a fraction culled at an instant | mass extinction | [Sp7](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:massext--> |

## Sampling

Sampling and fossils decide not how the tree grows, but how much of it you get to see.

By default you see every surviving species. **`sampling`** keeps a random fraction of the extant tips, so `sampling=0.5` gives you half [@stadler2009incomplete].

`n_extant` counts survivors, and sampling happens afterwards, so the two compose rather than cancel: `n_extant=20, sampling=0.5` grows to 20 survivors and then shows you about 10 of them. If you want 20 tips, ask for 40. The rest are not gone: they stay in the complete tree with the fate `unsampled`, which is why the run's summary counts them separately from the extinct.

## Fossils

**`fossils`** recovers lineages from the past [@heath2014fossilized; @gavryushkina2014sampledancestor]. Fossils are sampled along **every** branch of the complete tree at a rate you set, so `fossils=0.1` scatters dated samples through its history. They are a side output, written to `species_fossils.tsv` beside the trees (Appendix B); a fossil does not remove its lineage and does not appear in the extant tree ([Sp8](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:sampling-->).

## On the command line

The flags are the Python keywords: the base rates, the stop condition, and the sampling and fossil settings. Bare, `zombi2 species out/` runs the defaults the quickstart used: `birth` 1.0, stopping at 20 extant tips, every tip sampled. A varying rate goes on its flag in the written form, the same text Python takes.

```bash
# a birth–death tree of 20 surviving lineages
zombi2 species out/ --birth 1.0 --death 0.3 --n-extant 20 --seed 1

# grow to time 5; at time 3 each lineage dies with probability 0.75; half the survivors sampled
zombi2 species out/ --birth 1.0 --death 0.4 --total-time 5 \
    --mass-extinction 3 0.75 --sampling 0.5 --seed 2

# a skyline: full birth rate until time 3, a third of it after
zombi2 species out/ --total-time 5 --seed 3 \
    --birth "PerLineage(1.0).changing_at({0: 1.0, 3: 0.33})"
```
