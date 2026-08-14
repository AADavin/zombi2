# Species trees

The species tree is the backbone every other level runs on, so it is where almost every workflow begins.

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give a **speciation rate** and an **extinction rate**, and every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest.

The two rates set the tempo. Their difference fixes how fast diversity builds up. Their ratio fixes how much of the history is hidden, because a lineage that goes extinct takes its part of the tree with it. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process. The lineages that died are kept: they are in the complete tree, which is the tree the next level runs along — that is what lets a gene be transferred out of a lineage that later dies — while the extant tree is the one an analysis would be handed. Appendix B says which file holds each.

![A species tree grown by the birth–death process. Every lineage alive at a given moment has the same chance per unit time of splitting or of dying. The lineages that died are drawn dashed and stop where they died; the survivors reach the present. Both are in the complete tree, and only the solid ones are in the extant tree.](figures/species_tree.pdf){width=100%}

You also say when to stop: grow to a fixed **total time** (`total_time`), or to a fixed **number of surviving lineages** (`n_extant`).[^stopping]

`n_extant` bounds the run by construction; `total_time` does not, because standing diversity grows like exp((birth − death) · t), so a rate slightly too high or a time slightly too long is the difference between a thousand lineages and ten million. A run that passes **100,000 standing lineages** stops with an error rather than filling memory. Raise `max_lineages` if that is the size you want, or set it to `None` to lift the guard. It never truncates: a tree cut off at a size is no longer a sample from the process you asked for.

`total_time` is not conditioned on survival. A run that dies out raises an error rather than handing back a tree with no present, so if you loop over seeds and skip the failures, you are conditioning on survival yourself, and everything downstream inherits that.

[^stopping]: The two rules also give different tree shapes, which matters if you are going to estimate rates from the trees. `n_extant` stops the first moment that many lineages are alive together, then draws one more waiting time and puts the present where the next event would have fired, so the two newest tips have a real branch length instead of a zero-length one. With `death=0` that is exactly the standard way to sample a tree of a given size [@hartmann2010sampling]. With extinction it is not: the run stops the *first* time it touches the target, so the trees come out shallower than a birth–death process conditioned on that many tips: nothing measurable at `death=0`, about a tenth of the tree height at 10 tips with `death` at 0.4 of `birth`, a third to a half at 10 tips with `death` at 0.8 of `birth`, and back to nothing by 50 tips at moderate turnover. Use `total_time` if you need the conditioned distribution exactly; otherwise say in your methods which rule made the trees.

## Rates that vary

A birth or death rate need not be constant. It can depend on **time**, on **how crowded the tree is**, or on a lineage's **ancestry**. There are four of them, and the Literature table below names each one as the field does, beside the example that shows it.

- **On time.** The rates change at set points in time, which is the skyline, or episodic, tree: full rate until a breakpoint, a third of it after. Each entry in the schedule holds from its own time up to the next, and the earliest entry also applies *backwards* to the origin — so a schedule that starts at time 3 runs at that rate for the whole tree, not only after time 3. Start it at 0 whenever you mean "full rate until".
- **On total diversity.** The rate slows as the tree fills up, so diversity levels off at a carrying capacity instead of growing without bound.
- **On the parent's rate.** Each lineage inherits its parent's rate, nudged at every split, so rates wander across the tree and close relatives resemble each other.
- **By lineage.** Each lineage draws its own rate instead, with no memory of its parent. The distribution means a different thing in the two: inherited, it is the step taken at each split, which accumulates, so lineages deeper in the tree spread further apart; drawn afresh, it is the spread of rates across lineages and nothing accumulates.

## Other models

One model does not fit the modifier framework: a **mass extinction**, where at one instant only a fraction of the living lineages survive. Diversity crashes at the pulse and recovers after it.

## Literature

| What it does | Gallery | From the literature |
|---|---|---|
| rates change at set times | [Sp3](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:rateshift--> | skyline / episodic birth–death [@stadler2011mammalian] |
| rate slows as the tree fills | [Sp4](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:diversity--> | diversity-dependent diversification [@rabosky2008densitydependent; @etienne2012diversitydependence] |
| rates drift, inherited at each split | [Sp5](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:inherited--> | ClaDS [@maliet2019clads] |
| rates vary, drawn afresh per lineage | [Sp6](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:perlineage--> | uncorrelated ("relaxed") rates |
| a fraction culled at an instant | [Sp7](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:massext--> | mass extinction |

## Sampling

Two more choices decide not how the tree grows, but how much of it you get to see.

By default you see every surviving species. **`sampling`** keeps a random fraction of the extant tips, so `sampling=0.5` gives you half [@stadler2009incomplete]. It thins a tree that has already grown, so it costs nothing.

`n_extant` counts survivors, and sampling happens afterwards, so the two compose rather than cancel: `n_extant=20, sampling=0.5` grows to 20 survivors and then shows you about 10 of them. If you want 20 tips in hand, ask for 40. The rest are not gone: they stay in the complete tree with the fate `unsampled`, which is why the run's summary counts them separately from the extinct.

**`fossils`** does the opposite: it recovers lineages from the past [@heath2014fossilized; @gavryushkina2014sampledancestor]. Fossils are picked up along **every** branch of the complete tree at a rate you set, a surviving lineage's branch as readily as an extinct one, so `fossils=0.1` scatters dated samples through its history. They are a side output, reported alongside the trees; a fossil does not remove its lineage and does not appear in the extant tree ([Sp8](https://aadavin.github.io/zombi2/gallery.html#species)<!--gallery:sampling-->).

## On the command line

The command mirrors the Python call: the base rates, the stop condition, and the sampling and fossil settings each have a flag.

```bash
# a birth–death tree of 20 surviving lineages
zombi2 species out/ --birth 1.0 --death 0.3 --n-extant 20 --seed 1

# grow to time 5, with a mass extinction at time 3 and half the survivors sampled
zombi2 species out/ --birth 1.0 --death 0.4 --total-time 5 \
    --mass-extinction 3 0.75 --sampling 0.5 --seed 2
```
