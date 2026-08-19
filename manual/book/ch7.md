# Traits

The trait level evolves **phenotypes**: a body size, a habitat, the presence or absence of a structure. A trait evolves along the species tree like everything else here: the **complete** tree, extinct lineages included. Hand it the species run itself or its `.complete_tree`; every level takes either and runs on the complete tree. There are two kinds, continuous and discrete.

![The two kinds, on the same tree. A **continuous** trait starts at one value and drifts along every branch, so every node ends at a number of its own, with each branch painted with the value at its far end, and close relatives end up close. A **discrete** trait sits in one state and switches to another now and then; the dots mark the two switches, and every branch below a switch is in the new state until it switches again. The chips at the tips are the state each one ended in. The two are different functions, `simulate_continuous` and `simulate_discrete`, because they answer different questions.](figures/trait_kinds_print.png){width=100%}

## Continuous traits

A continuous trait **diffuses**. Give it a starting value and a rate, and it wanders down every branch, its variance growing in proportion to elapsed time. On its own that is **Brownian motion**, the null model of continuous trait evolution:

```python
from zombi2 import species, traits
from zombi2.params import PerLineage

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
```

```python
# BM: a body size diffusing from 0 at variance-rate σ² = 1.0
size = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)

size.values                # each extant tip's value, keyed by tip name: {'n19': -0.39, ...}
size.node_values           # every node of the complete tree, ancestors and extinct included
size.write("out/traits/")  # trait_values.tsv · trait_tree.nwk · trait_summary.json (Appendix B)
```

A discrete run gives the same shape back, states instead of numbers, and writes `trait_events.tsv` besides: the switch log, which is the file a conditioned run reads as its driver (Chapter 8).

Here `rate` is the Brownian variance-rate σ², the trait level's reading of "how fast". It is counted per lineage, the only scope here, so `PerLineage(1.0)` is the bare `1.0` in wrapper form, and like every rate in ZOMBI2 it takes the verbs. Two variations, the first by two more arguments, the second by a verb on the rate:

```python
# OU: the same diffusion, pulled toward an optimum value
traits.simulate_continuous(tree, start=0.0, rate=1.0,
                           reverts_to=2.0, pull=0.5, seed=1)

# early burst: the diffusion rate itself decays through time
traits.simulate_continuous(tree, start=0.0,
                           rate=PerLineage(1.0).changing_at({0: 1.0, 5: 0.2}), seed=1)
```

The **Ornstein–Uhlenbeck** process is Brownian motion with a rubber band: `reverts_to` is the optimum it is pulled back toward, and `pull` is how hard. **Early burst** (or ACDC) is a diffusion rate that decays as the tree ages, so most of the divergence happens near the root ([Tr3](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:early_burst-->); it is written with the same `changing_at` that gives the species tree its skyline, a schedule of steps, so a finer schedule is a finer approximation of the decay.

Two more arguments sit alongside `rate`. `regimes=` paints a multi-optimum OU, where clades pull toward different optima ([Tr6](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:regimes-->): a discrete run supplies the painting, and `reverts_to` becomes a dict keyed by its states, as in `regimes=habitat, reverts_to={"cave": -1.0, "surface": 3.0}`. And `at_speciation=` adds a jump *at* each split rather than along the branches, so change concentrates at branching ([Tr7](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:jumps-->). The value is the jump's variance: `at_speciation=0.5` draws each jump from a normal with variance 0.5, standard deviation about 0.7. None of these is a separate model with its own function and its own parameters, which is why they combine: a trait that bursts early *and* reverts to an optimum is one rate with one verb and two arguments.

`regimes=` is the one argument that asks you to give things up, and it refuses loudly rather than ignoring what you passed: with it the σ² is a plain number (not a modified rate), the jump variance is one number shared across regimes, and the run is one trait (no `correlation=`).

## Discrete traits

A discrete trait takes a finite set of states and switches between them along the branches: a continuous-time Markov chain, the field's **Mk model** [@lewis2001mk].

```python
# Mk: habitat flips between two states at rate 0.1
traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                         switch=0.1, start="marine", seed=1)
```

Leave `start` out and the root state is drawn uniformly from `states`. When the flips are not symmetric, replace the single rate with directed rates, one per ordered pair: the dict below, or a full `k × k` matrix whose entry at row *i*, column *j* is the rate *i* → *j* ([Tr9](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:asymmetric-->):

```python
# asymmetric: gains are commoner than losses
traits.simulate_discrete(tree, states=["absent", "present"],
                         switch={"absent->present": 0.2, "present->absent": 0.05},
                         seed=1)
```

Two more discrete models live on this same function. The **threshold** model reads the states off a continuous liability: `simulate_discrete(tree, states=["absent", "present"], liability=1.0, threshold=0.0)` diffuses a Brownian liability from `start` (a number here, 0.0 by default) and cuts it at the thresholds, `k − 1` increasing cuts for `k` states ([Tr10](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:threshold-->; `--liability` and `--threshold` on the command line). Give `liability` as a dict and add a `correlation={("wings", "eyes"): 0.6}` overlay and two discrete traits evolve **together**, their liabilities diffusing jointly and each cut by the shared thresholds, which is the correlated-binary model of the Literature table ([Tr13](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:dependent-->).

## Correlated traits

Two traits that evolve independently are two separate calls, in either order. Two traits that drift *together* cannot be simulated one before the other, because each is entangled with the other as it unfolds ([Tr11](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:correlated-->). Correlation is specified as per-trait rates plus a correlation overlay:

```python
traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8},        # one variance-rate per trait
    correlation={("size", "limb"): 0.6},    # the overlay, ∈ [−1, 1]
    seed=1)
```

The overlay is a dimensionless number in `[−1, 1]`, not a covariance matrix. Under `correlation=` the per-trait rates are plain numbers.

The overlay carries the reversion and the speciation jumps as well. Add `reverts_to` and `pull` and each trait reverts to its own optimum at its own strength; add `at_speciation` and the jump at each split is drawn under the same correlation the diffusion uses. Each argument takes one value shared across the traits, or a dict giving one per trait:

```python
traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8},
    correlation={("size", "limb"): 0.6},
    reverts_to={"size": 3.0, "limb": -1.0},   # each trait its own optimum
    pull={"size": 1.5, "limb": 0.4},          # and its own strength
    at_speciation=0.5, seed=1)
```

## Literature

Trait models arrive under a thicket of names, and a reader who wants "an OU model" or "a threshold model" should be able to find it. The names live here, in one table, each beside the example that shows it; the example carries the run that made it, so the table does not spell the call out a second time. It organises nothing else in the chapter.

| What it does | From the literature | Gallery |
|-------------------|--------------------------------|---|
| a value diffusing | Brownian motion (BM) [@felsenstein1985comparative] | [Tr1](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:bm--> |
| diffusion pulled to an optimum | Ornstein–Uhlenbeck (OU) [@hansen1997stabilizing; @butler2004phylogenetic] | [Tr2](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:ou--> |
| diffusion rate decays through time | Early burst (EB / ACDC) [@harmon2010earlyburst] | [Tr3](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:early_burst--> |
| diffusion rate drifts between lineages | Variable-rates BM [@maliet2019clads] | [Tr4](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:varying_rate--> |
| diffusion rate slows as the clade fills | Diversity-dependent / ecological limits [@etienne2012diversitydependence] | [Tr5](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:diversity_dependent--> |
| the optimum differs between painted clades | Multi-optimum OU (OUM) [@beaulieu2012ouwie] | [Tr6](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:regimes--> |
| the value jumps at each split | Cladogenetic / punctuational change | [Tr7](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:jumps--> |
| a discrete state switching | Mk (k-state Markov) | [Tr8](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:discrete--> |
| gains and losses at different rates | All-rates-different Mk (ARD) | [Tr9](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:asymmetric--> |
| discrete driven by continuous liability | Threshold / liability (Wright–Felsenstein) [@felsenstein2012threshold] | [Tr10](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:threshold--> |
| traits evolving together | Multivariate BM | [Tr11](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:correlated--> |
| traits reverting together, each to its own optimum | Multivariate OU, diagonal drift [@clavel2015mvmorph] | [Tr12](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:mv_ou--> |
| discrete traits evolving together | Correlated binary / Pagel [@pagel1994correlated] | [Tr13](https://aadavin.github.io/zombi2/gallery.html#traits)<!--gallery:dependent--> |

## On the command line

`--kind` is required, `continuous` or `discrete`, and it decides which of the other flags apply, since `--rate` and the OU flags belong to a continuous trait while `--states` and `--switch` belong to a discrete one. Any default would silently pick a model for you, so there is none:

```bash
# a continuous (Brownian) trait along a species tree
zombi2 traits out/ --kind continuous \
    --start 0.0 --rate 1.0 --seed 1

# a discrete two-state trait
zombi2 traits out/ --kind discrete \
    --states marine,terrestrial --switch 0.1 --seed 1

# the states the next command drives from — a discrete run writes trait_events.tsv by default,
# the driver file a conditioned genome, sequence or trait run reads (Chapter 8)
zombi2 traits out/ --kind discrete \
    --states cave,surface --switch 0.1 --seed 1

# a second trait, kept apart in out/traits/diet/ by --name, whose switch rate reads the first (Chapter 8)
zombi2 traits out/ --kind discrete --name diet --states plant,fish --seed 2 \
    --switch "PerLineage(0.2).scaled_by('out/traits/trait_events.tsv', {'cave': 5.0, 'surface': 1.0})"
```

An unnamed trait run writes to `out/traits/` and replaces whatever an earlier unnamed run left there, so the file the last command reads belongs to the cave/surface run just above it. `--name` gives a run its own directory, `out/traits/diet/`, which is what lets two traits sit side by side.

Every rate flag takes a rate in its written form, `--switch` as much as `--rate`, so the expression above is the same text the Python API takes. `--switch` also takes the keyword's other two shapes: the `{'a->b': rate}` dict and the `k × k` matrix.

Two keywords have no flag. `correlation=` grows several traits in one call, and the command line grows one trait per run. `regimes=` takes a discrete result object, the painting, handed to the call. Both stay in the Python API.

The trait evolves on the **complete** tree, extinct lineages included, so `species_complete.nwk` is the file to hand it. An external tree goes in with `--from` (a Newick file, or another run's directory); if it is not ultrametric you must declare each tip's fate with `--tip-fates`, a TSV of `tip<TAB>extant|extinct|unsampled`, because ZOMBI2 will not guess which early-ending tips are extinct.
