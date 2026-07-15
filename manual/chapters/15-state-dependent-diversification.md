# State-dependent diversification

The most studied coevolution is between the species tree and a **trait**: a lineage's character
sets how fast it speciates and goes extinct. This is the species–traits pair. When the trait drives
the tree (`traits:species`), the arrow points into S, so the tree becomes an **output** grown jointly
with the trait — the SSE family. The reverse edge (`species:traits`) instead lets speciation kick the
trait, an overlay on a given tree. Switch both on and you get ClaSSE.

## `traits:species` — the SSE family

A discrete or continuous trait drives speciation and extinction, and the tree is grown *jointly* with
the trait. Because the trait shapes the topology, this edge **produces** the tree (it takes no `-t`)
and a stopping condition (`n_tips` or `age`) instead. The driver is `simulate_sse`:

```python
from zombi2.coevolve import simulate_sse, BiSSE

# BiSSE: state 1 speciates 3x faster, so it dominates the standing tips
res = simulate_sse(
    BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2, q01=0.1, q10=0.1),
    n_tips=200, seed=1)
res.tree               # complete tree; prune() for the reconstructed one
res.labeled_values()   # the trait at the extant tips
```

![State-dependent diversification (BiSSE). **A**, the model: two states, the anagenetic transitions (arrow width = rate), and each state's speciation rate drawn as a fork (width = $\lambda$) — state 1 branches three times faster. **B**, one `simulate_sse` realization: branches are heavy where the lineage is in state 1 and light in state 0, lineages that go extinct end in a cross, and the extant tips carry chips. The fast (state-1) lineages proliferate and fill most of the standing tips — the diversification signal, written into the shape of the tree itself.](figures/sse.pdf){width=100%}

`BiSSE` is the **binary** state-dependent birth–death process [@maddison2007bisse]; `MuSSE` the
**k-state** variant [@fitzjohn2012diversitree]; `QuaSSE` the **continuous-trait** variant, whose
speciation and extinction are functions of the trait value [@fitzjohn2010quasse]; and `HiSSE` the
**hidden-state** model [@beaulieu2016hisse], where diversification is driven by an unobserved class
rather than the focal trait — the honest null against which a real association must be judged.

```python
import numpy as np
from zombi2.coevolve import MuSSE, QuaSSE, HiSSE

# k-state:
MuSSE(birth=[1, 3], death=[0.2, 0.2], Q=np.array([[-0.1, 0.1], [0.1, -0.1]]))
# continuous trait:
QuaSSE(speciation=lambda x: 1 + 2 / (1 + np.exp(-x)),
       extinction=lambda x: 0.2, sigma2=0.5, rate_bound=5.0, x0=0.0)
# hidden classes:
HiSSE(classes=[BiSSE(0.5, 0.7, 0.2, 0.2, 0.1, 0.1),
               BiSSE(2.0, 3.0, 0.2, 0.2, 0.1, 0.1)], hidden_transition=0.1)
```

From the command line, `--sse-model` picks the flavour:

```bash
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/
```

It writes `species_tree.nwk` (the tree the trait's rates shaped), `traits.tsv` (every node — tips
*and* ancestral states) and `trait_tree.nwk`. `--sse-model musse` takes `--birth`/`--death` vectors
and a `--q-matrix` file; `--sse-model quasse` takes a sigmoidal speciation (`--spec-low/high/center/slope`)
and Brownian `--diffusion`; and `--sse-model hisse` adds `--hidden-classes` diversification regimes
spanning the base rates up to `--hidden-scale`× faster (`--hidden-switch` between them). **What it
recovers:** the fast-speciating state accumulates lineages, so it dominates the standing tips — the
diversification signal is written into the tree shape itself.

![The continuous variant, QuaSSE. **A**, the model: the speciation rate is a rising function $\lambda(x)$ of the trait while extinction is flat, and the trait itself diffuses by Brownian motion (the axis is tinted with the viridis ramp used to paint the tree). **B**, one realization, each branch painted by its trait value: the high-value (yellow) lineages branch faster and proliferate, while the low-value (blue) lineages stay sparse and go extinct — the same "fast state fills the tips" signal as BiSSE, now on a continuous character.](figures/sse_quasse.pdf){width=100%}

## `species:traits` — cladogenetic trait jumps

The reverse arrow makes the trait jump *at* each speciation rather than (or in addition to) drifting
along the branches — speciational, or *cladogenetic*, evolution: sister species differ because
something happened at their split. On its own this edge has no arrow into S, so it runs on a **given**
tree, with a `Cladogenesis` kernel layered on an ordinary anagenetic trait model. Setting the
anagenetic rates to zero gives *purely* speciational change:

```python
import numpy as np
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2 import Cladogenesis, simulate_traits, Mk

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=60, age=6.0, seed=1)

# a purely speciational binary trait: no within-branch change (zero-rate
# Mk), a jump at each split
res = simulate_traits(tree, Mk(np.zeros((2, 2))),
                      cladogenesis=Cladogenesis(shift=0.4), seed=2)
```

`Cladogenesis(shift=…)` is the per-daughter state-hop probability for a discrete trait;
`Cladogenesis(jump_sigma2=…)` is the Gaussian jump variance for a continuous one. On the command line:

```bash
zombi2 coevolve --couple species:traits -t species_tree.nwk \
    --sse-model bisse --q01 0 --q10 0 --clado-shift 0.4 --seed 2 -o out/
```

**What it recovers:** change concentrated at the nodes — closely related tips can differ sharply while
long unbranched stretches stay constant, the signature a purely-gradual model cannot produce.

![Where trait change happens: anagenetic (as in BiSSE) vs cladogenetic (the ClaSSE addition), drawn on **one shared tree** so that only the *location* of change differs. **A**, anagenetic — the trait changes *along* the branches (open circles); at a speciation the daughters inherit the parent, so the amount of change scales with elapsed time and sister tips are usually alike. **B**, cladogenetic — the trait changes *at* the splits (filled diamonds), each daughter drawn as part of the speciation event; change scales with the number of speciations, so sister tips can differ sharply while long unbranched lineages stay constant. Same Gillespie, same tree — only the consequence of a speciation event differs.](figures/sse_cladogenetic.pdf){width=100%}

## Both arrows: ClaSSE

Turn on **both** `traits:species` and `species:traits` and you get the full ClaSSE feedback: the trait
shapes the tree *and* is kicked by its own branching. Because one arrow points into S, the tree is
again an output:

```python
res = simulate_sse(BiSSE(1, 3, 0.05, 0.05, 0.05, 0.05),
                   n_tips=200, cladogenesis=Cladogenesis(shift=0.3), seed=3)
```

```bash
zombi2 coevolve --couple traits:species --couple species:traits \
    --lambda0 1 --lambda1 3 --q01 0.05 --q10 0.05 --clado-shift 0.3 \
    --tips 200 --seed 3 -o out/
```
