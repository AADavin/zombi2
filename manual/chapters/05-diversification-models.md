# Diversification models

ZOMBI2's default species-tree sampler runs **backward in time**: it draws internal-node ages
i.i.d. from the reconstructed birth–death CDF and assembles a ranked ultrametric tree by uniform
coalescence. That machinery keeps the constant-rate `BirthDeath` and `Yule` models cheap, and it
extends cleanly to piecewise-constant rates. But a whole family of macroevolutionary models
cannot be written as "draw i.i.d. node ages": rates that depend on the standing diversity, rates
that drift lineage-by-lineage, scheduled shifts in one clade, instantaneous kills, and serially
sampled fossils all need the process to be simulated **forward in time**, retaining extinct
lineages. This chapter covers the diversification models ZOMBI2 ships today, flagging which run
backward and which are forward-only.

Every forward model is invoked through the same entry point, differing only in the model object
passed and the conditioning:

```python
import zombi2 as z

tree = z.simulate_species_tree(
    model,               # a diversification model object
    age=5.0,             # crown age (a fixed present); or n_tips=…
    direction="forward", # forward growth, extinct lineages retained
    seed=1,
)
```

A forward run returns the **complete** tree (extinct leaves carry `is_extant=False` at their death
times). `z.prune(tree)` extracts the survivors-only reconstructed counterpart. Any complete tree
feeds `simulate_genomes` unchanged, and extinct lineages become ghost transfer partners
automatically.

::: note
The forward growth loop is exact: because every model below has rates that are constant *between*
events, the next event time is drawn from the summed rate by an exact-Gillespie step — no thinning
bound needed. Mass-extinction pulses and incomplete sampling layer on top of any of them.
:::

## Episodic (skyline) birth–death

`EpisodicBirthDeath` lets speciation and extinction be **piecewise-constant through time** — the
model behind shifting diversification regimes and *gradual* extinction pulses
[@stadler2013birthdeathskyline; @hohna2016tess]. Rates are given one value per epoch, ordered from
the present backward, with the epoch boundaries supplied as strictly increasing **ages** before
the present:

```python
# a gradual extinction spike: normal extinction recently, high extinction older than age 1
epi = z.EpisodicBirthDeath(birth=[1.0, 1.0], death=[0.2, 3.0], shifts=[1.0])
tree = z.simulate_species_tree(epi, n_tips=30, age=4.0, seed=1)
```

`birth[i]` and `death[i]` apply to epoch `i`; `shifts` has one fewer entry than `birth` (the
boundaries between epochs). A single epoch (`shifts=[]`) reproduces the constant-rate `BirthDeath`.
Because the reconstructed tree is still a coalescent point process, ZOMBI2 samples i.i.d. node ages
from the numerically inverted piecewise CDF and assembles exactly as for the constant-rate case —
the tree stays ultrametric. This is the one time-varying model that still runs through the
**backward** sampler; it also runs forward in age mode.

![Piecewise-constant speciation and extinction rates across skyline epochs, with boundaries at fixed ages before the present.](figures/model_episodic.pdf)

### Incomplete extant sampling

Real phylogenies rarely include every living species. Pass `sampling_fraction`, the probability
$\rho$ that an extant species is sampled:

```python
z.EpisodicBirthDeath(birth=[1.0], death=[0.3], shifts=[], sampling_fraction=0.25)
```

In a forward run, unsampled extant tips are marked `is_extant=False`. Incomplete sampling keeps the
tree ultrametric, and is accepted as an overlay by every forward diversification model in this
chapter.

![Incomplete extant sampling: only a fraction of the living tips (filled) are retained; the rest (open) are dropped.](figures/model_sampling.pdf)

::: warning
Serial sampling *through time* — dated fossil tips, as in the fossilized birth–death process below
— is a different mechanism from incomplete extant sampling. It requires forward simulation with
retained extinct lineages, and is not available through the backward sampler.
:::

## Mass extinctions

Raising `death` over an epoch spreads extra extinction *smoothly* through a time window. A **mass
extinction** in the palaeobiological sense is instead an *instantaneous, tree-wide pulse*: at a
single instant a large fraction of the standing diversity is wiped out at once. Give any forward
model a `mass_extinctions` list of `(age, fraction)` pulses:

```python
# a radiation punctuated by two cataclysms (75% then 50% die), grown forward:
m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(1.0, 0.75), (2.5, 0.5)])
tree = z.simulate_species_tree(m, age=5.0, direction="forward", seed=1)

# pulses compose with an episodic background, too:
m = z.EpisodicBirthDeath(birth=[1.0, 1.4], death=[0.2, 0.3], shifts=[2.0],
                         mass_extinctions=[(1.0, 0.8)])
```

At each `age` before the present, every lineage then alive **independently** dies with probability
`fraction` (survives with `1 - fraction`) — the standard survival-pulse formulation of TreeSim and
TESS [@stadler2011mammalian; @hohna2016tess]. Killing is per-lineage Bernoulli, so the realized
fraction fluctuates around `fraction`. The victims become ordinary extinct leaves at the pulse
instant, so `simulate_genomes` treats them as ghost transfer partners — a mass extinction leaves a
*genomic* signature (families lost with the dead clades, transfers *from* the dead) with no extra
wiring.

Mass extinctions are a **forward** feature and need **age mode** (their times are ages before a
fixed present); the backward reconstructed sampler and `n_tips` mode reject them. A pulse of
`fraction=1.0` wipes the whole tree out, and the requirement of at least two surviving lineages
then rejects the run.

```bash
zombi2 species --mode forward --age 5 \
  --mass-extinction 1.0 0.75 --mass-extinction 2.5 0.5 -o out/
```

![A mass extinction: at one instant a fraction of the standing diversity is pruned, driving a sharp drop in the lineages-through-time curve.](figures/mass_extinction.pdf)

## Per-lineage rates: ClaDS

In the models so far, every lineage shares the same rates at any instant. **ClaDS**
[@maliet2019clads] instead gives each lineage its *own* speciation rate: at each speciation the two
daughters inherit the parent's rate times an independent lognormal jump,
$$\lambda_{\text{child}} = \lambda_{\text{parent}} \cdot \exp\!\big(\mathcal{N}(\log\alpha,\ \sigma^2)\big),$$
so rates drift lineage-by-lineage down the tree. It is the diversification counterpart of a relaxed
molecular clock, and captures the heavy among-clade rate variation real phylogenies show.

```python
# alpha<1 = speciation slows toward the present; sigma = jump spread; turnover = mu/lambda
m = z.ClaDS(lambda_0=1.0, alpha=0.9, sigma=0.2, turnover=0.1)
tree = z.simulate_species_tree(m, age=5.0, direction="forward", seed=1)   # or n_tips=…
```

`alpha` is the trend (`alpha<1` reproduces the empirically typical slow-down of speciation toward
the present), `sigma` the jump spread, and extinction is set by a constant turnover
$\varepsilon = \mu/\lambda$. `turnover=0` is ClaDS0 (pure birth with shifts); `turnover>0` adds
proportional extinction. ClaDS is **forward-only** — per-lineage rates have no closed-form
reconstructed CDF — and runs in either `age` or `n_tips` mode. It accepts `sampling_fraction` and
`mass_extinctions` overlays like `BirthDeath`.

```bash
zombi2 species --mode forward --diversification clads --birth 1.0 --age 5 -o out/
```

![A ClaDS tree, branches shaded by their per-lineage speciation rate: some clades radiate, others stall.](figures/clads.pdf)

## Diversity-dependent diversification

`DiversityDependent` [@rabosky2008densitydependent; @etienne2012diversitydependence] makes
speciation slow as the tree fills an ecological carrying capacity $K$,
$$\lambda(n) = \max\!\big(0,\ \lambda_0\,(1 - n/K)\big),$$
with constant extinction $\mu$. The tree radiates fast when small and saturates near $K$ — a
diversity brake, and the macroevolutionary analogue of the per-family `carrying_capacity` ZOMBI2
already offers for genes.

```python
m = z.DiversityDependent(lambda_0=2.0, death=0.2, carrying_capacity=50)
tree = z.simulate_species_tree(m, age=15.0, direction="forward", seed=1)   # or n_tips <= K
```

With $\mu=0$ the tree saturates at exactly $K$; with $\mu>0$ it settles near the equilibrium
$n^{*} = K\,(1 - \mu/\lambda_0)$. The model is **forward-only** but supports both `age` and `n_tips`
mode — in `n_tips` mode the target must be `<= K`. It accepts `sampling_fraction` and
`mass_extinctions` overlays.

```bash
zombi2 species --mode forward --diversification diversity-dependent \
  --birth 2 --death 0.2 -K 50 --age 15 -o out/
```

![Diversity-dependent diversification: a fast early radiation flattening into a plateau as the lineage count approaches the carrying capacity K.](figures/diversity_dependent.pdf)

## Clade-specific rate shifts

Where ClaDS shifts *every* lineage a little at *every* speciation, `CladeShiftBirthDeath` shifts
*one* clade a lot at a *scheduled* time — the discrete, hand-specified version of rate
heterogeneity. The tree runs at the background `(birth, death)` until, at each scheduled age before
the present, a uniformly chosen lineage then alive (and all its descendants) adopts a new
`(birth, death)` regime. It is the model for "a key innovation sparks a radiation in one clade."

```python
# a slow background; at age 3 one clade starts diversifying fast
m = z.CladeShiftBirthDeath(0.6, 0.4, clade_shifts=[(3.0, 2.0, 0.1)])
tree = z.simulate_species_tree(m, age=5.0, direction="forward", seed=1)
```

The shifted lineage is drawn at random: you cannot name an unborn clade in a forward run, and
contemporaneous lineages are exchangeable. Supply several `(age, birth, death)` shifts for several
clades. This model is **forward-only** and **age mode only** (the shifts are scheduled as ages
before a fixed present); it accepts `sampling_fraction` and `mass_extinctions` overlays.

```bash
zombi2 species --mode forward --age 5 --birth 0.6 --death 0.4 --clade-shift 3.0 2.0 0.1 -o out/
```

![A clade-specific shift: at a scheduled age one lineage and all its descendants adopt a new speciation/extinction regime, sparking a fast-diversifying clade.](figures/clade_shift.pdf)

## Fossilized birth–death

The same `BirthDeath` (or `EpisodicBirthDeath`) model gains **serial sampling through time** via
optional keywords [@gavryushkina2014sampledancestor; @heath2014fossilized]. Beyond speciation
($\lambda$) and extinction ($\mu$), lineages are sampled through time at rate $\psi$
(`fossilization`) — each a **dated fossil tip** — and extant lineages are sampled at the present
with probability $\rho$ (`sampling_fraction`). Because fossils are sampled points on extinct or
surviving lineages, this is a **forward** feature that requires retaining the extinct lineages.

```python
m = z.BirthDeath(birth=1.0, death=0.5, fossilization=0.5, sampling_fraction=0.9)
tree = z.simulate_species_tree(m, age=6.0, direction="forward", seed=1)  # complete tree + fossils
```

Fossil tips carry `sampled=True, is_extant=False` at their (past) sampling times; sampled extant
tips carry `sampled=True, is_extant=True`. With `removal` set below 1 (`--removal r`), a sampled
lineage *continues* with probability `1 - r` instead of being removed — a **sampled ancestor**,
represented as a degree-two node that the gene simulator passes genomes straight through, so DTL
simulation runs unchanged on such trees.

```bash
zombi2 species --mode forward --age 6 --fossilization 0.5 --sampling-fraction 0.9 --removal 0.5 -o out/
```

![Fossilized birth–death: dated fossil tips (sampled through time along extinct and surviving lineages, as diamonds) alongside sampled extant tips.](figures/model_fbd.pdf)

::: warning
The occurrence birth–death process [@andreoletti2022occurrence], which mixes fossil, occurrence,
and extant data types, builds on this same forward machinery but is not yet shipped.
:::

::: note
Every time-varying diversification model shares an identifiability caveat: distinct histories of
$\lambda(t)$ and $\mu(t)$ can produce the same distribution of reconstructed (extant-only) trees
[@louca2020identifiability]. Serial sampling — fossils — is one way to break that degeneracy, since
it constrains the extinct part of the process directly.
:::
