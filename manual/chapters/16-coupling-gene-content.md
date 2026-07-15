# Coupling gene content

Gene content can be coupled to the tree that carries it — genes that drive diversification, genes
reshuffled at each speciation — and to a **trait** — genes whose retention tracks a phenotype, a gene
that unlocks a phenotypic optimum. This chapter takes the species–genes and traits–genes pairs in
turn, and closes with the null models that let you tell a real coupling from the tree's own
heterogeneity.

The overlay examples below (every edge except the two that grow the tree) run on a species tree you
supply; they build one with

```python
from zombi2.species import BirthDeath, simulate_species_tree
tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=60, age=6.0, seed=1)
```

## Species and genes

### `genomes:species` — key-innovation diversification

Gene content can drive the tree: a small panel of binary **driver** ("key innovation") families whose
presence sets each lineage's speciation and extinction rate. Drivers are gained de novo (origination)
and — the interesting part — by **transfer**, which is *frequency-dependent*: a driver carried by more
of the live population is donated more often, so it spreads as the tree grows. That feedback is why the
tree and the gene content must grow together; this edge produces the tree. The neutral bulk genome does
not affect diversification, so it is overlaid afterward on the finished tree with the ordinary
`genomes` machinery (exact under independent families, not an approximation).

```python
from zombi2.coevolve import simulate_gene_diversification, GeneDiversification

res = simulate_gene_diversification(
    GeneDiversification(n_drivers=2, lambda0=1.0, mu0=0.2, driver_speciation=1.2,
                        transfer=0.8, loss=0.3, root_drivers=1),
    n_tips=200, seed=1)
res.tree               # the tree the drivers shaped
res.driver_names()     # ['D0', 'D1']
res.tip_prevalence()   # fraction of extant tips carrying each driver
```

From the command line, the neutral genome overlay is a second, ordinary `genomes` call on the tree the
first command wrote:

```bash
zombi2 coevolve --couple genomes:species --drivers 2 --root-drivers 1 \
    --lambda0 1 --mu0 0.2 --driver-speciation 1.2 \
    --driver-transfer 0.8 --driver-loss 0.3 \
    --tips 200 --seed 1 -o out/

# overlay the neutral genome on the resulting tree (factorisation made explicit):
zombi2 genomes -t out/species_tree.nwk --trans 1 --loss 0.5 \
    --write profiles trees -o out/
```

`coevolve` writes `species_tree.nwk`, `drivers.tsv` (per-node driver presence) and
`drivers_manifest.tsv` (the effect sizes and rates). `--driver-speciation`/`--driver-extinction` are
the per-driver log-rate effects; `--driver-transfer`/`--driver-loss`/`--driver-origination` the driver
dynamics. **What it recovers:** a successful key innovation both spreads across the tips
(high `tip_prevalence`) and leaves the clades that carry it more speciose — a genomic cause of a
diversification rate shift, rather than a trait one.

![Key-innovation diversification. **A**, the model: a lineage carrying the driver (D+) speciates faster than one without it (D−, the fork width is the speciation rate); the driver is gained by origination and by frequency-dependent transfer, and lost. **B**, one realization: the driver *originates* on a single branch (the +), and from there the carrier (heavy) clade radiates into a speciose group while the non-carrier (light) lineages stay sparse. Here it reaches 69% of the extant tips — a genomic cause of a diversification-rate shift.](figures/key_innovation.pdf){width=100%}

### `species:genomes` — punctuational genome

The reverse of `genomes:species`, and the genomic twin of `species:traits`: here gene content does *not*
affect diversification, so it is an overlay on a **given** tree. A genome is evolved down the tree with
a **cladogenetic burst** of gene loss and gain at each speciation (a founder-effect upheaval), on top
of the usual gradual along-branch change:

```python
from zombi2.coevolve import simulate_cladogenetic_genome, CladogeneticGenome

# purely punctuational: gene content changes ONLY at speciations
# (no anagenetic loss/origination)
res = simulate_cladogenetic_genome(
    tree,
    CladogeneticGenome(initial_families=30, loss=0.0, origination=0.0,
                       cladogenetic_loss=0.15, cladogenetic_gain=3.0),
    seed=2)
res.genome_sizes()     # {node: family count}
res.profile_matrix()   # families × extant tips
```

```bash
zombi2 coevolve --couple species:genomes -t species_tree.nwk \
    --genome-size 30 --clado-gene-loss 0.15 --clado-gene-gain 3 --seed 2 -o out/
```

`--clado-gene-loss` is the per-family drop probability at each speciation and `--clado-gene-gain` the
mean number of families gained (Poisson); `--gene-loss`/`--gene-origination` add gradual along-branch
change (both `0` = pure punctuation). It writes `profiles.tsv`/`presence.tsv` and `genome_sizes.tsv`.
**What it recovers:** the model's signature is that **sister tips differ** — change is injected at
their split, not spread evenly along the branches.

![The genomic twin of the cladogenetic-trait figure: a genome evolved down **one shared tree** two ways. **A**, gradual — families are lost and gained *along* the branches (circles at the branch midpoints), so gene-content turnover scales with time and sister genomes stay similar. **B**, punctuational — gene content changes only in a *burst at each speciation* (diamonds at the nodes), so sister tips can differ sharply in size and content. The tip bars are the extant genome sizes. Same marker grammar as the trait figure: change on a branch vs change at a node.](figures/punctuational_genome.pdf){width=100%}

### Both arrows: co-diversification

Turn on **both** `genomes:species` and `species:genomes` and the *same* driver families both set the
diversification rates *and* are reshuffled by a cladogenetic burst at every speciation — the genomic
twin of ClaSSE. Because a burst can hand one daughter a key innovation and not its sister, speciation
*itself* seeds the rate heterogeneity that then plays out along the branches. One arrow points into S,
so the tree is again an output:

```python
from zombi2.coevolve import GeneDiversification, simulate_co_diversification

res = simulate_co_diversification(
    GeneDiversification(3, lambda0=1.0, mu0=0.15, driver_speciation=1.0,
                        loss=0.0, origination=0.0, transfer=0.0, root_drivers=1,
                        cladogenetic_loss=0.15, cladogenetic_gain=0.2),
    n_tips=200, seed=5)
res.tip_prevalence()   # drivers still spread: the genomes:species signal survives
```

```bash
zombi2 coevolve --couple genomes:species --couple species:genomes \
    --drivers 3 --lambda0 1 --mu0 0.15 --driver-speciation 1.0 \
    --driver-loss 0 --driver-origination 0 --driver-transfer 0 --root-drivers 1 \
    --driver-clado-loss 0.15 --driver-clado-gain 0.2 --tips 200 --seed 5 -o out/
```

`--driver-clado-loss`/`--driver-clado-gain` are the per-driver drop/gain probabilities of the burst;
with both `0` the `species:genomes` arrow is off and this is plain `genomes:species`. **What it recovers:**
the diversification signal (a key innovation over-represented among the tips) *and* the punctuational
signal (sisters differing where a burst split them) at once.

## Traits and genes

### `traits:genomes` — trait-linked gene families

A lineage that becomes aerobic retains and acquires oxygen-using gene families; one that reverts sheds
them. This edge simulates that link: a trait is evolved down the tree, then a panel of gene families is
evolved whose **loss depends on the local trait value**. Writing the trait on a branch as $s$, a
*present* family with coupling weight $w_i$ is lost at rate

$$\text{loss}_i = \text{base\_loss} \cdot \exp(-\text{effect\_loss} \cdot w_i \cdot s),$$

so where the trait favours a responsive family it is retained, and where it does not it is purged
faster than baseline; inert families ($w_i = 0$) always lose at `base_loss`. Gain is a trait-blind
horizontal influx, and the trait-modulated loss then *selectively retains* it — the same
retention mechanism as the gene-family coupling model, which is what writes the trait–gene association
into the profiles. The driver takes the tree, a **pre-simulated trait**, and the coupling:

```python
import numpy as np
from zombi2 import simulate_traits, BrownianMotion
from zombi2.coevolve import simulate_trait_conditioned_genomes, TraitGeneCoupling

trait = simulate_traits(tree, BrownianMotion(sigma2=1.0), seed=1)
weights = np.zeros(40)
weights[::3] = 1.0             # every third family responds to the trait
res = simulate_trait_conditioned_genomes(
    tree, trait,
    TraitGeneCoupling(n_families=40, weights=weights,
                      effect_loss=2.0, base_loss=1.0, transfer=0.5),
    seed=1)
res.profiles           # ProfileMatrix: families x species (the linked data)
res.trait              # the TraitResult the genes were conditioned on
res.genomes()          # promote to a full Genomes (gene trees, reconciliations)
```

The command builds the trait for you from `--trait-model` (any `zombi2 traits` model) and selects the
responsive families with `--responsive` (a count, a fraction, or an id list):

```bash
zombi2 coevolve --couple traits:genomes -t species_tree.nwk \
    --trait-model bm --sigma2 1 --panel 40 --responsive 0.3 \
    --loss 0.5 --trans 1 --effect-loss 2 --seed 1 -o out/
```

It writes the usual gene-family output plus `traits.tsv`/`trait_tree.nwk` and `coupling.tsv` (the
per-family weights, so the exact trait–gene linkage that generated the profiles is on record for
downstream inference). **What it recovers:** responsive families are present where the trait favours
them and absent where it does not, while inert families do not distinguish the states — the signal is
entirely in the responsive panel, which is what an inference method should be able to pick out.

![Trait-linked gene families. **A**, the mechanism: a responsive family's loss rate falls with the trait ($\text{loss} = \text{base\_loss}\cdot e^{-\text{effect\_loss}\,w\,s}$), while an inert family loses at the flat baseline. **B**, one realization — the tree painted by the trait (viridis), a per-tip trait chip, then the gene-presence matrix. The responsive block is filled in the high-trait (yellow) tips and empty in the low-trait (blue) ones, tracking the trait almost perfectly; the inert block carries no such pattern. The signal lives entirely in the responsive families — the null is built in.](figures/trait_linked_genes.pdf){width=100%}

### `genomes:traits` — gene-conditioned trait

The reverse: here gene content conditions a **trait**. A binary *modifier* gene is gained and lost
along the tree, and its presence sets a continuous trait's **OU optimum** — a lineage that acquires the
gene is pulled toward a new adaptive peak (`theta_present`), one that loses it drifts back toward
`theta_absent`. "Gene presence enables a phenotypic shift." Also an overlay on a given tree:

```python
from zombi2.coevolve import simulate_gene_conditioned_trait, GeneConditionedTrait

res = simulate_gene_conditioned_trait(
    tree, GeneConditionedTrait(gene_gain=0.6, gene_loss=0.6, theta_absent=0.0,
                               theta_present=5.0, alpha=2.5, sigma2=0.4), seed=2)
res.gene_presence()    # {leaf: 0/1} — modifier presence at the tips
res.trait_values()     # {leaf: value} — carriers near 5, non-carriers near 0
```

```bash
zombi2 coevolve --couple genomes:traits -t species_tree.nwk \
    --modifier-gain 0.6 --modifier-loss 0.6 --theta-absent 0 --theta-present 5 \
    --trait-alpha 2.5 --trait-sigma2 0.4 --seed 2 -o out/
```

`--modifier-gain`/`--modifier-loss` set the binary modifier's dynamics (`--root-modifier` to start it
present); `--theta-absent`/`--theta-present` are the two OU optima; `--trait-alpha` is the
mean-reversion strength (`0` = Brownian) and `--trait-sigma2` the diffusion. **What it recovers:** tips
carrying the modifier sit near `theta_present`, those without near `theta_absent` — a discrete genomic
event reading out as a shift in a continuous phenotype.

![Gene-conditioned trait. **A**, the mechanism: one lineage's trait sits near `theta_absent` while the modifier is absent (light), then the gene is gained (the +) and the trait climbs to the new OU optimum `theta_present` (heavy). **B**, one realization — the tree drawn heavy where the modifier is present, light where absent, and each tip's trait value as a dot on a shared axis. Carriers (filled) sit at `theta_present`, non-carriers (open) at `theta_absent`: a discrete genomic event reading out as a shift in a continuous phenotype.](figures/gene_conditioned_trait.pdf){width=100%}

### Both arrows: trait–gene feedback

Turn on **both** `traits:genomes` and `genomes:traits` and the trait and a coupled gene panel modulate
*each other*: the trait sets the panel's retention while the panel sets the trait's optimum. Neither
arrow points into S, so this stays an **overlay** on a given tree — but because each depends on the
other's *current* value, the two are integrated together along each branch rather than one after the
other. The loop is self-reinforcing (carrying the panel pulls the trait up; a high trait keeps the
panel), so a lineage settles into a panel-rich/high-trait or a panel-poor/low-trait regime, and the
tips end up with the trait and the panel **correlated even though neither was imposed**:

```python
from zombi2.coevolve import TraitGeneFeedback, simulate_trait_gene_feedback

res = simulate_trait_gene_feedback(
    tree,
    TraitGeneFeedback(n_families=24, effect_loss=1.5, base_loss=1.0, gain=1.0,
                      theta_low=-3.0, theta_high=3.0, alpha=1.0, sigma2=0.5),
    seed=2)
res.trait_gene_correlation()   # the emergent trait-gene association
```

```bash
zombi2 coevolve --couple traits:genomes --couple genomes:traits -t species_tree.nwk \
    --panel 24 --effect-loss 1.5 --loss 1.0 --trans 1.0 \
    --theta-absent -3 --theta-present 3 --trait-alpha 1 --trait-sigma2 0.5 \
    --seed 2 -o out/
```

Here `--theta-absent`/`--theta-present` are the trait's optima at an empty/full panel and `--loss`/
`--trans` the panel's base loss/gain. Setting `--effect-loss 0` recovers pure `genomes:traits`, and
setting `--theta-present` equal to `--theta-absent` recovers pure `traits:genomes` — the joint model
contains both single edges as limits. **What it recovers:** a trait–gene-content association that is
*emergent* rather than built in; the decoupled control (`--effect-loss 0` with equal thetas) shows
none.

## Null models of coevolution

Every model in this part is a *claim*: this driver shapes that target. The hard part is not
simulating the claim but knowing, from data, whether it is true — and that is an inference problem
with a famous trap. Fit a trait-dependent speciation model (BiSSE) and it will almost always report
that your trait drives diversification, *even for a trait that does not* [@raboskygoldberg2015]. The
reason is that real trees always carry rate heterogeneity from causes unrelated to the trait, and
BiSSE's only way to describe a fast-diversifying clade is "the trait did it." Its naïve null — a
constant-rate tree with *no* heterogeneity — is a strawman that real data beat for reasons that have
nothing to do with the character. Imagine a clade of bacteria in which the aerobic lineages happen
to be more diverse: did aerobiosis *cause* the radiation, or does it merely sit in a part of the tree
that was diversifying quickly anyway?

The fix is a null model that is as flexible as the alternative, minus the causal link: a dataset with
the same amount of variation in the target, but where that variation is **not** produced by the
driver [@beaulieu2016hisse]. ZOMBI2 gives every edge in this part exactly such a null, so
"simulate coupled → simulate the matched null → run your detector on both → measure its
false-positive rate" is a one-command workflow.

![A null keeps the tree's variation and cuts only the trait's grip on it. Coupled — the trait fills
the fast-diversifying clade, so it looks causal. Neutral — a balanced tree with no fast clade, a weak
test. CID — the *same* fast clade as the coupled tree, but the trait is scattered across fast and slow
clades: the honest test of whether the trait tracks diversification. Panels A and C are the same
tree.](figures/coevolve_null_archetypes.pdf){width=100%}

Cutting the arrow honestly takes one of three forms:

- **`neutral`** — set the coupling strength to zero. The arrow is cut and nothing compensates, so the
  target loses its coupling-induced variation. This is the naïve null (the constant-rate strawman
  above): cheap, and the honest baseline for *"does my detector fire when there is truly no effect?"*
- **`cid`** (character-independent diversification) — re-introduce the *same* variation, but source
  it from a **hidden** driver uncorrelated with the observed one. This is the generalised HiSSE null:
  the tree really has fast and slow clades, but the observed character cannot explain them (panel C
  above). It is the *worthy opponent*, and it covers the four edges where a driver state sets a target
  rate.
- **`timing`** — for the two edges where change happens *at* speciation (`species:traits`,
  `species:genomes`) there is no hidden state to invoke; the honest null keeps the same amount of change
  but spreads it **along the branches** instead of piling it at the nodes — the
  punctuation-versus-gradual contrast [@pagel1999inferring]. The variance is matched analytically, from the
  tree's branch statistics.

![The timing null: the same change, moved off the speciations. Coupled — change happens at each
speciation, so sister tips differ sharply (punctuational). Timing null — the same amount of change
spread along the branches, so sisters differ only as much as their shared branch length
allows.](figures/coevolve_null_timing.pdf){width=100%}

For two of the four `cid` edges the null is almost free, because ZOMBI2 already produces a neutral
channel. In `genomes:species` the drivers shape a genuinely heterogeneous tree, while the **neutral bulk
genome** — the families that do not touch diversification — is a whole panel of real genes decoupled
from it. The null hands you that genome as the observed data and withholds the drivers as
ground-truth; `genomes:traits` reuses the same trick. Only `traits:genomes` needs one extra ingredient: a
second, independent neutral trait. (A neutral gene still carries a faint imprint of tree shape — bushy
clades have short branches — but that shared-tree confound is a *feature* of a good null: it is exactly
what a trustworthy detector must see through.)

Generate any null by adding `--null` to the edge's command:

```bash
# the coupled claim: a trait drives diversification
zombi2 coevolve --couple traits:species \
    --lambda0 1 --lambda1 3 --tips 200 --seed 1 -o out/alt

# its matched CID-2 null: same rate spread, no trait effect
zombi2 coevolve --couple traits:species \
    --lambda0 1 --lambda1 3 --tips 200 --seed 1 \
    --null cid --hidden 2 -o out/null

# the punctuational genome, spread along branches (species:genomes)
zombi2 coevolve --couple species:genomes -t species_tree.nwk \
    --clado-gene-loss 0.15 --clado-gene-gain 3 \
    --null timing -o out/null_punct
```

Every null run also writes a `null_manifest.tsv` recording which arrow was cut and how the target's
variance was preserved, so a downstream calibration is self-documenting. In Python the same nulls are
a coupling model's `.null(kind=...)` method, plus a `CID` factory for the `traits:species` case.
