# Coupled models

By default ZOMBI2 simulates in a **pipeline**: grow a species tree, then evolve a trait along it,
then evolve gene families along it. That works because the joint distribution *factorises* —
`P(tree)·P(trait | tree)·P(genes | tree, trait)` — so each stage runs on the frozen output of the
previous one. **Coupled models break that factorisation.** A coupling is a directed edge
`driver → target`: the driver's state modulates the target's rates. When the arrow points *into* the
tree (a trait or gene sets speciation/extinction) the tree can no longer be drawn first — it is an
**output**, grown forward together with the driver; when it points elsewhere the coupling is an
**overlay** on a tree you supply with `-t`. All of these live under one command,
`zombi2 coevolve --couple driver:target`, over the three nodes {species, traits, genes}.

| Model | Edge (`driver:target`) | Tree | Reach for it when |
| --- | --- | --- | --- |
| **GeneDiversification** | genes:species | output (forward) | gene content (key innovations, spread by HGT) drives the radiation |
| **CladogeneticGenome** | species:genes | given (`-t`) | gene gain/loss bursts *at* speciations — punctuational genome evolution |
| **GeneConditionedTrait** | genes:traits | given (`-t`) | a modifier gene unlocks a new phenotypic optimum |
| **TraitGeneCoupling** (trait-linked genomes) | traits:genes | given (`-t`) | a trait's history shapes which gene families are retained |
| **co-diversification** | genes:species + species:genes | output (forward) | speciation itself reshuffles the rate-setting drivers (genomic ClaSSE) |
| **TraitGeneFeedback** | traits:genes + genes:traits | given (`-t`) | a trait and a gene panel modulate each other, with no single imposed arrow |

SSE / ClaSSE (`traits:species`, `species:traits`) are the trait→species edges of the same command;
they are documented with the [diversification models](diversification.md).

## The models

### GeneDiversification (genes:species)

A small panel of binary **driver** ("key innovation") gene families whose *presence* sets each
lineage's speciation/extinction rate: a present driver scales λ by `exp(driver_speciation)` and μ by
`exp(driver_extinction)`. Drivers are gained de novo (`origination`) and — the interesting part — by
**frequency-dependent transfer** (`transfer`: a driver in more live genomes spreads faster), and lost
at `loss`. Because gain depends on the live population, the tree and the drivers must grow together,
so this arrow-into-S edge is **forward-only** and produces the tree (take `age` or `n_tips`, no `-t`).
The neutral bulk genome, which does not touch diversification, is overlaid afterward on the finished
tree with the ordinary [`genomes`](../cli.md) (exact under independent families). `root_drivers` seeds
the first *m* drivers at the root.

### CladogeneticGenome (species:genes)

The reverse of `genes:species`: gene content does **not** affect diversification, so this is an
overlay on a given tree — the genomic twin of cladogenetic trait evolution. A genome of
`initial_families` families is evolved down the tree with a **founder-effect burst** at every
speciation — a daughter drops each family it carries with probability `cladogenetic_loss` and gains a
Poisson(`cladogenetic_gain`) count of new families — on top of optional gradual along-branch change
(`loss`, `origination`). With both anagenetic rates 0 the change is **purely punctuational**, and its
signature is that *sister tips differ* because change is injected at their split rather than spread
along branches.

### GeneConditionedTrait (genes:traits)

The reverse of `traits:genes`: here gene content conditions a **trait**. A binary *modifier* gene
comes and goes along the tree (a two-state Markov chain, `gene_gain`/`gene_loss`, optionally
`root_gene` present at the root), and its presence sets a continuous trait's **OU optimum**: a lineage
carrying the gene is pulled toward `theta_present`, one without it drifts back to `theta_absent`, with
mean-reversion `alpha` (0 = Brownian) and diffusion `sigma2`. "Gene presence enables a trait shift."
Tips carrying the modifier end up near `theta_present`, those without near `theta_absent`.

### TraitGeneCoupling (traits:genes)

The trait-linked-genomes model: a trait is evolved down the tree, then a fixed **panel** of gene
families is evolved along it whose **loss depends on the local trait value**. A *responsive* family
carrying weight `w` is lost at `base_loss · exp(-effect_loss · w · s)` for trait value `s`, so it is
retained where the trait favours it and purged where it does not; inert families (`w = 0`) always lose
at `base_loss`. **Gain is field-blind transfer** — a constant influx that the trait-modulated loss
then selectively retains, so net gene content tracks the trait without the influx itself seeing it.
`TraitGeneCoupling.build(n_families, responsive, ...)` picks the responsive set (a count, a fraction,
or an explicit id/index list; `signed` randomises weight signs). The trait is followed exactly in
time — a discrete trait contributes its stochastic character map, a continuous one is sub-segmented
into `trait_steps` pieces per branch. Setting `effect_loss = 0` recovers plain, uncoupled gene-family
evolution as a null.

### co-diversification (genes:species + species:genes)

Both species↔genes arrows at once: the same driver panel **sets** the diversification rates
(`genes:species`) *and* is **reshuffled by a cladogenetic burst** at every speciation
(`species:genes`: a daughter drops each carried driver with probability `cladogenetic_loss` and gains
each absent one with probability `cladogenetic_gain`). Because a burst can hand one daughter a key
innovation and not its sister, speciation *itself* seeds rate heterogeneity — the genomic analogue of
ClaSSE. One arrow points into S, so the tree is an **output** (`simulate_co_diversification`, or
`--couple genes:species --couple species:genes`). It reduces to `GeneDiversification` when both
cladogenetic probabilities are 0.

### TraitGeneFeedback (traits:genes + genes:traits)

Both traits↔genes arrows at once: a continuous trait and a coupled panel of `n_families` modulate each
other, integrated jointly along each branch — the panel's present count sets the trait's OU optimum
(interpolated between `theta_low` at an empty panel and `theta_high` at a full one), while the trait
sets each responsive family's retention exactly as in `traits:genes` (`effect_loss`, `base_loss`,
`gain`). No single edge is imposed, yet the tips end up correlated. It is an overlay on a given tree
and contains its two single edges as limits. `root_fraction` seeds the panel at the root.

## Command line

`--couple driver:target` selects the edge(s); the order reads as the arrow. Edges into species grow
the tree (`--age`/`--tips`, no `-t`); the other edges overlay a tree you pass with `-t`. All runs
below reuse a tree from `zombi2 species`.

```bash
T=out/species_tree.nwk

# genes:species — key-innovation drivers spread by HGT drive the radiation (tree is an output)
zombi2 coevolve --couple genes:species --drivers 2 --root-drivers 1 \
    --lambda0 1 --mu0 0.2 --driver-speciation 1.2 --driver-transfer 0.8 --driver-loss 0.3 \
    --tips 60 --seed 1 -o keygene/
# overlay the neutral genome on the grown tree (the factorization, made explicit)
zombi2 genomes -t keygene/species_tree.nwk --trans 1 --loss 0.5 --write profiles trees -o keygene/

# species:genes — purely punctuational genome: change ONLY at speciations
zombi2 coevolve --couple species:genes -t $T \
    --genome-size 30 --clado-gene-loss 0.15 --clado-gene-gain 3 --seed 2 -o punct/

# genes:traits — a modifier gene unlocks a phenotypic optimum at 5 (vs 0 without it)
zombi2 coevolve --couple genes:traits -t $T \
    --modifier-gain 0.6 --modifier-loss 0.6 --theta-absent 0 --theta-present 5 \
    --trait-alpha 2.5 --trait-sigma2 0.4 --seed 2 -o genetrait/

# traits:genes — 30% of a 40-family panel respond to a binary trait
zombi2 coevolve --couple traits:genes -t $T \
    --trait-model mk --states 2 --rate 0.3 --trait-center \
    --panel 40 --responsive 0.3 --weight 1 --effect-loss 3 \
    --loss 0.4 --trans 1.0 --write all --seed 7 -o tgenes/

# co-diversification — both species<->genes arrows (tree is an output)
zombi2 coevolve --couple genes:species --couple species:genes --drivers 3 --root-drivers 1 \
    --lambda0 1 --mu0 0.2 --driver-speciation 1.2 --driver-clado-loss 0.2 --driver-clado-gain 0.1 \
    --tips 60 --seed 4 -o codiv/

# trait-gene feedback — both traits<->genes arrows on a given tree
zombi2 coevolve --couple traits:genes --couple genes:traits -t $T \
    --panel 30 --effect-loss 2 --theta-absent 0 --theta-present 5 \
    --panel-root-fraction 0.5 --seed 5 -o feedback/
```

## Python

The models and drivers are in `zombi2.coevolve` (each name also aliases at the top level, so
`zombi2.simulate_gene_diversification` works too):

```python
import zombi2 as z
from zombi2.coevolve import (
    GeneDiversification, simulate_gene_diversification, simulate_co_diversification,
    CladogeneticGenome, simulate_cladogenetic_genome,
    GeneConditionedTrait, simulate_gene_conditioned_trait,
    TraitGeneCoupling, simulate_trait_linked_genomes,
    TraitGeneFeedback, simulate_trait_gene_feedback,
)

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=40, age=5.0, seed=1)

# genes:species — the tree is grown jointly with the drivers (no -t; take age or n_tips)
gd = simulate_gene_diversification(
    GeneDiversification(2, lambda0=1.0, mu0=0.2, driver_speciation=1.2,
                        transfer=0.8, loss=0.3, root_drivers=1),
    n_tips=60, seed=1)
gd.tree                      # the complete tree the drivers shaped (extinct lineages kept)
gd.node_drivers              # per-node driver presence

# species:genes — an overlay on a given tree
cg = simulate_cladogenetic_genome(
    tree, CladogeneticGenome(30, cladogenetic_loss=0.15, cladogenetic_gain=3), seed=2)
cg.profile_matrix().presence()             # families x extant tips (0/1)

# genes:traits — a modifier gene sets a continuous trait's OU optimum
gct = simulate_gene_conditioned_trait(
    tree, GeneConditionedTrait(gene_gain=0.6, gene_loss=0.6,
                               theta_absent=0, theta_present=5, alpha=2.5, sigma2=0.4), seed=2)

# traits:genes — trait-linked gene families (center a binary trait for two-sided coupling)
coupling = TraitGeneCoupling.build(40, 0.3, weight=1.0, effect_loss=3.0,
                                   base_loss=0.5, transfer=1.0,
                                   state_values=[-1.0, 1.0], seed=1)
tl = simulate_trait_linked_genomes(tree, z.Mk.equal_rates(2, 0.4), coupling, seed=2)
tl.profiles.presence()       # panel families x extant species (the trait-linked data)
tl.trait                     # the trait the genes were conditioned on

# the joint models
simulate_co_diversification(
    GeneDiversification(3, lambda0=1.0, mu0=0.2, driver_speciation=1.2,
                        cladogenetic_loss=0.2, cladogenetic_gain=0.1, root_drivers=1),
    n_tips=60, seed=4)
simulate_trait_gene_feedback(
    tree, TraitGeneFeedback(n_families=30, effect_loss=2.0,
                            theta_low=-3, theta_high=3, root_fraction=0.5), seed=5)
```

## Output

Every run writes `coevolve.log` (the manifest) and, for overlay edges, echoes the `-t` tree as
`species_tree.nwk`. The into-species edges (`genes:species`, co-diversification) instead **write** the
grown `species_tree.nwk` (complete, extinct lineages kept) plus `drivers.tsv` (per-node driver
presence) and `drivers_manifest.tsv` (the effect sizes β and rates). `species:genes` writes
`Profiles.tsv`/`Presence.tsv` (families × extant tips) and `genome_sizes.tsv`. The trait-side edges
(`genes:traits`, `traits:genes`, feedback) always write `traits.tsv` (the trait at every node);
`genes:traits` and `traits:genes` also write `trait_tree.nwk`, and `traits:genes` additionally writes the usual gene-family files chosen with
`--write` ({profiles, trace, trees, events, transfers, summary}, or `all`) plus `coupling.tsv`, the
responsive-family manifest that records the exact trait↔gene linkage for downstream inference.

## Validation

- **GeneDiversification** — with two drivers, one a strong speciation driver, the extant tips are
  biased toward carrying it
  (`test_gene_diversification.py::test_speciation_driver_biases_tips_toward_it`).
- **CladogeneticGenome** — under pure punctuation, sister tips differ because the change is injected
  at their split, not along the branches
  (`test_cladogenetic_genome.py::test_sister_tips_differ_because_change_is_at_the_split`).
- **GeneConditionedTrait** — tips that carry the modifier track the present optimum `theta_present`,
  those without it the absent one
  (`test_gene_conditioned_trait.py::test_carriers_track_present_optimum`).
- **TraitGeneCoupling** — inject a strong coupling and recover it: responsive families track the
  trait across the tips while inert families do not
  (`test_trait_coupling.py::test_inject_recover_trait_tracks_responsive_families`).
- **co-diversification** — the cladogenetic burst differentiates sisters, seeding the rate
  heterogeneity (`test_co_diversification.py::test_burst_differentiates_sisters`).
- **TraitGeneFeedback** — the joint run writes a measurable trait↔gene-panel correlation into the
  tips (`test_trait_gene_feedback.py::test_feedback_writes_a_trait_gene_correlation`).

## References

- Maddison, Midford & Otto (2007). Estimating a binary character's effect on speciation and
  extinction. *Systematic Biology* 56(5): 701–710. (BiSSE — the state-dependent template.)
- Goldberg & Igić (2012). Tempo and mode in plant breeding system evolution. *Evolution* 66(12):
  3701–3709. (ClaSSE — cladogenetic state change, the analogue behind the joint models.)
- Bokma, F. (2008). Detection of "punctuated equilibrium" by Bayesian estimation. *Journal of
  Evolutionary Biology* 21(5): 1218–1227. (Change concentrated at speciation.)
- Davin, A. A. (2025). Timing the tree of bacteria with horizontally transferred genes and a trait
  linked to the Great Oxidation Event. (The trait-linked gene-family generator.)
