# Coevolution (coupled models)

By default ZOMBI2 simulates in a **pipeline**: grow a species tree, then evolve a trait along it,
then evolve gene families along it. That works because the joint distribution *factorises* —
`P(tree)·P(trait | tree)·P(genes | tree, trait)` — so each stage runs on the frozen output of the
previous one. **Coevolution breaks that factorisation.** A coupling is a directed edge
`driver → target`: the driver's state modulates the target's rates, and the two levels must be grown
**together** — one process, one Gillespie. All of them live under a single command,
`zombi2 coevolve --couple driver:target`, over the four levels {species, traits, genes}.

The direction of the arrow decides how the run behaves. When the arrow points *into* the tree — a
trait or a gene panel sets speciation/extinction — the tree can no longer be drawn first; it is an
**output**, grown forward together with the driver (give a stopping condition, `--age` or `--tips`,
no `-t`). When it points elsewhere, the coupling is an **overlay** on a tree you supply with `-t`.
State-dependent diversification (SSE) is simply the family of couplings whose arrow points into the
species tree.

| Model | Edge (driver→target) | Tree | What it does |
| --- | --- | --- | --- |
| **BiSSE** | traits→species | output (forward) | a binary trait sets each lineage's λ, μ; the fast state comes to dominate the tips |
| **MuSSE** | traits→species | output (forward) | a k-state character drives diversification (BiSSE is `k = 2`) |
| **HiSSE** | traits→species | output (forward) | binary observed + hidden classes; heterogeneity can live off the observed trait |
| **QuaSSE** | traits→species | output (forward) | a continuous trait sets a bounded λ(x) |
| **ClaSSE** | traits→species + species→traits | output (forward) | SSE *plus* a state jump *at* each speciation (cladogenetic) |
| **GeneDiversification** | genes→species | output (forward) | gene content (key innovations, spread by HGT) drives the radiation |
| **CladogeneticGenome** | species→genes | given (`-t`) | gene gain/loss bursts *at* speciations — punctuational genome evolution |
| **GeneConditionedTrait** | genes→traits | given (`-t`) | a modifier gene unlocks a new phenotypic optimum |
| **TraitGeneCoupling** (trait-linked genomes) | traits→genes | given (`-t`) | a trait's history shapes which gene families are retained |
| **co-diversification** | genes→species + species→genes | output (forward) | speciation itself reshuffles the rate-setting drivers (genomic ClaSSE) |
| **TraitGeneFeedback** | traits→genes + genes→traits | given (`-t`) | a trait and a gene panel modulate each other, with no single imposed arrow |

## The models

### State-dependent diversification (SSE)

State-dependent speciation and extinction (SSE) models let a **trait drive the shape of the tree**: a
lineage's character state sets its speciation and extinction rates, so the tree and the trait must be
grown **together** (Maddison, Midford & Otto 2007). This is the `traits:species` edge — an arrow
*into* species, so the run is **forward-only** and **produces** the tree rather than taking one. The
variants differ only in the kind of trait doing the driving: binary, k-state, hidden, continuous, or a
trait that also jumps *at* each speciation.

**BiSSE.** Two states (`0`, `1`), each with its own speciation rate (`lambda0`/`lambda1`) and
extinction rate (`mu0`/`mu1`), plus asymmetric anagenetic transitions (`q01`, `q10`). The classic
binary state-dependent model (Maddison, Midford & Otto 2007); the fast-speciating state comes to
dominate the standing tips. The default `--sse-model`.

**MuSSE.** The k-state generalisation: length-`k` `birth` and `death` rate vectors and a `k × k`
anagenetic `Q` matrix (off-diagonals ≥ 0; the diagonal is recomputed so rows sum to zero, exactly as
in [`Mk`](../guide/traits.md)). Use it when a multi-state character — not just a binary one — drives
diversification (FitzJohn 2012, *diversitree*). BiSSE is the `k = 2` special case.

**HiSSE.** Extends BiSSE with unobserved **hidden classes**: each observed state comes in `H`
variants, each its own diversification regime, with switch rates between classes (Beaulieu & O'Meara
2016). It is the honest null for SSE inference — rate heterogeneity that lives on a *hidden* class is
not falsely pinned on the observed character. Build it from one `BiSSE` per hidden class plus a
`hidden_transition` matrix (or a scalar for a symmetric rate); the tips report the **observed** state,
with the `(observed, hidden)` pair still available per node. **Python-only** (no `--sse-model hisse`).

**QuaSSE.** A **continuous** trait diffuses (Brownian motion, `sigma2`, optional `drift`) along every
lineage and the rates are functions of its current value (FitzJohn 2010). The rate functions must be
**bounded** — an unbounded λ(x) under a diffusing x has no valid thinning bound — so you pass a
`rate_bound` on λ(x) + μ(x); `QuaSSE.sigmoid(low, high, center, slope)` builds a convenient bounded
speciation curve. On the CLI the trait is a sigmoidal speciation (`--spec-low/high/center/slope`) plus
a constant extinction (`--qmu`).

**ClaSSE.** Not a separate class but the **both-arrows** combination: a discrete or continuous SSE
model *plus* a [`Cladogenesis`](../guide/traits.md) kernel that jumps each daughter's state **at**
speciation (Goldberg & Igić 2012). The trait both shapes the tree (`traits:species`) *and* is kicked
by its branching (`species:traits`), so change is concentrated at nodes rather than spread along
branches. `shift` is the per-daughter state-hop probability (discrete); `jump_sigma2` is the Gaussian
jump variance (continuous, `quasse`). With `Q = 0` the cladogenesis kernel supplies all the state
dynamics.

### GeneDiversification (genes→species)

A small panel of binary **driver** ("key innovation") gene families whose *presence* sets each
lineage's speciation/extinction rate: a present driver scales λ by `exp(driver_speciation)` and μ by
`exp(driver_extinction)`. Drivers are gained de novo (`origination`) and — the interesting part — by
**frequency-dependent transfer** (`transfer`: a driver in more live genomes spreads faster), and lost
at `loss`. Because gain depends on the live population, the tree and the drivers must grow together,
so this arrow-into-S edge is **forward-only** and produces the tree (take `age` or `n_tips`, no `-t`).
The neutral bulk genome, which does not touch diversification, is overlaid afterward on the finished
tree with the ordinary [`genomes`](../cli.md) (exact under independent families). `root_drivers` seeds
the first *m* drivers at the root.

### CladogeneticGenome (species→genes)

The reverse of `genes:species`: gene content does **not** affect diversification, so this is an
overlay on a given tree — the genomic twin of cladogenetic trait evolution. A genome of
`initial_families` families is evolved down the tree with a **founder-effect burst** at every
speciation — a daughter drops each family it carries with probability `cladogenetic_loss` and gains a
Poisson(`cladogenetic_gain`) count of new families — on top of optional gradual along-branch change
(`loss`, `origination`). With both anagenetic rates 0 the change is **purely punctuational**, and its
signature is that *sister tips differ* because change is injected at their split rather than spread
along branches.

### GeneConditionedTrait (genes→traits)

The reverse of `traits:genes`: here gene content conditions a **trait**. A binary *modifier* gene
comes and goes along the tree (a two-state Markov chain, `gene_gain`/`gene_loss`, optionally
`root_gene` present at the root), and its presence sets a continuous trait's **OU optimum**: a lineage
carrying the gene is pulled toward `theta_present`, one without it drifts back to `theta_absent`, with
mean-reversion `alpha` (0 = Brownian) and diffusion `sigma2`. "Gene presence enables a trait shift."
Tips carrying the modifier end up near `theta_present`, those without near `theta_absent`.

### TraitGeneCoupling (traits→genes)

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

### co-diversification (genes→species + species→genes)

Both species↔genes arrows at once: the same driver panel **sets** the diversification rates
(`genes:species`) *and* is **reshuffled by a cladogenetic burst** at every speciation
(`species:genes`: a daughter drops each carried driver with probability `cladogenetic_loss` and gains
each absent one with probability `cladogenetic_gain`). Because a burst can hand one daughter a key
innovation and not its sister, speciation *itself* seeds rate heterogeneity — the genomic analogue of
ClaSSE. One arrow points into S, so the tree is an **output** (`simulate_co_diversification`, or
`--couple genes:species --couple species:genes`). It reduces to `GeneDiversification` when both
cladogenetic probabilities are 0.

### TraitGeneFeedback (traits→genes + genes→traits)

Both traits↔genes arrows at once: a continuous trait and a coupled panel of `n_families` modulate each
other, integrated jointly along each branch — the panel's present count sets the trait's OU optimum
(interpolated between `theta_low` at an empty panel and `theta_high` at a full one), while the trait
sets each responsive family's retention exactly as in `traits:genes` (`effect_loss`, `base_loss`,
`gain`). No single edge is imposed, yet the tips end up correlated. It is an overlay on a given tree
and contains its two single edges as limits. `root_fraction` seeds the panel at the root.

## Command line

`--couple driver:target` selects the edge(s); the order reads as the arrow. Edges into species grow
the tree (`--age`/`--tips`, no `-t`); the other edges overlay a tree you pass with `-t`. For SSE,
`--sse-model` picks `bisse` (default), `musse`, or `quasse`, and `--root-state` sets the root state
index for `bisse`/`musse` (default: the character's stationary distribution). HiSSE is not exposed on
`--sse-model` — use the Python API for it. The overlay runs below reuse a tree from `zombi2 species`.

```bash
T=out/species_tree.nwk

# --- SSE: a trait grows the tree (traits:species; no -t, give --tips or --age) ---

# BiSSE: state 1 speciates 3x faster, so it comes to dominate the standing tips
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/bisse

# MuSSE: k-state — birth/death vectors + a k x k transition-rate matrix file
printf "0 0.2 0.2\n0.2 0 0.2\n0.2 0.2 0\n" > q3.txt
zombi2 coevolve --couple traits:species --sse-model musse \
    --birth 1 1 3 --death 0.2 0.2 0.2 --q-matrix q3.txt --tips 200 --seed 1 -o out/musse

# QuaSSE: continuous trait — sigmoidal speciation + constant extinction + Brownian diffusion
zombi2 coevolve --couple traits:species --sse-model quasse \
    --spec-low 0.4 --spec-high 3 --spec-center 0 --spec-slope 3 \
    --qmu 0.2 --diffusion 0.5 --root-value -1.5 --tips 200 --seed 1 -o out/quasse

# ClaSSE: both arrows — BiSSE rates + a cladogenetic state hop at each speciation
zombi2 coevolve --couple traits:species --couple species:traits \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.05 --q10 0.05 \
    --clado-shift 0.3 --tips 200 --seed 3 -o out/classe

# --- gene- and trait-coupled edges ---

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

Run the CLI as `python -m zombi2 coevolve ...` (not a bare `zombi2`) if the entry point is not on your
PATH.

## Python

The models and drivers live in `zombi2.coevolve`; the cladogenetic kernel in `zombi2.traits` (each
also re-exports at the top level, so `zombi2.BiSSE` / `zombi2.simulate_gene_diversification` /
`zombi2.Cladogenesis` work too):

```python
import zombi2 as z
from zombi2.coevolve import (
    BiSSE, MuSSE, HiSSE, QuaSSE, simulate_sse,
    GeneDiversification, simulate_gene_diversification, simulate_co_diversification,
    CladogeneticGenome, simulate_cladogenetic_genome,
    GeneConditionedTrait, simulate_gene_conditioned_trait,
    TraitGeneCoupling, simulate_trait_linked_genomes,
    TraitGeneFeedback, simulate_trait_gene_feedback,
)
from zombi2.traits import Cladogenesis

# --- SSE: the tree is grown jointly with the driving trait (no -t; take age or n_tips) ---

# BiSSE: state 1 speciates 3x faster -> tips biased toward state 1
res = simulate_sse(BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2, q01=0.1, q10=0.1),
                   n_tips=200, seed=1)
res.tree                    # complete tree (extinct lineages kept; z.prune() for the reconstructed one)
res.labeled_values()        # the observed trait at the extant tips

# MuSSE: k-state — birth/death vectors + a k x k Q
Q = [[0, 0.2, 0.2], [0.2, 0, 0.2], [0.2, 0.2, 0]]
musse = simulate_sse(MuSSE(birth=[1, 1, 3], death=[0.2, 0.2, 0.2], Q=Q), n_tips=200, seed=1)

# HiSSE: hidden classes drive the tree while the observed character stays neutral
fast, slow = BiSSE(2.5, 2.5, 0.2, 0.2, 0.3, 0.3), BiSSE(0.4, 0.4, 0.2, 0.2, 0.3, 0.3)
hisse = simulate_sse(HiSSE([fast, slow], hidden_transition=0.15), age=1.5, seed=0)

# QuaSSE: a bounded sigmoidal speciation on a Brownian continuous trait
spec = QuaSSE.sigmoid(low=0.4, high=3.0, center=0.0, slope=3.0)
quasse = simulate_sse(QuaSSE(spec, lambda x: 0.2, sigma2=0.5, rate_bound=3.2, x0=-1.5),
                      age=2.5, seed=0)

# ClaSSE: BiSSE + a cladogenetic jump at each speciation (both arrows)
classe = simulate_sse(BiSSE(1, 3, 0.2, 0.2, 0.05, 0.05),
                      cladogenesis=Cladogenesis(shift=0.3), n_tips=200, seed=3)

# --- gene- and trait-coupled edges ---

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

For the into-species edges (SSE, `genes:species`, co-diversification) provide exactly one stopping
condition: `age` (fixed crown age, random tip count) or `n_tips` (grow until this many extant tips
first coexist, random age); the run is conditioned on at least two extant survivors. `simulate_sse`
returns a `TraitResult`: `.tree` is the complete tree (extinct leaves carry `is_extant=False`),
`.values` are the extant tips' states, and `.history` is the realized character map (discrete models;
`None` for QuaSSE).

## Output

Every run writes `coevolve.log` (the manifest) and, for overlay edges, echoes the `-t` tree as
`species_tree.nwk`.

The **into-species** edges grow the tree, so they **write** `species_tree.nwk` (complete, extinct
lineages kept). SSE additionally writes `traits.tsv` (every node — tips *and* ancestral states) and
`trait_tree.nwk` (the trait annotated on every node); prune to the reconstructed, survivors-only tree
with `zombi2.prune(result.tree)` for downstream analysis. The gene-driven into-species edges
(`genes:species`, co-diversification) instead write `drivers.tsv` (per-node driver presence) and
`drivers_manifest.tsv` (the effect sizes β and rates).

The **overlay** edges write onto the tree you passed. `species:genes` writes
`Profiles.tsv`/`Presence.tsv` (families × extant tips) and `genome_sizes.tsv`. The trait-side edges
(`genes:traits`, `traits:genes`, feedback) always write `traits.tsv` (the trait at every node);
`genes:traits` and `traits:genes` also write `trait_tree.nwk`, and `traits:genes` additionally writes
the usual gene-family files chosen with `--write` ({profiles, trace, trees, events, transfers,
summary}, or `all`) plus `coupling.tsv`, the responsive-family manifest that records the exact
trait↔gene linkage for downstream inference.

## Validation

- **BiSSE** — a state that speciates 3× faster strongly biases the standing tips toward it
  (`test_sse.py::test_sse_faster_speciation_biases_tips`).
- **MuSSE** — with three states sharing extinction and a symmetric transition matrix, the
  fastest-speciating state (3× the others) is over-represented among the standing tips, far above the
  1/3 state-independent baseline
  (`test_sse.py::test_musse_fastest_speciation_state_over_represented_in_tips`).
- **HiSSE** — the fast *hidden* class dominates the tips while the *observed* character stays neutral
  (`test_sse.py::test_hisse_hidden_drives_diversification_not_observed`).
- **QuaSSE** — when speciation rises with the trait, surviving tips are biased to high values versus a
  constant-rate null (`test_sse.py::test_quasse_x_dependent_speciation_biases_the_trait`).
- **ClaSSE** — with anagenetic diffusion switched off, each parent→child step is exactly one
  cladogenetic jump, and those jumps are distributed `Normal(0, jump_sigma2)`: the empirical mean is
  ~0 and the empirical variance matches the `jump_sigma2` parameter to several sigma
  (`test_sse.py::test_classe_continuous_jumps_are_normal_zero_jump_sigma2`).
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

- Maddison, W. P., Midford, P. E. & Otto, S. P. (2007). Estimating a binary character's effect on
  speciation and extinction. *Systematic Biology* 56(5): 701–710. (BiSSE — the state-dependent
  template.)
- FitzJohn, R. G. (2010). Quantitative traits and diversification. *Systematic Biology* 59(6):
  619–633. (QuaSSE)
- FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R.
  *Methods in Ecology and Evolution* 3(6): 1084–1092. (MuSSE)
- Beaulieu, J. M. & O'Meara, B. C. (2016). Detecting hidden diversification shifts in models of
  trait-dependent speciation and extinction. *Systematic Biology* 65(4): 583–601. (HiSSE)
- Goldberg, E. E. & Igić, B. (2012). Tempo and mode in plant breeding system evolution. *Evolution*
  66(12): 3701–3709. (ClaSSE — cladogenetic state change + SSE, the analogue behind the joint models.)
- Bokma, F. (2008). Detection of "punctuated equilibrium" by Bayesian estimation. *Journal of
  Evolutionary Biology* 21(5): 1218–1227. (Change concentrated at speciation.)
- Davin, A. A. (2025). Timing the tree of bacteria with horizontally transferred genes and a trait
  linked to the Great Oxidation Event. (The trait-linked gene-family generator.)
