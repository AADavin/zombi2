# Coevolution

## The idea

ZOMBI2 normally simulates in a **pipeline**: it builds a species tree, then evolves a trait along it,
then evolves gene families along it. This works because the couplings between the three axes are
*one-directional* and the joint distribution factorises,

$$P(\text{tree}) \cdot P(\text{trait} \mid \text{tree}) \cdot P(\text{genes} \mid \text{tree}, \text{trait}),$$

so each stage can be drawn on the frozen output of the previous one. The moment a trait or a genome
feeds **back** into the process that generated the tree, this factorisation breaks: the tree, the
trait and the gene content must be grown *together*. The **`coevolve`** mode is for these coupled
scenarios.

There are three processes:

- **S** — species diversification (the birth–death process that grows the tree),
- **T** — a phenotypic trait (BM/OU/EB/Mk/threshold),
- **G** — gene-family content (the DTL process).

A **coupling is a directed edge** `driver -> target`: the driver's state modulates the target's
rates. A coevolution scenario is a **set of directed edges** on these three nodes. That single
abstraction covers everything from state-dependent diversification to full three-way feedback, and it
makes the one thing that matters — *direction* — explicit.

On the command line the coupling is one repeatable flag, `--couple driver:target`, where the order
reads as the arrow (driver first):

```bash
zombi2 coevolve --couple traits:species  ...   # T->S : the trait sets speciation/extinction (SSE)
zombi2 coevolve --couple species:traits  ...   # S->T : speciation drives the trait (cladogenetic)
zombi2 coevolve --couple traits:species --couple species:traits  ...   # both arrows = ClaSSE
```

So `--couple species:traits` and `--couple traits:species` are deliberately **different models**, and
a bidirectional coupling is simply *both* edges. The `:` (rather than `->`) keeps the flag shell-safe.

The one rule that governs their difficulty is: **does any active edge point into S?** If no edge
points into S, the tree is fixed — it is read from `-t/--tree` and every coupling is an *overlay* on a
frozen tree. If an edge does point into S (`traits:species` or `genes:species`), the tree topology
depends on the coupled state and cannot be drawn first: the tree becomes an **output**, and those runs
are forward-only and take no `-t`.

::: note
Each `--couple` edge is a distinct model, not a switch on one monolithic engine. Overlay edges (no
arrow into S) reuse the existing trait and gene-family engines unchanged; only the into-S edges leave
the fast paths.
:::

The remainder of this chapter documents the coupling with the most developed workflow — trait-linked
gene families — in depth, then summarises the state-dependent diversification edge.

## Trait-linked gene families

Gene families and phenotypic traits do not evolve independently. A lineage that becomes aerobic
retains and acquires oxygen-using gene families; one that reverts to anaeroby sheds them. ZOMBI2
simulates that link directly: it evolves a trait down the tree, then evolves a panel of gene families
whose **loss and gain depend on the local trait value**. The resulting phylogenetic profile carries a
*known*, trait-linked signal — the forward generator behind studies that read gene content as a
record of a trait's history, such as timing the bacterial tree from the Great Oxidation Event
[@davin2025goe].

This is the `traits:genes` edge. It runs on a **given** tree, and is exposed both as the Python driver
`simulate_trait_linked_genomes` and as the command `zombi2 coevolve-genetrait`.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=60, age=6.0, seed=1)

# a binary aerobic(1)/anaerobic(0) trait, then genes conditioned on it
coupling = z.TraitGeneCoupling.build(n_families=40, responsive=0.3, weight=1.0,
                                     effect_loss=3.0, base_loss=0.5, transfer=1.0, seed=1)
res = z.simulate_trait_linked_genomes(tree, z.Mk.equal_rates(2, 0.4), coupling, seed=2)

res.profiles.presence()        # panel families × extant species (0/1) — the trait-linked data
res.trait.labeled_values()     # the trait at the tips, from the same run
```

The trait models are those of the trait-evolution chapter and the coupled rate machinery is that of
the gene-family coupling model; only the family-side rate model is new, so the whole output pipeline
(profiles, gene trees, reconciliations) applies unchanged.

### The model

A fixed **panel** of $N$ gene families is seeded present at the root. Each family carries a coupling
**weight** $w_i$ ($0$ = inert). Writing the trait value on a branch as $s$, a *present* family is lost
at rate

$$\text{loss}_i = \text{base\_loss} \cdot \exp(-\text{effect\_loss} \cdot w_i \cdot s),$$

so where the trait favours a responsive family ($w_i \cdot s$ large and positive) it is retained, and
where it does not ($w_i \cdot s$ negative) it is purged faster than the baseline. Inert families
($w_i = 0$) always lose at `base_loss`.

**Gain is horizontal transfer** — a field-blind influx: a family flows into a lineage at a constant
rate, and the trait-modulated *loss* then selectively retains it, kept where the trait favours it and
purged where it does not. So the **net** gene content of a lineage tracks its trait even though the
influx itself is trait-blind. That differential retention is what writes the trait–gene association
into the profiles.

::: note
Coupling through *loss* is the mechanism that produces a clean, datable signal (it is the validated
device of the Potts coupling model). `effect_gain` optionally scales a lineage's transfer activity by
$\exp(\text{effect\_gain} \cdot s)$, but it is a donor-side effect and is **off by default** — the
retention channel already makes net gene content track the trait.
:::

### Choosing the responsive families

`TraitGeneCoupling.build(n_families, responsive, ...)` populates the weight vector. The `responsive`
selector is the flexible part:

```python
z.TraitGeneCoupling.build(50, 8)                 # 8 families, chosen at random
z.TraitGeneCoupling.build(50, 0.3)               # a random 30% of the panel
z.TraitGeneCoupling.build(50, ["F3", "F7", 12])  # exactly these families (id or index)
z.TraitGeneCoupling.build(50, 10, signed=True)   # half favoured by a high trait value,
                                                 # half by a low one
```

`weight` sets each responsive family's magnitude; `signed=True` randomises its sign so some families
co-occur with a high trait value and others with a low one. `effect_loss` is the overall coupling
strength ($0$ recovers plain, uncoupled gene-family evolution). The remaining rate parameters —
`base_loss`, `transfer`, `duplication`, `origination` — are the panel's base DTL rates.

### The trait as a covariate in time

The trait value varies *along* each branch, and the simulation follows it exactly:

- **Discrete traits** (`Mk`, threshold) contribute their exact stochastic character map — the
  per-branch (state, duration) segments — so a mid-branch state change is honoured to the instant it
  happens (it becomes a rate-refresh point in the Gillespie loop).
- **Continuous traits** (`BrownianMotion`, `OrnsteinUhlenbeck`) are sub-segmented into `trait_steps`
  pieces per branch (default 16), with the value interpolated between the node endpoints and held
  constant across each piece.

For a binary trait it is usually best to **center** the two states around zero
(`state_values=[-1.0, 1.0]`), so the trait pushes a responsive family's retention *up* in one state
and *down* in the other — a symmetric, two-sided coupling — rather than only lowering loss in the
"on" state:

```python
coupling = z.TraitGeneCoupling.build(40, 0.3, weight=1.0, effect_loss=3.0,
                                     base_loss=0.5, transfer=1.0,
                                     state_values=[-1.0, 1.0], seed=1)
```

`simulate_trait_linked_genomes` accepts either a trait **model** (evolved for you) or an
already-simulated `TraitResult`, so you can inspect or reuse the exact trait the genes were
conditioned on. It returns a `TraitLinkedResult` whose `.profiles` is the $N \times$ extant-species
panel matrix, `.trait` the `TraitResult` the genes were conditioned on, `.coupling` the weights and
effect sizes, and `.genomes()` a promotion to a standard `Genomes` for gene trees and reconciliations.

### From the command line

`zombi2 coevolve-genetrait` runs the whole thing on a species tree you provide. It simulates the
trait (`--trait-model`, reusing every `zombi2 trait` model), builds the coupling (`--panel`,
`--responsive`, `--weight`, `--effect-loss`), and writes the gene-family output plus the trait and a
coupling manifest:

```bash
T=species_tree.nwk

# a binary aerobic/anaerobic trait; 30% of a 40-family panel respond to it
zombi2 coevolve-genetrait -t $T \
    --trait-model mk --states 2 --rate 0.3 --trait-center \
    --panel 40 --responsive 0.3 --weight 1 --effect-loss 3 \
    --loss 0.4 --trans 1.0 --write all --seed 7 -o out/
```

Besides the usual gene-family files (chosen with `--write`, exactly as in `genomes`), it always writes
`traits.tsv` / `trait_tree.nwk` (the trait at every node) and `coupling.tsv` (the per-family coupling
weights and effect sizes, so the exact trait–gene linkage that generated the profiles is on record for
downstream inference). `--responsive` takes a count, a fraction, an id/index list, or `@file` of ids;
`--signed` randomises the weight signs; `--trait-center` centers a discrete trait's states;
`--trait-steps K` sets the within-branch resolution for a continuous trait; `--trait-file` reuses a
precomputed trait; and `--effect-gain` turns on the optional donor-side coupling.

### What it recovers

Inject a strong coupling and the trait shows up in the profiles: responsive families are present where
the trait favours them and absent where it does not, while inert families do not distinguish the
states. Concretely, with a two-clade trait (half the tips aerobic, half anaerobic) and
$\text{effect\_loss} = 3$, responsive families sit at roughly 0.7 prevalence in the aerobic tips and
0.1 in the anaerobic ones, whereas inert families are indistinguishable between the two — the signal
is entirely in the responsive panel, which is exactly what a downstream inference should be able to
pick out.

::: tip
Keep `base_loss` moderate relative to `transfer` so the *inert* families persist as a control. With an
over-large `base_loss`, an unprotected family — having only the field-blind influx to regain it — is
lost tree-wide and the inert rows go all-zero.
:::

## State-dependent diversification

The reverse direction of coupling — a trait shaping the **species tree** — is the `traits:species`
edge. A discrete or continuous trait drives speciation and extinction, and the tree is grown *jointly*
with the trait, so the command takes **no `-t` tree** (it produces one) and a stopping condition
instead:

```bash
# BiSSE: state 1 speciates faster, so it comes to dominate the standing tips
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/
```

This writes `species_tree.nwk` (the tree the trait's rates shaped), `traits.tsv` (every node — tips
*and* ancestral states), and `trait_tree.nwk`. The `--sse-model` selector chooses the flavour: `bisse`
is the **binary** state-dependent birth–death process [@maddison2007bisse]; `musse` the **k-state**
variant [@fitzjohn2012diversitree]; and `quasse` the **continuous-trait** variant [@fitzjohn2010quasse].

From Python the driver is `simulate_sse`:

```python
res = z.simulate_sse(z.BiSSE(1, 3, 0.2, 0.2, 0.1, 0.1), n_tips=200, seed=1)
res.tree              # complete tree (extinct lineages kept; z.prune() for the reconstructed one)
res.labeled_values()  # the trait at the extant tips
```

`z.BiSSE`, `z.MuSSE`, `z.QuaSSE` and `z.HiSSE` (hidden-state SSE) are all on the public API. `HiSSE`
is the honest null: it lets diversification-rate variation be driven by an unobserved state rather
than the focal trait [@beaulieu2016hisse]. The trait can also be made to jump *at* speciation
(cladogenetic evolution) with the `species:traits` edge; turning on **both** `traits:species` and
`species:traits` gives the full ClaSSE feedback, where the trait shapes the tree and is kicked by its
branching.

::: warning
State-dependent diversification models are notorious for high false-positive rates in *inference*: a
neutral trait can look strongly diversification-linked [@maddison2015unreplicated]. ZOMBI2 is a
**forward generator**, so it is not itself fooled — but the whole point of generating these scenarios
is to test downstream inference, so it is worth generating the null (a zero-effect coupling) as a
control, exactly as the trait-linked model keeps inert families as controls.
:::
