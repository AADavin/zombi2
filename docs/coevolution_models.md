# Coevolution: coupled species, traits and gene families

ZOMBI2 simulates in a **pipeline**: build a species tree, then evolve a trait along it, then
evolve gene families along it (optionally conditioned on the trait). That works because today's
couplings are *one-directional* and the joint distribution **factorises**,

```
P(tree) · P(trait | tree) · P(genes | tree, trait),
```

so each stage can be drawn on the frozen output of the previous one. The moment a trait or a
genome is allowed to feed **back** into the process that generated the tree, this factorisation
breaks and the stages can no longer be run in sequence — the tree, the trait and the gene content
must be grown *together*. The proposed **`coevolve`** mode is exactly for those non-factorising,
feedback scenarios.

This note fixes the model, the command-line shape, and a phased build order. All six single
directed edges now ship under `coevolve --couple driver:target` (see the status table below);
what remains a design is the fully *joint* model that combines several edges with feedback.

## The coupling graph

There are three processes:

- **S** — species diversification (the birth–death process that grows the tree),
- **T** — a phenotypic trait ([BM/OU/EB/Mk/threshold](guide/traits.md)),
- **G** — gene-family content (the [DTL](guide/gene-families.md) process).

A **coupling is a directed edge** `driver → target`: the *driver's* state modulates the
*target's* rates. A coevolution scenario is just a **set of directed edges** on these three nodes.
That single abstraction covers everything from SSE to full three-way feedback, and it makes the
one thing that actually matters — *direction* — explicit rather than implied.

**CLI convention.** One repeatable flag, `--couple driver:target`, where the order reads as the
arrow (driver first, left to right):

```bash
zombi2 coevolve --couple traits:species   ...   # T→S : the trait sets speciation/extinction (SSE)
zombi2 coevolve --couple species:traits   ...   # S→T : speciation drives the trait (cladogenetic)
zombi2 coevolve --couple traits:species --couple species:traits   ...   # both arrows = ClaSSE
zombi2 coevolve --all             ...           # every edge — the fully joint model
```

So `--couple species:traits` and `--couple traits:species` are deliberately **different models**,
and bidirectional coupling is simply *both* edges. `:` (not `->`) keeps it shell-safe.

## Summary — the six directed edges

| Edge (`driver:target`) | Reading | Model | Tree | Status |
|---|---|---|---|---|
| `traits:species` | trait sets speciation/extinction | **SSE** (BiSSE / MuSSE / QuaSSE / HiSSE) | **output** (forward) | **shipped** — `coevolve --couple traits:species` |
| `genes:species` | gene content sets diversification | key-innovation genes + HGT | **output** (forward) | **shipped** — `coevolve --couple genes:species` |
| `species:traits` | trait jumps *at* speciation | cladogenetic / speciational trait evolution | input (given tree) | **shipped** — `coevolve --couple species:traits` (both arrows = **ClaSSE**) |
| `species:genes` | gene gain/loss bursts at speciation | cladogenetic genome upheaval | input (given tree) | **shipped** — `coevolve --couple species:genes` |
| `traits:genes` | trait sets gene loss/gain | **trait-linked gene families** | input (given tree) | **shipped** — `coevolve --couple traits:genes` |
| `genes:traits` | gene presence enables a trait shift | gene-conditioned trait | input (given tree) | **shipped** — `coevolve --couple genes:traits` |

## The one rule: complexity = arrows *into* S

Everything about the difficulty (and about whether the command takes a tree or makes one) follows
from a single question: **does any active edge point into S?**

- **No edge into S** — the tree is fixed. It is read from `-t/--tree` (or a prior `species` run),
  and every coupling is an *overlay* on a frozen tree. This is a pipeline, not a joint simulation;
  even several axes at once stay a sequence of overlays.
- **An edge into S** (`traits:species` and/or `genes:species`) — the tree topology now *depends*
  on the coupled state, so it cannot be drawn first: **the tree is an output**. These runs are
  **forward-only** (they generate the complete tree, extinct lineages included) and take no
  `-t`. This is the same forward machinery as [`species --mode forward`](species_tree_models.md).

!!! note "`genes:species` was the complexity cliff — and how it was tamed"
    An arrow from **G into S** means the genome content must be known *as the tree grows*: the
    tree can't be built first. The v1 (**shipped**) tames this by observing that, under ZOMBI's
    **independent** families, only the handful of families that actually touch diversification
    need to ride in the forward loop. So a small panel of binary **driver** ("key innovation")
    families grows jointly with the tree, and the entire neutral genome — which does not affect
    the tree — is overlaid *afterward* by the ordinary [`genomes`](cli.md) / `simulate_genomes` on
    the finished tree (exact, not an approximation). The genuinely new machinery is the
    **frequency-dependent transfer** (a driver in more genomes is donated more often), which is
    why the loop must know the live population — and is what makes this more than a static SSE.

## Most three-axis scenarios still decompose

"Simulate species, traits and genes simultaneously" does **not** require one monolithic engine.
Because genes and traits usually do not feed back into S, even a scenario that touches all three
axes typically splits into "grow the `S` + its into-S drivers forward, then overlay the rest":

| Scenario | Edges | How it runs |
|---|---|---|
| A trait drives a radiation **and** shapes gene content | `traits:species`, `traits:genes` | forward-grow **(tree, trait)** jointly, then **overlay genes** conditioned on the trait — the grown tree is fed straight into today's [`simulate_trait_linked_genomes`](guide/trait-linked-genomes.md) |
| Full ClaSSE with trait-linked genes | `traits:species`, `species:traits`, `traits:genes` | same: grow the `S`+trait core forward, overlay genes downstream |
| Genome content itself drives diversification | `genes:species` (± anything) | **the merged loop** — tree and genome grow together in one interleaved event stream |

So true, all-in-one simultaneity is required in exactly one case: an arrow from **G into S**.
Otherwise "S + T + G together" is a forward core plus overlays, reusing pieces that already exist.

## Using it today (`traits:species`, `species:traits`, ClaSSE)

The first into-species edge is **shipped**. A discrete or continuous trait drives
speciation/extinction and the tree is grown *jointly* with it — so the command takes **no `-t`
tree** (it produces one) and a stopping condition instead:

```bash
# BiSSE: state 1 speciates faster, so it comes to dominate the standing tips
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/
```

This writes `species_tree.nwk` (the tree the trait's rates shaped), `traits.tsv` (every node —
tips *and* ancestral states), and `trait_tree.nwk`. `--sse-model musse` is the k-state variant
(`--birth`/`--death` vectors + a `--q-matrix` file); `--sse-model quasse` is the continuous-trait
variant (sigmoidal speciation via `--spec-low/high/center/slope` + Brownian `--diffusion`). From
Python the driver is `simulate_sse`:

```python
from zombi2.coevolve import simulate_sse, BiSSE
from zombi2.species import prune
res = simulate_sse(BiSSE(1, 3, 0.2, 0.2, 0.1, 0.1), n_tips=200, seed=1)
res.tree                 # complete tree (extinct lineages kept; prune() for the reconstructed one)
res.labeled_values()     # the trait at the extant tips
```

`BiSSE` / `MuSSE` / `QuaSSE` and `HiSSE` (hidden-state SSE, the honest null) are all on the
public API.

### The reverse arrow — `species:traits` and ClaSSE

The second arrow is also shipped: **`species:traits`** makes the trait jump *at* each speciation
(cladogenetic / speciational evolution). On its own it has no arrow into S, so it runs on a
**given** tree (`-t`), with a `Cladogenesis` kernel layered on an ordinary anagenetic model:

```bash
# a purely speciational binary trait on an existing tree (no within-branch change: --q01/--q10 0)
zombi2 coevolve --couple species:traits -t out/species_tree.nwk \
    --sse-model bisse --q01 0 --q10 0 --clado-shift 0.4 --seed 2 -o clado/
```

Turn on **both** arrows and you get the full **ClaSSE** feedback — the trait shapes the tree *and*
is kicked by its branching:

```bash
zombi2 coevolve --couple traits:species --couple species:traits \
    --lambda0 1 --lambda1 3 --q01 0.05 --q10 0.05 --clado-shift 0.3 \
    --tips 200 --seed 3 -o classe/
```

`--clado-shift` is the per-daughter state-hop probability (discrete traits); `--clado-jump` is the
Gaussian jump variance for a continuous (`quasse`) trait. In Python the kernel is
`Cladogenesis(shift=…, jump_sigma2=…)` (`from zombi2.traits import Cladogenesis`), accepted by both
`simulate_sse(..., cladogenesis=…)` (ClaSSE) and `simulate_traits(tree, model, cladogenesis=…)`
(`species:traits` on a fixed tree).

### Gene content drives the tree — `genes:species`

The second into-species edge is shipped too: a small panel of binary **driver** ("key innovation")
gene families whose presence sets each lineage's speciation/extinction rate. Drivers are gained de
novo (origination) and — the interesting part — by **transfer**, which is *frequency-dependent* (a
driver in more live genomes spreads faster), so the tree and the gene content grow together. The
neutral bulk genome is overlaid afterward on the grown tree:

```bash
# 2 drivers; one seeded at the root with a strong speciation effect, spread by HGT
zombi2 coevolve --couple genes:species --drivers 2 --root-drivers 1 \
    --lambda0 1 --mu0 0.2 --driver-speciation 1.2 --driver-transfer 0.8 --driver-loss 0.3 \
    --tips 200 --seed 1 -o keygene/

# overlay the neutral genome on the resulting tree (the factorization, made explicit)
zombi2 genomes -t keygene/species_tree.nwk --trans 1 --loss 0.5 --write profiles trees -o keygene/
```

`coevolve` writes `species_tree.nwk`, `drivers.tsv` (per-node driver presence) and
`drivers_manifest.tsv` (the effect sizes β and rates). `--driver-speciation`/`--driver-extinction`
are the per-driver log-rate effects; `--driver-loss`/`--driver-origination`/`--driver-transfer` the
gene dynamics. In Python: `simulate_gene_diversification(GeneDiversification(…), n_tips=…)`
(`from zombi2.coevolve import simulate_gene_diversification, GeneDiversification`).

### Speciation drives the genome — `species:genes`

The reverse of `genes:species`: here gene content does **not** affect diversification, so this is an
overlay on a **given** tree — the genomic twin of `species:traits`. A genome is evolved down the tree
with a *cladogenetic burst* of gene loss and gain at each speciation (founder-effect upheaval), on
top of the usual gradual (anagenetic) change — **punctuational genome evolution**:

```bash
# purely punctuational genomes: change ONLY at speciations (no anagenetic gene-loss/-origination)
zombi2 coevolve --couple species:genes -t species_tree.nwk \
    --genome-size 30 --clado-gene-loss 0.15 --clado-gene-gain 3 --seed 2 -o punctgenome/
```

`--clado-gene-loss` is the per-family drop probability at each speciation and `--clado-gene-gain` the
mean number of new families gained (Poisson); `--gene-loss`/`--gene-origination` add gradual
along-branch change (both 0 = pure punctuation). It writes `Profiles.tsv`/`Presence.tsv` (families ×
extant tips) and `genome_sizes.tsv`. In Python:
`simulate_cladogenetic_genome(tree, CladogeneticGenome(…))`
(`from zombi2.coevolve import simulate_cladogenetic_genome, CladogeneticGenome`). The signature of the model is that
**sister tips differ** — change is injected at their split, not spread along the branches.

### Gene content shapes a trait — `genes:traits`

The reverse of `traits:genes`: here gene content conditions a **trait**. A binary *modifier* gene comes
and goes along the tree (gain/loss), and its presence sets a continuous trait's **OU optimum** — so a
lineage that acquires the gene is pulled toward a new adaptive peak (`theta_present`), and one that
loses it drifts back to `theta_absent`. "Gene presence enables a trait shift." Also an overlay on a
given tree:

```bash
# a modifier gene that unlocks a phenotypic optimum at 5 (vs 0 without it)
zombi2 coevolve --couple genes:traits -t species_tree.nwk \
    --modifier-gain 0.6 --modifier-loss 0.6 --theta-absent 0 --theta-present 5 \
    --trait-alpha 2.5 --trait-sigma2 0.4 --seed 2 -o genetrait/
```

It writes `traits.tsv` (per node: modifier presence 0/1 and the trait value) and `trait_tree.nwk`.
`--trait-alpha` is the OU mean-reversion (0 = Brownian); the modifier's own dynamics are
`--modifier-gain`/`--modifier-loss`/`--root-modifier`. In Python:
`simulate_gene_conditioned_trait(tree, GeneConditionedTrait(…))`
(`from zombi2.coevolve import simulate_gene_conditioned_trait, GeneConditionedTrait`). The signal: **tips carrying the
modifier sit near `theta_present`**, those without near `theta_absent`.

Only the full three-way `--all` remains on the roadmap below.

## CLI reference

Every edge is one `zombi2 coevolve --couple driver:target` invocation; the examples throughout
this page are the authoritative per-edge usage, and `zombi2 coevolve -h` prints the full flag set.
The one edge whose flags differ most from `genomes`/`trait` is **`traits:genes`** (trait-conditioned
gene families), documented in full here.

### `traits:genes` — trait-conditioned gene families

The **`traits:genes`** edge links the two halves of the toolkit: it evolves a phenotypic trait
along the tree, then evolves a **panel** of gene families whose loss and gain **depend on the local
trait value**, so the resulting profile carries a known, trait-linked signal (the forward generator
behind reading gene content as a record of a trait's history — e.g. dating the tree from the Great
Oxidation Event). It simulates the trait with any [`trait`](cli.md#trait-a-phenotypic-trait) model
(`--trait-model`), builds the coupling, and writes the gene-family output alongside the trait and a
coupling manifest. (This was the standalone `coevolve-genetrait` command before it was folded into
`coevolve`.)

```bash
T=out/species_tree.nwk

# a binary aerobic(1)/anaerobic(0) trait; 30% of a 40-family panel respond to it
zombi2 coevolve --couple traits:genes -t $T \
    --trait-model mk --states 2 --rate 0.3 --trait-center \
    --panel 40 --responsive 0.3 --weight 1 --effect-loss 3 \
    --loss 0.4 --trans 1.0 --write all --seed 7 -o out/
```

A *responsive* family is retained where the trait favours it (loss scaled by
`exp(-effect_loss · weight · trait)`) and purged where it does not; gain is field-blind
horizontal transfer, so the **net** gene content of a lineage tracks its trait. `--responsive`
chooses which families respond — a count (`8`), a fraction (`0.3`), an id/index list
(`F3,F7,12`), or `@file` of ids — and `--signed` randomises the weight signs so some families
co-occur with a high trait value and others with a low one. `--trait-center` centers a discrete
trait's states (recommended for a binary character, giving a symmetric two-sided coupling), and
`--trait-steps K` sets the within-branch resolution for a continuous trait (discrete traits use
their exact stochastic map). `--effect-gain` optionally scales a lineage's transfer activity by
the trait too (off by default).

It writes the gene-family files selected by `--write` (as [`genomes`](cli.md#choosing-the-output-and-the-rust-engine)),
and always adds **`traits.tsv`** / **`trait_tree.nwk`** (the trait at every node) and
**`coupling.tsv`** (the per-family weights and effect sizes — the trait↔gene linkage on record
for downstream inference). Reuse a precomputed trait instead of simulating one with
`--trait-file traits.tsv` (a `node`/`value` table over **every** node — tips and ancestors —
with numeric values, as `zombi2 trait` writes). See
[Trait-linked gene families](guide/trait-linked-genomes.md) for the model.

| Option | Meaning |
| --- | --- |
| `--couple traits:genes` | select the trait-conditioned-genes edge (required) |
| `--tree` / `-t` | input species tree in Newick format (required) |
| `--trait-model {bm,ou,eb,mk,threshold}` | trait to evolve then couple to gene families (default `bm`); its parameters are the [`trait`](cli.md#trait-a-phenotypic-trait) flags (`--sigma2`, `--alpha`/`--theta`, `--rate`, `--states`/`--ordered`/`--q-matrix`, `--thresholds`, …) |
| `--trait-file TSV` | reuse a precomputed trait instead — a numeric `node`/`value` table over **every** node (as `zombi2 trait` writes); overrides `--trait-model` |
| `--trait-center` | [discrete] center the state values around their mean (two-sided coupling; recommended for a binary trait) |
| `--trait-steps K` | [continuous] within-branch resolution — sub-segment each branch into K pieces (default `16`; ignored for discrete traits) |
| `--panel` | number of gene families in the panel (default `50`) |
| `--loss` `--trans` `--dup` `--orig` | panel base rates — baseline per-copy loss (default `0.5`), transfer/HGT gain (default `1.0`), duplication, origination |
| `--responsive SPEC` | which families respond: a count, a fraction (e.g. `0.3`), an id/index list (`F3,F7,12`), or `@FILE` (default `0.3`) |
| `--weight` / `--signed` | coupling weight of each responsive family (default `1.0`) / randomise its sign |
| `--effect-loss` | retention coupling strength: loss scales by `exp(-effect_loss · weight · trait)` (default `2.0`; `0` = uncoupled) |
| `--effect-gain` | optional donor-side HGT-activity coupling: transfer scales by `exp(effect_gain · trait)` (default `0`) |
| `--write {profiles,trace,trees,events,transfers,summary,all}` | which gene-family files to write (default `profiles trees`); `traits.tsv` / `trait_tree.nwk` / `coupling.tsv` are always written too |
| `--sparse` / `--annotate-species` | sparse profile table / label internal gene-tree nodes (as in `genomes`) |
| `--seed` / `-o` / `--out` | RNG seed / output directory (required) |

For the into-species edges (`traits:species`, `genes:species`) and the other overlay edges, the
per-edge examples above list their flags; `zombi2 coevolve -h` is the authoritative list.

## The engine: one generic per-lineage state

The investment that unlocks the whole into-S family is a single generalisation of the forward
birth–death loop ([`species_forward.py`](species_tree_models.md#forward-simulation-implemented)):
today a growing lineage is an *unlabelled topology*; the coevolve engine gives each lineage a
**state bag** and lets the birth–death rates be a function of it.

```
lineage.state   # an open container of the active coupled processes: a trait value, a genome, or both
rate(lineage)   # speciation / extinction read off lineage.state
on_speciation   # daughters inherit (copy) the state, optionally with a change kernel applied
```

Crucially the loop does not care *what* the state is:

- put a **trait** in the bag and read it in `rate()` → `traits:species` (SSE);
- put a **genome** in the bag → `genes:species`;
- put **both** → three-way.

And the reverse arrows are nearly free in the same loop: applying a change kernel to the daughters
in `on_speciation` gives `species:traits` (cladogenetic trait jumps) and `species:genes` — so the
full `traits:species` **+** `species:traits` feedback (ClaSSE) needs only the trait→rate function
*and* the speciation→trait kernel, both small and additive. `--all` is then not a separate build
but the state where every edge's rate-function/kernel is wired on; it becomes a validation
milestone once the individual edges each work.

## Phased build order

- **Phase 0 — the umbrella. ✅ done.** `coevolve` with the `--couple driver:target` parser; the
  former standalone `coevolve-genetrait` command is folded in as `--couple traits:genes`
  (the standalone command was removed — a clean break, no alias).
- **Phase 1 — `traits:species` (SSE). ✅ done.** The forward joint tree+trait engine is in
  [`zombi2/sse.py`](https://github.com/AADavin/zombi2/blob/main/zombi2/sse.py) — `BiSSE` (binary),
  `MuSSE` (k-state), `QuaSSE` (continuous) and `HiSSE` (hidden-state), driven by `simulate_sse` and
  exposed as `coevolve --couple traits:species`. Next: fold the speciation→trait change kernel into
  the same loop for Phase 2.
- **Phase 2 — `species:traits` and full ClaSSE. ✅ done.** The `Cladogenesis` kernel
  (`zombi2/traits.py`) jumps the trait at speciation; it feeds both `simulate_traits` (the
  `species:traits` edge on a given tree) and the forward `simulate_sse` loop, so
  `--couple traits:species --couple species:traits` is the complete `traits↔species` ClaSSE feedback.
- **Phase 3 — `genes:species` (key innovations + HGT). ✅ done.** A forward joint loop
  ([`zombi2/gene_diversification.py`](https://github.com/AADavin/zombi2/blob/main/zombi2/gene_diversification.py))
  grows the tree together with a panel of binary driver families whose presence sets λ/μ; gain is by
  origination and **frequency-dependent transfer**. The neutral genome is overlaid afterward with
  `genomes` (exact under independent families), so the merged loop stays small. `GeneDiversification`
  / `simulate_gene_diversification`, exposed as `coevolve --couple genes:species`. **v1 scope:** binary
  drivers, driver *presence* profiles (not yet full driver gene trees), the edge runs on its own.
- **Phase 4a — `species:genes` (cladogenetic genome). ✅ done.** An overlay on a given tree
  ([`zombi2/cladogenetic_genome.py`](https://github.com/AADavin/zombi2/blob/main/zombi2/cladogenetic_genome.py)):
  a genome is evolved down the tree with a founder-effect burst of gene loss/gain at each speciation
  (`CladogeneticGenome` / `simulate_cladogenetic_genome`, `coevolve --couple species:genes`). v1 scope:
  presence/absence genome, gain by origination (no HGT yet), runs on its own.
- **Phase 4b — `genes:traits` (gene-conditioned trait). ✅ done.** An overlay on a given tree
  ([`zombi2/gene_conditioned_trait.py`](https://github.com/AADavin/zombi2/blob/main/zombi2/gene_conditioned_trait.py)):
  a binary modifier gene (a two-state Markov chain along the tree) switches a continuous trait's OU
  optimum, so acquiring the gene pulls the trait to a new peak (`GeneConditionedTrait` /
  `simulate_gene_conditioned_trait`, `coevolve --couple genes:traits`). **All six directed edges now
  exist individually.**
- **Phase 5 — the pairwise *joint* (both-arrow) models. ✅ done for all three pairs.** Each node-pair's
  two edges can be switched on together (`--couple A:B --couple B:A`), giving a bidirectional model
  where the *same* coupled object drives both directions — as ClaSSE (`traits:species` +
  `species:traits`) already did. Now also: **co-diversification** (`genes:species` + `species:genes`) —
  the driver panel both sets λ/μ *and* is reshuffled by a cladogenetic burst at each speciation, so
  speciation itself seeds rate heterogeneity (one arrow into S → tree is an output;
  `cladogenetic_loss`/`cladogenetic_gain` on `GeneDiversification`, `simulate_co_diversification`); and
  **trait–gene feedback** (`traits:genes` + `genes:traits`) — a trait and a coupled panel modulate each
  other, integrated jointly along each branch, so the tips end up correlated with no single edge imposed
  (overlay; [`zombi2/trait_gene_feedback.py`](https://github.com/AADavin/zombi2/blob/main/zombi2/trait_gene_feedback.py),
  `TraitGeneFeedback` / `simulate_trait_gene_feedback`). Each contains its two single edges as limits.
- **Phase 4 — `--all`.** The fully joint model: all six edges active at once (every pair
  bidirectional = maximal mutual feedback). There is no single causal direction — forward time
  resolves the mutual dependence (current states set next-event rates). Needs the generic
  per-lineage-state forward loop carrying both trait and genome, plus stability guards on the
  feedback loops (esp. trait↔gene). Best expressed by composing edges with `--couple`; `--all`
  would be sugar for "all six".

## Caveats

!!! note "SSE draws spurious associations easily"
    State-dependent diversification models are notorious for high false-positive rates in
    *inference*: a neutral trait can look strongly diversification-linked (Rabosky & Goldberg
    2015; Maddison & FitzJohn 2015). ZOMBI2 is a **forward generator**, so it is not itself
    fooled — but the whole point of generating these scenarios is to test downstream inference, so
    the simulator should make it easy to produce the null (`--couple` absent, or a zero-effect
    coupling) as a control, exactly as the `traits:genes` model keeps inert families as controls.

Fitting or even interpreting the fully-coupled `--all` model is hard (many entangled knobs); as a
*simulator* that is fine, but three-way feedback is best treated as a stress-test/showcase rather
than a routine analysis mode. Performance-wise, only the into-S edges leave the fast paths — every
overlay edge keeps the existing engines.

## Key references

- Maddison, Midford & Otto (2007), *Syst. Biol.* — BiSSE (binary state-dependent speciation/extinction).
- FitzJohn (2010), *Syst. Biol.* — QuaSSE (quantitative-trait SSE); FitzJohn (2012), *Methods Ecol. Evol.* — diversitree, MuSSE.
- Beaulieu & O'Meara (2016), *Syst. Biol.* — HiSSE (hidden-state SSE, the honest null).
- Goldberg & Igić (2012), *Evolution* — ClaSSE (cladogenetic state change + SSE): the model behind the reverse `species:traits` arrow.
- Bokma (2008), *J. Evol. Biol.*; Pagel (1999), *Nature* — speciational / punctuational trait evolution (change concentrated at branching).
- Rabosky & Goldberg (2015), *Syst. Biol.*; Maddison & FitzJohn (2015), *Syst. Biol.* — the false-positive caveat for SSE inference.
- Maliet, Hartig & Morlon (2019), *Nat. Ecol. Evol.* — ClaDS: the per-lineage-rate cousin already shipped, and a template for carrying per-lineage state forward.
