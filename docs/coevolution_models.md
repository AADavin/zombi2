# Coevolution: coupled species, traits and gene families — a design

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

This note fixes the model, the command-line shape, and a phased build order. Only one edge of the
model ships today (`traits:genes`, as [`coevolve-genetrait`](guide/trait-linked-genomes.md));
everything else here is a design for what `coevolve` becomes.

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
| `genes:species` | gene content sets diversification | gene-content-dependent diversification | **output** (forward) | proposed — Phase 3 (the merged engine) |
| `species:traits` | trait jumps *at* speciation | cladogenetic / speciational trait evolution | input (given tree) | **shipped** — `coevolve --couple species:traits` (both arrows = **ClaSSE**) |
| `species:genes` | gene gain/loss bursts at speciation | cladogenetic genome upheaval | input | proposed |
| `traits:genes` | trait sets gene loss/gain | **trait-linked gene families** | input | **shipped** — [`coevolve-genetrait`](guide/trait-linked-genomes.md) |
| `genes:traits` | gene presence enables a trait shift | gene-conditioned trait | input | proposed |

## The one rule: complexity = arrows *into* S

Everything about the difficulty (and about whether the command takes a tree or makes one) follows
from a single question: **does any active edge point into S?**

- **No edge into S** — the tree is fixed. It is read from `-t/--tree` (or a prior `species` run),
  and every coupling is an *overlay* on a frozen tree. This is a pipeline, not a joint simulation;
  even several axes at once stay a sequence of overlays.
- **An edge into S** (`traits:species` and/or `genes:species`) — the tree topology now *depends*
  on the coupled state, so it cannot be drawn first: **the tree is an output**. These runs are
  **forward-only** (they generate the complete tree, extinct lineages included) and take no
  `-t`. This is the same forward machinery as [`species --model forward`](species_tree_models.md).

!!! warning "`genes:species` is the one real complexity cliff"
    An arrow from **G into S** means the genome content must be known *as the tree grows*, so the
    species birth–death events and the gene-family DTL events interleave in **one** event stream.
    That merges the two Gillespie engines into a single loop and forfeits the fast
    [Rust genome engine](guide/rust-engine.md) (which assumes a frozen tree). It is expensive
    whether or not a trait is also present — so paying for `genes:species` once is what unlocks
    genuine three-way simulation. Every *other* edge is comparatively cheap.

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
import zombi2 as z
res = z.simulate_sse(z.BiSSE(1, 3, 0.2, 0.2, 0.1, 0.1), n_tips=200, seed=1)
res.tree                 # complete tree (extinct lineages kept; z.prune() for the reconstructed one)
res.labeled_values()     # the trait at the extant tips
```

`z.BiSSE` / `z.MuSSE` / `z.QuaSSE` and `z.HiSSE` (hidden-state SSE, the honest null) are all on the
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
`z.Cladogenesis(shift=…, jump_sigma2=…)`, accepted by both `z.simulate_sse(..., cladogenesis=…)`
(ClaSSE) and `z.simulate_traits(tree, model, cladogenesis=…)` (`species:traits` on a fixed tree).

The into-species `genes:species` edge (and the full three-way `--all`) remain on the roadmap below.

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

- **Phase 0 — the umbrella.** Add `coevolve` with the `--couple driver:target` parser; fold the
  existing [`coevolve-genetrait`](guide/trait-linked-genomes.md) in as `--couple traits:genes`
  (keep `coevolve-genetrait` as a deprecated alias). Ships immediately; no engine work.
- **Phase 1 — `traits:species` (SSE). ✅ done.** The forward joint tree+trait engine is in
  [`zombi2/sse.py`](https://github.com/AADavin/zombi2/blob/main/zombi2/sse.py) — `BiSSE` (binary),
  `MuSSE` (k-state), `QuaSSE` (continuous) and `HiSSE` (hidden-state), driven by `simulate_sse` and
  exposed as `coevolve --couple traits:species`. Next: fold the speciation→trait change kernel into
  the same loop for Phase 2.
- **Phase 2 — `species:traits` and full ClaSSE. ✅ done.** The `Cladogenesis` kernel
  (`zombi2/traits.py`) jumps the trait at speciation; it feeds both `simulate_traits` (the
  `species:traits` edge on a given tree) and the forward `simulate_sse` loop, so
  `--couple traits:species --couple species:traits` is the complete `traits↔species` ClaSSE feedback.
- **Phase 3 — `genes:species` (the merged loop).** The one hard build: interleave species and
  gene-family events in a single forward stream. This is also what makes genuine three-way
  (`--all`) real.
- **Phase 4 — `--all` and the remaining overlay edges** (`species:genes`, `genes:traits`) as
  additive rate-functions/kernels, validated against the single-edge cases.

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
