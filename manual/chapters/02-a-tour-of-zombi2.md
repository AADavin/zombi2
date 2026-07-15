# A tour of ZOMBI2

Before the tutorial chapters, this one lays out the ideas the whole book leans on: the **four
levels** ZOMBI2 simulates, the **two ways** it composes them, the single shape **every rate**
takes, and the **vocabulary** that recurs throughout. None of it is deep — it is the map you keep
in your head while reading the rest. Later chapters refer back here rather than re-deriving it.

## The four levels of evolution

Evolution is one process, but it leaves its trace at several levels at once: a lineage diversifies
into a species tree; along that tree, genomes gain and lose gene families; the organisms carry
traits that drift and adapt; and inside every gene, molecular sequences accumulate substitutions.
Reality draws no line between these levels. We do — because a tractable model usually captures only
one of them at a time.

ZOMBI2 is organised around that separation. It simulates four levels:

- **species** — the dated tree of lineages (a birth–death process);
- **genomes** — the gene families that live along that tree, gained and lost by duplication,
  transfer, loss and origination;
- **traits** — phenotypes evolving on the tree (Brownian motion, Ornstein–Uhlenbeck, Mk,
  threshold, biogeographic ranges);
- **sequences** — the nucleotides or amino acids inside each gene.

You always choose *which* levels to simulate — you need not run them all. A study of gene content
may never descend to sequences; a study of diversification may never leave the species level.

## Two ways to compose them: pipeline and coevolution

Once you run more than one level, the question is how they fit together. There are two answers.

![The four levels as a diamond — the same diamond the [coevolution](#coevolution-the-framework) part later fills with couplings. In the default *pipeline*, each level is simulated *along* (conditioned on) the one above: traits and gene content along the species tree, sequences along the gene trees. Coupling the levels instead of layering them is coevolution (Part VI). You always choose which levels to simulate — you need not run them all.](figures/levels_diamond.pdf){width=80%}

The usual way is **hierarchical**: simulate one level, then the next *conditioned on* it — a
species tree, then gene families along its branches, then sequences along the resulting gene trees.
Each level treats the one above as a fixed backbone. This keeps the levels independent and easy to
reason about, and lets a single species tree seed many genome, trait, or sequence runs, so one
backbone underlies a whole benchmark. It works because the joint distribution *factorises* — for
instance `P(tree)·P(trait | tree)·P(genes | tree, trait)` — so each stage runs on the frozen output
of the previous one.

Sometimes the levels are *not* independent: a trait may change the rate at which genes are gained
or lost; gene content may decide which lineages survive and diversify. When influence runs *between*
levels the factorisation breaks and a one-directional pipeline is no longer enough — the levels must
be grown **jointly**. ZOMBI2's [`coevolve`](#coevolution-the-framework) mode does this, coupling
species, genomes and traits along directed links you choose. A coupling is a directed edge
`driver → target`: the driver's state modulates the target's rates. **[Coevolution](#coevolution-the-framework)** (Part V)
is the whole story; the point here is just that it is the alternative to the pipeline, reached only
when a level feeds back on what produced it.

## How rates work: how many clocks, how fast

Every process in ZOMBI2 — a lineage speciating, a gene duplicating, a nucleotide mutating — is
driven by a **rate**, and they all share one shape. Meeting it once makes every model in the book
read the same way.

Think of every pending event as a **clock** that ticks at some rate; when it ticks, the event
fires. A whole model is then two answers:

1. **How many clocks?** One per lineage? one per gene copy? one shared by the whole tree? This is
   the *opportunity* — how many independent chances the event has right now.
2. **How fast does each tick?** The *base rate*, possibly rescaled by context — time, crowding,
   which lineage or family it is.

Put together, the instantaneous propensity of an event is

$$\text{propensity} \;=\; \underbrace{\text{base}}_{\text{how fast}} \;\times\; \underbrace{\text{opportunity}}_{\text{how many clocks}} \;\times\; \underbrace{\text{modifiers}}_{\text{context}}.$$

The three terms have distinct jobs and distinct units. The **base** carries the units — it is a
rate, $\text{time}^{-1}$. The **opportunity** is a *count* (dimensionless). The **modifiers** are
*multipliers* (also dimensionless).

![How many clocks, how fast. The opportunity — one clock per lineage, one per gene copy, or one shared by the whole process — sets how many independent chances an event has right now. A count that *tracks the growing quantity* (per lineage, per copy) compounds into exponential growth; a single shared clock keeps the total rate constant, so growth is linear. Both curves are ZOMBI2 runs at the *same* per-clock rate (`BirthDeath` vs `SharedBirthDeath`, mean over 40 seeds) — only the number of clocks differs.](figures/rate_clocks.pdf){width=100%}

### The opportunity — "per what?"

The opportunity answers *per what?* The choices form a short ladder, finest unit to coarsest, plus
one that stands apart:

| Opportunity | One clock per… | Rides events like |
|---|---|---|
| **`site`** | nucleotide | substitution, insertion, deletion |
| **`copy`** | gene copy | duplication, loss, transfer |
| **`lineage`** | genome / branch | speciation, extinction; genome-level loss |
| **`shared`** | *nothing — one clock for the whole thing* | a shared diversification budget |

The per-unit rungs nest — $\texttt{site} \subset \texttt{copy} \subset \texttt{lineage}$ — and
`shared` is the odd one out: not one clock per unit, but a single clock for the entire process.

**This choice, not the units, decides exponential versus linear growth** — because a count that
*tracks the growing quantity* compounds, and one that does not, cannot. You can see it directly in
the species models. `BirthDeath` puts one clock on each lineage (`lineage`); `SharedBirthDeath` puts
one clock on the whole tree (`shared`):

```python
import zombi2 as z

# one clock per lineage: births compound → exponential
per_lineage = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), age=8.0,
                                      direction="forward", seed=3)

# one shared clock: constant total birth rate → linear
shared = z.simulate_species_tree(z.SharedBirthDeath(1.0, 0.2), age=8.0,
                                 direction="forward", seed=3)
```

With the same rates and seed, the per-lineage tree reaches thousands of tips while the shared tree
reaches a few dozen: same *speed* per clock, a different *number of clocks*. The genome level offers
the same choice one rung down — gene-family rates counted **per gene copy** (`--rate-per copy`, the
default, so families grow exponentially) or **per lineage** (`--rate-per lineage`, size-independent,
so content grows linearly). "Per lineage" is exactly the size-independent measure speciation uses
one level up: a lineage carries one genome, so *per genome* **is** *per lineage*.

### The modifiers — context

The opportunity decides *how many* clocks; **modifiers** decide how fast a given clock ticks
relative to the base, by context (some gene families turn over faster than others; some branches run
hot — the relaxed-clock idea, which is also what a [molecular clock](#molecular-clocks) does to
substitution rates). Modifiers are dimensionless multipliers: they change *how fast*, never *how
many chances*.

The payoff is that the same shape holds everywhere:

| Level | Opportunity it uses | "How fast" set by |
|---|---|---|
| **species** | `lineage` (or `shared`) | the diversification model |
| **genome content** | `copy` (default) or `lineage` | the DTL rates, ± per-family/per-lineage modifiers |
| **sequences** | `site` | the substitution rate, ± a [molecular clock](#molecular-clocks) |

Read any rate in the book by asking the two questions: *how many clocks (per what?)*, and *how fast
does each tick?*

::: note
Today the opportunity is chosen by the *model you pick* — `BirthDeath` versus `SharedBirthDeath`,
per-copy versus per-lineage rates — or, for genomes, by `--rate-per`. A planned refinement makes it
a single named knob, `per=` (`--rate-per` / `--per`), available the same way at every level, with
those classes kept as friendly shorthands. When it lands, the code will read exactly the way this
section teaches; nothing here changes.
:::

## The ZOMBI2 vocabulary

A few terms recur throughout, some of them specific to ZOMBI2. Collected here for reference:

- **Level** — one of the four processes above (species, genomes, traits, sequences).
- **Opportunity** — *how many clocks* a rate is counted over: `site`, `copy`, `lineage`, or
  `shared` (see above). Per-copy diversifies exponentially, per-lineage/shared linearly.
- **Modifier** — a dimensionless multiplier that rescales a base rate by context (per family, per
  lineage), changing *how fast*, not *how many chances*.
- **Complete vs reconstructed tree** — the *complete* tree keeps every lineage, including those that
  went extinct; the *reconstructed* tree is pruned to the sampled survivors. `direction="forward"`
  gives the complete tree; the default backward sampler gives the reconstructed one, and
  [ghost lineages](#ghost-lineages) can graft the dead back on.
- **Extant / extinct / unsampled** — a tip that reached the present and was sampled is *extant*
  (`is_extant=True`); one that died before the present is *extinct*; one alive but not sampled is
  *unsampled*. Extinct and unsampled tips are collectively the "ghost" lineages.
- **Gene family / copy / profile matrix** — a *family* is a set of homologous gene *copies*; the
  *profile matrix* is families × extant species, holding copy numbers (or presence/absence).
- **Reconciliation** — the mapping of a gene tree into the species tree, labelling each gene-tree
  node with the duplication, transfer, loss or speciation event that produced it.
- **Event log** — the per-family genealogy of **O** (origination), **D** (duplication),
  **T** (transfer), **L** (loss) and **S** (speciation) events; every reconstructed gene tree is
  post-processing over it.
- **Genome representation** — how much structure a genome carries: *unordered* (a set of families),
  *ordered* (genes on a chromosome), or *nucleotide* (a real sequence). See
  [Genome evolution](#genome-evolution).
- **Coupling / driver / target / edge** — a coevolution coupling is a directed *edge*
  `driver → target`; the *driver*'s state sets the *target*'s rates.
- **Overlay vs output** — a coupling whose arrow does not point into the species tree is an
  *overlay* on a tree you supply; one that points *into* the tree makes the tree an *output*, grown
  jointly. See [Coevolution](#coevolution-the-framework).
