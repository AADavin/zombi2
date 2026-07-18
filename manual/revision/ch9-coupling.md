# Coupling levels

The four levels of the book have so far run past each other. A genome evolves along a species tree, a sequence down a gene tree, a trait along whichever tree you hand it, but nothing at one level has ever read anything at another: aquatic and terrestrial lineages lost genes at the same rate, and a body size drifted the same whether its owner was speciating fast or slow. This chapter is where the levels start to listen. It is the capstone of the book, and it rests on a single idea, already stated in Chapter 2 and made literal here:

> A **coupling** is a parameter that reads its value from another level, instead of being a number you type.

There is exactly one mechanism, a modifier:

```python
loss = 0.25 * mod.Driven(source, mapping)
```

`mod.Driven` reads the driver's value on each lineage and multiplies the base rate by the factor the `mapping` assigns to it. That is the whole of Part III. Everything else in this chapter is about *where the driver comes from*.

*[Draft — the `mod.Driven` modifier and the `joint` command are the design target of `docs/design/coupling-api.md`; they are not built yet. Today the same couplings ship under a single `coevolve` command with a `--couple driver:target` grammar, and as a class zoo in `zombi2.coevolve` (`BiSSE`, `MuSSE`, `QuaSSE`, `HiSSE`, `simulate_trait_conditioned_genomes`, and the rest). The chapter documents the target; the divergences are noted as they arise.]*

## One question splits the whole chapter: can the driver be grown first?

Chapter 2 gave three ways two levels can relate. **Independent** is the case with no coupling at all: the two levels never read each other, and they are two ordinary runs in any order. That leaves the two cases where a coupling exists, **conditioned** and **joint**, and it is tempting to file them as two different topics. They are not. They are the *same* mechanism, `mod.Driven`, separated by one question:

> **Can the driver be grown first, on its own, and handed over?**

- **Yes.** Grow the driver, write it to a file, and pass the file to the level it drives. The `source` is a **filename**. This is **conditioned**: two commands, in order.
- **No** — the driver is entangled with the very thing it drives, each depending on the other as it unfolds, so neither can finish before the other starts. The `source` is a **live level name**, and both grow together in one call. This is **joint**: one command.

Put the two side by side and the seam is a single argument. A habitat trait that drives gene loss can be grown first, because the habitat does not care how many genes a lineage has:

```python
# conditioned: the driver is a file, grown first and handed over
habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.1, seed=1)
habitat.write("habitat.tsv")
genomes.simulate_unordered(tree,
    loss = 0.25 * mod.Driven("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0}), seed=2)
```

A trait that drives *speciation* cannot, because the tree it would be written onto is the very tree its driving is still shaping. So the trait is named live and grown with the tree:

```python
# joint: the driver is a live level, grown alongside what it drives
joint.simulate(
    birth = 1.0 * mod.Driven("trait", {"small": 1.0, "large": 2.0}),   # trait drives speciation
    death = 0.2,
    trait = traits.discrete(states=["small", "large"], switch=0.1),     # grown WITH the tree
    n_tips = 100, seed = 1)
```

Same modifier, same `mapping`. The only difference is that `source` is `"habitat.tsv"` in the first and `"trait"` in the second: a file you already have versus a level you are still growing. This is why conditioning and joint models are one chapter and not two.

The rule underneath is the factorisation rule of Chapter 2, read one more time. *Every factor you can write on its own is a run you can do on its own.* When the driver's factor separates out, `P(Traits | Species) · P(Genomes | Species, Traits)`, you can run it on its own, which is conditioning. When it does not separate, `P(Species, Traits)`, you cannot, which is joint. "Can the driver be grown first?" is just that algebra asked in words.

### Why the framing is not "does it change the tree"

There is a tempting shortcut: joint models are the ones that grow the tree. For version 1 it is even true, because the only live-driver models ZOMBI2 ships are the ones where a level feeds back into the species tree (a trait or gene content driving speciation), so the tree does become an output every time. But that coincidence is a trap. The moment a **tree-fixed** joint model lands, a trait and gene content driving each other on a tree that stays put (SPEC §3), "joint means the tree changed" breaks, while "the driver could not be grown first" still holds without a word of reframing. So take *can the driver be grown first?* as the spine, and treat "v1's joint models all happen to grow the tree" as a passing scope note, not the definition.

## The mapping is the response

`mod.Driven` takes two things: a `source` and a `mapping`. The `source` we have just split into file versus live level. The `mapping` is the other half, and it answers a separate question: once you know the driver's value on a lineage, what factor does the rate get multiplied by? It comes in four shapes, the same four the coupling grammar has always used:

- **Table** — a discrete driver becomes a dict, one factor per state: `{"aquatic": 3.0, "terrestrial": 1.0}`.
- **Curve** — a continuous driver becomes a function: `lambda x: exp(0.5 * x)`, a rate that rises smoothly with body size.
- **Scalar** — a single multiplier, when the driver is already binary and there is nothing to look up.
- **Jump** — the response fires *at an event* rather than continuously: a burst of change at each speciation, not a steady pull along the branch.

The `source` says where the driver lives; the `mapping` says how its value is read. Every coupling in the rest of the chapter is one choice of each.

## Conditioned couplings: grow the driver, write a file, hand it over

A conditioned coupling is the easy half, because you already know how to do both steps. You have grown a driver level in its own chapter, and you have handed a file from one run to the next every time you passed a species tree to a genome run. Conditioning is nothing more than passing a *second* file, a column of driver values, into a rate.

Which pairs allow it is fixed by the geometry of the four levels (SPEC §3). A pair can be conditioned only when neither level lives on the other, that is, when they sit on separate branches of the layout. That is a short list, and **every entry involves Traits**, because traits are the only level off the main chain:

- **A trait drives gene gain or loss** (Traits → Genomes). The habitat example above is the canonical one: cave lineages lose their eye genes faster than surface lineages.

  ```python
  genomes.simulate_unordered(tree,
      loss = 0.25 * mod.Driven("habitat.tsv", {"cave": 4.0, "surface": 1.0}),
      duplication=0.2, transfer=0.1, seed=2)
  ```

- **Gene content drives a trait's optimum** (Genomes → Traits). Grow the genomes first, write a per-lineage gene count, and let it set the value an Ornstein–Uhlenbeck trait is pulled toward, so gene-rich lineages evolve toward a larger body size.

  ```python
  traits.simulate_continuous(tree, start=0.0, rate=1.0, pull=0.5,
      reverts_to = mod.Driven("gene_counts.tsv", lambda n: 0.1 * n), seed=2)
  ```

  *[Draft — this is the one conditioned coupling whose target is a **value** (the OU optimum), not a rate. `mod.Driven` as designed "multiplies the base rate"; driving an optimum needs it to set-or-scale a value instead. The spelling here is provisional and the question is flagged under *Still to design*.]*

- **A trait drives selection or clock speed on sequences** (Traits → Sequences). The trait sets either the substitution rate (a faster clock on some lineages) or the model's `dN/dS`, so a habitat shift relaxes purifying selection on a gene.

  ```python
  sequences.simulate_sequences(gene_trees, model=gy94(omega=0.2),
      substitution = 1.0 * mod.Driven("habitat.tsv", {"fast": 3.0, "slow": 1.0}),
      length=300, seed=2)
  ```

  A sequence rides a *gene* tree, not the species tree the trait was grown on, so each gene-tree branch reads the trait value of the species branch it sits inside. That mapping through the reconciliation is handled for you, but it is worth knowing it is there.

In every case the coupling **folds into the target level's own command**. There is no separate coupling step and no coupling object to build; you grow the driver, then call the ordinary genome, trait, or sequence run with one rate that happens to be `Driven` instead of a bare number. On the command line the same shape appears as a flag on the target subcommand, `--loss-driven-by habitat.tsv`, rather than a command of its own.

Two directions are conditionings you might expect and cannot have. A genome cannot be conditioned on the species tree, and a sequence cannot be conditioned on its genome, because in each case one level *lives on* the other rather than beside it, and a level always reads the tree it lives on already, for free. Those are not couplings at all. And the reverse of the last case, a **sequence** driving a trait or a gene, is real but deferred: sequences are target-only in version 1, because driving *out* of a sequence needs the substitution step broken up mid-branch (SPEC §10). The hole is genuine; it is not pretended shut.

## Joint models: grow the two together

When the driver cannot be grown first, there is nothing to write to a file, because the file would have to be written onto a tree that does not exist yet. The two levels are grown in one Gillespie that fires both kinds of event against each other: a speciation event uses the current trait value to set its rate, and a trait-change event evolves the trait on the tree as it has grown so far. Out comes both levels at once. This is the `joint` command, and it is the only place in the book that produces two levels from one call.

Version 1 ships two joint pairs, and both are cases where a level reaches back into the species tree, so the tree itself is grown as an output:

- **A trait drives speciation** — `P(Species, Traits)`. The worked example above is the whole of it: a body-size state that makes large lineages speciate twice as fast. Letting `death` be `Driven` too makes extinction state-dependent as well, which is the model the literature calls BiSSE.

- **Gene content drives speciation** — `P(Species, Genomes)`. The presence of a key gene, a toxin, a transporter, lifts a lineage's speciation rate; the genome and the tree grow together. The genome enters as a **process spec** rather than a finished run, the same way the trait did:

  ```python
  joint.simulate(
      birth = 1.0 * mod.Driven("genome", {"toxin_present": 1.8, "toxin_absent": 1.0}),
      death = 0.2,
      genome = genomes.unordered(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5),
      n_tips = 100, seed = 1)
  ```

The `mapping` shapes come back here with names the literature will recognise. A **Table** on a discrete trait is BiSSE (two states) or MuSSE (several); a **Curve** on a continuous trait is QuaSSE; a **Jump** response is a change that fires *at* each split rather than along the branch, which is how a trait both drives speciation and shifts at cladogenesis (ClaSSE), or how gene content drives speciation with a burst of gene gain and loss at every node. Same four shapes as the conditioned half; only the `source` moved from a file to a live level.

*[Scope note — the two pairs above both grow the tree. The third joint case in SPEC §3, a trait and gene content driving *each other* on a tree that stays fixed, is genuinely joint (neither can be grown first) but does not touch the tree, and it is deferred to `experimental`. It is the case that proves "joint" cannot mean "changes the tree".]*

### Not everything that looks like a connection is one

One distinction keeps this chapter from swallowing material that belongs elsewhere. A trait that jumps at a speciation event, or a genome that changes only at splits, looks like a coupling to the species level, but it is not one. It is the level reading the tree it *already lives on*, which every level does for free. A cladogenetic trait shift and a punctuational burst of gene change are options of the trait's and the genome's own models, and they stay in Chapters 8 and 6 respectively. A coupling, in this chapter's sense, always reads a *different* level. When in doubt, ask whether the rate reads a value that some *other* level produced; if it only reads the tree, it is not here.

## The literature → command bridge

The state-dependent and key-innovation models arrive under a wall of acronyms, and a reader who wants "a BiSSE model" or "QuaSSE" should be able to find the door. The names live here, in one table, and organise nothing else in the chapter.

| Literature name | What it does | ZOMBI2 |
|---|---|---|
| BiSSE | a binary trait drives speciation (and extinction) | `joint`: `birth = … * mod.Driven("trait", {…})`, a Table |
| MuSSE | a multi-state trait drives speciation | `joint`, a Table with more states |
| QuaSSE | a continuous trait drives speciation | `joint`, a Curve `lambda x: …` |
| HiSSE | a hidden state drives speciation | `joint`, a trait carrying hidden rate classes |
| ClaSSE | a trait drives speciation *and* shifts at each split | `joint` + a Jump response |
| Key innovation | the presence of a gene drives speciation | `joint`: `birth = … * mod.Driven("genome", {…})` |
| Trait–gene feedback / co-diversification | a trait and gene content drive each other, tree fixed | `joint`, tree-fixed — deferred to `experimental` |

As with every bridge table in the book, the acronyms sit in this one place so they can be searched for, and never in a section heading. The class names (`BiSSE`, `QuaSSE`, and the rest) survive in `zombi2.coevolve` for exactly the same reason: they are the field's search terms.

## Outputs

A conditioned run writes what any ordinary level run writes, its genome, trait, or sequence output, plus the **driver file** alongside it, so the pairing that produced the pattern is kept on disk. A joint run writes **both** levels: the grown species tree together with the trait history, or together with the genomes, each in the same format it would have had from its own command. Because a joint run grows the tree, the tree it writes is a *complete* tree in the sense of Chapter 4, with the extinct lineages that shaped the trait or gene distribution still in place. The full list of files lives in Appendix B.

## Nulls: a recipe, not a feature

The chapter closes on the question every coupling eventually raises: given a tree where the trait and the gene loss really do line up, how do you know the coupling is real and not an accident of the tree? This matters enough to end on, and the honest answer is short.

A tree manufactures associations on its own. Two things that both evolve down the same branches will look correlated at the tips whether or not either drives the other, because they share ancestry and they share the tree's timing. So a pattern is not evidence of a coupling until you know what *no* coupling produces **on that same tree**. The baseline has to be simulated, on the identical tree, from a model where the coupling is switched off, and then compared against what you saw.

ZOMBI2 gives you that baseline with the primitives already in the book. There is no null function and no `--null` flag, because a null is just a run with the coupling changed, and each of the three standard nulls is a one-line edit to the rate:

```python
loss = 0.25 * mod.Driven("habitat.tsv", {"cave": 4.0, "surface": 1.0})   # the coupling under test

loss = 0.25                                    # independent null — drop the coupling entirely
loss = 0.25 * mod.ByBranch(spread=0.5)         # CID null — rate varies across the tree, but NOT by the trait
loss = 0.25 * mod.Driven(shuffle("habitat.tsv"), {"cave": 4.0, "surface": 1.0})   # shuffle null — break the pairing
```

Each says something different. The **independent** null asks whether *any* rate variation is needed at all. The **CID** null is the sharper one: it lets the loss rate vary across the tree just as much as the real coupling would, but tied to nothing, so it separates "the trait matters" from "some branches simply run hot". It is literally the across-branch clock modifier of Chapter 7, `mod.ByBranch`, reused here as a null. The **shuffle** null keeps the trait's own distribution intact and only permutes which tip carries which value, so it destroys the pairing while preserving everything else. Wrap any of them in `for seed in range(100)` and you have a distribution to compare against.

What ZOMBI2 does **not** do is finish the test. Choosing the association measure and computing a p-value from that distribution is inference, not simulation, and it is the user's job. The simulator generates the baseline and stops there. The only utility that might earn a place in `tools` is a tiny `shuffle()` helper for permuting a tip-value file; even that is a convenience, not a null subsystem.

## Still to design

Part III is the newest of the design targets, and several of its corners are agreed in shape but not in spelling. The chapter should not pretend otherwise.

- **The name `Driven`.** `Driven` versus `DependsOn` is unsettled, as is the exact convention for naming a live driver inside `joint` (`"trait"` and `"genome"` as bare strings are placeholders for whatever referencing scheme the level names settle into).
- **Driving a value, not a rate.** The gene-content-drives-optimum coupling drives an OU optimum, which is a value, while `mod.Driven` is designed to multiply a *rate*. Whether `Driven` learns to set-or-scale a value, or that case gets its own spelling, is open.
- **The process-spec versus runner split.** Inside `joint`, a level enters as a *description* (`traits.discrete(...)`, `genomes.unordered(...)`) rather than a run (`traits.simulate_discrete(...)`, `genomes.simulate_unordered(...)`). This spec/run pair is the intended shape but is not yet built or named on either level.
- **State-dependent extinction.** Driving `death` as cleanly as `birth`, so that a trait raises the extinction rate as well as the speciation rate, needs its surface pinned down; the `death = … * mod.Driven(...)` form above is the intent, not a settled API.
