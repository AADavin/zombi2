# A tour of ZOMBI2

ZOMBI2 simulates evolution at four levels, and lets you either run them one after another or grow
them together. This chapter introduces the four levels, the three ways they can relate, and the
single shape every rate takes. It is the vocabulary the rest of the book uses.

## The four levels of ZOMBI2

Evolution leaves its trace at several levels at once. A lineage diversifies into a species tree;
along that tree, genomes gain and lose genes; inside every gene, sequences accumulate substitutions;
and the organisms carry traits that drift and adapt. ZOMBI2 simulates four:

- **Species** — the dated tree of lineages: a strictly bifurcating rooted tree.
- **Genomes** — the genes that exist in each lineage of the tree.
- **Sequences** — the nucleotides or amino acids inside each gene.
- **Traits** — phenotypes evolving along a tree: a body size, a habitat, the presence of a structure.

The first three form a chain of ancestry: a genome lives on the species tree, and a sequence lives
inside a gene. Traits are different. A trait can ride any tree, so it branches to the side rather
than sitting in the chain. You always choose which levels to run, and you need not run them all. A
study of gene content may never descend to sequences; a study of diversification may never leave the
species level.

![The four levels of ZOMBI2. Species, genomes and sequences form a chain of ancestry: a genome lives on the species tree, a sequence inside a gene. Traits branch to the side, because a trait can ride any tree.](figures/fig-2-1-four-levels.svg)

## How the levels relate: independent, conditioned, joint

Everything evolves on the tree, so every level is already conditioned on it. The real question is how
two levels relate to *each other*, and there are exactly three answers. Take a trait and a genome,
both evolving on the same species tree.

- They may ignore each other. Gene loss runs at the same rate everywhere, whatever the trait is
  doing. The two are **independent**.
- One may read the other. Aquatic lineages lose their olfactory genes faster, so the loss rate
  depends on a habitat trait. The genome is **conditioned** on the trait.
- Or neither can go first, because each depends on the other as it unfolds. The two are **joint**.

![Three ways two levels relate. The dashed arrows are the tree, always present; the only thing that changes is what happens between traits and genomes. Independent adds nothing, conditioned adds one arrow, joint grows the two together.](figures/composition.svg)

These are not degrees of one thing; they are three different runs. The notation makes it exact.
Writing `P(Genomes | Species)` for "genomes simulated on a species tree", the three cases are:

- **Independent:** `P(Traits | Species) · P(Genomes | Species)`
- **Conditioned:** `P(Traits | Species) · P(Genomes | Species, Traits)`
- **Joint:** `P(Traits, Genomes | Species)`

The factorisation is not decoration. *Every factor you can write on its own is a run you can do on
its own.* Independent levels are two separate runs, in any order. A conditioned level is still two
runs, but ordered: simulate the driver first, write it to a file, and hand that file to the second
run, exactly as you already hand a species tree to a genome run. A joint pair does not factorise, so
it cannot be split. It is a single run that produces both levels at once.

### When the tree itself is grown

So far the species tree has stayed fixed, sitting behind the bar as something you supply. But a
coupling can point *back into the tree*. If a trait changes the speciation rate, then
faster-speciating lineages leave more descendants, and the tree's own shape comes to depend on the
trait as it grows. The tree stops being an input and becomes an output, grown together with the
trait. In the notation it crosses the bar, into the joint term:

$$P(\text{Species}, \text{Traits}).$$

![Three runs of rising entanglement. A plain pipeline is all tree and no couplings (left). A trait that drives speciation is joint and grows the tree (middle). A trait that also conditions gene loss is joint with the tree and conditioned on by the genome (right).](figures/composition2.svg)

These same three ideas cover the whole book. A classic run, a species tree then genomes then
sequences, is all substrate and no couplings: three independent runs in order. A trait driving
speciation is joint, and the tree is grown. A trait that drives speciation *and* conditions gene loss
mixes both: one joint run, then a conditioned one that reads its output. A single study can combine
all three, because you classify each coupling on its own, not the run as a whole.

Underneath every case, a coupling is one simple thing: a parameter that reads its value from another
level, instead of being a number you type. Part III works through the couplings ZOMBI2 offers, in
that order of difficulty: **Conditioning** for the one-way case, **Joint models** for the inseparable
one, and **Null models** for telling a real coupling from the tree's own noise.

## How rates work

Every event in ZOMBI2, a lineage speciating, a gene duplicating, a site mutating, happens at a
**rate**, and every rate has the same three parts:

$$\text{rate} \;=\; \underbrace{\text{base}}_{\text{how fast}} \;\times\; \underbrace{\text{count}}_{\text{how many, per what}} \;\times\; \underbrace{\text{modifiers}}_{\text{context}}.$$

The **base** is the speed of a single event, in units of inverse time. The **count** is how many
independent opportunities the event has right now, and answering *per what?* is the crux. The
**modifiers** are context multipliers, dimensionless, that let one branch or one family run faster
than another.

**Per what?** Speciation is counted per lineage: every branch alive right now is one opportunity for
the tree to split. Gene gain and loss can be counted per gene copy, so a large family loses genes
faster, or per lineage, which is size-independent. Substitutions are counted per site. This choice,
not the base rate, decides whether a quantity grows exponentially or linearly: a count that tracks
the growing quantity compounds, and one that does not, cannot.

The count is set by a single knob, `per=`, the same at the species and genome levels:

```python
import zombi2 as z

# one opportunity per lineage (the default): births compound, so growth is exponential
per_lineage = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), age=8.0,
                                      direction="forward", seed=3)

# one shared opportunity for the whole tree: the total rate is constant, so growth is linear
shared = z.simulate_species_tree(z.BirthDeath(1.0, 0.2, per="shared"), age=8.0,
                                 direction="forward", seed=3)
```

With the same base rate and seed, the per-lineage tree reaches thousands of tips while the shared one
reaches a few dozen: the same speed, a different number of opportunities.

**Modifiers** change how fast a given event fires relative to the base, never how many opportunities
it has. Some gene families turn over faster than others; some branches simply run hot. The relaxed
molecular clock of the sequences level is one of these: a per-branch modifier on the substitution
rate. That is the one place the word *clock* belongs.

So every rate in the book reads the same way. Ask the two questions, *per what* and *how fast*:

| Level | Counted per… | Speed set by |
|---|---|---|
| **Species** | lineage | the diversification process |
| **Genomes** | copy or lineage | the gain and loss rates |
| **Sequences** | site | the substitution rate, times a clock |

## The ZOMBI2 vocabulary

A few terms recur throughout the book, collected here for reference:

- **Level** — one of the four: Species, Genomes, Sequences, Traits.
- **Rate** — how often an event fires: scope(base) × modifiers, as above.
- **Complete vs reconstructed tree** — the *complete* tree keeps every lineage, including the ones
  that went extinct; the *reconstructed* tree is pruned to the sampled survivors.
- **Extant / extinct / unsampled** — a tip that reached the present and was sampled is *extant*; one
  that died before the present is *extinct*; one alive but not sampled is *unsampled*.
