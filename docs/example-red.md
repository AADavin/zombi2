# An example: can you trust RED?

This page is here to show **the kind of question ZOMBI2 exists to answer**.

Phylogenetics is full of methods that cannot be checked on the data they are meant for, because that
data hides the answer. You can only grade a method where the truth is known — and the truth is known
in a simulation, by construction. So the recipe is always the same three moves: measure one honest
number on real data, reproduce that number in a simulation where the answer *is* known, and grade the
method there.

What follows is one worked instance of that recipe, start to finish. The full study, with its code
and data, lives in [`analyses/red/`](https://github.com/AADavin/zombi2/tree/main/analyses/red).

## The question

**Relative Evolutionary Divergence** (RED; Parks et al. 2018) turns a phylogram into a relative
divergence scale. Walking root to tip, it places every node along its path in proportion to
accumulated branch length, giving a number between 0 (root) and 1 (tips). GTDB uses it to normalise
taxonomic ranks across the tree of life, so that a phylum sits at a comparable RED whatever lineage
it belongs to.

That works only if RED is a faithful stand-in for relative divergence *time*. Under a strict clock it
is exact — branch length is proportional to time. Under rate variation it is not, and real lineages
do not share a clock. **Does RED still recover the ages?**

## Why real data cannot answer it

A phylogram measures substitutions, which are rate × time. From substitutions alone, rate and time
are jointly unidentifiable — so the true node ages you would grade RED against cannot be read off the
tree. Dating the tree first would assume a rate model, which is the thing in question, and the test
becomes circular.

## The one number we borrow

Every genome in the GTDB archaeal tree is extant, so every tip is the same amount of *time* from the
root. Any spread in root-to-tip *substitutions* can therefore only come from rate variation. The
coefficient of variation of those distances is a model-free summary of how ragged real trees are:

![The GTDB archaeal root-to-tip distribution](assets/red/observable.png)

**CV = 0.232** across 10,122 genomes — the fastest lineage has accumulated roughly four times the
substitutions of the slowest. That single number is all we take from the real world. No real branch
length enters the simulation.

## Grading the method where truth is known

In ZOMBI2 the dated tree *is* the truth, and the phylogram is what a real study would have seen:

```python
import numpy as np
from zombi2 import genomes, sequences, species
from zombi2.rates import modifiers as mod
from zombi2.sequences import substitution_models as sm
from zombi2.tree import read_newick
from zombi2.tree import relative_evolutionary_divergence as red_of

# a dated tree — we know its node ages, because we simulated them
sp = species.simulate_species_tree(birth=1.0, n_extant=200, seed=100)
truth = red_of(sp.extant_tree)                 # RED is exact on a dated tree

# evolve it under a relaxed clock and read the phylogram that leaves behind
g = genomes.simulate_genomes_family(sp.complete_tree, initial_families=5,
                                    duplication=0.01, loss=0.01, seed=100)
seq = sequences.simulate_sequences(
    g, model=sm.jc69(), length=1,
    substitution=1.0 * mod.ByLineage(spread=0.544, dist="lognormal"), seed=7)

phylogram, _ = read_newick(seq.species_phylogram["extant"])
estimate = red_of(phylogram)                   # RED, now on substitutions instead of time

nodes = [i for i in truth if i in estimate and sp.extant_tree.nodes[i].children is not None]
r = np.corrcoef([truth[i] for i in nodes], [estimate[i] for i in nodes])[0, 1]
print(f"RED recovers relative node age with r = {r:.3f} across {len(nodes)} nodes")
```

That is one tree and one draw of the clock. The study repeats it over eight trees and sweeps the
clock's heterogeneity σ, so the answer comes with an error bar rather than a seed.

## Calibrating, then reading off

Sweeping σ and finding where the simulated raggedness crosses the real one gives the clock that makes
a simulated tree as ragged as real archaea:

![Calibrating the clock](assets/red/clock_recovery.png)

One CV can say *how much* rate variation there is, but not how it is **arranged** — and arrangement
matters for a method that walks root to tip. So all three clocks ZOMBI2 wires at the sequence level
are carried through: **uncorrelated** (`ByLineage`, every lineage drawing independently, lognormal or
gamma tail) and **autocorrelated** (`FromParent`, each lineage inheriting its parent's rate, so
relatives evolve alike). They reach the same raggedness at very different σ, which is why the
comparison is made at matched CV.

![RED accuracy against raggedness](assets/red/red_bridge.png)

## The answer

At the raggedness real archaea show:

| clock | Pearson r | error (% of tree depth) |
|---|---|---|
| uncorrelated, lognormal | 0.953 | 5.8% |
| uncorrelated, gamma | 0.942 | 6.1% |
| autocorrelated | 0.993 | 2.3% |

**RED holds up at real archaeal raggedness**, and it holds up under every arrangement of rate
variation we can put it under — so the conclusion does not rest on guessing which one archaea have.
Uncorrelated rates are the harder case, so 0.94 is the conservative number to quote.

Two things worth being precise about. RED is an **ordinal** proxy: even at its best there is a
few-percent age error, so use it to order divergences and normalise ranks — its designed job — not to
read absolute times off. And it **does** break down, just past where real data sits: by CV ≈ 0.46 the
correlation falls to ≈ 0.79 and the error rises to ≈ 9%. Real archaea are on the safe side of that,
which is a quantitative version of the assumption GTDB relies on.

## The shape of the recipe

Nothing above is specific to RED. The move is general, and it is what a simulator is *for*:

> A method that cannot be graded on the data it is meant for can still be graded — if one honest
> number says how demanding the real case is, that number is reproduced where the answer is known,
> and the method is graded there.

The full write-up, including assumptions, limitations and every number's provenance, is in
[`analyses/red/REPORT.md`](https://github.com/AADavin/zombi2/blob/main/analyses/red/REPORT.md).
It regenerates from fixed seeds in three commands.

## References

- Drummond AJ, Ho SYW, Phillips MJ, Rambaut A (2006). *Relaxed phylogenetics and dating with
  confidence.* PLoS Biology 4:e88.
- Parks DH, Chuvochina M, Waite DW, et al. (2018). *A standardized bacterial taxonomy based on genome
  phylogeny substantially revises the tree of life.* Nature Biotechnology 36:996–1004.
- Rinke C, Chuvochina M, Mussig AJ, et al. (2021). *A standardized archaeal taxonomy for the Genome
  Taxonomy Database.* Nature Microbiology 6:946–959.
- Thorne JL, Kishino H, Painter IS (1998). *Estimating the rate of evolution of the rate of molecular
  evolution.* Molecular Biology and Evolution 15(12):1647–1657.
