# Incomplete lineage sorting

ZOMBI2's gene-family engine forces every gene lineage to coalesce *exactly at the speciation it passes
through*: at a species split the parent copy is cloned into both daughters, so a reconstructed gene tree
always matches the species tree. That is a modelling limit, not a fact of biology. Real gene trees
disagree with the species tree, and one of the fundamental reasons is **incomplete lineage sorting**
(ILS): looking backward in time, two lineages can fail to find a common ancestor within the branch that
separates their species, and instead sort — stochastically — in a deeper ancestral population. This
supplement describes an **experimental** ZOMBI2 feature that simulates gene trees under the
**multispecies coalescent**, where that is exactly what happens.

::: warning
Everything in this supplement lives in `zombi2.experimental`. It is shipped so you can use and iterate on
it, but it has **not** yet cleared ZOMBI2's core bar: the API may change, and the outputs are not yet
validated for publication. It needs **no optional dependencies** — it is
pure `numpy`. This is a separate, self-contained document; it is **not** part of the main manual.
:::

## The multispecies coalescent

Read time backward. Inside a branch, the $k$ gene lineages present coalesce as a **Kingman coalescent**:
each pair merges at rate $1/N$ per unit branch time, so the waiting time to the next of the
$\binom{k}{2}$ possible coalescences is exponential with total rate $k(k-1)/2N$. Any lineages that have
*not* coalesced when the branch reaches its ancestral end **escape** into the parent branch — where they
meet the lineages escaping the *sister* branch and can coalesce with those instead. When they do, the
gene tree disagrees with the species tree.

The single control is the population size $N$ (`population_size`), in the tree's own time units. The
amount of ILS is governed by $\text{branch length}/N$: large $N$ or short branches give slow coalescence
and much discordance; small $N$ or long branches recover the no-ILS limit, in which coalescence is forced
back to the nodes and the gene tree equals the species tree. Two facts pin the model down. For the rooted
triple $((A,B),C)$ with internal branch length $T$, the gene tree matches the species tree with
probability

$$P(\text{concordant}) \;=\; 1 - \tfrac{2}{3}\,e^{-T/N},$$

and the two discordant resolutions are equally likely; as $T/N \to 0$ all three approach $1/3$ — the
deep-ILS limit, and for four or more taxa the entrance to the *anomaly zone*, where the most probable
gene tree need not be the species tree.

Because a gene tree is laid out on the species tree's own time axis, a coalescence that predates the
species root simply sits at negative time — the honest signature of deep coalescence.

::: note
ZOMBI2 never builds gene trees during the forward simulation; it reconstructs them afterwards from the
event log. ILS is therefore a **separate backward pass** over that log — the forward Gillespie engines,
including the fast Rust path, are untouched.
:::

## Plain ILS: the coalescent in the species tree

The first mode runs the coalescent inside a species tree, sampling one allele per species (single-copy
orthologs). It answers a clean question: given this dated species tree and this $N$, how much do gene
trees disagree with it?

```
zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 --seed 1 -o out/
```

writes `out/gene_trees.nwk` (one gene tree per line) and reports the fraction whose topology matches the
species tree — a direct read-out of how much ILS the chosen $N$ produced. From Python:

```python
from zombi2.experimental import MultispeciesCoalescent
from zombi2.tree import read_newick

species = read_newick(open("species_tree.nwk").read())
gene_trees = MultispeciesCoalescent(population_size=0.5).sample_gene_trees(species, 1000)
```

## DTL + ILS: the coalescent in the locus tree

The second mode is the realistic combination: run the same coalescent inside a **locus tree** — a gene
family's duplication, transfer and loss history from a `zombi2 genomes` run. This is the three-level
model of the field (species tree, then locus tree, then gene tree; Rasmussen & Kellis 2012; SimPhy,
Mallo et al. 2016), with ZOMBI2's own event model — transfers included — supplying the locus tree.

The routing rests on a single idea. A locus is **founded by a single copy** at three kinds of event: the
family **origination**, a **duplication**'s new copy, and a **transferred** (or gene-converted) copy. Its
sampled alleles must therefore coalesce down to that one founder *by the event time* — a coalescent
**conditioned to reach a single lineage within the branch**, the *bounded coalescent* (its within-branch
probabilities are the classical coalescent transition probabilities of Tavaré 1984). A **speciation is
not a founding event** — the lineages existed in the ancestral population — so it allows deep
coalescence. That one distinction is what keeps the model causal: an allele may predate the *speciations*
it passes through (ILS), but never the birth of its own locus.

::: tip
The tightest self-check falls straight out of this. At the no-ILS limit ($N \to 0$) every coalescence is
forced back to its event time, and the model reproduces ZOMBI2's *own* deterministic gene-tree
reconstruction — exactly, across duplications, transfers, losses and originations.
:::

The feature is a post-process on a genomes run written with its event trace:

```
zombi2 genomes -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.1 --write trace -o run/
zombi2 experimental ils -t species_tree.nwk --events-trace run/Events_trace.tsv -N 0.5 -o out/
```

The second command writes `out/gene_trees/<family>.nwk`, one coalescent gene tree per surviving family
(or `--replicates R` gene trees per family). From Python the same is one call on a `Genomes` result:

```python
from zombi2.experimental import MultispeciesCoalescent
from zombi2.genomes.simulation import simulate_genomes

genomes = simulate_genomes(species_tree, duplication=0.1, transfer=0.05, loss=0.1, output="genomes")
gene_trees = MultispeciesCoalescent(population_size=0.5).sample_family_gene_trees(genomes)
```

::: note
Both modes share one engine: the censored coalescent inside a *container*. Plain ILS uses the species
tree; DTL + ILS uses each family's locus tree. Only the container and the per-event routing (the founder
bottleneck above) differ.
:::

## Status and limitations

- **Experimental.** In `zombi2.experimental`, reached on the command line only through
  `zombi2 experimental ils`. Pure `numpy`, no optional dependencies. APIs and outputs may change.
- **One global population size.** $N$ is constant across the tree; a per-branch (heterogeneous) $N$ is a
  planned refinement, as is a species origin/stem to bound the deepest coalescences instead of the
  standard unbounded ancestral population above the root.
- **Alleles, not copy number.** ILS here sorts the *sequence* of a copy that exists; gene **copy number**
  is still assumed to fix instantaneously, so each duplication or loss stays pinned to its node in the
  locus tree. Letting copy number itself sort across speciations — a duplicate segregating through a
  split, i.e. **hemiplasy of gene content**, in which one duplication reads as a gain plus a loss — is a
  further, largely open extension.
