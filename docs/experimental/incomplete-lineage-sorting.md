# Incomplete lineage sorting

ZOMBI2's core gene-family engine forces every gene lineage to coalesce **exactly at the speciation it
passes through**: at a species split the parent copy is cloned into both daughters, so the
reconstructed gene tree always matches the species tree. That is the *no-ILS* limit. Real gene trees
disagree with the species tree, and one reason is **incomplete lineage sorting** (ILS): looking
backward in time, two lineages can fail to coalesce within the branch that separates their species and
instead sort stochastically in a deeper, ancestral population. This experimental feature simulates
gene trees under the **multispecies coalescent**, where that is exactly what happens.

!!! warning "Experimental"
    Everything here lives in `zombi2.experimental`: usable, but not yet validated for publication, and
    the API may change. It is **pure `numpy`** — no optional dependencies to install. On the command
    line it is reached only through the `zombi2 experimental` group.

## The model: the multispecies coalescent

Read time backward. Inside a branch, the `k` gene lineages present coalesce as a **Kingman
coalescent**: each pair merges at rate `1 / N` per unit branch time, so the waiting time to the next
of the `C(k, 2)` possible coalescences is exponential with total rate `k(k−1) / (2N)`. Any lineages
that have *not* coalesced when the branch reaches its ancestral end **escape** into the parent branch —
where they meet the lineages escaping the *sister* branch and can coalesce with those instead. When
they do, the gene tree disagrees with the species tree. That is ILS.

The single knob is the population size `N` (`population_size`), in the **tree's own time units**. The
amount of ILS is governed by `branch_length / N`:

- **large `N`, or short branches** → coalescence is slow → lineages routinely escape → much discordance;
- **small `N`, or long branches** → coalescence is forced at the nodes → the no-ILS limit (gene tree = species tree).

For the classic rooted triple `((A,B),C)` with internal branch length `T`, the gene tree matches the
species tree with probability

```
P(concordant) = 1 − (2/3)·e^(−T/N)
```

(Hudson 1983; Nei 1987), and the two discordant resolutions `((A,C),B)` and `((B,C),A)` are **equally
likely**. As `T/N → 0` all three topologies approach `1/3` — the deep-ILS limit, and for four or more
taxa the entrance to the *anomaly zone*, where the most probable gene tree need not be the species tree.

Because the gene tree is a genealogy laid out on the **species tree's own time axis**, a coalescence
that predates the species root simply has `time < 0` — the signature of deep coalescence.

## Command line

```bash
# 1000 gene trees under the multispecies coalescent, one allele per species
zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 --seed 1 -o out/
```

| flag | meaning |
|------|---------|
| `-t/--tree` | the species-tree Newick — the coalescent *container* |
| `-N/--population-size` | population size in the tree's time units (larger → more ILS) |
| `-n/--replicates` | number of independent gene trees to draw |
| `-k/--samples` | gene copies sampled per species tip (default 1 = single-copy orthologs) |

It writes `out/gene_trees.nwk` (one gene tree per line), a copy of the species tree, and a run log.
With one copy per species it also reports the fraction of gene trees whose topology matches the
species tree — a direct read-out of how much ILS the chosen `N` produced.

## From the API

```python
from zombi2.experimental import MultispeciesCoalescent
from zombi2.tree import read_newick

species = read_newick(open("species_tree.nwk").read())
msc = MultispeciesCoalescent(population_size=0.5)
gene_trees = msc.sample_gene_trees(species, 1000)     # list[Tree]
```

## DTL + ILS: the coalescent within the locus tree

Run the same coalescent inside a **locus tree** — a gene family's duplication / transfer / loss history
from a `zombi2 genomes` run — and you get gene trees under *DTL + ILS*: the two things that make a gene
tree disagree with the species tree, together. This is the three-level model of SimPhy and DLCoal
(Rasmussen & Kellis 2012): species tree → locus tree (DTL) → gene tree (coalescent).

The routing has one governing idea. A locus is **founded by a single copy** at three kinds of event: the
family **origination**, a **duplication**'s new copy, and a **transferred** (or gene-converted) copy. Its
sampled alleles must therefore coalesce down to that one founder *by the event time* — a coalescent
**conditioned to reach a single lineage within the branch** (the *bounded coalescent*). A **speciation is
not a founding event** (the lineages existed in the ancestral population), so it allows deep coalescence.
That single distinction keeps the model causal: an allele may predate the *speciations* it passes through
(ILS), but never the birth of its own locus.

At the no-ILS limit (`N → 0`) the coalescent is forced back to the event times and reproduces ZOMBI2's own
deterministic gene-tree reconstruction exactly, across duplications, transfers and losses — the model's
tightest self-check.

```bash
# a gene-family run, written with the event trace...
zombi2 genomes -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.1 --write trace -o run/

# ...then a coalescent gene tree per family under DTL + ILS
zombi2 experimental ils -t species_tree.nwk --events-trace run/Events_trace.tsv -N 0.5 -o out/
```

It writes `out/gene_trees/<family>.nwk` (one gene tree per surviving family; `--replicates R` draws `R`
independent gene trees per family). From Python:

```python
from zombi2.experimental import MultispeciesCoalescent
from zombi2.genomes.simulation import simulate_genomes

genomes = simulate_genomes(species_tree, duplication=0.1, transfer=0.05, loss=0.1, output="genomes")
gene_trees = MultispeciesCoalescent(population_size=0.5).sample_family_gene_trees(genomes)  # {family: Tree}
```

!!! note "One coalescent, two containers"
    `sample_gene_tree` and `sample_family_gene_trees` share one engine: the censored coalescent inside a
    *container*. Plain ILS uses the **species tree**; DTL + ILS uses each family's **locus tree**. Only the
    container and the per-event routing (the founder bottleneck above) differ.

!!! note "What this does *not* yet model"
    ILS here sorts *alleles*; gene **copy number** is still assumed to fix instantaneously, so each
    duplication or loss stays pinned to its node in the locus tree. Letting copy number itself sort across
    speciations — a duplicate segregating through a split, i.e. hemiplasy of gene content — is a further
    (and largely open) extension.
