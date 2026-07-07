# Ordered chromosomes

By default a genome in ZOMBI2 is **order-free** — a multiset of gene families with copy numbers, all
you need for phylogenetic profiles (see [gene families](gene-family.md)). The **ordered** level adds
one more thing genes carry: **position**. Genes sit on an ordered, circular chromosome, each with a
strand orientation, and the chromosome evolves not just by gaining and losing genes but by *shuffling*
them — inversions and transpositions rearrange contiguous segments so that gene **order** itself
carries phylogenetic signal. Reach for it when synteny, operons, or gene-order evolution matter.
Selecting the level is `--genome-model ordered`.

| Model | Substrate | Rearrangements | Reach for it when |
| --- | --- | --- | --- |
| **Ordered** | circular chromosome of oriented genes (no intergenes) | inversion, transposition on gene segments (distance in genes) | synteny, operons, gene-order evolution matter |

## The models

### Ordered genome

A circular chromosome of genes, each carrying a strand orientation, with **no intergenic regions** —
the basic ZOMBI1 representation. On top of the usual duplication / transfer / loss / origination, two
rearrangements act on a **contiguous segment** whose length is drawn from an `extension` continuation
probability (`extension=None` → single genes; higher → longer segments): **inversion** reverses the
segment and flips every strand, and **transposition** cuts and pastes it elsewhere. Both change gene
order/orientation but **not** gene content, so they leave the profile matrix and the gene trees
unchanged and show up only in the event log and the final chromosome order. Rearrangement rates are
per gene copy; segment length is set by `--mean-length` (in genes). Rearrangements require the
`shared` rate model.

## Command line

`--genome-model ordered` selects the level; `--inversion`/`--transposition` are the rearrangement
rates (per gene copy), and `--mean-length` sets the segment length (in genes).

```bash
# ordered: gene-order rearrangements on a circular chromosome
zombi2 genomes -t species_tree.nwk --genome-model ordered \
    --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.4 \
    --inversion 0.3 --transposition 0.3 --mean-length 2 \
    --initial-families 30 --write profiles trees --seed 1 -o out/
```

## Python

Models live in `zombi2.genomes` (and re-export at the top level, so `zombi2.OrderedGenome` also
works). The ordered level runs through `simulate_genomes` with an `OrderedGenome` factory:

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import simulate_genomes, SharedRates, OrderedGenome

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

# rearrangements need OrderedGenome, which takes the extension knob
rates = SharedRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.4,
                    inversion=0.3, transposition=0.3)
genomes = simulate_genomes(tree, rates, initial_families=30, seed=1,
                           genome_factory=lambda ids: OrderedGenome(ids, extension=0.5))
leaf = next(iter(genomes.leaf_genomes.values()))
leaf.chromosome     # ordered list of OrderedGene(gid, family, orientation=±1)
```

## Output

The ordered level writes the usual gene-family output — `Profiles.tsv` / `Presence.tsv` (copy-number
and presence matrices over extant leaves), `species_nodes.tsv`, and per-family reconstructed gene
trees under `gene_trees/` when `trees` is requested; inversions and transpositions appear in the event
log and the final chromosome order, not in the profiles. `species_tree.nwk` is always written and
`genomes.log` is the run manifest.

## Validation

- **Ordered.** The mean inversion-event count equals `inversion_rate × initial_families ×
  total_branch_length` — a Poisson oracle: because rearrangements conserve genome content the size
  stays constant, so the integrated per-branch hazard is exact and the observed mean over many
  fixed-tree replicates lands within Monte-Carlo error of it
  (`test_ordered_genome.py::test_inversion_count_matches_poisson_mean`).

## References

- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a
  phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.
  *Bioinformatics* 36(4): 1286–1288.
