# Genome structure

By default a genome in ZOMBI2 is **order-free** — a multiset of gene families with copy numbers,
all you need for phylogenetic profiles. The models here add **structure**: where genes sit, how they
are oriented, and how structural events (inversions, transpositions, indels) rearrange them.
Two levels are available, selected with `--genome-model`. The **ordered** level places genes on a
circular chromosome and evolves gene-*order* by segment rearrangements, counting distance in genes.
The **nucleotide** level works one level down — a genome is a sequence of individual nucleotides and
events act on variable-length segments of them, so paralogy, xenology, gene order and orientation are
all resolved at nucleotide resolution and a gene tree is reconstructed for every stretch of shared
ancestry.

| Model | Substrate | Rearrangements | Reach for it when |
| --- | --- | --- | --- |
| **Ordered** | circular chromosome of oriented genes (no intergenes) | inversion, transposition on gene segments (distance in genes) | synteny, operons, gene-order evolution matter |
| **Nucleotide** | sequence of nucleotides; genes emerge as *blocks* | inversion, transposition, intergenic insertion/deletion (distance in nucleotides) | you need nucleotide-resolution structure, real genomes, or per-block gene trees |

## The models

### Ordered genome

A circular chromosome of genes, each carrying a strand orientation, with **no intergenic regions** —
the basic ZOMBI1 representation. On top of the usual duplication / transfer / loss / origination,
two rearrangements act on a **contiguous segment** whose length is drawn from an `extension`
continuation probability (`extension=None` → single genes; higher → longer segments): **inversion**
reverses the segment and flips every strand, and **transposition** cuts and pastes it elsewhere. Both
change gene order/orientation but **not** gene content, so they leave the profile matrix and the gene
trees unchanged and show up only in the event log and the final chromosome order. Rearrangement rates
are per gene copy; segment length is set by `--mean-length` (in genes). Rearrangements require the
`shared` rate model.

### Nucleotide genome

A genome is a sequence of individual nucleotides, starting from `initial_chromosomes` chromosome(s)
of `root_length` nucleotides at the root. Duplication, transfer, loss, inversion and transposition
act on variable-length segments and are **per-nucleotide** rates (the total genome rate is
`rate × current_length`, so longer genomes evolve faster); origination is **per branch**; segment
length is geometric with mean `1/(1 − extension)` nucleotides (`extension=0.99` → ~100 nt). Two extra
knobs edit only intergene positions: `insertion` lays down a run of novel nucleotides (a fresh block)
and `deletion` removes a run from within a single intergene, each with mean `indel_mean_length`. The
simulator partitions surviving material into **blocks** — maximal segments of unbroken shared
ancestry — and every block gets its own reconstructed gene tree; `profile_matrix`, `leaf_mosaic`
(ordered signed blocks) and `trace_back` (each nucleotide's ancestral origin) read out a leaf.
Declaring `gene_intervals` (or a real `--gff`) switches on *genic mode*, where genes are never split
and pseudogenization / homologous-replacement transfer become available (see
[nucleotide genomes](../guide/nucleotide-genomes.md)).

!!! warning "Keep gain ≤ loss"
    Duplication and additive transfer grow the genome without a cap. Over long ages keep them at or
    below `loss` to avoid runaway growth.

## Command line

`--genome-model ordered|nucleotide` selects the level; `--inversion`/`--transposition` are the
rearrangement rates (per gene copy for ordered, per nucleotide for nucleotide), and `--mean-length`
sets the segment length (genes for ordered, nucleotides for nucleotide).

```bash
# ordered: gene-order rearrangements on a circular chromosome
zombi2 genomes -t species_tree.nwk --genome-model ordered \
    --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.4 \
    --inversion 0.3 --transposition 0.3 --mean-length 2 \
    --initial-families 30 --write profiles trees --seed 1 -o out/

# nucleotide: structural events at nucleotide resolution, blocks + per-block trees
zombi2 genomes -t species_tree.nwk --genome-model nucleotide \
    --inversion 0.001 --transposition 5e-5 --loss 1.5e-4 --dup 1e-4 --trans 5e-5 --orig 0.2 \
    --root-length 1000 --write profiles trees --seed 1 -o out/
```

## Python

Models live in `zombi2.genomes` (and re-export at the top level, so `zombi2.OrderedGenome` also
works). The ordered level runs through `simulate_genomes` with an `OrderedGenome` factory; the
nucleotide level has its own entry point `simulate_nucleotide_genomes`:

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import (
    simulate_genomes, SharedRates, OrderedGenome, simulate_nucleotide_genomes)

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

# ordered: rearrangements need OrderedGenome, which takes the extension knob
rates = SharedRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.4,
                    inversion=0.3, transposition=0.3)
genomes = simulate_genomes(tree, rates, initial_families=30, seed=1,
                           genome_factory=lambda ids: OrderedGenome(ids, extension=0.5))
leaf = next(iter(genomes.leaf_genomes.values()))
leaf.chromosome     # ordered list of OrderedGene(gid, family, orientation=±1)

# nucleotide: variable-length structural events; genes emerge as blocks
result = simulate_nucleotide_genomes(
    tree, root_length=1000,
    duplication=1e-4, transfer=5e-5, loss=1.5e-4,
    inversion=1e-3, transposition=5e-5, origination=0.2,
    insertion=1e-4, deletion=1e-4, seed=1)
block_ids, species, matrix = result.profile_matrix()   # copy number per extant leaf
result.leaf_mosaic(tree.leaves()[0])   # the leaf as ordered, signed blocks
```

## Output

**Ordered** writes the usual gene-family output — `Profiles.tsv` / `Presence.tsv` (copy-number and
presence matrices over extant leaves), `species_nodes.tsv`, and per-family reconstructed gene trees
under `gene_trees/` when `trees` is requested; inversions and transpositions appear in the event log
and the final chromosome order, not in the profiles. **Nucleotide** additionally emits the block-based
architecture: `Profiles.tsv`/`Presence.tsv` are over **blocks**, `blocks.tsv` describes each block
(and its `kind`/`gene_id` in genic mode), `Mosaics.tsv` gives every leaf as an ordered signed block
sequence, `gene_trees/` holds one reconstructed tree per block, and `Reconciled_complete.nwk` /
`Reconciled_extant.nwk` / `Reconciliation_events.tsv` record the block reconciliations. `--write
ancestral` additionally simulates DNA and reconstructs the genome at every node. `species_tree.nwk`
is always written and `genomes.log` is the run manifest.

## Validation

- **Ordered.** The mean inversion-event count equals `inversion_rate × initial_families ×
  total_branch_length` — a Poisson oracle: because rearrangements conserve genome content the size
  stays constant, so the integrated per-branch hazard is exact and the observed mean over many
  fixed-tree replicates lands within Monte-Carlo error of it
  (`test_ordered_genome.py::test_inversion_count_matches_poisson_mean`).
- **Nucleotide — inversion oracle.** Applying random inversions to a nucleotide genome reproduces an
  independent array oracle cell-for-cell, and every inversion preserves length
  (`test_nucleotide_genome.py::test_inversion_matches_oracle_random`).
- **Nucleotide — all events oracle.** A full mixed stream of structural events (duplication,
  transfer, loss, inversion, transposition) matches the array oracle exactly
  (`test_nucleotide_genome.py::test_all_events_match_oracle`).
- **Nucleotide — indel lengths.** Intergenic insertion/deletion run lengths are geometric with mean
  equal to `indel_mean_length` (`test_nucleotide_indels.py::test_draw_indel_length_matches_mean`).

## References

- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a
  phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.
  *Bioinformatics* 36(4): 1286–1288.
