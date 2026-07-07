# Nucleotide genomes

The **nucleotide** level works one level down from [ordered chromosomes](ordered.md): a genome is a
sequence of individual nucleotides and events act on variable-length **segments** of them, so
paralogy, xenology, gene order and orientation are all resolved at nucleotide resolution. Surviving
material is partitioned into **blocks** — maximal stretches of unbroken shared ancestry — and a gene
tree is reconstructed for every block. Reach for it when you need nucleotide-resolution structure,
want to start from a real genome, or need per-block gene trees. Selecting the level is
`--genome-model nucleotide`.

| Model | Substrate | Rearrangements | Reach for it when |
| --- | --- | --- | --- |
| **Nucleotide** | sequence of nucleotides; genes emerge as *blocks* | inversion, transposition, intergenic insertion/deletion (distance in nucleotides) | you need nucleotide-resolution structure, real genomes, or per-block gene trees |

## The models

### Nucleotide genome

A genome is a sequence of individual nucleotides, starting from `initial_chromosomes` chromosome(s) of
`root_length` nucleotides at the root. Duplication, transfer, loss, inversion and transposition act on
variable-length segments and are **per-nucleotide** rates (the total genome rate is
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

`--genome-model nucleotide` selects the level; `--inversion`/`--transposition` are the rearrangement
rates (per nucleotide), and `--mean-length` sets the segment length (in nucleotides). `--root-length`
sets the root chromosome length; `--insertion`/`--deletion` with `--indel-mean-length` edit intergene
positions; and declaring genes with `--genes` or `--gff` switches on genic mode
(`--pseudogenization`, `--replacement`).

```bash
# nucleotide: structural events at nucleotide resolution, blocks + per-block trees
zombi2 genomes -t species_tree.nwk --genome-model nucleotide \
    --inversion 0.001 --transposition 5e-5 --loss 1.5e-4 --dup 1e-4 --trans 5e-5 --orig 0.2 \
    --root-length 1000 --write profiles trees --seed 1 -o out/
```

## Python

Models live in `zombi2.genomes` (and re-export at the top level). The nucleotide level has its own
entry point `simulate_nucleotide_genomes`:

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import simulate_nucleotide_genomes

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

# variable-length structural events; genes emerge as blocks
result = simulate_nucleotide_genomes(
    tree, root_length=1000,
    duplication=1e-4, transfer=5e-5, loss=1.5e-4,
    inversion=1e-3, transposition=5e-5, origination=0.2,
    insertion=1e-4, deletion=1e-4, seed=1)
block_ids, species, matrix = result.profile_matrix()   # copy number per extant leaf
result.leaf_mosaic(tree.leaves()[0])   # the leaf as ordered, signed blocks
```

## Output

The nucleotide level emits the block-based architecture: `Profiles.tsv`/`Presence.tsv` are over
**blocks**, `blocks.tsv` describes each block (and its `kind`/`gene_id` in genic mode), `Mosaics.tsv`
gives every leaf as an ordered signed block sequence, `gene_trees/` holds one reconstructed tree per
block, and `Reconciled_complete.nwk` / `Reconciled_extant.nwk` / `Reconciliation_events.tsv` record
the block reconciliations. `--write ancestral` additionally simulates DNA and reconstructs the genome
at every node. `species_tree.nwk` is always written and `genomes.log` is the run manifest.

## Validation

- **Inversion oracle.** Applying random inversions to a nucleotide genome reproduces an independent
  array oracle cell-for-cell, and every inversion preserves length
  (`test_nucleotide_genome.py::test_inversion_matches_oracle_random`).
- **All events oracle.** A full mixed stream of structural events (duplication, transfer, loss,
  inversion, transposition) matches the array oracle exactly
  (`test_nucleotide_genome.py::test_all_events_match_oracle`).
- **Indel lengths.** Intergenic insertion/deletion run lengths are geometric with mean equal to
  `indel_mean_length` (`test_nucleotide_indels.py::test_draw_indel_length_matches_mean`).

## References

- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a
  phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.
  *Bioinformatics* 36(4): 1286–1288.
