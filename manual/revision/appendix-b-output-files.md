# Appendix B — Output files

Every `simulate_*` returns a result; `result.write("out/", outputs=[...])` writes the files, and
omitting `outputs` writes the **default** set. Trees are Newick, tables and logs are TSV, sequences are
FASTA. The **Default** column says whether a file is written with no arguments (**yes**), only when you
name its token (**no**), or is available in Python but has no file yet (**Python**).

## Species trees — `simulate_species_tree`

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Complete tree | `species_complete.nwk` | Newick | yes | every lineage, including extinct and unsampled |
| Extant tree | `species_extant.nwk` | Newick | yes | only the sampled survivors |
| Event log | `species_events.tsv` | TSV | yes | every speciation/extinction — `time · kind · lineage · children` |
| Fossils | `species_fossils.tsv` | TSV | yes¹ | sampled fossil lineages — `lineage · time` |

¹ written only if fossil sampling recovered any.

## Genomes, unordered — `simulate_genomes_unordered`

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Event log | `genome_events.tsv` | TSV | yes | the source of truth — `time · kind · lineage · family · copy · parent · recipient` |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Gene trees | `.gene_trees` (`GeneTree.to_newick()`) | Newick | Python | each family's true genealogy (`.complete` and `.extant`) |

## Genomes, ordered — `simulate_genomes_ordered`

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Event log | `genome_events.tsv` | TSV | yes | as unordered |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Gene order | `gene_order.tsv` | TSV | yes | signed gene order of each leaf — `species · chromosome · position · strand · family · gene` |
| Rearrangements | `rearrangements.tsv` | TSV | no | inversions/transpositions/translocations — `time · kind · lineage · chromosome · start · length · dest_chromosome · dest_position · flipped` |
| Chromosome events | `chromosome_events.tsv` | TSV | no | chromosome-network edges — `time · kind · lineage · parents · children` |
| Gene trees | `.gene_trees` (`GeneTree.to_newick()`) | Newick | Python | as unordered |

## Sequences — `simulate_sequences`

One file per gene family (`<f>` = family number); a family with no surviving copy writes none.

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Alignments | `sequences_alignment_fam<f>.fasta` | FASTA | yes | one row per extant gene copy (labels `g<copy>_n<species>`) |
| Ancestral | `sequences_ancestral_fam<f>.fasta` | FASTA | no | reconstructed sequence at every internal gene-tree node |
| Phylograms | — | Newick | soon | gene trees in substitutions/site — forthcoming |

## Traits — `simulate_continuous` / `simulate_discrete`

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Values | `trait_values.tsv` | TSV | yes | value at each extant tip — `node · trait` |
| Changes | `trait_changes.tsv` | TSV | no | realized transitions per branch — `time · kind · lineage · from · to` (header-only for a continuous trait) |
| Trait tree | `trait_tree.nwk` | Newick | no | tree with every node annotated `[&trait=…]` (opens in FigTree / iTOL) |

## Tools

The scoring and reconciliation commands (tree-distance, reconciliation-accuracy, the undated simulator,
the ALE likelihood) are quarantined during the clean-core rebuild; their files are documented here as
each returns.
