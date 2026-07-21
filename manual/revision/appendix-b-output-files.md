# Appendix B — Output files

Every `simulate_*` returns a result; `result.write("out/", outputs=[...])` writes the files, and
omitting `outputs` writes the **default** set. Trees are Newick, tables and logs are TSV, sequences are
FASTA. Tree branch lengths are **time** everywhere except the sequence phylograms, whose lengths are in
**substitutions per site**. The **Default** column says whether a file is written with no arguments
(**yes**), only when you name its token (**no**), or is available in Python but has no file yet
(**Python**).

Every `zombi2` **command** also writes a run log next to its outputs (`species.log`, `genomes.log`,
`sequences.log`, `traits.log`): the version, the timestamp, the command line, and every resolved
parameter. Rates are recorded in their **written form** — `birth<TAB>1.0 * OnTime({0: 1, 3: 0.3})` —
so a line pastes straight back into the flag or a `--params` file. It is a CLI artifact, not a
`result.write()` output, so it has no row in the tables below.

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

The `zombi2 genomes` **command** also writes `genome_species_tree.nwk` — the complete species tree
canonicalised so its `n<id>` labels match the event log's `lineage` column — so `zombi2 sequences
--genomes DIR` can replay the gene genealogy from that directory alone. Like `names.tsv` (external
input trees) and the `.log`, it is a CLI artifact, not a `result.write()` output.

## Sequences — `simulate_sequences`

The `zombi2 sequences` command replays a prior `zombi2 genomes` output directory (`--genomes DIR`) —
its `genome_species_tree.nwk` and `genome_events.tsv`. Gene outputs are written **one file per gene
family** (`<f>` = family number); a family with no surviving copy writes none. Every node is labelled
`g<copy>`, so a phylogram's tips pair with its alignment and its internal nodes with the ancestral
sequences.

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Alignments | `sequences_alignment_fam<f>.fasta` | FASTA | yes | one row per extant gene copy — nucleotides or amino acids, following the model |
| Phylograms | `sequences_phylogram_fam<f>_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the gene tree each family's sequences were drawn along |
| Ancestral | `sequences_ancestral_fam<f>.fasta` | FASTA | no | reconstructed sequence at every internal node |
| Species phylogram | `sequences_species_phylogram_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | no | the species tree scaled by the molecular clock |

## Traits — `simulate_continuous` / `simulate_discrete`

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Values | `trait_values.tsv` | TSV | yes | value at each extant tip — `node · trait` |
| Changes | `trait_changes.tsv` | TSV | no | realized changes — `time · kind · lineage · from · to`, where `kind` is `on_branch` or `on_speciation` (a continuous trait diffuses, so it logs only its `at_speciation` jumps: header-only without them) |
| Trait tree | `trait_tree.nwk` | Newick | no | tree with every node annotated `[&trait=…]` (opens in FigTree / iTOL) |
| Driver | `trait_driver.tsv` | TSV | no | the conditioning driver file a genome/sequence run reads via `mod.DrivenBy(...)` (a discrete trait only) |

## Coupling — no new files

Coupling adds no formats. A **conditioned** run writes the target level's files plus the **driver file**
it read (above), keeping the pairing on disk; a **joint** run writes **both** levels, each in its own
format.

## Tools

The scoring and reconciliation commands (tree-distance, reconciliation-accuracy, the undated simulator,
the ALE likelihood) are quarantined during the clean-core rebuild; their files are documented here as
each returns.
