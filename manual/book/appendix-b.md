# Output files

One table per level, listing the files it writes. **Default** says whether a file is written with no
arguments (**yes**), only when you name its token (**no**), or is reachable from Python alone
(**Python**). Trees are Newick, tables and logs TSV, sequences FASTA; branch lengths are time, except
in the sequence phylograms, which are in substitutions per site.

## Where the files go

```
out/run.zombi2              the run report — one page for the whole run
out/species/                species_{complete,extant}.nwk · species_events.tsv · species_fates.tsv
out/genomes/                genome_events.tsv · profiles.tsv · genomes.tsv · initial_genome.tsv
out/genomes/                rearrangement_events.tsv · chromosome_events.tsv  (ordered, nucleotide)
out/genomes/                block_events.tsv                                  (nucleotide)
out/genomes/markers.tsv     one file for the run, not a directory        (zombi2 tools format)
out/genomes/gene_trees/     gene_tree_fam<f>_complete.nwk · …_extant.nwk
out/genomes/homology/       homology_fam<f>.tsv                          (zombi2 tools format)
out/genomes/recphylo/       recphylo_fam<f>.xml                          (zombi2 tools format)
out/sequences/              clock_species_tree_complete.nwk · …_extant.nwk
out/sequences/alignments/   fam<f>.fasta
out/sequences/phylograms/   phylogram_fam<f>_*.nwk
out/sequences/genomes/      genome_<lineage>.fasta                       (nucleotide runs)
out/traits/                 trait_values.tsv · trait_tree.nwk · trait_events.tsv
```

Each level directory also holds that command's log — `species.log`, `genomes.log`, `sequences.log`,
`traits.log` — with the version, the command line and every resolved parameter, rates in their written
form. Beside it, a summary of what came out: `species_summary.json`, `genome_summary.json`,
`sequences_summary.json`, `trait_summary.json`, written at every level but the ordered and nucleotide
genome resolutions. `--flat` writes every file straight into the run directory instead.

## Species trees — `simulate_species_tree`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Complete tree | `species_complete.nwk` | Newick | yes | every lineage, including extinct and unsampled |
| Extant tree | `species_extant.nwk` | Newick | yes | only the sampled survivors |
| Event log | `species_events.tsv` | TSV | yes | every speciation/extinction — `time` · `kind` · `parents` · `children`, one row per event. A `speciation` row is the lineage that ended and its two children, `;`-packed (`n0` → `n1;n2`); an `extinction` row is the dying lineage as the parent with no children. A lineage that died is written `e<id>` |
| Tip fates | `species_fates.tsv` | TSV | yes | each tip's resolved fate — `lineage` · `fate` (`extant` / `extinct` / `unsampled`) |
| Fossils | `species_fossils.tsv` | TSV | yes¹ | sampled fossil lineages — `lineage` · `time` |

¹ written only if fossil sampling recovered any.

## Genomes, family — `simulate_genomes_family`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | the source of truth — `time` · `kind` · `family` · `parents` · `children`, **one row per event**. The participants are gene copies written `n<species>_g<copy>`, `;`-packed where there are two, so the branch each one sat on rides inside the token and no `lineage`, `recipient` or `donor` column is needed: a transfer is one row naming the donor's copy and the copy that arrived. See [One row per event](#one-row-per-event) below |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Genomes | `genomes.tsv` | TSV | yes | every node's gene content, **ancestors included** — `lineage` · `family` · `copy`. **One row per gene copy**, so a lineage holding six genes has six rows; two rows sharing a `family` are two copies of it. `copy` is the same identifier the event log uses, so a gene can be followed from the genome it sits in back to the event that made it. `profiles.tsv` is the same information counted, and only for the extant tips |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `family` · `copy`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Conditioning | `conditioned_on` | text | conditioned | written **only when a rate was conditioned**: the run's levels this run read via `DrivenBy` (one per line, e.g. `traits`). It records the dependency so re-running a driver level (say the trait) refuses to leave this run silently stale, or clears it under `--force`. A run with no driven rate writes no such file |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | each family's true genealogy, in `genomes/gene_trees/`. Leaves are `n<species>_g<copy>`; **internal nodes are labelled `<event>_n<species>`** — `duplication_n45`, `transfer_n45`, `speciation_n45` — the event that ended that gene and the species branch it was on, which is what makes the tree readable on its own. A family with no surviving copy writes no `_extant` file |
| Tip names | `names.tsv` | TSV | external tree | written **only when the tree came from `--from` with its own tip labels** — `node` · `name`, mapping ZOMBI's `n<id>` back to the labels you supplied. Every other output names nodes `n<id>`, so this is the join back to your taxa; `profiles.tsv`'s columns key on the same ids |
| Family origination | `.gene_trees[f].origination` | float | Python | when the family was founded — where its gene tree's root branch begins |

### One row per event

Every event log ZOMBI2 writes opens with `time` and `kind`, then the payload columns for that file,
the same on every row. Here they are `family`, `parents` and `children` — the gene copies the event
ended and the ones it began, `;`-packed. A copy is written
`n<species>_g<copy>`, so the branch it lived on rides inside the token (`e<id>` where that lineage
died) and there is no `lineage`, `recipient` or `donor` column to keep in step with it.

| Kind | `parents` | `children` |
|------|-----------|------------|
| `origination` | — | the founding copy |
| `duplication` | the copy that ends | its two descendants, both on the same branch |
| `speciation` | the copy that ends at the split | one copy in each daughter species |
| `loss` | the copy that ends | — |
| `transfer_additive` | the donor's copy | the continuation on the donor branch, then the copy that arrived |
| `transfer_replacing` | the donor's copy, then the copy it overwrote | the same two |

One row of each kind, from a real run (the file is tab-separated; padded here to line the columns up):

```
time                 kind                family  parents        children
0.0                  origination         0                      n0_g0
0.281867             speciation          14      n0_g14         n1_g15;n2_g29
0.2959603504961474   duplication         4       n1_g19         n1_g43;n1_g44
0.11415202408218499  loss                0       n0_g0
1.2866954826434076   transfer_additive   8       n1_g23         n1_g255;n5_g256
0.36393613420080373  transfer_replacing  1       n2_g30;n1_g16  n2_g47;n1_g48
```

A transfer's children read **donor first, recipient second**, so the row alone says which way the
material went: the `transfer_additive` above left `n1` and landed on `n5`. A `transfer_replacing`
carries two parents — the donor's copy, then the copy it overwrote on the recipient branch — and
writes **no separate `loss` row** for the copy it displaced. That death is what the kind means.

One row is one event, so counting rows by kind counts events. `genome_summary.json` counts the
biology rather than the rows, and differs in three places: it reports `transfer` as one number over
both kinds, counts each replacing transfer's displaced copy under `loss`, and separates the families
the run started with (`initial`) from the ones the origination rate made (`origination`), where the
log writes both as `origination` rows and the initial ones sit at time 0.

It also reports `empty_genomes`: the extant genomes that came out with no genes at all. There is no
floor at this resolution — loss is counted per copy, and the last copy is a copy like any other — so
a lineage can lose everything, and an empty genome is otherwise invisible, having no row in
`profiles.tsv` and no gene tree.

Speciation is the largest single kind in the file, and the rows are not redundant: a gene tree's
internal nodes are labelled `speciation_n14` — the kind and the branch, no copy id — so this log is
the only record of the internal gene copies and their parentage. In Python an `Event` is still one
gene-tree **edge**; `zombi2.genomes.events.events_from_tsv` expands each row back into one per edge,
which is what the gene trees are derived from.

## Genomes, ordered — `simulate_genomes_ordered`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the gene genealogy with the place each event happened.** The family resolution's five columns, then five more — `time` · `kind` · `family` · `parents` · `children` · `chromosome` · `start` · `length` · `dest_chromosome` · `dest_position`¹. The kinds and the copy tokens are [the shared ones](#one-row-per-event); the arc is in that branch's own genome just before the event, as `gene_order.tsv` numbers it, and `chromosome` is a chromosome **id**, the one `chromosome_events.tsv` writes `c<id>`. `dest_position` is where the material landed — the tandem block for a duplication, the arriving block for a transfer — and `dest_chromosome` is the recipient's chromosome, the recipient *branch* already being inside the arriving copy's token. So one row carries a whole transfer, both ends of it. A speciation copies a genome whole and leaves all five coordinate cells empty |
| Rearrangements | `rearrangement_events.tsv` | TSV | yes | every inversion, transposition and translocation — `time` · `kind` · `lineage` · `chromosome` · `start` · `length` · `dest_chromosome` · `dest_position` · `flipped`¹. There is no token of its own: `events` writes both this and the log above, a run doing two different things to a genome. Rearrangements begin and end no gene lineage, so they have no `parents` and no `children` and a file of their own; a segment has no id either, which is why this is the one event log that still names its branch in a column. `dest_chromosome` is set by a translocation alone (a transposition lands on the chromosome it left) and `flipped` by the two that move a segment, which may land inverted |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Gene order | `gene_order.tsv` | TSV | yes | signed gene order of **every node**, ancestors included — `lineage` · `chromosome` · `position` · `strand` · `family` · `copy` |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `strand` · `family` · `copy`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Conditioning | `conditioned_on` | text | conditioned | as at the family resolution: written **only when a rate, an extent or `transfer_to` was conditioned** — the run's levels this run read via `DrivenBy`, one per line |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — `time` · `kind` · `parents` · `children`, chromosomes written `n<species>_c<id>` so the branch rides in the token here too. Kinds are `initial` (a replicon the run starts with, at time 0), `speciation`, `fission`, `fusion`, `origination` and `loss`; a fusion is the one row with two parents |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | as at the family resolution — position is orthogonal to genealogy |

¹ a segment is named by `start` (its first position, in the chromosome's frame just before the event)
and `length` (how many genes it covered), counted rightwards from `start` and **wrapping past position
0 on a circular chromosome** — so `start + length` greater than the chromosome's gene count means the
segment crossed the origin. `dest_position` is an index into what was left after it was excised. A
segmental event acts on a run of genes of several families at once and this table has one `family`
column, so `genome_events.tsv` writes one row per gene lineage, each repeating the same arc.

## Genomes, nucleotide — `simulate_genomes_nucleotide`

From `zombi2 genomes --resolution nucleotide` or `result.write(dir, outputs=[...])`.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the genealogy, in the format every resolution writes** — the same `time` · `kind` · `family` · `parents` · `children` as at the family resolution, [one row per event](#one-row-per-event). So one reader serves all three resolutions and `zombi2 tools` works here unchanged. `family` is the declared gene (else the recovered root-block), and the copy in each token is a **gene** id — the one the gene trees, alignments and homology tables use. It is derived onto the root-block partition, where a copy either covers a block in full or does not touch it, which is what makes a duplication a bifurcation here. A nucleotide transfer is always additive, so `transfer_replacing` cannot appear |
| Block events | `block_events.tsv` | TSV | yes | **this resolution's own record**, which has no counterpart elsewhere: the copy-lineage log over ancestral intervals — `time` · `kind` · `parents` · `children` · `chromosome` · `source` · `start` · `end`. Copy lineages are written `n<species>_g<copy>` and the chromosome `n<species>_c<id>`, as everywhere else. One row per **ancestral interval** an event touched, so an event spanning several blocks writes several rows sharing a `time` and `kind` — and a duplication mints a child *without* ending its parent, because a copy lineage covers an extent and the event covers a sub-extent. That is the right model for sequence and the wrong shape for a gene tree, which is why the genealogy above is a separate table rather than this one renamed. Kinds are `initial` · `origination` · `duplication` · `loss` · `transfer` · `speciation`: `initial` is the initial genome being laid down at time 0 (one row per replicon), kept distinct from `origination`, the de-novo births the `--origination` rate makes. A speciation re-mints one copy lineage into one per daughter without touching sequence, so it has no interval and writes a single row with both daughters and every coordinate cell empty; a transfer's blocks keep their source coordinates and leave `chromosome` empty. `source`/`start`/`end` are **ancestral** coordinates — the physical ones are a different frame and live in the rearrangement log. This is what `read_nucleotide_genomes` replays |
| Rearrangements | `rearrangement_events.tsv` | TSV | yes | every inversion, transposition and translocation, in **physical** coordinates along the chromosome as `blocks.tsv` numbers it — `time` · `kind` · `lineage` · `chromosome` · `start` · `length` · `dest_chromosome` · `dest_position` · `flipped`, the same columns the ordered resolution writes with base pairs in place of genes. No token of its own either: `events` writes all three tables here. They begin and end no lineage, so `parents` and `children` would be empty on every row and they get a file of their own. A nucleotide translocation records no `dest_position`: its blocks keep their source coordinates and the engine places the arc |
| Blocks | `blocks.tsv` | TSV | yes | every node's genome as its block mosaic, ancestors included — `lineage` · `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. The rows of one chromosome tile it end to end from 0. The largest file this level writes: blocks are not kept maximal during a run, so it grows with their number × every node |
| Genes | `genes.tsv` | TSV | yes | the declared genes in initial coordinates — `family` · `name` · `source` · `start` · `end` · `strand` (the **coding** strand). Header-only when none were declared |
| Initial sequence | `initial_sequence.fasta` | FASTA | yes¹ | the initial DNA the run was given (`--fasta`), one `>source<n>` record per replicon. Written only when a FASTA was supplied; it is what lets a separate `zombi2 sequences` run found its blocks from the real sequence |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Conditioning | `conditioned_on` | text | conditioned | as at the family resolution: written **only when a rate, an extent or `transfer_to` was conditioned** — the run's levels this run read via `DrivenBy`, one per line |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — same format and same kinds as ordered |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | one tree per declared gene (else per recovered root-block), in `gene_trees/` |
| GFF | `genome_<lineage>.gff` | GFF3 | yes | one file per node, in `gff/`: that genome's **genes** in its own coordinates — the annotation to read beside the sequence level's `genome_<lineage>.fasta`, which it names its sequences to match |
| BED | `genome_<lineage>.bed` | BED | yes | one file per node, in `bed/`: that genome's **blocks**, spacer included, each named by the ancestral interval it descends from — the ancestry as a browser track |

Three event files, because a nucleotide run records three different things: the genealogy in the one
format every resolution writes, the interval record only this resolution has, and what moved without
beginning or ending anything. Ancestral coordinates are in `block_events.tsv` and physical ones in
`rearrangement_events.tsv` — two frames, so two files, and no row carrying columns from both.

The events index against the species tree canonicalised so its `n<id>` labels match the branch named
inside each copy token, so a genomes run needs that exact tree to be replayable. A run grown by `zombi2 species`
already has it; a run whose tree came from `--from` gets a copy written to its own
`species/species_complete.nwk`, rather than a second file under another name. Either way `zombi2
sequences` can replay the gene genealogy from the run directory alone. Like `names.tsv` (external
input trees) and the `.log`, that copy is a CLI artifact, not a `result.write()` output.

## Sequences — `simulate_sequences`

The `zombi2 sequences` command replays a prior `zombi2 genomes` run — its own run directory, or `--from` another —
its species tree and its `genome_events.tsv`. Gene outputs are written **one file per gene
family** (`<f>` = family number); a family with no surviving copy writes none. Every node is labelled
`n<species>_g<copy>`, so a phylogram's tips pair with its alignment and its internal nodes with the
ancestral sequences.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Alignments | `fam<f>.fasta` | FASTA | yes | one row per extant gene copy — nucleotides or amino acids, following the model. They go in `alignments/`, which is what lets the name be this short |
| Phylograms | `phylogram_fam<f>_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the gene tree each family's sequences were drawn along, in `phylograms/`. Under `+Γ`/`+I` the branch length is the **mean** over sites, which is what the rate classes are normalised to |
| Ancestral | `sequences_ancestral_fam<f>.fasta` | FASTA | no | the sequence at every node that is not an extant tip: internal nodes, and the tips where a copy was lost or its species died. One per family, so they go in `ancestral/` |
| Founding | `sequences_founding.fasta` | FASTA | no | one record `fam<f>` per family — the sequence it originated with, where its phylogram's root branch begins |
| Clock species tree | `clock_species_tree_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the species tree with its branches in substitutions/site — the molecular clock made visible. The mean over sites under `+Γ`/`+I`, as for the phylograms. A driven substitution rate shows here too: a branch is the rate times the driver integrated along it, so this is where you read what the trait did |
| Conditioning | `conditioned_on` | text | conditioned | written **only when the substitution rate was conditioned**: the run's levels this run read via `DrivenBy` (one per line, e.g. `traits`). It records the dependency so re-running the trait refuses to leave this run silently stale, or clears it under `--force`. A run with no driven rate writes no such file |
| Genomes | `genome_<lineage>.fasta` | FASTA | yes | one file per **node** of the complete tree — extant, extinct and ancestral alike — one record `<lineage>_chr<c>` per chromosome: the assembled genome, its blocks concatenated in physical order, in `genomes/`. **Nucleotide genome runs only**: a family or ordered run has gene families, not coordinates, so there is nothing to lay out, and no `genomes/` is created. The biggest thing this level writes — a whole genome times every node |
| Initial genome | `genome_initial.fasta` | FASTA | yes | the genome the run **started** with, as sequence — the state the stem leads *from*, which is not any node's. In `genomes/` with the rest, being a whole-genome FASTA like they are. Nucleotide runs only |

On a **nucleotide** genome run every block evolves, spacer as well as gene, so a genome of *b* blocks
writes *b* alignments and *b* phylograms — that is what makes the genomes assemblable. The number in
those filenames is then a **root block index**, not a gene family id, and the files say so: `block6.fasta`
and `phylogram_block6_complete.nwk` in place of `fam6.…`. The two numbering schemes are different, so go
from a gene to its block with `genomes.block_of(family)` (Ch7). `zombi2 sequences` reads a nucleotide
handoff too — it recognises one by its `blocks.tsv` — so all of this is reachable from the command line.

## Joint — `simulate_joint` / `zombi2 joint`

A joint run grows two levels at once, so it writes both, each in the format its own command would
give it: the species files, and then the driver's — the trait's or the genomes'. There is no output
of its own beyond the run log.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Species tree | `species/species_complete.nwk` · `…_extant.nwk` · `species_events.tsv` | Newick, TSV | yes | the grown tree — complete, so the extinct lineages whose fate the driver decided are kept |
| The trait it grew | `traits/trait_values.tsv` · `trait_events.tsv` · `trait_tree.nwk` | TSV, Newick | yes¹ | as the traits level writes them |
| The genomes it grew | `genomes/genome_events.tsv` · `profiles.tsv` · `genomes.tsv` · `gene_trees/` | TSV, Newick | yes¹ | as the genomes level writes them |
| Run log | `species/joint.log` | TSV | yes | the resolved parameters, as every command writes |

¹ whichever driver the run used — one per run, never both.

## Traits — `simulate_continuous` / `simulate_discrete`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Values | `trait_values.tsv` | TSV | yes | value at every node (tips, extinct, internal) — `node` · `kind` · `trait`, where `kind` is the tip's fate (`extant` / `extinct` / `unsampled`) or `ancestor`, so `kind == "extant"` isolates the observed tips |
| Events | `trait_events.tsv` | TSV | yes (discrete) | the trait's whole history — an `initial` row giving the state at t=0, then every switch: `time` · `kind` · `lineage` · `from` · `to`, where `kind` is `initial` · `on_branch` · `on_speciation`. The one event log whose payload is a **state change** rather than a birth and a death, so it keeps its `lineage` and has no `parents` / `children`. Times are full precision (they drive a conditioned run's Gillespie). **This is also the conditioning file**: a genome/sequence run drives a rate with `mod.DrivenBy("trait_events.tsv", …)`, replaying it against the shared tree. A continuous trait carries only the `initial` row and any `at_speciation` jumps (a diffusion can't be rebuilt from events), and that holds for a multi-optimum (`regimes=`) run too. The exception is a **correlated** multi-trait run, which carries no event log at all: its value is one number per trait, and the single `from` / `to` columns do not hold a vector |
| Trait tree | `trait_tree.nwk` | Newick | no | tree with every node annotated `[&trait=…]` (opens in FigTree / iTOL) |
| Summary | `trait_summary.json` | JSON | yes | what came out, not what was asked for — `tips` · `nodes` · `events` (the `on_branch` and `on_speciation` counts), then `states` · `most_common_share` for a discrete trait, or `values` (min/mean/max) · `value_at_root_node` for a continuous one. The root node sits at the end of the stem, so that value is not the one the run started from |

## Conditioning and joining — no new files

Neither adds a format. A **conditioned** run writes the target level's files plus the **driver file**
it read (above), keeping the pairing on disk; a **joint** run writes **both** levels, each in its own
format.

The `zombi2 tools` commands write their own files — the homology matrix and the reconciliation/scoring
outputs — catalogued in Appendix C.
