# Output files

One table per level, listing the files it writes. **Default** says whether a file is written with no
arguments (**yes**), only when you name its token (**no**), or is reachable from Python alone
(**Python**). Trees are Newick, tables and logs TSV, sequences FASTA; branch lengths are time, except
in the sequence phylograms, which are in substitutions per site.

**Columns may be added; read them by name.** A later version may add a column to any of these tables —
a new kind of event needs somewhere to say what it did. What will not change is an existing column's
name or its meaning, and no column is removed within a major version. So a script that reads a column
by its header keeps working across upgrades, and one that reads by position, or assumes a column
count, does not. New columns go at the **end** of the row. The same file already differs by
resolution — `genome_events.tsv` has five columns at the family resolution and ten at the ordered
one — so matching by name is what these tables have always required.

## Where the files go

```
out/run.zombi2              the run report — one page for the whole run
out/species/                species_{complete,extant}.nwk · species_events.tsv · species_fates.tsv
out/genomes/                genome_events.tsv · initial_genome.tsv
out/genomes/                genomes.tsv                                       (family)
out/genomes/                profiles.tsv                                      (family, ordered)
out/genomes/                rearrangement_events.tsv · chromosome_events.tsv  (ordered, nucleotide)
out/genomes/                gene_order.tsv                                    (ordered)
out/genomes/                block_events.tsv · blocks.tsv · genes.tsv         (nucleotide)
out/genomes/markers.tsv     one file for the run, not a directory        (zombi2 tools format)
out/genomes/gene_trees/     gene_tree_fam<f>_complete.nwk · …_extant.nwk
out/genomes/gff/            genome_<lineage>.gff                              (nucleotide)
out/genomes/bed/            genome_<lineage>.bed                              (nucleotide)
out/genomes/homology/       homology_fam<f>.tsv                          (zombi2 tools format)
out/genomes/recphylo/       recphylo_fam<f>.xml                          (zombi2 tools format)
out/sequences/              clock_species_tree_complete.nwk · …_extant.nwk
out/sequences/alignments/   fam<f>.fasta
out/sequences/phylograms/   phylogram_fam<f>_*.nwk
out/sequences/genomes/      genome_<lineage>.fasta                       (nucleotide runs)
out/traits/                 trait_values.tsv · trait_tree.nwk · trait_events.tsv
out/traits/<name>/          the same files, when --name was given
```

Each level directory also holds that command's log — `species.log`, `genomes.log`, `sequences.log`,
`traits.log` — with the version, the command line and every resolved parameter, rates in their written
form. Beside it, a summary of what came out: `species_summary.json`, `genome_summary.json`,
`sequences_summary.json`, `trait_summary.json` — written at every level, and at all three genome
resolutions, except by `zombi2 genomes --stream`, which writes each family straight to disk and never
holds the whole run to summarise. A streamed run also keeps its own `species_complete.nwk` beside its
tables. It still reads the tree from `species/`, as any genomes run does; the copy is so that the
`genomes/` directory stands alone as the handoff a later level or another tool reads.

`--flat` writes every file straight into the run directory instead, except two that a flat run does
not write at all: `run.zombi2` and `conditioned_on`. Both are built from a grouped run's level
directories, which a flat run does not have.

Every directory a run fills with one file per family or per node — `gene_trees/`, `gff/` and `bed/` at
the genome level, `alignments/`, `ancestral/`, `phylograms/` and `genomes/` at the sequence level — is
**emptied before that run fills it**, so its contents describe one run and nothing else. Two do not:
`--flat`, which empties nothing, and the `homology/` and `recphylo/` directories `zombi2 tools format`
writes, which are left as they are.

## Species trees: `simulate_species_tree`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Complete tree | `species_complete.nwk` | Newick | yes | every lineage, including extinct and unsampled |
| Extant tree | `species_extant.nwk` | Newick | yes | only the sampled survivors |
| Event log | `species_events.tsv` | TSV | yes | every speciation/extinction — `time` · `kind` · `parents` · `children`, one row per event. A `speciation` row is the lineage that ended and its two children, `;`-packed (`n0` → `n1;n2`); an `extinction` row is the dying lineage as the parent with no children. A lineage that died is written `e<id>` |
| Tip fates | `species_fates.tsv` | TSV | yes | each tip's resolved fate — `lineage` · `fate` (`extant` / `extinct` / `unsampled`) |
| Fossils | `species_fossils.tsv` | TSV | yes¹ | sampled fossil lineages — `lineage` · `time` |
| Summary | `species_summary.json` | JSON | yes | what the run produced — counts by fate, tree height, stem, total branch length, the realised birth and death rates, and the seed |

¹ written only if fossil sampling recovered any.

## Genomes, family: `simulate_genomes_family`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | the run's true history — `time` · `kind` · `family` · `parents` · `children`, **one row per event**. See [One row per event](#one-row-per-event) below |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Genomes | `genomes.tsv` | TSV | yes | every node's gene content, ancestors included — `lineage` · `family` · `copy`, **one row per gene copy**. `copy` is the identifier the event log uses, so a gene can be traced back to the event that made it. `profiles.tsv` is the same information counted, for the extant tips only |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `family` · `copy`. It has no `lineage` column because it is no node's genome: every node sits at the *end* of its branch |
| Conditioning | `conditioned_on` | text | conditioned | the levels this run reads as a driver, one per line, whether for a rate or for `transfer_to`. Written only when something was conditioned. Re-running a driver level then refuses rather than leaving this run stale, unless you pass `--force` |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | each family's true genealogy, in `genomes/gene_trees/`. Leaves are `n<species>_g<copy>`; internal nodes are labelled `<event>_n<species>` (`duplication_n45`, `transfer_n45`), naming the event that ended that gene and the branch it was on. A family with no surviving copy writes no `_extant` file |
| Tip names | `names.tsv` | TSV | external tree | `node` · `name`, mapping ZOMBI2's `n<id>` back to the labels you supplied. Written only when the tree came from `--from` with its own tip labels; it is the join from every other output back to your taxa |
| Summary | `genome_summary.json` | JSON | yes | what the run produced — events by kind, families born, surviving and died out, genes per genome, `empty_genomes`, whether the family-size cap bound, and the seed |
| Species tree | `species_complete.nwk` | Newick | yes (Python) | the tree the run evolved along. Every other file here is indexed by its node labels, so the directory is not readable — by anyone or by `genomes.read_run()` — without it. `result.write()` writes it by default, so a directory written from Python stands alone; `zombi2 genomes` leaves it out — except under `--stream`, which is family-only — because a run already keeps one copy at `species/species_complete.nwk`, shared by every level |
| Tip fates | `species/species_fates.tsv` | TSV | external tree | each tip's resolved fate — `lineage` · `fate`, in the format the species level writes. Written beside the copied tree, only when the tree came from `--from` |
| Family origination | `.gene_trees[f].origination` | float | Python | when the family was founded — where its gene tree's root branch begins |
| Driver views | `.presence(name)` · `.completion(name)` | driver | Python | a named family's presence (`present` / `absent`) and a declared module's completion (a fraction) along every lineage, for use as a driver (Ch9). Read off the families' gene trees, so they change *inside* a branch; the **ordered** resolution gives the same two |

### One row per event

Every event log opens with `time` and `kind`, then the payload columns for that file, the same on
every row. Here they are `family`, `parents` and `children`: the gene copies the event ended and the
ones it began, `;`-packed.

A gene copy is written `n<species>_g<copy>` — copy 30 on branch `n2` is `n2_g30`, and `e2_g30` if
that lineage went extinct. The branch is part of the name, so the file needs no `lineage`,
`recipient` or `donor` column: a transfer is one row naming the donor's copy and the copy that
arrived.

| Kind | `parents` | `children` |
|------|-----------|------------|
| `origination` | — | the founding copy |
| `duplication` | the copy that ends | its two descendants, both on the same branch |
| `speciation` | the copy that ends at the split | one copy in each daughter species |
| `loss` | the copy that ends | — |
| `transfer_additive` | the donor's copy | the continuation on the donor branch, then the copy that arrived |
| `transfer_replacing` | the donor's copy, then the copy it overwrote | the same two |

One row of each kind, from a real run:

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
biology instead, and differs in three places: it reports `transfer` as one number over both kinds, it
counts each replacing transfer's displaced copy under `loss`, and it separates the families the run
started with (`initial`) from the ones the origination rate made. It also reports `empty_genomes`,
the extant genomes that came out with no genes at all, which are easy to miss: their column of
`profiles.tsv` is all zeros, they have no row in `genomes.tsv`, and no `_extant` gene tree names them.

Speciation is the largest kind in the file, and those rows are not redundant. A gene tree labels its
internal nodes `speciation_n14`, with no copy id, so this log is the only record of the internal gene
copies and their parentage. In Python an `Event` is one event; the per-edge record is a `GeneEdge`,
and `zombi2.genomes.events.edges_from_tsv` expands each row into one per edge, which is what the gene
trees are built from.

## Genomes, ordered: `simulate_genomes_ordered`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | the gene genealogy, with the place each event happened. The family resolution's five columns plus five more — `time` · `kind` · `family` · `parents` · `children` · `chromosome` · `start` · `length` · `dest_chromosome` · `dest_position`¹. Kinds and copy names are [the shared ones](#one-row-per-event). Coordinates are in that branch's own genome just before the event, as `gene_order.tsv` numbers it. `dest_position` and `dest_chromosome` are where the material landed, so one row carries both ends of a transfer. A speciation copies a genome whole and leaves all five coordinate cells empty |
| Rearrangements | `rearrangement_events.tsv` | TSV | yes | every inversion, transposition and translocation — `time` · `kind` · `lineage` · `chromosome` · `start` · `length` · `dest_chromosome` · `dest_position` · `flipped`¹. They begin and end no gene lineage, so they have no `parents` and `children` and get a file of their own; a segment has no name either, which is why this is the one event log that still puts its branch in a column. `dest_chromosome` is set only by a translocation, and `flipped` by the two events that move a segment, which may land inverted. The `events` token writes this file and the log above together |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Gene order | `gene_order.tsv` | TSV | yes | signed gene order of every node, ancestors included — `lineage` · `chromosome` · `topology` · `position` · `strand` · `family` · `copy`. `topology` is `circular` or `linear`, written beside every gene: it decides where a segmental event stops and which chromosomes may fuse. A chromosome with no genes has no rows here, and so no topology |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `topology` · `position` · `strand` · `family` · `copy`. No `lineage` column, for the reason given at the family resolution |
| Conditioning | `conditioned_on` | text | conditioned | as at the family resolution, and written when a rate or `transfer_to` was conditioned |
| Tip names | `names.tsv` | TSV | external tree | as at the family resolution |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — `time` · `kind` · `parents` · `children`, chromosomes named `n<species>_c<id>` on the same pattern as gene copies. Kinds are `initial` (a replicon the run starts with, at time 0), `speciation`, `fission`, `fusion`, `origination` and `loss`; a fusion is the one row with two parents |
| Summary | `genome_summary.json` | JSON | yes | `events` counted as biology rather than rows (see the family resolution), `families` born/surviving/died_out, genes and chromosomes per genome, and the `rearrangements` and `chromosome_events` by kind |
| Species tree | `species_complete.nwk` | Newick | yes (Python) | as at the family resolution |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | as at the family resolution — position is orthogonal to genealogy |

¹ a segment is named by `start` (its first position, in the chromosome's frame just before the event)
and `length` (how many genes it covered), counted rightwards from `start` and **wrapping past position
0 on a circular chromosome** — so `start + length` greater than the chromosome's gene count means the
segment crossed the origin. `dest_position` is an index into what was left after it was excised. A
segmental event acts on a run of genes of several families at once and this table has one `family`
column, so `genome_events.tsv` writes one row per gene lineage, each repeating the same arc.

## Genomes, nucleotide: `simulate_genomes_nucleotide`

From `zombi2 genomes --resolution nucleotide` or `result.write(dir, outputs=[...])`.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the genealogy, in the format every resolution writes** — the same `time` · `kind` · `family` · `parents` · `children` as at the family resolution, [one row per event](#one-row-per-event). So one reader serves all three resolutions and `zombi2 tools` works here unchanged. `family` is the declared gene (else the recovered root-block), and the copy in each token is a **gene** id — the one the gene trees, alignments and homology tables use. It is derived onto the root-block partition, where a copy either covers a block in full or does not touch it, which is what makes a duplication a bifurcation here. A nucleotide transfer is always additive, so `transfer_replacing` cannot appear |
| Block events | `block_events.tsv` | TSV | yes | this resolution's own record, with no counterpart elsewhere — `time` · `kind` · `parents` · `children` · `chromosome` · `source` · `start` · `end`. Copy lineages and chromosomes are named as everywhere else. There is **one row per ancestral interval an event touched**, not one row per event, so an event covering several blocks writes several rows sharing a `time` and `kind`, and a duplication starts a child without ending its parent (the parent still covers the rest of its extent). That is why the genealogy above is a separate table: this shape suits sequence and not a gene tree. Kinds are `initial` · `origination` · `duplication` · `loss` · `transfer` · `speciation`, where `initial` is the starting genome laid down at time 0, one row per replicon. A speciation splits one copy lineage into one per daughter without touching sequence, so it writes a single row with both daughters and every coordinate cell empty; a transfer's blocks keep their source coordinates and leave `chromosome` empty. `source`/`start`/`end` are **ancestral** coordinates; the physical ones live in the rearrangement log. This is what `read_nucleotide_genomes` replays |
| Rearrangements | `rearrangement_events.tsv` | TSV | yes | every inversion, transposition and translocation, in **physical** coordinates along the chromosome as `blocks.tsv` numbers it — `time` · `kind` · `lineage` · `chromosome` · `start` · `length` · `dest_chromosome` · `dest_position` · `flipped`, the same columns the ordered resolution writes with base pairs in place of genes. No token of its own either: `events` writes all three tables here. They begin and end no lineage, so `parents` and `children` would be empty on every row and they get a file of their own. A nucleotide translocation records no `dest_position`: its blocks keep their source coordinates and the engine places the arc |
| Blocks | `blocks.tsv` | TSV | yes | every node's genome as its block mosaic, ancestors included — `lineage` · `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. The rows of one chromosome tile it end to end from 0. The largest file this level writes: blocks are not kept maximal during a run, so it grows with their number × every node |
| Summary | `genome_summary.json` | JSON | yes | what came out: `events` counted as *biology* rather than rows (see the ordered table), `declared_genes`, and base pairs and chromosomes per genome. No phyletic profiles here — the unit is the base pair |
| Species tree | `species_complete.nwk` | Newick | yes (Python) | as at the family resolution |
| Genes | `genes.tsv` | TSV | yes | the declared genes in initial coordinates — `family` · `name` · `source` · `start` · `end` · `strand` (the **coding** strand). Header-only when none were declared |
| Initial sequence | `initial_sequence.fasta` | FASTA | yes¹ | the initial DNA the run was given (`--fasta`), one `>source<n>` record per replicon. Written only when a FASTA was supplied; it is what lets a separate `zombi2 sequences` run found its blocks from the real sequence |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Conditioning | `conditioned_on` | text | conditioned | as at the family resolution: written **only when a rate or `transfer_to` was conditioned** — the levels this run read as a driver, one per line |
| Tip names | `names.tsv` | TSV | external tree | as at the family resolution |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — same format and same kinds as ordered |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | one tree per declared gene (else per recovered root-block), in `gene_trees/` |
| GFF | `genome_<lineage>.gff` | GFF3 | yes | one file per node, plus `genome_initial.gff` for the genome the run **started** with, in `gff/`: that genome's **genes** in its own coordinates — the annotation to read beside the sequence level's `genome_<lineage>.fasta`, which it names its sequences to match |
| BED | `genome_<lineage>.bed` | BED | yes | one file per node, plus `genome_initial.bed` for the genome the run **started** with, in `bed/`: that genome's **blocks**, spacer included, each named by the ancestral interval it descends from — the ancestry as a browser track |
| Driver views | `.presence(name)` · `.completion(name)` | driver | Python | as at the family resolution, over **declared genes**: the name is the GFF's `ID` / `Name`, and what takes a gene away is an arc of DNA rather than a whole copy |

Three event files, because a nucleotide run records three different things: the genealogy in the one
format every resolution writes, the interval record only this resolution has, and what moved without
beginning or ending anything. Ancestral coordinates are in `block_events.tsv` and physical ones in
`rearrangement_events.tsv` — two frames, so two files, and no row carrying columns from both.

The events index against the species tree canonicalised so its `n<id>` labels match the branch named
inside each copy token, so a genomes run needs that exact tree to be replayable. A run grown by `zombi2 species`
already has it; a run whose tree came from `--from` gets a copy written to its own
`species/species_complete.nwk`, rather than a second file under another name, with
`species_fates.tsv` beside it so a later level on this run reads each tip's fate from the record
instead of guessing it from tip depth. Either way `zombi2 sequences` can replay the gene genealogy
from the run directory alone. Like `names.tsv` (external input trees) and the `.log`, those copies
are CLI artifacts, not `result.write()` outputs.

## Sequences: `simulate_sequences`

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
| Conditioning | `conditioned_on` | text | conditioned | written **only when the substitution rate was conditioned**: the levels this run read as a driver (one per line, e.g. `traits`). It records the dependency so re-running the trait refuses to leave this run silently stale, or clears it under `--force`. A run with no driven rate writes no such file |
| Genomes | `genome_<lineage>.fasta` | FASTA | yes | one file per **node** of the complete tree — extant, extinct and ancestral alike — one record `<lineage>_chr<c>` per chromosome: the assembled genome, its blocks concatenated in physical order, in `genomes/`. **Nucleotide genome runs only**: a family or ordered run has gene families, not coordinates, so there is nothing to lay out, and no `genomes/` is created. The biggest thing this level writes — a whole genome times every node |
| Initial genome | `genome_initial.fasta` | FASTA | yes | the genome the run **started** with, as sequence — the state the stem leads *from*, which is not any node's. In `genomes/` with the rest, being a whole-genome FASTA like they are. Nucleotide runs only |
| Driver views | `.gc()` · `.composition(letters)` | driver | Python | the share of a lineage's sequence that is those letters, pooled over all its families, for use as a driver (Ch9) — GC content, or any amino-acid frequency. A number, so it takes a `Curve` or a `Scalar`; it drives a trait or a further sequence run, never the genome its gene trees came from |
| Summary | `sequences_summary.json` | JSON | yes | what came out — `unit` (`family` or `block`), families with sequences, how many sequences, sites min/max, `mean_pairwise_identity` (the saturation check the command warns on), assembled genomes, and the seed |

On a **nucleotide** genome run every block evolves, spacer as well as gene, so a genome of *b* blocks
writes *b* alignments and *b* phylograms — that is what makes the genomes assemblable. The number in
those filenames is then a **root block index**, not a gene family id, and the files say so: `block6.fasta`,
`phylogram_block6_complete.nwk` and `sequences_ancestral_block6.fasta` in place of `fam6.…`, and
`sequences_founding.fasta`'s records are `block<b>` rather than `fam<f>`. The two numbering schemes are
different, so go from a gene to its block with `genomes.block_of(family)` (Ch6). `zombi2 sequences` reads a nucleotide
handoff too — it recognises one by its `blocks.tsv` — so all of this is reachable from the command line.

## Joint: `simulate_joint` / `zombi2 joint`

A joint run grows two levels at once, so it writes both, each in the format its own command would
give it: the species files, and then the driver's — the trait's or the genomes'. Its own output is one
file: `joint_summary.json`, at the run root beside the run report.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Species tree | `species/species_complete.nwk` · `…_extant.nwk` · `species_events.tsv` · `species_fates.tsv` · `species_summary.json` | Newick, TSV, JSON | yes | the grown tree — complete, so the extinct lineages whose fate the driver decided are kept |
| The trait it grew | `traits/trait_values.tsv` · `trait_events.tsv` · `trait_tree.nwk` · `trait_summary.json` | TSV, Newick, JSON | yes¹ | as the traits level writes them |
| The genomes it grew | `genomes/genome_events.tsv` · `profiles.tsv` · `genomes.tsv` · `initial_genome.tsv` · `gene_trees/` · `genome_summary.json` · `species_complete.nwk` | TSV, Newick, JSON | yes¹ | as the genomes level writes them — including the species tree, which a joint run grew and so keeps here as well as under `species/` |
| Summary | `joint_summary.json` | JSON | yes | both levels' summaries in one file — `seed` · `driver` (`trait` / `genome`) · `species`, whose realised rates are what the driver did · one of `trait` / `genome`, the driver's own. It sits at the **run root** as well as under each level, because neither level was grown first. The same payloads the two levels write alone, so there is nothing new to read |
| Run log | `species/joint.log` | TSV | yes | the resolved parameters, as every command writes |

¹ whichever driver the run used — one per run, never both.

## Traits: `simulate_continuous` / `simulate_discrete`

`zombi2 traits --name NAME` writes this trait's files to `traits/NAME/` instead of `traits/`, so one
run can hold several traits and one can drive another (Ch9). `--name` and `--flat` are refused
together.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Values | `trait_values.tsv` | TSV | yes | value at every node (tips, extinct, internal) — `node` · `kind` · `trait`, where `kind` is the tip's fate (`extant` / `extinct` / `unsampled`) or `ancestor`, so `kind == "extant"` isolates the observed tips |
| Events | `trait_events.tsv` | TSV | yes (CLI, discrete) | the trait's whole history — an `initial` row giving the state at t=0, then every switch: `time` · `kind` · `lineage` · `from` · `to`, where `kind` is `initial` · `on_branch` · `on_speciation`. The one event log whose payload is a **state change** rather than a birth and a death, so it keeps its `lineage` and has no `parents` / `children`. Times are full precision (they drive a conditioned run's Gillespie). **This is also the driver file**: a genome, sequence or trait run drives a rate with `scaled_by("trait_events.tsv", …)`, replaying it against the shared tree. A continuous trait carries only the `initial` row and any `at_speciation` jumps (a diffusion can't be rebuilt from events), and that holds for a multi-optimum (`regimes=`) run and a **correlated** multi-trait one alike. A correlated run **widens** the table instead of repeating a row per trait — `from:<trait>` · `to:<trait>`, one pair apiece, exactly as `trait_values.tsv` widens — because a correlated jump moves every trait at once and is one event |
| Trait tree | `trait_tree.nwk` | Newick | yes (CLI) | tree with every node annotated `[&trait=…]` (opens in FigTree / iTOL). `zombi2 traits` writes it by default; the Python `TraitsResult.write` default is `("values",)` alone |
| Summary | `trait_summary.json` | JSON | yes (CLI) | what came out, not what was asked for — `tips` · `nodes` · `events` (the `on_branch` and `on_speciation` counts), then `states` · `most_common_share` for a discrete trait, or `values` (min/mean/max) · `value_at_root_node` for a continuous one. The root node sits at the end of the stem, so that value is not the one the run started from |
| Conditioning | `conditioned_on` | text | conditioned | the levels this run reads as a driver, one per line, in the trait's own directory — a trait driven by a trait grown first records `traits` (Ch9). Written only when `--rate` or `--switch` was conditioned. Both sides sit under `traits/`, so the record is kept but re-running the driver trait does not invalidate this run |
| Tip names | `names.tsv` | TSV | external tree | `node` · `name`, mapping ZOMBI2's `n<id>` back to the labels you supplied. Written only when the tree came from `--from` with its own tip labels; it is the join from every other output back to your taxa |

## Conditioning and joining: no new files

Neither adds a format. A **conditioned** run writes the driven level's own files and one extra record,
`conditioned_on`, naming what it read (above), so the pairing is kept on disk; a **joint** run writes
**both** levels, each in its own format.

The `zombi2 tools` commands write their own files — the homology matrix and the reconciliation/scoring
outputs — catalogued in Appendix C.
