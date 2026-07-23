# Output files

Every `simulate_*` returns a result; `result.write("out/", outputs=[...])` writes the files, and
omitting `outputs` writes the **default** set. Trees are Newick, tables and logs are TSV, sequences are
FASTA. Tree branch lengths are **time** everywhere except the sequence phylograms, whose lengths are in
**substitutions per site**. The **Default** column says whether a file is written with no arguments
(**yes**), only when you name its token (**no**), or is available in Python but has no file yet
(**Python**).

A species-tree node is written `n<id>` everywhere it appears, and the column holding one is always
called `lineage` (or `donor` / `recipient` where a row names two), so a node reads the same in any
file of a run. A gene copy is always `copy`. In `blocks.tsv` alone, `gene` means something else — the
genic classification of a block, `0` for spacer or the family id for a declared gene.

## Where the files go

A `zombi2` command groups what it writes, one directory per level, and gives the outputs that run to
one file per gene family a directory of their own:

```
out/species/                species_complete.nwk · species_extant.nwk · species_events.tsv
out/genomes/                genome_events.tsv · profiles.tsv · genomes.tsv · genomes.log
out/genomes/gene_trees/     gene_tree_fam<f>_complete.nwk · …_extant.nwk
out/sequences/              clock_species_tree_complete.nwk · sequences.log
out/sequences/alignments/   fam<f>.fasta
out/sequences/phylograms/   phylogram_fam<f>_*.nwk
out/traits/                 trait_values.tsv · trait_tree.nwk · trait_events.tsv · traits.log
```

Filenames keep their prefix inside their directory — `species/species_complete.nwk` — so a file
still names itself once it has been moved or copied somewhere else.

`--flat` on any command writes its files straight into the output directory instead, for a tool that
expects one directory. The same files are written either way; only the directories differ. A run's
files being grouped is a CLI matter, not the library's: `result.write(dir)` writes into whatever
directory you hand it.

Every command that reads a prior level takes the **run directory**, and finds the file itself, in
either layout: `zombi2 genomes out/` and `zombi2 traits out/` pick up that run's complete species
tree, and `zombi2 sequences out/` picks up its handoff files. `--from` reads from somewhere else —
a Newick file, or another run — which is also how you write a run separate from the one you read.

Every `zombi2` **command** also writes a run log (`species.log`, `genomes.log`,
`sequences.log`, `traits.log`): the version, the timestamp, the command line, and every resolved
parameter. Rates are recorded in their **written form** — `birth<TAB>1.0 * OnTime({0: 1, 3: 0.3})` —
so a line pastes straight back into the flag or a `--params` file. It is a CLI artifact, not a
`result.write()` output, so it has no row in the tables below.

## Species trees — `simulate_species_tree`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Complete tree | `species_complete.nwk` | Newick | yes | every lineage, including extinct and unsampled |
| Extant tree | `species_extant.nwk` | Newick | yes | only the sampled survivors |
| Event log | `species_events.tsv` | TSV | yes | every speciation/extinction — `time` · `kind` · `lineage` · `children` |
| Fossils | `species_fossils.tsv` | TSV | yes¹ | sampled fossil lineages — `lineage` · `time` |

¹ written only if fossil sampling recovered any.

## Genomes, unordered — `simulate_genomes_unordered`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | the source of truth — `time` · `kind` · `lineage` · `family` · `copy` · `parent` · `recipient` · `donor`. A transfer writes two rows sharing a `parent`, and **both name `donor`**, the branch the material left, so either row gives the whole edge without pairing them. On the arriving row `lineage` and `recipient` are the same branch, so `donor` is what makes a self-transfer (`donor` == `recipient`) visible at all |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Genomes | `genomes.tsv` | TSV | yes | every node's gene content, **ancestors included** — `lineage` · `family` · `copy`. **One row per gene copy**, so a lineage holding six genes has six rows; two rows sharing a `family` are two copies of it. `copy` is the same identifier the event log uses, so a gene can be followed from the genome it sits in back to the event that made it. `profiles.tsv` is the same information counted, and only for the extant tips |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `family` · `copy`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | each family's true genealogy, in `genomes/gene_trees/`. A family with no surviving copy writes no `_extant` file |
| Family origination | `.gene_trees[f].origination` | float | Python | when the family was founded — where its gene tree's root branch begins |


## Genomes, ordered — `simulate_genomes_ordered`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the run's whole history, in one time-ordered table.** The gene genealogy as unordered, plus **where** each event happened and the ancestry-neutral rearrangements as kinds of their own — `time` · `kind` · `lineage` · `family` · `copy` · `parent` · `recipient` · `donor` · `dest_lineage` · `chromosome` · `position` · `length` · `dest_chromosome` · `dest_position` · `flipped`¹. `position` / `length` are coordinates in that branch's own genome just before the event, as `gene_order` numbers it, and are filled **once per event** — on its first row — because the arc belongs to the event, not to each copy it touched. Filtering on a non-empty `position` therefore gives one row per event that moved genes, which is what a replay walks. A transfer writes a row on each branch and each names the whole edge: `donor` is on both, and the departing row adds `dest_lineage` for where the material went |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Gene order | `gene_order.tsv` | TSV | yes | signed gene order of **every node**, ancestors included — `lineage` · `chromosome` · `position` · `strand` · `family` · `copy` |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `strand` · `family` · `copy`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — `time` · `kind` · `lineage` · `parents` · `children` |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | as unordered — position is orthogonal to genealogy |

¹ a run is named by `start` (its first position, in the chromosome's frame just before the event) and
`length` (how many genes it covered), counted rightwards from `start` and **wrapping past position 0 on
a circular chromosome** — so `start + length` greater than the chromosome's gene count means the run
crossed the origin. `dest_position` is an index into what was left after the run was excised.

## Genomes, nucleotide — `simulate_genomes_nucleotide`

From `zombi2 genomes --resolution nucleotide` or `result.write(dir, outputs=[...])`.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the run's whole history, in one time-ordered table**: the copy-lineage genealogy and the ancestry-neutral rearrangements — `time` · `kind` · `lineage` · `chromosome` · `copy` · `parent` · `recipient` · `source` · `start` · `end` · `position` · `length` · `dest_chromosome` · `dest_position` · `flipped`. One row per **ancestral interval** an event touched, so an event spanning several blocks writes several rows. `source`/`start`/`end` are **ancestral** coordinates — which stretch of which root sequence — while `position`/`length` are **physical** ones along the chromosome, as `blocks.tsv` numbers it: different frames, so different columns |
| Blocks | `blocks.tsv` | TSV | yes | every node's genome as its block mosaic, ancestors included — `lineage` · `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. The rows of one chromosome tile it end to end from 0. The largest file this level writes: blocks are not kept maximal during a run, so it grows with their number × every node |
| Genes | `genes.tsv` | TSV | yes | the declared genes in root coordinates — `family` · `name` · `source` · `start` · `end` · `strand` (the **coding** strand). Header-only when none were declared |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — same format as ordered |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | one tree per declared gene (else per recovered root-block) |

The nucleotide log needs no separate positions file: its events carry ancestral coordinates already.

The events index against the species tree canonicalised so its `n<id>` labels match the `lineage`
column, so a genomes run needs that exact tree to be replayable. A run grown by `zombi2 species`
already has it; a run whose tree came from `--from` gets a copy written to its own
`species/species_complete.nwk`, rather than a second file under another name. Either way `zombi2
sequences` can replay the gene genealogy from the run directory alone. Like `names.tsv` (external
input trees) and the `.log`, that copy is a CLI artifact, not a `result.write()` output.

## Sequences — `simulate_sequences`

The `zombi2 sequences` command replays a prior `zombi2 genomes` run — its own run directory, or `--from` another —
its species tree and its `genome_events.tsv`. Gene outputs are written **one file per gene
family** (`<f>` = family number); a family with no surviving copy writes none. Every node is labelled
`g<copy>`, so a phylogram's tips pair with its alignment and its internal nodes with the ancestral
sequences.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Alignments | `fam<f>.fasta` | FASTA | yes | one row per extant gene copy — nucleotides or amino acids, following the model. The command puts these in `sequences/alignments/`, which is what lets the name be this short |
| Phylograms | `phylogram_fam<f>_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the gene tree each family's sequences were drawn along, in `sequences/phylograms/` |
| Ancestral | `sequences_ancestral_fam<f>.fasta` | FASTA | no | the sequence at every node that is not an extant tip: internal nodes, and the tips where a copy was lost or its species died |
| Founding | `sequences_founding.fasta` | FASTA | no | one record `fam<f>` per family — the sequence it originated with, where its phylogram's root branch begins |
| Clock species tree | `clock_species_tree_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the species tree with its branches in substitutions/site — the molecular clock made visible |
| Genomes | `genome_<lineage>.fasta` | FASTA | yes | one file per **node** of the complete tree — extant, extinct and ancestral alike — one record `<lineage>_chr<c>` per chromosome: the assembled genome, its blocks concatenated in physical order. **Nucleotide genome runs only**: an unordered or ordered run has gene families, not coordinates, so there is nothing to lay out. The biggest thing this level writes — a whole genome times every node |
| Initial genome | `genome_initial.fasta` | FASTA | yes | the genome the run **started** with, as sequence — the state the stem leads *from*, which is not any node's. Nucleotide runs only |

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
| Species tree | `species/species_complete.nwk` · `…_extant.nwk` · `species_events.tsv` | Newick, TSV | yes | the grown tree — complete, so the extinct lineages the coupling decided the fate of are kept |
| The trait it grew | `traits/trait_values.tsv` · `trait_events.tsv` · `trait_tree.nwk` | TSV, Newick | yes¹ | as the traits level writes them |
| The genomes it grew | `genomes/genome_events.tsv` · `profiles.tsv` · `genomes.tsv` · `gene_trees/` | TSV, Newick | yes¹ | as the genomes level writes them |
| Run log | `species/joint.log` | TSV | yes | the resolved parameters, as every command writes |

¹ whichever driver the run used — one per run, never both.

## Traits — `simulate_continuous` / `simulate_discrete`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Values | `trait_values.tsv` | TSV | yes | value at each extant tip — `node` · `trait` |
| Events | `trait_events.tsv` | TSV | yes (discrete) | the trait's whole history — a `root` row giving the state at t=0, then every switch: `time` · `kind` · `lineage` · `from` · `to`, where `kind` is `root` · `on_branch` · `on_speciation`. Times are full precision (they drive a conditioned run's Gillespie). **This is also the conditioning file**: a genome/sequence run drives a rate with `mod.DrivenBy("trait_events.tsv", …)`, replaying it against the shared tree. A continuous trait carries only the `root` row and any `at_speciation` jumps (a diffusion can't be rebuilt from events) |
| Trait tree | `trait_tree.nwk` | Newick | no | tree with every node annotated `[&trait=…]` (opens in FigTree / iTOL) |

## Coupling — no new files

Coupling adds no formats. A **conditioned** run writes the target level's files plus the **driver file**
it read (above), keeping the pairing on disk; a **joint** run writes **both** levels, each in its own
format.

## Tools

The scoring and reconciliation commands (tree-distance, reconciliation-accuracy, the undated simulator,
the ALE likelihood) are quarantined during the clean-core rebuild; their files are documented here as
each returns.
