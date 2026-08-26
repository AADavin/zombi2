# Output files

One table per level, listing the files it writes: the file, and what it holds. Every file listed is
written with no arguments unless its line says otherwise; a few are written only when you name their
**token**, one of the short names `outputs=` and `--write` take, `events`, `alignments` or
`ancestral`. Those lines say so, naming the token where it is not simply the file's own name. What a
run gives back in **Python** follows each table: the objects and accessors, which are not files.
Trees are Newick, tables and logs TSV, sequences FASTA; branch lengths are time, except in the
sequence phylograms, which are in substitutions per site.

Where a file needs more than a phrase, it gets a paragraph under its table, in the table's own
order.

**Columns may be added; read them by name.** A later version may add a column to any of these
tables, because a new kind of event needs somewhere to say what it did. What will not change is an
existing column's name or its meaning, and no column is removed within a major version. So a script
that reads a column by its header keeps working across upgrades, and one that reads by position, or
assumes a column count, does not. New columns go at the **end** of the row. The same file already
differs by resolution, since `genome_events.tsv` has five columns at the family resolution and ten
at the ordered one, so matching by name is what these tables have always required.

## Where the files go

```
out/run.zombi2              the run report: one page for the whole run
out/joint_summary.json      the joint run's own file                          (joint runs)
out/species/                species_{complete,extant}.nwk · species_events.tsv · species_fates.tsv
out/species/                species_fossils.tsv                               (fossil sampling)
out/genomes/                genome_events.tsv · initial_genome.tsv
out/genomes/                genomes.tsv                                       (family)
out/genomes/                profiles.tsv                                      (family, ordered)
out/genomes/                rearrangement_events.tsv · chromosome_events.tsv  (ordered, nucleotide)
out/genomes/                gene_order.tsv                                    (ordered)
out/genomes/                block_events.tsv · blocks.tsv · genes.tsv         (nucleotide)
out/genomes/                initial_sequence.fasta                            (nucleotide + fasta)
out/genomes/                names.tsv                                    (--from with own labels)
out/genomes/markers.tsv     the marker table, one file for the run: Appendix D  (via zombi2 tools format)
out/genomes/gene_trees/     gene_tree_fam<f>_complete.nwk · …_extant.nwk
out/genomes/gff/            genome_<lineage>.gff · genome_initial.gff         (nucleotide)
out/genomes/bed/            genome_<lineage>.bed · genome_initial.bed         (nucleotide)
out/genomes/homology/       homology_fam<f>.tsv                          (via zombi2 tools format)
out/genomes/recphylo/       recphylo_fam<f>.xml                          (via zombi2 tools format)
out/sequences/              clock_species_tree_complete.nwk · …_extant.nwk
out/sequences/              sequence_events.tsv                          (record=True)
out/sequences/              sequences_founding.fasta                     (the founding token)
out/sequences/alignments/   fam<f>.fasta
out/sequences/ancestral/    sequences_ancestral_fam<f>.fasta             (the ancestral token)
out/sequences/phylograms/   phylogram_fam<f>_*.nwk
out/sequences/genomes/      genome_<lineage>.fasta · genome_initial.fasta  (nucleotide runs)
out/traits/                 trait_values.tsv · trait_tree.nwk · trait_events.tsv
out/traits/<name>/          the same files, when --name was given
out/<level>/conditioned_on  the driver record                            (conditioned runs)
```

Two rules hold everywhere before any table. A lineage is written `n<id>`, and `e<id>` once it goes
extinct, in every file including the trees, while an unsampled tip keeps its `n<id>`; a gene copy
is `n<species>_g<copy>` on the same convention. And a run whose tree came from `--from` keeps a copy
of that tree, with `species_fates.tsv` beside it, in the run's own `species/` directory, so the
directory stands alone for whatever reads it next.

**`run.zombi2`**: the run report, plain text for reading rather than parsing. One section for each
level in the directory, holding the seed, the parameters as resolved, what came out in numbers, and
each file with a one-line gloss; then a closing **TO REPRODUCE** with the exact commands, seeds filled
in, that rebuild the whole chain. Each command rewrites it, so it always describes the directory as
it stands.

Each level directory also holds that command's log (`species.log`, `genomes.log`, `sequences.log`,
`traits.log`) with the version, the command line and every resolved parameter, rates in their written
form. Beside it, a summary of what came out: `species_summary.json`, `genome_summary.json`,
`sequences_summary.json`, `trait_summary.json`, written at every level, and at all three genome
resolutions, except by `zombi2 genomes --stream`, which writes each family straight to disk and never
holds the whole run to summarise. A streamed run also keeps its own `species_complete.nwk` beside its
tables. It still reads the tree from `species/`, as any genomes run does; the copy is so that the
`genomes/` directory stands alone as the handoff a later level or another tool reads.

`--flat` writes every file straight into the run directory instead, except two that a flat run does
not write at all: `run.zombi2` and `conditioned_on`. Both are built from a grouped run's level
directories, which a flat run does not have.

The grouping above is the CLI's. **A Python run has no level directories**: `result.write("out/")`
fills the directory you name with that result's files, the per-family and per-node directories
included, so it is not `--flat`, which flattens those too. The same run written from Python and
from the CLI therefore leaves two different layouts, and `run.zombi2` and `conditioned_on` are
written by the CLI alone. To group a Python run, name the level directory yourself:
`result.write("out/species/")`.

Every directory a run fills with one file per family or per node (`gene_trees/`, `gff/` and `bed/` at
the genome level, `alignments/`, `ancestral/`, `phylograms/` and `genomes/` at the sequence level) is
**emptied before that run fills it**, so its contents describe one run and nothing else. Two
exceptions: a `--flat` run empties nothing, and the `homology/` and `recphylo/` directories that
`zombi2 tools format` writes are left as they are.

## Species trees: `simulate_species_tree`

| File | What it holds |
|---|---|
| `species_complete.nwk` | every lineage, extinct and unsampled included |
| `species_extant.nwk` | only the sampled survivors |
| `species_events.tsv` | every speciation and extinction, one row each |
| `species_fates.tsv` | each tip's resolved fate: `lineage` · `fate` (`extant` / `extinct` / `unsampled`) |
| `species_fossils.tsv` | sampled fossil lineages: `lineage` · `time`. Only if fossil sampling recovered any |
| `species_summary.json` | counts by fate, tree height, stem, total branch length, the realised birth and death rates, the seed |

From Python: `.complete_tree` · `.extant_tree` (trees), `.fossils` (the sampled fossils and their
times, present when the run asked for them), `.lineage_rates(kind)` (a dict).

**`species_events.tsv`**: `time` · `kind` · `parents` · `children`. A `speciation` row is the
lineage that ended and its two children, `;`-packed (`n0` → `n1;n2`); an `extinction` row is the
dying lineage as the parent with no children. A lineage that died is written `e<id>`.

**`.complete_tree` and `.extant_tree`**. The whole tree that grew, extinct lineages and all, is
what the next level runs along, which is what lets a gene be transferred out of a lineage that later
dies; the extant tree is the survivors', dated and bifurcating, the one an analysis would be handed.
Each holds every node, internal ones included, in `.nodes`, and answers `.leaves()`,
`.extant_leaves()`, `.extinct_leaves()` and `.unsampled_leaves()`.

**`.lineage_rates(kind)`**: `{lineage: rate}`, the birth (or `"death"`) rate each lineage itself ran
under, taken at the lineage's birth, which is when its factors are drawn. Every lineage of the
complete tree has one, extinct ones included. A rate that varies among lineages gives each a
different number; one that does not gives the same number for all of them. Under a rate that also
depends on time or diversity, this is its value at the start of the branch.

## Genomes, family: `simulate_genomes_family`

| File | What it holds |
|---|---|
| `genome_events.tsv` | the run's true history: `time` · `kind` · `family` · `parents` · `children`, [one row per event](#one-row-per-event) |
| `profiles.tsv` | family × extant-species copy counts |
| `genomes.tsv` | every node's gene content, ancestors included: `lineage` · `family` · `copy`, one row per gene copy |
| `initial_genome.tsv` | the genome the run **started** with, at the start of the root branch: `family` · `copy` |
| `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | each family's true genealogy, in `genomes/gene_trees/` |
| `genome_summary.json` | events by kind, families born, surviving and died out, genes per genome, `empty_genomes`, whether the family-size cap bound, the seed |
| `species_complete.nwk` | the tree the run evolved along, without which the directory cannot be read |
| `species/species_fates.tsv` | each tip's fate, in the format the species level writes. Only when the tree came from `--from` |
| `names.tsv` | `node` · `name`, mapping ZOMBI2's `n<id>` back to the labels you supplied. Only when the tree came from `--from` with its own tip labels |
| `conditioned_on` | the levels this run depends on as a driver, one per line. Only when something was conditioned |

From Python: `.genomes` · `.node_genomes` (the genomes), `.family_counts(node)` ·
`.has_family(node, name)` (a node's genome as `family → copies`, and whether a named family has a
copy there), `.gene_trees[f].origination` (when a family was founded, where its gene tree's root
branch begins), the driver views `.presence(name)` · `.completion(name)`, and `.gene_trees` ·
`.profiles` · `.initial_genome` · `.events` · `.seed`, the same objects the files hold before
`.write()` puts them on disk.

**`genomes.tsv`**: `copy` is the identifier the event log uses, so a gene can be traced back to the
event that made it. `profiles.tsv` is the same information counted, for the extant tips only.

**`initial_genome.tsv`** has no `lineage` column because it is no node's genome: every node
sits at the *end* of its branch.

**`gene_tree_fam<f>_*.nwk`**: leaves are `n<species>_g<copy>`; internal nodes are labelled
`<event>_n<species>` (`duplication_n45`, `transfer_n45`), naming the event that ended that gene and
the branch it was on. A family with no surviving copy writes no `_extant` file.

**`species_complete.nwk`**. Every other file here is indexed by its node labels, so the directory is
not readable, by anyone or by `genomes.read_run()`, without it. `result.write()` writes it by
default, so a directory written from Python stands alone; `zombi2 genomes` leaves it out, except
under `--stream`, which is family-only, because a run already keeps one copy at
`species/species_complete.nwk`, shared by every level. Its token in `outputs=` and `--write` is
`species_tree`, not the file's name.

**`names.tsv`**: the join from every other output back to your taxa.

**`conditioned_on`** is written whether the driver was used for a rate or for `transfer_to`.
Re-running a driver level then refuses rather than leaving this run stale, unless you pass `--force`.

**`.genomes` and `.node_genomes`**. `.genomes` is the observed dataset, the genome at each *extant*
tip, keyed by the tip name the tree writes (`n5`), so it joins to the tree and to a trait grown on the
same tree. `.node_genomes` is the run's own record, every node, extant and extinct and internal
alike, keyed by node id, for joining against `.complete_tree.nodes` or the event log. The pair is the
same at every resolution, and the distinction is the model's: one is what a dataset contains, the
other is what happened.

**`.presence(name)` and `.completion(name)`**: a named family's presence (`present` / `absent`) and
a declared module's completion (a fraction) along every lineage, for use as a driver (Ch8). Read off
the families' gene trees, so they change *inside* a branch; the **ordered** resolution gives the same
two.

### One row per event

Every event log opens with `time` and `kind`, then the payload columns for that file, the same on
every row. Here they are `family`, `parents` and `children`: the gene copies the event ended and the
ones it began, `;`-packed.

A gene copy is written `n<species>_g<copy>`, so copy 30 on branch `n2` is `n2_g30`, and `e2_g30` if
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
| `transfer_replacing` | the donor's copy, then the copy it overwrote | the continuations of those two lineages |

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
carries two parents, the donor's copy and then the copy it overwrote on the recipient branch, and
writes **no separate `loss` row** for the copy it displaced. That death is what the kind means.

The file is time-ordered, one row per event, so counting rows by kind counts events, except at the
**ordered** resolution, where a segmental event writes one row per gene it covers; dedupe on `time`
there before counting (the ordered section's Segment note). An absent value is an empty cell. `genome_summary.json` counts the
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

| File | What it holds |
|---|---|
| `genome_events.tsv` | the gene genealogy, with the place each event happened |
| `rearrangement_events.tsv` | every inversion, transposition and translocation |
| `profiles.tsv` | family × extant-species copy counts |
| `gene_order.tsv` | signed gene order of every node, ancestors included: `lineage` · `chromosome` · `topology` · `position` · `strand` · `family` · `copy` |
| `chromosome_events.tsv` | chromosome-network edges: `time` · `kind` · `parents` · `children` |
| `initial_genome.tsv` | the genome the run **started** with: `chromosome` · `topology` · `position` · `strand` · `family` · `copy` |
| `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | as at the family resolution, since position is orthogonal to genealogy |
| `genome_summary.json` | events as biology rather than rows, families born/surviving/died out, genes and chromosomes per genome, rearrangements and chromosome events by kind |
| `species_complete.nwk` | as at the family resolution |
| `names.tsv` | as at the family resolution |
| `conditioned_on` | as at the family resolution, and written when a rate or `transfer_to` was conditioned |

From Python: `.genomes` · `.node_genomes` (as at the family resolution, but each genome is a tuple of
**`Chromosome`** objects, each an `id`, a `topology`, and an ordered list of **`Gene`** objects with
`id`, `family` and `strand`), `.gene_order(node)` (one node's layout gene by gene, as
`(chromosome, position, strand, family, gene id)`), and `.rearrangements` · `.chromosome_events` ·
`.gene_trees` · `.profiles` in memory.

**`genome_events.tsv`**. The family resolution's five columns plus five more: `time` · `kind` ·
`family` · `parents` · `children` · `chromosome` · `start` · `length` · `dest_chromosome` ·
`dest_position`. Kinds and copy names are [the shared ones](#one-row-per-event). Coordinates are in
that branch's own genome just before the event, as `gene_order.tsv` numbers it. `dest_position` and
`dest_chromosome` are where the material landed, so one row carries both ends of a transfer. A
speciation copies a genome whole and leaves all five coordinate cells empty.

**`rearrangement_events.tsv`**: `time` · `kind` · `lineage` · `chromosome` · `start` · `length` ·
`dest_chromosome` · `dest_position` · `flipped`. These begin and end no gene lineage, so they have no
`parents` and `children` and get a file of their own; a segment has no name either, which is why this
is the one event log that still puts its branch in a column. `dest_chromosome` is set only by a
translocation; `flipped` is `1` or `0` on the two events that move a segment, which may land
inverted, and empty on an inversion, which flips by definition. The
`events` token writes this file and the log above together.

**A segment**, in both of those files, is named by `start` (its first position, in the chromosome's
frame just before the event) and `length` (how many genes it covered), counted rightwards from
`start` and **wrapping past position 0 on a circular chromosome**, so `start + length` greater than
the chromosome's gene count means the segment crossed the origin. `dest_position` is an index into
what was left after it was excised. A segmental event acts on a run of genes of several families at
once and `genome_events.tsv` has one `family` column, so it writes one row per gene lineage, each
repeating the same arc.

**`gene_order.tsv`**. `topology` is `circular` or `linear`, written beside every gene: it decides
where a segmental event stops and which chromosomes may fuse. A chromosome with no genes has no rows
here, and so no topology.

**`chromosome_events.tsv`**. Chromosomes are named `n<species>_c<id>`, on the same pattern as gene
copies. Kinds are `initial` (a replicon the run starts with, at time 0), `speciation`, `fission`,
`fusion`, `origination` and `loss`; a fusion is the one row with two parents.

**`initial_genome.tsv`** has no `lineage` column, for the reason given at the family resolution.


## Genomes, nucleotide: `simulate_genomes_nucleotide`

From `zombi2 genomes --resolution nucleotide` or `result.write(dir, outputs=[...])`.

| File | What it holds |
|---|---|
| `genome_events.tsv` | the genealogy, in the format every resolution writes |
| `block_events.tsv` | this resolution's own record: one row per ancestral interval an event touched |
| `rearrangement_events.tsv` | inversions, transpositions and translocations, in **physical** coordinates |
| `chromosome_events.tsv` | chromosome-network edges, same format and same kinds as ordered |
| `blocks.tsv` | every node's blocks, ancestors included: `lineage` · `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene` |
| `genes.tsv` | the declared genes in initial coordinates: `family` · `name` · `source` · `start` · `end` · `strand` (the **coding** strand). Header-only when none were declared |
| `initial_genome.tsv` | the genome the run **started** with: `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene` |
| `initial_sequence.fasta` | the initial DNA the run was given (`--fasta`), one `>source<n>` record per replicon. Only when a FASTA was supplied |
| `genome_summary.json` | what came out, counted three separate ways |
| `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | one tree per declared gene (else per recovered root-block), in `gene_trees/` |
| `gff/genome_<lineage>.gff` | each node's **genes** in its own coordinates, plus `genome_initial.gff` |
| `bed/genome_<lineage>.bed` | each node's **blocks**, spacer included, each named by the ancestral interval it descends from |
| `species_complete.nwk` | as at the family resolution |
| `names.tsv` | as at the family resolution |
| `conditioned_on` | as at the family resolution, written only when a rate or `transfer_to` was conditioned |

From Python: `.genomes` · `.node_genomes` (each a `NucleotideGenome`: a list of `Chromosome`s, each a
list of `Block`s), `.root_blocks` · `.block_trees` (the recovered root partition and a tree for every
interval), `.assembly(node)` · `.initial_assembly()`, `.mosaic(node)` · `.trace_back(node)` ·
`.describe(node)`,
`.deletions` (the indel log), `.gene_spans` · `.gene_names` · `.gene_strands` · `.block_of(family)`,
and the driver views `.presence(name)` · `.completion(name)`, which here are based on **declared
genes**: the name is the GFF's `ID` / `Name`, and what takes a gene away is an arc of DNA rather
than a whole copy.

**`genome_events.tsv`**: the same `time` · `kind` · `family` · `parents` · `children` as at the
family resolution, [one row per event](#one-row-per-event), so one reader serves all three
resolutions and `zombi2 tools` works here unchanged. `family` is the declared gene (else the
recovered root-block), and the copy in each token is a **gene** id, the one the gene trees,
alignments and homology tables use. It is derived onto the root-block partition, where a copy either
covers a block in full or does not touch it, which is what makes a duplication a bifurcation here. A
nucleotide transfer is always additive, so `transfer_replacing` cannot appear.

**`block_events.tsv`**: `time` · `kind` · `parents` · `children` · `chromosome` · `source` ·
`start` · `end` · `cuts`; copy lineages and chromosomes are named as everywhere else. There is one
row per ancestral interval an event touched, not one row per event, so an event covering several
blocks writes several rows sharing a `time` and `kind`, and a duplication starts a child without
ending its parent (the parent still covers the rest of its extent). That is why the genealogy above
is a separate table: this shape suits sequence and not a gene tree. Kinds are `initial` ·
`origination` · `insertion` · `duplication` · `deletion` · `loss` · `transfer` · `speciation`, where
`initial` is the starting genome laid down at time 0, one row per replicon, and `insertion` /
`deletion` are the **indels**: `insertion` is novel spacer on a fresh source, a root like an
origination but carrying no gene family, and `deletion` removes material without ending a copy
lineage. `cuts` gives the **ancestral** breakpoints that event used as `source:position` pairs,
semicolon separated: it is what tells an indel breakpoint from an ordinary one, and without it a run
with indels could not be read back. A speciation splits one copy lineage into one per daughter
without touching sequence, so it writes a single row with both daughters and every coordinate cell
empty; a transfer's blocks keep their source coordinates and leave `chromosome` empty.
`source`/`start`/`end` are **ancestral** coordinates; the physical ones live in the rearrangement
log. This is what `read_nucleotide_genomes` replays.

**`rearrangement_events.tsv`**: the same columns the ordered resolution writes with base pairs in
place of genes, along the chromosome as `blocks.tsv` numbers it, plus `cuts`, the **ancestral**
breakpoints the rearrangement used, which only this resolution can give (the ordered one leaves it
blank). No token of its own either: `events` writes all three tables here. A nucleotide translocation
records no `dest_position`: its blocks keep their source coordinates and the engine places the arc.

**`blocks.tsv`**: `position` is the block's physical offset, and the rows of one chromosome tile it
end to end from 0. `source` names the coordinate frame an interval lives in, either an initial replicon
or the fresh source an insertion or origination laid down, and `start`/`end` are the interval on
that source. `copy` is the copy-lineage token the event log uses for the block's carrier, and `gene`
numbers the declared gene a block is, 1 upward in declaration order, 0 for spacer. The largest file
this level writes: blocks are not kept maximal during a run, so it grows with their number × every
node.

**`genome_summary.json`**: `events` counted as *biology* rather than rows (see the ordered section),
`block_events` counting this resolution's own log the same way, `genes` (every gene the run holds,
those declared at the start and the de-novo ones `origination` added), and base pairs and chromosomes
per genome. No phyletic profiles here: the unit is the base pair. **The two counters do not bound each
other.** `events` counts gene-tree branchings and `block_events` counts events on the DNA: an arc
covering three genes is three gene events, one covering none is zero. A gene is never split, so an
event changes gene content only when its arc covers a gene end to end, which is why a run with small
extents against big genes reports block events by the dozen and gene events at zero, and warns that it
will. `deletions` and `base_pairs_deleted` count the indel log, which is neither of the other two: an
indel removes material without ending a copy lineage, so it appears in no event counter.
`block_events` also gains a kind `events` has no counterpart for: `insertion` brings sequence rather
than a gene family, so the gene-level counters cannot see it.

**`initial_sequence.fasta`** is what lets a separate `zombi2 sequences` run found its blocks from the
real sequence.

**`gff/` and `bed/`**: one file per node, plus `genome_initial.gff` / `.bed` for the genome the run
**started** with. The GFF is the annotation to read beside the sequence level's
`genome_<lineage>.fasta`, which names its sequences to match; the BED is the ancestry as a browser
track.

**`.describe(node)`** gives one node's genome written out block by block, for reading by eye:
`mosaic` as text, one line per block, each labelled with the gene it is or `intergene`.
Genes declared with a name are named; the rest are numbered by family.

**`.assembly(node)` and `.initial_assembly()`** say how a node's genome is built from the recovered
root blocks: `(block, gene, strand, lo, hi)` in physical order, `[lo, hi)` being the sub-range of that
block the node carries, the whole block unless an indel took a stretch out of it or opened a gap
inside it. `.mosaic(node)` and `.trace_back(node)` are the same genome at two grains, per
block and per nucleotide.

**`.root_blocks` and `.block_trees`**: the maximal never-cut intervals that some node still carries,
and a recovered tree for **every** one of them, spacer as well as gene. Cut at every node's
breakpoints rather than only the survivors', which is what lets any node's genome be rebuilt.

**`.deletions`** holds one `Deletion` per indel deletion: `time` · `lineage` · `chromosome` · `deleted`,
the last being `(copy, source, start, end)` rows in **ancestral** coordinates, the same payload a loss
carries. Kept out of `events` because a deletion ends no copy lineage, so the gene-tree recovery must
not read it. Only deletions are here: an **insertion** brings material that descends from nothing, so
it must begin a copy lineage and is recorded in `events` as a root of its own kind. Written into
`block_events.tsv` under kind `deletion`, which is also how it is read back.

**`.gene_spans` and `.block_of(family)`**: `{family: (source, start, end)}` in initial coordinates, a
named gene's family id and coding strand from a GFF, and the block index a family occupies. That last
one is the join between the two numbering schemes here: `.gene_spans` and `.gene_trees` are keyed by
family, `.root_blocks` and `.block_trees` by block index, both plain integers, so mixing them up is
otherwise silent.

Three event files, because a nucleotide run records three different things: the genealogy in the one
format every resolution writes, the interval record only this resolution has, and what moved without
beginning or ending anything. Ancestral coordinates are in `block_events.tsv` and physical ones in
`rearrangement_events.tsv`: two frames, so two files, and no row carrying columns from both.

The events index against the species tree canonicalised so its `n<id>` labels match the branch named
inside each copy token, so a genomes run needs that exact tree to be replayable. A run grown by `zombi2 species`
already has it; a run whose tree came from `--from` gets a copy written to its own
`species/species_complete.nwk`, rather than a second file under another name, with
`species_fates.tsv` beside it so a later level on this run reads each tip's fate from the record
instead of guessing it from tip depth. Either way `zombi2 sequences` can replay the gene genealogy
from the run directory alone. Like `names.tsv` (external input trees) and the `.log`, those copies
are CLI artifacts, not `result.write()` outputs.

## Sequences: `simulate_sequences`

The `zombi2 sequences` command replays a prior `zombi2 genomes` run, its own run directory or
`--from` another, reading its species tree and its `genome_events.tsv`. Gene outputs are written
**one file per gene family** (`<f>` = family number); a family with no surviving copy writes none.
Every node is labelled `n<species>_g<copy>`, so a phylogram's tips pair with its alignment and its
internal nodes with the ancestral sequences.

| File | What it holds |
|---|---|
| `fam<f>.fasta` | one row per extant gene copy: nucleotides or amino acids, following the model, in `alignments/` |
| `phylogram_fam<f>_complete.nwk` · `…_extant.nwk` | the gene tree each family's sequences were drawn along, in substitutions per site, in `phylograms/` |
| `clock_species_tree_complete.nwk` · `…_extant.nwk` | the species tree with its branches in substitutions per site, the molecular clock made visible |
| `sequences_ancestral_fam<f>.fasta` | the nodes the alignment leaves out, one file per family in `ancestral/`. Only when you name its token, `ancestral` |
| `sequences_founding.fasta` | one record `fam<f>` per family, the sequence it originated with. Only when you name its token, `founding` |
| `genome_<lineage>.fasta` | one file per node of the complete tree, one record `<lineage>_chr<c>` per chromosome, in `genomes/`. Nucleotide genome runs only |
| `genome_initial.fasta` | the genome the run **started** with, as sequence. Nucleotide runs only |
| `sequences_summary.json` | `unit` (`family` or `block`), families with sequences, how many, sites min/max, `mean_pairwise_identity`, assembled genomes, the seed |
| `conditioned_on` | the levels this run depended on as a driver. Only when the substitution rate was conditioned |

From Python: `.alignments` · `.ancestral` (the sequences), `.genomes` · `.node_genomes` ·
`.initial_genome` (the assembled genomes, present only when the run came from a **nucleotide**
genome), `.founding` · `.phylograms` · `.species_phylogram`, and the driver views `.gc()` ·
`.composition(letters)`.

**`phylogram_fam<f>_*.nwk`**: with rate variation across sites the branch length is the **mean** over sites, which is
what the rate classes are normalised to. Under a per-clade model set (`Models`) the lengths still mean
substitutions per site, since every model is normalised to one per unit length, but that normalisation
holds at stationarity, so on a branch whose composition is still relaxing toward its clade's
frequencies the realised count falls a little short of the length written here.

**`clock_species_tree_*.nwk`**: the mean over sites when rates vary across them, as for the phylograms. A driven
substitution rate shows here too: a branch is the rate times the driver integrated along it, so this
is where you read what the trait did. The token that writes the pair, in `outputs=` and `--write`,
is `species_phylogram`, not the file's name.

**`sequences_ancestral_fam<f>.fasta`**: internal nodes, and the tips where a copy was lost or its
species died.

**`genome_<lineage>.fasta`**: the assembled genome, its blocks concatenated in physical order, for
every node, extant, extinct and ancestral alike. A family or ordered run has gene families, not
coordinates, so there is nothing to lay out and no `genomes/` is created. It is the biggest thing this level
writes, a whole genome times every node. `genome_initial.fasta` is the state the stem leads *from*,
which is no node's.

**`conditioned_on`**: one level per line (`traits`, say). It records the dependency so re-running the
trait refuses to leave this run silently stale, or clears it under `--force`. A run with no driven
rate writes no such file.

**`sequences_summary.json`**. `mean_pairwise_identity` is the saturation check the command warns on:
every within-family pair, so it is the same number whether the run was held in memory or streamed.

**`.alignments` and `.ancestral`**. `.alignments` is the observable data, for each family the
sequence at every *extant* gene copy, the alignment a phylogenetic method would be handed;
`.ancestral` is every node the alignment leaves out. The run wrote a sequence at each node as it went,
so these are the exact ancestors, not estimates, and together they account for every node of the tree
exactly once.

**`.gc()` and `.composition(letters)`**: the share of a lineage's sequence that is those letters,
pooled over all its families, for use as a driver (Ch8): GC content, or any amino-acid frequency. A
number, so it takes a `Curve` or a `Scalar`; it drives a trait or a further sequence run, never the
genome its gene trees came from.

On a **nucleotide** genome run every block evolves, spacer as well as gene, so a genome of *b* blocks
writes *b* alignments and *b* phylograms, and that is what makes the genomes assemblable. The number in
those filenames is then a **root block index**, not a gene family id, and the files say so: `block6.fasta`,
`phylogram_block6_complete.nwk` and `sequences_ancestral_block6.fasta` in place of `fam6.…`, and
`sequences_founding.fasta`'s records are `block<b>` rather than `fam<f>`. The two numbering schemes are
different, so go from a gene to its block with `genomes.block_of(family)` (Ch5). `zombi2 sequences` reads a nucleotide
handoff too, recognising one by its `blocks.tsv`, so all of this is reachable from the command line.

## Joint: `joint.simulate` / `zombi2 joint`

A joint run simulates two levels at once, so it writes both, each in the format its own command would
give it: the species files, and then the driver's own, whether the trait's or the genomes'. Its own
output is one file: `joint_summary.json`, at the run root beside the run report.

| File | What it holds |
|---|---|
| `species/species_complete.nwk` · `…_extant.nwk` · `species_events.tsv` · `species_fates.tsv` · `species_summary.json` | the grown tree, complete, so the extinct lineages whose fate the driver decided are kept |
| `traits/trait_values.tsv` · `trait_events.tsv` · `trait_tree.nwk` · `trait_summary.json` | as the traits level writes them, when the trait was the driver |
| `genomes/genome_events.tsv` · `profiles.tsv` · `genomes.tsv` · `initial_genome.tsv` · `gene_trees/` · `genome_summary.json` · `species_complete.nwk` | as the genomes level writes them, when the genome was the driver |
| `joint_summary.json` | both levels' summaries in one file, at the **run root**, the one file that is the joint run's own |
| `species/joint.log` | the resolved parameters, as every command writes |

One driver per run, never both, so a run writes the trait files or the genome files, not each.

**`joint_summary.json`** holds four keys: `seed`; `driver`, saying which level it was (`trait` /
`genome`); `species`, whose realised rates are what the driver did; and a fourth named `trait` or
`genome`, the driver's own summary. It sits at the run root because neither level was grown first. The same payloads the two levels write alone, so there is nothing new
to read.

**The genome files** include the species tree, which a joint run grew and so keeps here as well as
under `species/`.


## Traits: `simulate_continuous` / `simulate_discrete`

`zombi2 traits --name NAME` writes this trait's files to `traits/NAME/` instead of `traits/`, so one
run can hold several traits and one can drive another (Ch8). `--name` and `--flat` are refused
together.

| File | What it holds |
|---|---|
| `trait_values.tsv` | the value at every node: `node` · `kind` · one column per trait, headed by the trait's name (`trait` when unnamed, side by side in a correlated run); `kind` is the tip's fate (`extant` / `extinct` / `unsampled`) or `ancestor` |
| `trait_events.tsv` | the trait's whole history: an `initial` row giving the state at t=0, then every switch |
| `trait_tree.nwk` | the tree with every node annotated `[&trait=…]`, for FigTree or iTOL |
| `trait_summary.json` | `tips` · `nodes` · `events`, then `states` · `most_common_share` for a discrete trait, or `values` (min/mean/max) · `value_at_root_node` for a continuous one |
| `names.tsv` | as at the genome level. Only when the tree came from `--from` with its own tip labels |
| `conditioned_on` | the levels this run depends on as a driver, in the trait's own directory. Only when `--rate` or `--switch` was conditioned |

`zombi2 traits` writes the values, the events, the tree and the summary; the Python
`TraitsResult.write` default matches the kind: a continuous run writes the values, the tree and the
summary, having no switches to log, and a discrete run writes everything, events included.

From Python: `.values` · `.node_values` (the values), `.values_by_id` (`.values` keyed by the bare
integer node id, `5` for tip `n5`, for joining against `.node_values`), and `.history`.

**`trait_events.tsv`**: `time` · `kind` · `lineage` · `from` · `to`, where `kind` is `initial` ·
`on_branch` · `on_speciation`. The one event log whose payload is a **state change** rather than a
birth and a death, so it keeps its `lineage` and has no `parents` / `children`. Times are full
precision, since they drive a conditioned run's Gillespie. **This is also the driver file**: a genome,
sequence or trait run drives a rate with `scaled_by("trait_events.tsv", …)`, replaying it against the
shared tree. A continuous trait carries only the `initial` row and any `at_speciation` jumps (a
diffusion cannot be rebuilt from events), and that holds for a multi-optimum (`regimes=`) run and a
**correlated** multi-trait one alike. A correlated run **widens** the table instead of repeating a row
per trait, giving `from:<trait>` · `to:<trait>`, one pair apiece, exactly as `trait_values.tsv`
widens, because a correlated jump moves every trait at once and is one event.

**`trait_summary.json`**: what came out, not what was asked for. The root node sits at the end of the
stem, so `value_at_root_node` is not the value the run started from.

**`conditioned_on`**. A trait driven by a trait grown first records `traits` (Ch8). Both sides sit
under `traits/`, so the record is kept but re-running the driver trait does not invalidate this run.

**`.values` and `.node_values`**. `.values` is the observable vector, the trait at each *extant* tip,
keyed by tip name (`n5`), the same names the Newick and `trait_values.tsv` use, so the dataset joins
the tree it came from. `.node_values` is every node, extant, extinct and internal alike: the true
ancestors at each split, from the same process that produced the tips. Discrete traits store the state
labels you gave, not integer indices.

**`.history`** is, for a **discrete** trait, the per-branch stochastic character map derived from the
event log: the ordered `(state, duration)` segments each branch passed through. `None` for a
continuous trait, which has no map, and for a threshold trait, whose liability crossings are un-timed.

## Conditioning and joining: no new files

Neither adds a format. A **conditioned** run writes the driven level's own files and one extra record,
`conditioned_on`, naming what it depended on (above), so the pairing is kept on disk; a **joint** run
writes **both** levels, each in its own format.

The `zombi2 tools` commands write their own files, the homology matrix and the reconciliation/scoring
outputs, catalogued in Appendix D.
