# Output files

Every `simulate_*` returns a result; `result.write("out/", outputs=[...])` writes the files, and
omitting `outputs` writes the **default** set. Trees are Newick, tables and logs are TSV, sequences are
FASTA. Tree branch lengths are **time** everywhere except the sequence phylograms, whose lengths are in
**substitutions per site**. The **Default** column says whether a file is written with no arguments
(**yes**), only when you name its token (**no**), or is available in Python but has no file yet
(**Python**).

A species-tree node is written `n<id>` everywhere it appears — or **`e<id>` when that lineage went
extinct**. The number is the identity and the letter an annotation: `5` names the lineage wherever it
appears, and a join can always strip the prefix. The column holding a node is always called `lineage`
(or `donor` / `recipient` where a row names two), so a node reads the same in any file of a run.

The letter marks the one fact about a branch you cannot recover from the tree's shape. It means a
complete tree **states its own extinctions**: the file survives being moved, copied or emailed, where
the sibling `species_fates.tsv` does not, and reading one back needs neither that table nor a guess
from tip depth. Two things deliberately stay `n`. Internal nodes: a speciation is not a fate. And
**unsampled** tips, which are alive — being unsampled is a property of the sampling you asked for,
not of the lineage, so the same branch would be named differently by two runs of the same tree.
`species_fates.tsv` remains the only place that tells an unsampled tip from an extant one.

Nothing that names only **extant** tips changed, which is most of what a downstream tool reads: the
extant tree, `profiles.tsv`, the alignments, the homology tables and the extant gene trees are all
`n<id>` throughout, so a pipeline built on them is unaffected. A gene copy is written `g<id>` the same way, and the column holding one is always
`copy` (or `parent`, the source copy of a duplication or transfer). Where a copy is named somewhere
with no column to say which species it sits in — a Newick leaf, a FASTA record, a homology table
header — it is written **`n<species>_g<copy>`**, both labels joined by a single `_`, so a tip or a
sequence says which genome it came from without a lookup, and splitting on the `_` recovers each
half. In
`blocks.tsv` alone, `gene` means something else — the genic classification of a block, `0` for spacer
or the family id for a declared gene.

## Where the files go

A `zombi2` command groups what it writes, one directory per level, and gives every many-files-per-run
output — one per gene family, or one per node — a directory of its own, so the level directory holds
only its handful of tables and logs:

```
out/species/                species_complete.nwk · species_extant.nwk · species_events.tsv · species_fates.tsv
out/genomes/                genome_events.tsv · profiles.tsv · genomes.tsv · genomes.log
out/genomes/gene_trees/     gene_tree_fam<f>_complete.nwk · …_extant.nwk
out/genomes/homology/       homology_fam<f>.tsv                          (zombi2 tools format)
out/genomes/recphylo/       recphylo_fam<f>.xml                          (zombi2 tools format)
out/sequences/              clock_species_tree_complete.nwk · sequences.log
out/sequences/alignments/   fam<f>.fasta
out/sequences/phylograms/   phylogram_fam<f>_*.nwk
out/sequences/genomes/      genome_<lineage>.fasta                       (nucleotide runs)
out/traits/                 trait_values.tsv · trait_tree.nwk · trait_events.tsv · traits.log
```

Filenames keep their prefix inside their directory — `species/species_complete.nwk` — so a file
still names itself once it has been moved or copied somewhere else.

The two groupings are decided in different places. Which **level** a file belongs to is the command's
business: `result.write(dir)` writes one level into the one directory it is handed, and a pipeline
gets `species/`, `genomes/`, `sequences/` because the command calls it once per level. The
subdirectory *inside* a level is the result's own: an output that is one file per gene family or per
node — the gene trees, the alignments, the phylograms, the assembled genome FASTAs, `gff` and `bed` —
gets a directory of its own wherever it is written, because a hundred families is two hundred files
and the handful of tables beside them would be lost. So `result.write("out/")` from Python produces
`out/gene_trees/` exactly as the command produces `out/genomes/gene_trees/`. An output the run has
none of writes nothing and leaves no empty directory behind.

`--flat` on any command writes every file straight into the output directory instead, for a tool that
expects one directory; in Python it is `result.write(dir, flat=True)`. The same files are written
either way — only the directories differ.

Every command that reads a prior level takes the **run directory**, and finds the file itself, in
either layout: `zombi2 genomes out/` and `zombi2 traits out/` pick up that run's complete species
tree, and `zombi2 sequences out/` picks up its handoff files. `--from` reads from somewhere else —
a Newick file, or another run — which is also how you write a run separate from the one you read.

Every `zombi2` **command** also writes a run log (`species.log`, `genomes.log`,
`sequences.log`, `traits.log`): the version, the timestamp, the command line, and every resolved
parameter. Rates are recorded in their **written form** — `birth<TAB>1.0 * OnTime({0: 1, 3: 0.3})` —
so a line pastes straight back into the flag or a `--params` file. It is a CLI artifact, not a
`result.write()` output, so it has no row in the tables below.

## The run summary

Beside each log sits a **summary**: `species_summary.json`, `genome_summary.json`,
`sequences_summary.json`. The log is what the run was *asked* for; the summary is what *happened*. They
are separate files so that a script reading outcomes does not have to step over parameters, and it is
JSON because the people who wanted it were writing collectors.

Unlike the log, it is a `result.write()` output — `outputs=("summary",)` from Python, `--write summary`
on the command line — so a run from either side has one, and `result.summary()` gives the same payload
as a dict without touching the disk.

**Every event count in it is deduplicated**: one number per event, not per row of the event log (see
the section below on why those differ). It also answers two questions the rest of the output leaves
open:

| Field | The question it answers |
|---|---|
| `events` | how many duplications, transfers, losses and speciations, counted **once each** — and `initial` separated from `origination`, since the starting genome is logged as origination too |
| `event_rows` | what a plain `wc -l` on the event log gives, so the difference is visible rather than a trap |
| `families.born` / `.surviving` | why `gene_trees/` holds more file pairs than the run reported families: the run reports the survivors |
| `family_size_cap` | whether `max_family_size` **bit**, and which families are sitting at it. A family at the ceiling had duplications and arriving transfers discarded, so its realised rates are below the ones you declared. This is the only place that says so |
| `realised_rates` (species) | speciations and extinctions per unit of branch length — the cheapest check there is that a tree came out at the rates you asked for |
| `mean_pairwise_identity` (sequences) | the number the run prints and warns on, in machine-readable form |

```python
from zombi2 import species

run = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
run.summary()["realised_rates"]      # {'birth': …, 'death': …}
run.write("out/species/", outputs=("complete", "summary"))
```

Every tree ZOMBI2 writes is ordinary Newick and opens in **FigTree**, **iTOL**, **Dendroscope** or
any other viewer; the species and gene trees carry internal-node labels, which most viewers will show
as node names. To load one in Python, `zombi2.tree.read_newick(text)` returns `(tree, names)` — the
`Tree` the levels take, and `{node id: your label}` for a tree that came with its own tip names. It is
the same reader `--from` uses, so an external tree behaves identically from Python and from the
command line:

```python
from zombi2.tree import read_newick

tree, names = read_newick(open("mytree.nwk").read())
```

A run that **read** a file records it too, one `input` line per file, holding its SHA-256 and its
path:

```
input	d6db05f110039fac…b52dfc	out/species/species_complete.nwk
input	cb2f514472da3958…3743f	out/species/species_fates.tsv
```

The species tree, a `--tip-fates` table, and the driver file of any conditioned rate are all inputs
in this sense. A path on its own does not pin a run down — two runs from two different trees can
carry the same `--from tree.nwk`, the same seed, and the same parameters — so the content is what
tells them apart, and what says a rerun is reading what it read before. See *Reproducing a run* in
Chapter 2.

## Species trees — `simulate_species_tree`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Complete tree | `species_complete.nwk` | Newick | yes | every lineage, including extinct and unsampled |
| Extant tree | `species_extant.nwk` | Newick | yes | only the sampled survivors |
| Event log | `species_events.tsv` | TSV | yes | every speciation/extinction — `time` · `kind` · `lineage` · `children` |
| Tip fates | `species_fates.tsv` | TSV | yes | each tip's resolved fate — `lineage` · `fate` (`extant` / `extinct` / `unsampled`) |
| Fossils | `species_fossils.tsv` | TSV | yes¹ | sampled fossil lineages — `lineage` · `time` |

¹ written only if fossil sampling recovered any.

## Genomes, family — `simulate_genomes_family`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | the source of truth — `time` · `kind` · `lineage` · `family` · `copy` · `parent` · `recipient` · `donor` · `event`. **One row per gene-tree edge, so duplications, transfers and speciations write two rows each** — see [One row per gene-tree edge](#one-row-per-gene-tree-edge-not-per-event) below before counting anything. A transfer's two rows share a `parent`, and **both name `donor`**, the branch the material left, so either row gives the whole edge without pairing them. On the arriving row `lineage` and `recipient` are the same branch, so `donor` is what makes a self-transfer (`donor` == `recipient`) visible at all |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Genomes | `genomes.tsv` | TSV | yes | every node's gene content, **ancestors included** — `lineage` · `family` · `copy`. **One row per gene copy**, so a lineage holding six genes has six rows; two rows sharing a `family` are two copies of it. `copy` is the same identifier the event log uses, so a gene can be followed from the genome it sits in back to the event that made it. `profiles.tsv` is the same information counted, and only for the extant tips |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `family` · `copy`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Conditioning | `conditioned_on` | text | conditioned | written **only when a rate was conditioned**: the run's levels this run read via `DrivenBy` (one per line, e.g. `traits`). It records the dependency so re-running a driver level (say the trait) refuses to leave this run silently stale, or clears it under `--force`. A run with no driven rate writes no such file |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | each family's true genealogy, in `genomes/gene_trees/`. Leaves are `n<species>_g<copy>`; **internal nodes are labelled `<event>_n<species>`** — `duplication_n45`, `transfer_n45`, `speciation_n45` — the event that ended that gene and the species branch it was on, which is what makes the tree readable on its own. A family with no surviving copy writes no `_extant` file |
| Tip names | `names.tsv` | TSV | external tree | written **only when the tree came from `--from` with its own tip labels** — `node` · `name`, mapping ZOMBI's `n<id>` back to the labels you supplied. Every other output names nodes `n<id>`, so this is the join back to your taxa; `profiles.tsv`'s columns key on the same ids |
| Family origination | `.gene_trees[f].origination` | float | Python | when the family was founded — where its gene tree's root branch begins |

### One row per gene-tree edge, not per event

`genome_events.tsv` holds **one row for every gene-tree edge that begins**, not one row per event.
An event that leaves two descendants therefore writes two rows, and an event that leaves one writes
one:

| Kind | Rows | Why |
|------|------|-----|
| `duplication` | **2** | the parent copy ends; two copies begin, both on `lineage` |
| `transfer` | **2** | the donor copy ends; the continuation begins on the donor branch, the transferred copy on `recipient` |
| `speciation` | **2** | the parent copy ends at the split; one copy begins in each daughter species |
| `loss` | 1 | the copy ends with no descendant |
| `origination` | 1 | a founding copy begins with no parent |

You do not have to do this arithmetic: `genome_summary.json` carries the deduplicated counts, and
`event_rows` beside them, so the two numbers are both named rather than one of them being a trap.

**Counting rows by kind therefore counts edges, not events, and doubles duplications, transfers and
speciations.** The `event` column is there so you never have to: it names the gene copy whose fate
the event is, so the rows of one event agree on it and, within a kind, it is exactly unique. Count
distinct values of it and the answer is right for every kind:

```bash
# transfers that happened — NOT `grep -c transfer`, which returns twice this
awk -F'\t' '$2=="transfer" {print $9}' genome_events.tsv | sort -u | wc -l
```

`event` repeats `parent` on the kinds that end a copy and `copy` on the two that do not, so one rule
covers all five and no reader has to know which column a given kind groups on. Group on it rather
than on `time`: times are floats, and pairing rows by float equality is both fragile and prone to
joining two events that happened to fire at once. Compare across kinds only as `(kind, event)` — a
copy that is originated and later lost names both of those events.

The edge is the unit because the file has to reconstruct the gene tree: every branch needs its own
row to carry the copy that starts it. Note that `chromosome_events.tsv`, in the same directory,
records **one row per event** — its rows are rearrangements, which begin no gene-tree edge.

## Genomes, ordered — `simulate_genomes_ordered`

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the run's whole history, in one time-ordered table.** The gene genealogy as at the family resolution, plus **where** each event happened and the ancestry-neutral rearrangements as kinds of their own — `time` · `kind` · `lineage` · `family` · `copy` · `parent` · `recipient` · `donor` · `event` · `dest_lineage` · `chromosome` · `position` · `length` · `dest_chromosome` · `dest_position` · `flipped`¹. `position` / `length` are coordinates in that branch's own genome just before the event, as `gene_order` numbers it, and are filled **once per event** — on its first row — because the arc belongs to the event, not to each copy it touched. Filtering on a non-empty `position` therefore gives one row per event that moved genes, which is what a replay walks. A transfer writes a row on each branch and each names the whole edge: `donor` is on both, and the departing row adds `dest_lineage` for where the material went |
| Profiles | `profiles.tsv` | TSV | yes | family × extant-species copy counts |
| Gene order | `gene_order.tsv` | TSV | yes | signed gene order of **every node**, ancestors included — `lineage` · `chromosome` · `position` · `strand` · `family` · `copy` |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `strand` · `family` · `copy`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — `time` · `kind` · `lineage` · `parents` · `children` |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | as at the family resolution — position is orthogonal to genealogy |

¹ a run is named by `start` (its first position, in the chromosome's frame just before the event) and
`length` (how many genes it covered), counted rightwards from `start` and **wrapping past position 0 on
a circular chromosome** — so `start + length` greater than the chromosome's gene count means the run
crossed the origin. `dest_position` is an index into what was left after the run was excised.

## Genomes, nucleotide — `simulate_genomes_nucleotide`

From `zombi2 genomes --resolution nucleotide` or `result.write(dir, outputs=[...])`.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Event log | `genome_events.tsv` | TSV | yes | **the genealogy, in the format every resolution writes** — the same columns as at the family resolution, `time` · `kind` · `lineage` · `family` · `copy` · `parent` · `recipient` · `donor` · `event`, one row per gene-tree edge. So one reader serves all three resolutions and `zombi2 tools` works here unchanged. `family` is the declared gene (else the recovered root-block), and `copy` is a **gene** id — the token the gene trees, alignments and homology tables use. It is derived onto the root-block partition, where a copy either covers a block in full or does not touch it, which is what makes a duplication a bifurcation here |
| Block events | `block_events.tsv` | TSV | yes | **this resolution's own record**, which has no counterpart elsewhere: the copy-lineage log over ancestral intervals, plus the ancestry-neutral rearrangements — `time` · `kind` · `lineage` · `chromosome` · `copy` · `parent` · `recipient` · `source` · `start` · `end` · `position` · `length` · `dest_chromosome` · `dest_position` · `flipped`. One row per **ancestral interval** an event touched, so an event spanning several blocks writes several rows — and a duplication mints a child *without* ending its parent, because a copy lineage covers an extent and the event covers a sub-extent. That is the right model for sequence and the wrong shape for a gene tree, which is why the genealogy above is a separate table rather than this one renamed. `kind` `initial` marks the initial genome being laid down at time 0 (one per replicon) — what the run *starts* with, kept distinct from `origination`, the de-novo births the `--origination` rate makes. `source`/`start`/`end` are **ancestral** coordinates while `position`/`length` are **physical** ones along the chromosome, as `blocks.tsv` numbers it: different frames, so different columns. This is what `read_nucleotide_genomes` replays |
| Blocks | `blocks.tsv` | TSV | yes | every node's genome as its block mosaic, ancestors included — `lineage` · `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. The rows of one chromosome tile it end to end from 0. The largest file this level writes: blocks are not kept maximal during a run, so it grows with their number × every node |
| Genes | `genes.tsv` | TSV | yes | the declared genes in initial coordinates — `family` · `name` · `source` · `start` · `end` · `strand` (the **coding** strand). Header-only when none were declared |
| Initial sequence | `initial_sequence.fasta` | FASTA | yes¹ | the initial DNA the run was given (`--fasta`), one `>source<n>` record per replicon. Written only when a FASTA was supplied; it is what lets a separate `zombi2 sequences` run found its blocks from the real sequence |
| Initial genome | `initial_genome.tsv` | TSV | yes | the genome the run **started** with, at the start of the root branch — `chromosome` · `position` · `source` · `start` · `end` · `strand` · `copy` · `gene`. Its own file, with no `lineage` column, because it belongs to no node: every `lineage` elsewhere is a node, and a node sits at the *end* of its branch |
| Chromosome events | `chromosome_events.tsv` | TSV | yes | chromosome-network edges — same format as ordered |
| Gene trees | `gene_tree_fam<f>_complete.nwk` · `…_extant.nwk` | Newick | yes | one tree per declared gene (else per recovered root-block), in `gene_trees/` |
| GFF | `genome_<lineage>.gff` | GFF3 | yes | one file per node, in `gff/`: that genome's **genes** in its own coordinates — the annotation to read beside the sequence level's `genome_<lineage>.fasta`, which it names its sequences to match |
| BED | `genome_<lineage>.bed` | BED | yes | one file per node, in `bed/`: that genome's **blocks**, spacer included, each named by the ancestral interval it descends from — the ancestry as a browser track |

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
`n<species>_g<copy>`, so a phylogram's tips pair with its alignment and its internal nodes with the
ancestral sequences.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Alignments | `fam<f>.fasta` | FASTA | yes | one row per extant gene copy — nucleotides or amino acids, following the model. They go in `alignments/`, which is what lets the name be this short |
| Phylograms | `phylogram_fam<f>_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the gene tree each family's sequences were drawn along, in `phylograms/` |
| Ancestral | `sequences_ancestral_fam<f>.fasta` | FASTA | no | the sequence at every node that is not an extant tip: internal nodes, and the tips where a copy was lost or its species died. One per family, so they go in `ancestral/` |
| Founding | `sequences_founding.fasta` | FASTA | no | one record `fam<f>` per family — the sequence it originated with, where its phylogram's root branch begins |
| Clock species tree | `clock_species_tree_complete.nwk` · `…_extant.nwk` | Newick (subs/site) | yes | the species tree with its branches in substitutions/site — the molecular clock made visible |
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

The `zombi2 tools` commands write their own files — the homology matrix and the reconciliation/scoring
outputs — catalogued in Appendix C.
