# Genomes III: nucleotide

The two previous chapters counted genes. This one counts base pairs. At the **nucleotide** resolution a chromosome is no longer a list of gene tokens but a **coordinate axis of DNA**, and an event is an arc on that axis: an inversion reverses 600 bp, a loss deletes 900 bp, a duplication copies 2 kb in tandem. Genes still exist and still get gene trees, but they are now stretches of that axis with a start and an end, and the DNA between them is simulated too.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_nucleotide

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=2)
g = simulate_genomes_nucleotide(
    tree, root_length=3000, genes=3, gene_length=400,
    inversion=1.0, inversion_length=600,
    duplication=0.3, duplication_length=300, loss=0.3, loss_length=300, seed=2)
```

This seeds a 3000 bp circular chromosome carrying three 400 bp genes, evenly spaced, and evolves it down the tree.

## Blocks, genes and intergenes

A chromosome is stored as an ordered list of **blocks**. A block is a run of DNA with one unbroken ancestry: the interval `[start, end)` on some ancestral **source**, read forward (`+`) or reverse-complemented (`−`). Events only ever *split* blocks, never merge them, so the block boundaries you see at a leaf are the accumulated breakpoints of its whole history.

Blocks come in two kinds, and the distinction runs through the rest of the chapter:

- a **gene** is a *declared, indivisible* block — one family, one id, **never split**. It carries a gene tree.
- an **intergene** is the spacer between genes. It fragments freely into as many blocks as events dictate. Its ancestry is tracked, but it gets no tree of its own.

A genome is therefore an alternating chain of intergenes and genes. Declare no genes and the whole chromosome is one big intergene, which is the uniform-sequence model.

Reading two leaves of the run above shows what the events did:

```
leaf n2, chromosome 2 (circular), 3000 bp        leaf n5, chromosome 5 (circular), 3000 bp
  [   0,  600) +   600 bp  intergene               [   0,  599) +   599 bp  intergene
  [ 600, 1000) +   400 bp  gene 1                  [1174, 1374) +   200 bp  intergene
  [1000, 1144) +   144 bp  intergene               [1000, 1144) −   144 bp  intergene
  [1144, 1374) −   230 bp  intergene               [ 600, 1000) −   400 bp  gene 1
  [1374, 1404) +    30 bp  intergene               [ 599,  600) −     1 bp  intergene
  [2000, 2262) −   262 bp  intergene               [1144, 1174) −    30 bp  intergene
  [1600, 2000) −   400 bp  gene 2                  [1374, 1600) +   226 bp  intergene
  [1404, 1600) −   196 bp  intergene               [1600, 2000) +   400 bp  gene 2
  [2262, 2600) +   338 bp  intergene               [2000, 2021) +    21 bp  intergene
  [2600, 3000) +   400 bp  gene 3                  [2021, 2319) −   298 bp  intergene
                                                   [2319, 2600) +   281 bp  intergene
                                                   [2600, 3000) +   400 bp  gene 3
```

Both leaves still carry all three genes. In `n2` an inversion covered gene 2, which now reads on the `−` strand with the spacer around it reversed; in `n5` a different inversion covered gene 1 instead. The coordinates are what make this readable: every block still names where it came from in the root, so `[600, 1000)` is gene 1 wherever it turns up and whichever way it points.

## Genes are never split

This is the rule the whole resolution is built on. An event either **engulfs a gene whole** or leaves it alone; a breakpoint never falls strictly inside one. So an event does not pick an arc and then clean up afterwards. It **nucleates at a random intergenic position**, and its far end is drawn from the positions where a breakpoint is actually legal, weighted toward the extent you asked for. Nothing is clipped, nothing is snapped, and no event is silently dropped.

Two consequences follow, and both surprise people.

**Gene turnover is emergent and size-dependent.** A gene changes copy number only when an event engulfs it entirely, and engulfing a big gene takes a big event. Take twenty genes in a 400 kb genome, fix the deletions at `loss=30.0, loss_length=2000`, and vary nothing but the gene size:

| gene length | genes left per leaf (of 20) |
|---|---|
| 100 bp | 14.0 |
| 500 bp | 14.8 |
| 2 000 bp | 16.8 |
| 5 000 bp | 18.4 |
| 10 000 bp | 19.8 |

Genes much longer than the deletions are essentially never lost. There is no gene-loss rate to set directly here: **to get more gene turnover, use larger events.**

**The realised extent is not the extent you asked for.** It is quantised to the legal breakpoints, and on a gene-dense genome the difference is large. Take thirty-one 3000 bp genes in a 100 kb genome — 93% genic, so the spacers are about 200 bp — and ask for an inversion of:

| asked | realised |
|---|---|
| 500 bp | 59 |
| 1 500 bp | 1 000 |
| 3 000 bp | 2 916 |
| 10 000 bp | 10 315 |

A 500 bp event has nowhere to go but inside a spacer, so it comes out at 59 bp. Long events land near what you asked for, because they can span whole genes. This conditioning *is* the model; what it does not do is quietly eat the event **rate** along with the length.

## The events

The event names are the ones you already know from Chapters 4 and 5, and they mean the same things. What changes is that each one now acts on an arc of DNA rather than a run of gene tokens, so each rate is paired with a mean length in base pairs, `<event>_length`:

| Event | Length knob | What it does |
|---|---|---|
| `inversion` | `inversion_length` | reverse-complement an arc in place |
| `transposition` | `transposition_length` | move an arc elsewhere on the same chromosome |
| `translocation` | `translocation_length` | move an arc to a different chromosome |
| `duplication` | `duplication_length` | copy an arc in tandem |
| `loss` | `loss_length` | delete an arc; never empties a chromosome |
| `transfer` | `transfer_length` | copy an arc into a contemporaneous recipient |
| `origination` | `origination_length` | lay down a **new gene** on fresh DNA |

Transposed and translocated arcs land inverted with probability `inversion_probability`. Transfer takes the same `transfer_to` (`"uniform"` or `"distance"`) and `self_transfer` arguments as Chapter 4, and is always additive — the donor keeps its copy, and replacing transfer is not yet implemented here. Origination always mints a gene, never plain spacer.

The first three events leave ancestry untouched and are recorded in `.rearrangements`; the last four change it and are recorded in `.events`, which is the log the gene trees are recovered from.

**One thing genuinely differs from the earlier chapters: these rates are per lineage, not per copy.** The rate says how often a lineage does the thing; the extent says how much DNA it touches. Chapter 5 counts the same events per gene, because a gene is a unit you can count chances on: a genome of 4 000 genes offers 4 000 places for a run to start. A nucleotide is not that unit. Counting per base pair would multiply the event count by the length of the genome, so the rate you type would need a string of leading zeros to mean anything — and how many zeros would depend on how long the genome happened to grow. Per lineage keeps the number you type readable and the event count flat as the genome grows; the same tree at `inversion=5.0`, with the genome a hundred times longer each row, gives:

```
   10 000 bp  ->  77 inversions
  100 000 bp  ->  77 inversions
1 000 000 bp  ->  77 inversions
```

The chromosome tier below is the exception: `fission`, `fusion` and `chromosome_loss` are counted **per chromosome**, and `chromosome_origination` per lineage.

Rates here are plain numbers. The `scope(base) × modifiers` grammar of the other levels is not yet wired at this resolution, so an `OnTime` skyline is not available.

## Chromosomes

The karyotype and the chromosome tier work exactly as in Chapter 5, and produce the same `chromosome_events` edge list — fission is a bifurcation, fusion the reticulation, origination a root, loss a leaf. What is new is that a chromosome now has a **size and a shape**: `chromosomes=3, root_length=4000` seeds three equal replicons, and a list of `(length, topology)` pairs seeds heterogeneous ones.

```python
g = simulate_genomes_nucleotide(
    tree, chromosomes=[(5_000_000, "circular"), (60_000, "circular")],
    genes=0, inversion=1.0, seed=1)                 # a chromosome and a plasmid
```

One rule is specific to this resolution: **a chromosome never exists without a gene.** A new replicon is born carrying one, and any event that would strip a chromosome of its last gene — a loss that would take it, a translocation carrying it away, a fission splitting off a geneless half — simply does not happen. This guarantees every lineage always holds some DNA; the rule is vacuous when no genes are declared.

## Seeding the root genome

The root genome is declared, in one of two ways.

**Evenly spaced genes** — `genes=N, gene_length=L` lays down `N` genes of `L` bp on each replicon, spreading the leftover DNA evenly between them. Good for controlled experiments, since gene density is then a number you set.

**A GFF file** — `gff="genome.gff"` takes exact coordinates from a real annotation. `##sequence-region` gives each replicon's extent, `gene` features give coordinates, strand and name, and other feature types are ignored. Names land in `result.gene_names`, so you can follow a named gene through the run. `gff=` and `genes=` are mutually exclusive; a GFF already declares the genes.

## The `NucleotideGenomesResult` object

`simulate_genomes_nucleotide` returns a **`NucleotideGenomesResult`**:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages included.
- `.genomes` — a dict from node id to that node's `NucleotideGenome`, a list of `Chromosome`s, each a list of `Block`s.
- `.events` — the copy-lineage genealogy: origination, loss, duplication, transfer, speciation.
- `.rearrangements` — the ancestry-neutral log: inversion, transposition, translocation.
- `.chromosome_events` — the chromosome network, as in Chapter 5.
- `.gene_spans` — `{family: (source, start, end)}`, where each declared gene sits in root coordinates.
- `.gene_names`, `.gene_strands` — a named gene's family id, and its coding strand, from a GFF.
- `.gene_trees` — one recovered gene tree per gene family, for the families that survive in at least one extant leaf.
- `.root_blocks` — the recovered root partition: the maximal never-cut intervals that some node still carries. Cut at **every** node's breakpoints, not just the survivors', which is what lets any node's genome be rebuilt.
- `.block_trees` — a recovered tree for **every** root block, spacer as well as gene, keyed by its index in `.root_blocks`.
- `.block_of(family)` — the block index a declared gene family occupies: the join between the two numbering schemes here, since `.gene_spans` and `.gene_trees` are keyed by family id while `.root_blocks` and `.block_trees` are keyed by block index. Both are plain integers, so mixing them up is silent; this is how you avoid it.
- `.seed`.

and four ways to read one node's genome, at four grains:

```python
g.mosaic(2)        # per block:      {chromosome: [(source, start, end, strand), ...]}
g.trace_back(2)    # per nucleotide: {chromosome: [(source, position, strand), ...]}
g.ancestry(2)      # the multiset of ancestral (source, position) it still carries
g.assembly(2)      # per piece:      {chromosome: [(block, gene, start, end, strand), ...]}
```

`ancestry` is the invariant worth knowing: rearrangements conserve it exactly, so an inversion-only run leaves every leaf holding a permutation of the root sequence. Loss makes it a subset, and duplication, transfer and origination add to it.

## Every block has a tree, so a genome can be rebuilt

A gene is a block that never splits, which is why it has a tree. But *no* block ever splits — a block is a run of one unbroken ancestry, and an event only ever cuts one in two, never merges two into one. So the recovery that builds a gene tree works on any block, and `.block_trees` runs it on all of them. The spacer between the genes gets a genealogy on exactly the same footing as the genes.

That is what makes a genome rebuildable rather than a list of loci. `.assembly(node)` says how a node's genome is put together out of those blocks: for each chromosome, the pieces in physical order, each one a stretch `[start, end)` of root block `block`, carrying gene id `gene` in that block's tree, read forward (`strand` `+1`) or reverse-complemented (`-1`). Pair each piece with a sequence evolved down its block's tree and you have the genome; that is precisely what Chapter 7 does.

A piece is always a whole block. That follows from cutting the partition at **every** node's breakpoints rather than only the survivors': a node's own boundaries are then all in the partition, so each of its blocks is a whole number of root blocks — one piece each, in physical order, running down the source where the block is inverted. It also means no node is ever missing a block, so `assembly` works at an extinct lineage and at the root as readily as at a surviving tip.

## Usage from Python

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_nucleotide

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=6, seed=4)

# rearrangement only — every leaf is a permutation of the root sequence
g = simulate_genomes_nucleotide(
    tree, root_length=100_000, inversion=5.0, inversion_length=1000, seed=4)

# a genic genome with real turnover: genes duplicate, are lost, transfer, and arise
g = simulate_genomes_nucleotide(
    tree, root_length=6000, genes=6, gene_length=400,
    duplication=2.0, duplication_length=900, loss=2.0, loss_length=900,
    transfer=1.0, transfer_length=900, transfer_to="distance",
    origination=0.5, origination_length=400, seed=4)

# a karyotype that splits and merges, with a plasmid
g = simulate_genomes_nucleotide(
    tree, chromosomes=3, root_length=4000, genes=6, gene_length=200,
    fission=0.2, fusion=0.2, chromosome_origination=0.05, chromosome_loss=0.05,
    inversion=1.0, seed=5)

# a real annotation as the root genome
g = simulate_genomes_nucleotide(
    tree, gff="ecoli.gff", inversion=2.0, inversion_length=5000,
    loss=1.0, loss_length=3000, seed=1)

# the outputs
g.gene_spans                            # where each gene sits, in root coordinates
family = min(g.gene_trees)              # gene_trees holds only the families that survive
g.gene_trees[family].to_newick("extant")
leaf = next(n.id for n in g.complete_tree.extant())
g.mosaic(leaf)                          # that leaf's genome, block by block
g.chromosome_events                     # the chromosome network
```

A gene family lost in **every** extant lineage still gets a tree — it is still history — but a complete one only: `.gene_trees[fam].extant` is `None`, and no `_extant.nwk` is written for it. Only a family deleted from every node in the tree has no tree at all, and it will still appear in `.gene_spans`, because it was declared. Ask `.gene_trees` what it holds rather than assuming a declared family is in it.

The second run starts from six genes and ends with each leaf carrying between nine and twelve families: some seeded genes duplicated, one was lost in every lineage, transfer moved copies sideways, and origination added new families along the way. Every gene tree it produces is a true genealogy recovered from the copy-lineage log, exactly as at the other two resolutions.

## Usage from the CLI

The nucleotide resolution is `--resolution nucleotide`. It takes the same event rates as Chapter 5 and adds two things: how to seed the genome, and how long an event is in base pairs.

```bash
# an evenly spaced seed: 5 kb, six genes of 300 bp, with inversions averaging 400 bp
zombi2 genomes out/ --resolution nucleotide \
  --root-length 5000 --genes 6 --gene-length 300 \
  --inversion 1.0 --inversion-length 400 --duplication 0.3 --loss 0.3 --seed 1

# or start from a real genome: the GFF declares the replicons and the genes
zombi2 genomes out/ --resolution nucleotide \
  --gff ecoli.gff --inversion 0.5 --loss 0.4 --loss-length 900 --seed 1
```

Every event kind has its own `--<event>-length`, the mean of a geometric draw in base pairs: `--inversion-length`, `--loss-length`, `--duplication-length`, `--transfer-length`, `--transposition-length`, `--translocation-length`, `--origination-length`.

Two flags from the other resolutions are refused here rather than ignored, because this engine has neither: `--initial-families` (it is seeded from a sequence, not a family count) and `--replacement` (its transfers are additive). Rates must be plain numbers — the `scope(base) × modifiers` grammar of SPEC §5 is not wired at this resolution, so a modifier expression is an error rather than something quietly dropped.

## Outputs

```bash
zombi2 genomes out/ --resolution nucleotide --root-length 5000
```

```
out/genome_events.tsv    the copy-lineage genealogy (the source of truth)
out/genes.tsv            where each gene sits in the root, and on which strand
out/blocks.tsv           every node's genome as its block mosaic
out/gene_trees/          one Newick per family, complete and extant
out/rearrangements.tsv   inversions, transpositions, translocations
out/chromosome_events.tsv  the chromosome network's edges
```

Everything is written by default, as at the ordered resolution — the point of this level is a history you can replay, and a history missing a table is not one. `blocks.tsv` is the big file: blocks are not kept maximal during a run, so a rearrangement-heavy genome carries far more of them than it has distinct ancestral runs, and the file grows with that number times every node. Narrow the set with `--write` when that matters — `--write events genes` for the two small ones.

`genes.tsv` is the gene declaration as the run saw it, de-novo originations included:

```
family  name  source  start  end   strand
1             0       534    834   1
2             0       1368   1668  1
3             0       2201   2501  1
```

The event log is wider than the one Chapters 4 and 5 write, because an event here is an arc rather than a gene: each row names the ancestral interval it touched, so one event that spanned several blocks writes several rows sharing a `time` and `kind`. That also means it is **not** the log `zombi2 sequences` replays — that command reads the unordered or ordered one, and says so if handed this.

Its width is also why there is no `genome_event_positions.tsv` here. Chapter 5 needs that companion file because its event log is the position-blind one Chapter 4 writes, so the coordinates have to live somewhere else; here an event *is* an interval, and `source`, `start` and `end` are columns of the log itself.

`rearrangements` and `chromosome_events` are the same tables as Chapter 5. Everything else stays in Python — `.trace_back(node)` for per-nucleotide ancestry, `.gene_trees` for the recovered genealogies, `.mosaic(node)` for a genome block by block. Appendix B catalogues the lot.
