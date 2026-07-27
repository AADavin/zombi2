# Genomes III: nucleotide

At the **nucleotide** resolution a chromosome is no longer a list of gene tokens but a **coordinate axis of DNA**, and an event is an arc on that axis: an inversion reverses 600 bp, a loss deletes 900 bp, a duplication copies 2 kb in tandem. Genes still exist and still get gene trees, but they are now stretches of that axis with a start and an end, and the DNA between them is simulated too.

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=2)
g = genomes.simulate_genomes_nucleotide(
    tree, root_length=3000, genes=3, gene_length=400,
    inversion=1.0, inversion_extent=600,
    duplication=0.3, duplication_extent=300, loss=0.3, loss_extent=300, seed=2)
```

This starts the run from a 3000 bp circular chromosome carrying three 400 bp genes, evenly spaced, and evolves it down the tree.

## Blocks: genes and intergenes

A chromosome is stored as an ordered list of **blocks**. A block is a run of DNA with one unbroken ancestry: the interval `[start, end)` on some ancestral **source**, read forward (`+`) or reverse-complemented (`−`). Events only ever *split* blocks, never merge them, so the block boundaries you see at a leaf are the accumulated breakpoints of its whole history.

Blocks come in two kinds:

- a **gene** is a *declared, indivisible* block — one family, one id, **never split**. It carries a gene tree.
- an **intergene** is the spacer between genes. It fragments freely into as many blocks as events dictate. 

A genome is therefore an alternating chain of intergenes and genes. Either extreme is legal. Declare no genes and the whole chromosome is one big intergene — the uniform-sequence model. Fill the whole replicon with genes and there is no spacer at all; that evolves too, because events break at the joins *between* genes (see *Genes are never split*, below).

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

An event either **engulfs a gene whole** or leaves it alone; a breakpoint never falls strictly inside one. So an event does not pick an arc and then clean up afterwards. Both of its ends are drawn **directly from the positions where a breakpoint is legal**. A genome can therefore be **all gene, with no spacer at all**: ten 100 bp genes in 1000 bp is a legal genome, and it evolves. Its breakpoints simply all fall at the joins between genes, so genes are inverted, moved, duplicated and lost whole. Genes may sit flush; they are not required to leave a gap.

## Extents

There is one important consequence: **the realised extent is not always the extent you asked for.** It is quantised to the legal breakpoints, and on a gene-dense genome the difference is large. Take thirty-one 3000 bp genes in a 100 kb genome — 93% genic, so the spacers are about 200 bp — and ask for an inversion of:

| asked | realised |
|---|---|
| 500 bp | 59 |
| 1 500 bp | 1 000 |
| 3 000 bp | 2 916 |
| 10 000 bp | 10 315 |

A 500 bp event has nowhere to go but inside a spacer, so it comes out at 59 bp. Long events land near what you asked for, because they can span whole genes.

**The correction runs both ways.** When every legal end lies *further* out than you asked — a long gene sitting just past the start — the arc snaps to the nearest legal breakpoint, which is **longer** than the extent you set. That is the `10 000 → 10 315` row above. When the genome cannot give what you asked at all, the arc is capped by the replicon and comes out shorter. Either way the event still fires: an extent is a request, and the genome answers it.

The one case that yields no event is degenerate — a replicon with no legal end within reach at all, such as one under 2 bp — where the event is skipped rather than forced.

## A note on rates

**Every rate here is per lineage.** The rate sets how often a lineage does the event; the extent (above) sets how much DNA it touches. Keeping the rate per lineage means the number you type reads the same whatever the genome's size: a rate counted per base pair would rise as the genome grew, so `inversion=5.0` would mean one thing at 10 kb and another at 1 Mb. Per lineage, the event count stays flat as the genome grows — the same tree at `inversion=5.0`, with the genome a hundred times longer each row, gives:

```
   10 000 bp  ->  77 inversions
  100 000 bp  ->  77 inversions
1 000 000 bp  ->  77 inversions
```

The chromosome tier below is the exception: `fission`, `fusion` and `chromosome_loss` are counted **per chromosome**, and `chromosome_origination` per lineage.

Rates here take the same written form as everywhere else — `scope(base) × modifiers` — and the scopes above are the defaults, so a bare number stays a bare number. The **skyline** works: `inversion = 5.0 * OnTime({0: 1.0, 3: 0.2})` drops the inversion rate fivefold at time 3, and the run re-reads its rates at each step rather than racing past it.

So does **conditioning**. Every rate here takes a `DrivenBy`, so a trait can drive how much DNA a lineage sheds — genome reduction as it is usually meant — and can drive the rearrangements too:

```python
loss = 0.8 * mod.DrivenBy(habitat, {"host": 20.0, "free": 0.5})
```

**The extent takes the same modifiers**, and that is a different statement:

```python
loss        = 0.8 * mod.DrivenBy(habitat, {"host": 20.0, "free": 0.5})   # deletes more often
loss_extent = 150 * mod.DrivenBy(habitat, {"host": 6.0,  "free": 1.0})   # deletes in bigger chunks
```

The first raises how often a host-restricted lineage deletes; the second raises how much each deletion takes. Set both and they multiply, so the DNA shed per unit time goes up by the product, not the sum. Which of the two you want is a modelling choice, and having only one knob would force you to describe genome reduction with the wrong one.

A modifier on an *extent* is read at the instant an event fires, so unlike the same modifier on a rate it changes no rate and adds no step to the run's clock. Chapter 9 covers what a driver is and how to grow one. Anything the engine has *not* wired raises rather than being quietly ignored.

## The initial genome

The **initial genome** — the genome the run starts from, at time 0, before any event — is declared in one of two ways. (It is not the same as the genome at the root *node*: the root branch is real simulated time, so that one already carries stem events. See `.initial_genome` below.)

**Evenly spaced genes** — `genes=N, gene_length=L` lays down `N` genes of `L` bp on each replicon, spreading the leftover DNA evenly between them. Good for controlled experiments, since gene density is then a number you set.

**A GFF file** — `gff="genome.gff"` takes exact coordinates from a real annotation. `##sequence-region` gives each replicon's extent, `gene` features give coordinates, strand and name, and other feature types are ignored. Names land in `result.gene_names`, so you can follow a named gene through the run. `gff=` and `genes=` are mutually exclusive; a GFF already declares the genes.

## The `NucleotideGenomesResult` object

`simulate_genomes_nucleotide` returns a **`NucleotideGenomesResult`**:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages included.
- `.genomes` — a dict from node id to that node's `NucleotideGenome`, a list of `Chromosome`s, each a list of `Block`s.
- `.initial_genome` — the genome the run **started** with, at the root lineage's origination. It is not `.genomes[root]`: a node sits at the **end** of its branch, and the root branch is real simulated time, so events happen along it. Written to its own `initial_genome.tsv`, with no `lineage` column, because it belongs to no node. It votes on the root partition like every other genome, so `.initial_assembly()` rebuilds it too.
- `.events` — the copy-lineage genealogy: origination, loss, duplication, transfer, speciation.
- `.rearrangements` — the ancestry-neutral log: inversion, transposition, translocation.
- `.chromosome_events` — the chromosome network, as in Chapter 5.
- `.gene_spans` — `{family: (source, start, end)}`, where each declared gene sits in initial coordinates.
- `.gene_names`, `.gene_strands` — a named gene's family id, and its coding strand, from a GFF.
- `.gene_trees` — one recovered gene tree per gene family, for the families that survive in at least one extant leaf.
- `.root_blocks` — the recovered root partition: the maximal never-cut intervals that some node still carries. Cut at **every** node's breakpoints, not just the survivors', which is what lets any node's genome be rebuilt.
- `.block_trees` — a recovered tree for **every** root block, spacer as well as gene, keyed by its index in `.root_blocks`.
- `.initial_assembly()` — the same for `.initial_genome`, as `(block, strand)` pairs. No gene id: the initial genome predates every event, so each block has exactly one sequence there.
- `.block_of(family)` — the block index a declared gene family occupies: the join between the two numbering schemes here, since `.gene_spans` and `.gene_trees` are keyed by family id while `.root_blocks` and `.block_trees` are keyed by block index. Both are plain integers, so mixing them up is silent; this is how you avoid it.
- `.seed`.

and four ways to read one node's genome, at four grains:

```python
g.mosaic(2)        # per block:      {chromosome: [(source, start, end, strand), ...]}
g.trace_back(2)    # per nucleotide: {chromosome: [(source, position, strand), ...]}
g.ancestry(2)      # the multiset of ancestral (source, position) it still carries
g.assembly(2)      # per piece:      {chromosome: [(block, gene, start, end, strand), ...]}
```

`ancestry` is the invariant worth knowing: rearrangements conserve it exactly, so an inversion-only run leaves every leaf holding a permutation of the initial sequence. Loss makes it a subset, and duplication, transfer and origination add to it.

## Usage from Python

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=6, seed=4)

# rearrangement only — every leaf is a permutation of the initial sequence
g = genomes.simulate_genomes_nucleotide(
    tree, root_length=100_000, inversion=5.0, inversion_extent=1000, seed=4)

# a genic genome with real turnover: genes duplicate, are lost, transfer, and arise
g = genomes.simulate_genomes_nucleotide(
    tree, root_length=6000, genes=6, gene_length=400,
    duplication=2.0, duplication_extent=900, loss=2.0, loss_extent=900,
    transfer=1.0, transfer_extent=900, transfer_to="distance",
    origination=0.5, origination_extent=400, seed=4)

# a karyotype that splits and merges, with a plasmid
g = genomes.simulate_genomes_nucleotide(
    tree, chromosomes=3, root_length=4000, genes=6, gene_length=200,
    fission=0.2, fusion=0.2, chromosome_origination=0.05, chromosome_loss=0.05,
    inversion=1.0, seed=5)

# a real annotation as the initial genome
g = genomes.simulate_genomes_nucleotide(
    tree, gff="ecoli.gff", inversion=2.0, inversion_extent=5000,
    loss=1.0, loss_extent=3000, seed=1)

# the outputs
g.gene_spans                            # where each gene sits, in initial coordinates
family = min(g.gene_trees)              # gene_trees holds only the families that survive
g.gene_trees[family].to_newick("extant")
leaf = next(n.id for n in g.complete_tree.extant())
g.mosaic(leaf)                          # that leaf's genome, block by block
g.chromosome_events                     # the chromosome network
```

A gene family lost in **every** extant lineage still gets a tree — it is still history — but a complete one only: `.gene_trees[fam].extant` is `None`, and no `_extant.nwk` is written for it. Only a family deleted from every node in the tree has no tree at all, and it will still appear in `.gene_spans`, because it was declared. Ask `.gene_trees` what it holds rather than assuming a declared family is in it.

## Usage from the CLI

The nucleotide resolution is `--resolution nucleotide`. It takes the same event rates as Chapter 5 and adds two things: how to set up the initial genome, and how long an event is in base pairs.

```bash
# an evenly spaced initial genome: 5 kb, six genes of 300 bp, with inversions averaging 400 bp
zombi2 genomes out/ --resolution nucleotide \
  --root-length 5000 --genes 6 --gene-length 300 \
  --inversion 1.0 --inversion-extent 400 --duplication 0.3 --loss 0.3 --seed 1

# or start from a real genome: the GFF declares the replicons and the genes,
# and a paired FASTA supplies the actual DNA those coordinates hold
zombi2 genomes out/ --resolution nucleotide \
  --gff ecoli.gff --fasta ecoli.fasta --inversion 0.5 --loss 0.4 --loss-extent 900 --seed 1
```

Every event kind has its own `--<event>-extent`, the mean of a geometric draw in base pairs: `--inversion-extent`, `--loss-extent`, `--duplication-extent`, `--transfer-extent`, `--transposition-extent`, `--translocation-extent`, `--origination-extent`.

## Outputs

```bash
zombi2 genomes out/ --resolution nucleotide --root-length 5000
```

```
out/genome_events.tsv    the whole history: the copy-lineage genealogy and the
                         rearrangements — in time order
out/genes.tsv            where each gene sits in the root, and on which strand
out/blocks.tsv           every node's genome as its block mosaic
out/initial_genome.tsv   the genome the run started with
out/gene_trees/          one Newick per family, complete and extant
out/genome_<lineage>.gff · .bed   every genome's genes and blocks
out/chromosome_events.tsv  the chromosome network's edges
```

Everything is written by default.
