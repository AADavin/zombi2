# Genomes III: nucleotide

At the **nucleotide** resolution a chromosome is a **coordinate axis of DNA** rather than a list of gene tokens, and events have an extension measured in base pairs: an inversion reverses 600 bp, a loss deletes 900 bp, a duplication copies 2 kb in tandem. Genes still exist and still get gene trees, but they are stretches of that axis with a start and an end, and the DNA between them is simulated too.

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=2)
g = genomes.simulate_genomes_nucleotide(
    tree, root_length=3000, genes=3, gene_length=400,
    inversion=1.0, inversion_extent=600,
    duplication=0.3, duplication_extent=300, loss=0.3, loss_extent=300, seed=2)
```

This starts the run from a 3000 bp circular chromosome carrying three 400 bp genes, evenly spaced, and evolves it down the tree.

## Genes and intergenes

A genome here is DNA, and its DNA is of two kinds:

- a **gene** is *declared* and *indivisible*: one family, one id, never cut in two. It carries a gene tree.
- an **intergene** is the spacer between genes. Nothing protects it, so events cut it wherever they land.

A genome is an alternating chain of intergenes and genes, and either extreme is legal. Declare no genes and the chromosome is one big intergene, the uniform-sequence model. Fill the replicon with genes and there is no spacer at all; that evolves too, because events break at the joins *between* genes (see *Genes are never split*, below).

## Extents

There is one important consequence: **the realised extent is not always the extent you asked for.** It is quantised to the legal breakpoints, and on a gene-dense genome the difference is large. Take thirty-one 3000 bp genes in a 100 kb genome, 93% genic, so the spacers are about 200 bp, and ask for an inversion of:

| asked | realised |
|---|---|
| 500 bp | 59 |
| 1 500 bp | 1 000 |
| 3 000 bp | 2 916 |
| 10 000 bp | 10 315 |

A 500 bp event has nowhere to go but inside a spacer, so it comes out at 59 bp. Long events land near what you asked for, because they can span whole genes.

**The correction runs both ways.** When every legal end lies *further* out than you asked, with a long gene sitting just past the start, the arc snaps to the nearest legal breakpoint, which is **longer** than the extent you set. That is the `10 000 → 10 315` row above. When the genome cannot give what you asked at all, the arc is capped by the replicon and comes out shorter. Either way the event still fires: an extent is a request, and the genome answers it.

The one case that yields no event is degenerate: a replicon with no legal end within reach at all, such as one under 2 bp, where the event is skipped rather than forced.

## A note on rates

**Every rate here is per lineage.** The rate sets how often a lineage does the event; the extent (above) sets how much DNA it touches. Keeping the rate per lineage means the number you type reads the same whatever the genome's size: a rate counted per base pair would rise as the genome grew, so `inversion=5.0` would mean one thing at 10 kb and another at 1 Mb. Per lineage, the event count stays flat as the genome grows. The same tree at `inversion=5.0`, with the genome a hundred times longer each row, gives:

```
   10 000 bp  ->  77 inversions
  100 000 bp  ->  77 inversions
1 000 000 bp  ->  77 inversions
```

The chromosome tier below is the exception: `fission`, `fusion` and `chromosome_loss` are counted **per chromosome**, and `chromosome_origination` per lineage. A fusion joins two chromosomes of the **same topology**, for the same reason it does at the ordered resolution: a ring and a molecule with two ends cannot become one molecule.

Rates here take the same written form as everywhere else, `scope(base) × modifiers`, and the scopes above are the defaults, so a bare number stays a bare number. The **skyline** works: `inversion = 5.0 * OnTime({0: 1.0, 3: 0.2})` drops the inversion rate fivefold at time 3, and the run re-reads its rates at each step rather than racing past it.

So does **conditioning**. Every rate here takes a `DrivenBy`, so a trait can drive how much DNA a lineage sheds, which is genome reduction as it is usually meant, and can drive the rearrangements too:

```python
from zombi2 import traits
from zombi2.rates import modifiers as mod

habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=2)
loss = 0.8 * mod.DrivenBy(habitat, {"host": 20.0, "free": 0.5})
```

**The extent takes the same modifiers**, and that is a different statement:

```python
loss        = 0.8 * mod.DrivenBy(habitat, {"host": 20.0, "free": 0.5})  # more often
loss_extent = 150 * mod.DrivenBy(habitat, {"host": 6.0,  "free": 1.0})  # bigger chunks
```

The first raises how often a host-restricted lineage deletes, the second how much each deletion takes. Set both and they multiply: the DNA shed per unit time goes up by the product, not the sum.

A modifier on an *extent* is read when an event fires, so unlike the same modifier on a rate it adds no step to the run's clock. Chapter 9 covers what a driver is and how to grow one; anything the level does not accept raises rather than being quietly ignored.

## Who receives a transfer

`transfer_to` is Chapter 4's recipient rule, and it works here unchanged: `"uniform"`, `"distance"` / `Distance(decay=)`, a `Clades(...)` kernel over named clades, or a `DrivenBy` weight read off a trait. It is not a rate: the numbers are weights normalised over the lineages alive when a transfer fires, so it says who receives and never how much transfer happens.

```python
tree6 = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=6)
flows = genomes.Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
g = genomes.simulate_genomes_nucleotide(
    tree6, root_length=6000, genes=6, gene_length=400,
    transfer=2.0, transfer_extent=900, seed=6,
    transfer_to=genomes.Clades({"A": ["n49", "n50"], "B": ["n30", "n36"]}, flows))
```

A transfer here is always **additive**, since the donor keeps its copy, so steering changes only which lineage the arc lands on. Nothing is taken from anyone, and a transfer whose every candidate weighs 0 simply does not fire.

## Genes are never split

An event either **engulfs a gene whole** or leaves it alone; a breakpoint never falls strictly inside one. So an event does not pick an arc and then clean up afterwards. Both of its ends are drawn **directly from the positions where a breakpoint is legal**. A genome can therefore be **all gene, with no spacer at all**: ten 100 bp genes in 1000 bp is a legal genome, and it evolves. Its breakpoints simply all fall at the joins between genes, so genes are inverted, moved, duplicated and lost whole. Genes may sit flush; they are not required to leave a gap.

Two leaves of the same run show what those events did. Each line is a **block**: a stretch of DNA with one unbroken ancestry, written as the interval it came from on the initial sequence, read forward (`+`) or reverse-complemented (`−`). A gene is always one block, because nothing may cut it. Spacer is not, so a run of several intergene lines in a row is simply spacer that has been cut apart and rearranged: the accumulated breakpoints of everything that happened to that lineage.

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

## The initial genome

The **initial genome**, the genome the run starts from at time 0 before any event, is declared in one of two ways. (It is not the same as the genome at the root *node*: the root branch is real simulated time, so that one already carries stem events. See `.initial_genome` below.)

**Evenly spaced genes.** `genes=N, gene_length=L` lays down `N` genes of `L` bp on each replicon, spreading the leftover DNA evenly between them. Good for controlled experiments, since gene density is then a number you set.

**A GFF file.** `gff="genome.gff"` takes exact coordinates from a real annotation. `##sequence-region` gives each replicon's extent, `gene` features give coordinates, strand and name, and other feature types are ignored. Names land in `result.gene_names`, so you can follow a named gene through the run. `gff=` and `genes=` are mutually exclusive; a GFF already declares the genes.

## The `NucleotideGenomesResult` object

`simulate_genomes_nucleotide` returns a **`NucleotideGenomesResult`**:

- `.complete_tree`, the species tree the genomes ran on, extinct lineages included.
- `.genomes`, a dict from node id to that node's `NucleotideGenome`, a list of `Chromosome`s, each a list of `Block`s.
- `.initial_genome`, the genome the run **started** with, at the root lineage's origination. It is not `.genomes[root]`: a node sits at the **end** of its branch, and the root branch is real simulated time, so events happen along it. Written to its own `initial_genome.tsv`, with no `lineage` column, because it belongs to no node. It votes on the root partition like every other genome, so `.initial_assembly()` rebuilds it too.
- `.events`, the copy-lineage genealogy: origination, loss, duplication, transfer, speciation.
- `.rearrangements`, the ancestry-neutral log: inversion, transposition, translocation.
- `.chromosome_events`, the chromosome network, as in Chapter 5.
- `.gene_spans`, `{family: (source, start, end)}`, where each declared gene sits in initial coordinates.
- `.gene_names`, `.gene_strands`, a named gene's family id and its coding strand, from a GFF.
- `.gene_trees`, one recovered gene tree per gene family: every family some node still carries, so a gene surviving only in lineages that died still gets a tree.
- `.root_blocks`, the recovered root partition: the maximal never-cut intervals that some node still carries. Cut at **every** node's breakpoints, not just the survivors', which is what lets any node's genome be rebuilt.
- `.block_trees`, a recovered tree for **every** root block, spacer as well as gene, keyed by its index in `.root_blocks`.
- `.initial_assembly()`, the same for `.initial_genome`, as `(block, strand)` pairs. No gene id: the initial genome predates every event, so each block has exactly one sequence there.
- `.block_of(family)`, the block index a declared gene family occupies: the join between the two numbering schemes here, since `.gene_spans` and `.gene_trees` are keyed by family id while `.root_blocks` and `.block_trees` are keyed by block index. Both are plain integers, so mixing them up is silent; this is how you avoid it.
- `.seed`.

and four ways to read one node's genome, at four grains:

```python
g.mosaic(2)      # per block:      {chromosome: [(source, start, end, strand), ...]}
g.trace_back(2)  # per nucleotide: {chromosome: [(source, position, strand), ...]}
g.ancestry(2)    # the multiset of ancestral (source, position) it still carries
g.assembly(2)    # per piece:      {chromosome: [(block, gene, strand), ...]}
```

`ancestry` is the invariant worth knowing: rearrangements conserve it exactly, so an inversion-only run leaves every leaf holding a permutation of the initial sequence. Loss makes it a subset, and duplication, transfer and origination add to it.

Everything in this chapter is one call. Here is the fullest version, a karyotype that splits and merges while the genes on it turn over:

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=6, seed=4)
g = genomes.simulate_genomes_nucleotide(
    tree, chromosomes=3, root_length=4000, genes=6, gene_length=200,
    duplication=2.0, duplication_extent=900, loss=2.0, loss_extent=900,
    transfer=1.0, transfer_extent=900, transfer_to="distance",
    fission=0.2, fusion=0.2, chromosome_origination=0.05, chromosome_loss=0.05,
    inversion=1.0, seed=5)
```

A family lost in **every** extant lineage still gets a complete tree, since it is still history, but `.gene_trees[fam].extant` is `None` and no `_extant.nwk` is written. Only a family deleted from every node has no tree at all, and it still appears in `.gene_spans` because it was declared. Ask `.gene_trees` what it holds rather than assuming a declared family is in it.

## On the command line

`--resolution nucleotide` takes the same event rates as Chapter 5 and adds two things: how to set up the initial genome, and how long an event is in base pairs.

```bash
# an evenly spaced initial genome: 5 kb, six 300 bp genes, inversions averaging 400 bp
zombi2 genomes out/ --resolution nucleotide \
  --root-length 5000 --genes 6 --gene-length 300 \
  --inversion 1.0 --inversion-extent 400 --duplication 0.3 --loss 0.3 --seed 1

# or start from a real genome: the GFF declares the replicons and the genes,
# and a paired FASTA supplies the actual DNA those coordinates hold
zombi2 genomes out/ --resolution nucleotide \
  --gff ecoli.gff --fasta ecoli.fasta \
  --inversion 0.5 --loss 0.4 --loss-extent 900 --seed 1
```

Every event kind has its own `--<event>-extent`, the mean of a geometric draw in base pairs: `--inversion-extent`, `--loss-extent`, `--duplication-extent`, `--transfer-extent`, `--transposition-extent`, `--translocation-extent`, `--origination-extent`.

## Outputs

| File | What it holds |
|---|---|
| `genome_events.tsv` | the gene genealogy, in the format every resolution writes |
| `block_events.tsv` | this resolution's own log: the copy lineages over ancestral intervals |
| `rearrangement_events.tsv` | the inversions, transpositions and translocations, in time order |
| `blocks.tsv` | every node's genome as its block mosaic |
| `genes.tsv` | where each declared gene sits in initial coordinates, and on which strand |
| `initial_genome.tsv` | the genome the run started with |
| `initial_sequence.fasta` | the initial DNA, when a `--fasta` supplied it |
| `chromosome_events.tsv` | the chromosome network, one row per edge |
| `gene_trees/` | one Newick per family, complete and extant |
| `gff/`, `bed/` | `genome_<lineage>.gff` (the genes) and `genome_<lineage>.bed` (the blocks), one file per node |

Everything is written by default. Appendix B gives the columns and the formats.
