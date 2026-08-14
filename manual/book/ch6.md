# Genomes III: nucleotide

At the **nucleotide** resolution a chromosome is a **coordinate axis of DNA** rather than a list of gene tokens, and events have an extension measured in base pairs: an inversion reverses 600 bp, a loss deletes 900 bp, a duplication copies 2 kb in tandem. Genes still exist and still get gene trees, but they are stretches of that axis with a start and an end, and the DNA between them is simulated too.

## Genes and intergenes

A genome here is DNA, and its DNA is of two kinds:

- a **gene** is *declared* and *indivisible*: one family, one id, never cut in two. It carries a gene tree.
- an **intergene** is the spacer between genes. Nothing protects it, so events cut it wherever they land.

## Extents

There is one important consequence: **the realised extent is not always the extent you asked for.** It is quantised to the legal breakpoints, and on a gene-dense genome the difference is large. Take thirty-one 3000 bp genes in a 100 kb genome, 93% genic, so the spacers are about 200 bp, and ask for an inversion of:

| asked | realised |
|---|---|
| 500 bp | 59 |
| 1 500 bp | 1 000 |
| 3 000 bp | 2 916 |
| 10 000 bp | 10 315 |

A 500 bp event has nowhere to go but inside a spacer, so it comes out at 59 bp. Long events land near what you asked for, because they can span whole genes.

**The correction can also run the other way.** When no legal end lies within reach of the extent you asked for, the arc snaps out to the nearest legal breakpoint, which is **longer** than the extent you set. Ask for a 5 bp inversion on a genome that is all gene — ten 100 bp genes in 1000 bp — and every event comes out at 100 bp, because the gene joins are the only places a breakpoint may fall. When the genome cannot give what you asked at all, the arc is capped by the replicon and comes out shorter. Either way the event still happens: an extent is a request, and the genome answers it.

The one case that yields no event is degenerate: a replicon with no legal end within reach at all, such as one under 2 bp, where the event is skipped rather than forced.

## A note on rates

**Every rate here is per lineage.** The rate sets how often a lineage does the event; the extent (above) sets how much DNA it touches. Keeping the rate per lineage means the number you type reads the same whatever the genome's size: a rate counted per base pair would rise as the genome grew, so `inversion=5.0` would mean one thing at 10 kb and another at 1 Mb. Per lineage, the event count stays flat as the genome grows. The same tree at `inversion=5.0`, with the genome ten times longer each row, gives:

```
   10 000 bp  ->  77 inversions
  100 000 bp  ->  77 inversions
1 000 000 bp  ->  77 inversions
```

The chromosome rates below are the exception: `fission`, `fusion` and `chromosome_loss` are counted **per chromosome**, and `chromosome_origination` per lineage. A fusion joins two chromosomes of the **same topology**, for the same reason it does at the ordered resolution: a ring and a molecule with two ends cannot become one molecule.

Rates here are written the same way as everywhere else — a scope with verbs chained onto it — and the scopes above are the defaults, so a bare number stays a bare number. The **skyline** works: `inversion = PerLineage(5.0).changing_at({0: 1.0, 3: 0.2})` drops the inversion rate fivefold at time 3, and the run re-reads its rates at each step rather than racing past it.

So does **conditioning**. Every rate here takes a `scaled_by`, so a trait can drive how much DNA a lineage sheds, which is genome reduction as it is usually meant, and can drive the rearrangements too:

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=56)
g = genomes.simulate_genomes_nucleotide(
    tree, root_length=3000, genes=3, gene_length=400,
    inversion=1.0, inversion_extent=600, seed=56)
```

```python
from zombi2 import traits
from zombi2.params import Extent, PerLineage

habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=2)
loss = PerLineage(0.8).scaled_by(habitat, {"host": 20.0, "free": 0.5})
```

**The extent takes the same verbs**, and that is a different statement:

```python
loss        = PerLineage(0.8).scaled_by(habitat, {"host": 20.0, "free": 0.5})  # more often
loss_extent = Extent(150).scaled_by(habitat, {"host": 6.0,  "free": 1.0})      # bigger chunks
```

The first raises how often a host-restricted lineage deletes, the second how much each deletion takes. Set both and they multiply: the DNA shed per unit time goes up by the product, not the sum.

A modifier on an *extent* is read when an event occurs, so unlike the same modifier on a rate it adds no step to the run's clock. Chapter 9 covers what a driver is and how to grow one; anything the level does not accept raises rather than being quietly ignored.

## Who receives a transfer

`transfer_to` is Chapter 4's recipient rule, and it works here unchanged: `"uniform"`, `"distance"` / `Distance(decay=)`, a `Clades(...)` kernel over named clades, or `Recipients().weighted_by(...)` read off a trait. It is not a rate: the numbers are weights normalised over the lineages alive when a transfer occurs, so it says who receives and never how much transfer happens.

```python
tree6 = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=6)
flows = genomes.Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
g6 = genomes.simulate_genomes_nucleotide(
    tree6, root_length=6000, genes=6, gene_length=400,
    transfer=2.0, transfer_extent=900, seed=6,
    transfer_to=genomes.Clades({"A": ["n49", "n50"], "B": ["n30", "n36"]}, flows))
```

A transfer here is always **additive**, since the donor keeps its copy, so steering changes only which lineage the arc lands on. Nothing is taken from anyone, and a transfer whose every candidate weighs 0 simply does not fire.

## Genes are never split

![Why a gene is never split. The strip beneath the replicon is the set of positions a breakpoint may take: every stretch of spacer, and no part of a gene. An event's two ends are drawn from that strip in the first place, rather than drawn anywhere and repaired afterwards, so the deletion above either takes gene 2 whole or does not touch it. A genome with no spacer at all is still legal: the strip is then just the joins between genes, and events move genes around whole.](figures/legal_breakpoints_print.png){width=100%}

An event either **engulfs a gene whole** or leaves it alone; a breakpoint never falls strictly inside one. So an event does not pick an arc and then clean up afterwards. Both of its ends are drawn **directly from the positions where a breakpoint is legal**. A genome can therefore be **all gene, with no spacer at all**: ten 100 bp genes in 1000 bp is a legal genome, and it evolves. Its breakpoints simply all fall at the joins between genes, so genes are inverted, moved, duplicated and lost whole. Genes may sit flush; they are not required to leave a gap.

Two leaves of the run at the top of the chapter show what those events did. Each line is a **block**: a stretch of DNA with one unbroken ancestry, written as the interval it came from on the initial sequence, read forward (`+`) or reverse-complemented (`−`). A gene is always one block, because nothing may cut it. Spacer is not, so a run of several intergene lines in a row is simply spacer that has been cut apart and rearranged: the accumulated breakpoints of everything that happened to that lineage.

```python
def show(node):
    for chromosome, blocks in g.mosaic(node).items():
        print(f"n{node}, chromosome {chromosome}")
        for source, start, end, strand in blocks:
            gene = next((f"gene {f}" for f, span in g.gene_spans.items()
                         if span == (source, start, end)), "intergene")
            print(f"  [{start:4d},{end:4d}) {'+' if strand == 1 else '−'} "
                  f"{end - start:4d} bp  {gene}")

show(2)
show(5)
```

```
n2, chromosome 2
  [   0, 600) +  600 bp  intergene
  [ 600,1000) +  400 bp  gene 1
  [1000,1140) +  140 bp  intergene
  [2000,2065) −   65 bp  intergene
  [1600,2000) −  400 bp  gene 2
  [1573,1600) −   27 bp  intergene
  [1410,1573) +  163 bp  intergene
  [1140,1410) −  270 bp  intergene
  [2065,2563) +  498 bp  intergene
  [2563,2596) −   33 bp  intergene
  [2596,2600) +    4 bp  intergene
  [2600,3000) +  400 bp  gene 3
n5, chromosome 5
  [   0, 512) −  512 bp  intergene
  [2600,3000) −  400 bp  gene 3
  [2532,2600) −   68 bp  intergene
  [ 512, 600) +   88 bp  intergene
  [ 600,1000) +  400 bp  gene 1
  [1000,1600) +  600 bp  intergene
  [1600,2000) +  400 bp  gene 2
  [2000,2532) +  532 bp  intergene
```

Both leaves still carry all three genes, and both are still 3000 bp. In `n2` an inversion covered gene 2, which now reads on the `−` strand with the spacer around it reversed; in `n5` a different inversion covered gene 3 instead, taking the spacer on either side of it with it. The coordinates are what make this readable: every block still names where it came from in the root, so `[600, 1000)` is gene 1 wherever it turns up and whichever way it points.

## The initial genome

The **initial genome**, the genome the run starts from at time 0 before any event, is declared in one of two ways. (It is not the same as the genome at the root *node*: the root branch is real simulated time, so that one already carries stem events. See `.initial_genome` below.)

**Evenly spaced genes.** `genes=N, gene_length=L` lays down `N` genes of `L` bp on each replicon, spreading the leftover DNA evenly between them. Good for controlled experiments, since gene density is then a number you set.

**A GFF file.** `gff="genome.gff"` takes exact coordinates from a real annotation. `##sequence-region` gives each replicon's extent, `gene` features give coordinates, strand and name, and other feature types are ignored. Names land in `result.gene_names`, so you can follow a named gene through the run. `gff=` and `genes=` are mutually exclusive; a GFF already declares the genes. Genes may touch but never overlap, since a gene is one indivisible block. Real annotations do overlap, usually by a base or two where genes abut in an operon, so an overlap is refused rather than guessed at: `trim_overlaps=True` (`--trim-overlaps`) pushes each overlapping gene's start to its neighbour's end instead, and drops any gene swallowed whole.

A GFF gives coordinates, not letters. `fasta="genome.fasta"` supplies the DNA those coordinates hold — one record per replicon, matched by id, each exactly its declared length — and a later `zombi2 sequences` run founds its blocks from that DNA (Chapter 7).

## What a run gives back

`simulate_genomes_nucleotide` returns a `NucleotideGenomesResult`. `.genomes` and `.node_genomes`
carry what they carry at every resolution — the observed tips, and every node — each genome now a
list of `Chromosome`s of `Block`s. What is particular to this resolution is that a genome is
*rebuilt* rather than stored: `.root_blocks` is the recovered root partition and `.assembly(node)`
says how one node's genome is made from it. Appendix B lists all of it.


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
