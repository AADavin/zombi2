# Gene-order export

`zombi2 tools export` turns a **nucleotide** genome simulation into the file formats that
gene-order / synteny studies consume — the analysis complement of the simulator's own outputs. It
reads a `zombi2 genomes --genome-resolution nucleotide` output directory and derives each format from
the reconstructed per-node gene orders, so it is a pure a-posteriori step: simulate once, export as
many views as you need.

It is the ZOMBI2 counterpart of the `zombiExporter` utility in Krister Swenson's ZOMBI fork.

## Inputs

The export reads what the genomes run already wrote — no re-simulation:

- **`bed/<node>.bed`** — each node's genes in genome order, with orientation (written by
  `--write bed`, which needs a gene annotation via `--genes` / `--gff`). This *is* the
  reconstructed gene order at every node, and it drives the order-based formats.
- **`species_tree.nwk`** — to walk parent→child edges.

So a run that produces everything the export needs looks like:

```bash
zombi2 species -b 1 -d 0.3 --tips 8 --age 3 --seed 1 -o run/S
zombi2 genomes -t run/S/species_tree.nwk --genome-resolution nucleotide \
    --genes genes.tsv --root-length 3000 \
    --inversion 0.01 --transposition 0.005 --write bed geneorder --seed 1 -o run/G
```

where `genes.tsv` is a BED-like `start end [name]` table (0-based, half-open), e.g.

```
0     100   g1
1000  1100  g2
2000  2100  g3
```

## Formats

```bash
zombi2 tools export run/G --format breakpoints gff posortho -o export/
```

Each `--format` writes one file into `-o DIR` (or prints to stdout if `-o` is omitted).

### `breakpoints` — adjacencies broken per tree edge

`breakpoints.tsv`: for every species-tree edge, the gene adjacencies that the rearrangements on
that branch broke. A genome is represented as its set of **circular signed-gene adjacencies** (each
gene has a head `_h` and a tail `_t`; consecutive genes meet at a pair of extremities), and the
broken set on an edge is `adjacencies(parent) − adjacencies(child)`.

```
parent  child  adjacency
root    i1     g1_h|g2_t
root    i1     g2_h|g3_t
i1      n5     g1_h|g3_h
```

An inversion of a single gene breaks exactly two adjacencies (the two that flanked it); a
transposition breaks the adjacencies at the cut and the paste. This is the format for measuring
synteny divergence along a tree.

### `gff` — every node's genes as one GFF3

`genes.gff`: a single GFF3 in which each node is a sequence (`seqid`) and every gene is a `gene`
feature (1-based inclusive coordinates; the gene family is the `Name` attribute, the `ID` is unique
per node so a family may recur after a duplication).

```
##gff-version 3
##sequence-region n1 1 3050
n1  zombi2  gene  116  215  .  -  .  ID=n1.0;Name=g1
n1  zombi2  gene  238  337  .  +  .  ID=n1.1;Name=g2
```

### `posortho` — positional ortholog sets

`positional_orthologs.tsv`: genes across the extant leaves grouped by ancestral family. Because a
gene keeps its identity across the whole tree, every occurrence of one family is an ortholog of the
others; the table lists, per family, its occurrences over the leaves.

```
family  leaf  strand  start
g1      n1    -       620
g1      n2    +       210
g2      n1    +       1240
```

## Exactness & caveats

`breakpoints` and `posortho` are **exact for content-conserving rearrangements** (inversion and
transposition, which permute and reorient genes but neither create nor destroy them). Under
**duplication or loss** the genomes on the two ends of an edge differ in gene *content*, and a gene
family recurs within one genome — so gene names stop being unique, a differing adjacency may reflect
a gained/lost gene rather than a broken one, and a positional-ortholog set becomes family-level
(all copies grouped). Interpret those runs accordingly.

## Planned

- **`dupinfo`** — duplications per gene family. The data lives in the block gene trees
  (`--write trees`), pending a small change so the reconciliation log carries the user gene name.
- **`ffgc`** — extant leaf genomes and sequences in the FFGC (Family-Free Genome Comparison) input
  format.

## Command line

```bash
# all three formats into export/
zombi2 tools export run/G --format breakpoints gff posortho -o export/

# one format to stdout
zombi2 tools export run/G --format breakpoints
```

The Python entry points live in `zombi2.tools.geneorder_export`
(`breakpoints_tsv`, `gff_text`, `posortho_tsv`, and the underlying `read_node_orders` /
`broken_adjacencies` / `positional_orthologs`).
