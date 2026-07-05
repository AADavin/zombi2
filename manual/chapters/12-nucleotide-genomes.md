# Nucleotide genomes

Chapter 7 introduced three models of genome evolution. The first two — unordered gene families
(Chapter 8) and ordered gene families (Chapter 11) — treat each gene as an indivisible token. This
chapter covers the third and most detailed: the **nucleotide genome** model, where a genome is a
real sequence of base pairs and structural events act on **variable-length segments** of
nucleotides. This resolves duplication, lateral gene transfer, and gene order and orientation at
nucleotide resolution, and reconstructs a gene tree for every stretch of shared ancestry — the
*blocks* of the sections below.

## The model

`simulate_nucleotide_genomes` evolves a genome forward along a fixed species tree. It starts from
`initial_chromosomes` chromosome(s) of `root_length` nucleotides at the root, and these events fire:

| Event | Effect |
|---|---|
| `duplication` | copy a segment elsewhere (a tandem duplicate) |
| `transfer` | copy a segment into another lineage (a laterally transferred copy) |
| `loss` | delete a segment |
| `inversion` | reverse a segment's orientation |
| `transposition` | move a segment |
| `origination` | insert a brand-new gene under a fresh source |
| `insertion` | add a run of novel intergenic nucleotides (an indel; off by default) |
| `deletion` | remove a run of intergenic nucleotides (an indel; off by default) |

Duplication, transfer, loss, inversion, transposition, and the two indels are **per-nucleotide**
rates — the total genome rate is `rate × current_length`, so longer genomes evolve faster — while
`origination` is **per branch**. Event lengths follow a geometric model with mean `1/(1 - extension)`
nucleotides (`extension=0.99` gives about 100 nt); indels have their own mean-length knob (see
*Intergenic indels* below). The two indels are off by default (rate 0).

```python
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

result = z.simulate_nucleotide_genomes(
    tree, root_length=1000,
    duplication=1e-4, transfer=5e-5, loss=1.5e-4,
    inversion=1e-3, transposition=5e-5, origination=0.2, seed=1)
```

![An inversion reverses a segment's orientation and a tandem duplication lengthens it, each acting on a variable-length stretch of nucleotides.](figures/nucleotide_events.pdf){width=90%}

## Blocks: units of shared ancestry

The simulator partitions the surviving material into **blocks** — maximal segments that share one
unbroken ancestry. Every event boundary splits blocks, so a block is the finest unit for which a
single gene tree is meaningful. Results are expressed over blocks:

```python
block_ids, species, matrix = result.profile_matrix()   # copy number of each block per extant leaf
```

![Each structural event carves out a segment; the same breakpoints partition every genome into shared blocks (left, panels A and B), and each block has its own reconstructed gene tree (right, panel C).](figures/nucleotide_tree.pdf){width=100%}

A duplication adds a tip to the block's tree, a loss prunes one, and an inversion leaves the genealogy
unchanged.

## Reading a leaf genome

```python
leaf = tree.leaves()[0]
result.leaf_mosaic(leaf)   # the genome as ordered, signed blocks: [(block_id, strand), ...]
result.trace_back(leaf)    # every nucleotide's ancestral origin: [(source, src_pos, strand), ...]
```

`leaf_mosaic` gives the leaf as a sequence of blocks with orientation; `trace_back` resolves each
nucleotide to where it came from.

![A leaf genome (bottom) traced back to its ancestral coordinates (top): collinear stretches keep the gradient, an inversion shows it reversed.](figures/nucleotide_segments.pdf)

## Per-block gene trees and reconciliation

With the default `output="genomes"` (the pure-Python engine), the result also carries the full event
log and a reconstructed gene tree per block:

```python
trees = result.block_gene_trees()        # {block_id: (complete_newick, extant_newick)}
result.write_reconciliations("out/")    # reconciled trees + the events table on disk
```

## Genes and intergenes

By default a genome is an unstructured sequence and "genes" are only recovered *post hoc* as blocks.
Pass `gene_intervals` — non-overlapping `(start, end)` (or `(start, end, name)`) intervals on the
root chromosome — to declare **genes** up front. Everything else is **intergene**. In this genic
mode:

- **Genes are never split.** Event breakpoints fall only in intergene positions, so every event
  moves, copies, inverts, or deletes a gene *as a whole*. Each gene is therefore exactly one block
  (one genealogy) wherever it survives; intergene stretches still fragment into many intergene blocks.
  A short event drawn entirely inside a gene is promoted to the whole gene.
- **Pseudogenization.** With probability `pseudogenization`, a loss that hits a gene *demotes* it to
  intergene — the sequence is retained, but the gene loses function. It is a state change on the
  continuing lineage (a `G` node in that gene's tree), not a deletion, and it is lineage-specific: the
  gene stays functional in sister lineages.
- **Homologous replacement transfer.** With probability `replacement`, a transfer replaces the
  recipient's syntenic copy instead of adding a new one. The homologous locus is found by the genes
  flanking the transferred segment; the recipient material between those flank genes is replaced (and
  logged as recipient losses). When the recipient has no such homolog, the transfer falls back to
  additive insertion.
- **Origination** mints a brand-new gene (its own gene tree), as in the base model.

```python
genes = [(100, 180, "dnaA"), (300, 360, "gyrB"), (500, 620, "rpoB")]
result = z.simulate_nucleotide_genomes(
    tree, inversion=1e-3, loss=8e-4, duplication=5e-4, transfer=5e-4,
    root_length=1000, extension=0.97, gene_intervals=genes,
    pseudogenization=0.3, replacement=0.4, seed=1)

result.gene_trees()          # {block_id: (complete, extant)} for the gene blocks
result.intergene_trees()     # ...and for the intergene blocks
result.pseudogenizations()   # [(block_id, gene_id, species_branch, time, gene_lineage), ...]
```

Blocks carry their classification (`block.kind` is `"gene"` or `"intergene"`, plus `block.gene_id`), so
`gene_blocks()` / `intergene_blocks()` partition the block set. Genic mode runs on the Python engine only
(the Rust `profiles` path does not model genes). On the CLI:

```bash
zombi2 genomes -t species_tree.nwk --genome-model nucleotide \
  --genes genes.tsv --pseudogenization 0.3 --replacement 0.4 \
  --inversion 0.001 --loss 0.0008 --write profiles trees -o out/
```

where `genes.tsv` is a BED/TSV of `start end [name]` lines. The run writes `genes.tsv` (the
annotation, including originated genes), gene and intergene trees under `Gene_trees/` and
`Intergene_trees/`, a `kind`/`gene_id` column in `blocks.tsv`, and `Pseudogenizations.tsv`.

### Starting from a real genome (GFF)

Instead of writing intervals by hand, point the model at a real annotation — e.g. a RefSeq bacterial
chromosome — and it copies the genome's **length** and **gene coordinates** (the intergenes are the
gaps). `read_gff` returns both; because bacterial genes sometimes overlap (shared start/stop codons,
nested ORFs) and the genic model forbids breakpoints inside a gene, overlaps are removed by trimming
— each gene's start is clipped to the previous gene's end, and a gene swallowed whole is dropped:

```python
g = z.read_gff("GCF_000005845.2_ASM584v2_genomic.gff")   # E. coli K-12 MG1655
g.length, len(g.genes), g.n_trimmed, g.n_dropped         # 4641652, 4480, 768, 26

result = z.simulate_nucleotide_genomes(
    tree, root_length=g.length, gene_intervals=g.genes,
    inversion=2e-6, loss=1.5e-6, extension=0.999, pseudogenization=0.3, seed=1)
```

On the CLI, `--gff` sets the length and genes in one step (superseding `--genes`/`--root-length`):

```bash
zombi2 genomes -t species_tree.nwk --genome-model nucleotide \
  --gff ecoli.gff --inversion 2e-6 --loss 1.5e-6 --pseudogenization 0.3 \
  --write profiles trees -o out/
```

The GFF may be gzipped. For a multi-sequence file (chromosome plus plasmids), the most-annotated
sequence is used by default; `--gff-seqid ID` (or `read_gff(..., seqid=...)`) picks another. The
genes keep their annotation names (locus tag / `Name`), so `genes.tsv` and the trees are labelled with
real gene ids.

::: tip
Real chromosomes are millions of nucleotides long, and every rate is per-nucleotide. Scale the
per-nucleotide rates down accordingly (the E. coli example uses `inversion=2e-6`, `loss=1.5e-6`) and
raise `extension` toward `0.999` for realistically long segments.
:::

## Intergenic indels

Two optional events change the *length* of the chromosome by editing intergenic DNA:

- **Insertion** lays down a run of novel nucleotides — a fresh source, its own block — inside an
  intergene stretch, lengthening it. It is an intergenic origination of new sequence, with its own
  genealogy rooted at the event.
- **Deletion** removes a run of nucleotides from *within a single* intergene stretch. It is clamped
  so it never reaches into, spans, or removes a gene, and never shrinks the chromosome below a small
  floor (`MIN_GENOME_LENGTH`, one nucleotide) — a deletion can never empty the genome.

Both are **per-nucleotide** rates, off by default (0). Unlike the structural events, an indel run's
length is geometric with mean `indel_mean_length` (default 10) — a knob *independent* of
`extension`, since indels and rearrangements have different length scales. In genic mode indels act
only in intergenes (a gene is never split, spanned, or deleted); with no genes declared they may act
anywhere.

```python
result = z.simulate_nucleotide_genomes(
    tree, root_length=3000, gene_intervals=genes,
    insertion=1e-3, deletion=1e-3, indel_mean_length=12, seed=5)
```

From the command line:

```bash
zombi2 genomes -t species_tree.nwk --genome-model nucleotide --root-length 3000 \
    --insertion 0.001 --deletion 0.001 --indel-mean-length 12 \
    --seed 5 --write profiles trees -o out/
```

Because an insertion roots a new block and a deletion terminates one, indels enter the block
genealogy exactly like an origination and a loss respectively. They require the Python engine
(`output="genomes"`); setting either rate routes the run there automatically.

## The Rust fast path

`output="profiles"` runs the compiled `zombi2_core` Rust engine over leaf segments only — much
faster, and enough for `profile_matrix()`, `leaf_mosaic()`, and `trace_back()`. It emits **no event
log**, so `block_gene_trees()` is unavailable, it does **not** model genes/intergenes or indels, and
it **requires** the built extension:

```python
result = z.simulate_nucleotide_genomes(tree, duplication=1e-4, loss=1.5e-4,
                                       inversion=1e-3, seed=1, output="profiles")
```
