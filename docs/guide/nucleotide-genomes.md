# Nucleotide genomes

The standard gene-family model ([gene families & rates](gene-families.md)) treats each gene
as an indivisible token. The **nucleotide genome** model works one level down: a genome is a
sequence of individual nucleotides, and structural events act on **variable-length segments**
of them. This resolves paralogy, xenology, and gene order/orientation at nucleotide
resolution, and reconstructs a gene tree for every stretch of shared ancestry.

## The model

`simulate_nucleotide_genomes` evolves a genome forward along a fixed species tree. It starts
from `initial_size` gene(s) of `root_length` nucleotides at the root, and these events fire:

| Event | Effect |
|---|---|
| `duplication` | copy a segment elsewhere (tandem / paralog) |
| `transfer` | copy a segment into another lineage (xenolog) |
| `loss` | delete a segment |
| `inversion` | reverse a segment's orientation |
| `transposition` | move a segment |
| `origination` | insert a brand-new gene under a fresh source |

Duplication/transfer/loss/inversion/transposition are **per-nucleotide** rates — the total
genome rate is `rate × current_length`, so longer genomes evolve faster — while `origination`
is **per branch**. Event lengths follow a geometric model with mean `1/(1 − extension)`
nucleotides (`extension=0.99` → ~100 nt).

```python
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

result = z.simulate_nucleotide_genomes(
    tree, root_length=1000,
    duplication=1e-4, transfer=5e-5, loss=1.5e-4,
    inversion=1e-3, transposition=5e-5, origination=0.2, seed=1)
```

!!! warning "Keep gain ≤ loss"
    Duplication and additive transfer grow the genome without a cap. Over long ages keep them
    at or below `loss` to avoid runaway growth.

<figure markdown="span">
![Structural events acting on a nucleotide genome](../img/nucleotide_events.svg)
<figcaption>Structural events on a nucleotide genome: an inversion reverses a segment's
orientation (the reversed colour gradient), a tandem duplication lengthens it — each acting on
a variable-length stretch of nucleotides.</figcaption>
</figure>

## Atoms: units of shared ancestry

The simulator partitions the surviving material into **atoms** — maximal segments that share
one unbroken ancestry. Every event boundary splits atoms, so an atom is the finest unit for
which a single gene tree is meaningful. Results are expressed over atoms:

```python
atom_ids, species, matrix = result.profile_matrix()   # copy number of each atom per extant leaf
```

<figure markdown="span">
![Events partition a genome into atoms, each with its own reconstructed tree](../img/nucleotide_tree.svg)
<figcaption>From events to atoms to gene trees: each structural event carves out a segment
(left); the same breakpoints partition every genome into shared <strong>atoms</strong>
(middle); and each atom has its own reconstructed gene tree (right) — a duplication adds a tip,
a loss prunes one, an inversion leaves the genealogy unchanged.</figcaption>
</figure>

## Reading a leaf genome

```python
leaf = tree.leaves()[0]
result.leaf_mosaic(leaf)   # the genome as ordered, signed atoms: [(atom_id, strand), ...]
result.trace_back(leaf)    # every nucleotide's ancestral origin: [(source, src_pos, strand), ...]
```

`leaf_mosaic` gives the leaf as a sequence of atoms with orientation; `trace_back` resolves
each nucleotide to where it came from.

<figure markdown="span">
![A leaf genome traced back to its ancestral coordinates](../img/nucleotide_segments.svg)
<figcaption>Tracing a leaf back to the ancestor: the top bar is the ancestral genome painted by
position; the bottom bar is a leaf, each nucleotide coloured by where it came from — collinear
stretches keep the gradient, an inversion shows it reversed.</figcaption>
</figure>

## Per-atom gene trees & reconciliation

With the default `output="genomes"` (pure-Python engine), the result also carries the full
event log and a reconstructed gene tree per atom:

```python
trees = result.atom_gene_trees()        # {atom_id: (complete_newick, extant_newick)}
result.write_reconciliations("out/")    # reconciled trees + the events table on disk
```

## Genes & intergenes

By default a genome is an unstructured sequence and "genes" are only recovered *post hoc* as
atoms. Pass `gene_intervals` — non-overlapping `(start, end)` (or `(start, end, name)`)
intervals on the root chromosome — to declare **genes** up front. Everything else is
**intergene**. In this *genic mode*:

- **Genes are never split.** Event breakpoints fall only in intergene positions, so every
  event moves, copies, inverts, or deletes a gene *as a whole*. Each gene is therefore exactly
  one atom (one genealogy) wherever it survives; intergene stretches still fragment into many
  intergene atoms. (A short event drawn entirely inside a gene is promoted to the whole gene.)
- **Pseudogenization.** With probability `pseudogenization`, a loss that hits a gene *demotes*
  it to intergene — the sequence is retained, but the gene loses function. It is a state change
  on the continuing lineage (a `G` node in that gene's tree), not a deletion, and it is
  lineage-specific: the gene stays functional in sister lineages.
- **Homologous replacement transfer.** With probability `replacement`, a transfer replaces the
  recipient's syntenic copy instead of adding a new one. The homologous locus is found by the
  genes flanking the transferred segment; the recipient material between those flank genes is
  replaced (and logged as recipient losses). When the recipient has no such homolog, the
  transfer falls back to additive insertion.
- **Origination** mints a brand-new gene (its own gene tree), as in the base model.

```python
genes = [(100, 180, "dnaA"), (300, 360, "gyrB"), (500, 620, "rpoB")]
result = z.simulate_nucleotide_genomes(
    tree, inversion=1e-3, loss=8e-4, duplication=5e-4, transfer=5e-4,
    root_length=1000, extension=0.97, gene_intervals=genes,
    pseudogenization=0.3, replacement=0.4, seed=1)

result.gene_trees()          # {atom_id: (complete, extant)} for the gene atoms
result.intergene_trees()     # …and for the intergene atoms
result.pseudogenizations()   # [(atom_id, gene_id, species_branch, time, gene_lineage), …]
```

Atoms carry their classification (`atom.kind` is `"gene"`/`"intergene"`, `atom.gene_id`), so
`gene_atoms()` / `intergene_atoms()` partition the atom set. Genic mode runs on the Python
engine only (the Rust `profiles` path does not model genes). On the CLI:

```bash
zombi2 genomes -t species_tree.nwk --rate-model nucleotide \
  --genes genes.tsv --pseudogenization 0.3 --replacement 0.4 \
  --inversion 0.001 --loss 0.0008 --output profiles trees -o out/
```

where `genes.tsv` is a BED/TSV of `start end [name]` lines. The run writes `genes.tsv` (the
annotation, including originated genes), gene/intergene trees under `Gene_trees/` and
`Intergene_trees/`, a `kind`/`gene_id` column in `atoms.tsv`, and `Pseudogenizations.tsv`.

### Starting from a real genome (GFF)

Instead of writing intervals by hand, point the model at a real annotation — e.g. a RefSeq
bacterial chromosome — and it copies the genome's **length** and **gene coordinates** (the
intergenes are the gaps). `read_gff` returns both; because bacterial genes sometimes overlap
(shared start/stop codons, nested ORFs) and the genic model forbids breakpoints inside a gene,
overlaps are removed by trimming — each gene's start is clipped to the previous gene's end, and a
gene swallowed whole is dropped:

```python
g = z.read_gff("GCF_000005845.2_ASM584v2_genomic.gff")   # E. coli K-12 MG1655
g.length, len(g.genes), g.n_trimmed, g.n_dropped         # 4641652, 4480, 768, 26

result = z.simulate_nucleotide_genomes(
    tree, root_length=g.length, gene_intervals=g.genes,
    inversion=2e-6, loss=1.5e-6, extension=0.999, pseudogenization=0.3, seed=1)
```

On the CLI, `--gff` sets the length and genes in one step (superseding `--genes`/`--root-length`):

```bash
zombi2 genomes -t species_tree.nwk --rate-model nucleotide \
  --gff ecoli.gff --inversion 2e-6 --loss 1.5e-6 --pseudogenization 0.3 \
  --output profiles trees -o out/
```

The GFF may be gzipped. For a multi-sequence file (chromosome + plasmids), the most-annotated
sequence is used by default; `--gff-seqid ID` (or `read_gff(..., seqid=...)`) picks another. The
genes keep their annotation names (locus tag / `Name`), so `genes.tsv` and the trees are labelled
with real gene ids.

## The Rust fast path

`output="profiles"` runs the compiled `zombi2_core` Rust engine over leaf segments only —
much faster, and enough for `profile_matrix()`, `leaf_mosaic()`, and `trace_back()`. It emits
**no event log**, so `atom_gene_trees()` / `atom_histories()` are unavailable, and it
**requires** the built extension (see [the Rust engine](rust-engine.md)):

```python
result = z.simulate_nucleotide_genomes(tree, duplication=1e-4, loss=1.5e-4,
                                       inversion=1e-3, seed=1, output="profiles")
```
