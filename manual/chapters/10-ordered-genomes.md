# Ordered genomes

Chapter 7 introduced three models of genome evolution. The first — unordered gene families
(Chapter 8) — treats a genome as an unordered *set*. This chapter covers the second, **ordered gene
families**: genes sit on a chromosome and their order matters, but distance is still counted in
genes, not nucleotides. (The third model, nucleotide genomes, is Chapter 11.)

By default a genome in ZOMBI2 is *order-free* (`UnorderedGenome`): a multiset of gene families with
copy numbers, which is all you need for phylogenetic profiles. When gene *order* matters — synteny,
operons, rearrangements — use `OrderedGenome`, the classic ZOMBI-1 model. An ordered genome is a
circular chromosome of genes, each carrying a strand orientation, with no intergenic regions.

![An ordered genome is a circular chromosome of genes, each with a strand orientation.](figures/genome_circular.pdf)

## Using an ordered genome

You select the representation by passing a `genome_factory` to `simulate_genomes`. Because
`OrderedGenome` needs an `extension` parameter (the segment-length knob), you pass a small factory
that builds one from the gene ids:

```python
from zombi2.genomes import SharedRates, OrderedGenome, simulate_genomes

rates = SharedRates(
    duplication=0.2, transfer=0.1, loss=0.2, origination=0.4,
    inversion=0.3, transposition=0.3,      # rearrangement rates (ordered genomes only)
)
genomes = simulate_genomes(
    tree, rates, initial_families=30, seed=1,
    genome_factory=lambda ids: OrderedGenome(ids, extension=0.5),   # continuation prob; mean 2 genes
)

leaf = next(iter(genomes.leaf_genomes.values()))
leaf.chromosome     # list of OrderedGene(gid, family, orientation=±1), in order
```

Each leaf genome exposes its `chromosome` as an ordered list of `OrderedGene` records, every one
tagged with its family and a strand orientation of $+1$ or $-1$.

### From the command line

The ordered level is also a `--genome-model`, so the same run is available from the `genomes`
command:

```bash
zombi2 genomes -t species_tree.nwk --genome-model ordered \
    --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.4 \
    --inversion 0.3 --transposition 0.2 --mean-length 2 \
    --initial-families 25 --seed 1 --write profiles trees events -o out/
```

Here `--inversion` and `--transposition` are the rearrangement rates (**per gene copy**, unlike the
per-nucleotide rates of the nucleotide model), and `--mean-length` is the segment-length knob: the
**mean number of genes** an event spans (`2` here; omit it for single-gene events), with the lengths
drawn geometrically around that mean. Because
rearrangements live on the shared per-copy rate model, they require the default `--rate-model
shared`; `--rate-model per-genome` runs an ordered genome with duplication, transfer and loss only.

## Segment events

Events act on a *contiguous segment* of the circular chromosome rather than a single gene. The
number of genes in a segment is **geometric**, set by its **mean length** — on the command line the
`--mean-length` knob (in genes), and in the Python API the equivalent per-step continuation
probability `extension` of `OrderedGenome`, the two related by
$\texttt{extension} = 1 - 1/\texttt{mean length}$. Every segment starts with one gene and each
further gene is added with probability `extension`, so a segment of $k$ genes has probability
$(1-\texttt{extension})\,\texttt{extension}^{\,k-1}$ and mean $1/(1-\texttt{extension})$. A mean of
one gene (`extension=None` or 0) gives single-gene events; larger means give longer segments with a
heavier tail.

![The number of genes an event affects is geometric, with the mean set by `--mean-length` (equivalently the continuation probability `extension` $= 1 - 1/\text{mean}$). Small means keep events local (mostly single genes); larger means give longer segments with a heavier tail. The dashed line marks the mean.](figures/segment_length.pdf){width=100%}

| Event | Effect on the segment |
|---|---|
| Duplication | tandem copy inserted after the segment |
| Loss | segment removed |
| Transfer | segment copied into a recipient at a chosen position |
| Inversion | segment reversed and every strand flipped |
| Transposition | segment cut and pasted elsewhere |

The first three events change gene *content* exactly as they do for unordered genomes; the ordered
representation simply also tracks where the affected genes sit and which way they point.

### Rearrangements: inversion and transposition

On top of duplication, transfer and loss, ordered genomes add two rearrangement events. An
*inversion* reverses a segment in place and flips the strand of every gene it contains.

![An inversion reverses a chromosomal segment and flips the orientation of each gene in it.](figures/genome_inversion.pdf){width=72%}

A *transposition* cuts a segment out and pastes it elsewhere on the chromosome.

![A transposition moves a chromosomal segment to a new location.](figures/genome_transposition.pdf){width=72%}

Inversions and transpositions change gene order and orientation but *not* gene content.

::: note
Because rearrangements preserve gene content, they leave the profile matrix and the gene trees
unchanged. They appear only in the event log and in the final chromosome order.
:::

## How events reach the genome

Inversion and transposition rates are emitted by `SharedRates` as candidate events, but a genome
only undergoes the events it declares in `supported_events()`:

- `UnorderedGenome` supports `{O, D, T, L}`. It silently ignores inversion and transposition rates.
- `OrderedGenome` supports `{O, D, T, L, I, P}`. It acts on them.

The same rate model therefore serves both representations, and the genome decides what applies.
Adding gene order required no change to the simulator, the sampler, the rate interface or the output
code — you swap the `genome_factory` and nothing else.
