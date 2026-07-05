# Ordered genomes

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
import zombi2 as z

rates = z.UniformRates(
    duplication=0.2, transfer=0.1, loss=0.2, origination=0.4,
    inversion=0.3, transposition=0.3,      # rearrangement rates (ordered genomes only)
)
genomes = z.simulate_genomes(
    tree, rates, initial_families=30, seed=1,
    genome_factory=lambda ids: z.OrderedGenome(ids, extension=0.5),
)

leaf = next(iter(genomes.leaf_genomes.values()))
leaf.chromosome     # list of OrderedGene(gid, family, orientation=±1), in order
```

Each leaf genome exposes its `chromosome` as an ordered list of `OrderedGene` records, every one
tagged with its family and a strand orientation of $+1$ or $-1$.

## Segment events

Events act on a *contiguous segment* of the circular chromosome rather than a single gene. The
segment length is drawn using the `extension` parameter, a per-step continuation probability:
`extension=None` produces single-gene events, while higher values produce longer segments.

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

![An inversion reverses a chromosomal segment and flips the orientation of each gene in it.](figures/genome_inversion.pdf)

A *transposition* cuts a segment out and pastes it elsewhere on the chromosome.

![A transposition moves a chromosomal segment to a new location.](figures/genome_transposition.pdf)

Inversions and transpositions change gene order and orientation but *not* gene content.

::: note
Because rearrangements preserve gene content, they leave the profile matrix and the gene trees
unchanged. They appear only in the event log and in the final chromosome order.
:::

## How events reach the genome

Inversion and transposition rates are emitted by `UniformRates` as candidate events, but a genome
only undergoes the events it declares in `supported_events()`:

- `UnorderedGenome` supports `{O, D, T, L}`. It silently ignores inversion and transposition rates.
- `OrderedGenome` supports `{O, D, T, L, I, P}`. It acts on them.

The same rate model therefore serves both representations, and the genome decides what applies.
Adding gene order required no change to the simulator, the sampler, the rate interface or the output
code — you swap the `genome_factory` and nothing else.
