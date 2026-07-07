# Ordered genomes

By default a genome is **order-free** (`UnorderedGenome`): a multiset of gene families with
copy numbers — all you need for phylogenetic profiles. When gene *order* matters (synteny,
operons, rearrangements) use **`OrderedGenome`**, the basic ZOMBI-1 model: a circular
chromosome of genes, each carrying a strand orientation, with no intergenic regions.

<figure markdown="span">
![A circular chromosome of oriented genes](../img/genome_circular.svg)
<figcaption>An ordered genome: a circular chromosome of genes, each with a strand
orientation — the substrate for segment duplications, inversions and transpositions.</figcaption>
</figure>

## Using it

`OrderedGenome` takes an `extension` parameter (the segment-length knob), so you pass a
small factory:

```python
from zombi2.genomes import SharedRates, simulate_genomes, OrderedGenome

rates = SharedRates(
    duplication=0.2, transfer=0.1, loss=0.2, origination=0.4,
    inversion=0.3, transposition=0.3,      # rearrangement rates (ordered genomes only)
)
genomes = simulate_genomes(
    tree, rates, initial_families=30, seed=1,
    genome_factory=lambda ids: OrderedGenome(ids, extension=0.5),
)

leaf = next(iter(genomes.leaf_genomes.values()))
leaf.chromosome     # list of OrderedGene(gid, family, orientation=±1), in order
```

## Segment events

Events act on a **contiguous segment** of the circular chromosome. Its length is drawn from
`extension` (a per-step continuation probability): `extension=None` → single genes;
higher values → longer segments.

| Event | Effect on the segment |
|---|---|
| Duplication | tandem copy inserted after the segment |
| Loss | segment removed |
| Transfer | segment copied into a recipient at a chosen position |
| **Inversion** | segment reversed and every strand flipped |
| **Transposition** | segment cut and pasted elsewhere |

Inversions and transpositions change gene order/orientation but **not** gene content, so
they leave the profile matrix and the gene trees unchanged — they show up in the event log
and in the final chromosome order.

## How events reach the genome

Inversion and transposition rates are emitted by `SharedRates` as candidate events. A
genome only undergoes the events it declares in `supported_events()`:

- `UnorderedGenome` supports `{O, D, T, L}` — it silently ignores inversion/transposition
  rates.
- `OrderedGenome` supports `{O, D, T, L, I, P}` — it acts on them.

So the same rate model works for both representations; the genome decides what applies. This
is the [extensibility](extending.md) design in action — adding gene order needed **no**
change to the simulator, sampler, rate interface or output code.
