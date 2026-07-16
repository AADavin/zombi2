```{=latex}
\part{Genomes}
```

# Genome evolution

Once ZOMBI2 has a species tree, it populates the branches with **genomes**. ZOMBI2 offers three
models of genome evolution, differing in how much structure a genome carries — from an unordered
*set* of gene families to a real nucleotide sequence.

- **Unordered gene families** — a genome is an unordered *set* of gene copies. Families have no
  position, no neighbours, no length; each evolves on its own. This is the fastest model, and the
  right one when what you care about is gene content and copy number: presence/absence profiles,
  reconciliations, family sizes.
- **Ordered gene families** — genes sit on an *ordered chromosome*. Position now matters:
  rearrangements shuffle neighbourhoods, and a transfer can land beside a syntenic locus. But
  distance is counted in *genes*, not nucleotides — length stays abstract.
- **Nucleotide genomes** — the chromosome is a real sequence of base pairs. Events have real
  lengths, genes and intergenes have real coordinates, and structural events act at nucleotide
  resolution.

![The three models of genome evolution, in increasing structure: an unordered set of gene families; genes on an ordered chromosome (order matters, length does not); and a real nucleotide sequence with genes and intergenes at true coordinates.](figures/genome_models.pdf){width=100%}

This chapter is a map of that space. The chapters that follow work through each level in turn; here
we lay out the two choices that define any genome run, and — most importantly — which events act at
which level.

## Two orthogonal choices

Two independent choices define a genome run.

The first is the **level** — unordered, ordered, or nucleotide. The level is what a genome *is*, and
therefore which events can act on it (see *Events, by level* below); each level *adds* events to the
one beneath it. On the command line, `--genome-resolution` selects the level: `unordered` (the default),
`ordered`, or `nucleotide`.

The second choice, orthogonal to the first, is the **opportunity** — how each rate is counted: **per
gene copy** (the default, so a bigger genome has proportionally more events) or **per lineage** (one
chance for the whole genome, size-independent), set by `--rate-per`. The [tour](#how-rates-work-how-many-clocks-how-fast) covers
the opportunity axis — and the per-family / per-lineage multipliers you can layer on top — in full;
here it matters only as the second of the two knobs.

| Level (`--genome-resolution`) | A genome is… | See |
|---|---|---|
| **unordered** | a *set* of gene families (presence/absence) | [Unordered genomes](#unordered-genomes) |
| **ordered** | genes on a chromosome (order matters, length does not) | [Ordered genomes](#ordered-genomes) |
| **nucleotide** | a real sequence (genes + intergenes at true coordinates) | [Nucleotide genomes](#nucleotide-genomes) |

## Events, by level

A genome model is defined as much by *what can happen to a genome* as by what a genome is. The three
levels form a hierarchy: each one inherits every event of the level below and adds its own, acting on
the extra structure that level introduces.

![The event sets nest — unordered $\subset$ ordered $\subset$ nucleotide. The four gene-family events act on any genome; a richer level only ever *adds* events, never removes them: gene order brings the two rearrangements, and real sequence brings the base-pair events.](figures/event_levels.pdf){width=80%}

- **Unordered** — the four gene-family events. **Origination (O)** seeds a brand-new family;
  **duplication (D)** copies a gene within a genome; **transfer (T)** copies a gene into another
  lineage alive at the same time; **loss (L)** removes a copy. These four act on gene *content*
  alone — which families are present, and in how many copies.
- **Ordered** — the four events above (now acting on contiguous *segments* of the chromosome) plus
  two **rearrangements**: **inversion** reverses a segment and flips its strands, and
  **transposition** cuts a segment and pastes it elsewhere. These change gene *order*, not content.
- **Nucleotide** — everything above, now at base-pair resolution, plus the events that exist only
  once a genome has sequence: **intergenic insertion and deletion** (indels that lengthen or shorten
  intergenic DNA), **pseudogenization** (a loss that demotes a gene to intergenic DNA rather than
  deleting it), **homologous replacement** (a transfer that overwrites the recipient's syntenic
  copy), and **substitution** (point mutation of the sequence itself).

| Level | Events it adds | These events act on… |
|---|---|---|
| **Unordered** | origination, duplication, transfer, loss | gene *content* (which families, how many copies) |
| **Ordered** | inversion, transposition | gene *order* (segments of the chromosome) |
| **Nucleotide** | insertion, deletion, pseudogenization, homologous replacement, substitution | the *sequence* (base pairs; gene/intergene coordinates) |

Each level is a strict superset: choosing a richer level never removes an event, it only adds the
ones its extra structure makes possible. A nucleotide genome still undergoes duplication, transfer,
loss and rearrangement — but now they have real lengths and can split or fuse at any base.

## Rates at the genome level

The opportunity axis of the [tour](#how-rates-work-how-many-clocks-how-fast) applies here directly,
as the `--rate-per` flag (`per=` in Python): **copy** (the default — `Rates(...)`, the total rate
scales with genome size, so families grow *exponentially*), **lineage** (`Rates(..., per="lineage")` —
one constant rate for the whole genome, size-independent, so content grows *linearly*), or **shared**
(`Rates(..., per="shared")` — a single clock for the whole family). On top of any of these,
**multipliers** make rates differ by context — per family (`FamilySampledRates`, or a `FamilyModifier`
overlay) or per lineage (a `LineageModifier`, the relaxed-clock analogue) — changing *how fast*, not
*how many chances*.

The level and the opportunity are independent: any rate model can drive any level, and the choice
never affects *which* events fire, only *how often* each family fires them.
[Unordered genomes](#unordered-genomes) works through the genome rate models in detail.
