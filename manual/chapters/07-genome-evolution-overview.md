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
one beneath it. On the command line, `--genome-model` selects the level: `unordered` (the default),
`ordered`, or `nucleotide`.

The second choice, orthogonal to the first and meaningful only within the unordered level, is how
the four gene-family rates *vary across families* — the **rate model** (`--rate-model`; see *Rate
models* below).

| Level (`--genome-model`) | A genome is… | Chapter |
|---|---|---|
| **unordered** | a *set* of gene families (presence/absence) | Chapter 8 |
| **ordered** | genes on a chromosome (order matters, length does not) | Chapter 11 |
| **nucleotide** | a real sequence (genes + intergenes at true coordinates) | Chapter 12 |

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

## Rate models

Within the unordered level, a second axis controls **how the four rates vary across gene families**.
This is the `--rate-model` flag (a `RateModel` object in Python). Two are available on the command
line — `shared` (the default) and `per-genome`; per-family rates are reached through the Python API.

| `--rate-model` | The rate is… | Family size grows | Object |
|---|---|---|---|
| **shared** (default) | one per-copy rate, the same for every family | exponentially | `SharedRates` |
| **per-genome** | one per-genome rate, size-independent | linearly | `PerGenomeRates` |
| *per-family* | each family draws its own rates | exponentially, per family | `FamilySampledRates` |

The level and the rate model are independent: any rate model can drive any unordered run, and the
choice never affects *which* events fire, only *how often* each family fires them. Chapter 8 works
through the rate models in detail. Gene-family *coupling* — where families gain and lose
non-independently (a Potts/Ising panel) — is not part of this release's CLI or manual; it is
available through the Python API (`zombi2.simulate_coupled`).
