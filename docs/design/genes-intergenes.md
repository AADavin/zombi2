# Genes & intergenes — design note

**Status: design to build.** The **genic layer** on the built nucleotide engine
(`zombi2/genomes/nucleotide.py`, merged in #196). Settled with Adrián on 2026-07-21. Read
`genome-api.md` and `SPEC.md §4` first.

---

## What it adds

Today a nucleotide genome is uniform sequence — every block is anonymous ancestry. The genic layer
**declares** some blocks to be **genes** (the units the gene trees are about) and leaves the rest as
**intergenes**, the spacer. Nothing else about the engine changes.

## Vocabulary (extends `SPEC.md §6`)

- **block** — *(unchanged)* a maximal run of one unbroken ancestry; the persistent unit.
- **gene** — a *declared, indivisible* block: one block, one **family**, one id, **never split**. It
  carries a gene tree.
- **intergene** — the free-fragmenting spacer between genes: it splits at any breakpoint into as many
  blocks as events dictate. We **track its block ancestry**, but do **not** build a per-intergene tree.

A genome is an alternating chain `I₀ G₁ I₁ G₂ … Gₙ Iₙ`. With **no genes declared** it is one big
intergene, and the engine is exactly today's uniform-sequence model.

---

## One method: extensions

There is **one** way an event picks its target — mixing two would make the model muddy. Every event
works the same way:

1. start at a breakpoint,
2. extend outward,
3. land when **both ends sit in an intergene** — never inside a gene. If an end would fall inside a
   gene, **redraw** (retry up to `N` times; then the event is a no-op).

So a gene is always **whole**: an event either engulfs it entirely (both breakpoints in the flanking
intergenes) or leaves it untouched. Genes are never split. In code this is just today's engine plus one
rule: **`_split_at` refuses to cut a genic block** — which is exactly what raises the redraw.

**Why redraw rather than snap the breakpoint out to the gene's edge:** length honesty. Snapping would
inflate a drawn extent to "whatever it takes to clear the gene" — with a 30 nt mean extent and a genome
`intergene(10) · gene(1 000 000) · intergene(20)`, a snap turns a 30 into a 1 000 000+ event. Redrawing
keeps the contract exact: **a successful event's span equals its drawn length, or it never happened.**

The extent is in **nucleotides**, but the **rate is per lineage** — *not* per nucleotide. With
extensions that matters: a per-nucleotide rate double-counts size, because a bigger genome would get
proportionally *more* events and each still spans an extent, so the churn per unit time grows with
length. Per lineage says "this lineage does N of these per unit time, wherever they land", which is both
the sane model and what keeps it tractable. Measured on an 8-tip tree, growing the genome 1 kb → 1 Mb:

| scope | 1 kb | 10 kb | 100 kb | 1 Mb |
|---|---|---|---|---|
| **per lineage** | 0.00 s | 0.01 s | 0.03 s | **0.20 s** |
| per nucleotide | 0.04 s | 0.33 s | **41 s** | (hopeless) |

Per lineage keeps the event count flat as the genome grows; per nucleotide explodes it. Block count then
tracks genome *structure* (how many genes were declared) rather than churn — which is exactly what makes
the block representation pay off.

**Which events:** all of them — `inversion`, `transposition`, `translocation`, `loss`, `duplication`,
`transfer`, `origination`, and the chromosome tier (`fission` / `fusion` cut-points are intergenic too).
**No new event names, no scopes, no per-gene variants.** Every event stays plain `(rate, extension)`.

### Consequence — gene turnover is emergent and size-dependent (document this clearly)

A gene changes copy number only when an event **engulfs it whole**, and engulfing a big gene needs a big
extent. So **large genes are rarely lost or duplicated**; small ones more easily. Gene turnover is a
*consequence* of event sizes, not a rate you set directly.

This is honest — a big gene really does need a big deletion — but it **must be documented prominently**,
because it surprises people: setting a small deletion rate on a genome of 1 000 nt genes leaves the genes
essentially untouched. The guide should say plainly: *to get more gene turnover, use larger events.*

---

## Seeding

The root genome is **declared**, in one of two forms:

- a **GFF** — exact gene coordinates (and names) on the seed replicon(s); or
- a **genes + intergene-distribution file** — a gene count (per-gene lengths, or a length distribution)
  plus an intergene-length distribution, laid down as the alternating chain.

Each seeded gene becomes a **family** with an id — the handle for its gene tree.

---

## Genealogy & recovery — almost nothing changes

A gene is a block that **never splits**, so:

- **each gene → one persistent block-lineage → one gene tree.** This falls out of the existing recovery
  for free (a gene is simply a root-block that is never cut).
- **intergenes fragment** into blocks whose ancestry the copy-lineage log already records. We keep that
  genealogy, but do **not** surface per-intergene trees.

So the recovery is the recovery we built. The only additions:

1. a **genic tag** on blocks (a gene id / family), set at seeding and inherited through every event;
2. **`_split_at` refuses to cut a genic block** — this raises the redraw and keeps genes atomic (an
   engulfed gene is wholly in or out, so the leaves-==-observed-copies cross-check still holds);
3. **`result.gene_trees` filters to genic families** — one tree per gene.

---

## Deferred

- **The per-copy dial** — an opt-in second mechanism that picks a *gene* directly (uniform over genes),
  making turnover size-blind and settable per gene. We start with extensions **only**. It is genuinely a
  *second selection* — size-blind control means picking genes, not nucleotides — and mixing two selection
  methods is what muddies the model, so it stays out until it earns its place.
- **Pseudogenization** — an explicit `gene → intergene` event.
- **Family heterogeneity** (`ByFamily` / `family_speed`).
- **Replacement transfer** — additive only for now.
- **Substitution-model routing** — the gene/intergene tag will select coding vs neutral models at the
  sequence level; out of scope here.

---

## Open — to settle while building

1. **`origination`** — is it specifically a new *gene* (+ flanking intergene), or new sequence that may
   be genic or not?
2. **Retry cap `N`** — default, and whether it is user-exposed.
3. **GFF / declaration-file schema** — the column set, and where multi-replicon karyotypes are expressed.
