# Rates: a primer

You have now met rates at four levels — on the [species tree](species-trees.md), in
[genomes](genomes.md), over [traits](traits.md), and along [sequences](sequences.md). They can look
like four unrelated systems. They are not: underneath all of them is a **single grammar**, and this
page is that grammar.

ZOMBI2 simulates **events**: a lineage splits, a gene is duplicated, a nucleotide mutates, a trait
shifts. Every event happens at some **rate**. Learn the grammar once here; it does not change from
the species tree down to the nucleotide.

## What a rate is

A rate is the **expected number of events per unit of time**. Time is measured in whatever your
input tree's branch lengths mean — millions of years, generations, expected substitutions. A loss
rate of `0.5` acting over a branch of length `2` gives, on average, one lost gene.

Events do not arrive on a timetable; they arrive at random, at that average pace (a Poisson
process). Two consequences are worth stating up front, because they trip people up:

- **A rate is not a probability.** It can be larger than 1. "Three events per unit time" is a rate
  of 3.
- **Runs vary.** Doubling a rate doubles the *expectation*, not the count you will see in any single
  simulation.

## The one shape

Every rate in ZOMBI2 has the same shape:

> **rate = base × opportunities × modifiers**

- **base** — how fast, *per opportunity, per unit time*. The number you set.
- **opportunities** — how many independent places the event can happen *right now*.
- **modifiers** — a set of multipliers, each attached to a specific context. Default 1, and they
  multiply together.

That is the whole grammar. Speciation, duplication, substitution, trait change — all of them read
this way. Everything below is just what each piece means at each level.

## Opportunities: "per what?"

The single most clarifying question you can ask about any rate is **"per what?"**

| The event acts… | opportunities = | so a bigger … means more events |
|---|---|---|
| per **lineage** | 1 per living lineage | more lineages |
| per **gene copy** | number of copies | bigger gene family / genome |
| per **nucleotide** | sequence length | longer sequence |

This is where the "per copy vs per lineage" choice actually lives — and "per lineage" is not special
to genomes (it is the same measure speciation uses). **Per copy** means *opportunities = number of
copies*, so a family with 100 copies is duplicated 100× as often as a family with one. **Per lineage**
means *opportunities = 1*, a fixed pace no matter how big the genome grows (one genome per lineage).
Choosing between them is simply choosing what counts as an opportunity.

## Modifiers: context that rescales the base

A modifier multiplies the base for one specific context. For example:

- *this branch* runs 2× faster (a relaxed clock),
- *this family* transfers 5× more readily,
- transfer *from lineage A to lineage B* is favoured (a modifier keyed on a **pair** of lineages),
- *this site* evolves at half speed (among-site rate variation).

Three rules make them predictable:

1. **Default 1** — a context with no modifier behaves normally.
2. **They multiply** — a family modifier of 3 and a branch modifier of 2 combine to 6. Order does
   not matter.
3. **Each keys on some context** — a family, a branch, a pair, a site, a trait value.

This is the flexibility dial. Want family-specific *and* branch-specific sequence evolution? Two
modifiers. Want to boost HGT between two clades? One modifier on transfer, keyed on the lineage
pair. You do not reach for a new model — you attach a modifier.

## Extent: how big, not how often

Some events touch a single thing (one copy lost, one nucleotide changed). Others cover a stretch —
an inversion flips a segment, an indel removes several nucleotides. For those there are **two
separate questions**:

- **how often does it start?** → the *rate* (everything above),
- **how big is it once it starts?** → the *extent*, drawn from its own distribution.

Keep these apart. The rate sets frequency; the extent sets size. A single-copy loss has
*opportunity = per copy* **and** *extent = 1* — the copy is both where it happens and what it hits.
An inversion has *opportunity = per copy* too (each gene is a place it can begin), but *extent > 1*
(it starts at one gene and flips a whole run). So an inversion is **per copy, exactly like
duplication and loss** — the *only* thing that sets it apart is its extent, not its opportunity.
Same event grammar, one extra dial.

!!! warning "The length² trap"
    If *both* the opportunities *and* the extent grow with length, the total material affected per
    unit time grows with length² — more places to start, times a bigger bite each time. Usually you
    want one or the other. Decide on purpose.

!!! note "Edges"
    A spanning event that runs off the end of a chromosome must *wrap*, *truncate*, or be *rejected*.
    That is a modelling choice, stated per model — never left implicit.

## Channels add, modifiers multiply

Two *different processes* that cause the same kind of event **add** their rates. Gene loss, for
instance, can come from pseudogenisation and from deletion; the total loss rate is the sum, and each
channel has its own base, opportunities and modifiers.

A modifier, by contrast, **multiplies** an existing rate. The test is one question: *is this a
separate cause (add), or a rescaling of an existing cause (multiply)?* Keep the two words apart and
the usual "is this rate a base, a modifier, or a sum?" confusion disappears.

## The same shape at every level

| Level | Example event | base | opportunities | typical modifiers | extent |
|---|---|---|---|---|---|
| **Species** | speciation | λ | per lineage | clade shift, diversity-dependence `(1 − N/K)`, trait state | point |
| **Species** | extinction | μ | per lineage | mass-extinction pulse, trait state | point |
| **Genomes** | duplication | d | per copy | family, branch, carrying capacity | point (segment on ordered/nt) |
| **Genomes** | transfer (HGT) | t | per copy | family, donor branch, recipient receptivity, lineage pair | point (segment on ordered/nt) |
| **Genomes** | loss | l | per copy | family, branch | point (segment on ordered/nt) |
| **Genomes** | inversion | i | per copy | family, branch | spans a run of genes / nt |
| **Trait** | discrete change (Mk) | qᵢⱼ | per lineage | branch, hidden state | point |
| **Trait** | continuous drift (BM) | σ² | per lineage · time | branch (relaxed), state | — (see note) |
| **Sequence** | substitution | μ × exchangeabilities | per nucleotide | branch (clock), site (Γ), family, selection (dN/dS) | point |

!!! note "Continuous traits"
    A drifting continuous trait has no countable events and no extent — it is the limit of infinitely
    many infinitesimal steps. But it still fits the skeleton: the **base** is the variance accrued per
    unit time (σ²), the **opportunity** is per lineage, and **modifiers** (a relaxed clock, a
    state-dependent rate) rescale it exactly as elsewhere. Same grammar, no discrete jumps.

## The three genome levels, side by side

The genome level is the one place where the **substrate** changes — a "copy" can be a whole gene or
a single nucleotide — so it is worth seeing the three genome models against the same vocabulary. The
rate grammar is identical across all three; only the opportunity unit and whether events carry extent
differ.

| | **Unordered** (default) | **Ordered** | **Nucleotide** |
|---|---|---|---|
| What it models | a bag of gene families with copy numbers | genes on an ordered, circular chromosome | a nucleotide-resolution genome; genes emerge as blocks |
| The "copy" (substrate unit) | a gene copy | a gene copy at a position | a nucleotide |
| **Opportunity unit** | per gene copy | per gene copy | per nucleotide |
| Positional structure | none (presence/absence + counts) | linear order + strand | full sequence coordinates |
| Core events | duplication, transfer, loss, origination | + inversion, transposition | structural events (inversion, transposition, indels) |
| **Extent** | point (extent = 1) | DTL point; rearrangements span a run of genes | events span a run of nucleotides |
| Extent measured in | — | genes | nucleotides |
| Seeded with | initial gene families | initial gene families | initial chromosomes |

!!! note "Gene conversion"
    Intra-genome gene conversion (one copy of a family overwriting another of the *same* family) is an
    unordered-genome event; it fires only on families holding two or more copies. See
    [Genomes](genomes.md).

Whatever the substrate, reading a genome rate is the same exercise: pick the event, ask *per what*
(a gene copy or a nucleotide), set the base, attach any modifiers, and — if the event spans ground —
give it an extent.

## How to read any rate

Faced with any rate in ZOMBI2, ask five questions in order:

1. **What event?**
2. **Per what?** (the opportunities)
3. **What is the base?** (the number you set)
4. **Any modifiers?** — and what does each key on?
5. **Point or spanning?** — if spanning, what is its extent distribution?

**Worked example — a transfer.** *Event:* horizontal transfer of a gene. *Per what:* per copy, so a
family with 8 copies donates 8× as often as a singleton. *Base:* the transfer rate `t`. *Modifiers:*
this family is a frequent mover (family modifier), the donor branch sits in a high-HGT clade (branch
modifier), and A→B transfer is favoured (pair modifier) — all three multiply onto `t`. *Extent:* on
an ordered genome the transfer can carry a segment, so it draws an extent; on an unordered genome it
is a single gene, extent 1.

## Common pitfalls

- **Treating a rate as a probability.** It is not; it can exceed 1.
- **Forgetting opportunities scale the total.** A modest per-copy rate on a 500-copy family is a
  *large* total rate. "Per copy" hides a multiplication.
- **The length² trap.** Scaling both opportunities and extent with sequence length (see the box
  above).
- **Adding what should multiply, or vice versa.** Separate causes add; rescalings multiply.

---

*See also:* [Species trees](species-trees.md) · [Genomes](genomes.md) · [Traits](traits.md) ·
[Sequences](sequences.md).
