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

Only the **base** carries units — it is a *rate*, `time⁻¹`. **Opportunities** (a count) and
**modifiers** (a multiplier) are both **dimensionless**: they say *how many* places and *how much
faster*, never touching the units. The product is again a `time⁻¹` rate — the actual propensity.
So every rate has two questions: **what are its units** (it is `time⁻¹`, always) and **per what**
(the opportunity). The first is fixed; the second is a modelling choice, and it is the interesting one.

That is the whole grammar. Speciation, duplication, substitution, trait change — all of them read
this way. Everything below is just what each piece means at each level.

## Opportunities: "per what?"

The single most clarifying question you can ask about any rate is **"per what?"** Every rate is per
unit time; the choice is whether it is *also* per some object, and which. The opportunity is a
**dimensionless count read off the current state** — it turns the base into the propensity:
`propensity = base × opportunities`.

| The base fires… | opportunities (the count) = | total propensity |
|---|---|---|
| **globally** (per the whole process) | 1 | just the base — constant |
| per **lineage** | number of living lineages | base × N |
| per **gene copy** | number of copies in the genome | base × copies |
| per **nucleotide** | sequence length | base × length |

**The count decides the dynamics** — specifically, whether it *tracks the quantity that is growing*:

- If it **tracks** (the count grows with the thing) → **multiplicative → exponential**. A per-copy
  duplication rate: more copies → more duplications → more copies. A per-lineage speciation rate:
  more lineages → more speciations → more lineages.
- If it is **fixed** (the count does not grow) → **additive → linear**. A per-genome (not per-copy)
  loss rate stays put as copies pile up; a global speciation rate stays put as lineages pile up.

This is the one place "per copy vs per lineage" actually lives, and it is the same axis up on the
species tree and down in the genome — but **mind the word "lineage"**, because it flips:

- **Species, per lineage** — each of `N` lineages speciates on its own. Opportunity = `N`, which *is*
  the growing quantity, so it tracks → the tree grows **exponentially**. (The standard birth–death.)
- **Genome, per copy** — each copy duplicates on its own. Opportunity = copies, the growing quantity
  → families grow **exponentially**.
- **Genome, per lineage** (`Rates(per="lineage")`) — here the *genome*, not the copy, is the unit;
  opportunity = 1 per family, fixed as copies pile up → families grow **linearly**.

So "per lineage" is exponential for speciation but linear for a gene family — *same words, opposite
scaling* — because in one case the lineages are what grows and in the other the copies are. The rule
is always the same: **does the opportunity count track the growing quantity?**

**And "per what" is really "independent, or shared?"** A per-lineage rate says the lineages act
**independently** — each one behaves as if the others were not there (the honest null, and the
default). A global rate says they **share one budget**: one event every `1/base` on average,
*whoever* it lands on, regardless of how many lineages exist — a strong coupling. That is the deep
reason per-lineage is standard.

!!! example "Two birth–death trees, same units"
    Speciate at rate `λ` **per lineage**: propensity `= λN`, so `E[N(t)] = N₀·e^{(λ−μ)t}` —
    exponential. Speciate at a constant **global** rate `Λ` (pick the lineage uniformly): propensity
    `= Λ`, so `E[N(t)] = N₀ + (Λ−M)t` — linear. Both `λ` and `Λ` are `time⁻¹`; the *only* difference
    is the opportunity count (`N` vs `1`). Equivalently, a global rate `Λ` shared among `N` lineages
    is a per-lineage rate `λ(N) = Λ/N` — diversity-dependence in disguise. The opportunity is a knob:
    `BirthDeath(per="lineage")` (the default) is exponential, `BirthDeath(per="shared")` (`--per shared`)
    is linear — the same model, one clock per lineage versus one for the whole tree. See
    [Species trees](species-trees.md).

## Modifiers: context that rescales the base

A modifier multiplies the base for one specific context. For example:

- *this lineage* runs 2× faster (a relaxed clock),
- *this family* transfers 5× more readily,
- transfer *from lineage A to lineage B* is favoured (a modifier keyed on a **pair** of lineages),
- *this site* evolves at half its base rate (among-site rate variation).

Three rules make them predictable:

1. **Default 1** — a context with no modifier behaves normally.
2. **They multiply** — a family modifier of 3 and a lineage modifier of 2 combine to 6. Order does
   not matter.
3. **Each keys on some context** — a family, a lineage, a pair, a site, a trait value.

This is the flexibility dial. Want family-specific *and* lineage-specific sequence evolution? Two
modifiers. Want to boost HGT between two clades? One modifier on transfer, keyed on the lineage
pair. You do not reach for a new model — you attach a modifier.

## Per-family rates: the base or a modifier (never an opportunity)

A common question: *how do I give each gene family its own rate?* Two ways — and the choice is exactly
the rate-vs-modifier distinction:

- **the base** — `FamilySampledRates` gives each family its own `(dup, transfer, loss)` **rate**, either
  drawn from distributions or fixed by name (`rates={"A": (0.8, 0.1, 0.2)}`, CLI `--family-rates FILE`).
  The rate *itself* varies per family.
- **a modifier** — `FamilyModifier` keeps one shared base rate and multiplies it by a per-family
  **factor** (`factors={"A": 1.6}`, or drawn `per_family=…`). A dimensionless multiplier varies.

They are interchangeable — `rate_f = base × factor_f`. Read the duplication rate off one genome where
family **A** holds 3 copies and family **B** holds 1:

| model | family A | family B | what varies |
|---|---|---|---|
| `Rates(dup=0.5)` | 0.5 × 3 = 1.5 | 0.5 × 1 = 0.5 | nothing per-family — one base, × copies |
| `FamilySampledRates` (A→0.8, B→0.2) | **0.8** × 3 = 2.4 | **0.2** × 1 = 0.2 | the **base rate** |
| `Rates(0.5)` + `FamilyModifier` (A→1.6, B→0.4) | 0.5 × **1.6** × 3 = 2.4 | 0.5 × **0.4** × 1 = 0.2 | a **modifier** |

The last two rows are identical (2.4, 0.2): same effective rate, two spellings. Pick by intent — you
have actual per-family rate *values* → the **base** (`FamilySampledRates`); you want to rescale a
baseline, or **stack** heterogeneities (family × branch, adding a `LineageModifier`) → a **modifier**
(`FamilyModifier`), because modifiers compose and a base model does not.

**The clarity guard:** *per-family is not an opportunity.* The opportunity (`per="copy" | "lineage" |
"shared"`) is *how the rate scales with counts* — the middle term of `base × opportunity × modifiers`.
Per-family heterogeneity is *which families get which values* — it lives in the base or the modifiers.
A family can be `per="copy"` **and** carry its own rate; they are orthogonal slots.

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
| **Genomes** | duplication | d | per copy | family, lineage, carrying capacity | point (segment on ordered/nt) |
| **Genomes** | transfer (HGT) | t | per copy | family, donor lineage, recipient receptivity, lineage pair | point (segment on ordered/nt) |
| **Genomes** | loss | l | per copy | family, lineage | point (segment on ordered/nt) |
| **Genomes** | inversion | i | per copy | family, lineage | spans a run of genes / nt |
| **Trait** | discrete change (Mk) | qᵢⱼ | per lineage | lineage, hidden state | point |
| **Trait** | continuous drift (BM) | σ² | per lineage · time | lineage (relaxed), state | — (see note) |
| **Sequence** | substitution | μ × exchangeabilities | per nucleotide | lineage (clock), site (Γ), family, selection (dN/dS) | point |

!!! note "Continuous traits"
    A drifting continuous trait has no countable events and no extent — it is the limit of infinitely
    many infinitesimal steps. But it still fits the skeleton: the **base** is the variance accrued per
    unit time (σ²), the **opportunity** is per lineage, and **modifiers** (a relaxed clock, a
    state-dependent rate) rescale it exactly as elsewhere. Same grammar, no discrete jumps.

## The three genome resolutions, side by side

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
this family is a frequent mover (family modifier), the donor lineage sits in a high-HGT clade (lineage
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
