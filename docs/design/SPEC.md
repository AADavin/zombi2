# ZOMBI2 — Model & Vocabulary Specification

This document is the single source of truth for how ZOMBI2 is organised, what it may and may not do,
and the exact words used to describe it. When any other file — code, docstring, CLI help, manual
chapter, design doc, README — disagrees with this document, that other file is a **fossil** and must be
changed to match. 

---

## 1. The four levels

ZOMBI2 simulates evolution at four **levels**. Always name them with **words**, in this order:

**Species · Genomes · Sequences · Traits.**

- Never abbreviate to single letters. Letters flipped meaning since ZOMBI1 (there they meant sequence
  and tree), so they misread; and Traits and Tree both start "Tr".
- Level names are **plural** as labels and in the notation. Prose may use the singular for one instance
  ("a trait", "a genome"). Write "Sequences", never "Sigma".

Three of the levels form a chain, and traits branch off it:

```
Species → Genomes → Sequences      (a genome lives on the species tree; a sequence lives inside a gene)
Species → Traits                   (a trait lives on the species tree)
```

These "lives-on" connections are always present; they are **not** couplings you add. Because a sequence
lives inside a gene, it sees the species tree only through its gene tree, so its notation conditions on
Genomes, not Species.

---

## 2. How levels relate: independent, conditioned, joint

Everything evolves on the tree, so every level is already conditioned on it. The real question is how
two levels relate **to each other**. There are exactly three answers, taught with
probability-factorisation notation, where `P(B | A)` reads "B simulated given A":

| Relation | Notation (a trait and a genome) | What you run |
|---|---|---|
| **Independent** | `P(Traits\|Species) · P(Genomes\|Species)` | two commands, any order |
| **Conditioned** | `P(Traits\|Species) · P(Genomes\|Species, Traits)` | two commands, driver first, passed as a **file** |
| **Joint** | `P(Traits, Genomes \| Species)`, or `P(Species, Traits)` when the tree is grown | **one** command |

The load-bearing rule: **every factor you can write on its own is a run you can do on its own.** When
two factors collapse into one term, so do their runs.

- **Independent** — neither reads the other (independent *of each other*, not of the tree).
- **Conditioned** — one reads the other; the driver can be grown first and held fixed, so it is two runs
  in order, the driver passed as a file.
- **Joint** — neither can go first, so one run makes both; when the coupling feeds back into the species
  tree, the tree becomes an **output** and crosses the bar (`P(Species, Traits)`).

---

## 3. Which pairs can be conditioned or joined

Not every pair of levels can carry every kind of coupling. What a pair allows depends only on whether
one level lives on the other or they sit on separate branches:

| Level pair | Can be **conditioned**? | Can be **joined**? |
|---|---|---|
| Species – Genomes | no (a genome lives on the species tree) | yes — gene content drives speciation (tree grown) |
| Species – Traits | no (a trait lives on the species tree) | yes — a trait drives speciation (tree grown) |
| Genomes – Sequences | no (a sequence lives inside a gene) | in principle yes, deferred |
| Genomes – Traits | yes, either direction | yes — the two drive each other |
| Traits – Sequences | yes (a trait drives a sequence; the reverse is deferred) | no |
| Species – Sequences | no | no — too far apart to connect |

The generating rule: a pair can be **conditioned** only on separate branches (Genomes–Traits,
Traits–Sequences); it can be **joined** either by a level feeding back into its own tree (Species–Genomes,
Species–Traits, tree grown as an output) or by two separate-branch levels driving each other
(Genomes–Traits). Species and Sequences are too far apart to connect.

---

## 4. Joint models and naming

**The models that must be simulated jointly** (one run produces both levels):

- a trait drives speciation (the tree is grown)
- gene content drives speciation (the tree is grown)
- a trait drives speciation and also changes at speciation events (the tree is grown)
- gene content drives speciation, with a burst of gene change at each split (the tree is grown)
- a trait and gene content drive each other (the tree stays fixed)

The literature calls these models by acronyms (BiSSE, MuSSE, QuaSSE, HiSSE, ClaSSE, key innovation,
co-diversification, trait–gene feedback). **Those names are deprecated as structure, not hidden:**
section headings and prose describe each model by what it does; the acronyms never organise a chapter.
The class names remain in the code as the field's search terms.

---

## 5. Rates

Every event fires at a **rate**, and every rate is written the same way: a **scope** wrapped around a
**base**, times **modifiers**.

```
effective rate  =  scope(base)  ×  modifiers
```

- **base** — the speed of one event (how fast), in units of inverse time (`time⁻¹`).
- **scope** — how many independent copies, lineages, or sites the event applies to right now (per
  what); answering **"per what?"** is the crux. It wraps the base and contributes a dimensionless
  factor.
- **modifiers** — dimensionless context multipliers (by lineage, by family). They change *how fast*,
  never *how many*. "per" is the scope word; the modifier word is a **preposition** — `On` / `By` /
  `From` (the families, below) — so `PerLineage` is a scope and `ByLineage` a modifier.

"Per what" by level:

| Level | Counted per | "How fast" set by |
|---|---|---|
| Species | lineage | the diversification process |
| Genomes, gene tier | copy — lineage for origination | duplication / transfer / loss / inversion / transposition / translocation |
| Genomes, chromosome tier | chromosome — lineage for origination | fission / fusion / chromosome loss / chromosome origination |
| Sequences | site | the substitution rate (× a clock) |
| Traits | lineage | the trait model |

One rule generates that column: **an event that acts on something already there is counted per that
thing; an event that makes something from nothing is counted per lineage.** Origination is not an
arbitrary exception — there is no existing gene for it to be "per", and likewise no parent replicon
for a de-novo plasmid.

Time is imposed by the species tree, at the beginning of the **stem**.

**How a rate is written (same at every level):** the scope wraps (`PerCopy(0.2)`, `PerLineage(0.5)`,
`Global(1.0)` — `Global` capitalised, since `global` is a Python keyword) and modifiers multiply
(`0.2 * ByFamily(...)`, `1.0 * OnTotalDiversity(cap=100)`); a bare number uses the rate's natural scope,
so the common case is just `birth=1.0`. There is **no `per=` argument** — the scope lives on each rate.
Two rules: (a) `*` composes only dimensionless modifiers onto one base (multiplying two rates is
`time⁻²`, impossible by construction); (b) **"per" is reserved for scopes** — a modifier never starts
with "per".

**One written form, everywhere.** That expression is not Python syntax that the CLI then translates —
it is *the* way a rate is written, and the CLI and the parameters file take it **verbatim**:

```
birth = 1.0 * OnTime({0: 1.0, 3: 0.3})          # Python
--birth "1.0 * OnTime({0: 1.0, 3: 0.3})"        # the command line
birth = "1.0 * OnTime({0: 1.0, 3: 0.3})"        # a --params TOML value
```

A bare number stays a bare number in all three (`--birth 1.0`, `birth = 1.0`). The `mod.` / `scope.`
qualifiers Python needs are optional in the other two, so a manual snippet pastes in unchanged. There
is **no second notation** — no per-modifier flags, no nested parameter tables; adding a modifier must
never add a flag. (Read by `rates/parse.py`; it parses the expression, it does not evaluate code.)

**A level rejects the modifiers it does not wire.** A modifier a level has not implemented must
**raise**, never be silently ignored — an unwired modifier that returns a factor of 1.0 is a run that
is quietly not the model the user asked for. Each level therefore declares what it wires
(`WIRED_MODIFIERS`), the CLI's help is **built from that declaration** rather than hand-listed, and
the engine's own gate may be stricter still where a rate is wired more narrowly than the level.

**A driver's number is not always a rate multiplier.** `DrivenBy(source, mapping)` is the one coupling
mechanism (§2), and the **slot** it sits in decides what the mapping's number means. In a rate it is an
ordinary modifier: dimensionless, multiplying, changing *how fast*. In a **choice slot** it is a
**weight**, normalised across the candidates the choice is made over, so it changes neither how fast
nor how many — only **who**. Today the one choice slot is the genome level's `transfer_to`, the "who
receives" of a horizontal transfer:

```
transfer    = 0.1 * DrivenBy(habitat, {"competent": 3.0, "normal": 1.0})   # a rate:   how much transfer
transfer_to =       DrivenBy(habitat, {"competent": 3.0, "normal": 1.0})   # a weight: where it lands
```

The first changes the total amount of transfer; the second redistributes the same transfers. A choice
slot takes the modifier **on its own**, never `base * modifier` — there is no base, because there is no
rate. A weight of 0 means "cannot receive"; when every candidate weighs 0 the transfer cannot happen at
all, so the event does not fire.

**Banned rate words:** "propensity" (say *rate*); "opportunity" as a noun (say **scope**, or ask **"per
what?"**); "clock" for the scope (reserve **clock** strictly for the by-lineage substitution-rate
modifier at the sequences level). **modifier** names the third factor only.

**The modifier families.** A modifier's name begins with the preposition that fixes its family:

| Preposition | Family | The factor is… | Examples |
|---|---|---|---|
| `On` | covariate | a deterministic function of a measured quantity | `OnTime`, `OnTotalDiversity` |
| `By` | independent | an i.i.d. draw, one per unit — **no memory** (uncorrelated) | `ByLineage`, `ByFamily` |
| `From` | inherited | inherited along a genealogical edge — **continuous memory** (autocorrelated) | `FromParent` |

So the uncorrelated / autocorrelated split is `ByLineage` vs `FromParent`, and one modifier —
`FromParent` — is ClaDS (species), the autocorrelated clock (sequences), and variable-rates BM
(traits). Four rules: **fully-qualify an `On` covariate** (`OnTotalDiversity`, since the preposition
does not fix its scope); **one memory-structure per axis** (`By…` none / `From…` continuous / `Markov`
discrete — never two at once; orthogonal axes compose); **named exceptions** — discrete-memory
switching is `Markov` (a mechanism name, no preposition), and any future prepositionless mechanism is
named likewise; **outside the grammar** — a modifier multiplies *one* rate, so a value-level process
(the OU trait's `reverts_to` / `pull`) is a function argument, not a modifier.

---

## 6. Canonical vocabulary 

Left column is correct; right column is a fossil to purge.

| Use this | Not this (fossil) |
|---|---|
| Species, Genomes, Sequences, Traits (words, plural) | single letters, "Sigma" |
| level (one of the four) | tier |
| resolution — unordered / ordered / nucleotide | "level" for the genome sub-axis; `--genome-model` |
| independent / conditioned / joint | pipeline / coevolution (as the framing) |
| conditioning; joining; a joint model | coevolution (as a category) |
| rate; effective rate = scope(base) × modifiers | propensity |
| scope; "per what?" | opportunity |
| clock (the sequences by-lineage rate modifier only) | clock (for the count) |
| the four levels of ZOMBI2 (the layout) | the diamond |
| complete tree / extant | "reconstructed" (only once, as Nee's synonym); "pruned" as a noun |
| ZOMBI2 | Zombi2 |
| ZOMBI1 | ZOMBI-1, ZOMBI 1, ZOMBI(1); "Zombi" except in citations/URLs |

Literature model names: deprecated in the manual (footnote at most), class names kept in the code.

---

## 7. Naming and branding

- The tool is **ZOMBI2** (already consistent; do not change to "Zombi2"). The package/CLI token is
  lowercase `zombi2`.
- Version 1 is **ZOMBI1** (no space). "Zombi" survives only in citations and URLs. Reject "ZOMBI(1)".
- Book subtitle: **"Simulating the Evolution of Species, Genomes, Sequences and Traits."**
- Trees: **complete** (keeps extinct lineages) vs **extant** (the sampled survivors); output filenames
  are frozen (`_extant.nwk` is the extant tree, kept from ZOMBI1).
- Every node in a written tree carries a branch length, **the root included** — no exceptions. A
  forward run starts from one lineage, so the root's branch is its **stem**: real simulated time in
  which events happen. Writing a crown-rooted `)n0;` would discard it. For a species tree the stem
  runs from the origin to the first split; for a gene tree, from the family's **origination** to the
  founding gene's first event; for a **phylogram**, it is that same stem in substitutions, because
  the founding sequence is drawn at origination and evolves across the stem like any other branch.

