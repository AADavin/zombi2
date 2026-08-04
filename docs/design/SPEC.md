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

These "lives-on" connections are always present; they are **not** a relation you add. Because a sequence
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
| **Conditioned** | `P(Traits\|Species) · P(Genomes\|Species, Traits)` | two runs, driver first |
| **Joint** | `P(Traits, Genomes \| Species)`, or `P(Species, Traits)` when the tree is grown | **one** run |

The load-bearing rule: **every factor you can write on its own is a run you can do on its own.** When
two factors collapse into one term, so do their runs.

- **Independent** — neither reads the other (independent *of each other*, not of the tree).
- **Conditioned** — one reads the other; the driver can be grown first and held fixed, so it is two runs
  in order.
  From Python the finished driver is handed over as the **object**; across two commands it is handed
  over as a **file**, because that is all a second process can read. The file is the CLI's way of
  passing it, not part of the model — do not define conditioning by it.
- **Joint** — neither can go first, so one run makes both; when a level feeds back into the species
  tree, the tree becomes an **output** and crosses the bar (`P(Species, Traits)`).

---

## 3. Which pairs can be conditioned or joined

Not every pair of levels can be conditioned or joined. What a pair allows depends only on whether
one level lives on the other or they sit on separate branches:

| Level pair | Can be **conditioned**? | Can be **joined**? |
|---|---|---|
| Species – Genomes | no (a genome lives on the species tree) | yes — gene content drives speciation (tree grown) |
| Species – Traits | no (a trait lives on the species tree) | yes — a trait drives speciation (tree grown) |
| Genomes – Sequences | no (a sequence lives inside a gene) | in principle yes, deferred |
| Genomes – Traits | yes, either direction | yes — the two drive each other |
| Traits – Sequences | yes, either direction | no |
| Species – Sequences | no | no — too far apart to connect |

The generating rule: a pair can be **conditioned** only on separate branches (Genomes–Traits,
Traits–Sequences); it can be **joined** either by a level feeding back into its own tree (Species–Genomes,
Species–Traits, tree grown as an output) or by two separate-branch levels driving each other
(Genomes–Traits). Species and Sequences are too far apart to connect.

This table is about the **model**, not about what is built. A pair marked *yes* is one the framework
permits; whether an engine implements it is a separate question, and the manual answers that one.

**A level can also drive itself.** Two traits on one tree, one gene family driving another: the
participants are at the same level, and nothing above changes. Driving is **one thing driving another**,
and the relation is fixed by the same question as always — *can the driver be finished before the target
starts?* If it can, that is conditioning, whichever levels the two sit at. So:

- **acyclic within a level** — trait A drives trait B, an earlier genome run drives a later one. Ordinary
  conditioning: `P(A|Species) · P(B|Species, A)`, two runs in order.
- **cyclic, or driven by an aggregate of the level itself** — every family's rate reading total genome
  size, sites evolving in each other's context. The level's units stop being independent, so this is a
  joint model and a new engine (§9), not a knob.

Cross-level and within-level are therefore not different mechanisms, and the vocabulary does not fork.

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

Every event answers three questions: **how often** it starts (its rate, this section), **where** it
starts (the copy, position or lineage it is drawn on), and **how much** it takes (its **extent**, §6).
At Species, at Traits, and at the genome level's family resolution the third answer is always "one", so
the rate is the whole story. Once an event acts on a **segment** the three come apart, and a process is
only described when both the rate and the extent are given.

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
  never *how many*. "per" is the scope word; a modifier is named for its family — `On` / `By` / `From`,
  plus `DrivenBy` (below) — so `PerLineage` is a scope and `ByLineage` a modifier.

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

**A level rejects the modifiers it does not support.** A modifier a level does not support must
**raise**, never be silently ignored — a modifier that returns a factor of 1.0 because nothing reads
it is a run that is quietly not the model the user asked for. Each level therefore declares what it
takes (`IMPLEMENTED_MODIFIERS`), the CLI's help is **built from that declaration** rather than hand-listed,
and the engine's own gate may be stricter still where a rate takes less than the level does.

Two different things get rejected, and the message must say which. A few combinations are
**meaningless** — `ByFamily` on a species or trait rate, where there are no gene families to draw a
factor per — and no implementation would make them mean anything; say so, and name the argument the
modifier does belong on. The rest are **not implemented yet**, which is a statement about the code
and not about the model; say that plainly and do not dress it up as a rule.

**`DrivenBy(source, mapping)` is the one mechanism** for both conditioning and joining (§2), within a
level as much as across two (§3). `source` says which thing is read, never how the run is organised: a
**finished result** (an object in Python, its written log across two commands) makes the run
conditioned; the **name of a level growing beside it** makes the run joint. A driver read from a file
and the same driver held in memory are the same model, so they are the same modifier.

**What the mapping's number means depends on what it is attached to.** On a rate it is an ordinary
modifier: dimensionless, multiplying, changing *how fast*. On a **choice** — an argument that decides
*who*, not how fast or how many — it is a **weight**, normalised across the candidates:

```
transfer    = 0.1 * DrivenBy(habitat, {"competent": 3.0, "normal": 1.0})   # a rate:   how much transfer
transfer_to =       DrivenBy(habitat, {"competent": 3.0, "normal": 1.0})   # a weight: where it lands
```

The genome level's `transfer_to`, the "who receives" of a horizontal transfer, is the only such
argument today. A choice takes the modifier **on its own**, never `base * modifier`, because there is
no rate to have a base. A weight of 0 means "cannot receive"; when every candidate weighs 0 the event
does not fire at all.

A weight may read **both** ends: a **kernel** over `(donor group, recipient group)` pairs
(`Between({...})`) steers transfer *between* groups rather than only *into* one. The groups come from
the tree (named clades, reading no other level) or from a trait (`DrivenBy(trait, Between(...))`). A
kernel only redistributes who receives, so one on a *rate* is refused: a rate has no donor to read it
on.

**Banned rate words:** "propensity" (say *rate*); "opportunity" as a noun (say **scope**, or ask **"per
what?"**); "clock" for the scope (reserve **clock** strictly for the by-lineage substitution-rate
modifier at the sequences level). **modifier** names the third factor only.

**The modifier families.** A modifier's name begins with the preposition that fixes its family:

| Preposition | Family | The factor is… | Examples |
|---|---|---|---|
| `On` | covariate | a deterministic function of a measured quantity | `OnTime`, `OnTotalDiversity` |
| `By` | independent | an i.i.d. draw, one per unit — **no memory** (uncorrelated) | `ByLineage`, `ByFamily` |
| `From` | inherited | inherited along a genealogical edge — **continuous memory** (autocorrelated) | `FromParent` |
| — | driver | the state of another simulated thing, read as the run walks the tree | `DrivenBy` |

`DrivenBy` sits outside the preposition scheme deliberately: the others say what kind of *function* the
factor is, while this one says the factor comes from somewhere else entirely. Naming it `On…` would
file it as a covariate and lose that.

So the uncorrelated / autocorrelated split is `ByLineage` vs `FromParent`, and one modifier —
`FromParent` — is ClaDS (species), the autocorrelated clock (sequences), and variable-rates BM
(traits). Three rules for naming the next one:

- **Fully qualify an `On` covariate** (`OnTotalDiversity`), since the preposition does not fix its scope.
- **One memory structure per axis**: `By…` none, `From…` continuous, and never two at once. Orthogonal
  axes compose. A discrete-memory mechanism would be named for the mechanism rather than a preposition
  (`Markov`); none is implemented.
- **A modifier multiplies one rate.** A process on the *value* rather than the rate — the OU trait's
  `reverts_to` / `pull` — is a function argument, not a modifier.

---

## 6. Extents

An event that acts on a **segment** — a run of genes at the ordered resolution, an arc of DNA at the
nucleotide one — carries an **extent**: how much it takes once it has started. An extent is written the
way a rate is, **minus the scope**, because it is already an absolute quantity and has no "per what?"
to answer:

```
extent  =  base × modifiers
```

- **base** — a number (the mean) or a distribution over sizes.
- **modifiers** — the same dimensionless multipliers a rate takes (`On` / `By` / `From`, `DrivenBy`).

One word, **extent**, at every resolution and for every event. Its **unit is set by the resolution** —
genes at ordered, base pairs at nucleotide — and needs no word of its own, because fixing that unit is
what a resolution *is*.

**A rate counts starts, not hits.** A rate says how often an event *begins*. A unit lying inside a
segment is therefore affected more often than the rate reads — about `rate × E[extent]` times — since it
is taken whenever any event begins on a segment that covers it. The two axes **multiply**, and quoting
one without the other describes nothing.

**A modifier attaches either to the lineage or to the contents.** `OnTime`, `ByLineage`, `FromParent`
and a trait-`DrivenBy` attach to the **lineage**: at any instant they are uniform across that lineage's
whole genome, so they compose with any extent unchanged. `ByFamily` attaches to the **contents**, and a
segment has several — so a content-attached modifier must weight the **segment, by what it covers**,
never the position the event started from. Weighting the start applies a family's own rate to its
*neighbours*, and the neighbourhood is reshuffled by every rearrangement, so the parameter would not
even mean a fixed thing over a run.

A resolution that does not support an extent, or a modifier on one, **raises** — the §5 rule, unchanged.

**Banned:** "extension" and "length" for this quantity (say **extent**); "size" for the same.

---

## 7. Canonical vocabulary 

Left column is correct; right column is a fossil to purge.

| Use this | Not this (fossil) |
|---|---|
| Species, Genomes, Sequences, Traits (words, plural) | single letters, "Sigma" |
| level (one of the four) | tier |
| resolution — family / ordered / nucleotide | "level" for the genome sub-axis; "unordered"; `--genome-model` |
| independent / conditioned / joint | pipeline / coevolution (as the framing) |
| conditioning; joining; a joint model | coevolution (as a category) |
| conditioning and joining (when the pair needs one name) | coupling (as the framing, a category or a level); the verb stays — "a transfer couples two lineages" |
| rate; effective rate = scope(base) × modifiers | propensity |
| scope; "per what?" | opportunity |
| extent — how much a segmental event takes | extension; length (for this quantity); size |
| clock (the sequences by-lineage rate modifier only) | clock (for the count) |
| the four levels of ZOMBI2 (the layout) | the diamond |
| complete tree / extant | "reconstructed" (only once, as Nee's synonym); "pruned" as a noun |
| ZOMBI2 | Zombi2 |
| ZOMBI1 | ZOMBI-1, ZOMBI 1, ZOMBI(1); "Zombi" except in citations/URLs |

Literature model names: deprecated in the manual (footnote at most), class names kept in the code.

---

## 8. Naming and branding

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

---

## 9. Adding a model

A new model **belongs to its level** (Species, Genomes, Sequences, Traits): a module in that level's
package and, when it ships, a flag (or a `--model` value) on that level's command — never a new command,
since the CLI is one command per level.

**Maturity is a tag, not a place.** A young model carries an `experimental` marker in its docstring — no
separate folder or command, no `gallery`/`sandbox`. It is Python-first and reaches the CLI only once you
promise to keep it; half-built work stays on a branch, not in the package.

**Graduation moves nothing:** drop the tag, add the flag, write the manual section, add the outputs to
Appendix B — the file was in the level's package all along.

A model that **breaks independence** — the level's units stop being independent (families that affect
each other, sites that evolve in context) — is a new engine, not a knob: it needs its own evolve path in
the level, and graduates by handing the core one *general* capability the specific model is then one use
of.

