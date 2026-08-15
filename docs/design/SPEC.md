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
- **Joint** — neither can go first, so **one run simulates both**. That is the whole of it. When the
  species tree is one of the two, it comes **out** of the run and crosses the bar
  (`P(Species, Traits)`); when it is not, the run takes a tree the way every other level does.

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
| Traits – Sequences | yes, either direction | yes — the two drive each other |
| Species – Sequences | no | no — too far apart to connect |

The generating rule: a pair can be **conditioned** only on separate branches (Genomes–Traits,
Traits–Sequences); it can be **joined** either by a level feeding back into its own tree (Species–Genomes,
Species–Traits, where the tree comes out of the run) or by two separate-branch levels driving each
other (Genomes–Traits, Traits–Sequences). Species and Sequences are too far apart to connect.

The same rule settles Traits–Sequences, which reads *yes* above because **both directions exist**: a
trait drives a gene's substitution rate, and that gene's composition drives the trait. Either on its
own is conditioning. Both at once is a cycle, and a cycle is joint.

This table is about the **model**, not about what is built. A pair marked *yes* is one the framework
permits; whether an engine implements it is a separate question, and the manual answers that one.

**A level can also drive itself.** Two traits on one tree, one gene family driving another: the
participants are at the same level, and nothing above changes. Driving is **one thing driving another**,
and the relation is fixed by the same question as always — *can the driver be finished before the driven level
starts?* If it can, that is conditioning, whichever levels the two sit at. So:

- **acyclic within a level** — trait A drives trait B, an earlier genome run drives a later one. Ordinary
  conditioning: `P(A|Species) · P(B|Species, A)`, two runs in order.
- **cyclic, or driven by an aggregate of the level itself** — every family's rate reading total genome
  size, sites evolving in each other's context. The level's units stop being independent, so this is a
  **joint** model rather than a knob. Whether it needs a new engine (§9) is a separate question, and
  the answer is per level: the family genome engine already races every family in one loop and takes
  this unchanged, so what it costs there is the parallel and streaming engines, which evolve one
  family per process. The sequences level, where a site's rate would read other sites, does need one.

Cross-level and within-level are therefore not different mechanisms, and the vocabulary does not fork.

---

## 4. Joint models and naming

**The models that must be simulated jointly** (one run produces both levels):

- a trait drives speciation (the tree is grown)
- gene content drives speciation (the tree is grown)
- a trait drives speciation and also changes at speciation events (the tree is grown)
- gene content drives speciation, with a burst of gene change at each split (the tree is grown)
- a trait and gene content drive each other (the tree stays fixed)
- a trait and a gene's sequence drive each other (the tree stays fixed)

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
- **modifiers** — dimensionless context multipliers (among lineages, among families). They change
  *how fast*, never *how many*. **Nobody writes a modifier**: a modifier is what a **verb** records,
  and the verbs are what a rate is written with (below). "per" is the scope word and nothing else:
  `PerLineage` is a scope, while a value that varies at random is written
  `varying_among('lineages', …)` — the unit **plural**, because a value varies *among* lineages
  rather than being counted *per* one.

**One verb replaces the base rather than multiplying it: `set_by`.** Two things are stated
absolutely in the literature rather than as a multiple — a driven rate ("the loss rate is 1.0 in
caves") and a trait's optimum — and writing those as a factor means inventing a background nobody
stated and dividing by it. So `loss = PerCopy().set_by(habitat, {"cave": 1.0, "surface": 0.25})`
takes **no base in front**: the driver supplies the whole number, in the rate's own units. The scope
is untouched — a per-copy rate set to 1.0 is still 1.0 *per copy* — because `set_by` answers *how
fast*, never *per what*. A rate carries **one** `set_by` and any number of multiplying verbs, and the
`set_by` is written first, because everything to its left is a base it would discard.

A `set_by` reads a driver like the other driven verbs, and a level must therefore declare it
**separately**: replacing a base is a capability an engine has or has not, and a gate that admitted it
alongside a `scaled_by` would accept it at four levels that cannot honour it.

"Per what" by level:

| Level | Counted per | "How fast" set by |
|---|---|---|
| Species | lineage | the diversification process |
| Genomes, gene rates | copy — lineage for origination | duplication / transfer / loss / inversion / transposition / translocation |
| Genomes, chromosome rates | chromosome — lineage for origination | fission / fusion / chromosome loss / chromosome origination |
| Sequences | site | the substitution rate (× a clock) |
| Traits | lineage | the trait model |

One rule generates that column: **an event that acts on something already there is counted per that
thing; an event that makes something from nothing is counted per lineage.** Origination is not an
arbitrary exception — there is no existing gene for it to be "per", and likewise no parent replicon
for a de-novo plasmid.

Time is imposed by the species tree, at the beginning of the **stem**.

**How a rate is written (same at every level):** the scope is the entry point (`PerCopy(0.2)`,
`PerLineage(0.5)`, `Global(1.0)` — `Global` capitalised, since `global` is a Python keyword) and the
**verbs chain onto it** (`PerCopy(0.2).varying_among('families', LogNormal(0.0, 0.5))`,
`PerLineage(1.0).scaled_by(TotalDiversity(cap=100))`); a bare number uses the rate's natural scope,
so the common case is just `birth=1.0`. There is **no `per=` argument** — the scope lives on each rate.
Two rules: (a) a verb composes only dimensionless factors onto one base, so two rates can never be
multiplied together (that would be `time⁻²`, impossible by construction); (b) **"per" is reserved for
scopes** — the unit a value varies among is written `among`, in the plural, never `per`.

**One written form, everywhere.** That expression is not Python syntax that the CLI then translates —
it is *the* way a rate is written, and the CLI and the parameters file take it **verbatim**:

```
birth = PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})          # Python
--birth "PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})"        # the command line
birth = "PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})"        # a --params TOML value
```

A bare number stays a bare number in all three (`--birth 1.0`, `birth = 1.0`). Every name the written
form may call is importable from `zombi2.params`, so Python needs no module qualifier; where one is
written (`scope.PerLineage(1.0)`) it is ignored by the other two, so a manual snippet pastes in
unchanged. There is **no second notation** — no per-modifier flags, no nested parameter tables;
adding a driver, a law or a mapping must never add a flag. (Read by `params/parse.py`; it parses the
expression, it does not evaluate code. Because the written form is a chain of method calls, the
parser whitelists the **verb names** as well as the class names.)

**A level rejects the modifiers it does not support.** A modifier a level does not support must
**raise**, never be silently ignored — a modifier that returns a factor of 1.0 because nothing reads
it is a run that is quietly not the model the user asked for. Each level therefore declares what it
takes (`IMPLEMENTED_MODIFIERS`), the CLI's help is **built from that declaration** rather than hand-listed,
and the engine's own gate may be stricter still where a rate takes less than the level does.

Two different things get rejected, and the message must say which. A few combinations are
**meaningless** — `varying_among('families', …)` on a species or trait rate, where there are no gene
families to draw a factor among — and no implementation would make them mean anything; say so, and
name the argument the modifier does belong on. The rest are **not implemented yet**, which is a
statement about the code and not about the model; say that plainly and do not dress it up as a rule.

**A driven parameter is the one mechanism** for both conditioning and joining (§2), within a level as
much as across two (§3). Whichever verb writes it, it takes a `driver` and a `mapping`. `driver` says
which thing is read, never how the run is organised: a **finished result** (an object in Python, its
written log across two commands) makes the run
conditioned; the **name of a level growing beside it** makes the run joint. A **`Clade`** is
conditioned too, and the limiting case of it: its value is a fact about the tree the run is already
walking, so there is nothing to grow first and no file — and a joint run refuses it, because a
growing tree has no clades yet. A driver read from a file
and the same driver held in memory are the same model, so they are the same modifier. `mapping` says
what the value becomes, and there are four shapes: a **`Table`** over named states, a **`Curve`** over
a number, a **`Scalar`** log-link `exp(strength · value)`, and a **`Between`**, a weight per (donor
group, recipient group) pair that only a choice takes (below).

**What the mapping's number means depends on what it is attached to.** On a rate or an **extent** (§6)
it is an ordinary modifier: dimensionless, multiplying, changing *how fast* or *how much* — unless the
verb is `set_by`, whose number carries the rate's own units and replaces the base. An extent takes no
`set_by`: it is already an absolute size drawn from a distribution, so there is no base to replace. On a
**choice** — an argument that decides *who*, not how fast or how many — it is a **weight**,
normalised across the candidates:

```
transfer    = PerCopy(0.1).scaled_by(habitat, {"competent": 3.0, "normal": 1.0})    # a rate:   how often transfer
transfer_to = Recipients().weighted_by(habitat, {"competent": 3.0, "normal": 1.0})  # a choice: which lineage receives
```

The genome level's `transfer_to`, the "who receives" of a horizontal transfer, is the only such
argument today. A choice is written from its own entry point, `Recipients()`, and never carries a
base, because there is no rate to have one. A weight of 0 means "cannot receive"; when every candidate weighs 0 the event
does not fire at all. A rate, an extent, a choice and a **model** are the four kinds of **target** —
the four things a driver can be attached to — and they are not the same as the three questions this section opens with:
*where* an event starts is drawn by the engine and takes no modifier, and a choice picks the lineage
that receives, not the segment.

The first three take a **factor**: the driver supplies a number and the target is multiplied by it. A
model takes none. A substitution model is an object rather than a quantity, so a driver on one
*selects* it instead of scaling it, and it is written with `set_by` — the verb that replaces a base
rather than multiplying one — `set_by(clade, {"endo": hky85(frequencies=...), "rest": hky85()})`.
Every branch shares one alphabet, so what varies along the tree is the process over the states, never
the states themselves.

A weight may read **both** ends: a **kernel** over `(donor group, recipient group)` pairs
(`Between({...})`) steers transfer *between* groups rather than only *into* one. The groups come from
the tree (named clades, reading no other level) or from a trait
(`Recipients().weighted_by(trait, Between(...))`). A
kernel only redistributes who receives, so one on a *rate* or an *extent* is refused: both are read on
one lineage and have no donor to condition on.

**Banned rate words:** "propensity" (say *rate*); "opportunity" as a noun (say **scope**, or ask **"per
what?"**); "clock" for the scope (reserve **clock** strictly for the by-lineage substitution-rate
modifier at the sequences level). **modifier** names the third factor only.

**A drawn value takes any distribution that can state its mean.** The **law** written beside the unit
says what the value is, and a bare distribution is the common case — the value itself, drawn once and
held: ``varying_among('families', LogNormal(0.0, 0.5))``, with any of the built-ins ``Fixed``,
``Exponential``, ``Gamma``, ``LogNormal``, ``Uniform``, ``Geometric``. ``Drift(dist)`` is the same
menu read as the per-split *step* instead. Each law owns and documents its own argument, so nobody
has to infer the role from the slot. A bare callable or a scipy frozen distribution is refused here,
though an **extent** takes either, because an extent is a size used as written rather than a
multiplier normalised to mean 1 (below).

Whatever the distribution, **the draw is normalised to mean 1**, by dividing by that distribution's
own mean. A drawn value is a *multiplier*, and one that does not average to 1 changes what the base
means — a base of 0.25 would stop being the average rate. So a distribution's **location is
normalised away** and what it contributes is its *shape*; ``Exponential(1.0)`` and ``Exponential(7.0)``
are one modifier. A distribution that cannot state its mean is refused rather than normalised by a
guess. A number that *is* the rate rather than a factor is a `set_by`, where nothing is normalised.

**The modifier families.** Four kinds, and a modifier's kind says who produces its number:

| Kind | The factor is… | Written |
|---|---|---|
| covariate | a deterministic function of a measured quantity | `changing_at({…})`, `scaled_by(TotalDiversity(cap=…))` |
| drawn | an i.i.d. draw, one per unit — **no memory** (uncorrelated) | `varying_among(unit, dist)` |
| inherited | the parent's, perturbed — **continuous memory** (autocorrelated) | `varying_among(unit, Drift(dist))` |
| driven | the state of another simulated thing, read as the run walks the tree | `scaled_by`, `weighted_by`, `set_by` |

**The unit is an argument, not a class.** A draw among families and a draw among lineages are one model at
two attachments, so a unit nobody has carried yet needs no name invented for it — which is what
`varying_among('chromosomes', …)` is. `ByFamily`, `ByLineage` and `FromParent` were names for three of
those cells and are **removed**: they were our coinages rather than the field's, and the field's own
names — the relaxed clock, ClaDS, rate heterogeneity across families — are what the prose uses.

**One object is one draw.** A drawn or inherited value is drawn once per unit by the *engine* and
kept. Which rates share that value is decided by **what you wrote, not by what the numbers
are**: one `Random` object read by several rates is one draw, shared between them, and two separately
built ones are two draws even when their arguments match.

```python
speed       = Random('families', LogNormal(0.0, 0.5))
duplication = PerCopy(0.2).varying_among(speed)                            # one draw: a fast family is fast at both
loss        = PerCopy(0.25).varying_among(speed)
duplication = PerCopy(0.2).varying_among('families', LogNormal(0.0, 0.5))  # two draws: the two rates vary independently
loss        = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))
```

That rule is the whole of it, and it replaces the separate family-wide argument this used to need.
Sharing is by **identity**, never by equality, so the question a reader has to answer is only ever
*did you write one thing or two?*. It follows that the **text form cannot express sharing** — two
flags parse to two objects — so a shared draw is Python-only until the written form can name a value.

So the uncorrelated / autocorrelated split is a bare distribution against `Drift(dist)`, and one
law — `varying_among('lineages', Drift(dist))` — is ClaDS at Species, the autocorrelated clock at
Sequences, and variable-rates BM at Traits. Three rules for the next one:

- **Fully qualify a measured driver** (`TotalDiversity`), since a bare word does not fix what is measured.
- **One memory structure per axis**: drawn has none, inherited has continuous memory, and a rate
  carries one or the other on a unit, never both. Several of the *same* kind compose and multiply, as
  any modifiers do. A discrete-memory mechanism would be named for the mechanism (`Markov`); none is
  implemented.
- **A verb multiplies one rate**, except `set_by`, which replaces its base. A process on the
  *value* rather than the rate — the OU trait's `reverts_to` / `pull` — is a function argument, not a
  modifier.

---

## 6. Extents

An event that acts on a **segment** — a run of genes at the ordered resolution, an arc of DNA at the
nucleotide one — carries an **extent**: how much it takes once it has started. An extent is written the
way a rate is, **minus the scope**, because it is already an absolute quantity and has no "per what?"
to answer:

```
extent  =  base × modifiers
```

- **base** — a number (the mean) or a distribution over sizes, written from `Extent(...)` when a verb
  is chained onto it.
- **modifiers** — the same dimensionless multipliers a rate takes, written with the same verbs
  (`changing_at`, `varying_among`, `scaled_by`).

One word, **extent**, at every resolution and for every event. Its **unit is set by the resolution** —
genes at ordered, base pairs at nucleotide — and needs no word of its own, because fixing that unit is
what a resolution *is*.

**A rate counts starts, not hits.** A rate says how often an event *begins*. A unit lying inside a
segment is therefore affected more often than the rate reads — about `rate × E[extent]` times — since it
is taken whenever any event begins on a segment that covers it. The two axes **multiply**, and quoting
one without the other describes nothing.

**A modifier attaches either to the lineage or to the contents.** `changing_at`,
`varying_among('lineages', …)` and a `scaled_by` reading a trait attach to the **lineage**: at any
instant they are uniform across that lineage's
whole genome, so they compose with any extent unchanged. `varying_among('families', …)` attaches to the **contents**, and a
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
| driver — the value a driven parameter reads, per lineage, as the run walks the tree (its first argument): grown by another level, or read off the tree itself (a clade) | source (for that argument); signal |
| target — what a factor is attached to: a rate, an **extent**, or a **choice** | "a target is a rate"; target (for the driven level — say *the driven level*) |
| choice — the target that decides who receives (`transfer_to`) | "which one"; slot |
| mapping — `Table` / `Curve` / `Scalar` / `Between` | response (the coevolve word) |
| weight — a `transfer_to` number, normalised across candidates | multiplier (there); a base in front of `Recipients()` (there) |
| rate; effective rate = scope(base) × modifiers | propensity |
| scope; "per what?" — name a group of rates by the scope they share: **the chromosome rates** (`PerChromosome`), which are not a resolution, since they exist at ordered and nucleotide alike | opportunity; "the chromosome tier" |
| verb — what reading a driver does to a parameter: `scaled_by` / `set_by` / `weighted_by`, plus the two shortcuts `varying_among` / `changing_at` | `ScaledBy` / `SetBy` / `Weights` / `OnTime` as names to call; `*` composing a rate |
| the unit a value varies among, plural: `'lineages'` / `'families'` / `'copies'` / `'sites'` / `'chromosomes'` | `per='family'`; `spread=` for a distribution |
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

