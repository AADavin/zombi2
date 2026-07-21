# ZOMBI2 manual — proposed detailed index

Round-1 restructure, drafted 2026-07-17; revised 2026-07-21. 9 chapters + 2 appendices, down from 17+1.
Provisional section titles; `⟨Dn⟩` marks a spot still gated by an open decision.

Conventions baked in: levels are **words** — Species, Genomes, Sequences, Traits — never single
letters (they flipped meaning since ZOMBI1). Levels relate as **independent / conditioned / joint**.
Every model chapter ends with **Outputs**; concept chapters are exempt essays.

---

## Part I — Getting started

### 1. Introduction  *(essay)*
- Why simulate: ground truth for methods
- What ZOMBI2 is: one platform, four levels, a fast engine
- What it can do — worked questions, biologically grounded, not one bullet per level ⟨C1.6; revised 2026-07-21⟩
- Installing it — brief; the rest lives in the README ⟨absorbed from the old ch3⟩
- For the impatient — the same short run in CLI **and** Python ⟨C1.11⟩

### 2. A tour of ZOMBI2  *(essay — the vocabulary chapter)*
- The four levels of ZOMBI2: Species, Genomes, Sequences, Traits — with the layout figure (Genomes **left**, Traits **right**; the layout is named nothing) ⟨D9, D10⟩
- The genome **resolution** dial: unordered → ordered → nucleotide (framing only; detail in Part II)
- How levels relate: **independent, conditioned, joint** — the composition figure + probability notation
- How rates work: how many, how fast — scope(base) × modifiers ⟨D3: section here, or its own chapter?⟩
- The ZOMBI2 vocabulary — a short glossary (level, rate, resolution, complete/extant) ⟨C2.11, D6⟩

*(A separate "Getting started" chapter was written and then **cut** on 2026-07-21 as redundant: its first run duplicated Ch1's, its output-folder tour duplicated the level chapters and Appendix B. Installation moved into Ch1.)*

---

## Part II — The four levels

> **Common structure ⟨D7⟩.** Every level chapter teaches its concept and models first, then closes
> the same way: **The objects** (the result type and its API) → **Usage from Python** →
> **Usage from the CLI** → **Outputs**. Only the chapter-specific tail notes are listed below.

### 3. Species trees  *(merges old 4 + 5)*
- **Forward and backward** — the two modes, with a figure; opens the chapter ⟨C4.2/C4.6⟩
- The birth–death process — the base (Yule as a note: death = 0, not a class) ⟨C4.3⟩
- **What the rate depends on** — one table, one process: nothing (constant) · time (episodic) · ancestry (the rate drifts as lineages split — "ClaDS" a footnote) · standing diversity · clade ⟨D5⟩
- **Scheduled shocks** — mass extinctions: part of the process, compose with any rate ⟨D5⟩
- **Sampling — what you observe** — the extant fraction, and fossils (FBD); the true tree is unchanged ⟨D5, C5.6⟩
- Ghost lineages — a tree transform, not a diversification model
- *(trait- and gene-driven diversification live in ch9 Coupling levels — one forward pointer)*
- **Tail:** the `Tree` object → Python → CLI → Outputs *(complete & extant Newick)* ⟨D6, D7⟩

### 4. Genomes I — Unordered  *(old 7 + the genealogy half of 8)*
- What an unordered genome is: a multiset of gene families — the base resolution
- The four events: origination, duplication, transfer, loss
- Rates: **per copy vs per lineage** (the scope), and **shared vs per-family** (heterogeneity) — kept distinct ⟨C7.5, D3⟩
- Transfers — recipient choice is the **mechanics** (how the event resolves): emission and receptivity ⟨C7.6, D5⟩, additive vs replacement, self-transfer
- Gene conversion — its own section, off by default ⟨C7.2⟩
- Bounding growth: hard cap and carrying capacity
- Gene trees and the event log — the event log is a full genealogy; reconstructing the gene trees; reconciliation ⟨from old ch8⟩
- **Tail:** the genome & profile objects → Python → CLI → Outputs *(profile matrix, sparse; event trace)* ⟨D7⟩

*(The old "Genomes II — Structured" was **split in two** on 2026-07-21: one chapter was carrying gene order, rearrangements, the karyotype and the chromosome network at once. The seam is order *within* a chromosome versus the chromosomes themselves.)*

### 5. Genomes II — Gene order
- From families to positions: what structure adds
- Events act on **segments** — the extension, and why neighbours share a history
- Rearrangements: inversion, transposition, translocation
- **Tail:** the `OrderedGenomesResult` object → Python → CLI → Outputs *(gene order, rearrangements)* ⟨D7⟩

### 6. Genomes III — Chromosomes
- The karyotype: number and topology (circular / linear)
- The chromosome tier: **fission, fusion, chromosome origination, chromosome loss**
- The **chromosome network** — species tree ⊃ chromosome network ⊃ gene trees; a reticulating graph, written as an edge list
- **Tail:** Python → CLI → Outputs *(chromosome events)* ⟨D7⟩

### 7. Sequence evolution  *(merges old 12 + 13)*
- The idea: one branch length = one expected substitution — the phylogram is the deliverable ⟨C12.3, was buried at the end⟩
- From a time-tree to a phylogram: what a relaxed clock does (it never changes topology)
- Rate variation, two sources: shared across genes in a lineage, and per family ⟨C12.2 reworded⟩
- The clocks — table first ⟨C13.3⟩, then uncorrelated and autocorrelated families
- Substitution models: DNA (JC/K80/HKY/GTR +Γ) and protein (LG/WAG…)
- One sentence: a clock is a trait (the Brownian identity) — cross-reference, no relocation ⟨C12.1⟩
- **Tail:** the sequence & alignment objects → Python → CLI → Outputs *(alignments, phylograms)* ⟨D7⟩

### 8. Trait evolution  *(old 11)*
- What a trait is, and that it rides **any** tree — species or gene ⟨C11.2⟩
- **Continuous vs discrete** — the real taxonomy ⟨C11.4⟩
- Continuous: Brownian motion, Ornstein–Uhlenbeck, early-burst, multivariate
- Discrete: Mk, correlated binary (Pagel), hidden rate classes, threshold
- Adaptation to regimes: multi-optimum OU
- Pagel's tree transforms
- Historical biogeography (DEC) — a discrete model with cladogenetic transitions ⟨C11.5⟩
- The change history — the within-branch stochastic map, exposed ⟨C11.6⟩
- **Tail:** the `TraitsResult` object → Python → CLI → Outputs *(trait values, change history, trait tree)* ⟨D7, C11.10⟩

---

## Part III — Coupling the levels

### 9. Conditioning and joining  *(merges old 9 Conditioning + 10 Joint + 11 Nulls into one chapter)*
- The one idea: a rate driven by another level — `mod.DrivenBy(source, mapping)`, a parameter that stops being a number you type
- Conditioned vs joint is **one** distinction: *can the driver be grown first?* — a file source (two commands, ordered) vs a live-level source (one command). NOT "does it change the tree" ⟨D14, D15⟩
- **Conditioned** — the driver is a file, grown first and handed over: a trait drives gene loss (Traits → Genomes); a gene drives a trait's optimum (Genomes → Traits); a trait drives selection / clock speed (Traits → Sequences)
- **Joint** — the driver is a live level, grown alongside (`zombi2 joint`): state-dependent diversification, a trait drives speciation (BiSSE / MuSSE / QuaSSE / HiSSE as footnotes) ⟨D11⟩; key-innovation, gene content drives speciation
- *(v1 ships only tree-changing joint models — a scope note, **not** the framing; tree-fixed trait–gene feedback is deferred to experimental)*
- *(model-flag couplings that only re-read the tree — cladogenetic traits, punctuational genomes — stay in their level chapter; pointers here)*
- Outputs — conditioned writes the driver alongside the target; a joint run writes both levels
- **Closing section — Nulls, a recipe not a feature:** the tree manufactures patterns, so a baseline must be simulated on the same tree; three one-line recipes (drop the coupling → independent; swap for `ByLineage` → CID; shuffle the driver file → permutation). ZOMBI2 generates the baseline; the user owns the statistic. *(No `--null` flag, no null subsystem.)*

---

## Appendices

### Appendix A — The Gillespie algorithm  *(unchanged)*
### Appendix B — Output files, in full  *(new: the shared write-parts reference/table)* ⟨D7, C10.6⟩
