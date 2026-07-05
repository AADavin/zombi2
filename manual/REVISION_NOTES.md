# Manual revision — working backlog

Adrián's feedback, organized so nothing is lost. Status: [ ] todo · [~] in progress ·
[x] done · [?] needs discussion/decision.

## ⭐ Outstanding questions for Adrián (when you're back)

Decisions I couldn't make for you (the rest I'm proceeding on with sensible defaults):

1. **"liability" term (§12.3.4).** It *is* the standard term — Wright's/Felsenstein's threshold
   model uses "liability" for the latent continuous variable. Keep it (with a one-line gloss),
   or use a plainer word? And: you dislike the phenogram figure — what should the replacement
   show (e.g. the liability trajectory + the threshold cut, tips coloured by state)?
2. **Pagel > 2 correlated traits (§12.3.2).** `CorrelatedBinary` is 2 binary traits today.
   Extending to k binary traits is a real feature (2^k state space). Want it built?
3. **Sequence module scope (§12.5 / Ch14).** Move Pagel-style tree transforms *out* of the
   sequence module into traits, so the sequence module only simulates sequences? (I'll check
   what it currently exposes and propose.)
4. **Ordered-genome `extension` (Ch10).** Same ugly knob as the nucleotide model — reparametrize
   to a **mean-length** knob there too, for consistency with the INDEL decision?
5. **genes+intergenes as default (§11.5).** My plan: fold it in (nucleotide model is genic when
   you supply gene coords, structural otherwise) — confirm?
6. **§12.6 DEC figure floats into §12.7** in the PDF (a LaTeX float-placement drift, not a
   markdown ordering bug). Want figures pinned in place globally?
7. Quick confirms I'll otherwise default on: homologous-replacement **orientation check** (add
   it), **GFF multi-contig** (use the most-annotated contig + warn).

Autonomous progress log: Ch10 done · Ch13 prose pass done · Ch12 overview table done ·
atoms→blocks rename running · INDELs next · Ch12 figures + namespace conversion pending a
zombi2-free window.

## Figure house style
- [x] Style agreed and written to `figures/STYLE.md`; 6 reference figures done.
- [~] Apply the locked style to **all remaining figures** (regenerate → copy SVGs to
  `docs/img/` → rebuild manual). Existing figures only; new figures are per-chapter below.
- [?] ClaDS colour: kept grayscale (Adrián: "the rest is perfect").

## Cross-cutting / design questions
- [x] **Species-tree chapters split** — Ch4 "Species trees (basic models)", Ch5
  "Species trees (advanced models)" (titles done; content already divided).
- [x] **API namespace** — DONE (1113 tests). Submodules `zombi2.species/genomes/traits/
  sequences/coevolve/abc/distributions` + top-level re-exports kept (back-compat). Ch4/Ch5
  examples converted to `from zombi2.species import ...`. Later chapters' examples still to convert
  (do per-chapter). Docs `reference/api.md` may want per-namespace sections eventually.
- [x] **"Atoms" → "blocks"** — DECIDED (Ch11 §11.2). Implement across code+docs+manual later.
- [?] **Sequence module scope** — move tree-transformation (Pagel-style) OUT of the
  sequence module into traits (Ch12 §12.5), so the sequence module ONLY simulates sequences.
- [?] **Traits on gene trees** contradicts the Ch1 levels figure (traits shown only on the
  species level). Reconcile the figure/framing.
- [ ] **Genes+intergenes default** (Ch11 §11.5) — should the nucleotide model be
  genes+intergenes by default/exclusively, dropping the separate section?

## Code features (agent-sized)
- [ ] **INDELs in the nucleotide model (Ch11)** — indels ONLY in intergene regions; must
  not break existing behaviour; test thoroughly (keep the E. coli test). DESIGN (settling):
  - Length parameterization: **CONFIRMED** — replace `--extension` (ugly 0.999) with a
    **mean event length** knob, **geometric under the hood** (`extension = 1 - 1/mean`).
  - Span vs gene content: **CONFIRMED** (Adrián: "you got it right") — measure the affected
    length in **FULL coordinates so genes in the span count** (a region enclosing 100 genes is
    "bigger" than one enclosing 1 gene, and is hit differently). Snap breakpoints to intergene
    boundaries (genes never split); whole genes inside the span are carried by the event.
    VERIFY whether the current engine measures full vs collapsed (rust `draw_length` clamps to
    available) and fix if collapsed.
  - Loss barrier: add a **minimum-genome-size floor** so deletions can't empty the genome.
- [ ] **Pagel correlated traits > 2** (Ch12 §12.3.2) — support >2 correlated binary traits
  if not already; else implement.
- [ ] **Homologous replacement transfer** (Ch11) — currently checks flanking genes; also
  check that orientation coincides?
- [ ] **GFF multiple contigs** (Ch11) — currently one seqid; decide behaviour.
- [ ] (carried over) banner-only-on-help · forward `species_tree_extant.nwk` · `u` vs `e`
  leaf naming — reverted earlier; redo as a deliberate pass.

## Ch4 — Species trees
- [x] Rename chapter "Species trees" (drop "birth-death").
- [x] "Every ZOMBI simulation begins with a species tree" → "every simulation has a species tree".

## Ch5 — (basic/advanced split TBD)
- Loved it; clear. (API-surface discussion above.)

## Ch7 — Gene-family evolution  (DONE this round)
- [x] Reopened with the **three genome-evolution models** (independent / ordered / nucleotide)
  + new figure `docs/img/genome_models.svg` (3-row B&W schematic).
- [x] §7.2.3 removed the "ZOMBI-1 style of heterogeneity" sentence.
- [x] §7.2 closed with a **rate-models overview table**.
- [x] Examples converted to `from zombi2.genomes import ...` / `from zombi2.distributions import ...`.

### Ch7 RE-SPLIT (2026-07-05, Adrián's correction) — DONE, verified, built
The clarity refactor wrongly fused the genome OVERVIEW into the unordered chapter. Re-split per
Adrián: a generalities overview chapter, then per-level chapters. **New numbering (everything below
the genome block shifted +1):**
- **Ch7 `07-genome-evolution-overview.md` "Genome evolution"** — generalities: 3 levels, the
  level×rate two-axes framing, levels table + rate-model table, and a new **"Events, by level"**
  section (Adrián's key point — explicit per-level events + the "each level adds to the one below"
  nesting, with a summary table). Fixed accuracy: **`--genome-model` exposes only `unordered` /
  `nucleotide`; ordered is Python-only** (`OrderedGenome` factory) — the levels table now says so.
- **Ch8 `08-unordered-genomes.md` "Unordered genomes"** — the four events (heading renamed
  **"The four unordered events"**), rates (shared/per-genome/family), transfers, growth. Rate-model
  table moved up to Ch7; Rates section back-references it.
- **Ch9** gene-trees-and-output (was 8) · **Ch10** coupling (was 9) · **Ch11** ordered (was 10) ·
  **Ch12** nucleotide (was 11) · **Ch13** trait (was 12) · **Ch14** coevolution (was 13) ·
  **Ch15** sequences (was 14).
- All inter-chapter refs updated; level naming standardized to **unordered/ordered/nucleotide**
  (dropped "independent gene families" in prose AND in the `genome_models.svg` row-1 label);
  Ch3 roadmap "uniform versus per-family" → "shared versus per-family" (post-rename). Full manual
  rebuilds = **87pp**. Level naming, output-chapter placement (after Unordered, before Coupling)
  both **confirmed by Adrián** via AskUserQuestion. NOT pushed — holding for Adrián's review.

## Ch9 — Coupling
- [?] Feedback deferred (Adrián doesn't understand the model yet).

## Ch10 — Ordered genomes
- [ ] Segment events: explain how the **number of genes affected** is chosen (geometric),
  parametrized. New figure: **PMF of 3 geometric distributions** (different p) with means marked.
- [ ] "How events reach the genome section" — needs the 3-model classification up front
  (1 unordered, 2 ordered).

## Ch11 — Nucleotide genomes
- [ ] INDELs (see code features).
- [?] §11.2 rename "atoms".
- [?] §11.5 genes+intergenes default vs separate section.
- [ ] Homologous replacement: orientation check.
- [ ] GFF multi-contig behaviour.

## Ch12 — Trait evolution (lots of figures)
Section §-numbers below are the OLD numbering (pre-restructure); mapped to live sections.
- [x] continuous vs discrete: chapter-opener two-panel figure (`fig_trait_overlay.py` →
  `trait_overlay`). At "The shared overlay engine". Full-width. NOT committed.
- [x] §12.2 overview table of models (committed: "Trait models at a glance").
- [x] **Pagel figure DONE** — `fig_trait_pagel.py` → `docs/img/trait_pagel.svg`, two-panel
  (state square with rate-width arrows + a realization on a tree). Adrián APPROVED; folds in
  BOTH the "Pagel matrix" and the "two coevolved traits" asks (he confirmed one figure is fine).
  Full-width in the "Correlated binary characters (Pagel)" section. NOT committed yet.
- [x] **Threshold figure DONE** — replaced the disliked phenogram: rewrote `fig_trait_threshold.py`
  as an OU-style painted tree (branches = liability viridis, threshold marked on the colour bar,
  tip chips = discrete 0/1). Adrián approved the proposal. Caption rewritten, full-width. NOT committed.
- [x] "liability" term: KEEP with a gloss (Adrián's call). Gloss still to add to the prose.
- [x] Mk: two-panel (`fig_trait_mk.py` rewritten) — 3-state chain diagram + grayscale stochastic
  map on a tree. 3 categorical states = 3 greys (dark/mid/light). Full-width. NOT committed.
- [x] corHMM / hidden rate classes (`fig_trait_hiddenmk.py` rewritten): Pagel-style two-panel —
  4-state (obs x hidden) width-coded diagram + realization (branches by hidden class, observed-
  change marks, obs tip chips). Full-width. NOT committed.
- [x] DEC (`fig_dec.py` restyled to house B&W): centered title, no subtitle, FS_* fonts, grey area
  cells, stacked legend. Full-width. NOT committed.
- **ALL 6 Ch12 figures DONE** (pagel, threshold, overlay, mk, hiddenmk, dec) + wired into
  ch12 + full manual rebuilt (78pp). Batch built in the locked style; Adrián to review the whole
  Ch12 PDF at once. Nothing committed.
- [ ] code-block overflows: re-verify against the current PDF (earlier overflow pass may cover).
- [?] §12.5 move tree-transformation out of sequence module — Adrián: DON'T move yet.
- [x] misplaced float (§12.6 fig under §12.7): fixed globally by the figure-float pin in preamble.tex.
- [x] Pagel >2 correlated traits (CODE): `CorrelatedBinaryK` DONE by agent on branch
  `feature/pagel-k-traits` (1152 tests). NOT merged — awaiting Adrián. Manual snippet drafted,
  held out of the manual until merge (every-snippet-must-run rule).

## Ch13 — Coevolution
- [ ] §13.1 code block overflow.
- [x] §13.1 remove the confusing note.
- [ ] general figure of the different modes + accompanying table.
- [ ] §13.2 rename section titles as directed edges: "Traits → Genes", "Genes → Traits", etc.
- [x] §13.2 remove the GOE self-citation.
- More feedback coming.

## Ch14 — Sequences
- (see §12.5 design: sequence module = simulate sequences only.)

---

# ROUND 4 FEEDBACK (2026-07-05, from Adrián reading the current full PDF)

Figure numbers below are the FULL-manual numbering he saw. Status: [ ] todo · [x] done.
**Cross-cutting:** many figures are too tall/narrow — need to be WIDER. Committed checkpoint = `be39ced`.

**LANDED (2026-07-05, 3 parallel figure agents + my prose lane; NOT yet committed):**
Ch12 all figure fixes + 3 new figures · Ch5 figs 5.1/5.3/5.4/5.5/5.6/5.7 widened (5.1 → ~210 tips)
+ "fast diversification" label + drop "ClaDS" from §5.3 title + 5.5/5.6 code overflows fixed ·
Ch6 fig 6.1 widened + §6.2 overflows · Ch7 fig 7.1 nucleotide-row gene letters (edited the hand-SVG
`docs/img/genome_models.svg` — no generator exists) · Ch10 fig 10.2 widened · Ch9 removed the dense
9.1 4-panel figure, split into 4 standalone figs (`potts_genome`/`potts_coupling`/`potts_lossrate`/
`potts_retention`) placed across The model / Building / Running sections · Ch3 & Ch8 "reconciled"→
"complete" gene tree (figure title + both captions) · Ch4 fig 4.1 caption (ten not twenty tips),
removed v1-age note + backward-mode tip box, extinct tips `x*`→`e*` in fig_species_tree_extinct.
Ch4 **age-vs-crown** two-panel figure DONE (`fig_age_crown.py` → `age_crown`, placed after the
age_type explanation). Full manual rebuilt = 83pp. STILL TODO: Ch13 rewrite, Ch14 clock code,
species-output `e`/`u` naming code pass.

## Ch3
- [ ] **Fig 3.2 caption is WRONG** — it shows a **complete** gene tree, NOT a *reconciled*
  gene tree. Fix the terminology (complete vs reconciled vs pruned). Same error recurs in Ch8.

## Ch4
- [ ] Section 4.2 deserves a NEW figure: the meaning of **age** vs **crown** when you set the
  age (imagine two side-by-side panels).
- [ ] **Fig 4.1 caption wrong** — the tree does NOT have 20 tips; fix the caption to match.
- [ ] Remove the note about **version 1 requiring an explicit age**.
- [ ] **Fig 4.2: dead lineages must be `e1, e2, …` not `x1, x2`** (extinct = e). Ties to the
  reverted `u` vs `e` leaf-naming code pass — do it for real (extinct `e<n>`, unsampled `u<n>`).
- [ ] Remove the tip box "Use the default backward mode when you want a clean reconstructed phylogeny."

## Ch5
- [ ] **Fig 5.1**: good, but WIDER (too tall/narrow) and much denser tree (~200 tips).
- [ ] **Fig 5.3**: much wider (full page width).
- [ ] Section 5.3: **remove "ClaDS" from the section title**.
- [ ] **Fig 5.4**: wider.
- [ ] **Fig 5.5**: wider.
- [ ] Section 5.5: **code-block overflow**.
- [ ] **Fig 5.6**: wider; rename the annotation "this clade diversifies now fast" → **"fast diversification"**.
- [ ] Section 5.6: **multiple code overflows**.
- [ ] **Fig 5.7**: too tall/narrow → wider.

## Ch6
- [ ] **Fig 6.1**: wider.
- [ ] Section 6.2: **code overflows**.

## Ch7
- [ ] **Fig 7.1**: great, but the **genes in the nucleotide-genome row should carry letters
  (A, B, C, D…)** like the other two rows.

## Ch8
- [ ] A figure is captioned/presented as a **reconciled gene tree but it is NOT** (same terminology
  bug as Fig 3.2). Fix.

## Ch9 — Coupling
- [ ] **Fig 9.1: remove it** — too small, barely readable, overkill.
- [ ] **Cut its panels into independent figures** and use them to illustrate the model in the
  corresponding sections. (Adrián now understands the coupling model.)

## Ch10
- [ ] **Fig 10.2** (genes-affected-per-event): very good — just WIDER.

## Ch12 (my new figures — feedback)  [DONE items landed this round; not yet committed]
- [x] **Fig 12.2 (BM)**: rewrote `fig_trait_bm.py` to house style (wide, FS_*, centered title,
  colour bar) + seed search for an even value spread → good contrast, chips match branches. Full-width.
- [x] **Fig 12.3 (OU)**: made full-width (already house-style). (Tips cluster at the optimum by
  nature — that IS the OU signal; left as-is otherwise.)
- [x] Section 12.3.3 **Early burst** figure: modernized `fig_trait_earlyburst.py` (was broken by the
  color_bar signature change) → house-style painted tree + rate-decay strip. Full-width.
- [x] **Fig 12.4 (multivariate)**: rewrote `fig_trait_multivariate.py` — B&W, larger, full-width,
  ASCII "rho", centered title.
- [x] **Fig 12.5 (Mk)**: removed the fine tip→chip lines; new tree seed (no polytomy-like branch).
- [x] **Fig 12.6 (Pagel)**: removed the fine tip→chip lines.
- [x] **Fig 12.7 (corHMM)**: removed the fine tip→chip lines.
- [x] **Fig 12.8 (threshold)**: state legend moved to the bottom-left.
- [x] **Multi-optimum OU** figure: new `fig_trait_multioptimum.py` — painted tree + regime tip bars;
  regime-A lineages pull to the low optimum, regime-B to the high. Full-width.
- [x] **DEC tree**: new `fig_dec_tree.py` — a DEC history on a tree (root range A, dispersal/extinction
  event marks on branches, tip range cells). Added after the DEC schematic. Full-width.
- **Ch12 FIGURES COMPLETE** (12 figures, all house style). All held uncommitted for Adrián's review.

## Ch13 — Coevolution  [REWRITTEN 2026-07-05]
- [x] Full rewrite (Adrián: 3 node-pair sections, equal depth). Kept the strong intro (pipeline-vs-
  coupled, six-edge table, "arrows into S" rule) + NEW `coevolve_modes` figure (S/T/G triangle, into-S
  edges heavy). Sections: 13.2 Species&Traits (traits:species SSE / species:traits cladogenetic /
  both=ClaSSE), 13.3 Species&Genes (genes:species key-innovation / species:genes punctuational),
  13.4 Traits&Genes (traits:genes trait-linked / genes:traits gene-conditioned), 13.5 inference caveat.
  ALL SIX edges documented at equal depth (model → Python → CLI → what-it-recovers). Every snippet
  agent-verified + smoke-tested (all 6 drivers run). Full manual = 85pp.

## Ch14 — Sequences (CODE + design)
- [ ] Add the **molecular-clock models from PhyloBayes** for phylogram creation: **white noise,
  log-normal, gamma multipliers** (+ any others). These relax the clock (chronogram → phylogram).
- [?] DESIGN: should the clock models live in the **trait-evolution module** (Adrián leans yes),
  perhaps as a **separate chapter**? Note: today this is `rate_variation.py` (see [[zombi2-rate-modules]]
  which said "leave alone" — now SUPERSEDED: Adrián wants it expanded).
