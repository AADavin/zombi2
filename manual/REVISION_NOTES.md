# Manual revision — working backlog

Adrián's feedback, organized so nothing is lost. Status: [ ] todo · [~] in progress ·
[x] done · [?] needs discussion/decision.

## 🧭 Readability & clarity revision (2026-07-15, ratified with Adrián)

A structural pass to improve readability, aligned to the programmed **rates** (opportunity /
"how many clocks") and **coevolution** (`--all`) work. Ratified decisions:

1. **Book grouped into Parts I–VI** (matches this README's plan) via raw-LaTeX `\part{}` at each
   part opener. Chapters keep their numbers, so every "Chapter N" prose ref stays valid.
2. **Ghost lineages demoted** from a chapter to a `## Ghost lineages` section inside Ch5 (advanced
   species trees) — its natural home (complete-vs-reconstructed). `06-ghost-lineages.md` removed;
   no inbound links existed.
3. **Rates taught once**: new short anchor `06-how-rates-work.md` ("how many clocks, how fast" —
   `base × opportunity × modifiers`), owning the rate-model detail. Ch7 overview keeps levels +
   events-by-level and points to Ch6; Ch8/Ch5/Ch15 reference it. Speaks `per=`/opportunity ahead of
   Part 3 (the unified `--rate-per`/`per=` knob), so it never teaches a soon-deprecated class name.
4. **Manual teaches the concept; the online guide keeps the exhaustive catalog** — esp. Ch13
   coevolution (trim to one worked example per pair; `--all` = one section when it lands).

**Pass 2 (Adrián's structural feedback, ratified via AskUserQuestion):** rates does not belong in
Part III after species trees; add a global-overview chapter to Part I; merge install + first-sim;
split the massive coevolution chapter. Decisions: the rate concept → a **section of the new Ch2
tour**; coevolution → **3 chapters**. FINAL TOC (16 chapters + App A):

- **I Getting started** — 01 Introduction · 02 **A tour of ZOMBI2** (four levels · pipeline vs
  coevolution · **how rates work** §2.3 · ZOMBI2 vocabulary) · 03 **Getting started** (install +
  first-sim MERGED)
- **II Species trees** — 04 basic · 05 advanced (+ §Ghost lineages)
- **III Genomes** — 06 Genome evolution (overview) · 07 unordered · 08 gene-trees/output · 09
  ordered · 10 nucleotide
- **IV Traits** — 11 Trait evolution
- **V Coevolution** — 12 **the framework** (+ `--all` roadmap note) · 13 **State-dependent
  diversification** (species–traits: SSE, cladogenetic, ClaSSE) · 14 **Coupling gene content**
  (species–genes + traits–genes + **Null models** §14.3)
- **VI Sequences** — 15 sequence evolution · 16 Relaxed molecular clocks
- **App A** Gillespie · ILS supplement

DONE + verified (uncommitted, held for review): full renumber via `git mv`; new
`02-a-tour-of-zombi2.md` + `03-getting-started.md`; Ch1 levels section removed (→ tour); old
`06-how-rates-work` folded into the tour; coevolution split `13`→`{12,13,14}`; all 6 Part dividers
placed; **every hard "Chapter N" prose ref converted to a stable heading-id link** (restructure-proof;
the clocks heading got an explicit `{#molecular-clocks}` id); `↔`→en-dash (lmroman lacks the glyph).
`make manual` clean (no undefined refs / missing glyphs), **128 pp**; all snippets run (including the
redistributed coevolution ones). STILL PENDING: [ ] vocabulary/`per=` sweep across the genome
chapters (the concept now lives in the tour) · [ ] Ch16 clocks reframed as a *rate modifier* tied to
the tour · [ ] Adrián reviews, then commit.

**Pass 3 (discussion, 2026-07-15) — further decisions, NOT yet executed:**
- **Sequences go BEFORE coevolution.** Parts become IV Traits → **V Sequences** → **VI Coevolution**
  (coevolution is the capstone that couples all four levels, so all four must be introduced first).
- **Coevolution is organised by what is TARGETED** (Adrián's call): framework, then chapters grouped
  by the arrow's target — *shaping the tree* (into-S: SSE, key innovation, the joints — these grow the
  tree), *shaping traits & gene content* (overlays), *shaping sequences* (the new node). This axis
  scales to the 4th node; the earlier pair-based split (12/13/14) will be re-cut this way. Exact chapter
  count held until the sequence-target content lands.
- **Sequence coevolution: sequences are a TARGET-only node** (a driver would be harder — deferred),
  and it **lives in the coevolution part** as the 4th node. An agent is building this now; the specific
  sequence-target edges/labels are TBD.
- **Traits stays ONE chapter** (no split).
- **Figures drafted** (`figures/scripts/`, rendered to `figures/<name>/`, NOT yet copied to `docs/img/`
  or wired into chapters — drafts for review): (1) `fig_coevolve_modes4.py` → the framework figure as a
  **4-node graph** (S/T/G triangle unchanged + sequences Q as a dashed target-only node); the two
  Q-edge labels (`selection`, `regime`) are **placeholders** pending the real model. (2)
  `fig_rate_clocks.py` → the tour's **"how many clocks"** figure: the three opportunities as clock
  schematics + a real-data exponential-vs-linear LTT (BirthDeath vs SharedBirthDeath, mean over 40
  seeds). NB figure scripts that import `zombi2` must `sys.path.insert(0, repo_root)` — the editable
  install points at the *main* checkout, which lacks `SharedBirthDeath`.

**EXECUTED for the review PDF (2026-07-15):** (1) **Sequences reordered before Coevolution** — Part V
= Sequence evolution (12 sequence, 13 clocks), Part VI = Coevolution (14 framework, 15 SDD, 16
coupling gene content); clean rotation via git mv/mv, all cross-refs survive (heading-id links). (2)
**`rate_clocks` figure wired into the tour** (`docs/img/rate_clocks.svg`, Figure 2.2 in §2.3 "How
rates work"), sits full-width and clean. `make manual` green, **128 pp**, `build/zombi2-manual.pdf`
ready for Adrián to check. STILL DEFERRED (needs the sequence-target model first): the **coevolution
by-target re-cut** (framework / shapes-the-tree / shapes-traits-and-genes / [future] shapes-sequences)
— chapters 14–16 are still the by-pair cut for now — and wiring the **4-node `coevolve_modes4`**
figure into the framework.

**FOUND the sequence-coevolution model (branch `claude/coevolve-grammar-migration`, design
`docs/design/coevolve-grammar.md`, 2026-07-15):** the whole coevolve subsystem is being reframed onto
ONE grammar — every edge is `driver → target-variable : response` — over the 4-node **diamond
{S, T, G, Σ}**. Sequences (Σ) ride *gene* trees, are **target-only**, and the **S–Σ edge is
FORBIDDEN** (Σ couples only to its tier neighbours T, G). Sequence edges: **T→Σ** = trait-driven
selection (ω = dN/dS) + substitution speed (`DriverClock`/`OmegaSelector`); **G→Σ** = post-duplication
relaxed selection (`GeneEventOmega`). Nulls become `response = 0` (uniform); structural names
(`traits:species`) become primary, literature names (SSE/ClaSSE) aliases. **This grammar IS Adrián's
"organise by what is targeted."** Sequence tier LANDED on the branch; the into-species reframe waits on
the `Rates(per=…)` rate rename merging to main.
- **`coevolve_modes4` figure FINALISED to the real model** (Σ node renders fine; T→Σ "selection
  (dN/dS)", G→Σ "relaxed selection", no S–Σ, subtitle explains why). It is the intended framework
  figure, to drop in **when the coevolution part is rewritten onto the grammar** — the right,
  non-throwaway home for the by-target re-cut + the sequence sections; best done as/after the grammar
  lands on main (every-snippet-runs).

**EXECUTED 2026-07-15 (later):** Framework (Ch14) reframed onto the **four-level diamond** + the
diamond figure `coevolve_modes4` wired in (Fig 14.1, p100); the tour's levels figure (Fig 2.1)
**replaced by a single-panel diamond** `fig_levels_diamond.py` → `docs/img/levels_diamond.svg` (four
levels + the pipeline arrows S→T, S→G, G→Σ). PDF rebuilt, **129 pp**, clean.

**⚠️ `per=` opportunity knob LANDED ON MAIN (Part 3, phases A–E) — worktree is ~13 commits behind
`origin/main`.** New on main: `BirthDeath(per=…)` (A), `Rates(per="copy"|"lineage"|"shared")` (B/C),
`Per(unit, rate)` per-event mixing (D), the finest rung named `site` (E); preset classes
(`SharedBirthDeath`, `PerCopyRates`, `PerLineageRates`) are now deprecated shorthands. **This makes the
tour's "How rates work" `::: note` STALE** (it calls `per=` a *planned* refinement). To include it: (1)
**sync the worktree to `origin/main`** (commit the manual WIP on the branch, then merge — needed so the
new snippets run); (2) rewrite the tour note to teach `per=` as real; (3) migrate snippets
`SharedBirthDeath(…)` → `BirthDeath(…, per="shared")` (incl. the `fig_rate_clocks.py` figure) and
`--rate-per`/`PerCopyRates` → `per=` in the genome chapters; optionally show `Per(unit, rate)` mixing.
Also on main: standalone **Tools PDF** (PR #148, `manual/tools_to_chapters.py`). PENDING Adrián's OK to
sync.

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
   it). **GFF multi-contig** — DONE, better than the old "most-annotated + warn" plan: a
   multi-sequence GFF now seeds one chromosome per sequence (`read_gff_all` / `--gff`), feeding the
   chromosome tier (Ch11 §*Multiple chromosomes and the chromosome tier*).

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
- [x] **INDELs in the nucleotide model (Ch12)** — DONE + folded into PR #10 (2026-07-05).
  Intergenic insertion (`N`) / deletion (`E`) events, **Python-engine only** (Rust profiles path
  refuses them; the WIP's Rust impl was dropped per Adrián — ship Python-only). Insertion lays a
  run of novel nucleotides (fresh block) in an intergene; deletion removes a run from *within one*
  intergene, clamped so it never touches a gene and never drops below `MIN_GENOME_LENGTH` (floor).
  Own **`indel_mean_length`** geometric knob (default 10), independent of `extension` — this is the
  "mean event length" idea, applied to indels (the broader `--extension` reparametrization for
  STRUCTURAL events remains a separate deferred polish). CLI `--insertion/--deletion/
  --indel-mean-length`; setting either rate auto-routes to the Python engine. New
  `tests/test_nucleotide_indels.py` (40 tests: grow/shrink, gene-integrity invariant, floor,
  determinism, default-off, gate). 1169 tests total. Ported cleanly onto the renamed `SharedRates`
  (was `UniformRates` on the stale feature/manual-revision WIP). Manual Ch12 §"Intergenic indels" +
  Ch7 events-by-level updated.
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
  nesting, with a summary table). Fixed accuracy of the levels table. (At the time of this pass
  `--genome-resolution` exposed only `unordered` / `nucleotide`; `--genome-resolution ordered` was added as a
  CLI level on 2026-07-05 — see below — so ordered is no longer Python-only.)
- **Ch8 `08-unordered-genomes.md` "Unordered genomes"** — the four events (heading renamed
  **"The four unordered events"**), rates (shared/per-genome/family), transfers, growth. Rate-model
  table moved up to Ch7; Rates section back-references it.
- Current chapter order (the former **coupling** chapter was removed with the Potts purge and the
  files renumbered to close the gap): **Ch9** gene-trees-and-output · **Ch10** ordered ·
  **Ch11** nucleotide · **Ch12** trait · **Ch13** coevolution · **Ch14** sequences ·
  **Ch15** molecular clocks.
- **Ordered is now a CLI level too (2026-07-05):** `--genome-resolution ordered` added (was Python-only).
  Wires `SharedRates(inversion=,transposition=)` + `genome_factory=OrderedGenome(ids, extension=)`;
  reuses `--inversion/--transposition` (per gene copy for ordered, per nt for nucleotide) and
  `--extension` (in genes; None=single-gene). Those three flags now default `None` and resolve
  per-level. Rearrangements need `--rate-model shared` (PerGenomeRates carries none → clean error).
  Ordered always uses the Python engine (Rust counts-only/trace paths gated off). Ch7 overview
  de-flagged ("Python API only" removed); Ch11 got a "From the command line" subsection. 1113 tests.
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
9.1 4-panel figure, split into 4 standalone coupling figs placed across The model / Building /
Running sections · Ch3 & Ch8 "reconciled"→
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
  both=ClaSSE), 13.3 Species&Genes (genomes:species key-innovation / species:genomes punctuational),
  13.4 Traits&Genes (traits:genomes trait-linked / genomes:traits gene-conditioned), 13.5 inference caveat.
  ALL SIX edges documented at equal depth (model → Python → CLI → what-it-recovers). Every snippet
  agent-verified + smoke-tested (all 6 drivers run). Full manual = 85pp.

## Ch14 — Sequences (CODE + design)
- [ ] Add the **molecular-clock models from PhyloBayes** for phylogram creation: **white noise,
  log-normal, gamma multipliers** (+ any others). These relax the clock (chronogram → phylogram).
- [?] DESIGN: should the clock models live in the **trait-evolution module** (Adrián leans yes),
  perhaps as a **separate chapter**? Note: today this is `rate_variation.py` (see [[zombi2-rate-modules]]
  which said "leave alone" — now SUPERSEDED: Adrián wants it expanded).
