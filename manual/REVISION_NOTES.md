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
- [ ] §12.1 two-panel figure: continuous vs discrete models.
- [ ] §12.2 overview table of models.
- [ ] §12.2.3 needs a figure.
- [ ] §12.2.4 code block overflow.
- [ ] figure: two panels, two coevolved traits.
- [ ] §12.3.1 code block overflows.
- [ ] §12.4 figure: two panels (Markov chain + simulated trait on phylogeny).
- [ ] §12.3.2 figure of the Pagel matrix; support >2 correlated traits (code).
- [ ] §12.3.3 illustrative two-panel figure.
- [ ] §12.3.4 the term "liability" — is it Felsenstein's? dislikes it; dislikes the
  phenogram — want a different figure.
- [ ] §12.4 code blocks overflow.
- [?] §12.5 move tree-transformation out of sequence module (design above).
- [ ] §12.6 needs figures; a figure is MISPLACED (the §12.6 figure is under §12.7) — fix.

## Ch13 — Coevolution
- [ ] §13.1 code block overflow.
- [x] §13.1 remove the confusing note.
- [ ] general figure of the different modes + accompanying table.
- [ ] §13.2 rename section titles as directed edges: "Traits → Genes", "Genes → Traits", etc.
- [x] §13.2 remove the GOE self-citation.
- More feedback coming.

## Ch14 — Sequences
- (see §12.5 design: sequence module = simulate sequences only.)
