# Coverage audit — 134 comments vs the current design (2026-07-18)

Cross-check of every dashboard comment against the six design docs + SPEC + the chapter drafts, run as a
19-agent workflow (partition → per-chapter audit → synthesis). Purpose: detect holes before building.

## Headline

| Status | Count |
|---|---|
| **Addressed** by the current design | 109 |
| **Partial** (touched, real gap remains) | 18 |
| **Hole** (nothing covers it) | 7 |
| **Total** | 134 |

Only **three** holes are sharp enough to matter before building; two are trivial prose edits the unwritten
Ch1 will subsume; two are a figure-content cluster that fell through decision D10's scope.

## The 7 holes

| id | ch | ask | why it's a hole |
|---|---|---|---|
| **C3.3** | 3 | Stop implying the seeded examples reproduce the printed figures | Text says "20 tips" over a 10-tip figure + a latent **family-7 no-survivor bug**; an unmade needs-decision, not a redesign gap. **Sharp.** |
| **C5.12** | 5 | Move the identifiability note into a collecting file of advanced notes | No advanced-notes mechanism (`OBSERVATIONS.md`) exists; the Louca/Pennell non-identifiability caveat — the key one for rate-inference users — has nowhere to live. **Sharp.** |
| **C4.5** | 4 | Verify the Stadler 2009 citation on the stem-vs-crown sentence | Stem/crown returns via `age_type`, but the misplaced-citation fix (drop Stadler 2009) is recorded nowhere and will silently return. **Sharp.** |
| **C7.4** | 7 | Change Fig 7.1 family legend from coloured squares to line segments (a square already means "duplication") | D10 gates only Figs 2.1/2.2/14.1; this legend fix falls through the gap. *(figure cluster)* |
| **C11.7** | 11 | Add a tree panel to the multivariate figure, both correlated traits painted on the phylogeny | Same D10 scope gap; the tree-panel request is in no draft or design doc. *(figure cluster)* |
| **C1.8** | 1 | Delete the "shortest possible taste" filler sentence | No Ch1 draft yet pins the deletion. *(trivial; Ch1 will subsume)* |
| **C1.9** | 1 | Make the lead-in an invitation to run something, not a minimality claim | No Ch1 draft yet pins the rephrase. *(trivial; Ch1 will subsume)* |

## The 18 partials (touched, but a real gap remains)

- **C1.10** — symmetry met (`sequences.simulate_sequences` added), but the redesign removes the top-level `z.simulate_*` this comment *praised*; **D8 gating still open**.
- **C2.7** — "you choose which levels" carried by P(…) notation, but the three-example figure is unbuilt (and its colour premise flipped to B&W in SPEC §7).
- **C3.1** — false demo-notebook line would drop on rewrite, but **"ship demo notebooks?" is an unresolved product decision**; `coevolve-grammar.md:370,553` still assume a notebook vehicle.
- **C3.2** — move unblocked (D10) but committed in no artifact; entangled with the still-open C3.3.
- **C4.2** — redesign deliberately *reverses* the ask (forward/backward demoted to an engine detail); index l.42 still promises the old opening.
- **C4.6** — opening reversed as C4.2; the two-mode figure never commissioned; complete-vs-reconstructed section still to design.
- **C5.4** — "one general class, no zoo" answered (`simulate_species_tree`), but whether `prune` is `Tree.prune` is still to design; **D8 open**.
- **C6.8** — word fix + figure conventions decided, but nothing relabels/redraws `event_levels.svg`; Ch6 carries zero figures.
- **C6.11** — chromosome-network section delivered, but the requested new figure is only *named* in the index, not drawn.
- **C7.2** — gene-conversion section exists, but the four-events opening still forward-references it ("A fifth, optional event… below") — the exact clause to move.
- **C9.4** — word-axis settled (resolution canonical, `--genome-model` a fossil), Rates zoo dissolved, but the **`Genomes-result/GenomeResult` collision + Chromosome export unresolved**; result API TBD; **D8 open**.
- **C10.7** — mean-length-in-nucleotides half addressed, but **choosing the length *distribution* (as ZOMBI1 did) + per-event-type is silently dropped**; design keeps geometric-only, one global knob.
- **C11.3** — per-level result objects kept, but the **cross-level symmetric-family decision unresolved** (D8), the **sequence result stays unnamed**, and **no "what each `simulate_*` returns" table** exists.
- **C12.3** — "why a phylogram" resolved, but the draft never states head-on that topology/ancestry is untouched; drops the evolve-down-the-timetree equivalence remark.
- **C13.3** — table-first clocks planned, but no captioned/numbered house standard is codified; the Ch7 draft's bridge table is a bare grid.
- **C14.5** — the `driver→target:response` formula is deleted (ask moot), but Ch9's surviving coupling definitions are blockquotes, not display math.
- **C15.5** — death-driving stated for the trait/BiSSE case, but the gene-content section never says a driver *gene* can also raise extinction (the poisoned-chalice point) — survives only as a "Still to design" note.
- **C15.6** — "blocks of code, no clarity" fixed (prose-led), but the "and maths" half is dropped: no equation after the exp-link log-rate model was replaced by a Table mapping.

## What this ADDS to the pre-build decision checklist

Beyond the manual/figure fixes, the partials surface a few genuine **decisions** not yet on the checklist:

1. **Result objects** — the `Genomes-result/GenomeResult` name collision, the unnamed sequence result, and a
   "what each `simulate_*` returns" table (C9.4, C11.3). *(→ naming pass + one small table.)*
2. **Length distribution** — do we support choosing the indel/segment length *distribution* (ZOMBI1 did) and
   per-event-type, or keep geometric-only, one global knob? (C10.7) *(→ scope decision.)*
3. **Demo notebooks** — ship them or not; `coevolve-grammar.md` still assumes they exist (C3.1). *(→ product decision.)*
4. **D8 gating** — the top-level `z.simulate_*` convenience the redesign removes was explicitly praised by a
   user (C1.10); confirm we override that. *(→ already decided "namespaced only"; just confirm.)*

The three sharp holes (C3.3, C5.12, C4.5) are **manual/figure/bug** work, not API blockers — queued for the
chapter pass. C3.3's latent **family-7 no-survivor bug** is worth a separate code look.
