# Build plan — turning the design into code

**Status: the executable plan.** The eight design docs (`species-api`, `genome-api`, `sequence-api`,
`trait-api`, `chromosome-network`, `coupling-api`, `result-api`, and `SPEC.md`) are the *concepts*. This
is how they become *code*, in a safe order. Principle: **concepts → code → chapter** — build a level,
then its chapter documents what was built. Agreed with Adrián on 2026-07-18.

---

## Strategy (the four agreed calls)

1. **Clean break, not a shim.** Delete the old zoos; target a major version. A thin deprecated alias for
   the old class names may live *one* release if cheap, but the design is a redesign, not an addition.
2. **Rewrite each level *and its tests* together**, in one PR. The old suite sometimes *encoded* old
   behaviour (and once a bug) — treat every surprising test change as a chance to catch that, not to
   preserve it.
3. **Redesign, not refactor → byte-identity is not a goal.** Old seeds are not expected to reproduce; a
   clean seed break is fine. (This is the opposite discipline from our refactor-era work, where
   byte-identity was sacred — do not confuse the two.)
4. **One level per PR, in dependency order** (below). Each PR is green (its rewritten tests + the chapter's
   worked examples run) before the next starts.

## Phase 0 — Foundations (build before any level)

Everything downstream imports these, so they come first:

- **`zombi2.modifiers` (`mod`)** — the shared modifier vocabulary: `OnTime`, `OnTotalDiversity`, `FromParent`
  (`spread`, `reverts_to`, `pull`), `ByFamily`, `ByLineage(dist=)`, `Markov`, `ByChromosomeSize`,
  and `DrivenBy` (Part III). Dimensionless, `*`-composable.
- **`zombi2.scope`** — the scope wrappers: `Global`, `PerCopy`, `PerLineage`, `PerSite`, `PerChromosome`.
  (`PerGenome` dropped — one genome per lineage.) These *wrap* a base; they do not multiply.
- **The `Rate`** — `scope(base) × modifiers` (`SPEC §5`); users never build one by hand.
- **The Result spine + event log** — `<Level>Result` base with `.events`, `.tree`, `.write(include=)`,
  `.seed`, and the **record / derive / write** memory model (`result-api.md`): the event log is the compact
  source of truth; rich views are derived lazily and streamed. Reuse the sparse `ProfileMatrix` (COO) and
  the event-trace lazy-replay work already in the tree.

## The level order (each = one PR)

| # | Level | Entry point → Result | Touches Rust? | Deletes |
|---|---|---|---|---|
| 1 | **Species** | `simulate_species_tree` → `SpeciesResult` | No (pure Python) | the 7 process classes |
| 2 | **Genomes — unordered** | `simulate_unordered` → `GenomesResult` | **Yes** (`zombi2_core`) | the `RateModel` hierarchy, `genome_factory`, `TransferModel` |
| 3 | **Sequences** | `simulate_sequences` → `SequencesResult` | reuses `exp(Qt)` | the 8-class `Clock` zoo |
| 4 | **Traits** | `simulate_continuous` / `simulate_discrete` → `TraitsResult` | No | the 13-class model zoo; `DEC` → experimental |
| 5 | **Genomes — structured** | `simulate_ordered` / `simulate_nucleotide` | **Yes** | — (extends #2) |
| 6 | **Coupling** | `mod.DrivenBy` + `joint.simulate` | via the driven level | the `coevolve` command |

**Why this order:** species is the foundation and pure Python (fastest, safest first); the shared grammar
it exercises (modifiers, scope, Result, record dial) is then proven for everyone. Genomes-unordered is next
because sequences and traits both consume its gene trees. Structured genomes (#5) is deferred behind
unordered because its chromosome-network machinery (re-mint, eNewick `#H`, per-gene-copy translocation,
length distributions) is the heaviest single piece. Coupling is last — it reaches into every level.

### Per-level definition of done

- The entry point returns the `<Level>Result` bundle; the old classes are gone.
- Rates use `scope(base) × modifiers`; the level's modifiers work and compose.
- The record dial scopes memory (`record=[...]`); rich views are lazy/streamed.
- The level's tests are rewritten to the new API and green; the chapter's worked examples run.
- The CLI command mirrors the Python API (`--write`/`record` selector included).

## Level-specific notes

- **Species (#1):** `birth`/`death = number × modifiers`; Yule = `death=0`; ClaDS = `FromParent`; skyline =
  `OnTime`; diversity-dependent = `OnTotalDiversity`. `mass_extinctions`/`sampling`(ρ)/`fossils`(v1 = side output,
  lineage not removed). Both `.complete_tree` + `.extant_tree`. Clade shift & ghost lineages are out (v1).
- **Genomes-unordered (#2):** D/T/L/O keyword rates; `ByFamily` (per-rate or family-wide slot); `TransferModel` → arguments
  (`transfer_to`/`replacement`/`self_transfer`; donor weight = modifier, recipient weight = mechanic).
  `.profiles` sparse + lazy. **Rust rebuild** (`maturin build --release`) after core changes.
- **Sequences (#3):** substitution model = a **menu** (`jc69`/`hky85`/`gtr`/`lg`/codon); clock = `ByLineage`/
  `FromParent`/`Markov` on the **species tree** (lineage clock), `ByFamily` = per-family speed, they
  compose; `+Γ` = `gamma=`. `.ancestral` from the recorded nodes.
- **Traits (#4):** BM native; OU = `reverts_to`+`pull`; EB = `× OnTime`; Mk = `switch=` (scalar / `{a->b}` /
  matrix); threshold via liability; correlated = one call + `correlation=` overlay. `DEC` → experimental;
  SSE leaves for #6.
- **Structured genomes (#5):** rearrangements (per gene copy, `inversion_length=Geometric(mean=)`);
  chromosomes (`scope.PerChromosome`; fission/fusion/loss; translocation **per gene copy**); the
  **chromosome network** — re-mint both children, eNewick `#H` + event table, recover-not-reconcile
  (Option A). Length distributions configurable per event type.
- **Coupling (#6):** `mod.DrivenBy(source, mapping)` — file source = conditioned (folds into the target
  level's command), live-level source = joint (`joint.simulate`, grows both). Gene-content source =
  presence (`"genes:toxin"`) or count (`"genes:count"`). Birth *and* death drivable. Nulls = recipes
  (drop / swap `ByLineage` / shuffle), no `.null()` API. Value-driving (OU optimum) deferred.

## Cross-cutting, alongside the levels

- **Appendix B (outputs catalogue)** — write as the Result objects land (which `--write`/`record` token
  produces which file, in what format).
- **The three audit issues** — fix **C3.3's latent family-7 no-survivor bug** during #1–#2 (it is a real
  code bug); C5.12 (identifiability note → an advanced-notes file) and C4.5 (citation) in the manual pass.
- **Part I chapters** (Ch1 Intro, Ch3 Getting started) — manual pass, parallel, not blocking.

## Risks

- **Rust parity (#2, #5).** The genome sim uses `zombi2_core`; keep the Python API a thin wrapper and
  rebuild + run parity checks after any `rust/` change.
- **The memory model is the one genuinely-new piece of infra** (event log as source of truth + lazy replay
  + streaming). It is the highest-risk build item; land it in Phase 0 and prove it on species before the
  data-heavy levels.
- **Seed breakage is expected** — audit any test that pins a seed and rewrite it, do not "restore" it.
