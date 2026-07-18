# Result objects — design target

**Status: the design to build.** What every `simulate_*` returns. Designed with Adrián on 2026-07-18.
The output half of the input grammar: as the inputs collapse to one `scope(base) × modifiers` grammar,
the outputs collapse to **one Result shape per level**. Parallels the level design docs; read `SPEC.md`
and them first.

---

## Principle: the result mirrors the input

Every `simulate_*` returns a **`<Level>Result` bundle** — never a bare object. The bundles share a common
spine and differ only in payload, exactly as the levels share the rate grammar and differ only in events.

**Confirmed:** `simulate_species_tree(...)` returns a **`SpeciesResult`** (not a bare `Tree`) — species
alone produces *two* trees plus the event log, which no single object can hold, and bundling keeps the
four levels symmetric.

## The common spine (on every result)

- `.events` — the event log (the **compact source of truth**; see memory model)
- `.tree` / `.complete` / `.reconstructed` — the tree(s) it ran on
- `.write(dir, [...])` — write the chosen outputs to disk (same **write** vocabulary as the CLI `--write`)
- `.seed`

## The bundles — names and payload (the "what each `simulate_*` returns" table)

Names are `<Level>Result` with the **plural level word**, which resolves the `GenomeResult`/`Genomes`
collision the coverage audit found (→ `GenomesResult`) and finally names the sequence result.

| Level | Returns | Payload accessors |
|---|---|---|
| Species | `SpeciesResult` | `.complete`, `.reconstructed`, `.fossils`, `.events` |
| Genomes | `GenomesResult` | `.genomes` (per-node content), `.gene_trees`, `.profiles` (sparse presence/count matrix), `.transfers`, `.events` (D/T/L/O); structured adds `.chromosome_network`, `.gene_order` (GFF/BED) |
| Sequences | `SequencesResult` | `.alignments`, `.ancestral`, `.events` |
| Traits | `TraitsResult` | `.values`, `.history`, `.events` |

## The memory model — three dials (declare what you want, pay only for that)

Richness must not explode memory even in-memory. Three dials, cheapest → richest:

1. **What's recorded** — `simulate_genomes(..., record=["profiles"])`.
   Declaring intent up front **scopes what the run computes and keeps.** Profiles-only makes the engine a
   *profile accumulator*: it tallies the sparse matrix as lineages evolve and **never builds the event log
   or the gene-tree objects at all.** Footprint = the matrix, streamable.
2. **What's derived** — if the event log *is* kept (the default), the rich views (`ancestral`,
   `gene_trees`, `profiles`) are **reconstructed lazily from the log on access**, never all-resident. Ask
   for one node's ancestral genome and it is replayed just-in-time; iterate and they stream one at a time.
3. **What's written** — `.write(dir, [...])` picks which of those hit disk, streamed node-by-node
   with bounded memory, however big the tree.

**The honest tradeoff:** narrow intent means the other views genuinely **are not there.** Accessing
`.gene_trees` on a profiles-only run raises `"not recorded — re-run with gene_trees"`, because deriving it
would need the event log that was not kept. Generality costs memory; the user decides how much.

**Default** = the friendly middle: keep the event log + trees (small, O(events)/O(tips)), derive the rest
on demand. **Scale** (millions of tips) = narrow `record=` and the footprint collapses to the sparse
matrix; the event log itself can be disk-backed so the in-memory bundle stays tiny.

Not new machinery: this is the sparse `ProfileMatrix` (COO) and the "event-trace + lazy replay" work
already in the codebase, promoted to *the* output model.

## Outputs catalogue

The full list of every file — which `--write` token writes what, in what format — is
**Appendix B** of the manual, the single cross-level reference (unblocked by this design).

## Still to design (naming only)

- **Decided: the record dial is `record=[...]`** (Python) — it scopes what the run computes and keeps.
  Asking to `write` an output implies recording it. (The concept was settled with the three dials.)
- **Decided (2026-07-19): disk output is `write` at every level** — the result method `.write(...)`, the
  CLI `--write`, and the existing code all use the one word; `include` is retired (it was a reactive dodge
  around a name collision — the collision should have unified on `write`, not invented a new word). The
  exact `.write()` argument form (positional list vs keyword) is a build detail.

## What to delete / change

- Rename to the `<Level>Result` scheme (plural level word); fix the `GenomeResult`/`Genomes` collision;
  name the currently-unnamed sequence result (`SequencesResult`).
- `simulate_species_tree` returns `SpeciesResult`, not a bare `Tree`.
- Rich accessors are lazy/streaming, gated by the record dial — not eager materialisation.
