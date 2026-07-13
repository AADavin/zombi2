# Design note: gene-order export via an enriched nucleotide event trace

**Status:** Phase 1 + 1.5 + 2a implemented — `Region` + `dest` persisted (`--write geneorder`),
and `zombi2 tools export --format breakpoints` (adjacencies broken per tree edge); Phase 2b
(gff / posortho / ffgc / dupinfo) pending · **Author:** Adrián (with Claude) · **Scope:**
`zombi2 genomes` + `zombi2 tools export`

**Decisions locked:** native coordinates only (no fork-compat layer); `Region` persisted for
**both** the nucleotide and ordered models; the structural-event log is a **distinct file**
`Geneorder_events.tsv` (not an overloaded `Events_trace.tsv`), written by `--write geneorder`.

## 1. What we want

Krister Swenson's fork exposes a gene-order **export surface** — `zombiExporter` with targets
`gff`, `bed`, `posortho` (positional orthologs), `ffgc`, `dupinfo` (family → #duplications), and
`breakpoints` (adjacencies broken per tree edge). We want the same for zombi2, ideally as an
**a-posteriori `zombi2 tools export`** command rather than something bolted into every simulation
run.

## 2. Why it can't be purely a-posteriori *today*

The gene-order state is built in memory during a run — the `SegmentRegistry` (every segment's
ancestral provenance) plus each event's `Region` (physical breakpoint / arc) — and then largely
**discarded**. Concretely, on disk today:

- **`Events_trace.tsv` has no positional columns.** Its header is
  `time · event · branch · donor · recipient · family · parent · child1 · child2` — event *type*
  and gene-tree topology, but no breakpoints, arc lengths, or coordinates. (It is also a
  *counts-model* artifact; see §5.)
- **The nucleotide model doesn't even write that trace.** It writes only *derived views*:
  `profiles`, `trees` (block gene trees), `reconciliations` (block histories), `ancestral`
  (ancestral DNA), `bed` (per-node gene coordinates).
- **The `Region` isn't persisted.** Each event's physical coordinates live on the input
  `Selection.region` (`Region(chromosome, start, length, strand)` in `events.py`), but the logged
  `EventRecord` keeps only `event, branch, time, genes[], donor, recipient, insertion` — no
  `region` field. So the breakpoint is available at event time and thrown away.

Consequence: the exporter targets split. `gff` / `bed` / `ffgc` / `posortho` are recoverable from
what's already written (BED coords, block reconciliations = positional orthology). But
**`breakpoints`** (and clean per-branch `dupinfo`) need the raw structural-event coordinates, which
are the one thing not on disk. That is the gap this note closes.

## 3. Proposal

Two small, additive changes plus one new tool:

1. **Persist the `Region`.** ✅ *done.* Added an optional `region` field to `EventRecord` and pass
   `region=selection.region` at the structural-event log site in
   `genome_sim.py::GenomeSimulator` (the `for group in genome.apply(...)` loop — `selection` is
   already in scope). The nucleotide model runs through this *unchanged* shared simulator, so this
   one change captures inversion / transposition / loss / duplication coordinates for it **and**
   for the ordered-genome model (which also fills `Region`) — verified by
   `tests/test_geneorder_events_output.py`.
2. **Write a distinct `Geneorder_events.tsv`.** ✅ *done.* A new `--write geneorder` target on the
   nucleotide model, serialised by `simulation.py::geneorder_events_from_log` (§4). A **separate
   file** rather than an overloaded `Events_trace.tsv`, because the counts model's trace is a
   *family*-event log consumed by `zombi2 sequence`; reusing the name for a *structural*-event log
   would be a same-name/two-schemas hazard. It sits alongside the existing `Karyotype_trace.tsv`,
   and — the payoff of a distinct file — the counts trace, its Rust writer, and `read_events_trace`
   are **untouched**.
3. **Add `zombi2 tools export`** — ⏳ *Phase 2.* Reads `Geneorder_events.tsv` + the root gene
   layout, **replays** the events to reconstruct per-node gene orders and block identities, and
   emits any target format (§6). The replay is the logic already exercised by
   `tests/test_geneorder_examples.py`, generalized from a hardcoded event list to reading rows.

## 4. `Geneorder_events.tsv` schema

One row per structural event (speciation / leaf markers filtered out), across every branch —
filter on `branch` for the per-branch view. Explicit typed columns rather than a packed string,
because zombi2 has a clean `Region`:

```
time  event  branch  family  chrom  start  length  strand  dest  donor  recipient
```

`chrom / start / length / strand` are the physical arc the event acted on; `dest` is the
paste / insert position where one applies. All empty for events with no region. Coordinates are
zombi2-native: half-open `[start, start+length)`, circular (wrapping arcs allowed).

Per-event coverage **as implemented (Phase 1 + 1.5)**:

| event (`.value`) | arc (chrom/start/length) | dest | fully described? |
|---|---|---|---|
| inversion `I` | ✓ | — | **yes** — the arc *is* the breakpoint |
| loss `L` | ✓ | — | **yes** — the removed arc |
| transposition `P` | ✓ source arc | ✓ paste position | **yes** |
| duplication `D` | ✓ source arc | — (tandem) | **yes** — the copy is tandem (immediately after the arc) |
| origination `O` | ✓ insert at `start`, gene of `length` | — | **yes** (novel gene; the root *seed* has no single position → empty) |
| transfer `T` | ✓ donor arc | ✓ recipient insert (`donor`/`recipient` cols name the branches) | **yes** |
| translocation `X` | ✓ source arc | — | partial — cross-chromosome dest needs a `dest_chrom` column (Phase 1.6) |

Only translocation (a multi-chromosome event, which also has `Karyotype_trace.tsv`) is not yet
fully described; every single-chromosome gene-order event is.

## 5. Why a distinct file (not an overloaded `Events_trace.tsv`)

`Events_trace.tsv` is a **counts-model** artifact: a *family*-event log (`O/D/T/L` + gene-tree
topology) with two writers (Python + a Rust fast path) and a strict header, consumed by
`zombi2 sequence`. `Geneorder_events.tsv` is a *structural*-event log for the nucleotide/ordered
models. Keeping them separate means:

- **Rust is never touched** (the new writer is Python-only; the nucleotide model has no Rust path
  for rearrangements anyway).
- No same-name/two-schemas ambiguity, and `read_events_trace` / `zombi2 sequence` keep working
  unchanged.
- It sits naturally beside the existing per-model `Karyotype_trace.tsv`.

## 6. `zombi2 tools export`

```
zombi2 tools export GENOMES_DIR --format breakpoints [-o OUT]     # phase 2a (implemented)
```

**Key simplification found while building:** a full event *replay* engine turned out to be
unnecessary for the order-based formats, because **`--write bed` already writes the reconstructed
per-node gene order + orientation** (`BED/<node>.bed`). So `zombi2 tools export` derives the
order-based formats by *reading and comparing those files*, not by re-simulating. Inputs by format:

- **`breakpoints`** ✅ — the adjacencies broken on each tree edge. Read each node's gene order from
  `BED/<node>.bed`, express a genome as its set of circular signed-gene adjacencies, and for each
  edge report `parent_adjacencies − child_adjacencies`. **Exact for content-conserving
  rearrangements** (inversion / transposition). With duplication/loss the two genomes differ in
  gene *content* (and gene names repeat, so extremity labels are no longer unique) — a differing
  adjacency then reflects gained/lost content, not a broken adjacency; those edges are approximate.
  `zombi2/tools/geneorder_export.py::breakpoints_tsv`.
- **`gff` / `posortho` / `ffgc`** ⏳ — also from the per-node BED orders (gff = reshape; posortho =
  group leaf genes by id; ffgc = leaf orders + `ancestral` FASTA).
- **`dupinfo`** ⏳ — NOT derivable cleanly from `Geneorder_events.tsv`: a nucleotide duplication is
  a *segment* event logged under its source, not a per-gene count. A gene-level dupinfo should come
  from the **block gene trees** (`--write trees`), which already reconcile per-block duplications.
- **`breakpoints` (coordinate-level)** — `Geneorder_events.tsv` remains the source for the raw
  rearrangement *coordinates* per branch, if a coordinate (rather than adjacency) view is wanted.

## 7. Interop with the fork's `zombiExporter` (native-only)

The fork's `_geneorderevents.tsv` uses `TIME, EVENT, BREAKPOINTS, LENGTH`, event codes
`INV/POS/TDUP/DUP/LOSS/ORIG/LFER(_F/_B)/AFER`, and **0-based inclusive** `"start,end"` breakpoint
strings. zombi2's encoding is deliberately native (half-open `[s, s+ℓ)`, single-letter
`EventType` codes); a **`--fork-compat` layer was considered and dropped** — `zombi2 tools export`
is the intended consumer, and a mapping to Krister's exact encoding can be added later if direct
`zombiExporter` interop is ever wanted.

## 8. Phasing

- **Phase 1** ✅ *done* — `Region → EventRecord` (both models) + `--write geneorder` →
  `Geneorder_events.tsv`. Tests in `tests/test_geneorder_events_output.py`.
- **Phase 1.5** ✅ *done* — added the `dest` column and captured the paste position (`P`), recipient
  insert (`T`) and novel-gene position (`O`); `D` is tandem by construction. Every
  single-chromosome event is now fully described. (Transposition `dest` is surfaced via a
  `self._event_region` stash set inside `apply` — zero new RNG draws, so no existing output moved.)
- **Phase 1.6** (small) — a `dest_chrom` column for translocation (`X`), the only remaining gap.
- **Phase 2a** ✅ *done* — `zombi2 tools export … --format breakpoints`: broken adjacencies per
  tree edge, from the per-node BED gene orders (no replay needed — see §6). Tests in
  `tests/test_geneorder_export.py`.
- **Phase 2b** — `gff` / `posortho` / `ffgc` (from the per-node BED orders) and a gene-level
  `dupinfo` (from the block gene trees).

## 9. Decisions & remaining questions

**Resolved:** coordinate convention = native only (§7); model scope = both nucleotide and ordered;
filename = distinct `Geneorder_events.tsv` (§5).

**Still open:**
1. **Replay-only vs. reuse derived views** — should `tools export` always replay from
   `Geneorder_events.tsv` (one canonical path), or short-circuit to the already-written
   `bed`/`reconciliations` where present?
2. **Priority order of export targets** — Krister flagged `breakpoints` and `posortho` as
   load-bearing for gene-order studies; confirm before building the rest.
3. **`dest` now or with the export tool** — do Phase 1.5 (`dest` columns) before or alongside the
   Phase 2 replay engine?
```
