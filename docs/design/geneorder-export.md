# Design note: gene-order export via an enriched nucleotide event trace

**Status:** Phase 1 implemented (`Region` persisted + `--write geneorder`); Phase 2 (`tools
export` replay) pending · **Author:** Adrián (with Claude) · **Scope:** `zombi2 genomes` +
a new `zombi2 tools export`

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
time  event  branch  family  chrom  start  length  strand  donor  recipient
```

`chrom / start / length / strand` are the physical arc the event acted on, empty for events with
no region. Coordinates are zombi2-native: half-open `[start, start+length)`, circular (wrapping
arcs allowed).

Per-event coverage **as implemented in Phase 1** (region = the acting `Selection.region`):

| event (`.value`) | chrom/start/length | strand | fully described? |
|---|---|---|---|
| inversion `I` | ✓ arc | 1 | **yes** — the arc *is* the breakpoint |
| loss `L` | ✓ arc | 1 | **yes** — the removed arc |
| transposition `P` | ✓ *source* arc | 1 | partial — paste **dest not yet logged** (⏳) |
| duplication `D` | ✓ *source* arc | 1 | partial — insert **dest not yet logged** (⏳) |
| origination `O` | — (empty) | — | ⏳ — separate log path, region not threaded yet |
| transfer `T` | — (empty) | — | ⏳ — separate log path; `donor`/`recipient` cols already filled |
| translocation `X` | ✓ arc | 1 | source arc; cross-chromosome dest ⏳ |

**Phase-1 gap to close next:** the paste/insert **`dest`** for `P`/`D`/`X` and the position for
`O`/`T` live on separate log paths (or are drawn inside `apply` and not returned). Threading them
adds a `dest` column and a couple more `region=` call sites. Inversions and losses — the headline
breakpoint case — are already complete.

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
zombi2 tools export GENOMES_DIR --format {gff,bed,posortho,ffgc,dupinfo,breakpoints} -o OUT
```

Engine = **replay** the enriched trace onto the root gene layout (from the run's `--genes`/`--gff`
input, recorded in the genomes output), reconstructing each node's ordered, signed gene list and
the block (= positional-ortholog) identities — exactly what `NucleotideResult.leaf_mosaic()` /
the `Block` decomposition give in memory, but rebuilt from disk. Then format:

- **`breakpoints`** — read directly from the trace (the broken adjacencies *are* the `I`/`P`/`X`
  rows); no replay needed.
- **`dupinfo`** — count `D` rows per family; no replay needed.
- **`gff` / `bed`** — per-node gene coordinates; from the existing BED writer or the replay.
- **`ffgc`** — extant leaf gene orders (replay) + sequences (`ancestral`/genome FASTA).
- **`posortho`** — block identity across leaves (replay, or the `reconciliations` output).

A single replay code path can serve all of them; the "read directly" cases are just shortcuts.

## 7. Interop with the fork's `zombiExporter` (native-only)

The fork's `_geneorderevents.tsv` uses `TIME, EVENT, BREAKPOINTS, LENGTH`, event codes
`INV/POS/TDUP/DUP/LOSS/ORIG/LFER(_F/_B)/AFER`, and **0-based inclusive** `"start,end"` breakpoint
strings. zombi2's encoding is deliberately native (half-open `[s, s+ℓ)`, single-letter
`EventType` codes); a **`--fork-compat` layer was considered and dropped** — `zombi2 tools export`
is the intended consumer, and a mapping to Krister's exact encoding can be added later if direct
`zombiExporter` interop is ever wanted.

## 8. Phasing

- **Phase 1** ✅ *done* — `Region → EventRecord` (both models) + `--write geneorder` →
  `Geneorder_events.tsv`. Inversions/losses fully described; the missing breakpoint data is now on
  disk. Tests in `tests/test_geneorder_events_output.py`.
- **Phase 1.5** — log the paste/insert `dest` for `P`/`D`/`X` and the position for `O`/`T` (a
  `dest` column + a few more `region=` sites), so every event is fully described.
- **Phase 2** — the replay engine (generalize `tests/test_geneorder_examples.py`'s replay) and
  `zombi2 tools export` for `breakpoints` / `dupinfo` (trace-direct) then
  `gff` / `bed` / `ffgc` / `posortho` (replay).

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
