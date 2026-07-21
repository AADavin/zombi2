# Porting the Swenson test suite (`thekswenson/Zombi`) to ZOMBI2

What was ported from Krister Swenson's fork test suite (`thekswenson/Zombi`, `root_genome` branch,
`tests/`), what ZOMBI2 already covered, and what was deliberately left — with reasons, so nothing
valuable is silently dropped. Read it alongside his tests.

> **History.** An earlier port was done against the *pre-reset* codebase. The clean-core rewrite
> quarantined it in `legacy/tests/`, where it targets names that no longer exist
> (`nucleotide_genome.SegmentRegistry`, `IdManager`, `TargetParams`). Those files are reference
> material, not code to revive. This note describes the port against the **clean core**.

## Why it is not a transcription

- **The fork** is the ZOMBI1 lineage: a `GenomeSimulator` monolith, the `zombi T/G/S` command line,
  `Filenames.py`, and a gene-order model of genes + intergenes with its own *specific* and *total*
  coordinate spaces.
- **ZOMBI2** is a ground-up rewrite: `species / genomes / sequences / traits`, the
  `zombi2 <level>` CLI, and three genome **resolutions** — unordered ⊂ ordered ⊂ nucleotide.

So a port re-expresses each scenario against a different model. Every concept survives; the numbers
have to be re-derived.

## Concept map

| Fork | ZOMBI2 |
| --- | --- |
| *division* | `Block` — a run of one unbroken ancestry (`genomes/nucleotide.py`) |
| *pieces* | the ordered block list; `blocks.tsv` |
| `ch.genes` (signed gene order) | the gene-carrying blocks in physical order; `gene_order.tsv` at the ordered resolution |
| `natural_cuts` / `init_divisions` | block boundaries |
| `make_inversion_intergenic(ch, bp1, bp2, dir)` | `Chromosome.invert(start, length)` — a half-open arc, not two breakpoints and a direction |
| `cut_and_paste` / `obtain_segment` | `_transpose` (ordered) / the transposition mutator (nucleotide) |
| `Inversion.afterToBeforeS` / `afterToBeforeT` | `Chromosome.trace_back()` — per-nucleotide ancestral origin, root-anchored rather than per-event |
| `Geneorder_events_per_branch/` | `genome_event_positions.tsv` + `rearrangements.tsv` |
| `All_genomes/` | `gene_order.tsv` / `blocks.tsv`, which now carry **every node**, ancestors included |
| `Leaves.tsv` | `species_extant.nwk` |

## Disposition, file by file

| Fork test | What it exercised | Where it landed |
| --- | --- | --- |
| `test_geneorder_events.py` (4) | replay the per-branch gene-order events → reconstruct the genomes | **Ported** → `tests/test_genomes_event_positions.py`. Needed new output first — see below. |
| `test_events.py` (26) | exact gene order / orientation / intergene lengths after a scripted inversion | **Ported (inversions)** → `tests/test_genomes_geneorder_examples.py`. Its D/T/L/O worked examples are **not** portable — see below. |
| `test_divisions.py` (8) | `natural_cuts` / `init_divisions` boundaries after inversions | **Ported** → `test_blocks_tile_the_chromosome_and_split_only_at_the_cuts`, plus the composing-inversions case. |
| `test_randomization.py` (5) | same seed → identical output directory | **Ported** → `tests/test_cli_determinism.py`, widened to the whole pipeline. |
| `test_genomes.py` (2) | `cut_and_paste`, coordinate exclusion | **Already covered** by `test_genomes_ordered.py` (`test_transpose_relocates_a_segment_within_the_chromosome_preserving_ids` and its flipping twin). Not duplicated. |
| `test_commandline.py` (5) | modes run; `All_genomes` vs `Genomes` crosscheck | **Already covered** by `test_cli.py` (files exist, per level) plus the replay test (the crosscheck's real content). |
| `test_pieces.py` (4) | divisions/pieces after events | **Not ported** — ~80% commented-out stubs in the fork; the live assertions are a subset of `test_divisions.py`. |

## What the port needed built first

Krister's best test could not run against the clean core as it stood, so three things changed. Each
is a feature in its own right, not test scaffolding:

1. **`gene_order.tsv` covers every node** (was extant tips only). A branch's rearrangements are
   meaningless without the genome it started from, which is its parent's row set.
2. **`genome_event_positions.tsv`** — where each D/T/L/O event happened. The gene-genealogy log is
   position-blind on purpose (identity and descent are resolution-blind), so the coordinates go in a
   companion table, the same split ZOMBI2 already makes for `rearrangements.tsv`. A transfer writes
   one row per branch, following the fork's leaving/arriving split, except that both rows name the
   edge outright instead of being matched by timestamp.
3. **A `write()` for the nucleotide resolution**, and `zombi2 genomes --resolution nucleotide`.
   Before this the whole resolution was in-memory only.

## What was deliberately not ported

- **The fork's per-nucleotide `afterToBeforeS`/`afterToBeforeT` arithmetic** (~6 cases). Its content
  is already proved for arbitrary inputs by the oracle in `test_genomes_nucleotide.py`, which checks
  random inversions against an independent brute-force array. One scenario is kept
  (`test_trace_back_maps_every_position_home`) to fix the *convention* — that is what a worked
  example is for; the rest would add bulk, not coverage.
- **The fork's tandem-duplication / loss / origination / transfer worked examples.** At the
  nucleotide resolution only inversion takes explicit coordinates; the others draw their position
  from the rng, so scripting one means reaching past the engine's public surface. The *properties*
  they assert are covered by the ordered-level tests and the oracle. **If it is worth making these
  mutators scriptable, this is the reason to.**
- **`test_pieces.py`**, which is mostly commented out upstream.

## The gap that is still open

The file-based replay (`tests/test_genomes_event_positions.py`) holds to **one chromosome per
genome**. The chromosome tier re-mints ids at speciation, fission and fusion, so replaying a
multi-chromosome run additionally needs `chromosome_events.tsv` to map a parent's chromosome onto
its daughters'. That is a separate question from whether the coordinates are right, which is what
the test is for.

## Running the ported tests

```
pytest tests/test_genomes_geneorder_examples.py \
       tests/test_genomes_event_positions.py \
       tests/test_cli_determinism.py
```
