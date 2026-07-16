# Porting the Swenson test suite (`thekswenson/Zombi`) to zombi2

This note records what was ported from Krister Swenson's fork test suite
(`thekswenson/Zombi`, `root_genome` branch, `tests/`), what was **already covered** by zombi2's
own suite, and what was **deliberately not ported** (and why). It is meant to be read alongside
Krister's original tests so nothing valuable is silently dropped.

## Why it isn't a mechanical copy

The two codebases have diverged architecturally:

- **The fork** is the ZOMBI‑1 lineage: a `GenomeSimulator` monolith, the `zombi T/G/S`
  command line, `Filenames.py`, and a gene‑order model built directly on genes + intergenes
  with its own *specific*/*total* coordinate spaces.
- **zombi2** is a ground‑up rewrite: `species / genomes / traits / sequences` subpackages over a
  Rust core, the `zombi2 species / genomes / sequence` CLI, and a **nucleotide/segment/block**
  genome model (`NucleotideGenome`, `Segment`, `SegmentRegistry`, `Block`).

So a port is a **re‑expression against a different internal API**, not a transcription. The good
news: every concept the fork tests rely on has a counterpart in zombi2 — it just lives in the
nucleotide layer, under different names.

## Concept map (fork → zombi2)

| Fork concept | zombi2 counterpart |
| --- | --- |
| *division* (maximal uncut interval) | `Block` — *"a maximal uncut interval `[start,end)` of one ancestral source"* (`nucleotide_sim.py`) |
| *pieces* (genome as ordered divisions) | ordered `Segment` list / `NucleotideResult.leaf_mosaic()` |
| coordinate map `Inversion.afterToBeforeS/T` | `SegmentRegistry.provenance` (`seg_id → (source, src_start, src_end)`) + `to_cells()` / `trace_back()` (root‑anchored, not per‑event) |
| `make_inversion_intergenic(ch, bp1, bp2, dir)` | `NucleotideGenome._apply_inversion(start, length)` |
| `cut_and_paste` / `obtain_segment` | `NucleotideGenome._apply_transposition(start, length, dest)` |
| `init_divisions` / `natural_cuts` | segment / block boundaries (the ancestral source‑intervals) |
| `Geneorder_events_per_branch/` files | `geneorder_events.tsv` (`--write geneorder`) — one file, one row per structural event, keyed by `branch` (plus the in‑memory `EventLog`) |

## Disposition of each fork test file

| Fork test | What it exercised | Disposition in zombi2 |
| --- | --- | --- |
| `test_events.py` | exact gene order / orientation / intergene lengths / coordinate maps after scripted inversions | **Ported (scenarios)** → `test_geneorder_examples.py`. Distinct scenarios kept as gene‑aware golden cases; the per‑nucleotide `afterToBeforeS/T` arithmetic is **already covered** (see below). |
| `test_genomes.py` | `cut_and_paste`, coordinate‑exclusion (transposition mechanics) | **Ported (scenarios)** → transposition cases in `test_geneorder_examples.py`. |
| `test_divisions.py` | `natural_cuts` / `init_divisions` boundaries after inversions | **Ported** → `test_inversion_subdivides_only_the_cut_intergene`. |
| `test_pieces.py` | divisions/pieces after events | **Not ported** — the live assertions are a subset of the above; the file was ~80% commented‑out stubs in the fork. |
| `test_commandline.py` | modes run + `All_genomes` vs `Genomes` crosscheck | **Already covered** by `test_cli.py` (file existence) + block/reconciliation invariants; not re‑ported. |
| `test_randomization.py` | same seed → identical output directory (T/G/S) | **Ported** → `test_pipeline_determinism.py` (CLI, `species → genomes → sequence`). |
| `test_geneorder_events.py` | replay per‑branch gene‑order events → reconstruct each genome | **Ported** → `test_geneorder_file_replay.py` (the file‑based replay). The in‑memory analogue also exists (`test_nucleotide_genome.py::test_mosaic_reassembles_each_leaf`), but it is *not* a substitute — see below. |

## What was added

- **`tests/test_geneorder_examples.py`** (9 tests) — gene‑aware inversion and transposition
  worked examples on the fork's `30_10.gff` / `30_6.gff` genomes, re‑seeded as gene‑annotated
  `NucleotideGenome`s. Each asserts the exact gene order, orientation, content conservation and
  division boundaries, **derived by hand** — correctness checks, not regression pins.
- **`tests/test_pipeline_determinism.py`** (1 test) — the whole‑pipeline CLI analogue of
  `test_randomization.py`: two same‑seed projects run through `species → genomes → sequence` and
  the output directories are compared by decompressed content (gzip headers carry an mtime;
  `.log` files record wall‑clock timing and are excluded).

## What was already covered (so it was **not** duplicated)

zombi2's `tests/test_nucleotide_genome.py` already stress‑tests the rearrangement geometry the
way the fork did in bulk: random `(s, ell)` inversions / transpositions / deletions / duplications
across 25 seeds, each checked against an **independent brute‑force array oracle**
(`ArrayGenome`), plus content‑conservation, bijection and involution invariants, block tiling,
mosaic reassembly, per‑block gene trees and reconciliation. Krister's *stress‑testing philosophy
is therefore already present* — the gap those tests leave is the **gene annotation**, which is
exactly what the new golden examples fill.

## What was **not** ported, and why

- **The fine‑grained `afterToBeforeS/T` coordinate‑arithmetic examples** (~6 in `test_events.py`):
  they assert per‑nucleotide coordinate mapping, which the random‑vs‑oracle tests already prove
  for arbitrary inputs. Porting them verbatim adds bulk, not coverage. The *distinct scenarios*
  they encode (wrapping arcs, enclosed genes, gene reordering) are preserved as golden cases.
- **`test_pieces.py`**: mostly commented‑out stubs in the fork.

## Correction: the file‑based replay *is* ported

An earlier revision of this note declined the file‑based `Geneorder_events_per_branch` replay
(`test_geneorder_events.py`) on the grounds that "zombi2 does not emit a per‑branch gene‑order event
file, so there is nothing to replay from disk", and that the in‑memory reconstruction already
covered it. **Both halves of that were wrong, and Krister was right to push back.**

- zombi2 *does* now emit the event log — `--write geneorder` → `geneorder_events.tsv`, one row per
  structural event with its breakpoints, keyed by `branch`. The stated blocker is gone.
- More importantly, the in‑memory test is **not a substitute**. zombi2's output is a set of *files*,
  and the files are what a user inferring rearrangements actually consumes. Verifying the in‑memory
  structures says nothing about whether the *written* breakpoints mean what they claim: a
  coordinate‑convention slip would yield a plausible file that replays to the wrong genome, and every
  in‑memory test would still pass.

`tests/test_geneorder_file_replay.py` closes this: it reads the written `geneorder_events.tsv`,
replays each branch's events onto the root genome with the simulator's own primitives, and asserts
the result reproduces the written genomes (`bed/<node>.bed`). A negative control (perturbing the
replay by a single base) confirms the test fails when the correspondence breaks, so it is not
vacuous. Scope today is a content‑conserving run (inversion + transposition) — the
rearrangement‑inference case; duplication / loss / transfer are the natural extension once their
file‑level semantics are pinned down.

## Notes for the merge

- The new gene‑order tests call the `_apply_inversion` / `_apply_transposition` primitives —
  matching the convention already used throughout `test_nucleotide_genome.py`. The
  multi‑chromosome refactor (`feat/nucleotide-translocation`, merged as PR #129) gave both
  primitives an **optional** `chrom=None` argument; the tests were re‑validated on that merged
  code and pass **unchanged**.

## Running

```
pytest tests/test_geneorder_examples.py tests/test_pipeline_determinism.py \
       tests/test_geneorder_file_replay.py
```
