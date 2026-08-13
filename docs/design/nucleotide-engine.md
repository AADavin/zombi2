# The nucleotide engine, mechanically — groundwork for the indel model

A walkthrough of `zombi2/genomes/nucleotide.py` and its bridge into
`zombi2/sequences/`, written before building indels
([`indels.md`](indels.md)) so the design is checked against what the code does rather
than against what its docstrings say. Every measured number here comes from a run; the scripts are
named where they matter and live outside the repo, under this session's scratchpad at
`/private/tmp/claude-503/-Users-aadria-Desktop-CLAUDE-ZOMBI2/eb40d4cf-e1ba-45fd-8e78-60e5fc83db99/scratchpad/`.
They are probes, not a study — if any of them earns a permanent home it is `analyses/`.

It is subordinate to [`SPEC.md`](SPEC.md).

---

## 1. The shape of a run

`simulate_genomes_nucleotide()` (line 2004) does four things before the loop starts.

**Resolves eleven rates.** A `_scoped` tuple pairs each rate name with its required scope —
`PerLineage` for the seven gene events, `PerChromosome` for fission, fusion and chromosome loss,
`PerLineage` for chromosome origination. `as_rate` fills the level's default where none was written,
and a rate carrying a *different* scope is refused rather than reinterpreted. Modifiers are
restricted to `IMPLEMENTED_MODIFIERS = (OnTime, Driven)`; anything else raises.

**Resolves seven extents.** `_as_bp_extent` requires a `Geometric` base. The stated reason matters
for indels: this engine draws each arc's far end directly out of the legal breakpoint set, so an
arbitrary shape would have to be re-weighted over that set instead of sampled.

**Lays down the initial karyotype.** Either `_replicon_specs` + `_even_gene_intervals`, or `read_gff`
(+ `read_fasta`). `_initial_blocks` makes the alternating gene/intergene chain, one source per
replicon, every block sharing the replicon's initial copy lineage.

**Resolves drivers.** Only a driver on a *rate* joins `trajs` and adds Gillespie breakpoints; a driver
on an extent is read when the event fires and deliberately stays out.

Everything then lands in `_Rates`, a frozen dataclass of eleven rates, seven extents and
`inversion_probability`.

---

## 2. The loop

One global-timeline Gillespie over the complete tree (lines 2310–2469), with `alive` / `gen` / `pos`
as the live-lineage set. Each turn:

1. build `ctx` (`copies`, `lineages`, `chromosomes`, `time`);
2. compute eleven rate totals through `_r`, which sums per lineage when the rate is driven and stays
   pooled when it is not;
3. compute `horizon` — the next species event, every rate's `next_change`, and every driver's
   next switch;
4. draw `t_ev`; if it beats the horizon, pick the class off a cumulative ladder (`b_trl`, `b_trp`, …),
   pick the lineage with `_pick`, read the extent with `_ext`, call the `_do_*`;
5. otherwise advance to the horizon, and if that is a species event, freeze the genome, retire the
   lineage and `_speciate` into the daughters.

`total_length` and `total_chromosomes` are maintained incrementally from what each `_do_*` returns.

**What this costs an indel model.** A new event class is not one edit. It touches `_Rates`, the
`_scoped` tuple, the `_as_bp_extent` calls, the `_extents` dict, the `r_*` block, `total`, `horizon`,
the `b_*` ladder, the dispatch chain and `_ext` — about ten sites inside one function, where a missed
one is silent rather than loud. Worth a single table of event classes instead, if the shape of the
loop is being opened anyway.

---

## 3. The cut set, and the mutators

`Chromosome._legal_cuts()` is the one place the rule lives: a gene contributes its leading edge and
nothing else, an intergene contributes its whole interior, a linear replicon adds its far end.
`_pick_legal_cut` picks uniformly from it; `_pick_arc_extent` picks the far end with weight
`exp(-d/mean)` restricted to that same set, by inverse CDF inside the chosen window.

`_arc_range` checks both ends *before* splitting either, so an illegal arc leaves the chromosome
untouched. `_split_at` → `_split_block` does the splitting and raises `_CutsGene` on a gene. The
mutators — `invert`, `duplicate`, `delete`, `originate`, `excise`, `place`, `transpose` — all take
explicit coordinates, so the engine and a scripted replay run the same code.

Measured (`exp3_semantics.py`):

- An all-gene chromosome has legal cuts `[(0, 0), (100, 100)]` — the two ends only. **Today an indel
  could only ever land in spacer**, which is the wrong biology and confirms the design note's point
  that indels need an exemption from this set.
- `delete(40, 5)` on one 100 bp spacer block leaves **two** blocks, `[0,40)` and `[45,100)`, both
  keeping copy 7, and returns the `Loss` row `(7, 0, 40, 45)`.

---

## 4. Three logs, two coordinate frames

| Log | Records | Frame |
|---|---|---|
| `events` | Origination · Loss · Duplication · Transfer · Speciation | **ancestral** (`source`, `start`, `end`) |
| `rearrangements` | Inversion · Translocation · Transposition | **physical** (chromosome position) |
| `chromosome_events` | the chromosome network | chromosome ids |

**This is the constraint the design note did not know about.** An inversion's breakpoints are
recorded as a physical position on a chromosome, not as an interval on a source. So the set of
ancestral positions that *genuine events* cut at cannot be reconstructed from the logs after the
fact. Any provenance scheme has to record the cut **at the moment it is made** — see §7.

---

## 5. The recovery, in two moves

**Move 1 — `_root_block_partition()`.** For every node's *final* genome and the initial genome, union
`(b.start, b.end)` per source; keep the intervals some node still covers. Bounds are plain integers.
Nothing on a `Block` says which event made a boundary — measured: an inversion of `[40,60)` and a
deletion of `[40,60)` leave the identical bound set `[0, 40, 60, 100]`.

**Move 2 — `_emit_block_events()`, per block.** Scans the whole origination / duplication / transfer /
loss log, BFS-walks the block-copy forest, and emits `GeneEdge` rows. `loss_of` uses an exact
`overlaps` test, so a loss ends a copy only on the blocks it actually covers.

Two passes exist: `_recover()` (declared genes, or every block when none are declared) and
`_recover_blocks()` (`every_block=True`). With no genes declared they do the same work twice.

`assembly()` / `_pieces()` then cut each node's blocks at the partition, with a guard that fires if
the partition fails to cover a block.

---

## 6. What the measurements say

`exp1_partition.py`, 19-node tree, 2 kb genome, `loss_extent=5`:

| `loss` | deletions | root blocks | mean block |
|---|---|---|---|
| 0 | 0 | 1 | 2000 bp |
| 2 | 17 | 35 | 57 bp |
| 5 | 37 | 75 | 27 bp |
| 20 | 156 | 288 | 6.9 bp |

Root blocks come out at **about two per deletion** — exactly `2N + 1` while breakpoints are sparse,
drifting below it as they start to coincide. The partition shatters exactly as predicted.

`exp2b_scaling.py`, 20 kb genome, same tree — the recovery cost against the deletion count:

| deletions | blocks | partition | replay |
|---|---|---|---|
| 406 | 795 | 0.008 s | 0.089 s |
| 1666 | 3050 | 0.057 s | 0.514 s |
| 6386 | 8546 | 0.572 s | 4.921 s |

Fitted exponent over successive pairs: **~1.7 for both moves** — superlinear, heading for quadratic,
because `_emit_block_events` rescans the whole log per block and the block count tracks the deletion
count. 6386 deletions on a 20 kb genome and a 19-node tree already costs 5.5 s. Scale either the
genome or the tree to something realistic and this is the wall, well before memory becomes the issue.

The genealogy grows as **blocks × tree nodes**, not as deletions: 1037 blocks on a 19-node tree gave
18 656 `GeneEdge` rows, about 18 per block.

`exp4_sequences.py` — the same 4000 bp of DNA, evolved whole and evolved shattered:

| | blocks | `_cdf_for` calls | time |
|---|---|---|---|
| no deletions | 1 | 19 | 0.005 s |
| 156 deletions | 302 | 5156 | 0.099 s |

**Twenty times the work for identical total sequence**, purely from fragmentation.

`exp5_cache.py` — a transition-matrix miss (`p_matrix` + `cumsum`) costs 8.4 µs against 0.1 µs for a
cache hit, so **a miss is ~90× a hit**. That is the number behind "splitting a branch at every indel
destroys the cache": each sub-interval is a fresh continuous length that never recurs, so every draw
becomes a miss.

---

## 7. Where [`indels.md`](indels.md) is wrong, and what it missed

**Wrong — the mechanism of the genealogy flood.** I said small deletions would flood the genealogy
with deaths. They do not. `_emit_block_events` uses an exact overlap test, so a partial deletion ends
the copy only on the block it covers and the copy survives on the flanks, which is correct. The flood
is real but its cause is the *partition*: every fragment becomes a root block and every root block
gets a full tree. The conclusion stands; the reason in the note needs replacing.

**Wrong — "give the bounds provenance" as a small change.** §4 shows the two logs are in different
coordinate frames, so this cannot be done as a post-hoc filter over the event log. The provenance has
to be recorded where the cut is made, which means threading a flag through `_arc_range` → `_split_at`
→ `_split_block`, having it survive inversion (where a block's left and right bounds swap meaning),
and travel with blocks through duplication, transfer and translocation. That is the crux of the whole
design and it is the part with no evidence behind it yet.

**Missed — the cost of a new event class in the loop.** Ten edit sites in one function (§2).

**Missed — the double recovery.** `_recover()` and `_recover_blocks()` are two full passes over a
structure whose cost is already ~E^1.7.

**Understated — the exemption from the legal cut set.** Measured in §3: with genes declared, indels
would be confined to spacer entirely. The exemption is not a refinement, it is required for the model
to be worth having.

---

## 8. What I would want settled before writing engine code

1. **How provenance is recorded.** A per-bound flag on `Block`, or a per-source set of event-made cuts
   recorded positively as they are made. The second is additive and cheaper to reason about; neither
   is free, and this decides whether the design works at all.
2. **Whether the event ladder in the loop becomes a table** before an eighth gene event is added to it.
3. **Whether indels get their own log** (`indel_events.tsv`) or ride in `block_events.tsv` with a new
   kind. This decides the read-back path in `read_nucleotide_genomes` and an Appendix B row.
4. **Whether `alignments` becomes gapped** at this resolution — unchanged from the note, still the
   decision with the widest blast radius.

A build order that follows from the above: provenance first, with the partition proved to stay coarse
under a high deletion rate *before* any insertion exists; then the deletion rate as its own event
class; then insertion; then the sequence level's presence mask; then the outputs.
