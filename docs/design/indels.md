# Indels at the nucleotide resolution — a design note

**Status: built, on `feat/indels`.** Breakpoint provenance, the `deletion` and `insertion` rates,
an indel's exemption from the legal cut set, and a root partition that stays coarse under all of it.
A sequence level that reassembles every node's genome and emits a true gapped alignment. A run that
round-trips through the files it wrote, and the command-line flags to drive one. The model is
complete; what remains is a manual chapter and a gallery example.

This note records where an indel model belongs, what stands in its way, and the design move that
clears it. It is subordinate to [`SPEC.md`](SPEC.md): where the two disagree, SPEC
wins.

Claims here were checked against runs before any of it was built — see
[`nucleotide-engine.md`](nucleotide-engine.md), which lists the two places this note was wrong and
the constraint it missed. Building it corrected two more, both marked below.

---

## The question

ZOMBI2 had no indel model. A sequence run evolves a fixed number of sites down each gene tree, so
every node holds the same count and the alignments are gap-free by construction. Adding insertions
and deletions of DNA means deciding two things: which level owns them, and what happens to the
alignment.

---

## Where they belong: Genomes, not Sequences

The obvious home is the sequence level, because that is where a simulator usually puts indels. It is
the wrong one.

`evolve_gene_tree()` samples only the *endpoints* of a branch, from `P(bl) = exp(Q·bl)`, vectorised
over a fixed number of sites. Indels there make a branch a process rather than two draws: the branch
splits at every event, the site indices shift at every node, `_site_classes()` stops working (it
hands out index arrays once per family and shares them across every node, which only holds while
indices stay put), and the transition-matrix cache stops paying — it works today because branch
lengths recur massively across blocks, and sub-interval lengths never recur.

The deeper objection is SPEC §9: a model in which the level's units stop being independent is a new
engine, not a knob. A deletion of several bases removes a run of sites together, so at the sequence
level an indel model breaks site independence.

The nucleotide resolution already does segmental events with an extent, so it does not break
anything there. It already has the deletion: `loss` deletes an arc of DNA. What is missing is the
insertion — novel DNA that is not a gene. `originate()` lays down new material on a fresh source, but
always as a gene, and that restriction is a choice rather than a constraint.

So an indel model is one new event and a second scale of an existing one, at a level built for both.

---

## What that fixes for free

`_laid_out()` walks each node's blocks and sums their sizes. If an indel changes a block, the GFF and
BED offsets, the gene spans, `assembly()` and `trace_back()` all come out right for that node with no
new code.

Had indels lived at the sequence level, the annotation would describe a genome the assembled FASTA no
longer matched — the nucleotide level writes coordinates in root-block sizes, and the sequence level
would have changed them afterwards.

---

## The obstacle: the root partition

`_root_block_partition()` reads the final genome at every node and unions `(b.start, b.end)` over all
blocks. Those bounds carry no provenance: an indel boundary is indistinguishable from an inversion
boundary — measured, the two leave the identical bound set.

Each deletion therefore adds **about two root blocks**, exactly `2N + 1` until breakpoints start to
coincide. On a 2 kb genome the mean root block falls from 2000 bp to 6.9 bp after 156 deletions. The
recovery cost grows as roughly **E^1.7** in the deletion count, because `_emit_block_events()`
rescans the whole event log once per block: 6386 deletions on a 20 kb genome and a 19-node tree
already costs 5.5 s. The genealogy grows as blocks × tree nodes. And the same 4000 bp of DNA costs
**20× more** to evolve at the sequence level shattered than whole.

This is structural, and giving indels their own event class and their own log — the way the
ancestry-neutral events already use `rearrangements` rather than `events` — does not touch it. That
separation is still wanted for deletion, for a different reason: a deletion changes a copy's size, it
does not end one, so it should not write a `GeneEdge`. But it is not the fix.

**The two indels are not symmetric, and this note first assumed they were.** A deletion ends no copy
lineage. An insertion necessarily *begins* one: novel DNA descends from nothing, so it can only be its
own ancestor, which means a fresh source and a fresh copy lineage — or the recovery finds no root and
builds no block tree for it. So an insertion belongs in `events` as a third kind of root beside
`initial` and `origination`, and is `origination` minus the gene: it brings **sequence**, where
origination brings a **gene family**. Built that way.

---

## The fix: an indel changes presence, not partition

Give the bounds provenance and have the partition ignore indel-made cuts.

An indel then becomes a presence fact **inside** a root block rather than a boundary between two of
them, and the partition stays exactly as coarse as it is today at any indel rate.

Measured, once built: 631 deletions leave the partition at **one** root block of 2000 bp, against 860
blocks averaging 2.3 bp for the same events as `loss`. With genes declared it stays exactly the
initial layout, and every gene still recovers as one root block.

**Insertion cannot be free in the same way, and the note should have said so.** Its breakpoint is
skipped like a deletion's, so the block it lands inside stays whole — measured, the original 2 kb
holds at one root block through 1252 insertions. But the material it brings is a new source, so it is
a root block of its own. Growth is therefore **one block per insertion**: linear, and it never
fragments what is already there, which is the difference that matters. The cost is real and is the
argument for the union-column design at the sequence level, where those runs become columns of the
block they landed in rather than blocks in their own right.

This is the crux, and it is not the small change this note first implied. It cannot be done as a
filter over the event log after the run, because the two logs are in **different coordinate frames**:
`events` records ancestral intervals, `rearrangements` records physical chromosome positions, so the
ancestral positions that genuine events cut at are not recoverable after the fact. The provenance has
to be recorded where the cut is made — threaded through `_arc_range()` → `_split_at()` →
`_split_block()`, surviving inversion (where a block's two bounds swap meaning) and travelling with
blocks through duplication, transfer and translocation. It holds: only bounds no genuine event ever
cut at are skippable, so a position an inversion also reached in another lineage stays in the
partition. Where that recording *lives* changed later — see below, it belongs on the events
themselves — but what is recorded did not.

One consequence falls out of the same change. Today no breakpoint may fall strictly inside a declared
gene, which would confine indels to spacer — wrong, since indels in coding DNA are most of the
biology. An indel may cut inside a gene precisely because its cut never reaches the partition, so the
gene still recovers as one root block with one tree. This has to be written as a property of the cut
set in `Chromosome._legal_cuts()`, not as a check at each mutator: a rule enforced per call site is a
rule some call site forgets.

**Built.** `_legal_cuts(indel=True)` opens every position, gene or not, and `_check_cut()` and
`_split_at()` take the same exemption; the two guards that refuse an ordinary event a gene's interior
are untouched. Measured on a genome that is 100% gene — ten 200 bp genes in 2000 bp, no spacer
anywhere, where before this an indel had nowhere at all to go: 499 deletions and 467 insertions land,
190 gene copies end up sitting in more than one block, and all ten genes still recover as exactly one
root block with one tree.

One thing it forced. `Chromosome.n_genes` counted *blocks*, and a gene cut by an indel is two of
them — so the guard that stops a chromosome losing its last gene would have counted two genes where
there is one and let the last one go. It counts distinct `(family, copy)` now.

---

## What the sequence level then evolves: the alignment

Each root block hands the sequence level a **union column set** — every position that ever existed
anywhere in the tree for that block — plus a presence mask per node. The sequence level evolves the
columns; each node's genome is the columns it carries.

The engine barely changes, because columns never move:

- `_site_classes()` keeps working as index arrays drawn once per family.
- The transition-matrix cache is untouched: branch lengths are still per branch, with no splitting.
- `_sample()` stays vectorised.

The one addition is that columns inserted on a branch are founded from the model's stationary
frequencies at that branch instead of inherited from the parent.

A run with no indels has a union set equal to the block, so it takes no extra draw and stays
bit-identical to what it is today.

---

## The cost, and how much of it is real

Evolving the union means evolving columns a node does not carry. A column deleted on some branch
still gets a state drawn at every node below the deletion. A column inserted on some branch is in the
array from the start, so it is evolved everywhere above and beside its insertion, and all of that is
discarded when the branch founds it afresh.

The union is about `L × (1 + per-lineage insertion × tips)`, since insertions in different lineages
land in different places and union rather than overlap. A hundred tips inserting 5% each gives about
a sixfold column set. That is not an artifact of the design: a true alignment of many sequences
really is mostly gaps, and any simulator that emits one pays for that width.

The compute half is avoidable. Carry the presence mask down the walk alongside `parent_states` and
sample only the columns the node holds — the same idiom the across-site-class path already uses, one
`_sample()` per group over its own index array. Presence changes only on branches where an indel
fired, so most nodes reuse the parent's index array unchanged and per-node work is what it is today.
Memory goes the same way: store each node at its own width, and keep presence as a few `(start,
end)` pairs, since indels are arcs.

---

## The one place this can be wrong

Two insertions that happen at the same offset in two different lineages are **not** the same column.
Keyed by coordinate they would merge into one, and the alignment would then assert homology between
material that was never related — the one error a simulator exists not to make.

Key columns by the insertion event that created them, which is the lesson `Block.copy` already
teaches one level up: identity, not coordinate. Fix a deterministic order for columns inserted at the
same offset, so a seeded run reproduces its alignment exactly.

---

## Rates and extents

Two named rates, `insertion` and `deletion`, per lineage like the other segmental rates, each with
its own extent in base pairs. Two are needed rather than one: no single geometric extent gives both
3 kb losses and 5 bp indels.

Because an indel is exempt from the legal cut set, its extent is not drawn from that set, so it is
free to take any distribution — unlike the existing segmental extents, which the nucleotide engine
restricts to geometric because it draws each far end directly from the legal breakpoints. That
matters: indel sizes in real data are closer to a power law than to a geometric.

**Built, and the restriction was load-bearing rather than inherited.** An extent reaches the mutator
as a *number*: `_ext()` calls `Extent.mean()`, and the drivers sampled a geometric around it. So
accepting a shape the engine could not express would have taken `Fixed(1)` and silently sampled it
geometrically — one nucleotide asked for, a tail out to five delivered. The fix is not to relax the
check alone but to stop collapsing the extent: on the indel path `_ext_draw()` calls
`Extent.sample()`, which already existed and already works with any base, scaling the draw rather
than the distribution's parameter. A deletion then stops going through `_pick_arc_extent()`
altogether, since that function exists to pick a far end out of a cut set an indel no longer has.

`Fixed(1)` now means exactly one nucleotide, measured; `Fixed(n)` removes exactly `n` per event; a
`scipy.stats.zipf` extent gives the power law. Deletion and insertion also stopped disagreeing at
small means, which they did while one drew through the cut set and the other did not.

One thing to keep in view: a heavy-tailed extent has no upper bound, so a single "deletion" can take
a large stretch — at which point it is doing `loss`'s job, and the vocabulary line between them (what
a lineage *has* against how much of a surviving copy it *carries*) is only as sharp as the extent
keeps it.

A deletion that removes a copy entirely is still recorded as a loss, because it is one. The
genealogy consequence stays emergent from the extent, which is the rule the level already uses for a
gene being engulfed whole rather than split.

---

## What the sequence level actually needed

Far less than this note assumed, and one thing it did not foresee.

The union column set is **already** the root block. The sequence level evolves each root block whole
— `per_block[i] = (b - a, …)` is its full ancestral extent — so the columns a node does not carry
are drawn and never read, which is exactly the waste described below, already paid. The presence
mask is not a new structure either: it is the sub-range each piece of `assembly()` carries. So
`evolve_gene_tree()`, `_site_classes()` and the transition-matrix cache are untouched, and the whole
change is that a piece became `(block, gene, strand, lo, hi)` and the assembly slices `[lo:hi]`
before flipping.

**What it did not foresee: an indel breakpoint is invisible to the partition but is still a real
boundary in the genome, and other events act on it.** A duplication can copy *part* of a root block;
a loss's arc is contiguous physically but spans several disjoint stretches of one, and only its two
ends become partition bounds. That broke the recovery three ways, each needing a different answer:

- **Duplication and transfer** now beget a block-copy when they *overlap* the block rather than
  cover it. The child carries the part it was given; which part is `assembly()`'s business, not the
  tree's.
- **A copy's death** cannot be decided from any one event's reach. `covers` says no death where
  there was one; `overlaps` says a death where the copy still holds the rest. The exact question is
  whether the copy carries any of the block at the end — `_carried_blocks()` — and the *time* is the
  last removal that touched it, not the first.
- **Attribution is of use, not creation.** An ordinary event that merely starts its arc where an
  indel once cut still needs that position in the partition; recording only splits that actually
  happen leaves it indel-attributed and skippable. Two corollaries: the initial layout's boundaries
  are structural and recorded as event-made at setup (or an indel landing on a gene edge lets the
  partition merge that gene away), and each source's origination bracket is unioned back past the
  skip.

Each of the three reads identically on a run without indels, which is why the suite kept passing
until one was wrong.

**The alignment is gapped**, and that answers the question this note left open. A block is evolved
over its whole ancestral extent, so an ungapped row hands back bases the lineage deleted — a 120 bp
row for a copy carrying one base of the block. Each row is now the block's full width with ``-``
where the lineage carries nothing, which is the true alignment and costs nothing to say: the
sub-ranges already record where the gaps go. It is written in place even though the assembly reads
the same strings, because a gap only ever falls *outside* a carried sub-range, so slicing one out
still yields unbroken sequence — and the assembled genome has no gap in it anywhere.

One thing that had to move with it: `_identity_counts()` walks the bytes actually present rather
than the model's alphabet, so two lineages that deleted the same stretch would have read as identical
across it. A gap is not a residue and two of them are not a match; pairs are counted per column over
the sequences that have a residue there, which on an ungapped alignment is every pair at every
column, exactly as before.

---

## The provenance had to move into the records

Found by building it. Kept as a log beside the run, the provenance **could not survive a write**: the
ancestral position an insertion opened its gap at was written nowhere, and deletions had no file at
all. Worse, it was not derivable after the fact either — `rearrangements` are in physical coordinates,
so there was no way to ask which ancestral positions an inversion cut at.

The fix is not a file for the side channel but to stop having one. **Every event record now carries
the breakpoints it used**, as ancestral `(source, position)` pairs: `Loss.cuts`, `Inversion.cuts` and
the rest, collected by `_arc_range()` as the splits happen and written as one column in
`block_events.tsv` and `rearrangement_events.tsv`. The indel log rides in the same table under kind
`deletion`. `_skippable_bounds()` then *reads* the logs rather than consulting bookkeeping kept
alongside them, and the distinction survives a write because it is in the record.

Three things fall out. The two coordinate frames the walkthrough complained about are now one, since
a rearrangement finally says where it cut in ancestral terms. The initial layout's structural
boundaries come from `initial_genome` rather than being recorded at setup. And a nucleotide run with
indels round-trips: partition, gene trees and assembly all identical after a write and a read, over
1762 indels.

The one thing this costs is a column on the **ordered** resolution's rearrangement file, which will
never fill it — that file is shared by both resolutions, and the ordered one has no ancestral
coordinates to give.

---

## What still needs deciding

- **Does `alignments` become gapped at this resolution?** The design forces yes — the presence mask
  is the gap pattern. That reaches Appendix B, the `zombi2 tools` readers and Phylustrator's
  alignment reader, so it is a decision, not a detail.
- **Frameshift.** The clean core has no codon model, so there is no frame to break today. If the
  codon models return from `legacy/`, an indel must be constrained to a multiple of three or
  explicitly allowed to frameshift.
- **`n_genes` must count distinct gene ids, not blocks.** Once an interior deletion can leave a gene
  as two physical blocks, the count behind "a chromosome never exists without a gene" miscounts
  silently.
- **Scope.** The family and ordered resolutions have no coordinates, so they get no indels under this
  design. Whether gapped alignments are wanted there is a separate question.

---

## A fossil to purge

The module docstring of `zombi2/genomes/ordered.py` already tells the reader that the nucleotide
resolution has indels. It does not. That line is wrong today and should be fixed whether or not this
note is ever built.
