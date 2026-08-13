# Indels at the nucleotide resolution — a design note

**Status: a proposal. Nothing in the code does this yet.** This note records where an indel model
belongs, what stands in its way, and the one design move that clears it. It is subordinate to
[`SPEC.md`](SPEC.md): where the two disagree, SPEC wins.

The claims here have since been checked against runs — see
[`nucleotide-engine.md`](nucleotide-engine.md), which also lists the two places this note was wrong
and the constraint it missed.

---

## The question

ZOMBI2 has no indel model. A sequence run evolves a fixed number of sites down each gene tree, so
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
separation is still wanted, for a different reason: an indel changes a copy's size, it does not end
one or start one, so it should not write a `GeneEdge`. But it is not the fix.

---

## The fix: an indel changes presence, not partition

Give the bounds provenance and have the partition ignore indel-made cuts.

An indel then becomes a presence fact **inside** a root block rather than a boundary between two of
them, and the partition stays exactly as coarse as it is today at any indel rate.

This is the crux, and it is not the small change this note first implied. It cannot be done as a
filter over the event log after the run, because the two logs are in **different coordinate frames**:
`events` records ancestral intervals, `rearrangements` records physical chromosome positions, so the
ancestral positions that genuine events cut at are not recoverable after the fact. The provenance has
to be recorded where the cut is made — threaded through `_arc_range()` → `_split_at()` →
`_split_block()`, surviving inversion (where a block's two bounds swap meaning) and travelling with
blocks through duplication, transfer and translocation. Nothing has been built or measured here yet;
it is the part of the design that could still fail.

One consequence falls out of the same change. Today no breakpoint may fall strictly inside a declared
gene, which would confine indels to spacer — wrong, since indels in coding DNA are most of the
biology. An indel may cut inside a gene precisely because its cut never reaches the partition, so the
gene still recovers as one root block with one tree. This has to be written as a property of the cut
set in `Chromosome._legal_cuts()`, not as a check at each mutator: a rule enforced per call site is a
rule some call site forgets.

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

A deletion that removes a copy entirely is still recorded as a loss, because it is one. The
genealogy consequence stays emergent from the extent, which is the rule the level already uses for a
gene being engulfed whole rather than split.

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
