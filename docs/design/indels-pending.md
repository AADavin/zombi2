# Indels — the user-facing changes this branch owes

**A staging area, not a document to publish.** While `feat/indels` is in flight, everything it changes
that a user could observe is written here instead of in `CHANGELOG.md` and Appendix B. When the work
is ready to publish, these entries are moved across in one edit and this file is deleted.

The reason to stage rather than edit in place: a branch that may be reshaped several times would
otherwise churn the changelog, and half-built entries under `[Unreleased]` read as shipped.

Rules for keeping this honest:

- **Add the entry in the same change as the code**, exactly as `CLAUDE.md` requires — the only
  difference is which file it lands in.
- **Write it in final form**, so publishing is a move and not a rewrite.
- **Delete an entry if the change is reverted.** A staged entry for work that no longer exists is
  worse than no entry.

---

## For `CHANGELOG.md`, under `[Unreleased]` → `### Added`

```markdown
- Indel deletions at the nucleotide resolution: `deletion` (per lineage) with `deletion_extent`
  removes an arc the way `loss` does, but ends no copy lineage, so it is recorded in `.deletions`
  and stays out of the genealogy. Its breakpoints do not cut the root partition, which is what lets
  a genome carry thousands of them and still recover the gene and block trees it would have had
  without them. `loss` changes what a lineage has, `deletion` how much of a surviving copy it
  carries. Python only for now, and `assembly()` — so `simulate_sequences` — refuses a run that used
  it. (#TBD)
```

```markdown
- Indel insertions at the nucleotide resolution: `insertion` (per lineage) with `insertion_extent`
  lays down a run of novel spacer at a legal position — `origination` without the gene. Novel DNA
  descends from nothing, so it arrives on a fresh source under a fresh copy lineage and is recorded
  in the event log as a root of its own kind, counted apart from `origination` and `initial` in
  `genome_summary.json` and written as `insertion` in `block_events.tsv`. The breakpoint it makes to
  open a gap does not cut the root partition, so the block it lands inside stays whole. Python only,
  and a run that used it can be neither assembled nor read back from files. (#TBD)
```

```markdown
- An indel may now fall **inside a gene**, where before it could only land in spacer: a breakpoint may
  never cut a gene, but an indel's breakpoint never reaches the root partition, so the gene still
  recovers as exactly one block with one tree and only how much of it a lineage carries has changed.
  A genome with no spacer at all evolves under indels rather than standing still. (#TBD)
```

```markdown
- An indel extent may take **any distribution**, where a segmental one must still be geometric:
  `deletion_extent=Fixed(1)` gives single-nucleotide indels, and a `scipy.stats.zipf` extent gives
  the power-law length distribution indels actually have. The segmental events sample an arc's far
  end out of the genome's legal breakpoints and so are parameterised by a mean; an indel may cut
  anywhere, so it draws its size outright and the shape survives. (#TBD)
```

```markdown
- **BREAKING:** a piece of `result.assembly(node)` is now `(block, gene, strand, lo, hi)` and one of
  `result.initial_assembly()` is `(block, strand, lo, hi)`, where `[lo, hi)` is the sub-range of that
  block the node carries. An indel leaves a lineage holding part of a block rather than all of it, so
  a piece can no longer be a whole block; without indels `lo` is 0 and `hi` the block's length, so
  the shape says what it always did at one extra pair of numbers. `simulate_sequences` reassembles
  every node's genome from an indel run accordingly. (#TBD)
```

```markdown
- A nucleotide run with indels writes a **true gapped alignment**: each row of a block is that
  block's full ancestral width, with `-` where that lineage carries nothing. The assembled genome
  FASTA is unchanged and gap-free — the alignment is what was evolved, the genome is what exists.
  Mean pairwise identity no longer counts a shared gap as a match. A run without indels is
  byte-for-byte what it was. (#TBD)
```

```markdown
- A nucleotide run with indels can be written and read back: every event record now carries the
  ancestral breakpoints it used, in a new `cuts` column of `block_events.tsv` and
  `rearrangement_events.tsv`, and the indel log rides in the first of those under kind `deletion`.
  What tells an indel breakpoint from an ordinary one is now in the record rather than in the run,
  so the root partition, the gene trees and the assembly come back exactly as they were. A
  rearrangement also finally states where it cut in **ancestral** coordinates, where before that log
  was physical-only and could not be joined to the block log at all. (#TBD)
- **BREAKING:** `rearrangement_events.tsv` gains a `cuts` column at both resolutions. The **ordered**
  resolution always leaves it empty — it has no ancestral coordinates to give — but the file is
  shared by both, so the column is there. (#TBD)
```

```markdown
- Indels reach the command line: `--deletion` and `--insertion` with `--deletion-extent` and
  `--insertion-extent`, at `--resolution nucleotide` only. The two extent flags take the written form
  rather than a bare number, so `--deletion-extent 'Fixed(1)'` is a single-nucleotide indel on the
  command line exactly as it is from Python. Either rate under `--resolution family` or `ordered` is
  refused rather than ignored. (#TBD)
```

```markdown
- `result.alignment_with_insertions(block)` reads a block's alignment with the material inserted into
  it spliced back in as real columns, so a lineage that gained a run shows its bases where the others
  show gaps. `alignments` holds what evolved — an insertion is a block of its own there, with its own
  tree — and a gene's alignment therefore showed every deletion and no insertion. Which lineage
  carries which run is read off the genome's layout, so a duplicated gene's two copies each get their
  own, and a run carried away from its host by a rearrangement is left where it now is. (#TBD)
```

```markdown
- A gallery example for indels: the real *M. genitalium* `MG_RS01985`, 161 bp, across five species —
  172 alignment columns, a 6 bp insertion carried by exactly one clade, a second insertion in one
  species and a 3 bp deletion in another. A gap means two opposite things and only the column says
  which. (#TBD)
```

Nothing else is left to stage for the model itself. Still outstanding, and not a changelog entry: a
manual chapter for Chapter 6 and a gallery example.

---

## For Appendix B (`manual/book/appendix-b.md`)

Both under **Genomes, nucleotide: `simulate_genomes_nucleotide`**.

**1. Extend the `genome_summary.json` row.** Append to the end of its Contents cell, after
"…and warns that it will":

> `deletions` and `base_pairs_deleted` count the indel log, which is neither of the other two: an
> indel removes material without ending a copy lineage, so it appears in no event counter

**2. Add a row**, after the `Driver views` row:

| Output | File | Format | Default | Contents |
|---|---|---|---|---|
| Assembly | `.assembly(node)` · `.initial_assembly()` | dict | Python | **changed shape:** a piece is `(block, gene, strand, lo, hi)`, `[lo, hi)` being the sub-range of that block the node carries — a whole block unless an indel took a stretch out of it |
| Indel log | `.deletions` | list | Python | one `Deletion` per indel — `time` · `lineage` · `chromosome` · `deleted`, the last being `(copy, source, start, end)` rows in **ancestral** coordinates, the same payload a loss carries. Kept out of `events` because a deletion ends no copy lineage, so the gene-tree recovery must not read it, and out of the root partition because cutting there would shatter it. Only deletions are here: an **insertion** brings material that descends from nothing, so it must begin a copy lineage and is recorded in `events` as a root of its own kind instead. Written into `block_events.tsv` under kind `deletion`, which is also how it is read back |

**3. Extend the `block_events.tsv` row.** Its Contents cell lists the kinds and now needs the new
column too:

> Kinds are `initial` · `origination` · `insertion` · `duplication` · `deletion` · `loss` ·
> `transfer` · `speciation`, where `insertion` and `deletion` are the indels — `insertion` is novel
> spacer on a fresh source, a root like an origination but carrying no gene family, and `deletion`
> removes material without ending a copy lineage. A last column, `cuts`, gives the **ancestral**
> breakpoints that event used as `source:position` pairs, semicolon separated: it is what tells an
> indel breakpoint from an ordinary one, and without it a run could not be read back

**4. Extend the `rearrangement_events.tsv` row** with the same column, and say that the ordered
resolution leaves it blank.

**5. Say that the indel flags exist** wherever Appendix B lists what `zombi2 genomes` takes:
`--deletion` · `--insertion` · `--deletion-extent` · `--insertion-extent`, nucleotide only.

**Also needs saying somewhere user-facing:** `genome_summary.json`'s `block_events` gains a kind that
its `events` counterpart does not have. `insertion` is DNA and nothing else — it brings sequence
rather than a gene family — so the gene-level counters have no such kind. Every other kind still
appears in both. Appendix B already tells the reader those two counters do not bound each other;
this is one more way in which they do not.

---

## Noticed while editing, not this branch's to fix

`[Unreleased]` in `CHANGELOG.md` carries **two `### Added` headings** (they are both on `main`). Keep
a Changelog expects one section per type, so they want merging — but by whoever owns those entries,
not as a side effect of this work.
