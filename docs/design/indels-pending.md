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

Still to write, as the remaining slices land: the legal-cut exemption (a deletion may fall inside a
gene), insertion, whatever the indel log is written to, the sequence level's gapped alignments, and
the CLI flags.

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
| Indel log | `.deletions` | list | Python | one `Deletion` per indel — `time` · `lineage` · `chromosome` · `deleted`, the last being `(copy, source, start, end)` rows in **ancestral** coordinates, the same payload a loss carries. Kept out of `events` because an indel ends no copy lineage and begins none, so the gene-tree recovery must not read it, and out of the root partition because cutting there would shatter it. Not written to a file yet, and `assembly()` — so a sequence run on top — refuses a genome that used it |

The `Default` column becomes a file rather than `Python` if the indel log is given one — an open
question in [`indels.md`](indels.md).

---

## Noticed while editing, not this branch's to fix

`[Unreleased]` in `CHANGELOG.md` carries **two `### Added` headings** (they are both on `main`). Keep
a Changelog expects one section per type, so they want merging — but by whoever owns those entries,
not as a side effect of this work.
