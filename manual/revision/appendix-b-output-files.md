# Appendix B — Output files, in full

Every level's simulator returns a **Result** object — `SpeciesResult`, `GenomesResult`,
`SequencesResult`, `TraitsResult` — and every Result writes to disk the same way:

```python
result.write("out/")                              # the level's default outputs
result.write("out/", outputs=["events", "extant"])  # or choose exactly which
```

The **write vocabulary** — the tokens you pass to `outputs=` — is per level, and each token maps to one
file (or, at the sequence level, one file *per gene family*). This appendix is the single place every
file, the token that writes it, and its format are listed. What a run *can* write is gated by what it
**recorded**: under a narrowed `record=` (Chapter 3) some views are never built, and asking to write one
of them raises rather than inventing it.

Three formats appear throughout:

- **`.nwk`** — a Newick tree. Branch lengths are **time** (a *chronogram*) everywhere except the sequence
  phylograms, whose lengths are in **substitutions per site** (a *phylogram*).
- **`.tsv`** — a tab-separated table with a header row.
- **`.fasta`** — FASTA sequences.

And two tree flavours recur (Chapter 4): the **complete** tree carries every lineage, the extinct and
unsampled included; the **extant** tree is pruned to the observed survivors. Files are named
`…_complete` / `…_extant` to match.

---

## Species — `SpeciesResult`

Default: every applicable output.

| `outputs=` token | File | Format | What it holds |
|---|---|---|---|
| `complete` | `species_complete.nwk` | Newick (time) | every lineage that ever lived — extinct and unsampled included |
| `extant` | `species_extant.nwk` | Newick (time) | the observed survivors (omitted if none survived) |
| `events` | `species_events.tsv` | TSV | the true history: every speciation and extinction with its time (always recorded) |
| `fossils` | `species_fossils.tsv` | TSV | recovered fossils as `lineage · time` — present only if a fossil rate was set |

## Genomes — `GenomesResult` (unordered) · `OrderedGenomesResult` (ordered)

Unordered default: `events`, `profiles`. Ordered adds `gene_order`.

| `outputs=` token | File | Format | What it holds |
|---|---|---|---|
| `events` | `genome_events.tsv` | TSV | the gene-family event log — origination / duplication / transfer / loss — the source of truth every gene tree is derived from |
| `profiles` | `profiles.tsv` | TSV | the family × extant-species copy-count matrix (the comparative-genomics table) |
| `gene_order` | `gene_order.tsv` | TSV | *(ordered)* the observed genomes' layout: one row per gene — `species · chromosome · position · strand · family · gene` |
| `rearrangements` | `rearrangements.tsv` | TSV | *(ordered)* the inversion / transposition / translocation log |
| `chromosome_events` | `chromosome_events.tsv` | TSV | *(ordered)* the chromosome genealogy — fission / fusion / loss edges |

## Sequences — `SequencesResult`

Default: `alignments`, `phylograms`. The gene outputs are written **one file per gene family**
(`_fam<f>`); tips and nodes are keyed by their gene id `g<copy>`.

| `outputs=` token | File | Format | What it holds |
|---|---|---|---|
| `alignments` | `sequences_alignment_fam<f>.fasta` | FASTA | the observable gene alignment — the extant gene sequences (skipped for a family with no survivor) |
| `ancestral` | `sequences_ancestral_fam<f>.fasta` | FASTA | the reconstructed sequence at every internal node — the truth against which to score an ancestral-reconstruction method |
| `phylograms` | `sequences_phylogram_fam<f>_complete.nwk` · `…_extant.nwk` | Newick (**subs/site**) | the gene tree the sequences were drawn along; every node labelled `g<copy>`, so the tips pair with the alignment and the internal nodes with the ancestral sequences |
| `species_phylogram` | `sequences_species_phylogram_complete.nwk` · `…_extant.nwk` | Newick (**subs/site**) | the species tree scaled by the molecular clock — which lineages ran fast or slow (written only when the run was given a `GenomesResult`, so a species tree exists to scale) |

## Traits — `TraitsResult`

Default: `values`.

| `outputs=` token | File | Format | What it holds |
|---|---|---|---|
| `values` | `trait_values.tsv` | TSV | the trait value at each extant tip — the observable comparative-data vector |
| `changes` | `trait_changes.tsv` | TSV | *(discrete)* the realized transition log — `time · kind · lineage · from · to` — the trait twin of the genome event log |
| `tree` | `trait_tree.nwk` | Newick (time) | the complete tree with every node's trait value annotated |

---

## Coupling — no new files

Coupling adds no file formats; it reuses the level outputs above.

A **conditioned** run writes exactly what the target level's run writes, plus the **driver file** it read,
so the pairing that produced the pattern is kept on disk. A **joint** run writes **both** levels — the
grown species tree together with the trait history, or together with the genomes — each in the format it
would have had from its own command. Because a joint run grows the tree, the tree it writes is a
*complete* tree (Chapter 4), with the extinct lineages that shaped the distribution still in place.

---

*This appendix is the as-built truth: when a level gains an output, its row lands here.*
