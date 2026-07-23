# Tools

`zombi2 tools` runs read-back analyses on a finished run. Where the level commands simulate, the tools
re-express what a run already wrote: they read its files and derive a new view of it. Each tool is a
sub-subcommand — `zombi2 tools <tool>` — and writes its output beside the run it read.

## `format` — analysis-ready tables

`zombi2 tools format DIR` reads a genomes run and writes tables derived from its gene trees, one
`--format` at a time, into a directory under `genomes/`. It works for every resolution: unordered and
ordered runs rebuild their gene trees from the event log; a nucleotide run recovers them from the
genome, and there one table is written per **declared gene** — the intergenic spacer is not a gene, so
it gets none. `--from PATH` reads a run that lives elsewhere; `--flat` writes the tables straight into
the output directory.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Homology matrix | `homology_fam<f>.tsv` | TSV | `--format homology` | one **n×n** table per family (n the extant leaves), in `genomes/homology/`. Row and column headers are the leaves `n<species>\|g<copy>`; each off-diagonal cell is the relation of that pair — `O` ortholog (their MRCA is a speciation), `P` paralog (a duplication), `X` xenolog (a transfer) — and the diagonal is `-`. Symmetric. A family with no surviving copy writes no table |

The homology matrix is exact, not inferred. ZOMBI simulated each gene tree's embedding in the species
tree, so the event at a leaf pair's most-recent common ancestor is **recorded** on the tree rather than
reconstructed from it: a speciation there makes the pair orthologs, a duplication paralogs, a transfer
xenologs. That is what makes these tables a ground-truth reference to score an orthology-inference
method against.

## The rest

The scoring and reconciliation commands (tree-distance, reconciliation-accuracy, the undated simulator,
the ALE likelihood) are quarantined during the clean-core rebuild; their files are documented here as
each returns.
