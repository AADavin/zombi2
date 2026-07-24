# Tools

`zombi2 tools` runs read-back analyses. Where the level commands simulate, the tools re-express what
already exists: they read files and derive a new view of them. Each tool is a sub-subcommand —
`zombi2 tools <tool>`. `format` reads a whole genomes run and writes its tables beside the run. `tree`
and `treedist` work on Newick trees instead: they read one or two `.nwk` files and write their result
to stdout by default, or to a file with `-o`.

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

## `tree` — one transform on a Newick tree

`zombi2 tools tree TREE` applies a single transform to a Newick tree and writes the result — Newick to
stdout, or to a file with `-o`. Exactly one action runs per call. `TREE` is a tree file, or `-` to read
the tree from stdin.

The actions split by whether they need each tip's fate. `--prune` does: it drops the dead and unsampled
lineages to leave the extant tree, so it reads the fates a ZOMBI complete tree carries (an ultrametric
tree, whose tips are all contemporaneous, counts as all-extant). A plain non-ultrametric tree with no
fates is refused. The rest — the geometric transforms and `--red` — ignore fates and load any tree, so
an inferred phylogram or a rounding-noisy dated tree goes through them unchanged.

| Action | What it writes |
|-----------------------------------|-------------------------------------------------|
| `--prune` | the extant tree: the dead and unsampled lineages dropped, and the unifurcations they leave behind suppressed so the tree stays bifurcating |
| `--round` | the tree snapped to exactly ultrametric, by extending the terminal branches to a common depth. `--tol` is the tolerance as a fraction of tree height (default `1e-3`); a tip-depth spread wider than that raises, because it is real tip-date signal — extinct lineages or serial samples — not rounding |
| `--stem LEN` / `--stem-add LEN` | the branch above the crown set to `LEN`, or extended by `LEN`; nothing below the crown moves |
| `--rescale-height H` / `--rescale-factor F` | every branch length scaled — so the root-to-tip height becomes `H`, or by a raw multiplier `F` |
| `--red` | the RED-rescaled tree: node depths become their Relative Evolutionary Divergence (Parks et al. 2018), ultrametric on `[0, 1]` with the root at 0 and every tip at 1. `--red --values` writes a two-column `node<TAB>RED` table instead of a tree |

```bash
# drop the extinct lineages, extant tree to stdout
zombi2 tools tree out/species/species_complete.nwk --prune

# snap a rounding-noisy dated tree to ultrametric, to a file
zombi2 tools tree dated.nwk --round -o dated_ultrametric.nwk

# the RED of every node, as a table
zombi2 tools tree out/species/species_extant.nwk --red --values
```

## `treedist` — distance between two trees

`zombi2 tools treedist TREE_A TREE_B` reports how far apart two rooted trees are over their shared
tips, printed as `<metric><TAB><value>` to stdout (or a file with `-o`). Pick the metric with
`--metric`:

| `--metric` | Distance |
|-----------------|--------------------------------------------------------------------|
| `rf` (default) | Robinson–Foulds — the number of clades present in one tree but not the other |
| `rf-normalized` | that count over the total number of non-trivial clades, so it lands in `[0, 1]` |
| `branch-score` | Kuhner–Felsenstein — the square root of the summed squared branch-length differences over every clade, terminal branches included; unlike RF it moves even when only the branch lengths differ |
| `all` | every metric above, one per line |

The tips are matched by node identity, and a ZOMBI tree labels every node `n<id>`, so two trees from
the same run are matched by those labels. The two trees must carry the same tip set, identically
labelled: a differing leaf set is an error, not a partial score. This is a current limitation.

```bash
# Robinson–Foulds between a true tree and an inferred one over the same tips
zombi2 tools treedist true.nwk inferred.nwk --metric rf

# every metric at once
zombi2 tools treedist true.nwk inferred.nwk --metric all
```

## The rest

The remaining scoring and reconciliation commands (reconciliation-accuracy, the undated simulator, the
ALE likelihood) are quarantined during the clean-core rebuild; their files are documented here as each
returns.
