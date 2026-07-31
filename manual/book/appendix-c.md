<!-- --8<-- [start:intro] -->

# Tools

`zombi2 tools` reads a finished run and derives a new view of it. Each tool is a sub-subcommand,
`zombi2 tools <tool>`. `format` reads a whole genomes run and writes its derived files beside the run;
`tree` and `treedist` read one or two Newick files and write to stdout, or to a file with `-o`.

<!-- --8<-- [end:intro] -->

## `format` — analysis-ready files

<!-- --8<-- [start:format] -->

`zombi2 tools format DIR --format FORMAT` reads a genomes run and writes files derived from its gene
trees, beside the run under `genomes/`. `--format` is required and takes several at once:
`--format homology recphylo`. Family and ordered runs rebuild their gene trees from the event log; a
nucleotide run recovers them from the genome, one file per **declared gene** (the intergenic spacer is
not a gene, so it gets none). `--from PATH` reads a run elsewhere; `--flat` writes straight into the
output directory.

| Output | File | Format | Ask for | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Homology matrix | `homology_fam<f>.tsv` | TSV | `--format homology` | one **n×n** table per family (n the extant leaves), in `genomes/homology/`. Row and column headers are the leaves `n<species>_g<copy>`; each off-diagonal cell says how that pair diverged and whether transfer is in its history since — `S`, `D` or `T`, each optionally with an `x` (see below) — and the diagonal is `-`. Symmetric. A family with no surviving copy writes no table |
| Marker table | `markers.tsv` | TSV | `--format markers` | **one row per family** for the whole run, a single file at `genomes/markers.tsv`: is it single-copy, is it universal, and does its true tree match the species tree |
| recPhyloXML | `recphylo_fam<f>.xml` | XML | `--format recphylo` | one file per family, in `genomes/recphylo/`: that family's **complete** gene tree written inside the **complete** species tree, in the recPhyloXML format. Written for every family, extinct ones included |
| extant-only reconciliation | `recphylo_fam<f>_{true,recoverable}.xml`, `family_origins.tsv` | XML + TSV | `--format recphylo --recphylo extant` | the same history projected onto what a dataset holds — the extant gene tree inside the extant species tree — written twice, and a table saying how each family entered |

These files are exact, not inferred: ZOMBI2 **recorded** each gene tree's embedding in the species tree
rather than reconstructing it, which is what makes them an answer key to score a method against.

### The homology matrix

A cell carries two independent facts about a pair of genes.

**How they diverged** — the event at their most-recent common ancestor. `S` a **speciation**, one gene
whose species split in two; `D` a **duplication**, one gene became two inside one genome; `T` a
**transfer**, one gene became two in two different lineages.

**Whether transfer is in their history since** — an `x` suffix, when a transfer sits on the path from
that common ancestor down to *either* gene. A pair that diverged **at** a transfer is already `T`, so
`Tx` would mean a *second* one further down.

The smallest interesting case: species `((a, b), c)`, where `c` donates a copy of its gene into `a`, so
`a` ends up holding two.

```
          n3_g1    n4_g2    n2_g3    n3_g4      n3 = a,  n4 = b,  n2 = c
  n3_g1       -        S        S       Sx      g4 = the copy c sent to a
  n4_g2       S        -        S       Sx      g3 = the copy c kept
  n2_g3       S        S        -        T
  n3_g4      Sx       Sx        T        -
```

Every pair involving the arrival `g4` carries the `x`; no other pair does. `n3_g1` against `n3_g4` —
two genes in the same genome — reads `Sx`: they diverged when `a` and `c` split, and one came back by
transfer.

The table is read off each family's **complete** gene tree. The pairs and the letters are the same
either way, but the `x` is a fact about the path: in a pruned tree a transfer whose donor-side copy
left no surviving descendant is suppressed, taking the record of the transfer with it. On an ordinary
run that is a fifth of all cells.

A cell is the event, not a label for it. "Ortholog" and "paralog" are readings laid over the event and
the published definitions disagree, so ZOMBI2 reports the event and the definition stays yours. For
genes to build a species tree from, use the marker table.

### The marker table

`--format markers` asks a question about a **family**, not about a pair: can I put this one in a
concatenation and trust the tree that comes out? One row per family that left a surviving copy:

| Column | Meaning |
|---|---|
| `family` | the family id |
| `genomes` · `copies` | how many extant genomes carry it, and how many copies in total |
| `single_copy` | every genome that has it has exactly one — so there is no choosing which copy to align |
| `universal` | every extant genome has it (the criterion a BUSCO-style marker set is built on) |
| `duplications` · `transfers` · `losses` | the family's own history, dead lineages included |
| `rf` | Robinson–Foulds distance between the family's true gene tree, each gene read as the genome it sits in, and the species tree **restricted to those genomes**. Empty where it would mean nothing: several copies in one genome (no one-to-one gene→genome map) or fewer than three genomes (no clade to disagree about) |
| `congruent` | `rf` is 0 — the family recovers the species tree exactly |

The last two are what the table is for. **A family can be single-copy *and* universal and still give
the wrong tree** — a duplication followed by loss of the other copy in each descendant, or a transfer
that replaced the resident gene. That is hidden paralogy, and in real data it passes every filter you
would apply. On a transfer-rich run with replacing transfers, 111 of 299 families came out single-copy
and universal, and 106 of those did **not** recover the species tree.

```bash
# every family that would make a trustworthy marker
zombi2 tools format out/ --format markers
awk -F'\t' 'NR==1 || ($4=="yes" && $5=="yes" && $10=="yes")' out/genomes/markers.tsv
```

`rf` is the same distance `zombi2 tools treedist` reports, on the same rooted-clade convention.

### recPhyloXML

recPhyloXML [@duchemin2018recphyloxml] is the community format for a gene tree embedded in a species
tree: every gene-tree node carries the species branch it sat on and the event that ended it, so a
viewer can draw one inside the other. It is normally the output of a reconciliation *method*; here
nothing is reconstructed, so the file is the true history.

The **complete** gene tree goes inside the **complete** species tree, because the events the format
exists to show are the losses, and a gene that died leaves nothing in the extant tree to hang a
`<loss>` on. Extinct and unsampled species are kept for the same reason: a transfer can arrive from a
lineage that later died, and that edge has to land somewhere.

| ZOMBI2 | recPhyloXML |
|--------|-------------|
| duplication | `<duplication speciesLocation="n<species>">` |
| speciation | `<speciation speciesLocation="n<species>">` — the *parent* species, the branch the gene was on when its species split |
| loss | `<loss speciesLocation="n<species>">` |
| gene at an extant tip | `<leaf speciesLocation="n<species>">` |
| gene at an extinct or unsampled tip | `<leaf …>` as well — the species tree says which fate that branch had |
| transfer | the format's own two steps: `<branchingOut speciesLocation="n<donor>">` on the node the copy left from, and `<transferBack destinationSpecies="n<recipient>">` opening the child that arrived |

Origination has no tag and needs none: a family founded mid-branch is a gene tree whose root starts
there. Branch lengths are left out, as in the format's own reference files; the dated trees are next
door in `genomes/gene_trees/` and `species/species_complete.nwk`.

In Python, `zombi2.tools.recphylo.recphylo_xml(gene_trees, tree)` returns the document as a string —
hand it every family for one file a viewer can draw them all in, or one family for that family's own.

### The extant-only reconciliation

`--recphylo extant` writes the same history projected onto the extant gene tree and the extant species
tree, which is what a reconciliation method is scored against, because that is all it ever sees.
`--recphylo both` writes both scopes.

The projection keeps what is observable. A speciation where the gene followed one daughter becomes a
single loss on the other, however many losses really happened inside that clade; it disappears
entirely when the abandoned daughter has no surviving descendant. A duplication whose second copy died
disappears too.

Transfers are the interesting case. When a copy arrives from a lineage that leaves no survivor, the
transfer node vanishes with the donor. If that copy has surviving relatives, the point where it rejoins
them is written as a transfer out of the branch they share — the strongest claim the extant tree can
express. If it has none, the family appears to begin where it landed.

| | rooted at | what it is for |
|---|---|---|
| `recphylo_fam<f>_true.xml` | where the family really originated | the answer key, with the losses that narrowed it — what an ancestral gene-content reconstruction is trying to recover |
| `recphylo_fam<f>_recoverable.xml` | the surviving copies' common ancestor | the ceiling. A family that left no trace above that point cannot be placed higher by any method |

Both are written so the gap between them is visible: that gap is the part of the history a perfect
method still cannot reach.

`family_origins.tsv` sits beside them, one row per family, with `family`, `entered_by`, `branch` and
`losses`. `entered_by` is `origination` when the family really began on the branch it is rooted at, and
`transfer` when it arrived there from a lineage nobody can see; `losses` counts the loss leaves in the
`true` file.

Unsampled survivors are treated exactly as extinct lineages: nothing in the data refers to either.

<!-- --8<-- [end:format] -->

## `tree` — one transform on a Newick tree

<!-- --8<-- [start:tree] -->

`zombi2 tools tree TREE` applies a single transform to a Newick tree and writes the result — Newick to
stdout, or to a file with `-o`. Exactly one action runs per call. `TREE` is a tree file, or `-` for
stdin.

`--prune` needs each tip's fate, so it reads the fates a ZOMBI complete tree carries; an ultrametric
tree counts as all-extant, and a plain non-ultrametric tree with no fates is refused. Every other
action ignores fates and loads any tree.

| Action | What it writes |
|-----------------------------------|-------------------------------------------------|
| `--prune` | the extant tree: the dead and unsampled lineages dropped, and the unifurcations they leave behind suppressed so the tree stays bifurcating |
| `--round` | the tree snapped to exactly ultrametric, by extending the terminal branches to a common depth. `--tol` is the tolerance as a fraction of tree height (default `1e-3`); a wider tip-depth spread raises, because it is real tip-date signal, not rounding |
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

<!-- --8<-- [end:tree] -->

## `treedist` — distance between two trees

<!-- --8<-- [start:treedist] -->

`zombi2 tools treedist TREE_A TREE_B` reports how far apart two rooted trees are over their shared
tips, printed as `<metric><TAB><value>` to stdout (or a file with `-o`). Pick the metric with
`--metric`:

| `--metric` | Distance |
|-----------------|--------------------------------------------------------------------|
| `rf` (default) | Robinson–Foulds — the number of clades present in one tree but not the other |
| `rf-normalized` | that count over the total number of non-trivial clades, so it lands in `[0, 1]` |
| `branch-score` | Kuhner–Felsenstein — the square root of the summed squared branch-length differences over every clade, terminal branches included; unlike RF it moves even when only the branch lengths differ |
| `all` | every metric above, one per line |

Tips are matched by **label** — the tip name for an external tree, or `n<id>` for a ZOMBI tree — so a
true tree and an inferred tree line up by taxon whatever order their files list them in. The two trees
must carry the same tip set: a differing leaf set is an error, not a partial score.

```bash
# Robinson–Foulds between a true tree and an inferred one over the same tips
zombi2 tools treedist true.nwk inferred.nwk --metric rf
```

### Comparing a gene tree to a species tree

A gene tree's tips are genes (`n<species>_g<copy>`) and a species tree's are species (`n<species>`), so
the two share no labels. `treedist` notices and compares them on **the species each gene sits in** —
"does this family's tree recover the species tree?" — saying so on stderr. It works only when the
family is **single-copy**, so the mapping is one-to-one; a family with two copies in some genome is
refused, naming the genomes at fault, rather than answered with a plausible number. Two trees of the
same kind are left alone and nothing is printed.

```bash
zombi2 tools treedist out/genomes/gene_trees/gene_tree_fam3_extant.nwk \
                      out/species/species_extant.nwk --metric all
```

<!-- --8<-- [end:treedist] -->
