# Tools

`zombi2 tools` runs read-back analyses. Where the level commands simulate, the tools re-express what
already exists: they read files and derive a new view of them. Each tool is a sub-subcommand —
`zombi2 tools <tool>`. `format` reads a whole genomes run and writes its derived files beside the run. `tree`
and `treedist` work on Newick trees instead: they read one or two `.nwk` files and write their result
to stdout by default, or to a file with `-o`.

## `format` — analysis-ready files

`zombi2 tools format DIR` reads a genomes run and writes files derived from its gene trees, into a
directory under `genomes/` per `--format` (several at once is fine: `--format homology recphylo`). It
works for every resolution: family and ordered runs rebuild their gene trees from the event log; a
nucleotide run recovers them from the genome, and there it is one file per **declared gene** — the
intergenic spacer is not a gene, so it gets none. `--from PATH` reads a run that lives elsewhere;
`--flat` writes straight into the output directory.

| Output | File | Format | Default | Contents |
|-----------|-----------------|-------|-----|------------------------|
| Homology matrix | `homology_fam<f>.tsv` | TSV | `--format homology` | one **n×n** table per family (n the extant leaves), in `genomes/homology/`. Row and column headers are the leaves `n<species>_g<copy>`; each off-diagonal cell says how that pair diverged and whether transfer is in its history since — `S`, `D` or `T`, each optionally with an `x` (see below) — and the diagonal is `-`. Symmetric. A family with no surviving copy writes no table |
| Marker table | `markers.tsv` | TSV | `--format markers` | **one row per family** for the whole run, in `genomes/markers/`: is it single-copy, is it universal, and does its true tree match the species tree. The answer to "which families can I build a species tree from?" — see below |
| recPhyloXML | `recphylo_fam<f>.xml` | XML | `--format recphylo` | one file per family, in `genomes/recphylo/`: that family's **complete** gene tree written inside the **complete** species tree, in the recPhyloXML format. Written for every family, extinct ones included |
| extant-only reconciliation | `recphylo_fam<f>_{true,recoverable}.xml`, `family_origins.tsv` | XML + TSV | `--format recphylo --recphylo extant` | the same history projected onto what a dataset holds — the extant gene tree inside the extant species tree — written twice, and a table saying how each family entered |

Both are exact, not inferred. ZOMBI simulated each gene tree's embedding in the species tree, so
every event is **recorded** on the tree rather than reconstructed from it. That is what makes these
files a ground-truth reference to score an inference method against.

### The homology matrix

Two genes are related on **two independent axes**, and a cell carries both.

**How they diverged** — the event at their most-recent common ancestor:

- `S` — a **speciation**: one gene, and the species carrying it split in two.
- `D` — a **duplication**: one gene became two inside one genome.
- `T` — a **transfer**: one gene became two, in two different lineages.

**Whether transfer is in their history since** — an `x` suffix, when a transfer sits on the path from
that common ancestor down to *either* gene. A pair that diverged **at** a transfer is already `T`, so
it takes no suffix for that transfer; `Tx` would mean a *second* one further down.

So a cell reads `S`, `D` or `T`, each optionally with an `x`.

The table states **the event, not an interpretation of it**. ZOMBI knows the event exactly, because
it recorded the history rather than inferring it; "ortholog" and "paralog" are a reading laid over
that event, and — as the section below sets out — several published definitions read it differently.
Giving you the event lets you apply whichever definition you work with; giving you a verdict would
quietly apply ours.

Here is the smallest interesting case: species `((a, b), c)`, where `c` donates a copy of its gene
into `a`, so `a` ends up holding two.

```
          n3_g1    n4_g2    n2_g3    n3_g4      n3 = a,  n4 = b,  n2 = c
  n3_g1       -        S        S       Sx      g4 = the copy c sent to a
  n4_g2       S        -        S       Sx      g3 = the copy c kept
  n2_g3       S        S        -        T
  n3_g4      Sx       Sx        T        -
```

Every pair involving `g4`, the arrival, carries the `x`; no other pair does. The copy `c` **kept**
(`g3`) never went anywhere, so it relates to `a` and `b` by ordinary vertical descent — only the gene
that moved has transfer in its history. And `n3_g1` against `n3_g4` — two genes **in the same
genome** — reads `Sx`: they diverged when `a` and `c` split, and one of them came back by transfer.

The table is read off each family's **complete** gene tree, not its pruned one. The pairs and the
divergence letters are the same either way, but the `x` is a fact about the path: a transfer whose
donor-side copy left no surviving descendant becomes a degree-two node in the pruned tree and is
suppressed, taking the record of the transfer with it. On an ordinary run that is a fifth of all
cells.

### If you came here looking for orthologs

Reasonably, and this is the section for you — because ZOMBI2 deliberately does not print the word.

When people ask for "the orthologs" they usually want one of three different things:

1. **Genes to build a species tree from** — one copy per genome, with a history that *is* the species
   history. This is the phylogenomics use, and the one ZOMBI2 answers directly: it is the **marker
   table** below, not the homology matrix.
2. **The same gene in another species**, for carrying annotation across. That is a claim about
   *function*. ZOMBI2 does not model function, so it cannot answer it and does not pretend to.
3. **Related by descent without duplication** — Fitch's definition, which is what the homology matrix
   reports, spelled as the event rather than the word.

The reason for spelling it as the event is that the definitions genuinely disagree, and most of them
need something a pair of genes does not carry. Fitch's relation is not one-to-one: after a
duplication in one lineage, *both* copies are orthologs of the single gene in the sister lineage — its
co-orthologs — so "the ortholog of gene X" is usually a set, not a gene. Paralogy is relative to a
reference speciation: the same two genes are in-paralogs or out-paralogs depending on which split you
are asking about. The label even depends on the *model* you fit — reconcile a gene tree under
duplication and loss alone and there is no transfer category at all, so a history that was a transfer
is explained as a duplication plus losses, and a xenologous pair comes out paralogous. And the
graph-based definitions the widely used tools implement (reciprocal best hits, clustering) are
approximations to the phylogenetic one that are known to drift from it.

Every one of those is the event at the common ancestor **plus a choice** — a reference speciation, a
clade, a model, a threshold. ZOMBI2 gives you the event, which is the part it knows exactly and the
part that is genuinely a property of the pair; the choice stays yours.

### The marker table

`--format markers` answers the first question above, and it is a question about a **family**, not
about a pair: *can I put this one in a concatenation and trust the tree that comes out?* One row per
family that left a surviving copy:

| Column | Meaning |
|---|---|
| `family` | the family id |
| `genomes` · `copies` | how many extant genomes carry it, and how many copies in total |
| `single_copy` | every genome that has it has exactly one — so there is no choosing which copy to align |
| `universal` | every extant genome has it (the criterion a BUSCO-style marker set is built on) |
| `duplications` · `transfers` · `losses` | the family's own history, dead lineages included |
| `rf` | Robinson–Foulds distance between the family's true gene tree, each gene read as the genome it sits in, and the species tree **restricted to those genomes** — so a family present in half the tree is judged against the half it occupies. Empty where it would mean nothing: several copies in one genome (no one-to-one gene→genome map) or fewer than three genomes (no clade to disagree about) |
| `congruent` | `rf` is 0 — the family recovers the species tree exactly |

The last two are what the table is for. **A family can be single-copy *and* universal and still give
the wrong tree** — a duplication followed by loss of the other copy in each descendant, or a transfer
that replaced the resident gene. That is hidden paralogy, and in real data it is invisible: it passes
every filter you would apply and quietly poisons the concatenation. Here it is a column. On a
transfer-rich run with replacing transfers, 111 of 299 families came out single-copy and universal —
and 106 of those did **not** recover the species tree.

```bash
# every family that would make a trustworthy marker
zombi2 tools format out/ --format markers
awk -F'\t' 'NR==1 || ($4=="yes" && $5=="yes" && $10=="yes")' out/genomes/markers/markers.tsv
```

`rf` is the same distance `zombi2 tools treedist` reports, on the same rooted-clade convention, so
the two agree — checked on real families in the test suite.

### recPhyloXML

recPhyloXML [@duchemin2018recphyloxml] is the community format for a gene tree embedded in a species tree: every
gene-tree node carries the species branch it sat on and the event that ended it, so a viewer can draw
one inside the other. It is normally the output of a reconciliation *method*; here nothing is
reconstructed, so the file is the true history and can be used as the answer key.

The **complete** gene tree goes in, inside the **complete** species tree, because the events the
format exists to show are the losses — a gene that died leaves nothing in the extant tree to hang a
`<loss>` on. Extinct and unsampled species are in the species tree for the same reason: a transfer can
arrive from a lineage that later died, and that edge has to land somewhere. The mapping from ZOMBI2's
event log is one gene-tree node to one `<clade>`:

| ZOMBI2 | recPhyloXML |
|--------|-------------|
| duplication | `<duplication speciesLocation="n<species>">` |
| speciation | `<speciation speciesLocation="n<species>">` — the *parent* species, the branch the gene was on when its species split |
| loss | `<loss speciesLocation="n<species>">` |
| gene at an extant tip | `<leaf speciesLocation="n<species>">` |
| gene at an extinct or unsampled tip | `<leaf …>` as well — the gene reached the end of a species branch that happened to die; the species tree says which fate that branch had |
| transfer | the format's own two steps: `<branchingOut speciesLocation="n<donor>">` on the node the copy left from, and `<transferBack destinationSpecies="n<recipient>">` opening the child that arrived |

Origination has no tag and needs none: a family founded mid-branch is simply a gene tree whose root
starts there. Branch lengths are left out, as they are in the format's own reference files; the dated
trees are next door in `genomes/gene_trees/` and `species/species_complete.nwk`.

In Python, `zombi2.tools.recphylo.recphylo_xml(gene_trees, tree)` returns the document as a string —
hand it every family for a single file a viewer can draw all of them in at once, or one family for
that family's own file, which is what the command writes.

### The extant-only reconciliation

The file above is the whole simulated history, and that is not what a reconciliation method should be
scored against, because a method never sees it. It sees survivors: an extant gene tree, and an extant
species tree. `--recphylo extant` writes the same history projected onto those, and `--recphylo both`
writes all of it.

The projection keeps what is observable and drops what is not. A speciation where the gene followed
one daughter becomes a single loss on the other — however many losses really happened down inside
that clade, one is all anyone can infer, and one is what you get. The same speciation disappears
entirely when the abandoned daughter has no surviving descendant, because the extant species tree does
not contain that split at all. A duplication whose second copy died disappears too: both copies sat on
one branch, so the data shows one copy and no method could say otherwise.

Transfers are the interesting case. When a copy arrives from a lineage that leaves no survivor, the
donor is simply not in the extant tree, and the transfer node vanishes with it. What is left is a copy
that appears from nowhere. If that copy has surviving relatives, the point where it rejoins them is a
real divergence and is written as a transfer out of the branch they share — the gene demonstrably was
on that branch, and the copy left from somewhere below it. That is a weaker claim than naming the
branch it truly left, and it is the strongest one the extant tree can express. If the copy has no
surviving relatives, the whole family appears to begin where it landed; `family_origins.tsv` records
that, since the arrival is invisible in the reconciliation itself.

Two files come out per family, and the difference between them is the point:

| | rooted at | what it is for |
|---|---|---|
| `recphylo_fam<f>_true.xml` | where the family really originated | the answer key. A family present in an ancestor and surviving in a scattered few genomes is recorded as present in that ancestor, with the losses that narrowed it — which is exactly what an ancestral gene-content reconstruction is trying to recover |
| `recphylo_fam<f>_recoverable.xml` | the surviving copies' common ancestor | the ceiling. A family that left no trace above that point cannot be placed higher by any method, so nothing here can mark one wrong for missing something invisible |

`true` contains `recoverable`: trimming the ancestral presence back to the surviving copies and
dropping the losses that go with it is mechanical, and the other direction is impossible. Both are
written so you can grade either way, and so the gap between them is visible — that gap is precisely
the part of the history a perfect method still cannot reach.

`family_origins.tsv` sits beside them, one row per family: `entered_by` is `origination` when the
family really began on the branch it is rooted at, and `transfer` when it arrived there from a lineage
nobody can see. `losses` counts the loss leaves in the `true` file.

Unsampled survivors are treated exactly as extinct lineages here. A lineage alive but not sampled is,
for a reconciliation, in the same position as one that died: nothing in the data refers to it.

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

The tips are matched by **label** — the tip name for an external tree, or `n<id>` for a ZOMBI tree —
so a true tree and an inferred tree line up by taxon, whatever order their files list the tips in. The
two trees must carry the same tip set: a differing leaf set is an error, not a partial score.

### Comparing a gene tree to a species tree

A gene tree's tips are genes (`n<species>_g<copy>`) and a species tree's are species (`n<species>`),
so the two are not the same kind of object and share no labels at all. `treedist` notices, and
compares them on **the species each gene sits in** — the question "does this family's tree recover the
species tree?" — saying on stderr that it has done so, because a distance between two different kinds
of tree should never appear without a word:

```bash
zombi2 tools treedist out/genomes/gene_trees/gene_tree_fam3_extant.nwk \
                      out/species/species_extant.nwk --metric all
```

This only works when the family is **single-copy**: one gene per species, so the mapping is
one-to-one. A family with two copies in some genome has no well-defined distance to a species tree,
and that is refused, naming the genomes at fault, rather than answered with a plausible number.
Two trees of the same kind — two gene trees, or two species trees — are left alone and nothing is
printed.

```bash
# Robinson–Foulds between a true tree and an inferred one over the same tips
zombi2 tools treedist true.nwk inferred.nwk --metric rf

# every metric at once
zombi2 tools treedist true.nwk inferred.nwk --metric all
```
