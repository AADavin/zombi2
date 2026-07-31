# Coming from ZOMBI v1

ZOMBI2 is a rewrite, not a new release of ZOMBI. Nothing carries over automatically.

## The commands

v1 took a mode letter, a parameters file and an output directory:

```
python3 Zombi.py T ./Parameters/SpeciesTreeParameters.tsv ./Output_folder
python3 Zombi.py G ./Parameters/GenomeParameters.tsv     ./Output_folder
python3 Zombi.py S ./Parameters/SequenceParameters.tsv   ./Output_folder
```

v2 names the *run directory* once per command, and each level finds the previous one in it. v1's three
`*Parameters.tsv` files become one TOML file, with a table named for each command. A key outside any
table is shared by every command; the command's own table overrides it.

```toml
# zombi2.toml
seed = 1

[species]
birth    = 1.0
death    = 0.2
n-extant = 20

[genomes]
duplication = 0.1
transfer    = 0.2
loss        = 0.3
origination = 0.2

[sequences]
model      = "wag"
length     = 100
divergence = 0.2
```

```bash
zombi2 species   out/ --params zombi2.toml
zombi2 genomes   out/ --params zombi2.toml
zombi2 sequences out/ --params zombi2.toml
```

Every key is a long option of that command, hyphens or underscores alike, so the same run can be
written as flags with no file at all. A flag on the command line overrides the file.

```bash
zombi2 species   out/ --birth 1.0 --death 0.2 --n-extant 20 --seed 1
zombi2 genomes   out/ --duplication 0.1 --transfer 0.2 --loss 0.3 --origination 0.2 --seed 1
zombi2 sequences out/ --model wag --length 100 --divergence 0.2 --seed 1
```

There are no mode letters and no aliases for them: `zombi2 T` is an error naming the six commands. The
levels are `species`, `genomes`, `sequences`, `traits`, plus `joint` and `tools`.

## Parameters: species trees

Every rate takes the same written form in the flag, in the file and in Python.

| v1 key | v2 |
|---|---|
| `SPECIATION` | `--birth` |
| `EXTINCTION` | `--death` |
| `STOPPING_RULE` | gone as a key: pass `--total-time` (v1's rule 0) **or** `--n-extant` (rule 1), never both |
| `TOTAL_TIME` | `--total-time` |
| `TOTAL_LINEAGES` | `--n-extant` — and it is *conditioned on survival*: a run that dies out is retried |
| `MAX_LINEAGES` | `--max-lineages` (default 100000). It **raises** rather than truncating: a tree cut off at a size is not a sample from the process you asked for |
| `MIN_LINEAGES` | gone. v1 called a run a failure below it; v2 refuses only the degenerate case, a run with nothing alive at the present |
| `MASS_EXTINCTION` | `--mass-extinction TIME FRACTION`, repeatable |
| `SCALE_TREE` | `zombi2 tools tree TREE --rescale-height H` or `--rescale-factor F` |
| `VERBOSE` | `--quiet` (inverted) |
| `SEED` | `--seed` |
| `SHIFT_SPECIATION_RATE_FREQUENCY`, `NUM_SPECIATION_RATE_CATEGORIES`, `BASE_SPECIATION` (and the `EXTINCTION` trio) | the rate grammar: `--birth "1.0 * FromParent(spread=0.2)"` for a rate inherited and nudged at each split, `OnTime` for a schedule, `OnTotalDiversity` for a ceiling. v2 has no *category* count — the drift is continuous |
| `TURNOVER` | **no equivalent.** v2 parameterises `birth` and `death` directly |
| `LINEAGE_PROFILE` | **no direct equivalent.** The nearest thing is `--birth "1.0 * OnTotalDiversity(cap=N)"`, which levels diversity off at a ceiling |

v1 encoded a distribution choice as a prefix on the value (`f:1`, `l:1`, `n:0.7;0.2`). v2 writes it as
part of the rate instead, in one notation shared by the flag, the `--params` file and Python:

```python
from zombi2 import species
from zombi2.rates import modifiers as mod

species.simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=1)                 # fixed
species.simulate_species_tree(birth=1.0 * mod.FromParent(spread=0.2), n_extant=20, seed=1)  # drifting
```

## Parameters: genomes

| v1 key | v2 |
|---|---|
| `DUPLICATION` · `TRANSFER` · `LOSS` · `ORIGINATION` | `--duplication` · `--transfer` · `--loss` · `--origination`. D/T/L are **per copy**, origination **per lineage** |
| `INVERSION` · `TRANSPOSITION` | `--inversion` · `--transposition`, with `--resolution ordered` or `nucleotide` |
| `DUPLICATION_EXTENSION` and the other `*_EXTENSION` keys | `--duplication-extent` and friends, on `--resolution nucleotide` (in base pairs). At the ordered resolution they are Python-only, `duplication_extent=` |
| `REPLACEMENT_TRANSFER` | `--replacement` |
| `ASSORTATIVE_TRANSFER` | `--transfer-to distance` — nearest equivalent, not identical: v2 weights each candidate by `exp(-decay × d / depth)` |
| `INITIAL_GENOME_SIZE` | `--initial-families` (default 100) |
| `EVENTS_PER_BRANCH` | gone. `genome_events.tsv` has a `lineage` column; group on it |
| `PROFILES` · `GENE_TREES` | `--write profiles gene_trees` |
| `RECONCILED_TREES` | `zombi2 tools format DIR --format recphylo` — **recPhyloXML, not Newick**, so this needs a new parser |
| `RATE_FILE` · `SCALE_RATES` | the rate grammar again: `--loss "0.25 * OnTime({0: 1.0, 3: 2.0})"`, or `DrivenBy` to read another level |
| `GENE_LENGTH` · `INTERGENE_LENGTH` | `--gene-length` · the spacer is what lies between genes on `--resolution nucleotide` |
| `MIN_GENOME_SIZE` | **not a setting.** On `--resolution ordered` and `nucleotide` a loss never takes a chromosome below its last gene: the run that would empty it does not fire. On `--resolution family` there is no floor, and a high loss rate empties a genome completely |
| `ALPHA` | **no equivalent** as a genome parameter |
| `PSEUDOGENIZATION` | **not in v2** |

## Parameters: sequences

| v1 key | v2 |
|---|---|
| `SEQUENCE_SIZE` | `--length` |
| `AA_MODEL` | `--model wag` (also `jtt`, `dayhoff`, `lg`, `poisson`) |
| the `AC`…`TG` and `A`/`C`/`G`/`T` keys | `--model gtr --gtr-rates … --frequencies …`, or `hky85` / `k80` / `jc69` |
| `KAPPA` | `--kappa` |
| `ST_RATE_MULTIPLIERS` · `GF_RATE_MULTIPLIERS` | the lineage clock: `--substitution "1.0 * ByLineage(spread=0.3)"` (uncorrelated) or `FromParent` (autocorrelated) |
| `SHIFT_SUBSTITUTION_RATE` · `SHIFT_CATEGORIES` · `BASE_RATE` | the same clock. No category count |
| `SCALING` · `SCALE_GENE_TREES` | `--divergence D` sets the rate from the height of the tree, which is what scaling was for |
| `SEQUENCE codon` · `CODON_MODEL` | **not in v2.** There are no codon models, so no dN/dS |
| `ALPHA` · `BETA` | **not in v2.** No across-site rate variation (`+Γ`), no invariant-sites class |

## Where the files went

| v1 | v2 |
|---|---|
| `T/CompleteTree.nwk` | `species/species_complete.nwk` |
| `T/ExtantTree.nwk` | `species/species_extant.nwk` |
| `T/Events.tsv` | `species/species_events.tsv` |
| — | `species/species_fates.tsv` — new, and the only place that tells an *unsampled* tip from an extant one |
| `G/Events.tsv` | `genomes/genome_events.tsv` |
| `G/Profiles.tsv` | `genomes/profiles.tsv` |
| `G/Gene_families/<n>_completetree.nwk` | `genomes/gene_trees/gene_tree_fam<n>_complete.nwk` |
| `G/Gene_families/<n>_prunedtree.nwk` | `genomes/gene_trees/gene_tree_fam<n>_extant.nwk` — **"pruned" is now "extant"** |
| `G/<node>_GENOME.tsv` | `genomes/genomes.tsv`, one table for every node rather than a file each |
| `S/SubstitutionScaledCompleteTree.nwk` | `sequences/clock_species_tree_complete.nwk` |
| `S/<n>_substitution_scaled.nwk` | `sequences/phylograms/phylogram_fam<n>_complete.nwk` |
| — | `sequences/alignments/fam<n>.fasta` |
| — | `<level>.log` — version, timestamp, command line, every resolved parameter, and the SHA-256 of each input file |
| — | `<level>_summary.json` — what came out: deduplicated event counts, families born and surviving, whether the family-size cap bit |

`--flat` writes everything into one directory, but the v2 names stay: it is a layout switch, not a
compatibility mode.

## What will break in your parsing code

This is the part that takes the time. Six changes each break a script silently or loudly:

**Event vocabulary.** v1 wrote letter codes; v2 writes words — `speciation`, `extinction`,
`duplication`, `loss`, `origination`, and **two** transfer kinds, `transfer_additive` and
`transfer_replacing`.

**One row per event.** Every v2 event log opens with `time` and `kind` and then carries payload
columns that are the same on every row of that file. v1's genome log was positional, with a `NODES`
field whose meaning depended on the event type; `genome_events.tsv` is now five columns, and one row
is one thing that happened (real rows from a real run, padded here to line the columns up):

```
time                 kind                family  parents        children
0.0                  origination         0                      n0_g0
0.281867             speciation          14      n0_g14         n1_g15;n2_g29
0.2959603504961474   duplication         4       n1_g19         n1_g43;n1_g44
0.11415202408218499  loss                0       n0_g0
1.2866954826434076   transfer_additive   8       n1_g23         n1_g255;n5_g256
0.36393613420080373  transfer_replacing  1       n2_g30;n1_g16  n2_g47;n1_g48
```

The file itself is tab-separated. `parents` is what the event ended, `children` what it began, packed
with `;` where there are two; an origination has no parents and a loss no children, so those cells
are empty. Counting rows by kind now counts events, so a `grep -c` gives the number of events of that
kind. `species_events.tsv` and `chromosome_events.tsv` are the same shape, with their own entities in
the two cells. One exception: a nucleotide run's `block_events.tsv` is keyed by *ancestral interval*,
so an event spanning several blocks writes a row per block.

**Participants carry their own lineage.** A gene copy is written `n<species>_g<copy>`: `n2_g30` is
copy 30 on branch `n2`. So the `lineage`, `recipient` and `donor` columns are gone — split a token
on its single `_` and you have both halves. `species_events.tsv` names lineages the same way (`n0` →
`n1;n2`) and `chromosome_events.tsv` names chromosomes `n<species>_c<id>`. A transfer's `children`
read **donor first, recipient second**, so one row says which way the material went without pairing
anything.

**A replacing transfer writes no loss row.** `transfer_replacing` carries two parents — the
donor's copy, then the copy it overwrote on the recipient branch — and writes **no separate `loss`
row** for the copy it displaced. A script that counts losses straight out of the log therefore comes
up short by exactly the number of replacing transfers. `genome_summary.json` counts the biology
instead: `transfer` as one number over both kinds, each displaced copy under `loss`, and the
families the run started with as `initial` rather than `origination`. This is the change that gives
a plausible wrong number instead of an error, so it is the one to check first.

**Rearrangements are their own file.** An inversion, a transposition or a translocation begins and
ends no gene lineage, so it has nothing to put in `parents` and `children`. At the ordered and
nucleotide resolutions they are in `rearrangement_events.tsv`, not in `genome_events.tsv`, and that
file is the one place a branch is still a column (`lineage`) — a segment has no id to carry one.

**Label prefixes.** A lineage that went extinct is `e<id>`, not `n<id>`, in the complete tree and the
event logs — so a `n\d+` pattern silently drops the extinct lineages, inside a copy token
(`e6_g138`) as well as on its own. Gene copies are `g<id>`, and a gene-tree leaf is
`n<species>_g<copy>`, the same token the event log packs into `parents` and `children`. Everything
that names only extant tips — the extant tree, `profiles.tsv`, the alignments — is `n<id>`
throughout.

The full column list for every file is in [Output files](output-files.md).
