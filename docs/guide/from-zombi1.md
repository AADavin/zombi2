# Coming from ZOMBI v1

ZOMBI2 is a rewrite, not a new release of ZOMBI. Nothing carries over automatically: the commands are
different, the parameter files are a different format, every output path has moved, and the event logs
use words where v1 used letter codes. This page is the mapping, so that porting a pipeline is
mechanical rather than archaeological.

Two things to know before you start.

**A v1 seed reproduces nothing in v2.** The two are different programs drawing from different streams
in a different order. Anything you published from v1 has to be regenerated, not re-derived — see
*Reproducibility* in Chapter 2 for what a seed does and does not fix.

**Pin the version you port against.** ZOMBI2 is pre-1.0 and a minor bump may rename an output or a
keyword. `pip install zombi2==0.16.0` and the run logs (which record the version) are what make a port
stay ported.

## The commands

v1 took a mode letter, a parameters file and an output directory:

```
python3 Zombi.py T ./Parameters/SpeciesTreeParameters.tsv ./Output_folder
python3 Zombi.py G ./Parameters/GenomeParameters.tsv     ./Output_folder
python3 Zombi.py S ./Parameters/SequenceParameters.tsv   ./Output_folder
```

v2 names the *run directory* once per command, and each level finds the previous one in it. There is no
parameters file to pass, though `--params` takes one if you want it:

```bash
zombi2 species   out/ --birth 1.0 --death 0.2 --n-extant 20 --seed 1
zombi2 genomes   out/ --duplication 0.1 --transfer 0.2 --loss 0.3 --origination 0.2 --seed 1
zombi2 sequences out/ --model wag --length 100 --seed 1
```

There are no mode letters and no aliases for them: `zombi2 T` is an error naming the six commands. The
levels are `species`, `genomes`, `sequences`, `traits`, plus `joint` and `tools`.

## Parameters: species trees

Every long option is also a `--params` key (hyphens or underscores), and every rate takes the same
written form in the flag, in the file and in Python.

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
| `MIN_GENOME_SIZE` | **no equivalent.** A high loss rate can empty a genome completely, and nothing stops it |
| `ALPHA` | **no equivalent** as a genome parameter |
| `PSEUDOGENIZATION` | **not in v2** |

## Parameters: sequences

| v1 key | v2 |
|---|---|
| `SEQUENCE_SIZE` | `--length` |
| `AA_MODEL` | `--model wag` (also `jtt`, `dayhoff`, `lg`, `poisson`) |
| the `AC`…`TG` and `A`/`C`/`G`/`T` keys | `--model gtr --rates … --frequencies …`, or `hky85` / `k80` / `jc69` |
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

This is the part that takes the time. Four changes each break a script silently or loudly:

**Event vocabulary.** v1 wrote letter codes; v2 writes words — `speciation`, `extinction`,
`duplication`, `transfer`, `loss`, `origination`.

**Column names.** v1's genome log was positional, with a `NODES` field whose meaning depended on the
event type. v2 has nine named columns: `time`, `kind`, `lineage`, `family`, `copy`, `parent`,
`recipient`, `donor`, `event`.

**One row per gene-tree edge, not per event.** A duplication, a transfer and a speciation each write
**two** rows, so counting rows inflates them exactly 2×. Group on the `event` column, or read the
deduplicated counts out of `<level>_summary.json`. This is the change that gives a plausible wrong
number instead of an error, so it is the one to fix first.

**Label prefixes.** A lineage that went extinct is `e<id>`, not `n<id>`, in the complete tree and the
event logs — so a `n\d+` pattern silently drops the extinct lineages. Gene copies are `g<id>`, and a
gene-tree leaf is `n<species>_g<copy>`. Everything that names only extant tips — the extant tree,
`profiles.tsv`, the alignments — is `n<id>` throughout.

## What v2 does that v1 did not

Worth knowing before you decide the port is only a cost: incomplete sampling and fossil recovery as
first-class options; **traits**, and rates that read another level (`DrivenBy`) so that gene loss can
depend on a habitat; the `joint` command for when neither level can be grown first; a nucleotide
resolution that starts from a real annotation; per-family marker tables and exact homology tables
(`zombi2 tools format`); run logs with input digests; a summary per run; and a staleness guard that
refuses to leave one level disagreeing with another.
