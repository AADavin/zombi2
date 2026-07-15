# Conventions

These are the conventions every ZOMBI2 model follows. They are what make a broad library
feel like **one tool** rather than a pile of scripts: the same names, the same outputs, the
same reproducibility guarantees, whichever level you simulate. A new model is not finished
until it follows them — see [Adding a model](adding-a-model.md).

The rule of thumb: **a user who has learned one command should be able to guess how the next
one behaves.** When in doubt, copy an existing model exactly.

## Naming

Names are assigned by the simulator, never invented per model.

| Kind | Prefix | Example | Meaning |
|---|---|---|---|
| Extant, sampled leaf | `n` | `n0`, `n12` | a species alive at the present and sampled |
| Internal node | `i` | `i0`, `i7` | a bifurcation |
| Extinct leaf | `e` | `e3` | a lineage that died before the present |
| Unsampled extant leaf | `u` | `u2` | alive at the present but not sampled (`sampling_fraction < 1`) |
| Gene lineage | `g` | `g0`, `g41` | one segment of a gene lineage between events (a gene-tree node) |
| Gene family | *(none)* | `1`, `2`, `37` | a family; a plain integer string |
| Nucleotide block | `block` | `+block0`, `-block3` | a syntenic segment, with strand sign |

Extinct (`e*`) and unsampled (`u*`) leaves are **distinct** — a mass extinction leaves `e*`
tips, incomplete sampling leaves `u*` tips — and both are kept in the complete tree so gene
families can transfer from them. User-supplied names from an input Newick are preserved.

## Reproducibility

Reproducibility is a first-class guarantee, not a nicety.

- **Every entry point takes `seed`.** Both the Python API and the CLI accept `--seed N` /
  `seed=N`. The API also accepts an explicit `rng=<numpy Generator>`.
- **One PRNG.** The pure-Python engines use `numpy.random.default_rng(seed)` (PCG64). The Rust
  engine takes the same integer seed; when only an `rng` is given, an integer seed is drawn
  from it and handed across.
- **Same seed → identical output.** For a fixed ZOMBI2 version, parameters, and seed, a run
  is byte-for-byte reproducible. Two documented exceptions: `--threads > 1` (Poisson-thinned
  parallel profiles are *statistically* identical, not byte-identical), and the Rust vs.
  pure-Python engines (statistically equivalent realizations, not bit-identical).

### The run manifest

**Every CLI command writes a run manifest** to the output directory — always, for
reproducibility (`species.log`, `genomes.log`, `traits.log`, `sequences.log`,
`coevolve.log`). It is a headed, tab-separated key–value file:

```
# ZOMBI2 run parameters
zombi2_version	0.2.0.dev0
timestamp	2026-07-06T14:23:45
command_line	zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o run/
age	5.0
birth	[1.0]
death	[0.3]
...
result	wrote run/species_tree.nwk (50 extant tips) in 0.021 s
```

The header fields (`zombi2_version`, `timestamp`, `command_line`) are followed by **every**
resolved argument, sorted, then a human-readable `result`. A model that adds a CLI flag gets
it recorded here for free — do not write a second, bespoke parameter file.

## Rates and time

- **Rates are per unit time.** Branch lengths are time; every rate is the parameter of an
  exponential (Gillespie) process in the same units.
- **Two words, and only two: `rate` and `modifier`.** A **rate** is a quantity *with units*
  (events, or substitutions/site, per unit time) — the base number you set. A **modifier** is a
  *dimensionless* multiplier (default 1) that scales a rate in some context (a family, a lineage, a
  lineage-pair, a site). Those are the only two "how fast" kinds anywhere the simulator runs forward.
- **Gene-family D/T/L are per gene copy** by default (`PerCopyRates`): the genome-wide rate is
  the per-copy rate times the number of copies. `PerLineageRates` instead fires at a constant
  per-lineage rate. `LineageRates` scales any base model by per-lineage factors.
- **Origination is per lineage**, not per copy — one new family per event, independent of
  genome size.
- State a model's rate semantics in its docstring in these terms (per copy / per lineage / per
  site), so a reader never has to guess what a number multiplies.
- **Odds, not rates, in the undated tools.** `zombi2 tools reconcile --model undated/reldated` and
  `zombi2 tools simulate` have no time, so their D/T/L parameters are dimensionless per-branch
  **odds** (each the ratio of that event to vertical descent), not rates. That is the one place the
  word "odds" applies; everywhere the simulator runs forward it is rates and modifiers.

## Outputs

- **Text, headed, tab-separated.** Tables are `.tsv` with a header row; trees are Newick
  (`.nwk`, timed unless a clock has rescaled them); sequences are FASTA.
- **Times are written `%.10g`.** Node times, event times.
- **The species tree always ships with its nodes.** Any command that emits a tree writes the
  `species_tree.nwk` + `species_nodes.tsv` pair (`name`, `time`, `is_leaf`, `is_extant`), so
  downstream tools can read node metadata without re-parsing Newick.

### Selecting output: `--write`

The `genomes` command writes a chosen subset of parts (`--write PART ...`); the default is
`profiles trees`.

| Part | Files |
|---|---|
| `profiles` | `profiles.tsv` + `presence.tsv` (or `profiles_sparse.tsv` with `--sparse`) |
| `trees` | `gene_trees/<family>_complete.nwk` + `_extant.nwk` |
| `trace` | `events_trace.tsv` (one compact, scalable event log) |
| `events` | `gene_family_events/<family>_events.tsv` (per-family detail) |
| `transfers` | `transfers.tsv` |
| `summary` | `gene_family_summary.tsv` |
| `branch_events` | `branch_events.tsv` (per-species-branch event counts, with `is_extant`) |
| `bed` | [nucleotide, genic] `genes.bed` + `bed/<node>.bed` (BED6 gene annotations) |
| `ancestral` | [nucleotide] `architecture/`, `genomes/<node>.fasta.gz`, `gene_alignments/` |

Representative schemas (quote these when adding a related output):

```
# profiles.tsv          family  n0  n1  n2 ...        (copy counts)
# events_trace.tsv      time  event  branch  donor  recipient  family  parent  child1  child2
# branch_events.tsv     branch  time  is_leaf  is_extant  origination  duplication  transfer_in  transfer_out  loss  inversion  transposition  total
# <node>.bed            chrom  chromStart  chromEnd  name  score  strand      (BED6, 0-based half-open)
# species_nodes.tsv     name  time  is_leaf  is_extant
```

Event codes are single letters: **O**rigination, **D**uplication, **T**ransfer, **L**oss,
**S**peciation.

## The Python API

- **Models are objects; simulators are functions.** You build a model (`BirthDeath(1.0, 0.3)`)
  and pass it to a `simulate_*` function (`simulate_species_tree(model, ...)`). No global
  state, no engine switch visible to the user.
- **Namespaces mirror the levels.** Everything is exported at the top level (`zombi2.BirthDeath`)
  and from a scikit-learn-style submodule (`zombi2.species.BirthDeath`) — the *same* objects.
  A new model is exported from both, and listed in `zombi2/__init__.py`'s `__all__`. The docs and
  examples use the **submodule form** (`from zombi2.species import BirthDeath, simulate_species_tree`)
  as canonical — it mirrors the levels and keeps `zombi2.<TAB>` legible as the catalog grows; the
  top-level alias stays for quick interactive use.
- **`seed=` / `rng=` on every simulator**, matching the reproducibility rule above.

## CLI grammar

The CLI is one subcommand per level: `species`, `genomes`, `traits`, `sequences`, `coevolve`.
Within a subcommand:

- A **`general`** argument group carries the shared flags: `--seed N` and the required
  `-o/--out DIR` (plus `-t/--tree FILE` where a tree is an input).
- **Model choice is a flag** (`--model`, `--mode`, `--diversification`, `--sse-model`, …), and
  each model's own parameters live in a **dedicated argument group whose description names the
  gating flag** — e.g. `add_argument_group("ornstein-uhlenbeck", "--model ou")`. This is what
  keeps `-h` legible as the library grows.
- **Output selection is `--write PART ...`** (and `--sparse` where profiles apply), never a
  bespoke per-command flag.
- Every subcommand ends with a runnable **`EXAMPLES`** epilog.

Flags are long and `--kebab-case`; a rate is named for the quantity it sets (`--dup`,
`--trans`, `--loss`, `--orig`). Match the neighbours.
