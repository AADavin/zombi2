# Design: naming consolidation — one word, one concept

**Status:** proposed (2026-07-14); the four open forks were **ratified the same day** — see
[Decisions](#decisions-ratified). Extends the rate-clarity line
([rate-vocabulary](rate-vocabulary.md), [rate-modifiers](rate-modifiers.md)) from *rates* to the
**whole vocabulary**. Where `rate-vocabulary` fixed one axis, this note fixes the handful of core
words that still carry two or three meanings — `level`, `model`, `rate`/`speed`, and the residue of
half-finished renames — and pins the **mechanism** that lets us rename a published API without
breaking anyone.

This is a clarity refactor: **no behaviour changes, no RNG-stream changes.** Every rename is a pure
alias, the trick already proven for `SharedRates → PerCopyRates`.

## The problem in one sentence

ZOMBI2 is one coherent tool whose *organising words* have drifted: the same word (`level`, `model`,
`rate`) names two or three different things, and three good renames
(`output→write`, `branch→lineage`, `shared→per-copy`) landed halfway, leaving two vocabularies
coexisting. A user who learned the first meaning trips on the second.

## Governing principles

- **P1 — One word, one concept.** Each core noun names exactly one thing. `model` = a stochastic
  *process*; `level` = one of the four simulation *domains*; `resolution` = the genome sub-axis
  (unordered/ordered/nucleotide); `rate` = a quantity per unit time; `modifier` = a dimensionless
  multiplier on a rate; `lineage` = the tree-entity (per [rate-vocabulary](rate-vocabulary.md)).
- **P2 — Canonical + a marked deprecation window.** Every rename keeps the old name working, but
  **marked** — a `DeprecationWarning` and *out of the public listing* — for one minor version, then
  it is dropped. Silent co-equal aliases (today's state) are the thing we are removing. The
  mechanism is specified once, below, and reused by every consolidation.
- **P3 — Guessability.** *A user who has learned one command guesses how the next behaves*
  ([conventions](../contributing/conventions.md)). Every decision here is judged by that test.
- **P4 — Zero behaviour change.** Pure renames + docs. Same seed → byte-identical output, before and
  after, at every step.

## The deprecation mechanism (specified once, reused everywhere)

Today's aliases are *silent and co-equal*: `SharedRates` is in `__all__`, `help(z.SharedRates)`
shows the `PerCopyRates` docstring with no notice, and `dir(zombi2)` lists both. That is what makes
"which is canonical?" unanswerable. The fix has three forms, one per surface:

- **Python names** — PEP-562 module `__getattr__`. Drop the deprecated name from `__all__` (so it
  leaves `dir()` and the mkdocstrings API reference), but resolve it in `__getattr__` with a
  `DeprecationWarning` that points at the canonical name. Old imports keep working; they warn once
  and vanish from the catalog.

  ```python
  _DEPRECATED = {"SharedRates": "PerCopyRates", "BranchRates": "LineageRates", ...}
  def __getattr__(name):                      # PEP 562
      if name in _DEPRECATED:
          new = _DEPRECATED[name]
          warnings.warn(f"zombi2.{name} was renamed to zombi2.{new}; "
                        f"the old name is removed in 0.4.", DeprecationWarning, stacklevel=2)
          return globals()[new]
      raise AttributeError(name)
  ```

- **CLI flags & value tokens** — keep them accepted but `help=argparse.SUPPRESS` (gone from `-h`),
  and print one stderr line on use: `zombi2: warning: --branch-rates is deprecated; use
  --lineage-rates`. A single `_deprecated_flag(old, new)` helper does this; several commands already
  half-do it.

- **Output filenames** — the one surface with *no* transparent alias (a file has one name). These
  are therefore batched into a **single** "output filenames v2" break (§C7), shipped once, with a
  one-version opt-in `--legacy-filenames` shim and a prominent CHANGELOG migration note — so
  downstream scripts absorb *one* rename event, not seven scattered ones.

**Version plan.** Aliases + deprecation warnings land in **0.3.0**; removal in **0.4.0**. The
output-filenames break lands in **0.3.0** behind the shim, shim removed in **0.4.0**.

---

# The consolidations

Ordered by clarity-per-effort. Each is an independently shippable PR.

## C1 — Retire the alias sprawl from the public surface *(ship first)*

**Confusion.** The public API exports both halves of five renames co-equally, with no runtime
signal: `PerCopyRates`+`SharedRates`, `PerLineageRates`+`PerGenomeRates`, `LineageRates`+`BranchRates`,
`LineageModifier`+`BranchModifier`, `read_lineage_rates`+`read_branch_rates` — all in `__all__`
([`zombi2/__init__.py:111`](../../zombi2/__init__.py)). On the CLI, `--rate-model`, the `--rate-per genome`
token, `--branch-rates`, and `--initial-chromosomes` are all still shown in `-h`; worse,
`--initial-chromosomes` *overrides* its own replacement `--n-chromosomes` when both are given.

**Decision.** Apply the P2 mechanism to what is *already renamed*: move the five alias names out of
`__all__` into the `__getattr__` deprecation map; `SUPPRESS` the deprecated CLI flags/tokens and warn
on use; fix the `--initial-chromosomes` precedence so the canonical flag wins (or errors on
conflict, matching how `--rate-per` + `--rate-model` already error).

**Back-compat.** Full. Every old name/flag keeps working; it just warns and leaves the catalog.
**Effort.** Small, mechanical. This PR also *builds* the reusable `__getattr__` map and
`_deprecated_flag` helper that C2–C9 reuse.

## C2 — Finish `branch → lineage` (decision already ratified)

**Confusion.** [rate-vocabulary](rate-vocabulary.md) *accepted* "retire branch, standardise on
lineage," but the rename only reached the genomes rate-file flag. Residue that still says *branch*:

| Surface | Today | → |
|---|---|---|
| sequence clock flags | `--branch-speed`, `--branch-bins`, `--branch-switch-rate`, `--branch-up-bias` | see note below |
| output part / file | `branch_events` / `Branch_events.tsv` | `lineage_events` / `Lineage_events.tsv` (§C7) |
| sequence output | `branch_rates.tsv` | `lineage_rates.tsv` (§C7) |
| API field | `GenePhylograms.branch_rate` | `lineage_rate` (alias) |
| **conventions.md itself** | lines 71–75, 99, 108: `PerGenomeRates`, "Origination is per branch", `BranchRates`, `branch` schema | rewrite to `lineage` |

The `conventions.md` rows are the tell: even the canonical doc still teaches the retired word.

**Decision.** Sweep *branch* out of every rate/lineage context. The four `--branch-*` **clock**
flags are legacy shorthands for `--clock autocorrelated-lognormal --clock-sigma …`
(see [`cli.py:2977`](../../zombi2/cli.py)); rather than rename them to `--lineage-*` we **deprecate
them toward the explicit `--clock` interface**, which is already `choices`-validated and modern —
killing the *branch* residue and a redundant parallel interface in one move. Reserve *branch* for
exactly one thing, unrelated to rates: `--branch-order` in `tools treedist` (the L_p order of the
branch-score distance), which is correct as-is.

**Back-compat.** Flags/fields: aliased + warned. Filenames: batched into §C7.
**Effort.** Small (mostly a prose + alias sweep).

## C3 — Two core terms (`rate`, `modifier`); `odds` is local to undated ALE

**Confusion.** The vocabulary drifted to *three* words for what, across the whole forward simulator,
is *two* things. "Rate" and "speed" are used as synonyms — `--subst-rate` (a base
substitutions/site/time) vs `--family-speed` (a multiplier) — and "modifier"
([rate-modifiers](rate-modifiers.md)) is *also* a multiplier. So "speed" is a redundant third word:
it means exactly what "modifier" already means.

**Decision — the core vocabulary is two terms, and it is dimensional.** Everywhere ZOMBI2 simulates
forward — species, genomes, traits, sequences — a "how fast" number is one of exactly two kinds
(new `conventions.md` section):

> **rate** = a quantity *with units* (events, or substitutions/site, *per unit time*) — the base
> number you set. **modifier** = a *dimensionless* multiplier, default 1, that scales a rate in a
> context (a family, a lineage, a lineage-pair, a site).

Dimensioned or dimensionless — there is no third dimension for a tempo number to occupy, which is
why two terms suffice. (*How-many* is the separate opportunity/count axis already named in
[rate-vocabulary](rate-vocabulary.md) — *rate = base × opportunities × modifiers* — a count, not a
synonym for a rate.)

**Retire "speed" — flag renames DEFERRED (needs a human decision).** On closer inspection the
substitution "speed" flags are subtler than "speed = modifier": `--family-speed` is the **σ of the
per-family multiplier distribution** (each family draws `~ LogNormal(0, σ)`), `--family-speeds` is a
*file* of per-family multipliers, and `--branch-speed` is the autocorrelated **clock σ**. So a
`--family-speeds` → `--family-modifiers` file rename is clean, but `--family-speed` / `--branch-speed`
are *spread* parameters, not modifiers — renaming them to `--*-modifier` would misdescribe a σ. These
five clock/spread flags (`--family-speed`, `--branch-speed`, `--branch-bins`, `--branch-switch-rate`,
`--branch-up-bias`) are therefore left **as-is** (byte-identical) pending Adrián's naming call — they
also carry the last of the `branch→lineage` residue on the CLI. Everything else in C3 ships: the
two-term rule, the odds scoping, and the docs. `--subst-rate` stays a genuine **rate**.

**One field exception.** `GammaRates` — among-**site** +Γ rate heterogeneity — keeps its name even
though such "rates" are technically mean-1 modifiers: *rate heterogeneity* / *+G rates* is bedrock
phylogenetics vocabulary, so field convention outranks the internal rule, exactly as the SSE
acronyms do in C8.

### Scoped sidebar: `odds` in the undated-ALE corner (not core)

One place is genuinely *not* rate-and-modifier, and it is small and self-contained: the **undated ALE
model** (`tools reconcile --model undated/reldated`, `tools simulate`). It has **no time**, so its
δ/τ/λ are not rates. Per Williams, Davín et al. 2023 (GBE), with vertical descent fixed to 1,
`P_D = δ/(1+δ+τ+λ)`, so **δ = P_D / P_S** — the *odds* of duplication against vertical descent:
dimensionless, unbounded, a ratio-to-reference. (The bounded, normalized `P_D … P_S ∈ [0,1]` are the
derived **probabilities** — real, but not the parameter a user sets.)

So **`odds` is introduced only where it lives** — a one-paragraph note in the undated-tools docs, and
their `--dup/--trans/--loss` help calling the numbers "per-branch odds," not "rates." It is **not**
added to the core primer or `conventions.md`'s main rule; the forward simulator never sees it. The
one consistency fix that rides along in this PR: those two tools currently call δτλ "rates" in
places — the single rate/odds conflation to correct. (Conceptually `odds` is the time-marginalised
shadow of a `rate`; `reldated` is the hybrid — but that is a sidebar, not the headline.)

**Effort.** Small (one conventions section + rename `--family-speed` with an alias + fix the two
undated tools' help text).

## C4 — `model` means a process; the genome axis is a *resolution*, not a "level"

**Confusion.** `--genome-model {unordered, ordered, nucleotide}` selects a genome *representation*
(its own metavar is even `LEVEL`), while every other `-model` flag selects a stochastic *process*
(`--subst-model`, `--sse-model`, `--omega-model`, `--trait-model`, `--score-model`). Meanwhile three
real process-selectors *drop* the suffix (`--clock`, `--diversification`, `--critic`). And internally
`dest="model"` is written by three unrelated flags — species `--mode` (`backward`/`forward`!), trait
`--model`, coevolve `--trait-model`.

**Decision.** The axis word must agree with C6, which reserves "level" for the four domains — so the
genome axis is **not** a "level":
- Rename `--genome-model → --genome-resolution` (alias kept, warned; metavar `RESOLUTION`). The
  values are a coarse→fine ladder — gene families ⊃ gene order ⊃ nucleotides — so "resolution" is
  literal and ties to the granularity ladder in [rate-vocabulary](rate-vocabulary.md). *(Word still
  open: `--genome-representation` / `--genome-type` are alternatives; whichever is chosen is used
  identically in the flag, the metavar `RESOLUTION`/…, and C6's prose — one word, three surfaces.)*
- Reserve "model" for a stochastic process: *a `--…-model` flag names a process and its `choices` are
  model names.* Fix the internal `dest="model"` collisions — species `--mode` stores `args.mode`;
  each `--*-model` gets its own dest. Internal only, no back-compat needed.

**Effort.** Small.

## C5 — CLI commands are plural nouns

**Confusion.** Commands mix number: `species`, `genomes` (plural) but `trait`, `sequence`
(singular). The packages *and* every guide page are uniformly plural (`zombi2.traits` / "Traits",
`zombi2.sequences` / "Sequences"). So `zombi2 trait` ↔ `zombi2.traits` ↔ guide "Traits" — three
forms for one level. `conventions.md:131` currently *enshrines* the singular.

**Decision.** Rename `trait → traits`, `sequence → sequences` (singular kept as a `SUPPRESS`-ed
accepted alias for the window). Now command == package == guide page == "level", everywhere. Leave
`coevolve` a verb: it is inherently relational (it couples two levels), and "coevolution" as a
command reads worse than the action. Update `conventions.md:131`.

**Back-compat.** Old singular commands still run (warned). **Effort.** Medium — a docs/examples
sweep (mechanical), plus every `_examples(...)` epilog.

## C6 — One domain word, and `level` reserved for the domain axis

**Confusion.** Two collisions on the tool's top-level metaphor:

1. **"Level" names two different axes.** The *domain* axis — "the four levels" = species, genomes,
   traits, sequences ([index.md](../index.md), [rates.md](../guide/rates.md)) — and the *resolution*
   axis inside genomes — "three levels" = unordered/ordered/nucleotide ([genomes.md](../guide/genomes.md)).
2. **The same domain has two names.** Coevolve's grammar is over `{species, traits, genes}`
   ([coevolution.md](../guide/coevolution.md)) — but "genes" *is* the domain the rest of the tool
   calls "genomes". A reader who learned "the four levels" meets "the three levels {species, traits,
   genes}" and must silently reconcile a rename (genomes→genes) and a disappearance (sequences).

**Decision.**
- **Reserve "level" for the four domains.** The genome sub-axis (unordered / ordered / nucleotide) is
  a **resolution** of the genomes level, never a "level" itself — named consistently in prose and in
  the `--genome-resolution` flag (C4). It is a granularity ladder (gene families ⊃ gene order ⊃
  nucleotides), not a fifth domain.
- **One domain word: `genomes`, not `genes`** *(ratified D1)*. Make the coevolve node name
  `genomes` too (`--couple traits:genomes`), with `genes` an accepted (warned) alias, so the domain
  word is identical everywhere. State explicitly in the coevolve guide: *the coevolve levels are a
  subset of the four — {species, genomes, traits} — because a coupling needs a driver and a target
  with rates; sequences are downstream and do not drive.* That single sentence closes the "why three
  not four" gap.
- **Level vs. unit.** `genomes` is the *level/node* word; **`gene family` stays the *unit* word.** So
  coevolve class names that describe the coupled *unit* (`GeneDiversification` = gene-family-content
  drives the tree, `TraitGeneCoupling` = trait↔gene-family) keep "Gene" and are *not* a violation —
  they name a gene family, not the genomes level. The `--couple` node token is the level; the class
  stem is the unit. C8 only unifies the *stem within* the `traits:genomes` edge, it does not
  genome-ify the unit names.

**Effort.** Small–medium (the node-name alias + a targeted doc reconciliation).

## C7 — Output filenames v2 (one break, one casing rule)

**Confusion.** Two casing conventions ship with *no rule*, and for two artifacts **both spellings
exist today**:

- `lowercase_snake` (species/trait/coevolve provenance) vs `Capitalized_Snake` (matrix/reconciliation/tools).
- Both ship: `gene_trees/` **and** `Gene_trees/`; `alignments/` **and** `Gene_alignments/`.
- `species_tree.log` is the lone log named for the artifact (all others are `<command>.log`).
- Arbitrary suffixes: `Events_trace.tsv` vs `Geneorder_events.tsv` both serialise the event log
  (`_trace` vs `_events`); `Reconciled_*.nwk` (adjective) vs `Reconciliation_*.tsv` (noun) in one
  trio; `genes.tsv`/`genes.bed` vs `Genes.gff`.
- `experimental ils` writes `gene_trees/` (a **dir**) in one branch and `gene_trees.nwk` (a **file**)
  in another.
- Opaque names: `RED.tsv`, `Mosaics.tsv`, `coupling.tsv`.

**Decision** *(ratified D2)*. **One casing rule: `lowercase_snake.ext` for every generated table,
tree, and directory** — the more standard *nix convention, and the one the most-referenced file
already uses (`species_tree.nwk`), so that pair needs no grandfathering. Applied in a single
versioned break behind the `--legacy-filenames` shim (P2). The visible churn falls on the currently
`Capitalized_Snake` products (the larger set), which is the deliberate cost of consistency.

| Fix | From → To |
|---|---|
| all data products lowercase | `Profiles.tsv`→`profiles.tsv`, `Presence.tsv`→`presence.tsv`, `Profiles_sparse.tsv`→`profiles_sparse.tsv`, `Events_trace.tsv`→`events_trace.tsv`, `Transfers.tsv`→`transfers.tsv`, `Gene_family_summary.tsv`→`gene_family_summary.tsv`, `Reconciled_*.nwk`→`reconciled_*.nwk`, `Reconciliation_*.tsv`→`reconciliation_*.tsv`, `Tree_distances.tsv`→`tree_distances.tsv`, `Karyotype_trace.tsv`→`karyotype_trace.tsv`, `Chromosomes.tsv`→`chromosomes.tsv`, `Mosaics.tsv`→`mosaics.tsv`, `Pseudogenizations.tsv`→`pseudogenizations.tsv`, `Geneorder_events.tsv`→`geneorder_events.tsv`, `Breakpoints.tsv`→`breakpoints.tsv`, `Positional_orthologs.tsv`→`positional_orthologs.tsv`, `Gene_family_profiles.tsv`→`gene_family_profiles.tsv` |
| all output dirs lowercase | `Gene_trees/`→`gene_trees/`, `Intergene_trees/`→`intergene_trees/`, `Gene_alignments/`→`gene_alignments/`, `Architecture/`→`architecture/`, `Genomes/`→`genomes/`, `BED/`→`bed/` |
| collisions (both ship) resolve *downward* | `gene_trees/`+`Gene_trees/`→`gene_trees/`; `alignments/`+`Gene_alignments/`+`gene_alignments/`→`gene_alignments/` — the lowercase choice fixes these for free |
| `genes.*` trio consistent for free | `Genes.gff`→`genes.gff` (now matches the already-lowercase `genes.tsv`/`genes.bed`) |
| log naming | `species_tree.log`→`species.log` (match `<command>.log`) |
| dir/file collision | `experimental ils` single-file `gene_trees.nwk`→`ils_gene_trees.nwk` (so it doesn't shadow the `gene_trees/` dir) |
| acronym | `RED.tsv`→`red.tsv` (lowercase rule wins; a header comment line names the acronym so it stays self-describing) |
| pair grouping | `presence.tsv` documented as the partner of `profiles.tsv` |
| opaque names | `mosaics.tsv`, `coupling.tsv` get a one-line self-describing header |
| branch→lineage (C2) | `Branch_events.tsv`→`lineage_events.tsv`; `branch_rates.tsv`→`lineage_rates.tsv` |

This also rewrites the `Capitalized` schema names in `conventions.md` (§Outputs table, lines
91–110). And it closes a gap the same convention already mandates: `tools` and the nucleotide
`genomes` path currently write **no** run-manifest log (`conventions.md` says *every* command must).

**Effort.** Medium. It touches many string literals but is mechanical, and the shim contains the
blast radius to one version.

## C8 — Coevolve classes follow the `driver:target` grammar

**Confusion.** The `--couple driver:target` grammar is elegant; the classes behind it are not
guessable from it.

- **`traits:genes` has two naming stems**: config `TraitGeneCoupling` vs engine/function/result
  `TraitLinkedRates` / `simulate_trait_linked_genomes` / `TraitLinkedResult` — both public.
- **`species:traits` is orphaned** as `Cladogenesis`, living in `zombi2.traits`, run via
  `simulate_traits`, with no `*Result`.
- **`species:genes`** (`CladogeneticGenome`) is muddied by the joint: [`__init__.py:149`](../../zombi2/__init__.py)
  says "species:genes = the co-diversification joint model", conflating the edge with
  `simulate_co_diversification`.
- SSE (`BiSSE`/`MuSSE`/…) doesn't map to the grammar — **but this one is left alone**: they are
  field-standard names, and forcing them into `traits:species` vocabulary would be *worse*. Document
  that SSE *is* the `traits:species` edge; don't rename it.

**Decision.**
- Collapse `traits:genes` to one stem: rename `TraitLinkedRates`/`…Result`/`simulate_trait_linked_genomes`
  to the `TraitGene*` stem (`TraitGeneRates`, `TraitGeneResult`, `simulate_trait_gene_genomes` — or
  keep `simulate_trait_linked_genomes` as the documented function and just align the *types*). Old
  names aliased.
- Make the **edge→class table** in the coevolve guide the canonical map, and cross-reference
  `Cladogenesis` from it so the `species:traits` edge is discoverable despite living in `traits`.
- Fix the `__init__.py:149` comment to distinguish the pure `species:genes` edge
  (`CladogeneticGenome`) from the `species↔genes` joint (`simulate_co_diversification`).

**Effort.** Small–medium (type aliases + a doc table).

## C9 — "Selection" disambiguated by *removing* the experimental tier

**Confusion.** The core codon models *are* selection — ω / dN/dS, "purifying/positive selection"
throughout [`sequences/codon_models.py`](../../zombi2/sequences/codon_models.py). The *experimental*
`selection` command is a different mechanism (ESM2 / a protein language model). They collide hard
enough that [`sequences/__init__.py:41`](../../zombi2/sequences/__init__.py) carries a comment
written only to dodge a `codon_selection.translate` name clash.

**Decision** *(ratified D3)*. **Remove the experimental protein-language-model selection family from
the published repo** — purge, not rename. The active PLM work already lives in the private
PLM-realism workspace; the code moves there (and stays in git history), where it can be revived as a
proper Extension. The collision then **resolves by removal**: "selection" becomes unambiguous — it
means the core codon ω / dN-dS models, full stop.

Scope (the ESM family; `experimental ils` is *not* ESM and stays):

- Delete `zombi2/experimental/{selection,codon_selection,genome_selection,nucleotide_selection,realism}.py`
  and their re-exports from `experimental/__init__.py` (`Critic`, `ESM2Critic`, `FixedProfileCritic`,
  `PLMSelection`, `CodonSelection`, `calibrate_beta`, `translate`, `CDS`, `GenomeSelection`,
  `read_cds_gff`, `frechet_esm_distance`, `BlockSelectionReport`, `NucleotideGenomeSelection`,
  `simulate_nucleotide_selection`).
- Remove the `experimental selection` CLI subcommand — `experimental ils` becomes the *sole*
  experimental model (simplifies that command's help to one entry).
- Drop the `zombi2[selection]` extra and its `torch`/`esm`/`fair-esm` deps from `pyproject.toml`.
- Remove the ESM tests (`test_selection.py`, `test_selection_live.py`, `test_codon_selection.py`,
  `test_nucleotide_selection.py`, `test_genome_selection*`, `test_realism.py`) and the docs
  (`docs/experimental/*selection*`, the ESM report under `CLAUDE/ZOMBI2_EMBEDDINGS_REPORT/` is
  already outside the repo).
- Delete the now-moot name-clash comment at `sequences/__init__.py:41`.

Because this *removes a shipped feature* (it was promoted in PR #100), it lands as its own PR with a
clear CHANGELOG "removed" entry, and the code is preserved in the private workspace first.

**Destination (decided).** The five modules **move into the private plm-realism workspace**
(`CLAUDE/ZOMBI2_workspaces/plm-realism/`) *first*, then are deleted from the repo — the live code is
preserved off-GitHub, not left to git history alone.

**Effort.** Small–medium (mechanical deletion + dep/test/doc cleanup).

## C-defer — Result-type regularisation (lower priority, do opportunistically)

The `simulate_*` return types use four conventions — `<Level>Result` (`GenomeResult`,
`NucleotideResult`, `TraitResult`, the coevolve `*Result`s), bare plural (`Genomes`,
`GenePhylograms`), descriptive (`RateScaledTree`), and reuse (`simulate_sse → TraitResult`). This is
real but low-traffic. **Cheap win now:** demote the engine-only `GenomeResult` from `__all__` (it is
a low-level dataclass redundant with the user-facing `Genomes`; users never construct it) and
publish a one-row-per-level "what `simulate_*` returns" table in the API reference so no one guesses.
**Full regularisation** (picking one convention and mass-renaming) is deferred — the payoff per rename
is small and the churn touches every level. Flagged so it is not forgotten, not scheduled.

---

# Decisions (ratified)

The four open forks, resolved 2026-07-14:

| # | Fork | Ratified outcome |
|---|---|---|
| D1 (C6) | coevolve domain node name | **`genomes` everywhere** — `--couple traits:genomes`; `genes` kept as a warned alias. One domain word tool-wide. |
| D2 (C7) | output-file casing | **`lowercase_snake` for all generated files** — no grandfathering; the `Capitalized` data products re-case. |
| D3 (C9) | experimental ESM selection | **Purge from the repo**, not rename — the PLM/ESM family moves to the private workspace; `experimental ils` remains the sole experimental model. |
| D4 (C5) | `coevolve` command | **Left as the one verb** (it is inherently relational). |

# Rollout

1. **PR 1 (C1)** — build the `__getattr__` map + `_deprecated_flag` helper; retire the existing five
   aliases + deprecated flags from the public listing. *No new renames yet — just the mechanism and
   the existing debt.* Independent of every ratified fork; can start immediately.
2. **PR 2 (C9)** — purge the experimental ESM/PLM selection family (preserve it in the private
   workspace first). Done early because it *shrinks* the surface every later PR has to sweep. The one
   PR that removes rather than renames — its own CHANGELOG "removed" entry.
3. **PR 3 (C2 + C3)** — finish `branch→lineage`; add the core two-term rate/modifier rule (retire
   "speed"); scope "odds" to the two undated-ALE tools; rewrite the stale `conventions.md` rows.
4. **PR 4 (C4 + C6)** — `--genome-model→--genome-resolution` (alias), `genomes` as the one domain
   word (coevolve node; `genes` warns), "level" reserved for the four domains. **C5 (plural commands
   `trait→traits`, `sequence→sequences`) is DEFERRED** — it uniquely changes an output file
   (`trait.log→traits.log`) and cascades into the params-log + a batch of test/doc updates, so it is
   left for a focused, human-reviewed follow-up rather than an autonomous rename. Likewise the
   sequence-clock `--branch-*` / `--family-speed` σ flags (C2/C3, above).
5. **PR 5 (C7)** — output filenames v2 (lowercase everywhere), behind `--legacy-filenames`, one
   CHANGELOG migration note; rewrites the `conventions.md` §Outputs schemas.
6. **PR 6 (C8)** — coevolve class/grammar alignment (`TraitGene*` stem, the edge→class table).
7. **0.4.0** — drop every alias and the `--legacy-filenames` shim.

Every PR except PR 2 (the deliberate removal) is byte-identical in simulation output (P4) and green
on the full suite; the only observable change is warnings on deprecated spellings.

## The rule, in one line

> **One word, one concept — `model` is a process, `level` is one of the four domains, `resolution`
> is the genome sub-axis, `rate` is a per-time quantity, `modifier` is a dimensionless multiplier,
> `lineage` is the tree-entity — and every rename keeps the old name working but *marked*, for
> exactly one version.**
