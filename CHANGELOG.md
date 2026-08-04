# Changelog

All notable changes to ZOMBI2 are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). While the project is pre-1.0, a **minor**
bump (0.x.0) carries new features or breaking changes and a **patch** bump (0.x.y) carries fixes.

A release is cut with `scripts/release.sh patch|minor|major` (the version is computed, not typed),
which moves the entries below from `[Unreleased]` into a dated version section.

## [Unreleased]

### Fixed
- **A run's output directory now describes that run and nothing else.** The per-unit directories —
  `gene_trees/`, `alignments/`, `phylograms/`, `genomes/`, `gff/`, `bed/` — hold one file per family,
  block or node, numbered by the run that made them, so a second run written into the same place
  interleaved two sets and left the leftovers indistinguishable from real output. A genome run that
  produced **zero** surviving families — announced clearly on the terminal — left the previous run's
  gene trees in place, and `zombi2 tools treedist` read one and printed `rf 0`, byte-identical to the
  earlier run's answer. They are emptied once, before a write fills them. `--flat` is untouched: there
  the directory is shared with every other output and every other level. (#316)
- **`--params` works on Python 3.10 again: `tomli` is now a declared dependency.** The package
  declares `requires-python = ">=3.10"` and CI tests 3.10, but the CLI's TOML reader falls back to
  `tomli` there and it was never in `dependencies` — supplied transitively by `pytest` and `mypy`,
  which both require it below 3.11. So the dev environment and CI had it and a plain
  `pip install zombi2` did not, and `--params` died with a raw `ModuleNotFoundError`. No test run
  inside this project could have found that; `tests/test_packaging.py` now compares the code's imports
  against the declared metadata instead. (#316)
- **`genome_summary.json` is written at the ordered and nucleotide resolutions too**, not only at
  family. It carries the corrected event counts — the ones that include a replacing transfer's
  displaced copy under `loss`, which the raw event log has no row for — and
  `docs/from-zombi1.md` names that undercount as the change most likely to hand a returning ZOMBI v1
  user a plausible wrong number, then points them at this file. It was missing at the two resolutions
  where the gap is *larger*: 64% at ordered, measured. All three resolutions now count events through
  one shared function, so they cannot drift. (#316)
- **`sequences_summary.json` no longer counts ancestral sequences it did not write.** They are
  reconstructed in memory either way but only land on disk when asked for, so a default run reported
  a count beside a directory that had none — and whoever inherits the folder cannot tell "never
  written" from "lost in transfer". `SequencesResult.summary()` is unchanged: it describes the run,
  where the written file describes the directory it sits in. (#316)
- **A run given no `--seed` now says which seed it drew**, on stderr, surviving `--quiet`. It was
  always recorded in the run report, but the report is a file the user has not opened: on screen an
  unseeded run looked exactly like a seeded one, and a class following a worksheet got 17, 19 and 15
  families where the sheet said 20 with no clue beyond the wrong answers. (#316)
- **The four gene-family rate defaults are stated in `--help`.** A pasted command that lost its tail
  still ran — at `--duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5`, none of which the
  help named, though `--initial-families` and `--max-family-size` named theirs. The wording carries
  the condition too: those apply only to a run given no rate at all. (#316)

### Added
- **A gallery example for `ByFamily`**: two genome runs on one species tree at the same mean rates,
  where only how much families differ from one another changes. With every family alike **no family
  is present in every genome**; a `ByFamily` draw gives 62 universal families and the bimodal
  gene-frequency spectrum real pangenomes show. `ByFamily` is the knob a comparative-genomics study
  leans on hardest and it had nothing showing what it does — a returning ZOMBI v1 user found it named
  nowhere in the migration guide and blank in the CLI help. Drawn with Phylustrator's new
  `genomes.grid`. (#316)

- **`zombi2 tools treedist --restrict` scores two trees on the taxa they share.** Differing leaf
  sets were an error, which refused the commonest comparison there is — a family's gene tree against
  the species tree, since only a universal single-copy family occupies every genome. A lecturer
  building a practical found not one of 22 families qualified, and an applied reviewer had none of
  his 21 single-copy families scored. `markers.tsv` already reported RF that way internally, so the
  capability existed and was simply not offered here. Opt-in, because silently scoring a different
  question than the one asked is worse than refusing; the refusal now names the flag. (#316)
- **The run report names the units of the rates it ran with** — what each is counted per, and that
  time is tree time in whatever unit the tree carries. A folder handed on is read long after the CLI
  help that answers this: `duplication 0.15` says nothing about per what. Only the slots a run used
  are listed, and the unit belongs to the slot rather than to the value — a modifier is a
  dimensionless multiplier by construction, so no rate expression can change it. Nothing
  dimensionless is listed, nor the continuous trait's `rate`, which is a variance rather than a rate.
  The `<level>.log` is untouched: it is `key<TAB>value` TSV that things parse. (#316)
- **The run report says where the software came from** — the project URL and a
  `pip install zombi2==<version>` line, pinned to the version that made the run. A folder handed on
  outlives the environment that produced it: a reviewer reconstructed an entire run from
  `run.zombi2` and still could not say where ZOMBI2 lived or what to install. (#316)
- **Images in the README resolve on PyPI.** They were repository-relative paths, and the README is
  the package's `long_description`, so the project page — a main discovery surface — showed six
  broken images. (#316)

### Changed
- **`tree.prune` takes `tips=` — a named set of leaves, whatever their fate**, beside the existing
  fate-based `keep="extant"`. The same operation on a different question: "the tree of the survivors"
  against "the tree of these taxa". `--restrict` prunes rather than intersecting clade sets, so
  branch lengths merge across the suppressed nodes and a length-aware metric still means something.
  (#316)
- **`zombi2[gallery]` asks for `phylustrator>=0.1.4`, not `==0.1.0`.** An exact pin on a companion
  library downgrades anyone who already has a newer one and makes any environment wanting one
  unresolvable. (#316)
- **A release is now gated on CI, and the README badge tells the truth.** 0.28.0 reached PyPI 43
  seconds after its CI run started and 19 minutes before it finished, so the artifact was published
  before a single test job had reported. `release.yml` now calls `ci.yml` as a reusable workflow and
  no publish job can start until the whole matrix is green. Separately, `cancel-in-progress` no longer
  cancels runs on `main`: a release push cancelled the run testing the commit before it, and GitHub
  renders cancelled as *failing*, so the badge told every visitor the build was broken while all six
  OS/Python jobs had passed. (#316)

## [0.28.0] - 2026-08-03

### Changed
- **BREAKING: each level now draws from its own random stream, so one seed on two levels no longer
  means the same numbers twice.** Every level opened `np.random.default_rng(seed)`, so two levels
  handed the same integer replayed the *same* PCG64 stream: a species tree's height and its genome's
  copy count at equal seeds came out correlated at Spearman −0.79 over 6000 seeds, with nothing in the
  output to say so. SPEC §2 calls two levels that do not read each other *independent*, and that is a
  claim about two random streams. Each level's generator is now spawned from a `SeedSequence` under a
  per-level key (`zombi2/rng.py`). **The species level is unchanged** — it keeps the unkeyed root
  stream, so every species tree ever grown under a seed still reproduces byte for byte — but
  **genomes, sequences, traits and joint runs at a given seed now produce a different (equally valid)
  realisation**. Pin the version for anything you need to rerun. (#315)
- **Branch lengths are written at full precision, so a tree round-trips through a file exactly.**
  `to_newick()` wrote 12 significant digits, which shifted every branch by ~2e−12 on the CLI's
  disk handoff between levels — enough to move every downstream waiting time, so
  `zombi2 genomes --seed 7` and `simulate_genomes_family(sp, seed=7)` produced *different* histories
  from the same tree and the same seed. Both were valid draws, but a seed that means one run through
  Python and another through the CLI is not a seed anyone can publish. `precision=` still takes a
  fixed digit count for a smaller file. (#315)
- **A genome run written with `stream_to=` now writes `species_complete.nwk` beside its outputs.**
  Every other file in that directory is indexed by the tree's node labels, and for a streamed run the
  directory is the only handoff there is. (#315)

### Added
- **`genomes.read_run(directory)` reopens a written genome run from Python**, and
  `simulate_sequences` accepts a directory or a `StreamedRun` wherever it accepts a result. The CLI
  has always reopened a run with `--from`; from Python the same handoff was a dead end, which mattered
  most for `stream_to=` — the feature whose whole point is that the run does not fit in memory, and
  whose handle could then not be passed to the next level. (#315)
- **A modifier of your own can declare the engines it was implemented for.** `Modifier` is public
  with a clean `factor(**context)` contract, and a subclass composed into a `Rate` correctly — and
  was then refused by every level, with no registry and no entry point, so extending the advertised
  grammar meant forking the package. Setting `implemented_for = ("species",)` on the subclass opens
  the gate for that engine and no other; everything undeclared is still refused by name. Seven
  engines take one; `sequences` does not and now says why, since it reads its modifiers itself rather
  than through the rate and could not honour one. Chapter 2 gains three worked examples — a rate
  following a measured curve, density dependence in the gene pool, and rearrangement scaling with the
  karyotype — each executed by the manual's doc tests. (#315)

### Fixed
- **`seed=None` from Python now records the seed it drew**, as the CLI already did.
  `docs/reproducibility.md` states there is no such thing as an unrepeatable ZOMBI2 run; that was true
  of the command line and false of the API, where `result.seed` stayed `None` and an interesting
  realisation found while exploring was gone. (#315)
- **A species run whose `species_fates.tsv` is missing is refused rather than silently consumed.**
  That is what SIGINT during the write phase leaves, and equally a partial copy or a tidied directory
  — and under `--sampling` every unsampled tip then read back as sampled, so the same seed and rates
  gave 30 extant genomes where the complete run gave 14, with exit 0 and nothing on stderr. An
  external tree still goes in as a file, which is what tells the two apart. (#315)
- **`--max-family-size` says so on stderr when it actually bound.** The cap is on by default, and when
  it binds it discards duplications and arriving transfers: at `--duplication 0.8` the default cap of
  10 has been measured pulling the realised rate to 0.32. It was recorded in `run.zombi2` and warned
  about in Chapter 4's prose, and nowhere a cluster job log would keep it — while the far less
  damaging all-genomes-empty case got a loud warning. (#315)
- **`FamilyGenomesResult.write(outputs=…)` and the ordered and nucleotide results now reject an
  unknown output token**, as species, sequences and traits always did. A typo wrote nothing and
  exited clean: silent data loss found three steps later, when the next tool has no input. (#315)
- **A Newick with duplicate tip labels is refused.** The simulation itself was fine, which is what
  made it dangerous: `names.tsv` — the documented join back to your own taxa — mapped two node ids to
  one name, and any downstream merge by taxon name silently duplicated or dropped rows. (#315)
- **A parallel run from an unguarded script now says what is wrong.** Workers re-import the caller's
  script; without `if __name__ == "__main__":` they ran the whole simulation again from the top, so
  the run "succeeded" having done N× the work and printed N copies of its output — or died with a
  bare `BrokenProcessPool`. Both paths now name the guard. A notebook and `python -c` are unaffected:
  they have no script to re-import, and still degrade to single-process. (#315)
- **A run reopened without its `genomes.tsv` says so** instead of raising `KeyError` from `.profiles`,
  and its `repr` no longer claims 0 nodes. (#315)
- **A discrete trait's `switch` rate honours any modifier it accepts, not only `DrivenBy`.** The
  engine chose between a constant generator and one rebuilt per stretch by asking "are there
  drivers?", which conflated two separate questions — *any* modifier makes the generator a function
  of the context, while only a `DrivenBy` names a source level to resolve. A non-driver modifier fell
  between the two and crashed the resolver. (#315)

### Changed
- **"wired" is gone from the vocabulary: `WIRED_MODIFIERS` is now `IMPLEMENTED_MODIFIERS`**, with
  `WIRED_SCOPES` and `WIRED_EXTENT_MODIFIERS` renamed to match. One word per concept — the levels
  declare what they *implement*, and the CLI's `RATES` help is still built from that declaration, so
  it cannot advertise what an engine refuses. These names are not in any `__all__`; nothing public
  changes. (#315)

## [0.27.0] - 2026-08-03

### Removed
- **Two optional-dependency extras that installed things for code that is not here.**
  `pip install zombi2[reconparser]` pulled in pandas for `zombi2/tools/reconparser/`, which does not
  exist, and `zombi2[bench]` pinned snakemake, its SLURM executor plugin, matplotlib and pyyaml for a
  Snakefile that is not in the tree — both advertised on PyPI. An extra is a promise about what is
  installable; it goes back in when the code does.

### Changed
- **The documentation says what the code does, in twenty-odd places it had drifted from.** The
  load-bearing one: Chapter 3 said `species_extant.nwk` was "what the next level reads", and it is
  the **complete** tree that the other levels run along — which is the whole point of keeping the
  dead lineages. Also: rearrangements are their own file at both structured resolutions (not rows in
  `genome_events.tsv`); `chromosome_events.tsv` has no lineage column; `assembly()` returns 3-tuples,
  not 5; every family some node still carries gets a gene tree, not only those surviving in an extant
  leaf; `Event`/`GeneEdge` and `edges_from_tsv` under their 0.26.0 names; the t=0 driver row is
  spelled `initial`, not `root`, so a hand-built file following the manual was rejected;
  `--gtr-rates` is `--exchangeabilities` (and the test that "covered" it was passing because the flag
  no longer parsed); `at_speciation` is refused for threshold traits; a `k×k` switch matrix cannot
  carry a modifier; `--at-speciation` takes a **variance**, not a width; `Tree`/`Node`/`prune` are
  documented under `zombi2.tree`; and the README no longer says ancestral sequences are written by
  default, since they are an opt-in `--write ancestral`.
- **Appendix A's modifier table gained the two rows it was missing**, and the test that guards it now
  covers them. The single "Traits" row described the *continuous* rate only, while a discrete
  `switch` rate takes `DrivenBy` and nothing else; the joint level had no row at all. The discrete
  engine now declares its `WIRED_MODIFIERS` like every other level, and `zombi2 joint -h` builds its
  RATES block from the joint engine's declaration instead of a hand-written tuple — it had been
  under-reporting `OnTime` and `OnTotalDiversity`, which a joint run does thread.

### Fixed
- **`--fasta` was accepted and ignored on a family or ordered run.** It was missing from the
  nucleotide-only list its own sibling `--gff` is in, so `zombi2 genomes DIR --fasta x.fasta` ran an
  ordinary family simulation, never opened the file, and wrote `fasta  x.fasta` into `genomes.log` as
  though it had applied — a run recorded as something it was not. It is refused now, like every other
  knob a resolution does not have.
- **The `Between` kernel guard reaches the last three slots.** The genome engines refuse a kernel in
  a rate; the two trait engines and the joint engine did not check, so a kernel there died mid-run
  with `AttributeError: 'Between' object has no attribute 'multiplier'` instead of naming the
  mistake. All three now go through the same `check_not_a_kernel`.
- **`correlation=` alongside `regimes=` was silently dropped.** `regimes=` dispatches before the
  correlated engine and threads no correlation, so the run was quietly the uncorrelated model: two
  runs differing only in a bogus correlation came back byte-identical. It raises now, and says why —
  multi-optimum OU evolves one trait, so there is no second trait for a correlation to be with.
- **`TraitsResult.summary()` on a correlated run.** A correlated run holds one value per trait at
  every node and `summary()` assumed one number: the continuous case raised `TypeError: float()
  argument must be … not 'dict'`, so `write(..., "summary")` failed outright for multi-trait
  continuous data, and the threshold case counted whole dicts as states, producing keys like
  `"{'a': 'b', 'b': 'b'}"`. Both now report per trait, in the shape `trait_values.tsv` already writes.
- **The ordered resolution is about 1.5× faster, with byte-identical output.** `_pick_gene` measured
  each lineage's genome with `_genome_size()` to decide whether the draw landed in it and then walked
  the chosen genome again through `_gene_in()` — every chromosome of every skipped lineage read once
  and the chosen lineage's read twice, per event, inside the Gillespie loop. Counting chromosome by
  chromosome finds the same gene in one pass. The draw is untouched, so a run is the same run: 1.70 s
  → 1.11 s on a 1,000-tip tree (medians of three interleaved rounds), and the event logs and genomes
  hash identically at three seeds. This recovers a regression that arrived in 0.21.0.
- **A `Between` kernel in a rate or extent slot now says so, everywhere.** A kernel weights a
  recipient by the (donor, recipient) pair, so it answers *who receives* and belongs in
  `transfer_to`; a rate is read on one lineage and has no donor. The family and ordered engines
  refused it on a rate, but the ordered engine did not check its **extents** and the nucleotide
  engine checked neither. Those two did not do the wrong thing quietly — `Between` implements no
  `multiplier`, so a run died part-way through with `AttributeError: 'Between' object has no
  attribute 'multiplier'`, a traceback from inside the engine naming neither the rate nor the
  mistake. All five slots now refuse before the run starts, through one shared guard so the message
  is the same wherever the kernel was put.
- **A joint gene-content mapping that can never fire is refused, as a joint trait's already was.**
  A typo'd `DrivenBy("genomes:toxin", {"presnt": 3.0, "absent": 1.0})` left every lineage at the
  default factor — birth was never driven, the run was the plain birth–death model — and it completed
  in silence, reporting a coupled run. Both gene-content drivers have an alphabet known before the
  race starts (a named family is `present` or `absent`; a count is a number, which a `{state: factor}`
  table can never equal), so the check is exhaustive.

- **Conditioning on a continuous trait from disk drove every lineage at one constant value.** A
  diffusion has no switches, so its `trait_events.tsv` holds only the `initial` row — and replaying
  that log built a driver frozen at the root value on every lineage, *accepted without a warning*. A
  run that looked conditioned was the undriven model with one constant factor, and the difference was
  large: on a 30-tip tree the same driven-loss run gave 79 losses from the file against 157 from the
  in-memory trait. The event log now **refuses** and names the file to use, and `trait_values.tsv` —
  which has always carried every node's value — is a driver in its own right, reproducing the
  in-memory result event for event. Discrete traits are untouched: their event log is still the exact
  stochastic map, and a genuinely never-switching discrete trait still loads as the constant it is.

### Changed
- **A continuous driver's resolution is cut per unit of time, and it is tunable.** It used to be a
  fixed eight stretches per branch, which makes the approximation as coarse as the branch is long:
  the error was worst exactly where the driver had had most time to move, and it could not be
  refined. `DrivenBy(source, mapping, step=…)` now takes a **duration** in the tree's own time units,
  so every stretch means the same thing wherever it sits — a branch twice as long gets twice as many,
  and halving the step doubles them everywhere. `step=None` takes 1% of the tree's height, a fraction
  rather than an absolute because a tree may be measured in substitutions or in millions of years.
  The same trait read at two resolutions stays two drivers rather than being silently resolved once
  and shared.

  Worth knowing when reading a driven run: the interpolation is the **straight line** between a
  branch's endpoint values, which is the mean of the Brownian bridge with the excursions dropped.
  Under a non-linear response curve those do not average out, so a smaller `step` is not only more
  precise, it removes a bias.


## [0.26.0] - 2026-08-02

### Changed
- **The conditioning diagram can draw the target's own Markov chain, and mark which arrow is
  driven.** A rate that reads a driver is usually **one transition**, not the whole chain, and a
  diagram that does not say so claims the model is symmetric when it is not. The `gene_drives_trait`
  figure drives `harmless → pathogenic` only — a toxin makes a lineage dangerous, it does not help it
  recover — and the driven arrow is drawn heavier than the one that is not.

### Added
- **Gene families can be grouped into modules, and a module's completion drives a rate.**
  `modules={"flagellum": ["flgA", …]}` names a group of declared families, and
  `g.completion("flagellum")` is a conditioning driver giving the fraction of it a lineage carries —
  a number in `[0, 1]`, read with a `Curve` like any continuous driver. A **fraction rather than a
  yes/no on purpose**: under independent loss the chance every family of a module survives falls off
  geometrically with its size (measured on a 200-tip tree, a module of three was complete at 189 tips
  and one of six at none), so a complete/incomplete driver would be a constant for anything but the
  smallest modules. A threshold is expressible where every other response shape already lives, in the
  curve. A module of one family is exactly that family's presence. Members must be named with
  `family_names=`; an anonymous family's id comes from the order events fired in.

### Fixed
- **Two gallery figures had silently lost their transfers.** They filtered a genome run's `.events`
  on `kind == "transfer"`, which after that attribute came to mean one row per event is never true —
  the file's vocabulary is `transfer_additive` — so `and e.recipient is not None` short-circuited
  before touching a field that no longer exists, and the list came back empty rather than raising.
  170 transfers read as 0 and the transfer-highway chart drew nothing, with no error anywhere. They
  read `.edges`, which is where a per-branch donor and recipient live.

### Changed
- **`result.events` now means what a row of `genome_events.tsv` means.** It used to hold one entry
  per gene-tree *edge*, so a duplication was two of them and a transfer likewise — counting
  duplications in Python gave twice the file's number, and a filter on `kind == "transfer"` matched
  everything in Python and nothing in the file, which is the mismatch that forced Chapter 4's
  wording to be corrected once already. The two now agree, kind for kind and count for count.
  **Breaking** for code reading the old shape: the edge list is `result.edges`, unchanged, and it is
  what a gene tree is built from. Family and ordered runs — traits and the nucleotide level already
  matched their files, which is how the gap was found.
- **The classes say which is which.** `Event` is one genome event, with the copies it ended in
  `parents` and the copies it began in `children`; the per-edge record is `GeneEdge`.
  `gene_trees_from_events` is `gene_trees_from_edges` and `events_from_tsv` is `edges_from_tsv`,
  because both take edges. `zombi2.species.Event` is untouched — a species-tree event is a different
  thing and keeps its name.
- **No run changed.** The reproducibility digests are byte-identical across this, because they hash
  `edges` — the finer record, and the one they always hashed.

## [0.25.0] - 2026-08-02

### Added
- **A gene family's presence can drive a rate.** `g.presence("tox")` is a conditioning driver like a
  grown trait, so `switch=0.1 * DrivenBy(g.presence("tox"), {"present": 8.0, "absent": 1.0})` makes a
  trait switch faster in lineages that carry the gene. Until now driving only ran one way — a trait
  could make gene loss faster, but a gene could not make a trait's rate faster — and the reason was a
  missing translator rather than a missing mechanism: the thing that answers *what state was this
  lineage in at time t* only knew how to be built from a trait. It reads the family's gene tree, so
  the signal changes **mid-branch** where a copy was actually gained or lost rather than only at the
  nodes, and a lineage that never held the family answers `absent` rather than raising. Family and
  ordered runs; only families named with `family_names=`. `presence(...).history(tree)` gives the
  per-branch map in the same shape `TraitsResult.history` has, so anything that draws a trait's
  history down a tree draws a gene's too — which is what the new gallery entry does.
- **A gallery entry for a module driving a trait through a step** — `module_drives_metabolism`.
  Four families make up aerobic respiration, and the response is *discontinuous*:
  `lambda f: 20.0 if f > 0.5 else 1.0`, so more than half the module makes a lineage aerobic and less
  makes it revert. Both directions read the module, oppositely, so the trait tracks gene content
  rather than accumulating — 97% of the tree's branch length has the trait on the side of the
  threshold its completion is. The curve diagram draws a step as a step rather than a ramp, which is
  the one thing a threshold is not.
- **A gallery entry for a gene driving a trait** — `gene_drives_trait`, the same tree painted twice:
  by the toxin family's presence, then by the pathogenicity whose switch rate reads it. The rates are
  chosen so the two panels have something to disagree about: the family covers 62% of the tree's
  branch length, and the realised switch rate is 1.1 per unit where it is present against 0.07 where
  it is not.

### Added
- **Site-specific amino-acid profiles.** `simulate_sequences(..., profiles={family: array})` gives a
  family one set of equilibrium frequencies **per position** instead of one shared by the whole gene,
  which is the difference between a buried hydrophobic site and the loop beside it. Each row of the
  `(L, K)` array becomes that site's own model over the base model's exchangeabilities — which pairs
  of residues interchange easily is chemistry and is kept; where each residue belongs is the
  profile's to say. Where the numbers come from is open: the manual shows one recipe from an
  alignment and one from a protein language model, whose output already *is* a distribution over
  amino acids at every position. Families without a profile are untouched, a flat profile is the
  model it was built from, and profiles compose with `+Γ` — a profile says which residues, a Gamma
  says how fast. Refused alongside `partitions` (both decide a family's per-site models), alongside
  `parallel`, and for a family the run does not have. An amino-acid profile needs a protein model and
  so belongs to a family or ordered run — a nucleotide genome is measured in base pairs and refuses
  protein models — where a profile is over the four bases instead, one row per base pair.

### Added
- **The sequence and trait levels are now validated against theory, not only against themselves.**
  `tests/test_validation.py` checks a run against a closed form rather than against an invariant, and
  it covered the species and genome levels only — the two levels a user is most likely to hand to a
  method for grading had no such check at all. Twelve more: JC69's `3/4·(1 − e^(−4d/3))` and Kimura's
  separate transition and transversion probabilities (written out from the 1969 and 1980 papers, not
  taken from `p_matrix()`); HKY85 and LG holding the frequencies they were given and keeping `i→j`
  and `j→i` balanced, which is reversibility; `+I+Γ` giving the Jukes–Cantor curve *averaged over* the
  rate classes rather than evaluated at their mean, so classes that are computed and then never
  applied are ruled out; both lineage clocks drawing mean-1 factors, where the historical lognormal
  bug lived; Brownian motion's per-branch `Normal(0, σ²·Δt)` and Felsenstein's tip covariance;
  Ornstein–Uhlenbeck's exact transition moments; an Mk trait firing each direction at the rate that
  direction was given (which catches a transposed rate dict) and saturating at ½ the way the chain
  says; correlated traits realising their ρ; and a trait driving the substitution rate, which is
  exact rather than statistical — a phylogram branch is the driven rate *integrated* over the branch,
  so it can be checked against the trait's own segments to floating point.

  Each test also rules out the plausible wrong model, not just confirms the right one: the
  uncorrected p-distance, κ = 1, a single rate across sites, an uncorrected lognormal clock, tips that
  ignore their shared ancestry, a diffusion with no pull, a chain whose switches never reverse,
  uncorrelated traits, and sampling a driver once per branch instead of integrating it. Resolution is
  measured by mutating the model rather than asserted: the checks catch a substitution rate off by
  2%, κ by 3%, a Brownian variance-rate by 5%, an OU pull by 2%, an Mk rate by 3%, a declared
  correlation by 5% and one stationary frequency shifted by 0.005 at fifteen standard errors. They
  add about four seconds to the suite.
- **Site profiles are validated too.** Over an equal-exchangeability model a profiled site is
  Felsenstein's F81, so `S·(1 − e^(−d/S))` with `S = 1 − Σπ²` gives each row its own closed form. A
  blocked profile — a flat half and a sharply peaked one — is checked half by half, which is what
  catches a profile applied to the wrong sites: the *total* over sites cannot, because summing over
  positions is permutation-invariant, so a shuffled profile has exactly the same total divergence as
  the right one. The formula also pins the renormalisation each rebuilt per-site model goes through,
  without which a peaked profile would quietly run at its own rate and a phylogram would stop meaning
  substitutions per site. A second check covers the other half of what a profile promises — that it
  says *where* residues belong, not *which pairs interchange* — by confirming changes still follow
  `π_i·π_j·S_ij` with the base model's exchangeabilities. That one is the loosest check in the file
  (one exchangeability must be off by about 20% to show, since pinning six of them needs branches
  short enough that a difference is still a substitution), but it rejects the failure the code path
  invites — rebuilding a site from its frequencies and dropping the chemistry, collapsing every model
  to F81 — at 270 standard errors.
- **The docs open on a worked study.** `analyses/red/` — does RED, the measure GTDB uses to normalise
  taxonomic ranks, still recover node ages once molecular rates vary? — was in the repository and
  linked from nowhere: not the README, not the site, not the manual. It is now a documentation page
  sitting directly after Home, framed as an example of the kind of question the tool exists to
  answer, and its figures are written into `docs/assets/red/` by the study's own `figures.py` so the
  two cannot drift apart.
- **The RED study now sweeps the autocorrelated clock too**, which closes the limitation it was
  written around. `FromParent` was a species-level modifier when the study was written and is wired
  at the sequence level now, so the claim that the uncorrelated case was all that could be tested had
  gone stale. At the raggedness real archaea show, RED recovers relative node ages with r = 0.993
  (nRMSE 2.3%) under autocorrelated rates against r = 0.94–0.95 (≈6%) uncorrelated — confirming the
  study's own prediction that autocorrelation is the easier case, and turning a caveat into a bound:
  the conclusion no longer depends on which arrangement real archaea have. The uncorrelated numbers
  are unchanged to four decimals.

### Removed
- **The yeast synteny/inversion study has left the repository** (`analyses/synteny_inversions/`). Its
  headline rate was stated in the analyses index where anyone browsing would read it, and it is not a
  number we are ready to stand behind — an unlinked folder is not the same as an unpublished result.
  The code, the data and the write-up are kept outside the repository and will come back when the
  estimate does. With one study left, the index is gone too: its conventions now sit in that study's
  own `REPORT.md`, under "Reproducing this recipe".

## [0.24.0] - 2026-08-02

### Added
- **A stated reproducibility contract, and a test that checks it.** What a seed guarantees was
  true but nowhere written down, and nothing verified it: every other test asserts a *property* of
  a run, so the whole suite would still pass if the random stream shifted by one draw and every
  number came out different. `tests/test_reproducibility.py` now hashes a run of each level against
  a recorded digest, and CI runs it on Linux, macOS and Windows across Python 3.10–3.13 — six
  configurations, six independently randomised hash seeds, one number. The contract itself is a new
  documentation page, including what is *not* promised (across versions, and between a serial and a
  `--parallel` run).
- **`mypy` is part of CI.** ZOMBI2 ships `py.typed`, so its annotations are read by the type checker
  of everyone who installs it — and 118 errors said they were not worth reading. All of them are
  fixed and the check is a CI job. Nothing about a run changed; the fixes are annotations, three
  latent `None` paths that could not be reached, and two stale signatures.
- **Coverage is measured in CI** and reported (95% of `zombi2` at the time of writing). Not
  enforced: a threshold turns a number worth reading into a number worth gaming.
- **A `CONTRIBUTING.md`** — how to set up, what to run before opening a pull request, and what a
  pull request is expected to carry.
- **Three gallery examples for continuous drivers.** A `Curve` maps a driver's value to a factor,
  and one exponential example did not show what that buys: a **saturating** curve, whose factor is
  bounded however high the driver goes; a **humped** one, fastest at an intermediate value, which no
  table of per-state multipliers can express; and **one trait driving another**, where the rate a
  body-size trait diffuses at reads a temperature trait, drawn as the same tree painted twice.

### Changed
- **The documentation is ordered the way it is read**: Home · Guide · Output files · Tools · API
  reference · Performance · Coming from ZOMBI v1 · FAQ.
- **`Chromosome.genes` is typed `Sequence[Gene]`, not `list[Gene]`.** A chromosome has two states —
  the engine's live one, which it mutates, and the frozen snapshot a finished run hands back, whose
  genes are a tuple. The old annotation described only the first and so was wrong about every
  chromosome a user ever touches.

## [0.23.0] - 2026-08-02

### Added
- **The README opens on one simulated dataset seen at all four levels** — the same 30-tip tree drawn
  with its extinct lineages, the gene order of every surviving genome with homologues linked, the
  alignment behind one gene family, and two coupled traits. The figure regenerates with
  `figures/scripts/fig_overview.py`.
- **A gallery entry for gene order along a tree** — `genome_synteny_tree`, ordered genomes drawn
  beside the tree that produced them with homologous genes linked across the tips.

### Changed
- **The CLI's flags and the Python keywords are the same words.** Three had drifted apart:
  `--frequencies` was `freqs`, `--gtr-rates` was `gtr(rates=)` and `zombi2 joint --trait-start` was
  `discrete(start=)`. **Breaking:** `hky85(freqs=…)` and `gtr(freqs=…)` are now `frequencies=`,
  `gtr(rates=…)` is `exchangeabilities=` (the name `reversible()` already used for the same six
  numbers), and the joint flag is `--start`, matching `zombi2 traits`. `--mass-extinction` keeps its
  singular name against a plural `mass_extinctions=`: the flag is repeatable and names one pulse
  where the parameter holds the list.

### Fixed
- **The gallery's conditioning examples read the trait by tip name.** Four of them still indexed
  `TraitsResult.values` by node id, which stopped working when that mapping was re-keyed — so the
  code printed beside those figures raised a `KeyError` for anyone who copied it.
- **`simulate_sequences` accepts an ordered genome run.** Its gate tested the result's class rather
  than what the level needs — `gene_trees` and a `complete_tree`, both of which an ordered run has —
  so `simulate_genomes_ordered(...)` into `simulate_sequences(...)` raised in Python while the
  identical two commands worked on the command line, whose directory handoff quietly rebuilds a
  family result on the way in. The docstring, Chapter 7 and Appendix B all said the Python route
  worked. A test now asserts that every CLI flag has a Python route, so the two front doors cannot
  drift apart again.
- **The extant tree keeps the complete tree's clade order.** `prune` rebuilt each node's children in
  ascending node id, which coincides with the original order only until a child is pruned away — the
  surviving descendant that replaces it can have a far larger id than its sibling, so the pair came
  out swapped. A figure showing both trees drew the same clade on opposite sides, and anything
  joining them by position disagreed with itself.

## [0.22.0] - 2026-08-01

### Fixed
- **Every tree ZOMBI2 writes is now ultrametric when read back.** Branch lengths were written at 7
  significant digits, and a tip's depth is a *sum* of them, so the rounding accumulated down the path
  and left `species_extant.nwk` about 1e-6 off — far above what `ape::is.ultrametric()` allows, on a
  tree that was ultrametric to 1e-16 in memory. The first thing anyone does with it in R therefore
  failed. Trees are written at 12 digits; `Tree.to_newick()` takes `precision=` for anything else.
- **Chapter 4 said `transfer`, the log writes `transfer_additive` / `transfer_replacing`** — a reader
  filtering on the documented kind got nothing back. The chapter now says what the log says.
- **The `"rest"` group is documented.** `Between({("rest", "A"): 8.0})` makes a named clade a transfer
  hotspot the whole rest of the tree donates into; it worked and appeared nowhere in the manual.
- **Chapter 9 no longer promises the staleness guard across run directories.** It records the
  dependency in the run directory, so a driver and a target written to two different directories with
  `--from` are not linked — which the chapter stated without qualification.
- **A trait dataset now joins the tree it came from.** `TraitsResult.values` was keyed by bare node
  ids (`5`) while every Newick label and every `trait_values.tsv` row says `n5` — so in Python the
  comparative vector and the tree beside it shared **no keys at all**, and nothing said so. It is now
  keyed by the tip's name. **Breaking** for code that indexed it by id: `values_by_id` is the old
  view, unchanged. The written files always agreed with each other; it was only the in-memory pair
  that did not.
- **`zombi2 tools tree --round` now produces a file that is actually ultrametric.** The snap was
  exact in memory (tip depths agreed to ~1e-16) and the writer undid it: a depth is a *sum* of branch
  lengths, so writing them at the usual 7 significant digits reintroduced a spread of ~1e-6 — well
  above what `ape::is.ultrametric()` allows, so the tree still came back rejected. `--round` writes
  at full precision; `Tree.to_newick()` takes a `precision=` for anything else that needs it.
- **`gene_order.tsv` and `initial_genome.tsv` record each chromosome's `topology`.** It decides where
  a segmental event stops and which chromosomes may fuse, and it is what a rearrangement format's
  per-chromosome terminator depends on — but it appeared in no output file, so a mixed circular and
  linear karyotype left nothing on disk saying which chromosome was a ring.
- **The `TO REPRODUCE` block runs.** It listed commands in pipeline order, which for a conditioned run
  is the wrong order: a rate driven by a trait must run after the trait that writes the file it
  reads, and traits come last in the pipeline. Copy-pasted, it failed on that line. A driver is now
  promoted above whatever is conditioned on it; a run with no conditioning is unchanged.
- **A `DrivenBy` mapping that names a state the driver never takes now says so.** One typo among
  otherwise-correct keys used to pass in total silence: the guard refused a mapping where *nothing*
  matched, but a mapping with one good key and one typo fires, so the run completed, reported itself
  as driven, and applied the factor the user cared about to nobody. It warns rather than raises,
  because when a driver is replayed from a file the only states known are the ones it actually
  reached, and a mapping may legitimately name one this realisation missed. (A joint run, whose
  alphabet is declared up front, already raised.)
- **The command line prints a library warning in its own voice** — `zombi2: warning: …` — instead of
  Python's file-and-line rendering, which for a CLI user wraps the one sentence that matters in noise
  and reads like a crash.
- **The conditioning and joining diagrams read backwards.** Their arrow ran from cause to effect —
  habitat to loss — but was labelled `DrivenBy`, which is passive, so reading along the arrow gave
  "habitat is driven by loss": the opposite of the model. The joining figure contained the proof, one
  arrow reading forwards (`creates`) and the other back. Both arrows now carry an active verb, and
  `DrivenBy` sits under the **target**, where it reads correctly and where it is actually typed. The
  gallery's five conditioning cards are drawn from the same helper and change with them; on the one
  card whose target is a choice slot (`transfer_to`) the expression now correctly shows no base.

## [0.21.0] - 2026-08-01

### Added
- **Extents on the command line at the `ordered` resolution** — `--inversion-extent`,
  `--duplication-extent`, `--loss-extent`, `--transfer-extent`, `--transposition-extent`,
  `--translocation-extent`, each the mean number of **genes** an event takes. They were reachable
  from Python only, so a command-line ordered run could make nothing but single-gene inversions —
  which flip one gene's strand and shuffle nothing, so the CLI produced something that looked like a
  rearrangement dataset and was not one. The same flags mean base pairs at `nucleotide`, as before.
- **`--topology` takes one label per chromosome** — `--topology circular,linear` for a mixed
  karyotype, which the manual has always shown for the Python argument and the flag could not express.
- **A correlated trait run carries the event log every other continuous run carries** — the `initial`
  row and one `on_speciation` row per jump, where it previously returned nothing and wrote a
  header-only `trait_events.tsv`. The table **widens** rather than repeating a row per trait
  (`from:<trait>` · `to:<trait>`, one pair apiece, as `trait_values.tsv` already does), because a
  correlated jump moves every trait at once and is one event.
- **`zombi2 traits` records `conditioned_on`** when its rate was driven, like `genomes` and
  `sequences` do, so a trait driven by another level is no longer an untracked dependency.
- **A substitution model can be built from your own matrix.**
  `substitution_models.reversible(exchangeabilities, freqs, name=…, alphabet=…)` takes a symmetric
  exchangeability matrix and stationary frequencies over any alphabet and normalises them like every
  model on the menu, so `gtr()` is visibly its four-state special case and `lg()` its twenty-state
  one. Python API: a K×K matrix is not a command-line-shaped thing, so there is no flag.
- **A family's sites can be split into partitions, each under its own model** —
  `simulate_sequences(genomes, partitions=[(hky85(kappa=2.0), 600), (jc69(), 400)])` in place of
  `model=` and `length=`. Every partition shares the run's one alphabet and one substitution rate,
  and every model is normalised the same way, so the family keeps one phylogram that is exact for all
  of them; each partition may carry its own `across_sites()` classes. Family and ordered runs only.
  Composes with `parallel=` and `stream_to=`. Experimental, Python-first.
- **Transfer steering at every resolution.** `transfer_to` now takes `Clades({…})` — weight by named
  clade — and `mod.DrivenBy(source, mapping)` — a trait-driven recipient weight, with a `Between`
  kernel that reads the donor too — at the **ordered** and **nucleotide** resolutions, not only at
  `family`. Same choice slot, same kernel: the numbers redistribute transfers without changing how
  many happen, a weight of 0 still means "cannot receive", and a transfer whose every candidate
  weighs 0 does not fire. `--transfer-to` already accepted the written form everywhere, so the flag
  works unchanged.
- **A trait can drive the substitution rate.** `DrivenBy` now works on `substitution`, so a lineage's
  habitat or lifestyle sets how fast the sequences inside its genes evolve — the same modifier that
  drives a genome rate, composing with either lineage clock. A driver that switches **mid-branch** is
  integrated across the switch rather than sampled once for the branch, so the phylograms and the
  clock species tree are the trees the alignments were really drawn along. Reachable from the command
  line in the rate's written form, as everywhere else:
  `--substitution "0.05 * DrivenBy('out/traits/trait_events.tsv', {'cave': 0.5, 'surface': 1.0})"`.
  `zombi2 sequences` writes a `conditioned_on` marker when its rate was driven, so re-running the
  trait beneath it refuses rather than leaving the sequences silently stale.
- **A trait can drive an ordered genome run.** `DrivenBy` now works at the `ordered` resolution — on
  the gene-family four, on `inversion` / `transposition` / `translocation`, on the chromosome tier
  and on every extent — so a lineage can rearrange its gene order faster, or in longer runs of genes,
  where a trait says so. Read wherever it switches mid-branch, as at the other two resolutions. The
  ordered engine now declares its own `WIRED_MODIFIERS` and `WIRED_EXTENT_MODIFIERS`, and
  `zombi2 genomes -h` builds its per-resolution sentence from them instead of a hand-written list.
- **Across-site rate variation — `+Γ` and invariant sites.**
  `hky85(kappa=2.0).across_sites(gamma_shape=0.5, invariant=0.1)`, or `--gamma-shape` ·
  `--invariant` · `--rate-categories` on the command line, gives every site one of a
  discretised-Gamma set of rate classes (four by default) plus an optional class that never changes.
  It decorates the **model**, the way the field spells it — the run reports `HKY85+I+G4` — and not
  the rate: across-site variation is not a modifier, and `substitution` still takes a lineage clock
  and nothing else. The classes are normalised to mean one, so a phylogram's branch lengths are
  still substitutions per site, now the mean over them, and a run with variation and one without are
  directly comparable. Works on every model on the menu, nucleotide and protein alike.
- **Ornstein–Uhlenbeck takes a modified variance-rate.** `reverts_to` / `pull` now compose with
  `OnTime`, `FromParent`, `OnTotalDiversity` and `DrivenBy` on `rate`, so a trait that bursts early
  *and* reverts to an optimum is one rate with one modifier and two arguments — which is what
  Chapter 8 has said all along and the code refused. The per-branch variance is the exact weighted
  integral `∫ e^{−2α(t₁−s)} σ²(s) ds`, stepping where the schedule, the standing diversity or the
  driver steps; the Brownian `∫ σ²(s) ds` is a different number, larger by an order of magnitude on
  a typical branch.
- **Multivariate Ornstein–Uhlenbeck** — `reverts_to` and `pull` alongside `correlation=`, each
  taking one value shared across the traits or one per trait. Each trait reverts to its own optimum
  at its own strength and the correlation rides in the diffusion. The drift is **diagonal**: one
  trait's deviation does not pull another, and a full drift matrix is refused by name rather than
  quietly approximated.
- **Jumps at speciation combine with correlated traits and with regimes.** A correlated jump is
  drawn under the same `correlation=` overlay the diffusion uses, and a multi-optimum OU (`regimes=`)
  run now jumps at each split like any other continuous trait.
- **A genome run says when it emptied a genome.** There is no floor at the family resolution — loss
  is counted per copy, and the last copy is a copy like any other — so a high loss rate can strip a
  lineage of every gene. The run now reports how many extant genomes came out empty, on stderr and
  as `empty_genomes` in `genome_summary.json`, instead of leaving a reader to work it out from
  `profiles.tsv` having no rows.
- **Relaxed (per-lineage) diversification rates** — `birth = 1.0 * mod.ByLineage(spread=σ)`. The
  species level already took the *inherited* form of rate variation (`FromParent`, ClaDS) and
  refused the independent one, which left the model with no null: `ByLineage` spreads lineage rates
  the same way and inherits none of it, so the tree-shape signature heritability leaves — fast
  clades hoarding the tips — can be told apart from the rate variation itself. It works on `birth`
  and on `death`, independently of each other. A rate carrying both `FromParent` and `ByLineage` is
  refused rather than letting one of them silently win: they are two answers to the same question.

### Fixed
- **`chromosome_loss` never takes an ordered genome's last genes.** It refused only the genome's last
  *chromosome*, so a lineage holding one gene-bearing chromosome beside an empty replicon — one
  `chromosome_origination` minted, or one a translocation emptied — could lose everything it had.
  The same floor `_lose_at` enforces one tier down, for the same reason.
- **`simulate_joint` declares what it supports and rejects the rest**, as every other level does. Its
  gate named two modifiers and let everything else through, so `ByFamily` on a joint `birth` was
  accepted and then silently returned a factor of 1. (`OnTime` and `OnTotalDiversity` were never
  affected — the loop threads both and steps at their breakpoints.)
- **A rate matrix that is not time-reversible is refused instead of being evolved under wrong
  transition probabilities.** `SubstitutionModel` is public, and `p_matrix()` computes `exp(Qt)` by
  eigendecomposing the symmetric `diag(√π)·Q·diag(1/√π)`, which is similar to `Q` only under detailed
  balance — so a hand-built non-reversible `Q` was quietly replaced by a different matrix and
  produced plausible, wrong sequences with no error at all. It now raises, naming the violated
  identity and pointing at `reversible()`. Also refused: rows that do not sum to 0, non-positive or
  unnormalised frequencies, an alphabet that does not name every state exactly once, and an
  exchangeability matrix with a non-zero diagonal (which was silently discarded).
- **A transfer the ordered engine drops because no candidate can receive now leaves the donor
  untouched.** The recipient is chosen before the donor's chromosome is anchored, so a run that wraps
  position 0 no longer rotates the donor's gene list for an event that did not happen. Byte-identical
  for every existing run: the pick consumes the random stream and anchoring does not.
- **`family_speed` beside a driven rate is refused instead of running a mismatched model.** It is a
  `ByFamily` draw and was missing from the guard that refuses `ByFamily` beside `DrivenBy`, so the
  run was accepted and then summed the total *without* the per-family multipliers while drawing the
  copy *with* them — a total saying one thing and a pick doing another.
- **`zombi2 genomes` records `conditioned_on` and the driver's SHA-256 for a `DrivenBy` on any rate**,
  not only `--duplication` / `--transfer` / `--loss` / `--origination` / `--transfer-to`. A run driven
  through `--inversion` or `--fission` — already legal at the nucleotide resolution — left no
  conditioning marker, so re-running the trait beneath it did not know it had orphaned the genome run,
  and pinned no digest of the driver file in the log.
- **A circular chromosome no longer fuses with a linear one.** At the ordered resolution `fusion`
  drew its partner from every other chromosome and gave the fused child whichever topology was
  picked first, so a ring and a molecule with two ends silently became one molecule — on a
  two-topology karyotype it collapsed almost every genome. The partner is drawn from the
  same-topology chromosomes only, the rule the nucleotide resolution already enforced; a genome of a
  single topology draws exactly as it did before.
- **`zombi2 sequences` no longer crashes on a run with no gene families.** A run whose genomes
  emptied, or one started with `--initial-families 0`, has no alignments, so mean pairwise identity
  is undefined; the run report tried to render it as a percentage and raised.
- **A joint run no longer accepts a per-lineage rate it does not thread.** `simulate_joint` rejected
  `FromParent` and checked `DrivenBy`, but let every other modifier through and ignored it, so
  `birth = 1.0 * ByLineage(...)` returned a tree grown without the rate variation asked for. It now
  raises, which matters more since the same expression started working on `zombi2 species`.

### Changed
- **The sequences level refuses a `divergence` given alongside a driven `substitution`.** The base is
  solved for by assuming the modifiers average to 1 along a root-to-tip path, which the two lineage
  clocks are mean-corrected to do and a driver deliberately is not — so allowing it would log a
  divergence the run does not realise.
- **A `regimes=` trait run writes the event log every other continuous run writes** — the `initial`
  row at t=0 and one `on_speciation` row per jump. It previously returned an empty log and wrote a
  header-only `trait_events.tsv`.
- **Clearer refusals at the traits level.** The combinations that stay blocked — a modified
  variance-rate with `regimes=`, per-trait modifiers with `correlation=`, a per-regime jump size, a
  modified liability variance-rate, a full drift matrix — now say plainly that they are not
  implemented yet and name the modifier actually given, instead of referring to an internal "slice".
- **The API reference is one page per level**, and each entry is stamped with what kind of thing it
  is — module, class, function, method, attribute — beside its heading and again in the contents
  column, so the reference reads as an index of each package rather than one long scroll. The
  overview page gains a table of the entry points: which function starts a run at each level and
  each genome resolution, what it returns, and which chapter covers it.
- **The docs site has the favicon the main site has been using**, rather than the theme's default
  mark.

## [0.20.0] - 2026-07-31

### Added
- **`zombi2 traits --name NAME`** — a run directory holds one slot per level, which is wrong for
  traits: a tree can carry several, and now one can drive another. `--name` writes each to
  `traits/NAME/`, so a driver's files are still there when the trait that reads them is written, and
  the run report gives each its own section. Without a name the plain `traits/` slot is unchanged.
- **A trait can drive another trait.** `DrivenBy` now works on a trait rate, in both the continuous and
  the discrete engine, so one trait's state can set another's variance-rate or switch rate. Driving is
  one thing driving another: with both participants at the same level it is still conditioning, and
  still the same modifier. A driver that switches **mid-branch** breaks the integral where it switches
  rather than being sampled once per branch, which is the difference between the model asked for and a
  2.8× wrong answer on a typical branch.
- **A discrete-bin (rate-category) clock** — `FromParent(spread=σ, bins=N)`. The rate takes one of `N`
  values on a geometric ladder and a daughter moves to a neighbouring rung, which is what a
  rate-category model assumes. It is a knob on `FromParent` rather than a modifier of its own because
  the model is `FromParent`'s: a daughter starts from its parent and is perturbed. `bins` defaults to
  `None`, the continuous form, so a run written before this draws exactly as it did.
- **Four gallery cards comparing the relaxed clocks** — uncorrelated lognormal, uncorrelated gamma,
  autocorrelated and discrete-bin, all down one species tree at one calibrated divergence and on one
  shared colour scale, so the only thing that differs between them is the pattern of rate variation.
  They replace the two phylogram examples.

### Changed
- **Conditioning and joining are two chapters**, not one. They answer the same question — can the driver
  be grown first? — but a reader following one does not need the other, and the single chapter kept
  putting the two side by side where they had to be told apart.

### Fixed
- **`--write` offers everything a result can write.** Its choices were hand-copied from each
  `write()`, with a comment saying so, and they had drifted: `species_tree` and `initial_sequence`
  were writable from Python and unnameable on the command line. Each result now declares its write
  vocabulary and the CLI reads it, so the two cannot disagree.
- **A run with an autocorrelated clock reported itself as `strict clock`.** The summary line only
  recognised `ByLineage`, so `FromParent` — and now its binned form — fell through to the default.

## [0.19.0] - 2026-07-31

### Added
- **A gallery example that is not the species→genomes→sequences pipeline** — four thousand trees under
  constant-rate and diversity-dependent birth, compared by their γ statistic, which is the other way to
  use a level: run it many times and measure the output.
- **`rearrangement_events.tsv`** — inversions, transpositions and translocations get a file of their
  own, at both the ordered and the nucleotide resolution, written with the event log. They begin and
  end no lineage, so in the genealogy table they were rows with nine columns empty; here every column
  is one they have. `position` is now `start`, because a `length` travels with it.

### Changed
- **The ordered `genome_events.tsv` is one row per event**, in the shared genealogy format —
  `time kind family parents children` with the coordinates of the event beside them. Participants are
  written `n<species>_g<copy>`, so the `lineage`, `recipient`, `donor` and `dest_lineage` columns are
  gone: a transfer is one row that names the donor's copy, the copy that arrived, the arc that left
  and where the block landed. **This changes the file's columns**, and a replacing transfer is now
  `transfer_replacing` carrying the copy it overwrote as a second parent, with no separate `loss` row.
- **`block_events.tsv` and `chromosome_events.tsv` carry `parents` and `children` too**, with copies
  written `n<species>_g<copy>` and chromosomes `n<species>_c<chromosome>`. A nucleotide speciation is
  one row naming both daughters rather than one row each; the chromosome log has no `lineage` column.
- **A replicon or a copy the run starts with is `initial`, not `origination`**, in
  `chromosome_events.tsv` as it already was in `block_events.tsv` — so counting `origination` gives
  what arose during the run.
- **"Wired" is gone from the error messages.** A level that refuses a modifier now says whether the
  combination is meaningless (a species tree has no gene families for `ByFamily` to vary over) or simply
  **not implemented yet** — where "not wired" meant the second and sounded like the first. The manual,
  the specification and the `RATES` help block use the same wording.
- **Long lines in the manual's code blocks wrap instead of vanishing.** Past about 104 characters a line
  ran off the paper and the remainder was silently not printed.
- **The gallery's sections fold**, with species open, and **conditioning and joining are now two
  sections** rather than one: conditioning grows the driver first and holds it fixed, joining grows both
  in one run. Deep links (`gallery.html#conditioning`) open a closed section.
- **A loss never takes an ordered chromosome below its last gene.** A run covering every gene still on
  the chromosome no longer fires, matching the floor the nucleotide resolution already enforced, so the
  two resolutions agree on what a chromosome is. Emptying the karyotype remains chromosome loss's own
  event. The family resolution has no floor: a high loss rate still empties a genome there. **This
  changes results** for ordered runs whose loss extents reach whole chromosomes — the random stream is
  untouched, so a run in which no loss is ever declined is byte-identical.
- **`--help` is orientation again, not a manual** — 836 rendered lines across the ten screens down to
  541. Every option states what it does and stops: the run directory is "the run directory", `--quiet`
  is "no progress bar". Units, defaults, valid values and genuinely surprising interactions stay; the
  reasoning moves to the manual, which already carried all of it.
- **`zombi2 tools format` requires `--format`.** It used to default to `homology` and write a table per
  family, which answers a question nobody asked. Running it bare now names the three choices.
- **A single-file output no longer gets a directory of its own** — the marker table is
  `genomes/markers.tsv`, not `genomes/markers/markers.tsv`. `homology` and `recphylo` write one file per
  family, so they keep theirs.
- **Every "wrote" line prints one complete path** you can copy — `wrote 1 table in out/genomes/markers.tsv`
  rather than a base directory plus a fragment to compose yourself.
- **The saturation warning is a sentence** (it was a paragraph): what happened, the measured identity,
  and the flag that fixes it.
- **`zombi2 -h` no longer has a "Coupling" heading.** The commands group as **Levels / Joint / Tools**,
  and the trailing paragraph about conditioning is gone. ZOMBI2 simulates evolution at four levels, and
  `joint` is the command for one of the three ways those levels can relate — independent, conditioned,
  joint — rather than a level of its own. Conditioning has no command: a driven rate is written on the
  level command that reads it, like any other rate. The word "coupling" is retired throughout, in the
  CLI, the manual, the specification and Appendix B's `Conditioning and joining — no new files`.

### Fixed
- **Four `--help` examples did not work.** Two `zombi2 species` examples used a seed whose run goes
  extinct before the present, so they exited 1; the only worked rate example on `zombi2 joint` used a
  filename driver, which `joint` always rejects; and `zombi2 tools -h` showed `tools format` without a
  `--format`.
- **`--genes` claimed it was ignored alongside `--gff`.** Passing both is an error, and always was.
- **`--write summary` was undiscoverable** on `genomes`, `sequences` and `traits` — the summary is
  written by default but was missing from all three help listings. The lists are now generated from the
  output tables, so they cannot drift again.
- **`zombi2 genomes` described two resolutions**, omitting `nucleotide`, in both the command list and
  its own description.
- **An error from a `zombi2 tools` subcommand named itself with the parent's whole usage string** —
  `zombi2 tools <tool> DIR [options] format: error: …`.
- **The `zombi2 joint` example in `--help` named a flag that does not exist** — `--families` for what is
  spelled `--family-names`, so the example failed if you copied it out of the help.
- **Gallery: extinct lineages are drawn dashed again.** Since lineages that go extinct became `e<id>`,
  the gallery still built its set of names as `n<id>`, so it matched nothing and four figures — extinct
  lineages, mass extinction, state-dependent extinction and MuSSE — drew every branch solid; the two
  state-dependent figures also left those branches uncoloured. The names now come from the tree's own
  `labels()`, and the copy-paste snippets on the detail views no longer teach the broken line.

## [0.18.0] - 2026-07-30

### Added
- A **FAQ** (in the docs) answering the handful of things that surprise newcomers — the `n`/`e`/`g`
  ids, the skipping family numbers, `_complete` vs `_extant` trees, how to view a tree, and the
  saturation warning.

### Changed
- **Each command now ends by listing every file it wrote** — one line per output directory naming the
  files in it (a per-family directory named as itself), then the run report. It was a two-line pointer
  at a couple of the files; the rest you had to find. The file-by-file descriptions live in
  `run.zombi2`, which the terminal names from the same source, so the two never disagree.
- **`trait_values.tsv` gains a `kind` column** — a tip's fate (`extant` / `extinct`, and `unsampled`
  under incomplete sampling) or `ancestor` for an internal node — so the observed tips a comparative
  method wants are one `kind == "extant"` filter away from the all-nodes table.
- **A bare nucleotide genome now starts from a 10 kb replicon with ten 500 bp genes** (it was a 1 kb
  all-intergenic replicon), so a first `--resolution nucleotide` run has real genes — and so gene
  trees, GFF/BED and a downstream `sequences` run — without having to know `--genes` and `--gene-length`
  first. The count is capped to what fits, so a smaller `--root-length` still runs.
- **A bare run no longer warns about its illustrative defaults.** The values it fills in are still
  recorded — in the `.log` and `run.zombi2` — but the paragraph-long stderr warning on every default
  run is gone.
- **The TO REPRODUCE block lists only the flags you actually typed** (plus the seed), so a run left on
  its defaults reproduces as the short command it was — and a family/ordered sequences run no longer
  prints nucleotide-only knobs (`--intergene-speed`) that never applied.

### Fixed
- **The bundled E. coli example ran with stale option names.** `examples/parameters/ecoli.toml` used
  the old `*-length` keys for the seven event extents (the CLI renamed them `*-extent`), so
  `zombi2 genomes --params parameters/ecoli.toml` refused the file with an unknown-parameter error.
  Its `[sequences]` step also gained `divergence = 0.2`, so the example's alignments no longer saturate.
- **`run.zombi2` claimed "machine-readable stats" for genome resolutions that write none.** The
  ordered and nucleotide resolutions drop only a `.log`, not a `_summary.json`, so their `records:`
  line now reads `(run parameters)`; the label appears only when a summary is actually present.
- **`--write summary` was a silent no-op on `--resolution ordered`** (it writes no summary file). It is
  no longer an accepted output there, so asking for it is a clear error rather than a file that never
  appears.
- **The log's `command_line` is now shell-quoted**, so a recorded rate expression
  (`--birth "1.0 * OnTime({0: 1.0})"`) pastes back into a shell intact rather than globbing on the `*`.

## [0.17.0] - 2026-07-30

### Added
- **A run report, `run.zombi2`.** Every command now writes one human-readable page at the run root,
  in the spirit of IQ-TREE's `.iqtree`: the parameters each level ran with, the input files it was
  computed on (by content hash), what came out, and every file it wrote — ending with a
  `TO REPRODUCE` block, one command per level rebuilt from the resolved parameters, that reruns the
  pipeline byte-for-byte. Because it records inputs by hash, it warns when a downstream level was
  computed on an upstream one that has since changed. It is a derived view of the per-level logs and
  summaries, rewritten on every run, so it always reflects the run as it stands; a `--flat` run keeps
  no per-level records and so has none. The run log now also records the Python, NumPy and platform
  versions — the environment a reproduction or a bug report turns on.

### Changed
- **`trait_values.tsv` now carries every node** — the extant tips, the extinct lineages and the
  internal nodes alike, each with its exact value, the same ones `trait_tree.nwk` annotates. It was
  the extant tips only. (The Python `TraitsResult.values` property is unchanged: still the extant-tip
  comparative dataset.)

### Fixed
- **Three manual tour/reference examples** a first-time reader hit now run as written (#287).

## [0.16.0] - 2026-07-29

### Added
- **The extant-only reconciliation** — `zombi2 tools format DIR --format recphylo --recphylo extant`
  writes the simulated history projected onto what a *dataset* holds: the extant gene tree inside the
  extant species tree. The complete recPhyloXML is the truth and is not what a reconciliation method
  should be scored against, because a method never sees it.

  The projection keeps what is observable and drops what is not. A speciation where the gene followed
  one daughter becomes **one** loss on the other, standing for however many really happened down
  inside it, and disappears entirely when that daughter has no surviving descendant. A duplication
  whose second copy died disappears: both copies sat on one branch. A transfer whose donor left no
  survivor leaves the arriving copy appearing from nowhere; where it rejoins its surviving relatives
  *is* the transfer, attributed to the branch they share.

  Two files come out per family and the gap between them is the point. `_true` is rooted where the
  family really originated, ancestral presence and the losses that narrowed it included — the answer
  key for ancestral gene content. `_recoverable` is that trimmed back to the surviving copies' common
  ancestor — the ceiling any method could reach. The first *contains* the second, so a scoring script
  can grade either way from it, and only one way from the other. `family_origins.tsv` records how each
  family entered (`origination`, or `transfer` from a lineage nobody can see), which the
  reconciliation itself cannot show and recPhyloXML has no tag for.

  Verified with a checker written against the trees rather than the builder — every node on a branch
  the extant tree names, every speciation and duplication child at or below its parent, leaf set
  identical to the extant gene tree — over three regimes and 71,000 nodes. It earned its keep: it
  caught the transfer rule twice.
- **`zombi2 sequences --stream`** writes each family's files as it finishes and keeps nothing. This is
  the level where a run's memory goes — every alignment and every ancestral sequence live at once — so
  what you can run was bounded by families × copies × sites rather than by time. At 200 species and
  640 families: **466 MB in memory against 329 MB streamed at 300 sites, and 1008 MB against 392 MB at
  1500**. The in-memory run grows with the sequences produced; the streamed one barely moves. It is a
  memory choice and not a modelling one — the same seed leaves byte-identical files either way, on the
  serial engine and on a real worker pool. Not available on a nucleotide run, which reassembles whole
  genomes and so needs every block at once.
- **`--strict`** on `species`, `genomes` and `sequences` refuses to fall back on an illustrative
  default. A bare run still works and still warns, because that is what teaches a newcomer the shape
  of a command — but in a pipeline that warning goes to a log nobody opens, so a dropped config key or
  a `--params` file that failed to apply gave a *successful* run whose science was not the one asked
  for.
- **A genome result writes the species tree it ran along.** `result.write()` wrote gene trees,
  profiles and events — every one of them indexed by that tree's node labels — and not the tree, so
  the Python quickstart handed back a dataset with no ground truth in it and said nothing. The CLI
  still keeps one canonical copy under `species/` rather than one per level.

### Changed
- **`max_family_size` is a plain number of copies in one genome** (breaking). It limits copies of one
  family in a *single* genome, but was written with a scope, and `PerLineage(n)` resolved to `n` ×
  the node count of the *species tree* — so the shipped default of `PerLineage(10)` was **1470 copies
  in one genome on a 147-node tree**, and more on a bigger one. A per-genome bound that moves when you
  add species is not a bound the person who set it can predict. Both old spellings are now refused
  with the arithmetic each implied, so a stale script fails loudly instead of running at a cap
  different from the one it reads; a float is refused rather than rounded. The default is **10**,
  which changes nothing for realistic rates — a tuned bacterial run of 3000 families is byte-identical
  at 10 and at the old effective 1470 — and constrains exactly the runs that were producing genomes
  with eleven copies of a family.
- **The `n_extant` stopping rule is written down.** The run stops the first moment that many lineages
  are alive together, then advances the clock by one more waiting time. Under pure birth that is
  exactly the general sampling approach (Hartmann, Wong & Stadler 2010, whose entry had been sitting
  in `references.bib` uncited); **with extinction it is a first-hitting rule and the trees are
  shallower** than the birth–death process conditioned on that many tips — measured at about a tenth
  of the tree height at 10 tips with `death/birth = 0.4`, a third to a half at 0.8, and back in the
  noise by 50 tips. Anyone estimating rates from these trees needs to know that, so the docstring and
  Chapter 3 now say it.
- **"Time runs forward from the origin", not "from the crown"** — sixteen places in the code said
  crown where they meant the origin, contradicting the book's own definition two chapters earlier and
  the specification. The crown is the first split; `t = 0` is the origin, at the top of the stem.
- Appendix A no longer claims that mean-correcting a random modifier leaves the tree alone. It fixes
  the per-split *factor*, not `E[N(t)]`: a branching process is convex in its rate, so standing
  diversity **rises** with `spread` and a run can grow explosively.
- `zombi2 tools format` prints its summary line under `--quiet`, like every other command. It alone
  printed nothing, so a scripted run had to go looking on disk to learn whether it had written
  anything.
- `numpy` and `tqdm` gained upper version bounds, so a future major release cannot break an install
  of a pinned ZOMBI2.

### Fixed
- **`--stream` crashed on every genome run.** The log's `write` line read a variable only the
  in-memory branch set, so a streamed run printed its success line and *then* died with
  `UnboundLocalError`, exit 1, having written no `genomes.log` — a `set -e` pipeline stopped after
  apparently succeeding. It shipped broken because nothing in the CLI tests mentioned the flag.
- **A run that sampled none of its survivors exited 0.** `sampling` can observe none of them, and the
  run then has no present: `zombi2 species` wrote no `species_extant.nwk` and exited 0, after which
  `zombi2 genomes` produced "0 gene families across 0 extant genomes", also at exit 0. In a batch of
  replicates that is a silently empty dataset (about 1% of seeds at `--sampling 0.2 --n-extant 20`).
  It is the same dead end the extinction guard already refuses, so it is refused the same way.
- **`--force` deleted the downstream run in silence.** Its notice was suppressed by `--quiet` and went
  to stdout. Every pipeline passes `--force`, because a re-run re-enters a directory that holds its
  downstream, and most pass `--quiet` too — so the one command that deletes work you already have was
  muted for exactly the caller who needed to hear it. It is a warning now: stderr, and not silenceable.
- **Extant gene trees named an extinct lineage `n<id>`.** A transfer node sits on the *donor's* branch
  and survives into the extant tree whenever the copy that moved has extant descendants, so a transfer
  out of a lineage that later died was written `n30` where the complete tree says `e30` — a label
  resolving in neither species tree file, on precisely the ghost-donor transfers a reconciliation
  benchmark is about.
- Ctrl-C printed a twenty-line traceback where every other error path gets a sentence. It is
  `zombi2: interrupted` and exit 130.
- The species run summary counts unsampled tips, so its parts add up to the tip count printed beside
  them. `0 extant + 8 extinct (28 tips)` was arithmetic that did not work.
- `zombi2 tools -h` advertised the homology table as `O`/`P` after 0.15.0 renamed it to the event, and
  `tools format -h` ran its three format descriptions together into one sentence. The same help said
  "Two `--format` choices" while listing three.
- `--help` cited "SPEC §5" at users. The specification does not ship with the package, so the RATES
  block explains itself instead.

## [0.15.0] - 2026-07-29

### Added
- **A marker table** — `zombi2 tools format DIR --format markers` writes one row per family saying
  whether it can be trusted to recover the species tree: `single_copy`, `universal`, the family's own
  event counts, and an `rf` distance between its true gene tree and the species tree **restricted to
  the genomes it occupies** (`rf = 0` is a perfect marker). This is the answer to the question most
  people are actually asking when they ask for "the orthologs", and it is a question about a *family*,
  which no pairwise label adds up to.

  Its value is the case that is invisible in real data: **a family can be single-copy *and* universal
  and still give the wrong tree** — hidden paralogy from a duplication plus reciprocal loss, or a
  transfer that replaced the resident gene. On a transfer-rich run, 111 of 299 families came out
  single-copy and universal and 106 of those did *not* recover the species tree. `rf` agrees with
  `zombi2 tools treedist`, checked against the command itself on real families.

### Changed
- **The homology table states the event, not a reading of it**: `S` (speciation), `D` (duplication)
  or `T` (transfer) at the pair's common ancestor, each optionally with `x` for a transfer *since*.
  It used to say `O`/`P` — ortholog and paralog — and that was worth undoing. The readings disagree
  with each other and most of them need something a pair of genes does not carry: Fitch's relation is
  not one-to-one, paralogy is relative to a chosen speciation, and the label even depends on the
  reconciliation model, since a duplication–loss model has no transfer category and explains one as a
  duplication plus losses. Each is *this event plus a choice*; ZOMBI knows the event exactly and the
  choice belongs to the reader. It also retires the absurdity of two genes in one genome being
  labelled orthologs — that pair now reads `Sx`, which is simply true.
- The Tools appendix gains an **"If you came here looking for orthologs"** section setting out the
  three different things people mean by the word and which of them ZOMBI2 answers.

## [0.14.0] - 2026-07-29

### Changed
- **A lineage that went extinct is now named `e<id>`** — in the complete species tree, the event
  logs, `genomes.tsv`, the complete gene trees, the trait tree and event log, the ancestral FASTA and
  phylograms, recPhyloXML, and the gff/bed and assembled-genome filenames. The number is the identity
  and the letter an annotation: `5` names the lineage wherever it appears and a join can always strip
  the prefix, which is what keeps this a marking rather than an identifier change.

  It marks the one fact about a branch you cannot recover from the tree's shape, and it means a
  complete tree **states its own extinctions**: the file survives being moved, copied or emailed,
  where the sibling `species_fates.tsv` does not, and `read_newick` no longer needs that table or a
  guess from tip depth to recover them.

  Two things deliberately stay `n`. Internal nodes — a speciation is not a fate. And **unsampled**
  tips, which are alive: being unsampled is a property of the sampling you asked for, not of the
  lineage, so the same branch would be named differently by two runs of the same tree.
  `species_fates.tsv` remains the only thing that tells an unsampled tip from an extant one.

  **Nothing that names only extant tips changed** — the extant tree, `profiles.tsv`, the alignments,
  the homology tables and the extant gene trees are `n<id>` throughout, so the ALE-ready trees, the
  FASTA joins and `treedist` are all unaffected.

### Added
- **`--parallel` runs conditioned rates.** A `DrivenBy` rate, a driven `transfer_to` and a `Clades`
  recipient rule used to make the parallel engine announce a fallback and hand the run to the serial
  loop — so anyone iterating on a coupled model was capped at one core, which is exactly the kind of
  run that is worth parallelising. **Conditioning does not couple families**, which is what makes the
  per-family decomposition survive it: the driver was grown *before* this run and is an input to it,
  so a lineage's factor at a given moment is the same number whichever family is asking, and no
  family can reach another through it. The workers now thread the driver trajectories with the rest
  of the run context, and the driven rate is summed over the family's own footprint rather than the
  whole tree. Measured at 20,000 families with a driven loss: **4.3 s serial → 2.7 s on four cores,
  2.5 s on eight** — a smaller margin than it would have been a day ago, because the family-cap fix
  above took most of the serial cost out first. Streaming (`stream_to=`) takes a driven rate too,
  where it used to raise.

  Verified by compensator, not by comparison: the parallel engine originates and loses at the rates
  it was *declared* with, to within a standard error, on the same standard the serial engine is held
  to. Worker-count invariance survives the driver.

### Changed
- **The genome level is no longer superlinear in genome size.** `max_family_size` asks "does this
  family already fill its quota here?" on every duplication and every arriving transfer, and it was
  answered by scanning the lineage's whole genome — an O(genome) step in the inner loop of a run
  whose genome is the thing that grows. It was 70% of a 4000-family run, and it was the *entire*
  superlinear term: with the cap removed the same run was linear. A per-lineage counter answers it by
  lookup, and the cost per event is now flat (2.0 µs at 1000 families and at 8000). Byte-identical —
  integers in, integers out, so the cap binds exactly where it bound before; verified across 48
  configurations of cap, replacement, transfer and self-transfer.

  | families | before | after |
  |---:|---:|---:|
  | 1,000 | 0.59 s | 0.36 s |
  | 4,000 | 6.44 s | 1.24 s |
  | 8,000 | 29.33 s | 2.45 s |

### Fixed
- **`--params` on a missing file printed a raw `[Errno 2]`**, and a broken one printed tomllib's
  position with no mention of the file. Both now name the file, and a backslash inside an ordinary
  TOML `"…"` string is called out for what it is — a path TOML has eaten before ZOMBI2 saw it.
- **An unknown flag now suggests the real one** — `--transfers` → `did you mean --transfer?` — and
  only when the guess is close, so a flag unlike anything gets no invented suggestion.
- **The run log records the options the run's resolution actually has.** A family run logged
  `root_length`, `gene_length`, `inversion` and the rest of the nucleotide and ordered knobs at their
  defaults, which reads as though it had them and chose those values.
- **Mean alignment identity is reported to a tenth of a percent**, not rounded to a whole one.

### Added
- Appendix B documents the gene-tree **internal node labels** (`duplication_n45` — the event that
  ended that gene and the branch it was on), names the viewers a `.nwk` opens in, and gives
  `zombi2.tree.read_newick`, which had appeared in no user-facing page. `read_newick` and `as_tree`
  are now on the API reference too.

## [0.13.0] - 2026-07-28

### Changed
- **The homology table says whether transfer is in a pair's history, not only whether it was the
  divergence.** A cell now carries two independent axes: how the pair diverged (`O` a speciation, `P`
  a copying — a duplication or a transfer, both of which turn one gene into two) and whether a
  transfer sits anywhere on the path since (`x`). So `O`, `P`, `Ox`, `Px`. Reading xenology off the
  common-ancestor event alone caught only the copy left behind against the copy that left, a sixfold
  undercount — and called the case that matters most a plain ortholog: two genes in the same genome,
  one an arrival from a relative, which no orthology method can reproduce, so the answer key scored
  the method wrong. All 650 such pairs in a measured run now read `Ox`. The table is also now read
  off each family's **complete** gene tree: pruning suppresses a transfer whose donor-side copy left
  no survivor, which was a fifth of all cells.
- **`zombi2 tools treedist` compares a gene tree to a species tree**, on the species each gene sits
  in, and says on stderr that it has — the two are not the same kind of object and used to fail as
  "different leaf sets". Refused, naming the genomes at fault, when the family is not single-copy.

### Added
- **CI runs the suite on Windows and macOS**, not only Linux, so the install instructions can say so.

### Fixed
- **A driver path could not be written on Windows.** The path sits inside a rate expression, read by
  Python's own parser, so `DrivenBy('C:\Users\me\trait.tsv', …)` failed as a truncated `\UXXXXXXXX`
  escape — and `C:\temp` was worse, parsing silently into a tab. A path now means itself, however it
  was written, escaped or not.
- **Every file was opened in the platform's default encoding**, which on Windows is cp1252 — so a
  UTF-8 tip label in an input tree was a decode error there and nowhere else. All 215 text reads and
  writes are now explicitly UTF-8.
- **348 Sphinx roles published verbatim to the docs site.** `docs/reference/api.md` is generated from
  the docstrings by mkdocstrings, which renders markdown and knows nothing about `:func:` — so they
  appeared as raw markup on the page a Python user is most likely to read. Two dead module names went
  with them (`zombi2.species_tree`, `zombi2.modifiers`).
- The manual's index said two appendices and listed two; there are three.

### Added
- **recPhyloXML output** — `zombi2 tools format DIR --format recphylo` writes each family's complete
  gene tree inside the complete species tree, in the community format for a gene tree embedded in a
  species tree. It is normally what a reconciliation *method* produces; here nothing is reconstructed,
  so the file is the true history and can be used as the answer key an inference is scored against.
  Losses, extinct species and transfers-from-the-dead are all in it, which is why the *complete* trees
  go in. Output validated by round-tripping through the format's own reference parser. In Python,
  `zombi2.tools.recphylo.recphylo_xml(gene_trees, tree)`.
- **Every run records the files it read, by content.** One `input<TAB><sha256><TAB><path>` line per
  input in the run log — the species tree, a `--tip-fates` table, the driver file of a conditioned
  rate. A path alone does not identify an input: two runs from two different trees could log the same
  `--from tree.nwk`, the same seed and the same parameters, and nothing told them apart.
- **A stated reproducibility policy**, in the manual's new *Reproducing a run* (Chapter 2): what a
  seed pins down within a version, why it is not promised across versions, and why a parallel run
  differs from a serial one.

### Changed
- **`result.write()` groups its many-files-per-run outputs, as the commands already did.** A genome
  run written from Python dumped 200 gene-tree Newicks into the directory alongside its four tables;
  the same run written by `zombi2 genomes` put them in `gene_trees/`. The subdirectory is now the
  *result's* decision, not the command's, so both layouts agree: `gene_trees/`, `gff/`, `bed/`,
  `alignments/`, `phylograms/`, `ancestral/`, and the assembled genome FASTAs in `genomes/`.
  `result.write(dir, flat=True)` is the old behaviour, and is what `--flat` now passes through. An
  output a run has none of writes nothing and leaves no empty directory. **A script reading files
  written by `.write()` will need the subdirectory in its path**; nothing written by a `zombi2`
  command moved.
- **`--family-speed` and `ByFamily` are no longer quadratic in genome size.** The per-lineage sums of
  the per-family multipliers were rebuilt from scratch on every event — the whole live gene pool, per
  event — where one event changes one lineage by one copy. They are now carried across events and
  only the touched lineage is rebuilt. At 1000 families a `--family-speed` run went from 157× a plain
  run to 7.9× (86.9 s to 4.0 s); at 300 families, 62× to 3.4×. Byte-identical: every seed gives the
  run it gave before, verified across seven configurations.

### Fixed
- The saturation warning said `--substitution` "is currently None" on a run that had defaulted it. It
  now names the rate the run actually used.
- **A result printed its whole contents.** Typing a result's name in a session or a notebook rendered
  the dataclass repr — 4.5 MB for a 40-tip genome run. Every result and the tree now repr as a
  one-line summary of what the run produced.
- **Handing a level the wrong thing raised an `AttributeError` from inside an engine.** Passing a
  Newick *string* or a *path* where a tree belongs — the easy mistake, since both are a tree to the
  person holding one — now says so, and says how to parse it. One shared `as_tree` at every level's
  entry, so what a level accepts is decided in one place.
- **The package shipped no `py.typed` marker**, so a type checker in someone else's project ignored
  every annotation in ZOMBI2 and treated the whole package as `Any` (PEP 561).

## [0.12.0] - 2026-07-28

### Added
- **`divergence` on the sequence level: state the outcome, and the rate is solved for.** The
  substitution rate is per unit *time*, so what it produces depends on the height of the tree it runs
  down — `1.0` is reasonable on a short tree and pure noise on a tall one, which is why no default can
  be right. `divergence=0.2` instead asks for 0.2 substitutions per site from root to tip, and the
  base falls out as `divergence / height`. Verified across a tenfold difference in tree height: the
  same `divergence` gives alignments within 5% identity of each other, where the same *rate* would
  not be close. It composes with `substitution`, which keeps saying what *kind* of clock: give the
  shape alone (`ByLineage(spread=0.3)`) and `divergence` sets its scale. A base number alongside is
  refused rather than overridden, and the run log records the resolved rate, so a run stays
  reproducible from its own command line.

### Fixed
- **A relaxed clock reported itself as a strict clock** in the `zombi2 sequences` summary line
  whenever the modifier was given without a base.

### Added
- **`zombi2 sequences` reports the alignment divergence it actually produced, and warns when the
  sequences are saturated.** The summary line gains `mean identity NN%`, measured over a bounded
  random sample of within-family pairs. When that identity sits within 15% of the model's own random
  floor (`Σπ²` — 25% for equal-frequency DNA, ~6% for a protein model) the run also prints a warning
  naming the realised identity, the floor, and `--substitution`. The rate is per unit *time*, so
  whether it yields a usable alignment depends on the height of the tree it runs down — which no flag
  reveals. Five independent first-time users were given the default rate on the documented quickstart
  tree and two of them nearly abandoned the tool believing the sequence level was broken. The warning
  goes to stderr and survives `--quiet`, so a scripted batch hears it; a healthy run prints nothing.
- **A "Genome reduction" gallery example — a trait *conditions* the genome.** In the renamed
  *Joining and conditioning* level: an irreversible endosymbiont lifestyle drives fast gene loss and
  near-zero gene gain through the same `DrivenBy` mechanism that couples a trait to speciation, so
  those lineages' genomes collapse — shown as per-tip genome-size bars beside a tree coloured by
  lifestyle. Uses Phylustrator 0.1.2's new `bars` panel (pinned via `zombi2[gallery]`). (#260)

### Added
- **`zombi2 species out/`, `zombi2 genomes out/` and `zombi2 sequences out/` now run with no other
  arguments.** Each fills what it was not given (`--birth 1.0 --n-extant 20`; `--duplication 0.2
  --transfer 0.1 --loss 0.25 --origination 0.5`; `--model jc69`) and says on stderr exactly which
  values it chose and that they are illustrative rather than estimates. Refusing to run taught a
  newcomer nothing about the shape of a command, and the rates are precisely the part nobody can
  guess on a first run — "every number was a guess" was among the loudest first-time complaints. A
  default is announced rather than silent for the same reason a saturated alignment now is: the run
  should never imply a number was chosen when it was not, and the `.log` records it either way.
  `genomes` fills in only when **no** rate of any kind was given: a run that names `--duplication
  0.3` or `--inversion 1.0` has had its model described, and quietly adding gene turnover to it
  would be the surprise.

### Changed
- **`max_family_size` is written with a scope**: `PerLineage(10)` (the new default, ten copies for
  every lineage in the complete tree) or `Global(50)` (fifty copies whatever the tree looks like);
  `None` still removes the cap. It used to be the `int`/`float` distinction alone — `10` absolute
  against `10.0` relative — which put a factor of the tree's size between two values Python calls
  equal, in a spelling no `--params` file or reader could be expected to notice. A bare number is now
  refused and names both forms. This is the same vocabulary rates already use for the same question,
  so it needs no new notation: the run log records it in written form (`Global(4.0)`) and it pastes
  straight back into the flag or a TOML. The scope wrapper also rejects a negative base at
  construction, so a validation branch disappeared with it.
- **A gene copy is named `n<species>_g<copy>` wherever there is no column to say which species it
  sits in** — gene-tree and phylogram Newick leaves, alignment and ancestral FASTA records, and
  homology table headers (which used `n<species>|g<copy>`; `_` replaces `|`, which aligners and tree
  builders treat as a field separator). Previously a FASTA record was `>g2179` and could not say
  which genome it came from, so anyone benchmarking orthology had to join the alignments back to
  `genomes.tsv` themselves — the first thing two independent first-time users each wrote by hand.
  The `g<copy>` half is unchanged and the tables still carry their own `lineage` column, so every
  file still joins on the same token: verified as 1254 FASTA records matching their gene-tree tips
  exactly, none mismatched.
- **Every resolution now writes the same `genome_events.tsv`.** One filename meant two different
  tables: at the nucleotide resolution the rows were ancestral intervals, a duplication wrote one row
  instead of two, `initial` replaced `origination`, and `lineage` on a transfer named the *donor*
  rather than the branch the copy is born on. The columns were close enough that making them match
  would have let the family reader accept one and build silently wrong gene trees. So the genealogy
  is now written in the shared format — derived onto the root-block partition, where a copy either
  covers a block in full or does not touch it, which is what makes a duplication a bifurcation — and
  the interval record keeps its own name, **`block_events.tsv`**. `zombi2 tools` and the sequence
  level now work identically across all three resolutions, and a new test asserts the header,
  per-kind row arity and transfer contract agree between them. The engine, the copy ids and the
  simulation output are unchanged: the genealogy is what the recovery already built to derive the
  gene trees, and was being discarded. `read_nucleotide_genomes` reads `block_events.tsv`.
- **Trees are written with 7 significant figures instead of 6.** Rounding accumulated along a
  root-to-tip path put the extant species tree outside the `1e-6 × height` ultrametricity that
  `read_newick` documents — measured at 19 of 20 seeds, so third-party dating and strict-clock tools
  could reject ZOMBI2's own output. At 7 figures the drift is ~2–3e-7 and 0 of 20 seeds exceed the
  criterion, verified up to 5000 tips. Applies to the species tree, gene trees, phylograms and trait
  trees; event logs were already full precision.
- **The quickstart is now one command sequence, identical in the README, on the landing page and in
  `zombi2 --help`.** The three had drifted apart — different `--n-extant`, missing rates, and no
  sequence step on the landing page. It now also passes `--substitution 0.05`, so the documented
  first run produces a usable alignment (83% mean identity) rather than a saturated one.
- **The landing page's version badge is stamped from `zombi2.__version__` at deploy time.** It read
  `v0.2.0` while PyPI served 0.11.0. The Pages workflow now rewrites it on every deploy and fails
  loudly if the badge is missing, so it cannot drift behind a release again.

### Fixed
- **Documented that `genome_events.tsv` holds one row per gene-tree *edge*, not per event.**
  Duplications, transfers and speciations each write two rows; losses and originations write one.
  Appendix B (which the output-files reference includes verbatim) previously stated this for
  transfers only, so counting rows by kind silently doubled the other two. Three of five first-time
  users hit this and two of them diagnosed it as a factor-of-two rate bug before finding the cause.
  The appendix now gives the per-kind row counts and a verified counting recipe — including why that
  recipe must *not* be applied to originations, which share `time` 0 and carry no `parent`.
- **`--stream` no longer implies it reproduces a serial run.** Its help said "the files are the
  same", meaning the same file *layout*; it reads as "the same content". It is a separate engine like
  `--parallel` and draws families in a different order, so the same seed gives a different (equally
  valid) run — verified as 173 families/3684 rows serial against 191/3996 streamed. The help now says
  so, in the wording `--parallel` already used.

### Changed
- **`zombi2[gallery]` now installs Phylustrator from PyPI.** The plotting library is published, so the
  gallery extra depends on it directly (`phylustrator>=0.1.1` — the release that draws transfers as
  visible arrows) instead of requiring a separate `pip install` from git. (#258)

## [0.11.0] - 2026-07-27

### Added
- **An examples gallery, published at `/gallery.html` and linked from the landing page.** Five levels
  — species, genomes, sequences, traits, and the joining that couples them — each a short, runnable
  recipe that simulates with ZOMBI2 and plots with Phylustrator (trees via `ph.trees`, genomes,
  synteny and alignments via `ph.genomes`). Click any figure for the exact code that reproduces it.
  The source lives in `gallery/` (self-contained, like `analyses/`) and a `[gallery]` optional extra
  regenerates the page; the landing page's four level cards now link to its per-level sections, with a
  "See the gallery" button beside "Read the docs". (#257)

## [0.10.0] - 2026-07-27

### Added
- **A joint trait can jump at the split as well as along branches.** `traits.discrete(at_speciation=…)`
  already worked inside `simulate_joint`, so a trait could drive speciation *and* change at each
  speciation — but nothing said so, nothing tested it, and Chapter 9 documented only the along-branch
  case, which made the combination look unavailable. It is now documented, with an example, and pinned
  by tests: a jump at a split is logged as `on_speciation`, a change along a branch as `on_branch`,
  so the two can be told apart when scoring a method. (#256)

## [0.9.0] - 2026-07-27

### Added
- **`--max-family-size` and `--family-speed` on `zombi2 genomes`.** Both existed in the Python API at
  the family and ordered resolutions but had no flag, which mattered once the growth cap became
  on-by-default at ordered: a command-line run had a bound it could neither adjust nor switch off. The
  cap keeps the model's int/float distinction — an int is an absolute copy count, a float that
  multiple of the tree's lineages — and takes `none` for unbounded growth. `--family-speed` takes a
  `ByFamily` draw in the same written form as every rate flag. Both are refused under
  `--resolution nucleotide`, each with its own reason rather than one blanket explanation. (#251)
- **An extent now takes modifiers, as SPEC §6 says it should: `base × modifiers`, no scope.** So
  `loss_extent = 150 * DrivenBy(habitat, {"host": 6.0, "free": 1.0})` makes host-restricted lineages
  delete in **bigger chunks** — a different model from `loss = … * DrivenBy(…)`, which makes them
  delete more **often**. Set both and they multiply. An extent takes the modifiers its resolution
  wires on a rate: `OnTime` at ordered, `OnTime` and `DrivenBy` at nucleotide; anything else raises.
  Unlike a rate's, an extent's modifier is read when an event fires, so it changes no rate and adds no
  step to the Gillespie clock. The concept has its own module, `zombi2.rates.extent` (`Extent`,
  `as_extent`), parallel to `rates.rate`; `as_extent` moved there from `rates.distributions` and now
  returns an `Extent` rather than a bare distribution. A scope on an extent is refused — it is already
  an absolute size, so there is no "per what?" to answer. (#250)
- **Conditioning at the nucleotide resolution: every rate there now takes a `DrivenBy`.** A trait can
  drive how much DNA a lineage sheds — genome reduction as it is usually meant, on a genome measured
  in base pairs rather than in family tokens — and can drive the **rearrangements** as well, which
  nothing could do before: Chapter 9 named "let a trait speed up inversion" as the limitation people
  hit first, and it is now expressible. Same mechanism and same written form as the family
  resolution, from a trait file or an in-memory result:
  `loss = 0.8 * DrivenBy(habitat, {"host": 20.0, "free": 0.5})`. A driven rate becomes per lineage, so
  the affected lineage is drawn with the same weights the total was summed with, and a driver
  switching mid-branch is a step the Gillespie stops at. `transfer_to` — where a transfer lands — stays
  family-resolution only; a nucleotide transfer's *rate* can be driven, its recipient rule cannot.
  (#249)
- **The nucleotide resolution now speaks the rate grammar.** Its rates were bare floats with the
  scope hardcoded; each is now a `scope(base) × modifiers` expression like every other level's, with
  the defaults stated rather than implied — **per lineage** for the gene events (the rate says how
  often a lineage does it, the extent how much DNA it touches, so the number reads the same whatever
  the genome's size) and **per chromosome** for the number-changing tier. A bare number stays a bare
  number, so nothing existing changes. The **skyline** now works there:
  `inversion = 5.0 * OnTime({0: 1.0, 3: 0.2})`, in Python, on the command line and in a `--params`
  file; the engine re-reads its rates at each step instead of racing past it. The resolution declares
  what it wires and refuses the rest — `DrivenBy` still raises, with a message saying so, rather than
  being silently dropped. (#248)
- **Per-family rate heterogeneity at the ordered resolution.** `ByFamily` on `duplication`,
  `transfer`, `loss`, `inversion`, `transposition` or `translocation`, and the family-wide
  `family_speed=`, now work at `--resolution ordered`; both were previously refused. The weight lands
  on the **segment an event covers**, not on the gene it started from (SPEC §6) — weighting the start
  would apply a family's own rate to its *neighbours*, and the neighbourhood is reshuffled by every
  rearrangement, so the parameter would not name a fixed thing over a run. A run's weight is the
  **mean** over the genes it covers, which is what makes a run with no weights set behave exactly as
  before. `ByFamily` remains refused on `origination` (there is no family yet to have drawn a factor)
  and on the chromosome tier (a fission acts on a whole replicon). (#247)
- **`max_family_size` now works at the ordered resolution too, and is on by default.** A segment may
  carry several families, and several copies of one, so a run is refused when it would take *any* of
  them past the quota — the whole run, never part of it, since clipping a run to the genes still under
  quota would quietly shorten runs exactly where the genome is crowded and so reshape the extent
  distribution. It reduces to the family resolution's condition when a run is a single gene. (#247)

### Changed
- **An ordered run is now bounded by default** (`max_family_size=10.0`, the same multiple of the
  tree's lineages the family resolution has always used), where it previously had **no cap at all**.
  Duplication compounds, so a family whose duplication rate sits above its loss rate — or one that
  drew a high `ByFamily` factor — could grow without bound; at `ByFamily(spread=1.0)` an ordered run
  produced roughly twice the duplications of the equivalent, already-capped family run. The two now
  agree. Pass `max_family_size=None` for the old unbounded behaviour. This **changes results** for any
  ordered run that would have exceeded the cap. (#247)
- **One word for how much a segmental event takes: `extent`.** The ordered resolution called it
  `<event>_extension` and the nucleotide one `<event>_length`; both are now `<event>_extent`, in
  Python and on the command line (`--inversion-extent`, `--loss-extent`, …). The unit is still set by
  the resolution — genes at ordered, base pairs at nucleotide — which is what a resolution is for.
  `--root-length` and `--gene-length` are **unchanged**: they size the initial genome, they are not
  extents. (#246)
- **A bare extent number is now the mean, at every resolution.** `duplication_extent=3` at the ordered
  resolution used to mean *exactly* three genes every time and now means runs averaging three, which
  is what the same number has always meant at the nucleotide resolution. This is a **change of
  results**, not just of spelling: an ordered run that passed a bare number will now produce a spread
  of segment sizes rather than one size. Write `Fixed(3)` for the old behaviour. `None` remains a
  single unit. The new `zombi2.rates.as_extent` coerces this; `as_distribution` is untouched, since a
  bare number is rightly a fixed value for the per-family rate specs that use it. (#246)
- **The nucleotide resolution refuses an extent shape it cannot honour.** It draws each arc's far end
  directly from the genome's legal breakpoints, so anything other than a geometric extent would have
  to be re-weighted over that set rather than drawn; passing one now raises instead of being silently
  approximated. An extent below 1 bp is also refused. (#246)
- **The manual's title page now names the ZOMBI2 version it documents** — `Version 0.8.0 — July 2026`,
  under the author. The version is read from the package at build time rather than typed into the
  manual, so the book cannot claim a release the code does not have, and the build stops if it cannot
  be read. (#246)

### Fixed
- Chapter 6 said an event whose extent the genome cannot supply "still fires, just shorter", which its
  own table contradicts two lines above — an arc whose nearest legal breakpoint lies beyond the extent
  you asked for comes out **longer**. The correction runs both ways, and the chapter now says so, along
  with the one degenerate case (a replicon with no legal end in reach) where the event is skipped. (#246)

## [0.8.0] - 2026-07-26

### Changed
- **The genome level's first resolution is now `family`, not `unordered`.** It is the default
  resolution and the first thing a new user meets, and it was named for the feature it lacks rather
  than for what it is. The other two resolutions keep their names, and the model itself is untouched —
  the three still share one engine, one rate grammar and one event log. This is a **breaking rename
  with no deprecation aliases**:

  | was | is |
  |---|---|
  | `simulate_genomes_unordered(…)` | `simulate_genomes_family(…)` |
  | `GenomesResult` | `FamilyGenomesResult` |
  | `UnorderedGenome` | `FamilyGenome` |
  | `genomes.unordered(…)` | `genomes.family(…)` |
  | `zombi2.genomes.unordered` | `zombi2.genomes.family` |
  | `--resolution unordered`, and the same value in a `--params` file | `--resolution family` |

  Renaming the result class also removes an asymmetry: the first resolution returned the unprefixed
  `GenomesResult` while its two siblings carried a prefix, which read as a default plus two variants
  rather than three resolutions of one model. The run log's `resolution` row and the run-completion
  summary now say `family`; no output filename, directory or column header changes, so nothing that
  reads a run's tables needs updating. Chapter 4 is now *Genomes I: gene families*. (#245)
- **`families=` is now `family_names=`**, on `simulate_genomes_family`, `simulate_genomes_ordered` and
  `genomes.family(…)`; `zombi2 joint --families` is likewise now `--family-names`. The argument takes a
  list of *names* and always did — the attribute it fills has been `result.family_names` all along — so
  input and output now agree, and `initial_families=100` (a count) no longer reads as a sibling of the
  name list beside it. (#245)

## [0.7.0] - 2026-07-26

### Added
- **Cross-level staleness guard.** A level refuses to re-run in place when a later level built from it
  is already in the run directory — re-running would leave that downstream output silently mismatched.
  This covers both the pipeline chain (re-running `genomes` would orphan the `sequences` under it) and
  `DrivenBy` **conditioning** (re-running a `traits` run that a `genomes` rate was conditioned on orphans
  that genomes, and the sequences beneath it — the dependency is recorded in a small `conditioned_on`
  marker). `--force` re-runs anyway and removes the now-stale downstream, so a run's levels can never
  quietly disagree. The forward pipeline is unaffected (each level is run once); applies to the default
  grouped layout — `--flat` commingles the levels and is left to the user. (#243)
- **Onboarding nudges.** After a level runs, a one-line pointer names the output file(s) worth looking
  at — a species run names both its complete and its extant tree — so a newcomer can find what they
  made. It says only *what was written and where*, never what to do next (a run is not one road);
  suppressed by `--quiet`. The top-level `--help` now leads with the levels and the plain quickstart,
  with the `DrivenBy` / joint coupling note moved below them. (#244)

### Fixed
- A non-nucleotide `sequences` run no longer leaves an empty `sequences/genomes/` directory behind (the
  assembled-genome FASTAs are a nucleotide-run output). (#244)
- The `sequences` run log now records the substitution model's **effective** parameters (e.g. `kappa 2.0`,
  `frequencies [0.25, 0.25, 0.25, 0.25]`) instead of the bare `None` args, so a run reproduces from its
  log alone, without the reader having to know each model's defaults. (#244)

## [0.6.0] - 2026-07-26

### Added
- `--quiet` on `zombi2 tools format`, matching the level commands — it suppresses the summary line,
  for a log file or a batch of runs. (#242)

### Changed
- **Gene-copy ids are now written `g<id>` in every genome table** — `genomes.tsv`, `genome_events.tsv`,
  `gene_order.tsv`, `blocks.tsv` and `initial_genome.tsv`, in both the `copy` and `parent` columns.
  This is the same token the gene-tree Newick leaves, the alignment FASTA headers and the homology
  tables already use, so a gene copy now joins across every file of a run with no translation step
  (previously the tables held a bare integer while the trees/alignments held `g<id>`). The column
  names are unchanged; this is a breaking change to the *values* those columns hold. (#242)
- **A joint (`zombi2 joint`) driver mapping that names a state outside the trait's declared alphabet is
  now refused.** A key that can never occur — a typo such as `{"caev": 4.0}` for a cave/surface trait —
  is caught up front instead of silently applying to nothing. The conditioned (file-driven) path is
  unchanged: it still refuses only a mapping that matches *none* of the driver's observed states. (#242)

### Fixed
- The run-completion line no longer doubles a trailing slash: a directory given as `out/` (as the
  quickstart shows) now echoes as `wrote out/`, not `wrote out//`. (#242)
- A conditioned rate whose `DrivenBy` points at a missing file now reports it as a missing driver file
  (with the path and how to fix it), instead of a bare `[Errno 2] No such file or directory`. (#242)
- `parallel=` no longer crashes with a raw `BrokenProcessPool` when called from a notebook, `python -c`,
  or a stdin heredoc (where worker processes cannot re-import your program). It falls back to
  single-process with a one-line note; a `.py` script or the CLI still uses every core. (#242)

## [0.5.0] - 2026-07-24

### Added
- **Opt-in parallelism** for the unordered-genome and sequence engines (`parallel=` in Python,
  `--parallel` on the CLI): a separate, worker-count-invariant engine that evolves independent units
  (gene families, gene trees) across processes. It covers D/T/L/O with `uniform`/`distance` transfer;
  a driven or clade-based `transfer_to` falls back to the serial engine (loudly). Also `stream_to=`
  to write each family straight to disk for the many-families regime. The serial path is the default
  and is unchanged. (#233)
- **`Clades` transfer_to rule + the `Between` kernel** — weight horizontal transfers by the donor's
  and the recipient's **named clade**, so transfer can be steered to run *between* two clades rather
  than within them (a topological rule, sibling of `"distance"`). Name a clade by a few of its tips or
  a node id. The trait-driven form is `mod.DrivenBy(trait, Between({...}))`. (#235)

### Changed
- **Faster nucleotide sequence runs at lower memory**: a run-wide CDF cache and lazy genome assembly,
  plus batched decoding and a streaming FASTA writer. (#234)

## [0.4.0] - 2026-07-24

### Added
- `analyses/` — a home for self-contained validation studies, each regenerating from fixed seeds:
  **RED** node-age validation (does Relative Evolutionary Divergence recover node ages under a
  realistic clock?) and the **yeast inversion-rate** study (inferring the genome inversion rate from
  synteny with the nucleotide genome model). (#227)
- `zombi2.tree` — one home for the `Tree` datatype and its toolkit, all free functions: `with_stem`,
  `make_ultrametric`, `rescale`, `relative_evolutionary_divergence`, `red_scaled`, `distance`
  (Robinson–Foulds / branch-score), and `read_newick(assume_extant=)`. (#228)
- CLI: `zombi2 tools tree` (prune / round / stem / rescale / RED) and `zombi2 tools treedist`. (#228)
- The **autocorrelated** molecular clock (`FromParent`) at the sequence level, alongside the
  uncorrelated `ByLineage`. (#228)
- Release tooling: this `CHANGELOG.md` and `scripts/release.sh`.

### Changed
- `Tree`, `Node`, `read_newick`, `to_newick`, `prune` moved from `zombi2.species` to `zombi2.tree`
  — import them from `zombi2.tree` now (`zombi2.species` keeps the simulator and `Event`). (#228)
- `docs/design/MAP.md` correctness pass: the quarantine now points at the sibling `ZOMBI2_LEGACY/`,
  the `joint` CLI is listed, sequences/traits are marked built, and stale signatures/sections fixed.
  (#228)

### Fixed
- `zombi2 tools treedist` matches tips by label, not by parse-order node ids, so two external trees
  are compared by taxon rather than by structural position. (#228)

## [0.3.0] - 2026-07-23

### Added
- The clean, pure-Python core grown from `docs/design/SPEC.md`: species trees, genomes at three
  resolutions (unordered ⊂ ordered ⊂ nucleotide), sequences, traits, and the joint engine — with the
  cross-level rate grammar (`zombi2.rates`) and a CLI. First release of the rewrite; `pip install`
  needs no build step.

### Removed
- The Rust engine (`zombi2_core`) and the old codebase, retired to the sibling `ZOMBI2_LEGACY/`.
  (0.2.0 was yanked on PyPI, superseded by this clean core.)
