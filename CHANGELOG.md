# Changelog

All notable changes to ZOMBI2 are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). While the project is pre-1.0, a **minor**
bump (0.x.0) carries new features or breaking changes and a **patch** bump (0.x.y) carries fixes.

A release is cut with `scripts/release.sh patch|minor|major` (the version is computed, not typed),
which moves the entries below from `[Unreleased]` into a dated version section.

## [Unreleased]

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
