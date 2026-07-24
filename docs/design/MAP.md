# ZOMBI2 — Code Map

**Status: AUTHORITATIVE. The code companion to [`SPEC.md`](SPEC.md).**

`SPEC.md` fixes the **model and the words**. This fixes the **files and the names**: every module, every
public name, and its one canonical home. When code, a docstring, the CLI, or a chapter names something
differently, this document wins and that other file is a **fossil**.

Read `SPEC.md` first (the model), then this (the shape). Two principles throughout:

- **One canonical path per name** — no top-level re-exports, no aliases, so nothing drifts. Reach
  everything through its level package (`from zombi2 import species`, `from zombi2.rates import scope`).
- **Concepts → code → chapter**, and this map is the **as-built truth**: where an older design doc
  disagrees, this map is right — the built engine uses
  `n_extant` / `total_time`, and time is always **forward from the crown** (SPEC §5).

---

## The clean core and the quarantine

The rewrite is a **clean core grown from `SPEC.md`**, not a migration of the old codebase. Everything not
yet rebuilt is **kept out of the clean core**, read-only, in the sibling **`../ZOMBI2_LEGACY/`** (the old
codebase, moved out of the repo — not importable, not wired to anything). Features are ported **in** from
there deliberately, one level at a time, renamed to this map — never carried, never re-imported. That is
what keeps the active tree small enough to hold in one head, and in one context window.

Legend:  ✅ built · 🔨 to build · 📦 not in the clean core (reference in the sibling `../ZOMBI2_LEGACY/`)

## The tree

```
zombi2/
  __init__.py          thin — no top-level re-exports (one canonical path per name)
  rates/               the cross-level rate grammar (SPEC §5): effective rate = scope(base) × modifiers
    scope.py        ✅ PerCopy · PerLineage · PerSite · PerChromosome · Global
    modifiers.py    ✅ OnTime · OnTotalDiversity · FromParent · ByLineage · ByFamily · DrivenBy   (Markov 🔨)
    rate.py         ✅ Rate · as_rate       (internal plumbing; users never build a Rate directly)
    parse.py        ✅ parse_rate           (the written form → a rate spec: the CLI and --params read the same expression as Python, SPEC §5)
    distributions.py✅ Fixed · Exponential · Gamma · LogNormal · Uniform · Geometric   (value / length distributions)
    mapping.py      ✅ Table · Curve · Scalar · as_mapping   (a driver value → a factor; DrivenBy's response — SPEC §2)
    driver.py       ✅ DriverTrajectory · load_driver   (a conditioned DrivenBy's file-backing: value/next-switch per lineage)
  tree.py            ✅ Tree · Node · read_newick · prune  (the shared dated-tree datatype + its toolkit — every level rides on it)
  species/           ✅ simulate_species_tree → SpeciesResult ;  Event  (the recorded true history)
  genomes/           ✅ simulate_genomes_unordered · simulate_genomes_ordered · simulate_genomes_nucleotide   (all three on the CLI via --resolution) — unordered: all four D/T/L/O DrivenBy-conditionable ✅ plus a driven `transfer_to` recipient weight ✅; ordered wires OnTime only; nucleotide takes constant rates
  sequences/         ✅ simulate_sequences → SequencesResult
  traits/            ✅ simulate_continuous · simulate_discrete → TraitsResult ;  discrete(...) process spec (for joint)
  joint/             ✅ simulate_joint → JointResult   (the FUSE engine; SPEC §2–4). Conditioned needs no engine — it folds into the target level via DrivenBy + rates/driver.py. "Coupling" is the concept (SPEC/manual), not a package.
  tools/             ~  read-back analyses on a finished run (the levels simulate; the tools re-express). `homology` ✅ — the true ortholog/paralog/xenolog matrix per gene tree (the event at each leaf pair's MRCA). (RED and tree distance are not here — they live in `zombi2.tree`, exposed on the CLI as `tools tree --red` / `tools treedist`.) The rest (reconciliation, recon-accuracy, …) 📦
  cli/               ✅ species · genomes · sequences · traits · joint · tools (clean, one flag per API keyword, rates in the written form)
```
(The old codebase is no longer in the repo — it lives in the sibling `../ZOMBI2_LEGACY/`; nothing under
`zombi2/` imports it.)

The cross-level primitives live in **`zombi2.rates`** (scope, modifiers, rate, distributions) because the
rate grammar is one thing shared by all four levels. The **`Tree`** and its toolkit live in
**`zombi2.tree`** — every level rides on the shared dated tree (the species engine and `read_newick`
produce it; genomes, sequences and traits consume it), so the datatype and everything you do to a tree
share one home (`from zombi2 import tree`).

## The names (as built)

| Canonical home | Public names |
|---|---|
| `zombi2.rates` | `from zombi2.rates import scope, modifiers` → `scope.Global`, `modifiers.OnTime({...})`. Scopes: `PerCopy · PerLineage · PerSite · PerChromosome · Global`. Modifiers: `OnTime · OnTotalDiversity · FromParent · ByLineage · ByFamily · DrivenBy`. Also `parse_rate("1.0 * OnTime({0: 1.0, 3: 0.3})")` — the same expression the CLI and `--params` take. |
| `zombi2.tree` | `from zombi2 import tree` (or `from zombi2.tree import Tree, read_newick, …`) — the shared dated-tree datatype and its toolkit. **`Tree`** (`nodes: {id: Node}`, `root`) is a lean dataclass; its **methods are structural self-queries only** — `leaves()`, `extant()`, `extinct()`, `unsampled()`, `to_newick()`. **`Node`** (`id, parent, birth_time, end_time, children, fate`). Everything that transforms or analyses a tree is a **free function**: `read_newick(newick, *, tip_fates=None)` → `(Tree, {id: label})`; the transforms `prune(tree, keep="extant")`, `with_stem(tree, length, *, mode="set")`, `make_ultrametric(tree, *, tol=1e-3)` (snaps rounding, raises on real tip-date signal), `rescale(tree, *, height=None, factor=None)`, `red_scaled(tree)` (all → `Tree`); the analyses `relative_evolutionary_divergence(tree)` → `{id: RED}` and `distance(a, b, *, metric="rf"|"rf-normalized"|"branch-score")` → `float`. (The toolkit grows by adding free functions, never by growing the class.) |
| `zombi2.species` | `simulate_species_tree(birth, death=0, *, n_extant=None, total_time=None, mass_extinctions=None, sampling=1.0, fossils=0.0, seed=None, max_lineages=100_000)` → `SpeciesResult(.complete_tree, .extant_tree, .fossils, .events, .seed)`. Also `Event` (the recorded true history). Trees are `zombi2.tree.Tree`. |
| `zombi2.genomes` | `simulate_genomes_unordered(tree, *, duplication=0, transfer=0, loss=0, origination=0, transfer_to="uniform", replacement=False, self_transfer=False, initial_families=100, families=None, family_speed=None, max_family_size=10.0, seed=None)` → `GenomesResult(.complete_tree, .genomes, .events, .seed, .family_counts())`. Also `GeneCopy(id, family)`, `Distance(decay=1.0)`. `transfer_to` takes `"uniform"` · `"distance"` / `Distance(decay=)` · `mod.DrivenBy(source, mapping)` — the **choice slot** of SPEC §5, where the mapping's numbers are per-candidate weights, not rate multipliers (unordered engine only). |
| `zombi2.genomes` (ordered) | `simulate_genomes_ordered(tree, *, duplication=0, transfer=0, loss=0, origination=0, inversion=0, transposition=0, translocation=0, chromosomes=1, topology="circular", fission=0, fusion=0, chromosome_origination=0, chromosome_loss=0, <event>_extension=None (→ Geometric(mean=1)), inversion_probability=0, transfer_to=…, replacement=…, self_transfer=…, initial_families=100, families=None, family_speed=None, max_family_size=None, seed=None)` → `OrderedGenomesResult(.complete_tree, .genomes, .events, .rearrangements, .chromosome_events, .seed, .family_counts(), .gene_order(), .event_positions)`. Every gene-level event acts on an **extension** (a run of consecutive genes, length ~ `<event>_extension`; origination is single). `topology` decides where a run stops: on a `"circular"` chromosome it **wraps** past position 0 and is capped only by the whole chromosome, on a `"linear"` one it stops at the last gene. Also `Gene(id, family, strand)`, `Chromosome(id, topology, genes)`, `Inversion` · `Transposition` · `Translocation` (identity-preserving, in `.rearrangements`), `EventPosition` (where each gene-genealogy `Event` fired — the positional companion to the position-blind log, in `.event_positions`), `ChromosomeEvent` (kinds: origination · speciation · fission · fusion · loss — the reticulating chromosome network's edge list). Shared spine (`Event`, live-set, transfer mechanics) lives in `genomes/{events,_live,_transfer}.py`; `ChromosomeEvent` and its `chromosome_events_tsv` writer live in `genomes/chromosomes.py`, one home for both the ordered and the nucleotide engine. |
| `zombi2.genomes` (nucleotide) | `simulate_genomes_nucleotide(tree, *, inversion=0, translocation=0, transposition=0, loss=0, duplication=0, transfer=0, origination=0, <event>_length=50.0, inversion_probability=0, fission=0, fusion=0, chromosome_origination=0, chromosome_loss=0, chromosomes=1, root_length=1000, topology="circular", genes=0, gene_length=100, gff=None, fasta=None, trim_overlaps=False, transfer_to=…, self_transfer=…, seed=None)` → `NucleotideGenomesResult(.complete_tree, .genomes, .events, .rearrangements, .chromosome_events, .seed, .gene_spans, .gene_names, .gene_strands, .mosaic(), .trace_back(), .ancestry(), .root_blocks, .gene_trees)`. The genome is a nucleotide sequence of `Block`s (a run of one unbroken ancestry) on `Chromosome`s; a **declared gene is indivisible** (an event that would cut one redraws). Its `Origination` · `Loss` · `Duplication` · `Transfer` · `Speciation` records are **positional already** (each names ancestral intervals), so there is no `EventPosition` companion here. Rates are **constants** — no `scope × modifiers` grammar yet. |

| `zombi2.tools` | `homology_table(root)` → `(labels, matrix)` and `homology_tsv(root)` classify one gene tree's leaves — the event at each pair's MRCA is a **speciation** → `O` (ortholog), **duplication** → `P` (paralog), or **transfer** → `X` (xenolog); `write_homology({family: GeneTree}, dir)` writes one `homology_fam<f>.tsv` per surviving family. Exact, not inferred: ZOMBI recorded the embedding (see `genomes/gene_trees.py`). On the CLI as `zombi2 tools format DIR` (`--format` defaults to `homology`). |

Every level returns a `<Level>Result` bundle sharing the spine `.events` / tree(s) / `.seed` /
`.write(dir, [...])`, with the `record=[...]` memory dial.

## The rebuild order

Level by level, each with its chapter written alongside:

1. **Species** ✅ — the forward birth–death engine. Chapter 4 written.
2. **Genomes** — unordered D/T/L/O ✅ (Ch5); **ordered** ✅ slices 1–3 — chromosomes as
   identity-bearing containers; **segmental** D/T/L/O + inversion/transposition/translocation (every
   gene-level event acts on an *extension* of consecutive genes, the ZOMBI1 model); and the
   number-changing tier (fission/fusion/origination/loss) forming a reticulating chromosome network,
   recorded as the `chromosome_events` edge-list ground truth.
   **Nucleotide** ✅ — blocks of unbroken ancestry, indivisible declared genes (from `--gff` or an
   even layout), the same event set in base pairs, its own outputs and `--resolution nucleotide`
   (unordered ⊂ ordered ⊂ nucleotide). Its rates are still constants: the
   `scope × modifiers` grammar 🔨 is not wired there.
3. **Sequences** ✅ — substitution models + both lineage clocks on the gene trees: `ByLineage` (uncorrelated / relaxed) and `FromParent` (autocorrelated, drifting parent→child down the species tree).
4. **Traits** ✅ — the continuous / discrete overlay models on the species tree.
5. **Coupling** — the one mechanism `mod.DrivenBy(source, mapping)` (SPEC §2–4). **Conditioned** ✅
   (source = a file: `rates/driver.py` + the target level runs it — e.g. genome loss driven by a trait;
   all four unordered D/T/L/O rates, plus the `transfer_to` **choice slot**, where the same modifier's
   numbers are per-candidate weights instead of rate multipliers — SPEC §5);
   **joint** ✅ (source = a live level: `joint/simulate_joint` grows both — a discrete trait drives
   speciation, BiSSE/MuSSE). There is **no `coupling` package**: conditioned folds into the target
   level, so the only engine is `joint`. Remaining: joint gene-content→speciation, conditioned
   trait→sequence-clock, continuous drivers (QuaSSE — needs thinning). "Coupling" stays the level's
   name in SPEC and the manual's Part III.

## Docs & manual

Same strategy, so the docs are **watched growing** instead of retrofitted:

- **Manual (the book)** — Chapters 1–9 written (`manual/book/ch1.md … ch9.md`) plus appendices A–C;
  the genome resolution ladder is Ch4–Ch6 (unordered ⊂ ordered ⊂ nucleotide).
- **Docs site** — each chapter publishes as a snippet include (Ch1–Ch9 under `docs/guide/`, the
  appendices under `docs/reference/`); the site shows only what the clean core supports, built with
  `mkdocs --strict`.

---
*When a convention needs to change, change `SPEC.md` (the words) or this map (the shape) **first**, then propagate.*
