# ZOMBI2 — Code Map

**Status: AUTHORITATIVE. The code companion to [`SPEC.md`](SPEC.md).**

`SPEC.md` fixes the **model and the words**. This fixes the **files and the names**: every module, every
public name, and its one canonical home. When code, a docstring, the CLI, or a chapter names something
differently, this document wins and that other file is a **fossil**.

Read `SPEC.md` first (the model), then this (the shape). Two principles throughout:

- **One canonical path per name** — no top-level re-exports, no aliases, so nothing drifts. Reach
  everything through its level package (`from zombi2 import species`, `from zombi2.rates import scope`).
- **Concepts → code → chapter**, and this map is the **as-built truth**: where an older design doc
  disagrees (e.g. `species-api.md` still says `n_tips` / `age`), this map is right — the built engine uses
  `n_extant` / `total_time`, and time is always **forward from the crown** (SPEC §5).

---

## The clean core and the quarantine

The rewrite is a **clean core grown from `SPEC.md`**, not a migration of the old codebase. Everything not
yet rebuilt is **quarantined** in `legacy/` at the repo root — read-only reference, **not importable**, not
wired to anything. Features are ported out of `legacy/` **deliberately, one level at a time, renamed to
this map** — never carried, never re-imported. That is what keeps the active tree small enough to hold in
one head, and in one context window.

Legend:  ✅ built · 🔨 to build · 📦 quarantined in `legacy/`

## The tree

```
zombi2/
  __init__.py          thin — no top-level re-exports (one canonical path per name)
  rates/               the cross-level rate grammar (SPEC §5): effective rate = scope(base) × modifiers
    scope.py        ✅ PerCopy · PerLineage · PerSite · Global
    modifiers.py    ✅ OnTime · OnTotalDiversity · FromParent · ByLineage   (ByFamily · Markov 🔨)
    rate.py         ✅ Rate · as_rate       (internal plumbing; users never build a Rate directly)
    distributions.py✅ Fixed · Exponential · Gamma · LogNormal · Uniform · Geometric   (value / length distributions)
  species/           ✅ simulate_species_tree → SpeciesResult ;  Tree · Node  (the shared dated tree)
  genomes/           ✅ simulate_genomes_unordered · simulate_genomes_ordered   (ordered = chromosomes + segmental D/T/L/O + inversion/transposition/translocation + tier ✅; nucleotide 🔨)
  sequences/         🔨 simulate_sequences → SequencesResult
  traits/            🔨 simulate_traits → TraitsResult
  coupling/          🔨 conditioned · joint       (what the old "coevolve" becomes; SPEC §2–4)
  cli/               ~  species · genomes (clean); every other subcommand 📦
legacy/              📦 repo root, not importable — the old code, kept only to port from
```

The cross-level primitives live in **`zombi2.rates`** (scope, modifiers, rate, distributions) because the
rate grammar is one thing shared by all four levels. The **`Tree`** belongs to **`zombi2.species`**: the
species level produces it and the other levels read it (`from zombi2.species import Tree`).

## The names (as built)

| Canonical home | Public names |
|---|---|
| `zombi2.rates` | `from zombi2.rates import scope, modifiers` → `scope.Global`, `modifiers.OnTime({...})`. Scopes: `PerCopy · PerLineage · PerSite · Global`. Modifiers: `OnTime · OnTotalDiversity · FromParent · ByLineage`. |
| `zombi2.species` | `simulate_species_tree(birth, death=0, *, n_extant=None, total_time=None, mass_extinctions=None, sampling=1.0, fossils=0.0, seed=None)` → `SpeciesResult(.complete_tree, .extant_tree, .fossils, .events, .seed)`. Also `Tree`, `Node`, `prune(tree, keep="extant")`. |
| `zombi2.genomes` | `simulate_genomes_unordered(tree, *, duplication=0, transfer=0, loss=0, origination=0, transfer_to="uniform", replacement=False, self_transfer=False, initial_families=0, seed=None)` → `GenomesResult(.complete_tree, .genomes, .events, .seed, .family_counts())`. Also `GeneCopy(id, family)`, `Distance(decay=1.0)`. |
| `zombi2.genomes` (ordered) | `simulate_genomes_ordered(tree, *, duplication=0, transfer=0, loss=0, origination=0, inversion=0, transposition=0, translocation=0, chromosomes=1, topology="circular", fission=0, fusion=0, chromosome_origination=0, chromosome_loss=0, <event>_extension=Geometric(mean=1), inversion_probability=0, transfer_to=…, replacement=…, self_transfer=…, initial_families=0, seed=None)` → `OrderedGenomesResult(.complete_tree, .genomes, .events, .rearrangements, .chromosome_events, .seed, .family_counts(), .gene_order())`. Every gene-level event acts on an **extension** (a run of consecutive genes, length ~ `<event>_extension`; origination is single). Also `Gene(id, family, strand)`, `Chromosome(id, topology, genes)`, `Inversion` · `Transposition` · `Translocation` (identity-preserving, in `.rearrangements`), `ChromosomeEvent` (kinds: origination · speciation · fission · fusion · loss — the reticulating chromosome network's edge list). Shared spine (`Event`, live-set, transfer mechanics) lives in `genomes/{events,_live,_transfer}.py`. |

Every level returns a `<Level>Result` bundle sharing the spine `.events` / tree(s) / `.seed` /
`.write(dir, [...])`, with the `record=[...]` memory dial (see [`result-api.md`](result-api.md)).

## The rebuild order

Level by level, each with its chapter written alongside:

1. **Species** ✅ — the forward birth–death engine. Chapter 4 written.
2. **Genomes** — unordered D/T/L/O ✅ (Ch5 drafted); **ordered** ✅ slices 1–3 — chromosomes as
   identity-bearing containers; **segmental** D/T/L/O + inversion/transposition/translocation (every
   gene-level event acts on an *extension* of consecutive genes, the ZOMBI1 model); and the
   number-changing tier (fission/fusion/origination/loss) forming a reticulating chromosome network,
   recorded as the `chromosome_events` edge-list ground truth (`chromosome-network.md`). Only
   **nucleotide** (genes, indels) 🔨 remains (`genome-api.md`: unordered ⊂ ordered ⊂ nucleotide).
3. **Sequences** 🔨 — substitution + clocks on the gene trees.
4. **Traits** 🔨 — the overlay models on the species tree.
5. **Coupling** 🔨 — conditioned and joint (SPEC §2–4).

## The move (species, first)

The first act of the clean core: `zombi2/species_tree.py` becomes the `zombi2/species/` package (the
engine + `Tree`); the old model-zoo `zombi2/species/` (`sim`, `model`, `forward`, `ghosts`, `_caps`) goes to
`legacy/`; `zombi2/genomes_unordered.py` becomes `zombi2/genomes/`; and everything downstream still built on
the old engines (coevolve, tools, most of the CLI, the old sequence/trait code) goes to `legacy/` until its
level is rebuilt. The public import becomes `from zombi2 import species; species.simulate_species_tree(...)`
— real at last.

## Docs & manual

Same strategy, so the docs are **watched growing** instead of retrofitted:

- **Manual (the book)** — the 11-chapter index stands (SPEC §9). Ch4 (Species) written, Ch5 (Genomes I)
  drafted; the rest follow their level.
- **Docs site** — rewrite **Home** and the **Guide → Species Tree** pages fresh; **park "Get started"** and
  every other page in `legacy/`. The site shows only what the clean core supports, and grows a page as each
  level lands.

---
*When a convention needs to change, change `SPEC.md` (the words) or this map (the shape) **first**, then propagate.*
