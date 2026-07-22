# Sequence API — design target

**Status: the design to build.** The target for rewriting `zombi2/sequences`, the detailed consequence
of `SPEC.md` for the sequence level. Designed with Adrián on 2026-07-18. **Not built yet** — today's code
ships an eight-class `Clock` hierarchy (§ *What to delete*). Parallels `species-api.md` and
`genome-api.md`; read those first.

---

## The problem it fixes

`zombi2/sequences/clocks.py` ships **eight** clock classes — `Clock`, `StrictClock`,
`UncorrelatedLogNormalClock`, `UncorrelatedGammaClock`, `WhiteNoiseClock`,
`AutocorrelatedLogNormalClock`, `CIRClock`, `RateVariation` — a zoo exactly like the seven species
processes and the `RateModel` hierarchy. But a clock is not a special kind of object: it is a **modifier
on the substitution rate**, the same `scope(base) × modifiers` grammar as every other level (`SPEC §5`).
"That is the one place the word *clock* belongs" (Ch2). Once clocks are modifiers, the zoo collapses to
the same small shared vocabulary the rest of ZOMBI2 uses.

## The entry point

One function, taking the **genome run** whose gene trees the sequences evolve along:

```python
genomes_run = genomes.simulate_genomes_unordered(species_run, duplication=0.2, loss=0.2, seed=1)
sequences.simulate_sequences(genomes_run, model=hky85(kappa=2.0), length=1000, seed=1)
```

(Today the surface is a `SequenceEvolution` class plus `evolve_on_tree`; the function is the target.)

**The whole run, not just its `.gene_trees`** (decided with Adrián, 2026-07-21). A bare
`{family: GeneTree}` mapping would run, but silently degraded: the lineage clock below is drawn *per
species branch* and shared across families, so without the species tree there is no clock and no
`species_phylogram`. A level reads the level above it whole; a mapping is rejected with the loud
error the other levels use.

## Two things vary, and they are different axes

A substitution rate can vary **across branches** (the clock) and **across sites** (rate heterogeneity).
They are orthogonal and must not be conflated:

- **Across branches — the clock.** A modifier on the substitution rate, exactly like a modifier at any
  other level. This is where the zoo lived.
- **Across sites — +Γ.** Some sites evolve faster than others, drawn from a Gamma. This is the classic
  `+Γ` of phylogenetics and stays its own argument (`gamma=0.5`, the shape α), not a clock.

```python
sequences.simulate_sequences(genomes_run,
    model=gtr(...),                       # the substitution model (a menu, see below)
    substitution=1.0 * mod.ByLineage(spread=0.3),   # the clock: across-lineage variation
    gamma=0.5,                            # +Γ: across-site variation (shape α)
    length=1000, seed=1)
```

## The substitution model is a MENU, not a zoo

Not everything collapses, and that is correct. JC69 / HKY / GTR / LG are genuinely different chemistry —
different rate matrices — so they stay a **menu of constructors**, each taking its own physical
parameters:

```python
model=jc69()                       # no free parameters
model=hky85(kappa=2.0)             # transition/transversion bias
model=gtr(rates=..., freqs=...)    # six exchangeabilities + base frequencies
model=lg()                         # empirical amino-acid matrix
model=gy94(omega=0.2)              # codon model, dN/dS
```

Faking a grammar over the matrices would be worse than a menu. The menu is the honest shape here.

## The clock collapses to three modifiers + the strict default

Every clock in the literature is one of three things happening to the rate along the tree, expressed with
the **same modifiers the other levels use**:

```python
# strict clock — no across-branch variation (the default; write nothing)
substitution = 1.0

# uncorrelated / relaxed — each lineage draws its rate independently (i.i.d.)
substitution = 1.0 * mod.ByLineage(spread=0.3)                 # lognormal (default)
substitution = 1.0 * mod.ByLineage(spread=0.3, dist="gamma")   # gamma; white-noise is another dist

# autocorrelated — the rate drifts continuously along the tree (geometric Brownian)
substitution = 1.0 * mod.FromParent(spread=0.3)

# CIR — the same, but mean-reverting (Ornstein–Uhlenbeck)
substitution = 1.0 * mod.FromParent(spread=0.3, reverts_to=1.0)

# the Markov clock — the rate hops between discrete categories along the branches
substitution = 1.0 * mod.Markov(rates=[0.5, 1.0, 2.0], switch=0.1)
```

Three modifiers, grouped by what memory the rate has:

- **`ByLineage`** — *no memory*: each lineage independent. (The uncorrelated / relaxed family.)
- **`FromParent`** — *continuous memory*: the rate drifts, parent to child. (Autocorrelated; CIR is this
  with `reverts_to`.)
- **`Markov`** — *discrete memory*: the rate switches between a few states via a CTMC on rate categories.

Two of the three are **shared across levels**, which is the whole point of the grammar:

- `ByLineage` is the lineage-twin of the genome level's `ByFamily` — the same i.i.d.-heterogeneity idea,
  by lineage instead of by family.
- `FromParent` is **literally the species `FromParent`** (ClaDS): a rate that drifts along the tree. The
  autocorrelated molecular clock and ClaDS diversification are the same modifier at two levels.
- `Markov` is new to sequences, but even it echoes species: a clade shift is one discrete rate jump;
  `Markov` is that, happening repeatedly at a rate.

## Which tree the clock rides — two axes, two modifiers (decided)

A substitution rate can vary two ways across branches, and they are **different modifiers we already
have** (decided with Adrián, 2026-07-18):

- **The lineage clock rides the *species* tree.** A clock is a property of a *lineage* — a whole species
  runs hot or cold, and every gene passing through that branch feels it. So `ByLineage` / `FromParent` give
  **one clock value per species lineage**, shared by all its genes. Each gene-tree branch reads the clock
  of the species branch it is reconciled to — ZOMBI2 knows that reconciliation exactly, so it is automatic.
- **Per-family variation is `ByFamily`** — the *same modifier as the genome level*: some families evolve
  faster than others whatever the lineage, each a constant speed.

They **compose**, reproducing today's lineage-clock × per-family-speed (`R_b · s_g`):

```python
substitution = 1.0 * mod.ByLineage(spread=0.3) * mod.ByFamily(spread=0.5)
```

**Deferred:** a fully *per-gene-tree-branch* clock (a single family fluctuating independently branch by
branch *within* a lineage) — exotic, rarely what "a clock" means; that one would be **`ByBranch`** (per
gene-tree branch), a name reserved for it. `ByLineage` is the species-lineage clock: one value per
species lineage, shared by every gene passing through it.

## The literature → command bridge (goes in the chapter)

The deprecated model names survive **only** in this table — a reader who knows "I want a CIR clock" finds
the command, and no one has to memorise the acronyms anywhere else. Every chapter carries one of these;
this is the sequence chapter's.

| Literature name | What it does | ZOMBI2 |
|---|---|---|
| Strict / global clock | one rate everywhere | `substitution = 1.0` (default) |
| Uncorrelated lognormal (UCLN) | each lineage i.i.d. lognormal | `1.0 * mod.ByLineage(spread=…)` |
| Uncorrelated gamma (UGAM) | each lineage i.i.d. gamma | `1.0 * mod.ByLineage(spread=…, dist="gamma")` |
| White-noise clock | each lineage i.i.d., short memory | `1.0 * mod.ByLineage(spread=…, dist=…)` |
| Autocorrelated lognormal (Thorne–Kishino) | rate drifts along the tree | `1.0 * mod.FromParent(spread=…)` |
| CIR clock | drift with mean-reversion | `1.0 * mod.FromParent(spread=…, reverts_to=…)` |
| Discrete-category / random local clock | rate hops between categories | `1.0 * mod.Markov(rates=[…], switch=…)` |
| +Γ rate heterogeneity | variation across sites | `gamma=α` (not a clock) |

## Still to design

- **Decided:** `ByLineage(spread=, dist=)` exposes the distribution — `dist="lognormal"` (default) or
  `"gamma"`. No separate "white-noise" label (per-branch i.i.d. *is* white-noise).
- **Decided:** CIR is `FromParent(spread=, reverts_to=, pull=)` — mean-reversion is `reverts_to` (target) +
  `pull` (strength) on `FromParent`, the *same two knobs* as the OU trait. Plain `FromParent(spread=)` = pure
  drift (autocorrelated clock); add `reverts_to`+`pull` = CIR. One modifier across species/sequences/traits.
- `Markov`'s exact signature (`rates=` as multipliers vs absolute; `switch=` one global rate vs a matrix).
- Where a **trait- or driver-conditioned** clock lives — that is `DriverClock` today, and it belongs to
  Part III (a `traits → sequences` conditioning), not to this menu.
- Codon-model surface (`gy94`/`mg94`, the M-series site models) — already shipped; confirm they read as
  menu constructors here.

## What to delete / change in `zombi2/sequences`

- Delete the `Clock` hierarchy (`StrictClock`, `UncorrelatedLogNormalClock`, `UncorrelatedGammaClock`,
  `WhiteNoiseClock`, `AutocorrelatedLogNormalClock`, `CIRClock`, `RateVariation`). Clocks become the
  shared `mod.ByLineage` / `mod.FromParent` / `mod.Markov` modifiers on the `substitution` rate.
- Keep the substitution models as a **menu** of constructors (`jc69`, `k80`, `hky85`, `gtr`, `lg`, …);
  they are genuinely different matrices, not a zoo.
- Keep `+Γ` as its own `gamma=` argument (across-site variation ≠ across-branch clock).
- Add the `sequences.simulate_sequences(genomes_run, …)` entry point over the existing evolution core.
- `DriverClock` → a Part III conditioning (`traits`/level → sequences), not a menu member.

## Sequences on a nucleotide genome (decided 2026-07-22, not yet built)

Today the sequence level runs only on an unordered or ordered genome run. It refuses a
`NucleotideGenomesResult` — a type check, not a missing capability: that result already carries the
`complete_tree` and the recovered trees the level needs. Wiring it is what makes a **real gene layout
from a GFF** and **evolving residues** the same run, which is the one combination the two paths
cannot express apart.

**What it evolves.** Every recovered **root block**, not just the declared genes — genes *and*
intergenic spacer — so the run reconstructs the whole genome rather than a handful of loci.

**The intergenic trees are built** — `NucleotideGenomesResult.block_trees`, one per recovered root
block, spacer as well as genes. It was a small change, not the large one first scoped here: a block
never splits, so its size is fixed and its whole genealogy is already in the event log exactly as a
gene's is, and `_recover_gene_trees` differed between the genic and uniform cases *only* in which
blocks it was pointed at. A gene's tree comes back with the same topology and branch lengths either
way; the `g<id>` labels differ, because segment ids are handed out as the recovery walks its targets.

**Two models, keyed on what is already recorded.** `Block.gene` is the family id for a gene and `0`
for spacer, so no new classification is needed. Genes take `model=`; spacer takes `intergene_model=`,
defaulting to `jc69` — flat, equal frequencies, no free parameters, which is the right null for
unconstrained sequence. Spacer also runs faster: `intergene_speed=3.0` multiplies the gene rate, so
the ratio (relaxed constraint) stays fixed when the overall rate changes.

**`length=` is rejected here, not ignored.** Each block brings its own length in bp from the genome,
so a single `length=` would contradict the coordinates the genomes run recorded — the same discipline
that rejects `--kappa` for `jc69`.

**Output.** The per-block alignments, plus one FASTA per extant lineage with the blocks concatenated
in genome order — `SequencesResult.genomes`, written as `genome_<lineage>.fasta`, one record per
chromosome.

**The assembly is split at the level boundary.** `NucleotideGenomesResult.assembly(node)` gives, per
chromosome, the pieces in physical order as `(block, gene, start, end, strand)`; the sequence level
slices, flips the `-1`s and concatenates. Nothing in the genomes level knows about letters, and the
layout is testable on its own — expanded back to one entry per nucleotide it must equal
`trace_back(node)`, which is how it is pinned, at every node rather than only the leaves.

Two facts the first sketch had wrong, and only the **extant leaves** vote on where the partition is
cut, which is what settles both. A working block is **not** a slice of one root block: at a leaf, that
leaf's own breakpoints are all in the partition, so the partition is at least as fine as its blocks and
a block spans one or more of them, cut into a piece each — every piece a *whole* block, and on a `-1`
block they come out in descending coordinate order. The sub-block case is real only at an **ancestor**,
and only where the ancestor has a breakpoint no survivor has: the descendants that inherited it died out
or lost the material, while another lineage kept that stretch unbroken. (An earlier note here blamed
transfer. It is not transfer — sub-block pieces appear just the same with transfer switched off; plain
extinction and loss are enough.)

The copy→gene-id join is what makes it work: `_emit_block_events` mints a fresh segment id at every
event, so a copy that duplicated twice is three genes in a row, and the one a genome still carries is
the **last** rung of that ladder. That map (`{(block, copy): gene id}`) now comes back from the
recovery beside the trees.

**Every node, not just the survivors.** `.genomes` covers the whole tree in one map — extant tips,
ancestors and extinct lineages — the same coverage the genome level's own `.genomes` has, written as
one `genome_<lineage>.fasta` apiece. (It briefly had a second `.ancestral_genomes` map and an
`ancestral` filename prefix; both went, because no node is a special case and an extinct tip is not an
ancestor.) Two changes made it possible, and both were forced by a stress run (91 nodes, three replicons, every event kind, 128 losses) that the small
runs could not produce:

1. **The root partition is cut at every node's breakpoints, not the extant leaves'.** Cutting from the
   survivors alone leaves two holes: material no survivor kept has no block at all, and an ancestor can
   hold a *fragment* of a block whose genealogy — being the survivors' — has no lineage for it (that one
   surfaced as a `KeyError` deep in the assembly). Counting every node closes both, because a node's own
   breakpoints are then all in the partition, so each of its blocks is a whole number of root blocks.
   Reading each node's **final** genome is enough: a breakpoint matters only where material survives on
   one side and not the other, and the surviving side carries it to the end of its branch.
2. **`_split` puts the dead tips in `.ancestral`.** A copy ended by a loss, or one whose species went
   extinct, is a node of the tree with a sequence at it. Leaving them out gave a phylogram whose tips
   named sequences that existed nowhere, and made an extinct lineage's genome unreconstructable.

Consequences worth knowing. `assembly` returns `(block, gene, strand)` — a piece is now *always* a whole
block, so the slice fields went away with the case that needed them. A gene surviving only in lineages
that died now gets a (complete-only) tree, where before it got none. The partition costs 1.34× the
blocks on that stress run, and recovers the 2 738 bp no survivor kept. All 91 genomes rebuild exactly,
3.26 Mb checked base by base.

Care needed: `_recover` is pinned by the Swenson port (`tests/test_nucleotide_model_krister.py`),
whose whole argument is that the *written files* must mean what they claim.
