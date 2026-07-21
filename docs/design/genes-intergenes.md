# Genes & intergenes — design note

**Status: design to revise.** The **genic layer** on top of the built nucleotide engine
(`zombi2/genomes/nucleotide.py`, merged in #196). Designed with Adrián on 2026-07-21 (discussion). **Not
built yet.** Read `genome-api.md` (the layering) and `SPEC.md §4` first. This note is the model; the
chapter is its exposition (concepts → code → chapter).

---

## What it adds

Today a nucleotide genome is uniform sequence — every block is anonymous ancestry. The genic layer
**declares** some blocks to be **genes** (the coherent evolutionary units the gene trees are actually
about) and leaves the rest as **intergenes**, the spacer. Nothing else about the engine changes: genes
and intergenes are both just blocks, and the recovery is the same recovery.

## Vocabulary (extends `SPEC.md §6`)

- **block** — *(unchanged)* a maximal run of one unbroken ancestry; the persistent unit of the
  representation and the ancestry.
- **gene** — a *declared, indivisible* block: **one** block, one **family**, one id, **never split**.
  It carries a gene tree. Atomic under every event.
- **intergene** — the free-fragmenting spacer between genes: it splits at any breakpoint into as many
  blocks as events dictate. Its blocks carry ancestry like any other (we track it), but we do **not**
  surface a per-intergene "gene tree".
- **segment** — *(unchanged)* the transient extent an event acts on.

A genome is therefore an **alternating chain** `I₀ G₁ I₁ G₂ … Gₙ Iₙ` — genes atomic, intergenes
fragmenting. With **no genes declared** it is one big intergene, and the engine is exactly today's
uniform-sequence model (the degenerate case). So the genic layer is: a **classification** of blocks +
one **hard rule** (a gene never splits) + (later) which substitution model applies downstream.

---

## The two-track event model (the core)

Two complementary ways an event finds its target. They are chosen so that **per-family rates and
multi-gene *nucleotide* extents never coexist** — which is what makes the whole thing well-defined.

| | **Track 1 — content** | **Track 2 — structure** |
|---|---|---|
| **anchors** | on a **gene** | in an **intergene** |
| **rate scales with** | family rates (`Σ rₘ`) — *per copy* | intergenic length — *per intergenic nucleotide* |
| **extent measured in** | **genes** (a run of adjacent genes) | **nucleotides** |
| **breakpoints** | in the flanking intergenes, by construction | drawn; **rejected** if inside a gene |
| **on a gene** | size-**blind** (targets the gene directly) | size-**dependent** (must engulf it whole) |
| **events** | `duplication`, `loss`, `transfer`, `origination` | `inversion`, `transposition`, `translocation`; `deletion`, `segmental_duplication`; `insertion` |

Both tracks keep genes **whole** — Track 1 because it is defined on genes, Track 2 because it rejects any
breakpoint that would fall inside one. So a gene is always atomic: wholly affected or wholly untouched.

### Track 2 — the breakpoint rule (rejection, not snapping)

A Track-2 event nucleates at a nucleotide inside an intergene and extends by a drawn length (in both
directions). If **either breakpoint lands strictly inside a gene, the draw fails and is retried** (up to
`N` times, default ~10; then the event is a no-op). *Engulfing* a whole gene is fine — both breakpoints
then sit in the flanking intergenes.

Why reject and not snap the breakpoint out to the gene edge: **length honesty.** Snapping would inflate a
drawn extent to "whatever it takes to clear the gene". With a 30-nt mean extent and a genome
`intergene(10) · gene(1 000 000) · intergene(20)`, a snap turns a 30 into a 1 000 000+ event. Rejection
keeps the contract exact: **a successful Track-2 event's span equals its drawn length — or it never
happened.** A consequence, and a desirable one: to reach *past* a big gene you must draw a big extent, so
**large genes are rarely swept up** (they are structurally stable); small ones are swept more easily.

Because the anchor is length-weighted over *intergenic* nucleotides, the rate `= rate × (total
intergenic length)`, so it auto-scales with available room — a gene-dense genome rearranges slowly, a
spacer-rich one freely — and the anchor itself never fails; only the extent can.

### Track 1 — nucleation owns the rate (per-family, with dragging)

A Track-1 event nucleates on **one gene**, then may **drag its neighbours** along. Model it as a single
Poisson process; the rate of "nucleate on copy `g`, extent `m` genes" factorises:

```
rate(nucleate at g, extent m) = r_g · P(extent = m)      r_g = g's family rate; P is family-blind
Σ_m …                         = r_g                       (the extent distribution integrates to 1)
total Gillespie rate          = Σ_g r_g
```

So the Gillespie needs nothing new:

1. total rate `= Σ_g r_g`; draw the waiting time from it;
2. pick the **anchor** gene with probability `r_g / Σ r_g`;
3. draw the extent `m` from a **family-blind** distribution → a run of `m` consecutive genes;
4. the two breakpoints fall at drawn positions in the two **flanking intergenes** (which fragment; a
   fraction of each rides along). Apply (copy / delete / move / send) the whole run.

The genes the extent *reaches* never enter the rate — they are passengers, dragged because they are
physically adjacent, not because of their own rate. So `r_family` is a **nucleation** rate, and a gene's
**realised** turnover is `r_g + (chance a hot neighbour drags it) + (chance a Track-2 event engulfs it)`.

### Two pathways for copy number — on purpose

A gene's copy number changes **both** ways:

- **Track 1** — a targeted, size-blind `duplication`/`loss`/`transfer` nucleated on the gene;
- **Track 2** — a size-dependent `segmental_duplication`/`deletion` that **engulfs** it whole.

Both are real biology (gene-level turnover *and* large structural mutations), and both are wanted. The
gene tree cannot tell which mechanism duplicated a gene — a duplication is a duplication — so the
recovery is agnostic (below).

### The default: uniform per copy

Family heterogeneity (`ByFamily` / `family_speed`, `genome-api.md`) stays **parked**. By default every
copy shares one rate, so `Σ_g r_g = rate × (gene count)` and the anchor pick collapses to **uniform**.
The nucleation-weighted pick is simply the forward hook: turning heterogeneity on later swaps the uniform
pick for a weighted one and changes nothing else.

---

## Seeding

The root genome is **declared**, two accepted forms:

- a **GFF** — exact gene coordinates (and names) on the seed replicon(s); or
- a **genes + intergene-distribution file** — a gene count (and per-gene lengths, or a length
  distribution) plus an intergene-length distribution, laid down as the alternating chain.

Each seeded gene becomes a **family** with an id (the handle for its gene tree and, later, a named
family — cf. `family_names` in the unordered result).

---

## Genealogy & recovery — almost nothing changes

A gene is a block that **never splits**, so:

- **each gene → exactly one persistent block-lineage → exactly one gene tree.** This is the whole point,
  and it drops out of the existing recovery for free — a gene is just a root-block that is never cut.
- **intergenes fragment** into many blocks whose ancestry the copy-lineage log already records. We
  **track their genealogy** (they are ordinary blocks) but do **not** surface per-intergene gene trees.

So the recovery is the recovery we built (root partition → per-block replay → cross-check). The only
additions:

1. **a genic tag** on blocks (a gene id / family) — set at seeding, inherited through every event;
2. **`_split_at` refuses to cut a genic block** — this is exactly what *raises* the Track-2 rejection and
   guarantees genes stay atomic (an engulfed gene is wholly in or out, so the cross-check still holds);
3. **`result.gene_trees` filters to genic families** — one tree per gene; intergene block-genealogies
   stay available in the log but are not built into trees;
4. **fission / fusion cut-points are intergenic too** — a chromosome fission must cut in an intergene
   (same reject-in-gene rule), so a gene is never split by the tier either.

---

## Naming & the one migration this implies

The layering contract (`genome-api.md`: *unordered ⊂ ordered ⊂ nucleotide*, shared params mean the same
thing) fixes the names:

- **`duplication` / `loss` / `transfer` / `origination`** are the **Track-1, per-family gene events** —
  the same meaning as in ordered/unordered. (With no genes declared they have nothing to act on, so they
  are inert; a gene-free run uses the Track-2 events.)
- The nucleotide-specific **Track-2** events are their own names: `inversion` / `transposition` /
  `translocation` (already built) + **`deletion`** and **`segmental_duplication`** (the size-dependent,
  gene-engulfing content changers) + **`insertion`** (intergenic growth / indels).

**Migration:** today's nucleotide `duplication` / `loss` are per-nucleotide segmental — i.e. they are the
Track-2 events. Aligning to the layering renames them **`segmental_duplication` / `deletion`**, and frees
`duplication` / `loss` for the Track-1 gene meaning. This is the one breaking rename; it is worth it for
the cross-level consistency the layering promises. *(Open — see below.)*

---

## Deferred

- **Pseudogenization** — an explicit `gene → intergene` event (drops indivisibility, and flips the block
  to neutral evolution downstream). Not now; later.
- **Family heterogeneity** (`ByFamily` / `family_speed`) — parked; the nucleation pick is the hook.
- **Replacement transfer** — additive only for now (as in the built engine).
- **Segmental `transfer`** (a Track-2 transfer that engulfs genes into a recipient) — start with Track-1
  `transfer` only; add later if wanted.
- **Substitution-model routing** — the gene/intergene tag will pick coding vs neutral models at the
  sequence level; out of scope here.

---

## Open — to settle in revision

1. **The rename.** `duplication`/`loss` → `segmental_duplication`/`deletion` for the Track-2 events,
   `duplication`/`loss` re-assigned to Track-1 gene events (recommended, per layering) — or keep the
   current names on the per-nt events and call the gene events `gene_*` (breaks the layering)?
2. **Track-1 extent shape** — run direction (one-sided from the nucleated gene vs bidirectional) and the
   distribution (Geometric in genes, mean `*_extension`, mirroring ordered).
3. **`insertion` / intergenic length** — is de-novo intergenic DNA its own event, or folded into
   `origination` (which currently births a whole new source)? And is `origination` now specifically a new
   **gene** (+ flanking intergene)?
4. **Retry cap `N`** — default and whether it is user-exposed.
5. **GFF/declaration file schema** — column set and where multi-replicon karyotypes are expressed.
