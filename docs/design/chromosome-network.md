# Chromosome-network API — design target

**Status: the design to build.** The target for the chromosome tier of `zombi2/genomes`: recover the
**chromosome network** — chromosomes as identity-bearing entities with a genealogy, not just events
that reshape a genome. Designed with Adrián on 2026-07-18. The tier events already exist
(`FISSION`, `FUSION`, `CHROMOSOME_ORIGINATION`, `CHROMOSOME_LOSS`, `TRANSLOCATION`, plus within-chromosome
`INVERSION` / `TRANSPOSITION`); what is missing is the **genealogy across the whole run** and a
**network output**. Parallels `species-api.md` and `genome-api.md`; read those first.

Principle: **concepts → code → chapter.** Build the network, then the chapter documents it.

Realism is **not** the goal. Events are abstract — split or merge a run of elements. `circular` is
just a label (topology gates which fusions/fissions are legal, nothing more). No origins, centromeres,
telomeres, or biological legality constraints.

---

## The problem it fixes

The tier events fire and are logged, but the genealogy they imply is **only reconstructable within a
single species branch**. Two gaps:

1. **Identity dies at speciation.** `OrderedGenome.clone_reminting` (`genome.py:801`) mints a *fresh*
   `chrom_id` for every daughter chromosome (`new_chromosome()` in `__init__`, `genome.py:606`). The
   parent→daughter correspondence exists transiently — the `zip(self.chromosomes, new.chromosomes)` at
   `genome.py:805` walks them in lockstep — but only the **gene** mapping is returned (`genome.py:810`);
   the chromosome mapping is thrown away. So a chromosome in a leaf cannot be traced back past its own
   species branch. The docstring even says chromosomes are "re-minted at speciation (like gene copies)"
   (`genome.py:58`) — but unlike gene copies, **no `SPECIATION` edge is recorded for them.**

2. **The output is a flat event list, not a network.** `karyotype_trace.tsv` (`simulation.py:481`,
   `cli/genomes.py:697`) writes one row per `ChromosomeEvent` with `parents` / `children` chrom_ids —
   the raw edges of fission/fusion/origination/loss *within branches* — but nothing assembles them into
   a connected genealogy, and the speciation edges are absent, so the rows do not form one graph.

The data model is ready for this. `ChromosomeEvent` already carries `parents: tuple[int,...]` and
`children: tuple[int,...]` and its docstring states its purpose outright: "so a **karyotype genealogy**
can be reconstructed" (`events.py:210`). We are cashing that in.

---

## 1. Chromosome identity — the chromosome lineage

A **chromosome lineage** is the analogue of a gene lineage: a maximal run of the *same* chromosome
identity between two tier events. It is born by one event, ends at the next, and re-mints on every
event that touches it — exactly the discipline gene lineages already follow (`reconciliation.py:1`,
"every event re-mints lineage ids ... the per-family event log is a complete genealogy").

Concretely, a `chrom_id` is **stable between events** and **fresh after** any event that reshapes,
splits, merges, or copies it:

| Event on a chromosome | New `chrom_id`(s)? |
|---|---|
| nothing (flows down a branch) | same id, longer lineage |
| **speciation** (genome cloned) | each daughter chromosome re-mints → 2 children per parent |
| **fission** | source keeps its id + 1 new id (today: `genome.py:833`) — see §2 |
| **fusion** | kept id survives, absorbed id ends (today: `genome.py:540`) |
| **origination** | 1 new id, no parent |
| **loss** | id ends, no child |
| inversion / transposition / translocation | **same id** — identity persists (§2) |

The lineage is what makes the genealogy a *network*: a node is one chromosome lineage; an edge joins a
lineage to the lineage(s) it descends from. Because a chromosome carries an ordered set of *gene*
lineages, the chromosome network is a genuine middle tier — coarser than gene trees, finer than the
species tree.

**The one required code change for identity:** record the speciation correspondence. `clone_reminting`
already computes it (the lockstep `zip`); it must emit a `ChromosomeEvent(SPECIATION, branch, t,
parents=(pchrom.chrom_id,), children=(d1_cid, d2_cid))` — one per parent chromosome, mirroring the gene
`SPECIATION` record. Without this the network has no vertical (time-descent) edges across species and is
merely a per-branch sketch.

---

## 2. Each event as a network operation

Classify every tier event by its shape in the network. Three primitives: **bifurcation** (one parent,
two children — a tree split), **reticulation** (two parents, one child — a merge; what makes this a
network and not a tree), **birth/death** (a root / a leaf).

| Event | Parents → children | Network shape | Today |
|---|---|---|---|
| **speciation** | 1 → 2 | **bifurcation** (both daughter species inherit the chromosome) | edge **not** recorded (gap §1) |
| **fission** | 1 → 2 | **bifurcation** (split a run into two lineages) | `genome.py:824`, logged `parents=(src,) children=(src,new)` |
| **fusion** | **2 → 1** | **reticulation** (two lineages join; the crux) | `genome.py:521`, logged `parents=(keep,absorbed) children=(keep,)` |
| **origination** | 0 → 1 | **birth / root** (a de-novo replicon; a plasmid) | `genome.py:503`, logged `children=(new,)` |
| **chromosome loss** | 1 → 0 | **death / leaf** (the lineage ends; its genes die as LOSS) | `genome.py:511`, logged `parents=(lost,)` |
| **translocation** | 1 → 1 (×2 chroms) | **identity persists on both** (genes move, replicons stay) | `genome.py:745`, `nucleotide_genome.py:890` |
| inversion / transposition | 1 → 1 | **identity persists** (intra-chromosome reshuffle) | logged in the *gene* stream, no `ChromosomeEvent` |

Notes that need to be exact:

- **Fission is a bifurcation, but the source id is reused.** Today the source chromosome keeps its
  `chrom_id` and only the excised arc gets a new one (`children=(src, new_cid)`, `genome.py:838`). For a
  clean genealogy the *surviving* segment is a **new lineage** too — the pre-fission chromosome ended.
  **Decided (Adrián, 2026-07-18): at fission, re-mint both children** (`parents=(src,), children=(a, b)`
  with two fresh ids), so a fission node looks like a speciation node and no id spans an event. The
  byte-identity break on multi-chromosome seeds is **accepted** — the redesign does not preserve old
  seeds. *(Alternative reuse-and-special-case rejected — it makes the source lineage's "birth" ambiguous.)*

- **Fusion is the reticulation and the whole reason this is a network.** Two parent lineages
  (`keep`, `absorbed`) end; one child lineage begins. Today the child reuses `keep`'s id
  (`children=(keep,)`, `genome.py:547`) — same objection as fission. **Decided: the same rule** — re-mint
  the child (`parents=(keep, absorbed), children=(fused,)`). A reticulation node has **in-degree 2** — this
  is the node that breaks tree tooling (§6).

- **Translocation does not touch chromosome identity.** Genes move from chromosome A to chromosome B;
  both replicons persist with their ids (`genome.py:745`, "gids are untouched -> lineage-neutral"). In
  the *chromosome* network it is a **no-op edge**. But it is not invisible: it is the tier's analogue of
  transfer, and it means **a gene lineage can cross from one chromosome lineage to another** (§4). So
  translocation is recorded in the gene stream (as today, role `"translocated"`) *and* should annotate
  which chromosome lineages it bridged, so the gene-→-chromosome path stays recoverable.

- **Inversion / transposition** are pure intra-chromosome reshuffles: identity persists, no network
  edge, no gene-lineage change. They belong to gene order, not the chromosome network.

---

## 3. Output format — eNewick + the event table (both)

The species tree is one Newick string; a gene family is one Newick string; the chromosome genealogy is
**one network**. A network is not a tree, so plain Newick cannot hold it. Two artefacts, because they
answer different questions:

**`karyotype_network.enewick` — the topology.** Extended Newick (eNewick) is the standard serialisation
for phylogenetic networks with reticulations: a reticulation node is written once in full and then
referenced by a shared label `#H<k>` everywhere else it appears, so a two-parent node is expressed in a
one-parent-per-occurrence string. A **fusion** is exactly this — the fused child hangs under *both*
parents, tagged `#H1`:

```
((chrA#H1)chrP1, (chrB, chrA#H1)chrP2) ... ;
      └ fusion of chrA and chrB: the child #H1 appears under both parents
```

Fission and speciation are ordinary bifurcations; origination is a root; loss is a leaf. eNewick is
parseable by existing network tools (PhyloNetworks, Dendroscope, ape/evobiR), which is the point of
choosing it over a bespoke format.

**`karyotype_events.tsv` — the audited edge list.** Keep today's table (`simulation.py:485`), extended:

```
time    event   branch    parents        children      genome
3.14    FI      n7        c12            c40;c41        g_root      # fission (both re-minted)
4.02    FU      n7        c40;c9         c55            g_root      # fusion  → reticulation
5.00    S       n7        c55            c88;c89        g_root      # speciation (NEW: now recorded)
```

`event` uses the existing `EventType` chars (`FI` / `FU` / `CO` / `CL`, plus `S`). `parents` /
`children` are `;`-joined chrom_ids (already the format). The table is the **ground truth**; the
eNewick is derived from it. Emit the eNewick only when the network is non-trivial (any fission / fusion
/ origination / loss / >1 seed chromosome), exactly the condition that already triggers the karyotype
outputs (`cli/genomes.py:114`).

**Reconciling into the species tree.** Every `ChromosomeEvent` already records its `branch` — the
species-tree node it fired on (`events.py:216`). That *is* the reconciliation: each chromosome node is
stamped with a species branch, precisely as gene reconciliation stamps each gene node with a species
branch (`reconcile`, `reconciliation.py:402`; it is exact annotation, never LCA inference, because the
simulator knows the truth). So the chromosome network reconciled against the species tree is the same
eNewick with each internal node labelled `#H../branch` — no new machinery, just carry `branch` through.
What "reconciling a *network*" means when a fusion joins two branches is the open question in §6.

---

## 4. Relation to gene trees — species tree ⊃ chromosome network ⊃ gene trees

The three levels nest. A gene lives *on* a chromosome; a chromosome lives *in* a species. So a gene's
recorded history should let you recover **which chromosome lineage carried it, when**.

Today the hook already exists but is under-used: every `EventRecord` carries a `region: Region | None`
(`events.py:202`) and `Region.chromosome` is the `chrom_id` the event landed on (`events.py:95`,
`genome.py:696`). So the log *already* stamps each gene event with a chromosome id. What is missing is
that (a) the id is not re-minted consistently (§1), so it is not a stable lineage key across the run,
and (b) there is no reader that assembles a gene's chromosome path.

**A gene's chromosome path** is then reconstructable, and should be a first-class output:

- At each gene event, read `region.chromosome` → the chromosome lineage the gene sat on at that time.
- A **translocation** (§2) is the only gene-level event that *changes* a gene's chromosome lineage
  without a duplication/transfer — record both the source and dest chromosome ids on that gene op so the
  hand-off is explicit (`Region.dest` exists for this, `events.py:99`).
- A **transfer** moves a gene between *genomes*; the recipient's `region.chromosome` names which
  chromosome lineage received it.

Output: extend `gene_order.tsv` (already `species / chromosome / position / family / gid`,
`simulation.py:469`) with, per gene, the **chromosome lineage id** rather than the leaf-local
`chrom_id`, so a gene tree's tips can be coloured by chromosome lineage and a gene's walk across
chromosome lineages is queryable. The chromosome network is thus the **connective tissue**: gene trees
say *which genes are related*; the species tree says *which organisms*; the chromosome network says
*which replicon carried each gene at each moment*, and it is the only one of the three that reticulates.

---

## 5. API surface

Following `genome-api.md`: no rate/model objects, keyword rates in the `count(base) × modifiers`
grammar, chromosomes are ordered/nucleotide-only so they live on `simulate_ordered` /
`simulate_nucleotide` (never `simulate_unordered` — a multiset has no chromosomes).

**What the user sets** — the initial karyotype and the four tier rates:

```python
genomes.simulate_ordered(
    tree,
    duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,   # gene tier (genome-api.md)
    # --- initial karyotype ---
    chromosomes=8,                 # number of chromosomes seeded at the root
    topology="linear",             # "circular" (default, bacteria) | "linear" | a list per chromosome
    # --- chromosome tier: count(base) × modifiers ---
    fission=0.02,                  # per chromosome (a chromosome splits in two)
    fusion=0.02,                   # per chromosome (two chromosomes merge — the reticulation)
    translocation=0.05,            # per gene copy (a gene moves to another chromosome; a rearrangement, like transposition)
    chromosome_origination=0.01,   # per genome (a de-novo replicon / plasmid appears)
    chromosome_loss=0.01,          # per chromosome (a whole chromosome and its genes die)
    seed=1,
)
```

- **`chromosomes=`** replaces today's `n_chromosomes=` (`cli/genomes.py:174`); **`topology=`** replaces
  the `--linear-chromosomes` boolean (`cli/genomes.py:180`) and accepts a per-chromosome list
  (`["circular", "linear"]`) for a mixed karyotype, which the engine already supports via
  `circular: bool | Sequence[bool]` (`genome.py:583`).
- **The tier rates follow the same grammar as every other rate.** Defaults answer *per what*: fission /
  fusion / loss are **per chromosome**; **translocation is per gene copy** (a rearrangement, like
  transposition — a per-chromosome translocation rate is a possible future addition, deferred); origination is **per genome** (matching the
  code's own weighting — `_choose_chromosome_weighted` vs `_choose_chromosome_uniform`, `genome.py:471`;
  and today's help text, `cli/genomes.py:188`). Override with a count wrapper, e.g.
  `fission = PerGenome(0.02)` for one fission budget per genome regardless of chromosome count, or bend
  with a modifier: `loss = 0.01 * mod.ByChromosomeSize(...)` (bigger replicons die faster — the code
  already size-weights the *pick*; a modifier would make the *rate* depend on size).

**What the user gets** — the network, alongside the gene trees and species tree:

```python
result = genomes.simulate_ordered(tree, chromosomes=8, fission=0.02, fusion=0.02, seed=1)

result.chromosome_network            # the eNewick string (§3)
result.karyotype_events              # the ChromosomeEvent edge list (§3)
result.gene_chromosome_path("g42")   # which chromosome lineage carried gene g42, over time (§4)
result.write("out/", write=["karyotype", "layout"])   # writes .enewick + events.tsv + gene_order.tsv
```

The `write=["karyotype"]` selector already exists (`simulation.py:481`); it gains the `.enewick`
artefact and the speciation edges. Single-chromosome, no-tier-event runs write nothing new — the
network output is opt-in exactly when the karyotype is non-trivial (`cli/genomes.py:114`), preserving
byte-identity for the common case.

---

## 6. Decisions locked (2026-07-18) — was: open questions

All three questions that needed Adrián are decided:

- **Fusion reticulation → recover it descriptively (Option A).** The chromosome network is output as
  eNewick `#H` + an event table, every node **branch-stamped** to its species branch; the reticulations
  live *below* the species level and are **not** projected onto species branches. There is **no formal
  reconciliation engine** for the reticulation in v1 — the network is *recovered, not reconciled*.
  (Adrián: recover the network; chromosome realism is not the point.)
- **Fission/fusion id → re-mint both (§2); the byte-identity break is accepted.** The redesign does not
  preserve old seeds, so the previously-guarded byte-identity invariant does not apply here.
- **Translocation → per gene copy (§5)** — a rearrangement, like transposition. It is surfaced in the
  gene layer with a chromosome cross-reference (§4), not as a horizontal edge in the network. A
  *per-chromosome* translocation rate is a possible future addition — **deferred, not now**.

Still deferred (not v1):

- **eNewick reader.** eNewick `#H` is the chosen output (not a bespoke edge-list); a network-aware reader
  is needed because tree tooling (Newick, degree-2 suppression, LCA) does not apply.
- **Extant-pruned network.** Ships the *complete* network only for v1; the extant-pruned network (keep a
  reticulation iff either parent path survives) is deferred.

- **One network per genome, or one global network?** A chromosome lineage lives inside a single evolving
  genome that itself splits at speciation, so the natural object is **one network per gene-family
  simulation run**, rooted at the seed karyotype and branching with the species tree — *not* one per
  species and *not* one per chromosome. This mirrors "one event log per run" (`events.py:223`). But
  transfer moves genes between genomes without moving chromosomes, so the network stays
  within-genome-lineage — good. Confirm this is the intended scope (vs a per-leaf snapshot of "which
  chromosomes exist now").

- **Speciation fan-out cost.** Recording a `SPECIATION` chromosome edge per parent chromosome per
  speciation adds `O(n_chrom × n_speciations)` rows. For a multi-chromosome eukaryote-scale run this is
  the dominant term of the karyotype stream. Likely fine (it is what makes identity work), but worth a
  `write=` opt-out that falls back to today's within-branch-only trace.

- **Nucleotide vs ordered parity.** Both models have the tier (`nucleotide_genome.py:1079` fission,
  `:890` translocation) and both emit `ChromosomeEvent`s through the shared `genome_sim.py:523` path, so
  the network is model-agnostic — but the speciation-edge fix must land in **both** `clone_reminting`s
  (ordered `genome.py:801` and the nucleotide equivalent). Confirm the nucleotide clone re-mints
  chromosomes the same way.

---

## What to build

1. **Record the speciation edge.** In `clone_reminting` (`genome.py:801` + nucleotide twin), emit a
   `ChromosomeEvent(SPECIATION, branch, t, parents=(parent_cid,), children=(d1, d2))` per parent
   chromosome — the lockstep `zip` already has the mapping; stop discarding it. This is the single change
   that gives chromosomes a genealogy across the whole run.
2. **Re-mint both children at fission and fusion** (§2), so no `chrom_id` spans an event and every node
   has a clean birth. (Gated on Adrián accepting the byte-identity break.)
3. **Assemble + serialise the network.** A `chromosome_network` builder over
   `EventLog.chromosome_records` → eNewick (`karyotype_network.enewick`) + the extended
   `karyotype_events.tsv`. Node labels carry the species `branch` (reconciliation-by-annotation, free).
4. **Gene-→-chromosome path.** Use `Region.chromosome` / `Region.dest` (already logged) with stable
   lineage ids to expose `gene_chromosome_path(gid)` and colour `gene_order.tsv` by chromosome lineage.
5. **API surface.** On `simulate_ordered` / `simulate_nucleotide`: `chromosomes=`, `topology=`, and the
   four tier rates in the `count(base) × modifiers` grammar (retire `n_chromosomes`,
   `--linear-chromosomes`, and the bare-float tier flags into the keyword-rate grammar).

## Still to design

- **Reconciling a reticulation into the species tree** — what the artefact *is* when a fusion node has
  two parents (§6). The single hardest open question; needs Adrián.
- **The extant-pruned network** — definition and whether it ships (a reticulation survives iff either
  parent path survives).
- **Fission/fusion id-reuse vs re-mint** — the byte-identity trade (§2, §6). Needs Adrián.
- **Translocation's place** — chromosome-network horizontal annotation vs gene-layer-only (§6).
- **Network-aware tooling** — reader/validator for the eNewick (existing gene-tree code cannot be
  reused); whether to lean on PhyloNetworks/Dendroscope conventions or ship a minimal own reader.
- **Names** — the count wrappers for the tier (`PerChromosome` / `PerGenome`), and the modifier for
  size-dependent rates (`ByChromosomeSize`, not "PerChromosomeSize" — "per" is reserved for counts,
  `SPEC §5`).
