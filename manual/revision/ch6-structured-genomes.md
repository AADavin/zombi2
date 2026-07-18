# Genomes II — Structured

In Chapter 5 a genome was a bag of genes: which families a lineage carries, and in how many copies, with no sense of where anything sits. That is the **unordered** resolution, and it is all you need when families evolve independently of one another. This chapter adds the one thing it left out: **position**. Give the genes an order along a chromosome and a whole class of events opens up that were meaningless before — reordering genes without changing which genes are there, splitting and merging chromosomes, and, at the finest grain, the actual nucleotides that genes and the spaces between them are made of.

There are two structured resolutions, and they are layers, not rival models:

```python
from zombi2 import genomes
genomes.simulate_unordered(tree, …)     # a multiset of gene families (Chapter 5)
genomes.simulate_ordered(tree, …)       # + position: rearrangements, chromosomes
genomes.simulate_nucleotide(tree, …)    # + DNA: genes, intergenes, indels
```

**unordered ⊂ ordered ⊂ nucleotide.** Everything from Chapter 5 (the four events, the rate grammar, transfers, conversion, growth caps, gene trees and the event log) is still here and still works the same way; each resolution just *adds* to the one below it. `simulate_ordered` runs the unordered core and then reorders its genes on chromosomes; `simulate_nucleotide` runs that and then gives every gene and gap a length in base pairs. Because they share one core, the shared behaviour cannot drift between them.

> *[Draft — the three-function surface `simulate_unordered` / `simulate_ordered` / `simulate_nucleotide` is the design target of `docs/design/genome-api.md`; it is not built yet. Today the same jobs are done by `simulate_genomes(..., genome_resolution="ordered")` and by a separate `simulate_nucleotide_genomes`. The ordered and nucleotide argument sets in particular are still under design (`genome-api.md`, "Still to design"), so the specific spellings below are the intended shape, flagged where they are not yet settled. This chapter documents the design; today's code is noted where it differs.]*

## What order adds

An ordered genome keeps its genes in a line (or a ring) with neighbours: gene *B* sits between *A* and *C*, and that adjacency is now part of the state. Nothing about the *content* changes (the same families in the same copy numbers), but the genome now has a geography, and geography is what synteny studies read. Two lineages that share a run of neighbouring genes share a **syntenic block**; the events below are exactly the ones that make and break those blocks.

The gene events you already know still fire, and they respect the order now: a duplication drops its new copy in **tandem**, right beside the original, rather than anywhere in the bag; a transfer lands its segment at a position; a loss removes a contiguous run. What is genuinely new are the events that move genes around without adding or removing any.

## Rearrangements

Two events reorder genes on a chromosome and leave its content untouched:

- **Inversion** picks a contiguous segment, reverses it, and flips the strand of every gene in it — the segment is now read back-to-front, as if the DNA had been cut out, turned over, and pasted back. This is the workhorse of bacterial gene-order evolution.
- **Transposition** cuts a segment out and pastes it somewhere else on the *same* chromosome. By default it keeps its orientation; a `transposition_flip` probability lets it reinsert reverse-complemented.

Both are counted the same way as duplication, transfer and loss, **per copy** (so a larger genome rearranges proportionally more often), and both take the same `scope(base) × modifiers` grammar as every other rate in the book:

```python
genomes.simulate_ordered(tree,
    duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,   # the Chapter 5 events
    inversion=3e-4, transposition=1e-4,                          # + rearrangements
    seed=1)
```

An event does not have to touch a single gene. A **per-event length distribution** sets how many neighbouring genes a rearrangement carries at once — `inversion_length=Geometric(mean=5)` flips runs of about five genes at a time — and any distribution constructor works (`Geometric`, `Poisson`, …), chosen separately per event type. A mean length of one recovers single-gene events, which is the default.

> *[Draft — the event-rate names (`inversion`, `transposition`) are the intended shape but not yet pinned in `genome-api.md` (the resolution-specific arguments are still to design); today's CLI spells them `--inversion`, `--transposition` and `--transposition-flip`, as bare floats requiring the per-copy count. The per-event length distribution `inversion_length=Geometric(mean=5)` **is** decided (`genome-api.md`, 2026-07-18); today's CLI still has only a single global `--mean-length` (in genes). The `scope(base) × modifiers` spelling is the intended target.]*

## Chromosomes

Up to here every gene lived on one chromosome. An ordered genome can carry several, and a chromosome is a thing in its own right — it has an identity, a **topology**, and events that act on whole chromosomes rather than on the genes inside them. You seed a karyotype with a count and a topology:

```python
genomes.simulate_ordered(tree,
    chromosomes=8,             # eight chromosomes at the root
    topology="linear",         # "circular" (default, bacteria) | "linear" | a list per chromosome
    …)
```

Topology is either shared by every chromosome or given per chromosome as a list (`["circular", "linear"]` for a mixed karyotype). It is just a label that decides which operations are legal: a **circular** chromosome is a ring with no privileged origin, so a segment may wrap past the origin; a **linear** one has two ends and its segments are clamped to them. ZOMBI2 makes no claim about origins, centromeres or telomeres — the topology is the whole of the biology it models here.

Four events act at the chromosome layer (the code calls it the *chromosome tier*), each with a natural scope:

- **Fission** splits one chromosome into two. A linear chromosome is cut at one breakpoint; a circular one at two, excising the arc between them into a new ring. Counted **per chromosome**, size-weighted, so a bigger chromosome splits more often.
- **Fusion** merges two chromosomes into one, appending the genes of the second onto the first. Only two chromosomes of the *same* topology may fuse — joining a ring to a line is ill-defined — and it too is counted **per chromosome**.
- **Chromosome origination** adds a new, empty replicon: a de-novo plasmid that genes reach later by origination or transfer. Counted **per lineage** (one genome per lineage), since it is the genome as a whole that acquires one.
- **Chromosome loss** deletes an entire chromosome and every gene on it; the genes die as ordinary losses in the gene log. Counted **per chromosome**.

```python
genomes.simulate_ordered(tree, chromosomes=8, topology="linear",
    fission=0.02, fusion=0.02,
    chromosome_origination=0.01, chromosome_loss=0.01, seed=1)
```

A fifth event, **translocation**, moves a segment of genes from one chromosome to another within the same genome — the cross-chromosome cousin of transposition, and, like transposition, counted **per gene copy**. It needs at least two chromosomes, and it is deliberately *identity-neutral*: both chromosomes persist, only the genes change address. It matters because it is the one way a gene can jump from one chromosome to another without being duplicated or transferred, which is exactly what makes the chromosome layer connect to the gene layer (below).

Each rate takes the same grammar as everything else; the defaults answer "per what?" so you rarely write a wrapper. To override, wrap the base (`fission = scope.PerLineage(0.02)` gives one fission budget per lineage — one genome per lineage — regardless of how many chromosomes it has), or bend it with a modifier (`loss = 0.01 * mod.ByChromosomeSize(...)` makes bigger replicons die faster).

> *[Draft — `chromosomes=` and `topology=` replace today's `--n-chromosomes` and `--linear-chromosomes`; the four tier rates are already real flags (`--fission`, `--fusion`, `--chromosome-origination`, `--chromosome-loss`) but as bare floats, not the wrapped grammar. The scope wrapper `scope.PerChromosome` and the modifier `mod.ByChromosomeSize` are decided but not yet built (`PerGenome` was dropped in favour of `PerLineage` — one genome per lineage). Translocation's scope is now decided: **per gene copy** (a rearrangement, like transposition; matching today's code).]*

## The chromosome network

Speciation, fission, origination and loss each have one parent and at most two children, so on their own the chromosomes of a genome would trace out a **tree** — a karyotype genealogy, one node per chromosome lineage, branching as chromosomes split and speciate. Fusion breaks that. A fusion takes *two* parent chromosomes and yields *one* child, and a node with two parents is a **reticulation** — the mark of a network, not a tree. So the true genealogy of a genome's chromosomes is a **chromosome network**, and fusion is the single event that makes it one.

This network is the missing middle of the three-way nesting the simulator records:

```
species tree  ⊃  chromosome network  ⊃  gene trees
```

A gene lives on a chromosome; a chromosome lives in a species. The gene trees say *which genes are related*; the species tree says *which organisms*; the chromosome network says *which replicon carried each gene, and when*. It is the only one of the three that reticulates, and it is the connective tissue between the other two: because translocation and transfer can move a gene from one chromosome lineage to another, a gene's walk across chromosomes is itself a recoverable history.

A network cannot be written as one Newick string, so ZOMBI2 records it two ways. The ground truth is an **event table**, one row per chromosome event with its time, its species branch, and the parent and child chromosome ids — the audited edge list the whole genealogy is assembled from. From it a topology is serialised as **extended Newick** (eNewick), the standard format for networks: a reticulation node is written once and then referred to by a shared label `#H1` wherever else it appears, so a fusion's two-parent child hangs under both its parents, tagged the same. eNewick is what phylogenetic-network tools read, which is the reason for choosing it over a bespoke format.

Every chromosome event already carries the species branch it fired on, so every node in the network can be **stamped** with the branch it happened on — an exact annotation, never inferred, because the simulator knows the truth. That branch-stamp is how the network is placed against the species tree: the same topology with each node labelled by its species branch, exactly as a reconciled gene tree stamps each node with the event that made it. What v1 does *not* do is formally **reconcile** the reticulation — a fusion node has two parent branches, and projecting that two-parent join onto the species tree is left open. So the chromosome network is **recovered, not reconciled**: branch-stamped and readable, with no reconciliation engine behind the fusion nodes.

> *[Draft — the chromosome network is a design target (`docs/design/chromosome-network.md`), **not yet built**. Today the events are logged as `ChromosomeEvent` records (`events.py`) and written to a flat `karyotype_trace.tsv`, but they are not assembled into a connected genealogy, and one edge is missing entirely: at speciation each daughter chromosome is re-minted with a fresh id and the parent→daughter correspondence is discarded, so a chromosome cannot yet be traced back past its own species branch. Recording that speciation edge is the one change that gives chromosomes a genealogy across the whole run. The three questions that were open are now **decided (2026-07-18)**: the reticulation is **recovered, not reconciled** (branch-stamped eNewick, no formal reconciliation engine for a two-parent fusion node); v1 ships the **complete** network only (the extant-pruned network is deferred); and fission and fusion **re-mint both children** for a clean, uniform genealogy, accepting the byte-identity break on multi-chromosome runs.]*

## The nucleotide resolution

The finest resolution gives the genome an actual sequence. `simulate_nucleotide` evolves a **circular string of nucleotides**, carved into **segments**, each one anchored to the interval of the ancestral genome it descends from. Tracing any present-day base back to its origin is therefore immediate (the segment already names its ancestral coordinates), which is what lets the simulator reconstruct the genome of every ancestral node, not only the tips.

The organising idea is the **block**: a maximal stretch with a single history. Two kinds of block matter, and the difference is what makes this resolution biological rather than a bare string.

- **Genes** are declared intervals, and they are **indivisible**: no rearrangement breakpoint, insertion or deletion is ever allowed to fall *inside* a gene. A gene is therefore exactly one block wherever it survives — one gene, one genealogy — and that invariant is what ties the nucleotide genome back to the gene trees of Chapter 5.
- **Intergenes** are the stretches between genes, and they are freely cut. Breakpoints land in them, and it is there that small insertions and deletions accumulate.

You either declare the genes or let there be none. With `--gff` (a real annotation, e.g. a RefSeq bacterial chromosome) or `--genes` (a BED file), the simulation starts from a genuine genome's architecture; with neither, the whole chromosome is intergene and events act anywhere.

Two events are new at this resolution, both confined to intergenes:

- **Insertion** lays down a run of fresh, random nucleotides — a new block from a new source.
- **Deletion** removes a run from within a single intergene (never spanning or touching a gene).

Both are off by default and counted **per nucleotide**, each with its own per-event length distribution (`insertion_length=` / `deletion_length=`), separate from the segment length that rearrangements use. One further event lives here too: **pseudogenization**, a sub-outcome of loss in which a gene loses its function but keeps its sequence and carries on as intergene, its history unbroken — the natural way to model a gene decaying rather than vanishing.

The rearrangements and chromosome events of the ordered resolution are all here as well, now measured in base pairs rather than gene counts: inversions and transpositions act on nucleotide arcs, and the same fission, fusion, translocation, origination and loss reshape the chromosomes, so a nucleotide run produces the very same chromosome network as an ordered one.

This chapter stops at the structured genome — the segments, blocks, genes and intergenes, and their coordinates. The *letters* themselves, the substitutions that turn one nucleotide into another down each gene tree, are the subject of Chapter 7. A nucleotide genome run can hand its gene trees and its ancestral genome coordinates straight to that level.

> *[Draft — `simulate_nucleotide` is the design target; today the entry point is `simulate_nucleotide_genomes`, and the events above are flags (`--insertion`, `--deletion`, `--indel-mean-length`, `--pseudogenization`, `--gff`, `--genes`). The nucleotide argument set is in `genome-api.md`'s "Still to design".]*

## The literature → command bridge

Readers arriving from the rearrangement and comparative-genomics literature already have names for these events. As in every chapter, the names live in one table and organise nothing.

| From the literature | What it does | Here |
|---|---|---|
| Inversion / reversal (GRIMM, MGR) | reverse a gene segment, strands flipped | `inversion=…` |
| Transposition | cut a segment and paste it elsewhere on the chromosome | `transposition=…` |
| Reciprocal translocation | a segment moves to another chromosome | `translocation=…` |
| Chromosome fission / fusion | split one chromosome, or merge two | `fission=…`, `fusion=…` |
| Plasmid / replicon gain | a de-novo chromosome appears | `chromosome_origination=…` |
| Karyotype / chromosome-number evolution | the chromosome network as a whole | fission + fusion → `karyotype_network.enewick` |
| Indels | insert or delete a run in an intergene | `insertion=…`, `deletion=…` |
| Pseudogene formation | a gene decays to non-functional, sequence kept | `pseudogenization=…` |
| Starting from a real genome | seed the architecture from an annotation | `--gff` / `--genes` |

## The objects

*[Draft — the result API is decided (the `GenomesResult` bundle, `docs/design/result-api.md`); the network accessors below are design targets, not yet built.]*

`simulate_ordered` and `simulate_nucleotide` return a **`GenomesResult`**, a superset of the unordered one's payload: the per-lineage genomes (now with their gene order and chromosomes), the gene trees and reconciliations, and the profile matrix, plus the structured extras. Like every result bundle it shares the common spine (`.events`, `.tree`, `.write(include=[...])`, `.seed`). From an ordered or nucleotide result you can read a leaf's gene order, ask for the **chromosome network** as eNewick, pull the **chromosome event table**, and trace a single gene's path across chromosome lineages:

```python
result = genomes.simulate_ordered(tree, chromosomes=8, fission=0.02, fusion=0.02, seed=1)
result.chromosome_network            # the eNewick string
result.karyotype_events              # the chromosome event edge list
result.gene_chromosome_path("g42")   # which chromosome lineage carried gene g42, over time
```

## Usage from Python

```python
from zombi2 import genomes, modifiers as mod
from zombi2 import scope        # scope wrappers: Global, PerCopy, PerLineage, PerChromosome, …

# an ordered bacterial genome: the Chapter 5 events, plus inversions
genomes.simulate_ordered(tree,
    duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
    inversion=3e-4, seed=1)

# several linear chromosomes that split, merge, and gain a plasmid — the chromosome network
genomes.simulate_ordered(tree,
    chromosomes=8, topology="linear",
    fission=0.02, fusion=0.02,
    chromosome_origination=0.01, chromosome_loss=0.01, seed=1)

# a nucleotide genome started from a real annotation, with intergenic indels
genomes.simulate_nucleotide(tree,
    gff="ecoli.gff",
    inversion=1e-3, insertion=1e-4, deletion=1e-4,
    pseudogenization=0.1, seed=1)
```

### A worked example

A small end-to-end run: grow a species tree, evolve an ordered genome of two circular chromosomes along it with a background of inversions and the odd fission and fusion, then read back the karyotype the run produced.

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_tips=12, seed=1)

result = genomes.simulate_ordered(tree,
    duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
    chromosomes=2, topology="circular",
    inversion=3e-4, fission=0.02, fusion=0.02, seed=1)

print(result.chromosome_network)          # eNewick, with #H1 at each fusion
for e in result.karyotype_events:         # time, event, branch, parents → children
    print(e)
result.write("out/", include=["karyotype", "layout"])
```

The fusions are the reticulations: each one shows up in the event table as two parent chromosome ids collapsing to one child, and in the eNewick as a `#H` node hanging under both its parents. A run with a single chromosome and no fission or fusion has a trivial genealogy and writes nothing new — the network output is opt-in, exactly when the karyotype is non-trivial.

> *[Draft — the worked example uses the target API; today the run is `simulate_genomes(..., genome_resolution="ordered", n_chromosomes=2)` and the network is a flat `karyotype_trace.tsv`, not eNewick.]*

## Usage from the CLI

*[Draft — the CLI re-fit to this API is still to be designed; today's `zombi2 genomes` selects the resolution with `--genome-resolution {unordered,ordered,nucleotide}` and takes the events as bare-float flags rather than the wrapped rate grammar.]*

```bash
# ordered genome with inversions and two chromosomes
zombi2 genomes --genome-resolution ordered -t species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 \
    --inversion 3e-4 --n-chromosomes 2 --fission 0.02 --fusion 0.02 \
    --seed 1 -o my_genomes

# nucleotide genome from a real annotation, with intergenic indels
zombi2 genomes --genome-resolution nucleotide -t species_tree.nwk \
    --gff ecoli.gff --inversion 1e-3 --insertion 1e-4 --deletion 1e-4 \
    --pseudogenization 0.1 --seed 1 -o my_genomes
```

## Outputs

*[Draft — to finalise with Appendix B; the eNewick network is a design target.]*

A structured run writes everything the unordered run does (gene trees, reconciliations, event log, profile matrix) and adds the structure:

- the **gene layout** (`gene_order.tsv`), which chromosome each gene sits on and where, so a leaf genome can be read back in order;
- the **chromosome event table** (`karyotype_events.tsv` in the target; `karyotype_trace.tsv` today), the fission / fusion / origination / loss genealogy, plus the speciation edges once identity is carried across speciation;
- the **chromosome network** (`karyotype_network.enewick`), emitted whenever the karyotype is non-trivial;
- for a nucleotide run, the genes and intergenes as intervals in **GFF/BED**, and, when you ask for it, the reconstructed DNA of every ancestral genome — the point where this chapter hands off to Chapter 7.

The full list of files lives in Appendix B.
