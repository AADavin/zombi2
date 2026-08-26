# Genomes I: gene families

Genomes live inside the species tree, and they can be simulated at three **resolutions**, one per chapter: **family** here, **ordered** in Chapter 4, **nucleotide** in Chapter 5. The simplest is a gene family evolving along the tree; the most detailed tracks every nucleotide across several chromosomes.

The **family** resolution is genomes made of gene families and nothing more: no position along a chromosome, no DNA sequence. Genes can originate anywhere in the tree, and then they evolve by duplication, transfers and losses.

## The four events

A genome at the family resolution evolves by four kinds of event, applied to every lineage as it runs down the species tree:

- **Duplication.** A gene copy is copied, so its family gains a member in that lineage.
- **Transfer.** A copy jumps from one lineage to another that is alive at the same moment.
- **Loss.** A gene copy is deleted; a family that loses its last copy is gone from that lineage.
- **Origination.** A brand-new family appears in a lineage, with one copy.

![The four events, one to a panel, labelled beneath. Each grey band is a lineage of the species tree, each black line one gene copy, and time reads left to right: a duplication forks a copy inside its lineage, a transfer carries one into another lineage, a loss ends one, and an origination starts a family from nothing.](figures/four_events_print.png){width=100%}

The initial genome, at the beginning of the stem, starts with `initial_families` families of one copy each (100 gene families by default); from there the four rates determine how they evolve.

### Four things to know

- **Families can evolve at different rates.** This can be used to simulate for example genomes with a core and an accessory part ([Ge4](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_pangenome_by_family-->). Appendix A gives the spelling.
- **A family's copies within one genome are capped**, at `max_family_size=10` by default, because a duplication rate above the loss rate can grow exponentially.
- **There is no minimum number of genes.** Loss is counted per copy and the last copy is a copy like any other, so a high loss rate can leave a lineage with nothing. `genome_summary.json` reports `empty_genomes` and the command warns, because an empty genome is otherwise invisible.[^chromosome-minimum]
- **ZOMBI2 produces both complete and extant trees for the gene trees too.** The complete tree keeps every gene lineage; the extant tree keeps only the copies that survive to the present. Both are described below.

[^chromosome-minimum]: The chromosome-based resolutions do keep a minimum: one gene per chromosome. A loss never takes a chromosome below its last gene, because a chromosome with nothing on it is not a chromosome. The minimum is per chromosome, not a bound on the genome's size: a whole chromosome, genes and all, still dies by `chromosome_loss` (Chapter 4).

## What the rate depends on

By default, duplication, transfer and loss are counted **per copy**: a family with ten copies is ten times as likely to duplicate or lose one as a family with a single copy, which is usually what you want: more genes, more chances. Origination is counted **per lineage** (per branch of the species tree): acquiring a wholly new family is a property of the lineage, not of any gene it already has.

The three events (D, T, L) can also be expressed as `PerLineage`, but that implies a different model of evolution. For example, `PerCopy(0.25)` means that every copy can undergo an event independently, so a big genome loses often; `PerLineage(0.25)` is different: the lineage loses at that rate whatever it holds, so deletion-biased genomes lose at their own pace and shrinking does not slow them down. In a genome of a hundred genes the same `0.25` is a total loss rate of 25 one way and of 0.25 the other, a hundredfold apart.

A rate can also depend on **time**, on **the clade where the gene is**, or on a different level. Using `changing_at` on a rate makes it change at set moments: the skyline genome, fast early and slow later; `Clade` multiplies a rate inside a named subtree of the species tree, and the two compose ([Ge5](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_clade_transition-->); and `scaled_by` makes a rate depend on a driver, which is Chapter 8. See Appendix A for more details.[^origins]

[^origins]: `families=[family("toxin", origin=("n1", 0.4))]` declares a family and places it instead of leaving it to the rate: it originates on the lineage you name, at the time you give, keeps its name, and **adds to** whatever `initial_families` and `origination` produce. With both of those at zero the tree carries exactly the families you declared.

## Lateral gene transfers

When a transfer event occurs, the copy that triggers the events is one of that family's live copies, anywhere on the tree, and it is copied to another lineage that is **alive at that same instant**.

Three arguments are important for transfers:

- **`transfer_to`**, who receives the copy. `"uniform"` (the default) picks any other contemporaneous lineage with equal chance; `"distance"` makes closer relatives likelier: a recipient at distance *d* from the donor gets weight exp(−decay · *d* / depth), the depth being the tree's mean root-to-tip time, so the same `decay` means the same thing on trees of different timescales. `"distance"` is `Distance(decay=1.0)`, and `Distance(decay=)` sets its own. A third rule, `Clades(...)`, weights recipients by **named clades of the tree**, described below.
- **`replacement`**, what happens on arrival. By default the incoming copy is **additive**: the recipient simply gains a copy. With `replacement=True` it **overwrites** a copy of the same family already present, picked at random when the recipient holds several, and falls back to additive when it has none.
- **`self_transfer`**, whether a lineage may donate to itself. Off by default. With additive arrival the lineage gains a copy, so the gene content changes as it would under a duplication, but the event is recorded as a transfer. Under `replacement` the arrival never overwrites its own continuation: it displaces one of the lineage's *other* copies, or falls back to additive when there is none.

Transfers can arrive **from a lineage that later goes extinct** [@szollosi2013lgtdead]. A genome run happens on the complete tree, dead branches included, so a gene can enter a survivor from a donor that leaves no other trace. That is the feature the software was originally named for: transfers from the dead, from zombie lineages.

### Transfer between named clades

`"distance"` biases transfer by relatedness, but sometimes you want to name the groups yourself: "let genes flow between these two clades, and nowhere else." `Clades` does that. You name each clade, by a few of its tips (the clade is the subtree below their MRCA) or by the id of its root node. Either way a lineage is written as the `n<id>` label the tree files use, read off `species_complete.nwk` or `species_events.tsv`, so the example below names each clade by two of its tips. The weights go in a `Between` kernel, a table with one weight for each ordered **(donor clade, recipient clade)** pair ([Ge3](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_transfer_highway-->).

```python
from zombi2 import species, genomes

sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=16, seed=1)
# genes flow only BETWEEN clade A and clade B, never within either, never to the rest
flows = genomes.Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
g = genomes.simulate_genomes_family(
    sp, transfer=1.0, initial_families=20, seed=2,
    transfer_to=genomes.Clades({"A": ["n27", "n28"], "B": ["n21", "n26"]}, flows))
```

There is one group you do not have to name: **`"rest"`** is every lineage outside the clades you did name, so `Between({("rest", "A"): 8.0})` makes clade A a transfer hotspot that the whole rest of the tree donates into, without having to enumerate the rest of the tree.

Each entry is a weight, read the same way `"distance"`'s weights are: normalised over the lineages alive at the instant a transfer occurs. The weight belongs to each candidate lineage, not to the clade as a whole, so at the same weight a clade of twenty lineages receives ten times the flow of a clade of two. Naming only `("A", "B")` and `("B", "A")` and setting `default=0.0` means every other pairing weighs 0: a clade-A donor can reach clade B but not another clade-A lineage, and the rest of the tree neither sends nor receives. Drop the `default=0.0` and unlisted pairs return to weight 1 (baseline), so `Between({("A", "B"): 5.0})` *enriches* A→B fivefold while leaving everything else to happen normally. A weight of 0 means "cannot receive", exactly as in Chapter 8: when a donor's every candidate weighs 0, the transfer has nowhere to land and does not happen.

## Profiles and gene trees

A run hands back two views of the same history, and the pair is worth keeping straight: the genomes at the *extant* tips are the observed dataset, and the genomes at **every** node, extinct and internal alike, are the run's own record of what happened. Appendix B names them and everything else.

**Profiles** are the classic comparative-genomics view [@pellegrini1999profiles]: how many copies of each gene family sit in each extant species.

```python
g.profiles.matrix        # families × extant-species copy counts, a NumPy array
g.profiles.presence      # the same as 0/1 presence/absence
g.profiles.to_tsv()      # the table as text
```

**Gene trees.** Every family has a `.complete` tree with every gene lineage, and, when at least one copy survives, an `.extant` tree pruned to the genes that do. A family that died out has no extant tree: `.extant` is `None`, `to_newick("extant")` returns `None`, and the run writes no `_extant` file for it.

```python
gt = g.gene_trees[7]                 # the gene tree of family 7
gt.to_newick("extant")               # the surviving copies as Newick ...
gt.to_newick("complete")             # ... or the whole genealogy
gt.origination                       # when the family began
```

The root of a gene tree carries a branch length, as the species tree's does. A family starts at its origination, and the founding gene lives a while before its first duplication, transfer or speciation; that wait is the root's branch, the family's stem. A gene that originated and never split at all is a one-node tree, written as its own lifespan, `n26_g545:0.6614353;`.

## Evolving families in parallel

Because families are independent, a run can evolve them **concurrently**, one family per worker process. It is off by default; `parallel` turns it on.

`parallel=True` uses every core and an integer sets the worker count; on the command line it is `--parallel` for all cores or `--parallel 8` for eight. It is a **separate engine**, not a faster path through the default one: each family draws from its own random stream, so the result is identical for any worker count, but it differs from a serial run of the same seed. Both are valid draws of the same process.

For very large runs of hundreds of thousands of families, or a million, the difficulty stops being speed and becomes memory. `stream_to` writes each family to a directory the moment it is done, and hands back a light handle holding a path rather than the ordinary result object (a `FamilyGenomesResult`) holding everything. Memory then stays flat however many families you run, so a run that would have held 2 GB in memory streams in about 40 MB, and the sequence level reads the families back off the disk afterwards. On the command line this is `--stream`.

```python
run = genomes.simulate_genomes_family(
    tree, origination=2.0, initial_families=5000, seed=1,
    parallel=8, stream_to="out/", outputs=("events", "profiles", "species_tree"))
run.path("events")            # out/genome_events.tsv — the log the sequence level replays
```

## On the command line

`zombi2 genomes` evolves gene families along the species tree already in the run directory, or one given with `--from` (a Newick file, or another run). The family resolution is the default, and a rate is a plain number or a quoted expression in the written form of Appendix A.

```bash
# duplication–loss–origination along a species tree
zombi2 genomes out/ --duplication 0.2 --loss 0.25 --origination 0.5 --seed 1

# horizontal transfer biased toward close relatives, overwriting resident copies
zombi2 genomes out/ --transfer 0.5 --transfer-to distance --replacement \
    --origination 0.4 --seed 3
```
