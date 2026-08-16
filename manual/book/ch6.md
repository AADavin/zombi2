# Sequences

The sequence level does two main things:

* It rescales the gene trees and the species tree from time into substitutions per site (**phylograms**).
* It evolves the residues that sit inside every gene, so each family ends with an alignment.

The alignment is the **observable** half: the sequence at every extant gene copy, which is what a phylogenetic method would be handed. The other half is every node the alignment leaves out — the internal ones, and the tips where a copy was lost or its species died. The run wrote a sequence at each node as it went, so those are the exact ancestors rather than estimates, and together the two account for every node of the gene tree exactly once ([Sq9](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:seq_ancestral-->).

![Where a sequence lives. A sequence is not evolved along the species tree: it is evolved along the **gene tree**, which the genome run produced and which runs inside the species tree. The two forks are different events. The first is a speciation, which hands the gene to both daughters; the second, marked with a square, is a duplication, and it happens inside one lineage — which is why one lineage can hold two tips of the same family. A run gives one sequence per node of that gene tree, so the tips come out as an alignment.](figures/sequence_nesting_print.png){width=88%}

## Creating phylograms

A gene tree in ZOMBI2 is by default a **chronogram**, its branch lengths measure time. What a sequence actually accumulates along a branch is not time but a number of *substitutions per site*, and that is time multiplied by an evolutionary rate. Turning one into the other is the whole job of the sequence level. Applying a rate to every branch rescales the tree from time into expected substitutions and yields a phylogram.

Two things therefore have to be chosen: *what* changes (the substitution model, the chemistry of which residue turns into which) and *how fast* it changes along each branch, which is the clock.

## The substitution models

ZOMBI2 implements nine standard models of sequence evolution.

| Model | Alphabet | Parameters | From the literature |
|---|---|---|---|
| `jc69` | bases | none | [@jukes1969evolution] |
| `k80` | bases | `kappa` | [@kimura1980simple] |
| `hky85` | bases | `kappa`, `frequencies` | [@hasegawa1985dating] |
| `gtr` | bases | `exchangeabilities`, `frequencies` | [@tavare1986some] |
| `poisson` | residues | none | equal-rate protein null |
| `dayhoff` | residues | none | [@dayhoff1978model] |
| `jtt` | residues | none | [@jones1992rapid] |
| `wag` | residues | none | [@whelan2001general] |
| `lg` | residues | none | [@le2008improved] |

The model decides the alphabet, and `length` counts whatever that alphabet holds. The nucleotide models are four different rate matrices rather than one model with four settings, but they nest in the order written — `jc69` is `k80` with `kappa=1`, and `k80` is `hky85` with equal base frequencies — so each step down the table adds free parameters. The protein matrices are **empirical**, each estimated once from a large set of real alignments and then used as a fixed table, which is why they take no parameters of their own ([Sq3](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:seq_protein-->).

### Your own matrix

If none of the nine is the model you want, write the matrix yourself. `reversible` takes a symmetric **exchangeability** matrix $S$ and stationary frequencies $\pi$, and returns a model like any other:

```python
from zombi2 import species, genomes, sequences
from zombi2.sequences.substitution_models import hky85, lg

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
my_genomes = genomes.simulate_genomes_family(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)
```

```python
import numpy as np
from zombi2.sequences.substitution_models import reversible

kappa = 2.0
S = np.array([[0, 1, kappa, 1],      # A ↔ C, A ↔ G, A ↔ T
              [1, 0, 1, kappa],      # C ↔ …
              [kappa, 1, 0, 1],
              [1, kappa, 1, 0]], dtype=float)
mine = reversible(S, frequencies=(0.3, 0.2, 0.2, 0.3), name="mine")
custom = sequences.simulate_sequences(my_genomes, model=mine, length=1000, seed=1)
```

The rate of $i \to j$ is $Q_{ij} = S_{ij}\,\pi_j$, the diagonal is minus the rest of its row, and the whole matrix is then scaled so that one unit of branch length is one expected substitution per site, the same scaling every model on the menu gets. So a phylogram from your matrix is comparable with one from `hky85` without converting anything. This is the constructor the menu itself uses: the matrix above *is* HKY85, so `mine` and `hky85(kappa=2.0, frequencies=(0.3, 0.2, 0.2, 0.3))` are the same model. Pass `alphabet=AMINO_ACIDS` for a twenty-state matrix of your own.

You give $S$ and $\pi$ rather than $Q$ directly, and the restriction is deliberate: the engine's matrix exponential is only valid for a **time-reversible** model, which a symmetric $S$ times $\pi$ satisfies by construction. A general $Q$ is refused with an error rather than run — a wrong transition matrix produces plausible sequences that are not the model you asked for. There is no command-line flag for a custom matrix: a twenty-state matrix is 190 numbers, which is a file format rather than an argument.

## Rate variation across sites

So far every site of a gene evolves at the same speed, which is a model no real gene obeys: some positions are held nearly fixed by what the protein has to do, and others drift freely. Two settings say so, and both go on the **model** rather than on the rate. `gamma_shape` gives each site a multiplier drawn from a Gamma with mean 1, so one number sets how unequal the sites are — a small shape means a few fast sites among many slow ones, a large one is nearly flat. `invariant` sets aside a fraction of sites that never change, which a Gamma alone fits badly. The field writes them as suffixes on the model's name and so does ZOMBI2: `+G` [@yang1994variable], `+I`, and together `+I+G` [@gu1995maximum]. On the command line they are `--gamma-shape` and `--invariant`, on any model of the menu.

```python
both = hky85(kappa=2.0).across_sites(gamma_shape=0.5, invariant=0.1)
print(both.name)                  # HKY85+I+G4
sequences.simulate_sequences(my_genomes, model=both, length=1000, seed=1)
```

**Branch lengths do not change.** The Gamma is cut into equal-probability classes, four by default and set with `rate_categories`, and the classes are normalised so the mean rate over all sites is exactly 1, invariant sites included. A branch length is therefore still substitutions per site, now the mean over them, and a run with rate variation and one without, at the same rate, have the *same tree* — what differs is how the change is spread across the columns. One visible consequence: the mean pairwise identity a run reports goes **up** at the same divergence, because the slow and invariant sites keep their matches while the fast ones saturate.

## Site-specific amino-acid profiles

Every model on the menu gives a gene one set of amino-acid frequencies, shared by all its sites. Real
proteins do not work that way: a buried position is some flavour of hydrophobic and essentially never
charged, while the loop next to it takes almost anything. A **profile** says so directly: one set of
frequencies per position, rather than one per gene.

You supply a table with a row per site and a column per amino acid, for whichever families you have
one for:

```python
import numpy as np

# a profile is an (L, 20) array: row i is the amino-acid frequencies at position i,
# with the columns in the model's own alphabet order, model.alphabet
protein = lg()
constrained = np.full((300, 20), 0.002)      # every residue possible, none likely...
constrained[:, protein.alphabet.index("L")] = 1.0   # ...except leucine, at every position

seqs = sequences.simulate_sequences(my_genomes, model=protein, length=300,
                                    profiles={0: constrained}, seed=1)
```

Each row is normalised to sum to 1, and the exchangeabilities (which pairs of amino acids swap
easily, the chemistry the model already encodes) are left alone. Families you leave out evolve under
`model` untouched. A row of exact zeros is refused, because it makes that site's matrix degenerate: a
real profile says a residue is unlikely at a position, not that it is impossible, so add a
pseudocount. A **flat** profile — every row the model's own frequencies — is statistically the model
without one.

Where the numbers come from is your business. From an alignment of the real family, column
frequencies with a pseudocount so nothing is impossible:

```python
def profile_from_alignment(columns, alphabet, pseudocount=0.1):
    """One row per alignment column, from the residues seen in it."""
    counts = np.array([[c.count(a) for a in alphabet] for c in columns], dtype=float)
    counts += pseudocount                       # nothing is impossible, only unlikely
    return counts / counts.sum(axis=1, keepdims=True)

columns = ["LLLIL", "GGGGA", "KKRKK"]           # three columns of a five-sequence alignment
print(profile_from_alignment(columns, protein.alphabet).shape)      # (3, 20)
```

Or from a **protein language model**, which already produces a distribution over amino acids at
every position — exactly the table this section wants. Compute it once, save the array beside your
script, and load it: calling a model at simulation time would put its version and its hardware inside
your run's reproducibility, whereas a saved array keeps the simulation bit-identical from the seed as
usual.

Profiles and `+Γ` are about different things and compose: a profile says **which** amino acids belong
at a site, a Gamma says **how fast** sites change. `partitions=` is the same idea at block grain — the
gene evolves as consecutive blocks of sites, each under its own model, `partitions=[(lg(), 100),
(wag(), 200)]` — and it replaces `model=` and `length=`, which it already answers.

**An amino-acid profile needs a protein model, so it belongs to a family or ordered run.** A
nucleotide genome is measured in base pairs and its blocks are read on either strand, so that
resolution refuses protein models altogether: there is no complement of an amino acid. Profiles there
are over the four bases, one row per base pair of the block, the row count matching the block's
length.

## Insertions and deletions

`insertion` and `deletion` give the sequences **indels**, so an alignment gains gaps: a lineage loses
sites it had, or gains sites the others never did. The rates are **relative to the substitution
rate** — `deletion=0.05` is five deletions for every hundred substitutions a site expects — so a
lineage's clock reaches its indels too, and the number means the same on a tree of any height.
`insertion_extent` and `deletion_extent` say how many sites an event takes, a bare number the mean of
a geometric draw (3 by default); `Fixed(1)` is a single-site indel, and any distribution works.

```python
gapped = sequences.simulate_sequences(my_genomes, model=hky85(2.0), length=300,
                                      insertion=0.05, deletion=0.05, seed=4)
```

The alignment columns are the union of every site any lineage ever had, and a lineage that carries
none of one shows a gap. Two refusals, each because something else already owns the sites: a
**nucleotide** genome brings its own indels (Chapter 5), where a base pair has a position these rates
do not know; and `partitions` and `profiles` are written against a site count an indel changes.

## Relaxed molecular clocks

Rate variation across sites says which *positions* change fast. A clock says which *lineages* do. The two are orthogonal and compose.

The rate is `substitution`, counted **per site**: a gene-tree branch of Δ*t* time accrues `substitution · Δt` substitutions at every site. Left alone it is `1.0` everywhere — the **strict clock**, one tempo for the whole tree. Let it change from lineage to lineage and it is what the field calls a **relaxed clock** [@lepage2007general], written by chaining a verb onto the rate exactly as at every other level. Whichever distribution you give is **normalised to mean 1**, so it contributes its *shape* and the base keeps meaning the average rate.

| What it does | From the literature | Gallery |
|-------------------|--------------------------------|---|
| one rate everywhere | Strict / global clock | [Sq1](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:seq_alignment--> |
| each lineage i.i.d. lognormal | Uncorrelated lognormal (UCLN) | [Sq5](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:clock_ucln--> |
| each lineage i.i.d. gamma | Uncorrelated gamma (UGAM) | [Sq6](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:clock_ugam--> |
| the rate drifts parent to daughter | Autocorrelated lognormal | [Sq7](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:clock_autocorrelated--> |
| the rate takes one of a few values, inherited in steps | Discrete rate categories | [Sq8](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:clock_discrete_bin--> |
| the rate reads another level | Trait-dependent rate of molecular evolution | [Co8](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_substitution--> |

Row by row, the written forms: the uncorrelated lognormal clock is `PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))`, and swapping the law for a `Gamma` gives UGAM; wrapping the law in `Drift(...)` gives the autocorrelated clock, and `Drift(..., bins=6)` the discrete categories (Appendix A spells each out); the driven rate is Chapter 8's. A bare distribution has **no memory** — each lineage is an independent draw — where `Drift` starts a daughter at its parent's rate and takes one step from it, so close relatives evolve at similar rates. And whichever you use, **the clock belongs to the species tree, not to the gene trees**: ZOMBI2 draws one rate per species branch, and every gene passing through that branch evolves at it, so a fast species is fast in all of its genes at once.

### Setting the rate, or letting the divergence set it

`substitution` is the number people most often get wrong, because the right value depends on the height of the tree it runs down, so no default can suit every tree. On a tree 20 time units tall the default `1.0` puts **20** substitutions on every site from origin to tip. Those are twenty *events*, not twenty visible differences: a site has four states, so it soon lands back on bases it already held, and past about one substitution per site two sequences are as different as two random ones. The history is not missing — it has been overwritten. Read the rate off backwards instead, from the divergence you want:

$$\text{substitution} = \frac{\text{substitutions per site you want, origin to tip}}{\text{height of the tree}}$$

On that 20-unit tree, sequences around 80% identical need `0.2 / 20 = 0.01`. The difference is not subtle: `0.01` gives tips about 77% identical and the default `1.0` gives **25%**, which is the score of two sequences with no shared history at all.

| Substitutions per site, origin to tip | Mean identity between two tips |
|---|---|
| 0.02 | 98% |
| 0.09 | 90% |
| 0.17 | 81% |
| 0.34 | 68% |
| 0.86 | 46% |
| 1.72 | 34% |

Measured on simulated JC69 alignments. The first column carries over between trees, being rate times height; the second is a mean over pairs of tips, so it depends on where the tree's splits sit, and at the deep end two ordinary trees can differ by twenty points. Read it as a guide, not a lookup. Two unrelated DNA sequences already match at 25% by chance, so the last row is nearly no signal; every run reports the identity it produced and warns when it lands that close to the floor.

`divergence` states the first column directly and lets ZOMBI2 do the division against the height of the tree — `divergence=0.2` in Python, `--divergence 0.2` on the command line. The height is the **species tree's**: the clock is drawn per species branch, so a family born late simply accrues less.

```python
sequences.simulate_sequences(my_genomes, model=hky85(), length=1000, divergence=0.2)
```

The two settings say different things and can be given together: `substitution` says what *kind* of clock, `divergence` how much change it must deliver. A relaxed clock calibrated to a divergence is written with the shape alone — `PerSite().varying_among('lineages', LogNormal(0.0, 0.3))` — and `divergence=0.2` beside it. Giving that rate a base number as well is an error rather than an override, since the base is exactly what `divergence` solves for; the resolved rate goes into the run log either way. One limit: `divergence` is refused alongside a **driven** rate. The clocks are normalised to mean 1, so the division still lands; a driver's factors average to whatever the driver does, so the solved base would deliver a divergence off by the driver's mean while the log claimed otherwise. Set the base yourself there — a driven rate is conditioning, two ordinary runs in order, and Chapter 8 covers it, including how a driver that switches partway along a branch is integrated across the branch rather than read once for it.

### A clade can evolve differently, not only faster

Everything above changes how *fast* a lineage evolves. `Models` changes what the change *looks like*: which residues turn into which, and what composition the sequence settles at.

```python
from zombi2.params import Clade
from zombi2.sequences import Models
from zombi2.sequences.substitution_models import hky85

at_rich = hky85(kappa=2.0, frequencies=(0.40, 0.10, 0.10, 0.40))   # equilibrium A+T = 0.80
model = Models().set_by(Clade({"endo": ["n76", "n112"]}),
                        {"endo": at_rich, "rest": hky85(kappa=2.0)})
```

Every group the clade paints needs a model, `"rest"` included — a lineage in no named clade is in `"rest"`, and a missing one is refused rather than filled in.

This is what an endosymbiont study needs. Its clade evolves faster *and* drifts toward AT, and the two mislead a tree-builder differently: a fast branch causes long-branch attraction, an AT-rich one causes compositional attraction, where unrelated AT-rich lineages group together because they look alike. Scoping the rate to the same clade gives both at once: `substitution=PerSite(1.0).scaled_by(clade, {"endo": 3.0, "rest": 1.0})` beside the model set.

Two limits. The alphabet and the across-site rate classes are shared, because one gene copy's sequence is one string and a site's class is drawn once for the family — asking for either to differ by clade is refused with the reason. And the driver of a model set must be a `Clade`: this level samples one transition matrix per branch, so a model that switched mid-branch, as a trait can, would need the branch cut at the switch — that is not built. A trait can still drive the *rate*.

One caveat before you read the output: every model is normalised to one expected substitution per site per unit branch length, and that holds **at stationarity**. A lineage that has just entered the AT-rich clade is still relaxing toward its frequencies, so while its composition catches up it accrues somewhat fewer substitutions than its branch length claims. That transient is usually the very thing being studied — the ordinary price of a model that varies along the tree, not a defect.

`Models` is Python only, and there is no flag for it: a model is a K×K matrix, which has no written form a flag can carry — the same reason `reversible()` has none.

### Where a sequence starts

A family does not begin at the first branching of its gene tree. It begins when it originates, and the founding gene then lives for a while — its **stem** — before anything splits it. So that is where the sequence starts: one draw from the model's stationary frequencies at the origination, which then evolves across the stem in the ordinary way and arrives at the root gene as a sequence that has already changed. `.founding` is that first draw; `.ancestral` holds what the root gene ended up with, and the two differ by however much the stem allowed.

This is why a phylogram's root carries a branch length: it is the stem in substitutions per site, exactly as every other branch is its own stretch of time converted by the rate. Under a strict clock of rate 1 a phylogram is its chronogram, root branch included.

## Keeping the history

Every other level records what happened to it. This one records the letters at every node and no events, and that is deliberate: a substitution log is bigger than the alignment it explains — a hundred families down an ordinary tree is close to a million rows, where the genome log for the same run holds a few thousand. `record=True` asks for it:

```python
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import jc69, simulate_sequences
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, n_extant=10, seed=1).complete_tree
g = simulate_genomes_family(ct, initial_families=4, duplication=0.05, loss=0.1, seed=2)

r = simulate_sequences(g, model=jc69(), length=300, seed=1, record=True)
r.events[0]                 # time · kind · lineage · gene · site · strand · after · from · to
r.write("out/", outputs=("alignments", "events"))    # sequence_events.tsv
```

`kind` is `substitution`, `insertion` or `deletion` — the indels above get rows too, one per site — and `from` and `to` are the states. A site is named by an **id**, not by a position: a position shifts with every insertion above it, and the alignment column index is stable only once the run has finished, so neither can go in a row written while the run is going. An insertion row carries the id it follows, and a deletion row the ids dropped; with those two a reader can rebuild any lineage's column list at any moment. The `strand` column is there because the nucleotide level's blocks can sit reverse-complemented; at this level it is always `+1`.

Recording changes the **sampler**. An ordinary run draws each branch's end in one step and never simulates the path between, because nothing asks what it was; a recorded run has to walk that path, event by event. The two give the same distribution at a branch's end, so a recorded run is a valid run — it is a *different realisation* for the same seed, the price of asking what happened rather than only where it ended.

What it will not walk, it refuses: `partitions` and `profiles`, which give a family several models; a nucleotide genome run, which evolves blocks rather than one sequence per family; and the parallel and streaming engines, which hand a family off before its rows could be collected. A **joint** run records as well — `joint=True` here, or `joint.simulate(..., record=True)` for a trait and a gene together (Chapter 9); those runs advance in slices of species time, and a slice is exactly the interval their rate is constant over, so a row's time is as exact there as on an ordinary branch.

## Large runs

This is the level where a run's memory goes. Every family's alignment and every ancestral sequence are held at once, so what you can run is bounded by families × copies × sites rather than by time. `stream_to` writes each family's files the moment it is finished and keeps nothing, handing back a light handle with a path instead of a `SequencesResult` holding everything. On the command line it is `--stream`.

Memory then stops growing with the sequences: raise the sites fivefold and an in-memory run doubles while a streamed one barely moves. What is left is the genome run being read, which is the remaining cost. It is a memory choice and not a modelling one: the same seed writes the same files either way, so a streamed run and an in-memory one are the same dataset. `outputs=` picks which files, exactly as `.write` takes them, and it composes with `parallel`. A **nucleotide** run cannot stream: it puts whole genomes back together, which needs every block's sequence at once.

## Running on a nucleotide genome

Hand the level a **nucleotide** genome run — one ZOMBI2 drew, or a real annotation you supply — and you get whole assembled genomes in FASTA, at every tip and at every ancestor. Genes and spacer get their own models: `model` evolves the genes, and `intergene_model` evolves the spacer at `intergene_speed` times the rate (3× by default), under `jc69` unless you hand it a model of its own — decorated or plain, like any other.

It is the same two steps as any other run — a genomes run at `--resolution nucleotide`, then a sequences run over it. What comes back is every node's assembled genome: `.genomes` at the tips, `.node_genomes` at the ancestors, reconstructed rather than estimated, and `.initial_genome` for the state the run started from.

### Starting from a real sequence

So far the founding sequence of each block is *drawn* from the model's frequencies, random ACGT. Hand the genomes run a **FASTA** alongside the GFF (`fasta=`, or `--fasta`) and it starts from the sequence you supply instead — the layout from the annotation, the letters from the DNA.

The FASTA has one `>seqid` record per GFF `##sequence-region`, each exactly its declared length. Every block is then founded from the real DNA at its own initial coordinates, so an assembled genome descends from exactly what you gave. A gene that origination invents mid-run has no supplied DNA (it did not exist initially), so its block still draws from the model.

## On the command line

The positional directory is the run being written, and it doubles as the input: `zombi2 sequences out/` runs **in place**, reading the genome run already in `out/` — its species tree and event log, replaying the gene genealogy from them — and writing beside it. `zombi2 sequences seqs/ --from out/` reads one run and writes into a fresh directory. A nucleotide run given a `--fasta` also hands its initial DNA across (in `initial_sequence.fasta`), so the sequences descend from your real sequence without you naming it twice. `--write` picks the outputs, the command line's `outputs=`.

```bash
# 1. genomes along a species tree (from the previous chapters)
zombi2 genomes out/ \
    --duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5 --seed 1

# 2. HKY85, 1000 sites, strict clock
zombi2 sequences seqs/ --from out/ --model hky85 --kappa 2.0 \
    --length 1000 --seed 1

# GTR with unequal frequencies under a relaxed clock, also writing the ancestral sequences
zombi2 sequences seqs/ --from out/ --model gtr \
    --frequencies 0.3 0.2 0.2 0.3 \
    --substitution "PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))" \
    --seed 1 --write alignments phylograms species_phylogram summary ancestral
```

A protein model is the same command with a different `--model`:

```bash
# proteins under LG, 300 residues per gene
zombi2 sequences seqs/ --from out/ --model lg --length 300 --seed 1
```

Because a protein model has no parameters, passing one is an error rather than a flag that gets quietly ignored: `--model lg --kappa 2.0` stops with *"these options don't apply to --model lg: --kappa"*.
