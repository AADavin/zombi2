# Sequences

The sequence level does two main things:

* It rescales the gene trees and the species tree from time into substitutions per site (**phylograms**).
* It evolves the residues that sit inside every gene, so each family ends with an alignment.

![Where a sequence lives. A sequence is not evolved along the species tree: it is evolved along the **gene tree**, which the genome run produced and which runs inside the species tree. The two forks are different events. The first is a speciation, which hands the gene to both daughters; the second, marked with a square, is a duplication, and it happens inside one lineage — which is why one lineage can hold two tips of the same family. A run gives one sequence per node of that gene tree, so the tips come out as an alignment.](figures/sequence_nesting_print.png){width=88%}

## Creating phylograms

A gene tree in ZOMBI2 is by default a **chronogram**, its branch lengths measure time. What a sequence actually accumulates along a branch is not time but a number of *substitutions per site*, and that is time multiplied by an evolutionary rate. Turning one into the other is the whole job of the sequence level. Applying a rate to every branch rescales the tree from time into expected substitutions and yields a phylogram.

Two things therefore have to be chosen: *what* changes (the substitution model, the chemistry of which residue turns into which) and *how fast* it changes along each branch, which is the clock.

## The substitution models

ZOMBI2 implements different standard models of sequence evolution:

The four nucleotide matrices are the standard published ones [@jukes1969evolution; @kimura1980simple; @hasegawa1985dating; @tavare1986some], as are the four protein ones [@dayhoff1978model; @jones1992rapid; @whelan2001general; @le2008improved]. The model decides the alphabet, and `length` counts whatever that alphabet holds: bases for a nucleotide model, residues for a protein one. The nucleotide models are four different rate matrices, not one model with four settings, but they nest in the order written (`jc69` is `k80` with `kappa=1`, and `k80` is `hky85` with equal base frequencies) so each step adds free parameters. The protein matrices are **empirical**, each estimated once from a large set of real alignments and then used as a fixed table, which is why they take no parameters.

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

You give $S$ and $\pi$ rather than $Q$ directly, and that restriction is deliberate. ZOMBI2 computes $P(t) = e^{Qt}$ through a symmetric eigendecomposition that is only valid for a **time-reversible** model, one where $\pi_i Q_{ij} = \pi_j Q_{ji}$ for every pair. A symmetric $S$ times $\pi$ satisfies that by construction, so there is no way to write a model here that the engine would evaluate wrongly. A general $Q$ handed straight to `SubstitutionModel` could be, and is refused with an error rather than run: non-reversible models such as UNREST are not implemented, and a wrong transition matrix produces plausible sequences that are not the model you asked for.

There is no command-line flag for a custom matrix. A twenty-state matrix is 190 numbers, which is a file format rather than an argument; `--model` offers the menu, and your own matrix is a Python constructor.

## Rate variation across sites

So far every site of a gene evolves at the same speed, which is a model no real gene obeys. Some positions are held nearly fixed by what the protein has to do; others drift freely. The standard way to say so is a **Gamma distribution of rates across sites**: each site gets a multiplier drawn from a Gamma with mean 1, so one number, its **shape**, sets how unequal the sites are. A small shape means a few fast sites among many slow ones; a large one is nearly flat.

You add it to the model, not to the rate:

```python
gamma_model = hky85(kappa=2.0).across_sites(gamma_shape=0.5)
varied = sequences.simulate_sequences(my_genomes, model=gamma_model, length=1000, seed=1)
```

The Gamma is cut into a small number of equal-probability classes, each represented by its mean [@yang1994variable], four by default and changed with `rate_categories`. Cutting it is not an approximation made for tidiness: a site's rate is what its branch length is computed from, and a continuous draw would give every site its own branch length and so its own transition matrix. With classes, the sites sharing a class share the work.

A second setting adds a class of sites that **never** change:

```python
both = hky85(kappa=2.0).across_sites(gamma_shape=0.5, invariant=0.1)
print(both.name)                  # HKY85+I+G4
sequences.simulate_sequences(my_genomes, model=both, length=1000, seed=1)
```

`invariant=0.1` sets aside a tenth of the sites as unchangeable. Real alignments have columns that are constant because the site cannot change rather than because it happened not to, and a Gamma alone fits those badly. Either works alone, and the model's name records what you chose, `HKY85+I+G4`, which is what the run prints and logs.

From the command line the same three settings are flags, and they apply to any model on the menu:

```bash
zombi2 sequences seqs/ --from out/ --model hky85 --gamma-shape 0.5 --invariant 0.1 --seed 1
```

**Branch lengths do not change.** The classes are normalised so the mean rate over all sites is exactly 1, invariant sites included, so a branch length in the phylograms is still substitutions per site, now the mean over them. A run with rate variation and a run without, at the same rate, have the *same tree*; what differs is how the change is spread across the columns. That is what makes the two comparable, and it is why the mean-1 normalisation is checked rather than assumed.

Two consequences worth knowing. The mean pairwise identity a run reports goes **up** under `+Γ` at the same divergence, because the slow and invariant sites keep their matches while the fast ones saturate. And on a nucleotide run the spacer between genes keeps its own model, which is flat by default: `model` does not reach `intergene_model`, since the spacer's job is to be the unconstrained null.

The field writes these as suffixes on the model's name, and so does ZOMBI2: `gamma_shape=` is `+G` [@yang1994variable], `invariant=` is `+I`, and the two together are `+I+G` [@gu1995maximum].

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
pseudocount.

Where the numbers come from is your business. Two ways.

**From an alignment of the real family.** Column frequencies, with a pseudocount so nothing is
impossible:

```python
def profile_from_alignment(columns, alphabet, pseudocount=0.1):
    """One row per alignment column, from the residues seen in it."""
    counts = np.array([[c.count(a) for a in alphabet] for c in columns], dtype=float)
    counts += pseudocount                       # nothing is impossible, only unlikely
    return counts / counts.sum(axis=1, keepdims=True)

columns = ["LLLIL", "GGGGA", "KKRKK"]           # three columns of a five-sequence alignment
print(profile_from_alignment(columns, protein.alphabet).shape)      # (3, 20)
```

**From a protein language model.** These already produce a distribution over amino acids at every
position, which is the table this section wants, so nothing has to be derived from it:

<!-- doc-test: skip — transformers and torch are not ZOMBI2 dependencies -->

```python
import torch
from transformers import T5ForConditionalGeneration, T5Tokenizer

def profile_from_prott5(sequence, alphabet, model_name="Rostlab/prot_t5_xl_uniref50"):
    tok = T5Tokenizer.from_pretrained(model_name)
    mdl = T5ForConditionalGeneration.from_pretrained(model_name).eval()
    ids = tok(" ".join(sequence), return_tensors="pt")
    with torch.no_grad():
        logits = mdl(**ids, decoder_input_ids=ids["input_ids"]).logits[0]
    p = torch.softmax(logits, dim=-1).numpy()
    cols = [tok.convert_tokens_to_ids(a) for a in alphabet]      # the model's own token ids
    prof = p[:len(sequence), cols]
    return prof / prof.sum(axis=1, keepdims=True)
```

**Compute the table once and keep it.** A language model's output depends on its version and on the
hardware it ran on, so calling one at simulation time would put something outside this manual's
reproducibility promise inside your run. Save the array beside your script and load it. The
simulation is then bit-identical from the seed as usual, and a reader can inspect the profile without
owning a GPU.

A **flat** profile, every row the model's own frequencies, is the model without one. Statistically,
not byte for byte: a hundred single-site models walk the random stream differently from one
hundred-site model, so the same seed gives a different draw from the same distribution.

Profiles and `+Γ` are about different things and compose. A profile says **which** amino acids belong
at a site; a Gamma says **how fast** sites change. Decorate the model as usual and you get both.

**An amino-acid profile needs a protein model, so it belongs to a family or ordered run.** A
nucleotide genome is measured in base pairs and its blocks are read on either strand, so that
resolution refuses protein models altogether: there is no complement of an amino acid. Profiles are
still accepted there, but over the four bases: a row per base pair, saying what belongs at that
coordinate. The row count then has to match the block's length in bp, since the genome run already
fixed it.

## Relaxed molecular clocks

Rate variation across sites says which *positions* change fast. A clock says which *lineages* do. The two are orthogonal and compose.

The rate itself is `substitution`, and it is counted **per site**: a gene-tree branch of Δ*t* time accrues `substitution · Δt` substitutions at every site. Leave it alone and it is `1.0` everywhere: the **strict clock**, one tempo for the whole tree.

This is the one number people most often get wrong, so it is worth doing an example in full. Suppose your tree runs 20 time units from the start of the stem to the leaves, and you leave the rate at `1.0`. Then every site accumulates `1.0 × 20 = 20` substitutions on the way from the origin to a tip.

Twenty substitutions per site is a great many, and this is where the intuition usually breaks: those are twenty *events*, not twenty visible differences. A site has only four states, so once it has been hit a couple of times it starts landing back on bases it already held, and the differences you can actually see stop accumulating long before the events do. Past about one substitution per site, two sequences are as different as two random ones, and the alignment no longer records where they came from. The rate has not stopped working; the history has simply been overwritten.

So the rate that suits a run depends on the height of the tree it runs down, which is why no default can be right for every tree. Read it off backwards instead, from the divergence you want:

$$\text{substitution} = \frac{\text{substitutions per site you want, origin to tip}}{\text{height of the tree}}$$

On the 20-unit tree above, sequences around 80% identical need roughly `0.2 / 20 = 0.01`, not `1.0`. The difference is not subtle: simulated down that tree, `0.01` gives tips about 77% identical, while the default `1.0` gives **25%**, precisely the score of two sequences with no shared history at all. The table below is measured the same way, on simulated JC69 alignments. Its first column is rate times height, so it transfers between trees of different heights. The second column does not transfer as cleanly: the identity is a mean over pairs of tips, and how much of the tree's height separates the average pair depends on where its splits sit. A tree whose splits are close to the tips keeps more identity than one whose splits are close to the root, and at the deep end of the table two ordinary trees can differ by twenty points. Read the second column as a guide, not a lookup:

| Substitutions per site, origin to tip | Mean identity between two tips |
|---|---|
| 0.02 | 98% |
| 0.09 | 90% |
| 0.17 | 81% |
| 0.34 | 68% |
| 0.86 | 46% |
| 1.72 | 34% |

The last row is close to the floor: two unrelated DNA sequences already match at 25% by chance, so 34% is almost no signal at all. Every run reports the identity it actually produced in its summary line, and warns when it comes out this close to the floor, so you never have to work it out from the flags alone.

You can also state the divergence and let ZOMBI2 do the division. `divergence` is that first column, substitutions per site from the root to a tip, and the rate is solved for from the height of the tree the run is about to use:

```python
sequences.simulate_sequences(my_genomes, model=hky85(), length=1000, divergence=0.2)
```

```bash
zombi2 sequences out/ --model hky85 --length 1000 --divergence 0.2
```

The two say different things and can be given together. `substitution` says what *kind* of clock, strict or relaxed by a modifier, and `divergence` says how far it drifts, so a relaxed clock calibrated to a divergence is written with the shape alone and the scale beside it:

```python
from zombi2.params import LogNormal, PerSite
# the shape: an uncorrelated clock, and the scale: 0.2 substitutions per site
substitution = PerSite().varying_among('lineages', LogNormal(0.0, 0.3))
divergence   = 0.2
```

Giving `substitution` a base *number* as well is an error rather than an override, because the base is precisely what `divergence` solves for, and a run whose rate came from somewhere other than its own command line is a run you cannot reproduce from it. The resolved rate is written into the run log either way.

A substitution rate that changes from lineage to lineage is what the field calls a **relaxed clock** [@lepage2007general]. It is not a new kind of object here: you chain a verb onto the rate, exactly as at every other level.

```python
from zombi2.params import Drift, Gamma, LogNormal, PerSite

# strict clock: one rate everywhere; the default, so write nothing
substitution = 1.0

# relaxed: each lineage draws its own rate, independently of its neighbours
substitution = PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))           # lognormal
substitution = PerSite(1.0).varying_among('lineages', Gamma(shape=4.0, scale=0.25))  # or any other

# relaxed: each lineage inherits its parent's rate and drifts from it
substitution = PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.3)))
```

The **law** is the second argument, written out: `LogNormal(0.0, sigma)` for the common case, or any
other built-in distribution — and whichever you give, the draw is **normalised to mean 1**, so what a
distribution contributes is its *shape* and the base keeps meaning the average rate.

**A bare distribution** has *no memory*: each lineage is an independent draw, so a lineage's rate tells you nothing about its neighbours'. It is the value each lineage gets: `LogNormal(0.0, sigma)` for the usual case, or any other built-in shape.

**`Drift`** has memory: a daughter starts at its parent's rate and multiplies it by one lognormal step, so close relatives evolve at similar rates. That is the **autocorrelated** clock. Both draws are mean-corrected, so widening the distribution spreads the lineages apart without moving the average rate off the number you typed. Rate variation across sites is not a modifier, and does not belong in the rate at all: it is part of the model, as above.

One important point: **the clock belongs to the species tree, not to the gene trees.**

ZOMBI2 draws one rate for each species branch. Every gene that passes through that branch then evolves at that rate. Each gene-tree branch looks up the species branch it sits inside, which the genome run already recorded. The consequence is that if a species evolves quickly, all of its genes evolve quickly together.

A reference table that can be handy to people who want to implement a specific model from the literature:

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| one rate everywhere | `substitution = 1.0` (default) | Strict / global clock |
| each lineage i.i.d. lognormal | `PerSite(1.0).varying_among('lineages', LogNormal(0.0, …))` | Uncorrelated lognormal (UCLN) |
| each lineage i.i.d. gamma | `PerSite(1.0).varying_among('lineages', Gamma(...))` | Uncorrelated gamma (UGAM) |
| the rate drifts parent to daughter | `PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, …)))` | Autocorrelated lognormal |
| the rate reads another level | `PerSite(1.0).scaled_by(trait, {…})` | Trait-dependent rate of molecular evolution |

### A trait can drive the rate

The two clocks above make a lineage fast or slow at random. A third verb makes it fast or slow for a *reason*: `scaled_by` reads a trait grown first and looks the factor up from that lineage's state.

```python
from zombi2 import traits

habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.3, seed=1)

result = sequences.simulate_sequences(my_genomes, model=hky85(), length=1000, seed=2,
    substitution = PerSite(0.05).scaled_by(habitat, {"cave": 0.5, "surface": 1.0}))
```

Cave lineages now evolve at half the rate of surface ones. The driver is the grown trait, or the path to the `trait_events.tsv` it wrote, the same two spellings every driven rate takes. This is conditioning, so it is two ordinary runs in order, and Chapter 9 covers the whole mechanism.

A clock and a driver **compose**, because verbs chain and their factors multiply. Written together, a lineage's branch length is the base rate, times the tempo it was dealt, times the factor its state gives:

```python
sequences.simulate_sequences(my_genomes, model=hky85(), length=1000, seed=2,
    substitution = PerSite(0.05).varying_among('lineages', LogNormal(0.0, 0.3))
                                .scaled_by(habitat, {"cave": 0.5, "surface": 1.0}))
```

A discrete trait switches partway along a branch, and ZOMBI2 does not read the driver once per branch. It integrates the rate across the branch, breaking at each switch. A lineage that leaves the cave halfway down a branch of length 2 accrues `0.05 × 0.5 × 1` substitutions per site before the move and `0.05 × 1.0 × 1` after it, so the branch is `0.075` long rather than `0.05` or `0.1`. The gene phylograms and the clock species tree carry that same number, so the tree a run writes is the tree its alignments were drawn along.

The reverse direction runs too: `result.gc()` makes a finished run's GC content drive a trait grown after it, or a further sequence run, and `result.composition(letters)` does the same for any letters of the run's alphabet, an amino-acid frequency say (Chapter 9). What the pair cannot be is **joined**, because a sequence lives inside a gene and never feeds back into the trait, so there is nothing for the two to decide together. Naming a live level (`scaled_by("trait", …)`) says so rather than looking for a file. One other limit here: `divergence` is refused alongside a driven rate, because it solves for the base by assuming the modifiers average to 1, which the two clocks are corrected to do and a driver is not. Set the base yourself there.

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

This is what an endosymbiont study needs. Its clade evolves faster *and* drifts toward AT, and the two mislead a tree-builder differently: a fast branch causes long-branch attraction, an AT-rich one causes compositional attraction, where unrelated AT-rich lineages group together because they look alike. Scoping the rate to the same clade gives both at once.

Two limits, and one caveat worth knowing before you read the output.

The alphabet is shared, because one gene copy's sequence is one string. So are the across-site rate classes: a site's class is drawn once for the family and holds all the way down the tree, so two clades cannot sort the same site into different classes. Both are refused with the reason rather than accepted and half applied.

The driver must be a clade. A trait switches partway along a branch, and this level samples one transition matrix per branch, so a model that changed mid-branch would need the branch cut at the switch — that is not built. A trait can still drive the *rate*.

The caveat: every model is normalised to one expected substitution per site per unit branch length, and that holds **at stationarity**. A lineage that has just entered the AT-rich clade is not yet at its frequencies, so while its composition is still relaxing it accrues somewhat fewer substitutions than its branch length claims. That transient is usually the thing being studied, and it is the ordinary price of a model that varies along the tree, not a defect. It shows up plainly if you raise the divergence: the clade's A+T climbs from 0.54 toward 0.80 as it has more time to get there.

`Models` is Python only, and there is no flag for it: a model is a K×K matrix, which has no written form a flag can carry — the same reason `reversible()` has none.

## The objects

`simulate_sequences` returns a **`SequencesResult`**, which carries:

- `.alignments`, the observable data: for each family, the sequence at every **extant** gene copy. This is the alignment a phylogenetic method would be handed.
- `.ancestral`, : internal nodes, and the tips where a copy was lost or its species died. The run wrote a sequence at each node as it went, so these are the exact ancestors, not estimates. With `.alignments` it accounts for every node of the tree exactly once, so every label in a complete phylogram names a sequence.
- `.founding`, for each family, the sequence it began with, at its origination.
- `.phylograms`, for each family, its gene tree with branch lengths converted from time into substitutions per site: the tree the sequences were drawn along.
- `.species_phylogram`, the same conversion applied to the species tree, so the clock is visible as branch lengths.
- `.genomes`, the assembled genome of each **extant** tip, keyed by tip name — the genomes you would have sequenced. `.node_genomes` is the same for every node, ancestors and extinct lineages included, and `.initial_genome` is the run's starting point. All three are present only when the run came from a **nucleotide** genome. See below.

As with every level, the bundle also carries `.seed` and `.write(directory, outputs=[...])` to put the chosen outputs on disk.

### Where a sequence starts

## Large runs

This is the level where a run's memory goes. Every family's alignment and every ancestral sequence are held at once, so what you can run is bounded by families × copies × sites rather than by time. `stream_to` writes each family's files the moment it is finished and keeps nothing, handing back a light handle with a path instead of a `SequencesResult` holding everything. On the command line it is `--stream`.

Memory then stops growing with the sequences. Two hundred species and 640 families at 300 sites costs 466 MB in memory and 329 MB streamed; raise the sites to 1500 and the in-memory run goes to 1008 MB while the streamed one barely moves, to 392 MB. What is left is the genome run being read, which is the floor.

It is a memory choice and not a modelling one: the same seed writes the same files either way, so a streamed run and an in-memory one are the same dataset. `outputs=` picks which files, exactly as `.write` takes them, and it composes with `parallel`. A **nucleotide** run cannot stream: it puts whole genomes back together, which needs every block's sequence at once.

## Running on a nucleotide genome

Hand the level a **nucleotide** genome run and you get whole assembled genomes in FASTA. The genome
can be one ZOMBI2 drew, as below, or a real annotation you supply:

<!-- doc-test: skip — needs an annotation and its FASTA, which the reader supplies -->
```python
my_genomes = genomes.simulate_genomes_nucleotide(
    tree, gff="ecoli.gff", fasta="ecoli.fasta", inversion=1.0, inversion_extent=5000,
    duplication=0.3, loss=0.3, duplication_extent=3000, loss_extent=3000, seed=1)
```

Genes and spacer get their own models. `model` evolves the genes; `intergene_model` evolves the spacer, at `intergene_speed` times the rate (3× by default), under `jc69` by default, which is flat and has no free parameters.

```python
from zombi2 import species, genomes, sequences
from zombi2.sequences.substitution_models import hky85

tree = species.simulate_species_tree(
    birth=1.0, death=0.2, n_extant=5, seed=1).complete_tree
my_genomes = genomes.simulate_genomes_nucleotide(
    tree, root_length=6000, genes=6, gene_length=400,
    inversion=1.0, inversion_extent=500, duplication=1.0, loss=1.0,
    duplication_extent=1200, loss_extent=1200, seed=1)

result = sequences.simulate_sequences(my_genomes, model=hky85(kappa=3.0),
                                      intergene_speed=3.0, substitution=0.05, seed=1)

result.genomes["n5"]             # {chromosome: sequence}, a whole assembled genome
result.node_genomes["n0"]        # the same at an ancestor, reconstructed not estimated
result.initial_genome            # the genome the run started with
```

From the command line it is the same two commands as any other run:

```bash
zombi2 genomes out/ --resolution nucleotide --gff ecoli.gff --trim-overlaps \
  --inversion 5.0 --inversion-extent 50000 --loss 2.0 --loss-extent 8000 --seed 7

zombi2 sequences out/ --model hky85 --kappa 3.0 --substitution 0.02 \
  --intergene-speed 3.0 --seed 7
```

### Starting from a real sequence

So far the founding sequence of each block is *drawn* from the model's frequencies, random ACGT. Hand the genomes run a **FASTA** alongside the GFF and it starts from the sequence you supply instead:

<!-- doc-test: skip — needs an annotation and its FASTA, which the reader supplies -->
```python
my_genomes = genomes.simulate_genomes_nucleotide(
    tree, gff="ecoli.gff", fasta="ecoli.fasta",     # layout AND letters
    inversion=1.0, loss=0.3, loss_extent=3000, seed=1)
result = sequences.simulate_sequences(
    my_genomes, model=hky85(kappa=3.0), substitution=0.05, seed=1)
```

The FASTA has one `>seqid` record per GFF `##sequence-region`, each exactly its declared length. Every block is then founded from the real DNA at its own initial coordinates, so an assembled genome descends from exactly what you gave. A gene that origination invents mid-run has no supplied DNA (it did not exist initially), so its block still draws from the model.

## On the command line

On the command line the genome run is handed over as a **directory**, the run directory itself, which by then holds the genomes. `zombi2 sequences out/` reads that run's species tree and event log and replays the gene genealogy from them, so the two commands chain without anything else passing between them. A nucleotide run given a `--fasta` also hands its initial DNA across (in `initial_sequence.fasta`), so the sequences descend from your real sequence without you naming it twice. Point `--from` at another run to read one and write somewhere else.

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
