# Inferring the inversion rate from synteny

**What we infer:** the **inversion rate** — how often inversions reshuffle gene order — in
inversions per gene per Myr. **How:** match a small set of synteny *observables* between real
genomes and a ZOMBI2 forward simulation run down a **dated** species tree.

## The idea

Over time, an **inversion** reverses a chromosomal segment (and flips the strand of the genes in
it), so gene order slowly diverges between species. The *amount* of divergence between two genomes,
plotted against *how long ago they split*, carries the inversion rate. The one ingredient that makes
this a **rate** rather than a bare count is a **dated tree**: it turns "how much gene order has been
scrambled" into "inversions per gene per Myr". We summarise gene order with a few **observables**,
simulate genome evolution with ZOMBI2 under a proposed inversion rate down the real dated tree, and
accept the rate whose observables match the data (Approximate Bayesian Computation).

## Data and the timescale window

We use **two well-scaled clades to measure the rate**, at two depths inside the usable window:

| Clade | Genomes | Crown age | Role |
|---|---|---|---|
| *Lachancea* | GRYC/Génolevures, chromosome-level + strand (9 species) | 76 Myr | rate |
| *Kluyveromyces* + *Eremothecium* | NCBI RefSeq (6 species) | 94 Myr (pairs to 188 Myr) | rate (corroboration) |

Gene order comes from the annotations. Of the ~5,000 protein-coding genes in each genome, the ones
that are single-copy and present in **every** species of the clade (1:1 orthologs, mmseqs2) form the
core we analyse — **3,306 genes** in *Lachancea* and **2,593** in *Kluyveromyces*+*Eremothecium*.
Every pairwise comparison below is over that core. Both dated trees are pruned from the
time-calibrated phylogeny of Shen et al. (2018) and rescaled to Myr — which is where the *rate* units
come from. The *Kluyveromyces*+*Eremothecium* clade is the classic pre-WGD synteny group: its deepest
pair, *K. lactis* vs *Ashbya* (*Eremothecium*) *gossypii*, is exactly the one Keogh et al. (2000)
measured.

<figure markdown="span">
  ![](../img/synteny/lachancea_tree.png){ width="80%" }
  <figcaption markdown="span">Figure 1. The dated *Lachancea* tree (9 species, pruned from Shen et al. 2018), time axis in millions of years ago; 76-Myr crown.</figcaption>
</figure>

<figure markdown="span">
  ![](../img/synteny/kle_tree.png){ width="80%" }
  <figcaption markdown="span">Figure 2. The dated *Kluyveromyces* + *Eremothecium* tree (6 species); 94-Myr crown, with cross-genus pairs reaching ~188 Myr, so a stronger inversion signal.</figcaption>
</figure>

**Why these clades — the timescale window.** Inversions are only readable from synteny in a window
of roughly **20–200 Myr**: younger than ~20 Myr, too few inversions have happened (no signal); older
than ~200 Myr, local gene order is already scrambled (nothing left to read). This is also why the
subphylum-wide (400 Myr) tree is *not* used here — its gene order is saturated and its all-species
core is too sparse to see local inversions. Both *Lachancea* and *Kluyveromyces*+*Eremothecium* sit
comfortably inside the window, where both synteny observables fall off measurably with divergence
time — the decay the rate is read from (Figure 3).

<figure markdown="span">
  ![](../img/synteny/observables_real.png){ width="100%" }
  <figcaption markdown="span">Figure 3. The two synteny observables on the real genomes, one column per clade: gene-order pairwise conservation (top) and mean conserved block size (bottom) versus divergence time. Both decay measurably with time — that decay is what the fit reads. Each point is a pair of genomes.</figcaption>
</figure>

## The observables

On the core single-copy gene order of each genome, for every pair of genomes we compute two
observables (Figure 5):

| Observable | What it measures |
|---|---|
| **Gene-order conservation** — fraction of neighbouring gene-pairs still shared, vs divergence time | how fast synteny erodes → the rate |
| **Conserved block size** — mean length of the runs that stay together, vs divergence time | breakpoint spacing (corroborates conservation) |

Conservation is the workhorse. A reversal is a **two-break** operation — an inversion cuts exactly
two adjacencies, at its two ends (Figure 4) — so the rate at which conservation decays with time *is*
the inversion rate. Conserved block size is the *same* signal seen the other way round: the colinear
runs *between* those breaks. It corroborates conservation rather than adding independent information —
a fact that turns out to matter when we ask what the data can and cannot pin down (Results).

<figure markdown="span">
  ![](../img/genome_inversion.svg){ width="72%" }
  <figcaption markdown="span">Figure 4. An inversion reverses a chromosomal segment and flips the strand of every gene inside it — the elementary event whose rate we measure.</figcaption>
</figure>

<figure markdown="span">
  ![](../img/synteny/observables_toy.png){ width="95%" }
  <figcaption markdown="span">Figure 5. The two observables on a toy pair of genomes related by an inversion of the {4,5} segment. **Gene-order conservation** is the fraction of neighbouring gene-pairs that survive (here 3 of 5). **Conserved block sizes** are the lengths of the runs that stay together (3, 2, 1); collected across all genome pairs, their distribution is the second observable.</figcaption>
</figure>

## Fitting the rate

The goal is to find the inversion rate (and length) under which ZOMBI2 produces genomes whose synteny
observables match the real ones. There are two general strategies for this.

**1. Maximum likelihood.** Write the probability of the observed gene orders as a function of the
rate and maximise it. This is the gold standard when it is available — but for inversions the
likelihood is intractable. The probability of turning one gene order into another is a sum over
*every* inversion history that could connect them, and the number of such histories grows
combinatorially with genome size, so it cannot be evaluated. What *is* computable is the **minimum**
number of inversions between two signed gene orders — the inversion distance, solved in polynomial
time by Hannenhalli–Pevzner theory (Hannenhalli & Pevzner 1999) — but the minimum is not the
likelihood: it ignores the many longer histories real evolution also takes, and it saturates once
genomes are well diverged. So a full ML treatment of inversions is not practical here.

**2. Approximate Bayesian Computation (ABC).** When the likelihood is out of reach but *simulating*
the process is easy, ABC replaces the likelihood with the simulator: propose parameter values,
simulate a genome down the dated tree, reduce it to the same observables as the data, and keep the
proposals whose observables land closest to the real ones. Because ZOMBI2 simulates inversions on an
ordered genome directly and cheaply, and the model we need to capture this behaviour has only one or
two parameters to tune, we can simply **sweep a grid** over (inversion rate × inversion length)
instead of sampling a prior: every cell is one simulation, and its distance to the real observables
says how well that (rate, length) fits. Mapping the *whole* grid — rather than reporting only the
best cell — is what lets us read off not just the estimate but *what the data can and cannot pin
down* (Figure 6). This is the strategy we use.

We are **not flying blind**. Two anchors frame the answer before the sweep runs:

1. **Literature.** Fischer et al. (2006) measured yeast inversion rates of ~3×10⁻⁴–2×10⁻³ per
   gene·Myr; this sets the grid's range.
2. **A moment estimate.** The observed conservation decays roughly as exp(−2·r·t), so
   r ≈ −ln(conservation) / (2·t) gives a first-pass rate directly from the data.

The moment estimate gives the ballpark; the grid sweep refines it and, crucially, shows how tightly
each parameter is actually constrained.

## Results

**With the inversion length fixed at a few genes (small, as yeast inversions are — Keogh et al.
2000), the inversion rate is ≈ 3–4×10⁻⁴ per gene·Myr in both clades** — *Lachancea* **2.7×10⁻⁴**
and *Kluyveromyces*+*Eremothecium* **3.9×10⁻⁴** (moment estimates 1.8 and 6.7×10⁻⁴). Both sit inside
the range of Fischer et al. (2006) (3×10⁻⁴–2×10⁻³). That is the number the recipe delivers.

Figure 6 maps the misfit to the real data over the whole (inversion rate × inversion length) grid. The
bright band is the result — a **diagonal ridge** of equally good fits, because rate and length trade
off: more, shorter inversions scramble gene order about as much as fewer, longer ones (a reversal
breaks two adjacencies whatever its size). The ridge is **narrow across** — fix a length and the rate
is pinned — but runs the **full height** — any length fits, with its own matching rate. So gene order
constrains the *rate*, not the *length*. It does fix one edge: single-gene inversions merely flip a
strand without separating neighbours, so they cannot erode gene order — the fit collapses below two
genes.

<figure markdown="span">
  ![](../img/synteny/heatmap.png){ width="100%" }
  <figcaption markdown="span">Figure 6. The fit to the real data over the (inversion rate × inversion length) grid, one panel per clade — brighter (yellow) is a closer fit; the star is the best cell. The bright band is a rate–length trade-off: it is narrow in rate (fixing a length pins the rate) but spans all lengths (gene order does not constrain the inversion size). The two bright ends are the same ridge, not two answers. Reading the rate at a fixed literature length (a few genes) gives ~3–4×10⁻⁴ per gene·Myr.</figcaption>
</figure>

So we **fix the length from the literature** (a few genes) rather than fit it, and read the rate off
the ridge. The dependence is mild above two genes: the rate holds to within a factor of ~2 for lengths
from ~4 to ~200 genes, rising only for the very shortest inversions. The two clades agree to within
their spread, at the low end of the Fischer et al. (2006) range, consistent with these being
genomically stable pre-WGD yeasts. The simulations are chromosome-faithful — real karyotypes of 8 and
7 linear chromosomes — via ZOMBI2's multichromosome ordered genome.

!!! tip "Match the parameter to the data"
    Gene order tells you the *rate* of rearrangement, not the *size* of individual events — a
    reversal breaks two adjacencies whatever its length, so the two are averaged away. Measuring
    event size needs a different view: the events one at a time, from a reversal-distance
    reconstruction or a sequence alignment. Choosing a recipe is partly knowing which parameter each
    kind of data can actually constrain — synteny is the right tool for the rate, not for the size.

## Assumptions and limitations

- **Dated tree and gene order are both real** — the tree from Shen et al. (2018), the annotations
  from GRYC/Génolevures and NCBI.
- **Inversions-only model.** We treat *all* local rearrangement as inversions. Small inversions
  dominate yeast micro-synteny (Fischer et al. 2006; Keogh et al. 2000), so ignoring the occasional
  translocation makes the rate at most a slight upper bound.
- **Inversion length is fixed, not fitted.** Gene order constrains the rate, not event size (a
  reversal breaks two adjacencies whatever its length; Figure 6), so we set the inversion length from
  the literature (Keogh et al. 2000) rather than infer it. The rate barely moves with the assumed
  length above a few genes, so this costs the headline little; measuring the size itself would need a
  reversal-distance reconstruction or a sequence alignment (see the note above).
- **Genes are dimensionless tokens.** The model has gene *order* but no nucleotide lengths or
  intergenic spacing, so "inversion length" is a gene count, not a physical span. A nucleotide-level
  model would make length physical, but the rate — the quantity the data actually constrain — is
  unchanged, so this is not pursued here.
- **Single-copy core only**, so gene gain, loss and duplication are factored out. The core (3,306 /
  2,593 genes) is density-matched: the simulation evolves the full ~5,000-gene genome and observes an
  evenly-spread subset of the same size, so the fitted rate is a genome-wide per-gene rate.

## Reproducing this recipe

The pipeline is four steps:

1. **Dated trees** — prune the time-calibrated phylogeny of Shen et al. (2018) to each clade and
   rescale to Myr.
2. **Real gene order** — take the annotations, keep the single-copy core (1:1 orthologs via mmseqs2),
   and record each gene's chromosome, family, and strand.
3. **Grid sweep** — for each (inversion rate × inversion length) cell, simulate a chromosome-faithful
   genome down the dated tree with ZOMBI2's multichromosome ordered genome, reduce it to the two
   observables, and score its distance to the real data.
4. **Read the rate** — fix the inversion length from the literature and read the best-fitting rate
   off the ridge.

!!! note "Runnable code"
    This page documents the method and the result. The scripts are being ported onto ZOMBI2's
    multichromosome ordered genome (now part of the core engine) and will ship alongside the recipe.

## References

- Fischer G, Rocha EPC, Brunet F, Vergassola M, Dujon B (2006). *Highly variable rates of genome
  rearrangements between hemiascomycetous yeast lineages.* PLoS Genetics 2:e32.
- Hannenhalli S, Pevzner PA (1999). *Transforming cabbage into turnip: polynomial algorithm for
  sorting signed permutations by reversals.* Journal of the ACM 46(1):1–27.
- Keogh RS, Seoighe C, Wolfe KH (2000). *Prevalence of small inversions in yeast gene order
  evolution.* Yeast 16:1009–1020.
- Shen X-X, Opulente DA, Kominek J, et al. (2018). *Tempo and mode of genome evolution in the budding
  yeast subphylum.* Cell 175:1533–1545.
