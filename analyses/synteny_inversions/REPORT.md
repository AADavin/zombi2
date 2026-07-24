# Inferring the inversion rate from synteny

**What we infer:** the **inversion rate** — how often inversions reshuffle gene order — in inversions
per gene per Myr. **How:** match a small set of synteny *observables* between real yeast genomes and a
ZOMBI2 forward simulation run down a **dated** species tree, using the **nucleotide** genome model, so
inversions act at real base-pair coordinates.

## The idea

An **inversion** reverses a chromosomal segment and flips the strand of the genes in it, so gene order
slowly diverges between species. The *amount* of divergence between two genomes, plotted against *how
long ago they split*, carries the inversion rate. The ingredient that makes this a **rate** rather than
a bare count is a **dated tree**: it turns "how much gene order has been scrambled" into "inversions
per gene per Myr". We summarise gene order with a few **observables**, simulate genome evolution under a
proposed inversion rate down the real dated tree, and accept the rate whose observables match the data
(Approximate Bayesian Computation). In rearrangement vocabulary: ABC with the **breakpoint distance**
(equivalently, the gene-order conservation below) and the **conserved segment length** as summary
statistics.

## Data and the timescale window

Two well-scaled clades measure the rate at two depths inside the usable window:

| Clade | Genomes | Crown age | Role |
|---|---|---|---|
| *Lachancea* | GRYC/Génolevures, chromosome-level + strand (9 species) | 76 Myr | rate |
| *Kluyveromyces* + *Eremothecium* | NCBI RefSeq (6 species) | 94 Myr (pairs to 188 Myr) | rate (corroboration) |

Gene order comes from the annotations. Of the ~5,000 protein-coding genes in each genome, the
single-copy genes present in **every** species of the clade (1:1 orthologs, mmseqs2) form the core we
analyse — **3,306 genes** in *Lachancea* and **2,593** in *Kluyveromyces*+*Eremothecium*. Every pairwise
comparison is over that core. Both dated trees are pruned from the time-calibrated phylogeny of Shen et
al. (2018) and rescaled to Myr — which is where the *rate* units come from.

**Why these clades — the timescale window.** Inversions are only readable from synteny in a window of
roughly **20–200 Myr**: younger, too few inversions have happened (no signal); older, local gene order
is already scrambled (nothing left to read). Both clades sit comfortably inside the window, where both
observables fall off measurably with divergence time — the decay the rate is read from (Figure 1).

## The observables

On the core single-copy gene order of each genome, for every pair of genomes we compute two observables:

| Observable | What it measures |
|---|---|
| **Gene-order conservation** — fraction of neighbouring gene-pairs still shared, vs divergence time | how fast synteny erodes → the rate |
| **Conserved block size** — mean length of the runs that stay together, vs divergence time | breakpoint spacing (corroborates conservation) |

Conservation is the workhorse. A reversal is a **two-break** operation — an inversion cuts exactly two
adjacencies, at its two ends — so the rate at which conservation decays with time *is* the inversion
rate. Conserved block size is the same signal seen the other way round: the colinear runs *between*
those breaks; it corroborates conservation rather than adding independent information. Both are standard
objects: gene-order conservation is the complement of the **breakpoint distance**, and conserved block
size is the **conserved segment length** (Caprara & Lancia 2000). We use the breakpoint distance rather
than the reversal distance because its expectation under a rate is analytically clean (it gives the
exp(−2·r·t) moment estimate below), it is cheap, and it is *event-agnostic* — a translocation breaks
adjacencies too, so it stays the right currency once more event types are added.

![Figure 1](figures/observables.png)

*Figure 1. Gene-order conservation (top) and conserved block size (bottom, log scale) versus divergence
time, one column per clade. Open circles are real genome pairs; filled circles are the best-fit ZOMBI2
simulation. Both decay measurably with time — that decay is what the fit reads, and the simulation at
the fitted rate reproduces it.*

## Fitting the rate — ABC over a grid

We find the inversion rate (and length) under which ZOMBI2 produces genomes whose synteny observables
match the real ones. The exact likelihood of turning one gene order into another sums over *every*
inversion history connecting them and has no closed form (York, Durrett & Nielsen 2002, 2007;
Hannenhalli & Pevzner 1999 solve only the *minimum* distance, which saturates once genomes diverge). So
we use **ABC**: propose a (rate, length), simulate genome evolution down the dated tree — a genome at
every tip, as in the data — reduce to the same pairwise observables, and keep the proposals whose
observables land closest to the real ones. Because ZOMBI2 simulates inversions directly and the model
has only two parameters, we **sweep a grid** over (inversion rate × inversion length): every cell is one
simulation, and mapping the *whole* grid — not just the best cell — is what shows *what the data can and
cannot pin down* (Figure 2).

Two anchors frame the answer before the sweep: Fischer et al. (2006) measured yeast inversion rates of
~3×10⁻⁴–2×10⁻³ per gene·Myr (setting the grid's range), and the observed conservation decays roughly as
exp(−2·r·t), so r ≈ −ln(conservation)/(2·t) gives a first-pass moment estimate.

**The nucleotide model.** Unlike the earlier version of this recipe — which used an ordered genome of
*dimensionless* gene tokens — the simulation here is ZOMBI2's **nucleotide** genome: genes have real
base-pair spans on linear chromosomes, and an inversion reverses an arc at nucleotide coordinates, so a
breakpoint that falls inside a gene cuts it, exactly as a real inversion would. Inversion length is
therefore a **physical span**, converted to a gene count only for reporting. No DNA is simulated — only
the gene layout descends the tree — so a run stays fast.

## Results

**With the inversion length fixed at a few genes (small, as yeast inversions are — Keogh et al. 2000),
the inversion rate is ≈ 3–5×10⁻⁴ per gene·Myr in both clades** — *Lachancea* **2.7×10⁻⁴** and
*Kluyveromyces*+*Eremothecium* **4.6×10⁻⁴** (read at a length of 4 genes; the length itself is not
identifiable, so it is fixed, not fitted). Both sit inside the range of Fischer et al. (2006)
(3×10⁻⁴–2×10⁻³) and reproduce the earlier ordered-model result (*Lachancea* 2.7×10⁻⁴ there too),
confirming that making the genome nucleotide-explicit does not move the rate — the quantity the data
actually constrain.

Figure 2 maps the misfit over the whole (inversion rate × inversion length) grid. The result is a
**diagonal ridge** of equally good fits, because rate and length trade off: more, shorter inversions
scramble gene order about as much as fewer, longer ones (a reversal breaks two adjacencies whatever its
size). The ridge is **narrow across** — fix a length and the rate is pinned (Figure 3, left) — but runs
the **full height** — any length fits, with its own matching rate (Figure 3, right). So gene order
constrains the *rate*, not the *length*.

The short-length edge is where the nucleotide model departs from the dimensionless one, informatively.
In an ordered model of gene tokens a single-gene inversion merely flips a strand without separating
neighbours, so it cannot erode gene order and the fit collapses below two genes. Here an inversion has a
physical span: a ~1-gene-length inversion is a few kilobases, so it can still reverse a gene *pair* or
cut a gene at a breakpoint, and thus still erodes order. That softening is real but clade-dependent — it
shows up as a residual low-misfit patch at length 1 in the milder *Lachancea* clade (Figure 3, right,
where the *Lachancea* curve does not rise at length 1), while the more-eroded *Kluyveromyces* clade still
shows the collapse. It does not touch the reported rate, which is read at a few genes.

![Figure 2](figures/ridge.png)

*Figure 2. The fit to the real data over the (inversion rate × inversion length) grid, one panel per
clade — darker is a closer fit, the white contour outlines the best-fitting region, the star is the best
cell. The bright band is a rate–length trade-off: narrow in rate (fixing a length pins the rate) but
spanning all lengths (gene order does not constrain the inversion size).*

![Figure 3](figures/constraints.png)

*Figure 3. The identifiability statement. Left: misfit versus inversion rate (best over length) has a
clear minimum — the rate is constrained. Right: misfit versus inversion length (best over rate) is flat
above a few genes — the length is not. Both clades agree.*

> **Match the parameter to the data.** Gene order tells you the *rate* of rearrangement, not the *size*
> of individual events — a reversal breaks two adjacencies whatever its length, so the two are averaged
> away. Measuring event size needs a different view (a reversal-distance reconstruction or a sequence
> alignment). Choosing a recipe is partly knowing which parameter each kind of data can constrain —
> synteny is the right tool for the rate, not the size. This is the same lesson as the
> [RED recipe](../red/REPORT.md): a summary pins one thing cleanly and hands what it cannot pin to a
> stated modelling choice (here, the length, fixed from the literature).

## Assumptions and limitations

- **Empirical inputs, simulated process.** The two inputs are real — the dated tree (Shen et al. 2018)
  and the observed gene orders (GRYC/Génolevures, NCBI) — and only the rearrangement process is
  simulated; the inferred rate is the one that makes the simulated genomes consistent with those inputs.
- **Inversions-only model.** Yeast gene orders also see occasional translocations (Fischer et al. 2006;
  Keogh et al. 2000), but small inversions dominate the micro-synteny signal, so treating all local
  rearrangement as inversions makes the rate at most a slight upper bound. Because the summary statistic
  is the (event-agnostic) breakpoint distance, the natural extension is to add a translocation rate as a
  second axis — the multi-event regime where ABC earns its keep. ZOMBI2's nucleotide model supports
  cross-chromosome translocations directly, so this is one more axis in the sweep.
- **Inversion length is fixed, not fitted.** Gene order constrains the rate, not event size (Figure 2),
  so we set the length from the literature (Keogh et al. 2000). The rate barely moves with the assumed
  length above a few genes.
- **Length is now physical.** The earlier ordered-model version of this recipe noted that its genes were
  dimensionless tokens, so "inversion length" was a gene count with no physical span — and predicted that
  a nucleotide-level model would make length physical while leaving the rate unchanged. This port uses
  that nucleotide model and confirms the prediction: length is a base-pair span (breakpoints can even cut
  genes), and the fitted rate is the same.
- **Single-copy core only**, so gene gain, loss and duplication are factored out. The core is
  density-matched: the simulation evolves the full ~5,000-gene genome and observes an evenly-spread
  subset of the real core's size, so the fitted rate is a genome-wide per-gene rate.

## Reproducing this recipe

```bash
cd analyses/synteny_inversions
python fit.py        # grid sweep over (inversion rate x length) for both clades -> results/*.npz,*.json
python figures.py    # Figures 1-3 from the results
```

The real gene orders (`data/<clade>/observed_signed_order.json`, single-copy core, chromosome + strand)
and dated trees (`data/<clade>/tree_myr.nwk`, Myr) ship with the recipe; regenerating them from the raw
annotations (GRYC/Génolevures, NCBI) needs the mmseqs2 ortholog pipeline, which is out of scope here. The
simulation is ZOMBI2's `zombi2.genomes.simulate_genomes_nucleotide` (inversions only), on the real
karyotype (8 / 6 linear chromosomes). Every number regenerates from fixed seeds.

## References

- Caprara A, Lancia G (2000). *Experimental and statistical analysis of sorting by reversals.* In:
  Sankoff D, Nadeau JH (eds), *Comparative Genomics.* Springer, 171–183.
- Fischer G, Rocha EPC, Brunet F, Vergassola M, Dujon B (2006). *Highly variable rates of genome
  rearrangements between hemiascomycetous yeast lineages.* PLoS Genetics 2:e32.
- Hannenhalli S, Pevzner PA (1999). *Transforming cabbage into turnip: polynomial algorithm for sorting
  signed permutations by reversals.* Journal of the ACM 46(1):1–27.
- Keogh RS, Seoighe C, Wolfe KH (2000). *Prevalence of small inversions in yeast gene order evolution.*
  Yeast 16:1009–1020.
- Kurtzman CP (2003). *Phylogenetic circumscription of Saccharomyces, Kluyveromyces and other members of
  the Saccharomycetaceae…* FEMS Yeast Research 4(3):233–245.
- Shen X-X, Opulente DA, Kominek J, et al. (2018). *Tempo and mode of genome evolution in the budding
  yeast subphylum.* Cell 175:1533–1545.
- York TL, Durrett R, Nielsen R (2002). *Bayesian estimation of the number of inversions in the history
  of two chromosomes.* Journal of Computational Biology 9(6):805–818.
- York TL, Durrett R, Nielsen R (2007). *Dependence of paracentric inversion rate on tract length.* BMC
  Bioinformatics 8:115.
