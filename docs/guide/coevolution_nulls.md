# Null models of coevolution

!!! tip "In one line"
    Generate any null with `zombi2 coevolve --couple <edge> --null {neutral,cid,timing}`, or in
    Python with a coupling model's [`.null(kind=...)`](#python-api) method (plus the `CID` factory
    for the `traits:species` case). Every null run also writes a `null_manifest.tsv` recording
    exactly which arrow was cut and how the target's variance was preserved. The rest of this page
    is the *why* and the *how*.

Every coevolution model in ZOMBI2 is a directed arrow `driver → target`: the driver's
state modulates the target's rates. Once you can **simulate** such a coupling, the next
question is always inferential — *given data, is the coupling real?* Answering that needs
a **null**: a matched dataset generated with the arrow **cut**, against which a detector's
false-positive rate can be measured. This page defines one null (in fact a small, closed
set of nulls) for every edge, so that "simulate coupled → simulate the matched null →
calibrate your detector" becomes a first-class workflow rather than something each user
hand-rolls.

## Why a naive null is not enough

The cautionary tale is BiSSE. A trait-dependent speciation model asks *"does the observed
character set λ and μ?"*, and the obvious null is a **constant-rate** birth–death process:
no rate variation at all. That null is inadequate, and famously so (Rabosky & Goldberg
2015): real trees **always** carry diversification-rate heterogeneity, from causes that
have nothing to do with the character being tested. If the only two hypotheses on offer
are *"flat rates"* and *"the trait drives rates"*, then **any** heterogeneity — however it
arose — is attributed to the trait. BiSSE fits reject the constant-rate null on trees where
the trait is demonstrably irrelevant, because the null was never the right one to reject.

The fix (Beaulieu & O'Meara 2016) is to make the null *as flexible as the alternative,
minus the causal link*. Diversification is allowed to vary across the tree — via an
**unobserved** ("hidden") factor — but that variation is **independent of the observed
character**. This is a **CID** (Character-Independent Diversification) model. A CID null
has the same number of free rate parameters as the alternative, so rejecting it is
evidence for *the trait specifically*, not merely for *some* rate variation. `CID-2` uses
two hidden rate classes (matched to BiSSE's complexity); `CID-4` uses four (matched to a
full HiSSE).

The principle generalises verbatim to every ZOMBI2 coupling: **a good null keeps the
target's variance and cuts only the arrow.**

## Three archetypes

Cutting the arrow honestly takes exactly one of three forms, and each edge falls into a
single one, decided by *what the driver is*.

**① Neutral null (`neutral`).** Set the coupling strength to its no-effect value (β = 0,
`theta_present = theta_absent`, λ₀ = λ₁, equal panel rates). The arrow is cut and *nothing
compensates* — the target loses its coupling-induced heterogeneity. This is the naive,
constant-rate-style null. It is the weakest of the three (it is exactly the null
Rabosky & Goldberg warn about), but it is cheap and it is the honest baseline for *"does my
detector fire when there is literally no effect and no confound?"* Every edge supports it.

**② Character-independent null (`cid`).** Re-introduce the *same amount* of target
heterogeneity, but source it from an **unobserved** driver that is uncorrelated with the
observed one. This is the generalised HiSSE/CID: the target genuinely varies, the observed
driver genuinely does not explain it. This is the honest, matched-complexity null and the
headline of this layer. It is defined for the four **state→rate** edges (a driver *state*
modulates a target *rate*).

**③ Decoupled-timing null (`timing`).** For the two **event-driven** edges — change that
happens *at speciation* — there is no driver *state* to hide; the driver is the branching
process itself. The honest null keeps the same saltational *amount* of change but
redistributes it **along branches** instead of concentrating it at nodes. This is the classic
punctuation-vs-gradual / node-density contrast (Bokma 2008; Pagel 1999): same variance,
decoupled from the cladogenesis events.

The variance is matched **analytically, in expectation** — not by measuring the realized jumps
of a particular coupled run. The coupled model implies a per-branch expected number of
cladogenetic events (its speciation nodes) and a per-event variance; their product, spread over
total tree length, gives the anagenetic rate the null uses. So the null is a deterministic
function of the coupling *parameters* alone — the same parameters produce the same null
regardless of which coupled realization it is being compared against, which is what makes a
calibration reproducible. (An empirical, per-realization match would couple the null back to the
very dataset it is meant to be independent of.)

![Coupled model versus its neutral and CID nulls: three schematic trees. Coupled — the trait fills
the fast-diversifying clade, so it looks causal. Neutral — a balanced tree with no fast clade, a
weak test. CID — the *same* fast clade as the coupled tree, but the trait is scattered across fast
and slow clades, so it is the honest test of whether the trait tracks
diversification.](../img/coevolve_null_archetypes.svg)

## The layer, edge by edge

| Edge (driver → target) | The arrow | `neutral` | `cid` | `timing` |
| --- | --- | --- | --- | --- |
| `traits:species` (SSE) | trait state → λ, μ | λ₀=λ₁, μ₀=μ₁ (constant-rate BD) | **HiSSE CID-2/-4**: hidden classes carry λ, μ; equal within each class | — |
| `genes:species` (key innov.) | driver-gene presence → λ, μ | β_spec = β_ext = 0 (genes neutral) | *(free)* drivers hidden; analyse the **neutral overlay genome** | — |
| `genes:traits` (modifier → OU) | gene presence → OU optimum θ | θ_present = θ_absent (plain OU) | *(free)* modifier hidden; observe a **neutral overlay genome** | — |
| `traits:genes` (trait → panel) | trait state → panel loss/gain | effect_loss = 0 (uncoupled panel) | hidden trait drives the panel; observe a **second neutral trait** | — |
| `species:traits` (cladogenetic) | speciation event → trait jump | shift = 0 / jump_sigma2 = 0 (anagenetic only) | — | same jump variance, Poisson **along** branches |
| `species:genes` (clado genome) | speciation event → gene burst | clado_loss = clado_gain = 0 (anagenetic only) | — | same turnover, Poisson **along** branches |

The symmetry is the point: the four **state→rate** edges all take a **HiSSE-shaped** null
(hidden, uncorrelated driver); the two **→at-speciation** edges take a **punctuation-anywhere**
null. That closes the set — there is no edge without a defined honest null.

The three **joint** models compose their edges' nulls: cutting *both* arrows of ClaSSE,
co-diversification, or trait–gene feedback yields the fully decoupled null, and cutting one
arrow isolates the contribution of the other. The CLI (below) treats `--null` on a joint run
as "cut every requested arrow with its archetype".

### Reading each `cid` null concretely

One recipe drives all four, and only the ``traits:species`` case needs bespoke machinery:

> **Let a *hidden* driver shape the target (so its heterogeneity is real), then hand the
> analyst a *neutral channel of the observed type, evolving on the same tree*, and keep the
> hidden driver as ground-truth.**

The neutral channel is not a synthetic decoy — it is a genuine evolutionary process (a gene
family, a trait) with real phylogenetic autocorrelation, simply one with **no causal effect**
on the target. And for two of the four edges it already exists in the toolbox for free.

- **`traits:species` → CID.** The one edge with a native hidden driver: `H` hidden classes,
  each a `BiSSE` with **equal** λ across the two observed states; a `hidden_transition` matrix
  moves lineages between them. The tree gets real fast/slow clades from the hidden class; the
  observed character is spread across them and cannot explain the diversity. This is exactly
  what `CID` (a constrained `HiSSE`) builds, and
  what the figure script `fig_sse_hisse.py` draws — *"the honest null a raw BiSSE fit would
  wrongly read as a trait effect."*
- **`genes:species` → CID (free).** The standard workflow **already** produces the null: the
  drivers shape a genuinely heterogeneous tree, and the **neutral bulk genome** overlaid
  afterward with `zombi2 genomes` is a whole panel of real families that vary across that tree
  *without causing any of it*. The null hands the analyst `{tree + neutral Profiles.tsv}` and
  **withholds** the drivers (kept only as ground-truth). No new simulation — the decoupled
  observations were there all along.
- **`genes:traits` → CID (free).** Same trick: the modifier gene shapes a trait with real
  optimum shifts, then a **neutral genome** is overlaid and presented as the observed gene
  content, with the modifier withheld. The trait varies; the genes the analyst tests do not
  explain it.
- **`traits:genes` → CID (one extra trait).** The only edge whose observable is a *trait*, so
  there is no free genome to reuse: a **hidden** trait drives the panel's retention while a
  **second, independent neutral trait** — one extra `simulate_traits` call on the same tree — is
  presented as the observed trait. The panel carries real, trait-shaped heterogeneity that the
  observed trait cannot account for.

!!! note "The shared-tree imprint is a feature"
    On a rate-heterogeneous tree even a *neutral* gene's presence/absence carries a faint imprint
    of tree shape (bushy clades have short branches → slightly less loss). That is not a flaw in
    the null — it is the confound a trustworthy detector must see through, and it makes a neutral
    genome a *stronger* CID null than an abstract hidden class. A good null contains the
    confound; it does not erase it.

### Reading each `timing` null concretely

![The timing null: two copies of the same tree, each teal tick one unit of change. Coupled — change
happens at each speciation, so the ticks cluster at the nodes and sister tips differ sharply. Timing
null — the same number of ticks spread along the branches, so sisters differ only as much as their
shared branch length allows.](../img/coevolve_null_timing.svg)

- **`species:traits` → timing.** Drop the `Cladogenesis` kernel and add the **same** expected
  amount of change as an anagenetic process spread over branch length. For a continuous trait,
  the coupled model adds variance `jump_sigma2` at each of the `E[n]` speciation events on a
  lineage's path; the null replaces it with Brownian diffusion of rate
  `sigma² = jump_sigma2 · E[n] / L` over path length `L`, giving the same expected tip variance.
  For a discrete trait, `shift` per node over `E[n]` nodes maps to a matched anagenetic `Q`.
  Both are computed from the parameters and the tree's expected node count — no realized run is
  inspected. Sister tips now differ *no more than their shared branch length allows*, so a
  detector keyed on node-concentrated change should not fire.
- **`species:genes` → timing.** Replace the founder-effect burst with matched anagenetic
  `loss`/`origination` — the per-event `cladogenetic_loss`/`cladogenetic_gain` × `E[n]` nodes,
  spread over branch length as a constant rate — so the *same* expected total gene turnover is
  applied gradually. The punctuational signature — *sister tips differ because change was
  injected at their split* — disappears, while the marginal amount of turnover is held fixed.

## Python API

The layer follows the idiom ZOMBI2 already uses for its hand-rolled nulls —
`CorrelatedBinary.independent(...)` in [traits](traits.md) and the Potts `J = 0` null — and
generalises it. Two entry points:

**A `.null(kind=...)` method on a coupling model**, returning a *model* with the arrow cut, so
a coupled run and its matched null share every non-coupling parameter automatically. This covers
the nulls that are a pure reparameterisation: `"neutral"` on **any** SSE/gene/trait model (a
`HiSSE`'s collapses its hidden classes to a constant-rate `BiSSE`; a `QuaSSE`'s is a constant
speciation), `"cid"` on a discrete character — `BiSSE` (returns a `CID`) or `MuSSE` (the k-state
generalisation) — and `"timing"` on the two at-speciation models:

```python
from zombi2.coevolve import BiSSE, simulate_sse

alt  = BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2, q01=0.1, q10=0.1)  # trait drives λ
null = alt.null(kind="cid", n_hidden=2)   # CID-2: same rate spread, cut from the trait

coupled  = simulate_sse(alt,  n_tips=200, seed=1)
baseline = simulate_sse(null, n_tips=200, seed=1)   # matched null dataset
```

`kind` accepts `"neutral"`, `"cid"`, or `"timing"`; asking an edge for an archetype it does
not define (e.g. `timing` on a state→rate edge, or `cid` on a non-binary model) raises with a
message pointing here.

**The three gene/trait `cid` nulls are a *workflow*, not a model method** — their honest null is
not a reparameterised model but "run the coupled model, then observe a **neutral channel on the
same tree** while withholding the driver" (see above). For `genes:species` and `genes:traits`
that channel is the neutral genome you already overlay with [`zombi2 genomes`](../cli.md); for
`traits:genes` it is one extra `simulate_traits` call. The CLI `--null cid` (below) orchestrates
this and records the hidden driver as ground-truth; in Python it is the ordinary `simulate_*`
call plus the neutral overlay.

**A named `CID` factory** for the canonical case, so users do not have to know it is "HiSSE with
equal within-class rates":

```python
from zombi2.coevolve import CID

CID.two(lambda_slow=0.5, lambda_fast=2.0, mu=0.2, switch=0.15)   # CID-2
CID.four([0.4, 0.8, 1.6, 2.4], switch=0.1)                       # CID-4
```

`CID.two(...)` is sugar over a `HiSSE` whose two hidden classes carry the
two diversification rates with the **observed-state rates forced equal** — i.e. it *is* a HiSSE,
constrained to be character-independent. Observed-character transitions default to
`q01 = q10 = 0.1` (override per call).

## Command line

One flag on `coevolve` reuses every edge parameter — pass `--null` to generate the matched null
instead of the coupled model:

```
--null {none,neutral,cid,timing}   default: none — cut the arrow this way
--hidden N                         [cid, traits:species] hidden rate classes (2 = CID-2, 4 = CID-4)
```

```bash
# Coupled: a trait drives diversification (BiSSE)
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/alt

# Matched CID-2 null: same rate spread, no trait effect — feed both to your detector
zombi2 coevolve --couple traits:species --sse-model bisse --null cid --hidden 2 \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/null

# Timing null for the punctuational genome: same turnover, spread along branches
zombi2 coevolve --couple species:genes -t species_tree.nwk --null timing \
    --genome-size 30 --clado-gene-loss 0.15 --clado-gene-gain 3 --seed 2 -o out/null_punct

# CID workflow (genes:species): the drivers shape the tree, a NEUTRAL overlay genome is the
# observed data (Profiles.tsv), the drivers are withheld in drivers_ground_truth.tsv
zombi2 coevolve --couple genes:species --drivers 2 --driver-speciation 1.4 \
    --tips 60 --null cid --genome-size 30 --seed 1 -o out/cid_genes
```

The valid `--null` values depend on the edge (see the table): the state→rate edges accept
`neutral`/`cid`, the at-speciation edges accept `neutral`/`timing`, and `none` (the default)
runs the ordinary coupled model. An invalid combination errors early, pointing back here.

The three gene/trait `cid` workflows overlay a neutral genome (a plain `transfer=1.0`, `loss=0.5`
panel of `--genome-size` families) or, for `traits:genes`, a second independent neutral trait — the
neutral rates are recorded in the manifest; re-run [`zombi2 genomes`](../cli.md) for other rates.

### Provenance

A null run writes the ordinary edge outputs **plus** a `null_manifest.tsv` (and a line in
`coevolve.log`) recording (i) which arrow was cut, (ii) the archetype used, and (iii) how the
target's variance was preserved — for `cid`, the hidden-class rates and switch matrix; for
`timing`, the matched anagenetic rate and the variance it targets. A downstream calibration is
then self-documenting: the null dataset carries a machine-readable record of *exactly what was
decoupled*.

## Related nulls elsewhere in ZOMBI2

This layer unifies the "cut the arrow, keep the variance" idea that already appeared, hand-rolled,
in a couple of places:

- **`CorrelatedBinary.independent(...)`** ([traits](traits.md)) — the null for correlated
  binary-trait evolution, the naming precedent this layer's `.null()` follows.
- **The `effect_loss = 0` limit** of [trait-linked gene families](coevolution.md#trait-conditioned-gene-families-traitsgenes)
  is exactly `TraitGeneCoupling.null("neutral")`.

And `CID` (the `traits:species` `cid` null) is a constrained [`HiSSE`](coevolution.md) — the same
model the figure script `fig_sse_hisse.py` draws as "the honest null" — now with a first-class
constructor so you no longer set λ₀ = λ₁ by hand.

## Validation intent

Each null should ship with an inject-and-fail-to-recover test — the mirror image of the
coevolution [validation](coevolution.md) suite:

- **`cid` nulls** — the *observed* driver is statistically independent of the target's
  heterogeneity, by construction. For `traits:species`, the tip-state fraction is uncorrelated
  with clade size even though clade sizes vary widely (the existing
  `test_hisse_hidden_drives_diversification_not_observed` is the template). For `genes:species`,
  observed-driver prevalence is uncorrelated with the rate the tree actually experienced.
- **`timing` nulls** — the punctuational signature is absent: sister-tip divergence is
  explained by shared branch length alone (no excess at nodes), while the marginal amount of
  change matches the coupled run to within sampling error.
- **`neutral` nulls** — the target reduces to its uncoupled model exactly (a regression guard
  that β = 0 recovers plain birth–death / plain OU / an uncoupled panel).

A second, integration-level check is the actual payoff: run an external detector (e.g. a BiSSE
fit) on many `cid` replicates and confirm its rejection rate sits at the nominal α, not the
inflated rate a constant-rate null would show — reproducing Rabosky & Goldberg's result inside
ZOMBI2, for every edge.

## References

- Rabosky, D. L. & Goldberg, E. E. (2015). Model inadequacy and mistaken inferences of
  trait-dependent speciation and extinction. *Systematic Biology* 64(2): 340–355. (Why the
  constant-rate null inflates false positives — the motivation for this whole layer.)
- Beaulieu, J. M. & O'Meara, B. C. (2016). Detecting hidden diversification shifts in models of
  trait-dependent speciation and extinction. *Systematic Biology* 65(4): 583–601. (HiSSE and the
  CID null.)
- Maddison, W. P. & FitzJohn, R. G. (2015). The unsolved challenge to phylogenetic correlation
  tests for categorical characters. *Systematic Biology* 64(1): 127–136. (Proper nulls for
  trait-correlation tests — the same problem for `traits:genes`/`genes:traits`.)
- Uyeda, J. C., Zenil-Ferguson, R. & Pennell, M. W. (2018). Rethinking phylogenetic comparative
  methods. *Systematic Biology* 67(6): 1091–1109. (Unreplicated bursts and the limits of what a
  single tree can distinguish — why the `timing` nulls matter.)
- Pagel, M. (1994). Detecting correlated evolution on phylogenies. *Proc. R. Soc. B* 255:
  37–45. (The independent-evolution null; the `CorrelatedBinary.independent` precedent.)
- Bokma, F. (2008). Detection of "punctuated equilibrium" by Bayesian estimation. *Journal of
  Evolutionary Biology* 21(5): 1218–1227. (Change at speciation vs. along branches — the
  `timing` contrast.)
- Pagel, M. (1999). Inferring the historical patterns of biological evolution. *Nature* 401:
  877–884. (Speciational / punctuational change.)
