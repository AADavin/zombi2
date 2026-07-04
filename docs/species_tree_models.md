# Species-tree models: a roadmap for ZOMBI2

ZOMBI2 currently simulates species trees **backward in time** as a *constant-rate
reconstructed birth–death* process conditioned on the number of extant tips: it draws the
internal-node ages i.i.d. from the reconstructed-process CDF and assembles a ranked tree by
uniform coalescence (Hartmann, Wong & Stadler 2010; Stadler 2009). It also runs the process
**forward** for the complete tree (extinct lineages retained), which unlocks the whole family of
models that need retained extinct lineages or lineage-heterogeneous rates. Shipped today:
constant-rate `BirthDeath`/`Yule`, the episodic (skyline) `EpisodicBirthDeath`, the fossilized
birth–death extensions, explicit **mass-extinction pulses**, per-lineage **`ClaDS`**,
**`DiversityDependent`** (density-dependent) diversification, and **`CladeShiftBirthDeath`**
(scheduled clade-specific rate shifts) (all detailed below).

This note proposes a prioritized menu of additional models. Candidates are judged on
**scientific value** and, crucially, on **how well they fit the backward sampler**: models
that keep the "draw i.i.d. node ages → assemble" shape (possibly with a numerically
inverted CDF) are cheap; models that need a *forward* simulation retaining extinct lineages
(fossils, sampled ancestors) or a type-/trait-structured process (SSE-family) need new
machinery and are deferred.

## Summary

| Model | Fwd/Bwd | Conditions on | Fit to ZOMBI2 backward sampler | Priority |
|---|---|---|---|---|
| Constant-rate BD / Yule | backward | N tips (+age) | **shipped** | — |
| **Episodic (skyline) BD-sampling** | backward | N tips, time, epochs | piecewise-constant CDF, numerically inverted per epoch | **shipped** |
| **Explicit mass extinctions (survival pulses)** | forward | age + pulse (time, fraction) | **shipped** — instantaneous tree-wide kill, layered on any BD/episodic background | — |
| **Time-varying sampling (BDST)** | backward | N sampled tips | same machinery as episodic, ψ(t) | **High** |
| **Diversity-dependent (logistic) BD** | forward | age or N tips (≤K) | **shipped** — λ(n)=λ₀(1−n/K), exact-Gillespie forward (trivial once the loop knows n) | — |
| **Birth–death-shift / ClaDS (per-lineage rates)** | forward | age or N tips | **shipped** — each lineage carries its own λ; daughters jump lognormally at speciation | — |
| **Clade-specific rate shifts** | forward | age + shifts (time, λ, μ) | **shipped** — at a scheduled age a random lineage + its descendants adopt new (λ, μ) | — |
| Multi-type BD (BiSSE/MTBD) | forward | N tips, tip states | needs latent trait history → forward/augmented | Low |
| Continuous trait-dependent (QuaSSE) | forward | N tips, traits | joint tree+trait → very hard | Low |
| **Fossilized birth–death (FBD)** + sampled ancestors | **forward** | extant + fossils | **requires retained extinct lineages** | Low now / High post-ghosts |
| Skyline FBD | **forward** | tips, epochs, fossils | FBD + episodic | Low now / High post-ghosts |
| Occurrence BD (OBDP) | **forward** | mixed data types | forward + occurrence sampling | Low |

\* diversity-dependence *can* be done backward but only with rejection/importance sampling.

## Near-term (fit the current backward sampler)

**1. Episodic / skyline birth–death-sampling — High.** Piecewise-constant `λ(t), μ(t)`
(and sampling `ψ(t)`) over user-defined epochs. Captures mass extinctions and changing
diversification regimes. *Fit:* within each epoch the reconstructed-process CDF keeps its
closed form; sample node ages by numerically inverting the piecewise CDF (a generalization
of the critical-rate branch we already special-case). This is the recommended first
addition and the foundation for several others. (Stadler & Bonhoeffer 2013; Höhna et al.
TESS 2016.)

**2. Time-varying sampling (BDST) — High.** As above but the knob that varies is the
sampling rate `ψ(t)` and we condition on *sampled* tips. Important for phylodynamics
(epidemic trees) and paleontology (variable preservation). *Fit:* reuses the episodic
machinery. (Stadler et al. 2013.)

**3. Diversity-dependent diversification — shipped (forward).** Rates depend on current lineage
count, `λ(n) = λ₀·max(0, 1 − n/K)` (ecological carrying capacity). The CDF becomes path-dependent
so pure i.i.d. *backward* sampling breaks — but *forward* the loop already knows `n`, so it is a
one-line rate function (see [below](#per-lineage-and-diversity-dependent-rates)). Nice conceptual
parallel to our per-family `carrying_capacity` for genes.

**4. Clade-specific rate shifts — shipped (forward).** A finite set of scheduled shifts: at a
given age a lineage (and its descendants) adopts a new (λ, μ) regime — a radiating or collapsing
clade. Forward-natural (the shifted lineage just carries new rates); see
[below](#per-lineage-and-diversity-dependent-rates).

## Deferred — need forward simulation with extinct lineages

These are scientifically valuable but fundamentally require **retaining extinct lineages**
(the "ghost lineage" extension already on ZOMBI2's roadmap), because sampled fossils are
extinct-lineage tips:

- **Fossilized birth–death (FBD)** and **sampled-ancestor** trees (Gavryushkina et al.
  2014; Heath et al. 2014) — the standard model for fossil-calibrated dating.
- **Skyline FBD** (episodic FBD) and the **occurrence birth–death process** (Andréoletti
  et al. 2022, Stadler group) for combined-evidence macroevolution/epidemiology.

Recommendation: build the forward simulator with extinct-lineage retention first (it also
unlocks ghost-lineage HGT for gene families), then FBD becomes a natural addition.

## Deferred — need type/trait-structured processes

The state-dependent speciation–extinction (SSE) family — **BiSSE/MuSSE/MTBD**, **QuaSSE**,
**ClaDS** — requires jointly sampling the tree and a latent trait/rate history. This is a
different machinery from our backward sampler; for these we recommend users lean on
diversitree, RevBayes, or BEAST2 unless there is strong demand to integrate them natively.

## How a new model plugs in

A species-tree model in ZOMBI2 is an object with a `sample_internal_age(u, A)` method (its
inverse-CDF sampler); `simulate_species_tree` supplies conditioning (N tips, age) and does
the ranked assembly, which is model-agnostic. So an episodic model is essentially a new
`sample_internal_age` that inverts a piecewise CDF numerically — the assembly, I/O and the
whole downstream gene-family machinery are untouched.

## Forward simulation (implemented)

`simulate_species_tree` runs *backward* and returns the reconstructed tree (survivors only).
`simulate_species_tree(..., direction="forward")` runs the birth–death process *forward* and returns the
**complete** tree — extinct lineages included natively (`is_extant=False` leaves at their death
times). It is the second route to a complete tree, alongside grafting ghosts onto a backward
tree with [`add_ghost_lineages`](ghost_lineages.md); pass either to `simulate_genomes` and
transfers use the dead lineages automatically.

```python
# grow for a fixed crown age (number of extant tips is random):
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.4), age=5.0, direction="forward", seed=1)

# ...or grow until N extant lineages coexist (age is random):
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.5), n_tips=50, direction="forward", seed=1)

recon = z.prune(tree)   # the reconstructed (survivors-only) counterpart
```

Conventions match the backward crown tree: rooted at the crown (`time == 0`), present at
`total_age`, `age` = crown age; conditioned on ≥2 sampled survivors. Verified against theory
(Yule: mean extant ≈ `2·e^{λ·age}`).

`EpisodicBirthDeath` is supported in **age mode** (time-varying λ/μ and incomplete sampling
`ρ<1`, which marks extant-but-unsampled tips `is_extant=False`). In age mode the present is
fixed at `age`, so the model's ages-before-present map to tree-time `age − t`; a recent
mass-extinction epoch sharply thins the extant tips (e.g. mean extant 100 → 12 when the last
epoch's μ jumps). `n_tips` mode is constant-rate only (the present must be fixed for episodic
rates).

### Explicit mass extinctions (instantaneous survival pulses)

Raising μ over an epoch spreads extra extinction *smoothly* through a time window. A **mass
extinction** in the palaeobiological sense is instead an *instantaneous, tree-wide pulse* — at a
single instant a large fraction of the standing diversity is wiped out at once (a bolide, a
Snowball Earth). Both `BirthDeath`/`Yule` and `EpisodicBirthDeath` take a `mass_extinctions`
list of `(age, fraction)` pulses: at each `age` before the present, every lineage then alive
independently dies with probability `fraction` (equivalently survives with probability
`1 − fraction`). This is the standard mass-extinction formulation of TreeSim (`sim.rateshift.taxa`)
and TESS's explicit mass-extinction birth–death (Stadler 2011; Höhna et al. 2016), here layered
on top of whatever background diversification model you use.

```python
# a constant-rate radiation punctuated by two cataclysms (75% then 50% die):
m = z.BirthDeath(1.0, 0.3, mass_extinctions=[(1.0, 0.75), (2.5, 0.5)])
tree = z.simulate_species_tree(m, age=5.0, direction="forward", seed=1)

# pulses compose with an episodic background, too:
m = z.EpisodicBirthDeath([1.0, 1.4], [0.2, 0.3], [2.0], mass_extinctions=[(1.0, 0.8)])
```

Because a pulse's time is an age before the present, mass extinctions require **age mode** (a
fixed present), like the episodic models; `n_tips` mode and the backward reconstructed sampler
reject them (the killed lineages are not part of the survivors-only tree). The victims become
ordinary extinct leaves at the pulse instant, so `simulate_genomes` treats them as ghost transfer
partners automatically — a mass extinction therefore leaves a genomic signature (families lost
with the dead clades, transfers *from* the dead) with no extra wiring. Killing is per-lineage
Bernoulli, so the realized fraction fluctuates around `fraction` (exactly it in expectation); a
pulse of `fraction=1.0` wipes the tree out entirely and is rejected by the ≥2-survivor
conditioning.

CLI: `zombi2 species --mode forward --age 5 --mass-extinction 1.0 0.75 --mass-extinction 2.5 0.5`.

### Fossilized birth–death (serial / through-time sampling)

The same `BirthDeath` (or `EpisodicBirthDeath`) model gains serial sampling through optional
kwargs: beyond speciation (λ) and extinction (μ), lineages are sampled *through time* at rate ψ
(`fossilization`) — each a **dated fossil tip** — and extant lineages are sampled at the present
with probability ρ (`sampling_fraction`). Sampling removes the lineage by default (`removal=1`),
so every sample is a terminal tip and the tree stays binary (the gene-family machinery is
unaffected).

```python
m = z.BirthDeath(birth=1.0, death=0.5, fossilization=0.5, sampling_fraction=0.9)
tree = z.simulate_species_tree(m, age=6.0, direction="forward", seed=1)   # complete tree + fossils
fbd = z.prune(tree, keep="sampled")   # the sampled tree: dated fossil tips + sampled extant tips
```

Fossil tips carry `sampled=True, is_extant=False` at their (past) sampling times; sampled extant
tips carry `sampled=True, is_extant=True`. `prune(..., keep="sampled")` extracts the FBD sampled tree
(the dated-tip tree used in total-evidence dating); `prune` still gives the extant-only
reconstructed tree. Verified: fossil count scales with ψ (0 at ψ=0), fossils are dated before the
present, and the sampled tree has one tip per sample.

**Sampled ancestors.** With `removal=r<1`, a sampled lineage *continues* with probability `1−r`
instead of being removed — a **sampled ancestor** (the SA-FBD model), represented as a degree-two
node (`sampled=True`, one child). `prune(..., keep="sampled")` keeps these as degree-two nodes; the gene
simulator passes genomes straight through them (they are not gene events), so DTL simulation runs
unchanged on SA trees.

**Episodic FBD.** `EpisodicBirthDeath(birth[], death[], shifts[], *, fossilization=[...],
sampling_fraction=…, removal=…)` composes time-varying λ/μ/ψ with fossil sampling (age mode, like
the other episodic models). So a mass-extinction epoch and a changing fossilization rate can be
combined in one forward run.

### Per-lineage and diversity-dependent rates

The models above are *lineage-homogeneous*: at any instant every lineage shares the same rates
(possibly varying in time). Three further shipped models break that — each lineage (or the tree as
a whole) carries rates that vary with the branching itself. All are **forward-only** and, because
their rates are constant *between* events, are grown by an **exact-Gillespie** loop (draw the next
event time from the summed rate — no thinning bound needed), which also carries mass extinctions
and incomplete sampling.

**ClaDS — per-lineage rate shifts** (`ClaDS`; Maliet, Hartig & Morlon 2019). Every lineage has its
own speciation rate; at each speciation the two daughters inherit the parent's rate times an
independent lognormal jump, `λ_child = λ_parent · exp(N(log α, σ²))`. `α` is the trend (`α<1` =
the empirically typical slow-down of speciation toward the present), `σ` the jump spread, and
extinction is set by a constant turnover `ε = μ/λ` (`ε=0` = ClaDS0, pure birth with shifts). This
generates the heavy among-lineage rate variation ClaDS was designed to capture, and it is the
forward cousin of the relaxed-clock `RateVariation` ZOMBI2 already applies to gene rates.

```python
m = z.ClaDS(lambda_0=1.0, alpha=0.9, sigma=0.2, turnover=0.1)
tree = z.simulate_species_tree(m, age=5.0, direction="forward", seed=1)   # or n_tips=…
```

**Diversity-dependent (density-dependent) birth–death** (`DiversityDependent`; Rabosky & Lovette
2008; Etienne et al. 2012). The speciation rate declines as the tree fills its carrying capacity,
`λ(n) = max(0, λ₀·(1 − n/K))`, with constant `μ`. The tree radiates fast when small and saturates
near `K` (with `μ=0`) or near the equilibrium `n* = K·(1 − μ/λ₀)` — a diversity-brake, the
macroevolutionary analogue of the per-family `carrying_capacity` ZOMBI2 offers for genes.

```python
m = z.DiversityDependent(lambda_0=2.0, death=0.2, carrying_capacity=50)
tree = z.simulate_species_tree(m, age=15.0, direction="forward", seed=1)   # or n_tips ≤ K
```

Both support `age` **and** `n_tips` mode (their rates don't reference age-before-present), reject
the backward sampler (no closed-form reconstructed CDF), and accept `sampling_fraction` and
`mass_extinctions` overlays. For `DiversityDependent`, `n_tips` must be `≤ K`. CLI: `zombi2 species
--mode forward --diversification clads|diversity-dependent …`.

**Clade-specific rate shifts** (`CladeShiftBirthDeath`). Where ClaDS shifts *every* lineage a
little at *every* speciation, this shifts *one* clade a lot at a *scheduled* time — the discrete,
hand-specified version of rate heterogeneity. Diversification runs at the background `(birth,
death)` until, at each scheduled age before the present, a uniformly chosen lineage then alive (and
all of its descendants) adopts a new `(birth, death)` regime. It's the model for "a key innovation
sparks a radiation in one clade" or "a clade enters a slow-down," and the direct forward analogue
of the node/time-specific rate shifts on the roadmap (you can't name an unborn clade in a forward
run, so the shifted lineage is drawn at random — contemporaneous lineages are exchangeable).

```python
# a slow background; at age 3 one clade starts diversifying fast
m = z.CladeShiftBirthDeath(0.6, 0.4, clade_shifts=[(3.0, 2.0, 0.1)])
tree = z.simulate_species_tree(m, age=5.0, direction="forward", seed=1)
```

Because the shift schedule is in ages before the present, this model is **age mode only** (unlike
ClaDS/DD); it too rejects the backward sampler and accepts `sampling_fraction`/`mass_extinctions`.
CLI: `zombi2 species --mode forward --age 5 --clade-shift AGE BIRTH DEATH` (repeatable).

## Key references

- Stadler (2009), *J. Theor. Biol.* — reconstructed birth–death process.
- Hartmann, Wong & Stadler (2010), *Syst. Biol.* — sampling trees / backward sampler.
- Stadler & Bonhoeffer (2013), *Phil. Trans. R. Soc. B* — birth–death skyline.
- Höhna et al. (2016), *Bioinformatics* — TESS (episodic BD simulation/inference).
- Gavryushkina et al. (2014), *PLoS Comput. Biol.* — sampled-ancestor FBD.
- Maliet, Hartig & Morlon (2019), *Nat. Ecol. Evol.* — ClaDS (per-lineage rate shifts).
- Rabosky & Lovette (2008), *Proc. R. Soc. B*; Etienne et al. (2012), *Proc. R. Soc. B* —
  diversity-dependent diversification.
- Andréoletti et al. (2022), *Syst. Biol.* — occurrence birth–death process.
- Louca & Pennell (2020), *Nature* — identifiability limits of time-varying BD (a caveat
  for how far to push time-varying models).
- Louca (2020), *Bioinformatics* — castor, large-scale tree simulation.
