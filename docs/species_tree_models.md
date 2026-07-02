# Species-tree models: a roadmap for ZOMBI2

ZOMBI2 currently simulates species trees **backward in time** as a *constant-rate
reconstructed birth–death* process conditioned on the number of extant tips: it draws the
internal-node ages i.i.d. from the reconstructed-process CDF and assembles a ranked tree by
uniform coalescence (Hartmann, Wong & Stadler 2010; Stadler 2009). `BirthDeath` and `Yule`
are the two models shipped today.

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
| **Episodic (skyline) BD-sampling** | backward | N tips, time, epochs | piecewise-constant CDF, numerically inverted per epoch | **High** |
| **Time-varying sampling (BDST)** | backward | N sampled tips | same machinery as episodic, ψ(t) | **High** |
| Diversity-dependent (logistic) BD | backward* | N tips | CDF becomes path-dependent (couples to topology) → rejection/sequential sampling | Medium |
| Node/time-specific rate shifts | backward | N tips | extension of episodic; shifts tied to nodes | Medium |
| Birth–death-shift / ClaDS (per-lineage rates) | forward | N tips | needs joint tree+rate history → forward | Low |
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

**3. Diversity-dependent diversification — Medium.** Rates depend on current lineage count,
e.g. `λ(n) = λ₀·max(0, 1 − n/K)` (ecological carrying capacity). *Fit:* the CDF becomes
path-dependent (the rate at a node depends on how many lineages exist there), so pure i.i.d.
sampling breaks; feasible via rejection or sequential/importance sampling on top of the
current framework. Nice conceptual parallel to our per-family `carrying_capacity` for genes.

**4. Node/time-specific rate shifts — Medium.** A finite set of shifts at named times or
nodes (e.g. a K–Pg extinction spike). An extension of the episodic model.

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
`simulate_species_tree_forward` runs the birth–death process *forward* and returns the
**complete** tree — extinct lineages included natively (`is_extant=False` leaves at their death
times). It is the second route to a complete tree, alongside grafting ghosts onto a backward
tree with [`add_ghost_lineages`](ghost_lineages.md); pass either to `simulate_genomes` and
transfers use the dead lineages automatically.

```python
# grow for a fixed crown age (number of extant tips is random):
tree = z.simulate_species_tree_forward(z.BirthDeath(1.0, 0.4), age=5.0, seed=1)

# ...or grow until N extant lineages coexist (age is random):
tree = z.simulate_species_tree_forward(z.BirthDeath(1.0, 0.5), n_tips=50, seed=1)

recon = z.prune_to_extant(tree)   # the reconstructed (survivors-only) counterpart
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

### Fossilized birth–death (serial / through-time sampling)

`FossilizedBirthDeath(birth, death, fossilization, sampling)` adds serial sampling: beyond
speciation (λ) and extinction (μ), lineages are sampled *through time* at rate ψ
(`fossilization`) — each a **dated fossil tip** — and extant lineages are sampled at the present
with probability ρ (`sampling`). Sampling removes the lineage (removal probability 1), so every
sample is a terminal tip and the tree stays binary (the gene-family machinery is unaffected).

```python
m = z.FossilizedBirthDeath(birth=1.0, death=0.5, fossilization=0.5, sampling=0.9)
tree = z.simulate_species_tree_forward(m, age=6.0, seed=1)   # complete tree + fossils
fbd = z.prune_to_sampled(tree)   # the sampled tree: dated fossil tips + sampled extant tips
```

Fossil tips carry `sampled=True, is_extant=False` at their (past) sampling times; sampled extant
tips carry `sampled=True, is_extant=True`. `prune_to_sampled` extracts the FBD sampled tree
(the dated-tip tree used in total-evidence dating); `prune_to_extant` still gives the extant-only
reconstructed tree. Verified: fossil count scales with ψ (0 at ψ=0), fossils are dated before the
present, and the sampled tree has one tip per sample.

**Sampled ancestors.** With `removal=r<1`, a sampled lineage *continues* with probability `1−r`
instead of being removed — a **sampled ancestor** (the SA-FBD model), represented as a degree-two
node (`sampled=True`, one child). `prune_to_sampled` keeps these as degree-two nodes; the gene
simulator passes genomes straight through them (they are not gene events), so DTL simulation runs
unchanged on SA trees.

**Episodic FBD.** `EpisodicFossilizedBirthDeath(birth[], death[], fossilization[], shifts[], *,
sampling, removal)` composes time-varying λ/μ/ψ with fossil sampling (age mode, like the other
episodic models). So a mass-extinction epoch and a changing fossilization rate can be combined in
one forward run.

## Key references

- Stadler (2009), *J. Theor. Biol.* — reconstructed birth–death process.
- Hartmann, Wong & Stadler (2010), *Syst. Biol.* — sampling trees / backward sampler.
- Stadler & Bonhoeffer (2013), *Phil. Trans. R. Soc. B* — birth–death skyline.
- Höhna et al. (2016), *Bioinformatics* — TESS (episodic BD simulation/inference).
- Gavryushkina et al. (2014), *PLoS Comput. Biol.* — sampled-ancestor FBD.
- Andréoletti et al. (2022), *Syst. Biol.* — occurrence birth–death process.
- Louca & Pennell (2020), *Nature* — identifiability limits of time-varying BD (a caveat
  for how far to push time-varying models).
- Louca (2020), *Bioinformatics* — castor, large-scale tree simulation.
