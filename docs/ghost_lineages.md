# Ghost (dead) lineages — design

ZOMBI2 samples the species tree **backward** as a *reconstructed* birth–death tree: only
lineages with sampled extant descendants survive. That is efficient and conditions cleanly
on `N` tips, but it drops every **dead lineage** — species that went extinct, and (under
incomplete sampling) lineages alive today but not sampled. Those matter for gene families:
in reality a gene can be transferred from a lineage that later vanished — a **transfer from
the dead** (Szöllősi, Tannier, Lartillot & Daubin 2013). Without ghosts, every simulated
transfer is between sampled co-existing lineages, which biases transfer-inference benchmarks.

This page describes a principled way to put the dead lineages back **without abandoning the
backward tree**: *un-prune* it.

## The construction

The reconstructed tree is the complete birth–death tree with the dead branches pruned off.
Un-pruning inverts that, using the exact conditional law of the complete tree given the
reconstructed one (Nee, May & Harvey 1994; Stadler 2009; Lambert & Stadler 2013):

> **Along each edge of the reconstructed tree, dead lineages attach as an inhomogeneous
> Poisson process with intensity `λ(t)·E(t)`**, where `E(t)` is the probability that a lineage
> present at time `t` leaves **no sampled descendant** at the present.

Intuition: a surviving lineage still speciates at rate `λ`; each side-branch independently
leaves a sampled descendant (prob `1−E`) or nothing (prob `E`). The reconstructed tree kept
only the branchings where *both* sides survived — those are its nodes. So *between* its nodes,
every branching had a dying sibling, and those occur at rate `λ·E(t)`.

`E(t)` is exactly the survival quantity ZOMBI2 already solves: the ODE
`dE/dτ = μ − (λ+μ)E + λE²` with `E(0) = 1−ρ` (τ = time before present). Constant-rate
`BirthDeath` has the closed form; `EpisodicBirthDeath` integrates it on a grid.

### Two kinds of dead, one quantity

* **Extinct** — arose and died before the present (needs `μ > 0`).
* **Unsampled extant** — alive today but not sampled (needs `ρ < 1`).

Both are just "leaves no *sampled* descendant," unified by `E(t)` through its `E(0)=1−ρ`
boundary condition.

### Growing each ghost

An attachment at time `t` roots a birth–death subtree grown forward to the present,
**conditioned on leaving no sampled descendant** (the event of probability `E(t)`):

* **Rejection** (`method="rejection"`, default — simple, provably correct): grow a normal BD
  subtree from `t`; sample each present-day tip with prob `ρ`; reject and retry if any tip is
  sampled. Ghosts are small (conditioned to die out), and the expected total work is `O(λ ·
  tree length)`, so it stays cheap. A size guard rejects runaway supercritical attempts early
  (negligible bias).
* **Direct h-transform** (`method="htransform"` — rejection-free): the conditioned process is
  itself a birth–death with per-lineage birth `λ·E(τ)` and death `μ/E(τ)`, so each subtree is
  drawn in one pass. The death rate diverges as `E→0` (i.e. `ρ=1`, `τ→0`), which is exactly
  what forces extinction before the present; we sample the next-event age by inverting the
  cumulative hazard `T(τ)=∫_τ^A g` tabulated on a grid, handling the log-singular first cell
  analytically (`g ≈ 1/τ` there). Statistically equivalent to rejection (verified), and faster
  when heavy extinction/incomplete sampling would make rejection retry a lot.

### Why it is exact

Pruning the ghosts back off returns the original reconstructed tree unchanged (ghosts have no
sampled descendants), so the sampled-tree statistics are untouched — this is an exact
augmentation, not an approximation. Sanity check: Yule with complete sampling (`μ=0, ρ=1`)
gives `E≡0` → zero ghosts → reconstructed = complete.

## Fit with ZOMBI2

The architecture already supports this as a single additive step:

* **New tree step** `add_ghost_lineages(tree, model, *, seed)`: walk each reconstructed edge,
  draw attachments from the `λ(t)E(t)` Poisson process, and graft a conditioned BD subtree at
  each (inserting a binary junction that splits the edge). Ghost tips get `is_extant=False`.
* **Forward gene sim — unchanged**: ghosts are just branches. `branches_alive_at(t)` includes
  them, speciation-cloning gives them a genome, and `leaf_genomes` records only extant tips.
  So transfers can now draw **ghost donors/recipients** automatically.
* **Reconciliation — unchanged**: gene trees already prune to sampled tips, so a gene from a
  ghost donor surfaces as a transfer from a lineage absent from the species tree.

## API

```python
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.5), n_tips=20, age=5.0, seed=1)
z.add_ghost_lineages(tree, z.BirthDeath(1.0, 0.5), seed=1)   # augments in place
genomes = z.simulate_genomes(tree, duplication=0.1, transfer=0.2, loss=0.15,
                             origination=0.5, seed=42)        # transfers can use ghosts

# episodic rates + incomplete sampling (adds unsampled-extant ghosts as well as extinct ones):
m = z.EpisodicBirthDeath(birth=[1.0, 1.6], death=[0.3, 0.8], shifts=[3.0], sampling_fraction=0.6)
tree = z.simulate_species_tree(m, n_tips=50, age=6.0, seed=1)
z.add_ghost_lineages(tree, m, seed=2)
```

## Scope & roadmap

* **Supported:** constant-rate `BirthDeath`/`Yule`, and `EpisodicBirthDeath` with time-varying
  `λ(t)`/`μ(t)` (`E(t)` from the model's ODE grid) and incomplete sampling `ρ<1`. With `ρ=1`
  ghosts are extinct-before-present; with `ρ<1` they also include lineages alive today but
  unsampled. Subtrees are grown by rejection (default) or the direct **h-transform**
  (`method="htransform"`) — both statistically equivalent and exhaustively sanity-checked.
* **Next:** none outstanding for ghosts.

## Sanity checks (tests)

1. **Yule / no extinction** → zero ghosts; the tree is unchanged.
2. **Pruning invariant** → pruning the augmented tree back to its extant leaves reproduces the
   original reconstructed tree (identical Newick).
3. **Extant tips untouched** → the set and times of sampled leaves are identical before/after.
4. **Density scaling** → mean ghost count rises with `μ` (and is 0 at `μ=0`).
5. **Reproducibility** → a fixed seed gives an identical augmented tree.
