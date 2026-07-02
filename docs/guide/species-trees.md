# Species trees

ZOMBI2 simulates the species tree **backward in time** as a *reconstructed birth–death
process* conditioned on the number of extant tips. Concretely, it draws the internal-node
ages i.i.d. from the reconstructed-process CDF and assembles a ranked tree by uniform
coalescence (Hartmann, Wong & Stadler 2010).

## Models

```python
import zombi2 as z

z.BirthDeath(birth=1.0, death=0.3)   # speciation λ, extinction μ
z.Yule(birth=1.0)                    # pure birth == BirthDeath(birth, death=0)
```

## Simulating

```python
tree = z.simulate_species_tree(
    z.BirthDeath(1.0, 0.3),
    n_tips=20,          # condition on the number of extant species (>= 2)
    age=5.0,            # tree age
    age_type="crown",   # "crown": age of the root; "stem": time of origin
    seed=1,             # or rng=<numpy Generator>
)
```

- **`n_tips`** — the tree has exactly this many extant leaves.
- **`age` / `age_type`** — with `"crown"` the root sits at time 0 and every extant leaf at
  `age`; with `"stem"` the age is the origin time and a stem precedes the crown. (v1
  requires an explicit `age`; conditioning on `n_tips` alone is on the roadmap.)

## The `Tree` object

```python
tree.to_newick()          # timed Newick (branch lengths from node times)
tree.leaves()             # extant leaves
tree.internal_nodes()
tree.branches_alive_at(t) # lineages crossing time t (used by the gene-family loop)
tree.total_age
```

Node times increase forward from the root (time 0) to the extant leaves (`total_age`); a
branch is identified by its child node and spans `(parent.time, node.time]`.

## Episodic (skyline) birth–death

`EpisodicBirthDeath` lets speciation and extinction rates be **piecewise-constant through
time** — the model behind mass extinctions and shifting diversification regimes. Rates are
given one value per epoch (ordered from the present backward), with the epoch boundaries as
strictly increasing **ages** before the present:

```python
# a mass extinction: normal extinction recently, a spike older than age 1
epi = z.EpisodicBirthDeath(birth=[1.0, 1.0], death=[0.2, 3.0], shifts=[1.0])
tree = z.simulate_species_tree(epi, n_tips=30, age=4.0, seed=1)
```

- `birth[i]`, `death[i]` apply to epoch `i`; `shifts` has one fewer entry than `birth`
  (the boundaries). One epoch (`shifts=[]`) reproduces the constant-rate `BirthDeath`.
- The reconstructed tree is still a coalescent point process, so ZOMBI2 samples i.i.d.
  node ages from the (numerically inverted) CDF and assembles exactly as before — the
  tree stays ultrametric.

### Incomplete extant sampling

Pass `sampling_fraction=ρ` (probability an extant species is sampled):

```python
z.EpisodicBirthDeath(birth=[1.0], death=[0.3], shifts=[], sampling_fraction=0.25)
```

!!! note "Scope"
    This covers episodic *diversification* and incomplete *extant* sampling — both keep
    the tree ultrametric. Serial sampling *through time* (dated tips / fossils, as in the
    fossilized birth–death process) needs forward simulation with retained extinct
    lineages and is on the [roadmap](../species_tree_models.md), not implemented yet.

See the [`episodic_species_trees.ipynb`](https://github.com/AADavin/zombi2/blob/main/examples/episodic_species_trees.ipynb)
notebook for worked examples.

## What's coming

Further models — diversity-dependent diversification, node-specific rate shifts, and
(via forward simulation) fossilized birth–death — are laid out in the
[species-tree roadmap](../species_tree_models.md). A new model is just an object with a
`sample_internal_age(u, A)` inverse-CDF sampler; assembly and everything downstream is
model-agnostic.
