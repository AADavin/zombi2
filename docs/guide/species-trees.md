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

## What's coming

The current menu is constant-rate birth–death (and Yule). A prioritized set of further
models — episodic/skyline birth–death, diversity-dependent diversification, fossilized
birth–death — is laid out in the [species-tree roadmap](../species_tree_models.md). A new
model is just an object with a `sample_internal_age(u, A)` inverse-CDF sampler; assembly
and everything downstream is model-agnostic.
