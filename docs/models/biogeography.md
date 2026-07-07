# Biogeography (DEC)

Biogeography models evolve a lineage's **geographic range** — a set of discrete areas — down a
species tree. ZOMBI2 ships the classic **DEC** process (Ree & Smith 2008): along a branch the range
changes anagenetically (a lineage *disperses* into new areas and goes *locally extinct* in old
ones), and at every speciation the ancestral range is split between the two daughters by
**cladogenesis**. The anagenetic part is an [`Mk`](../guide/traits.md#the-mk-model) chain over the
enumerated ranges, so DEC sits alongside the other discrete-trait models but adds the area-set state
space and the cladogenetic split at nodes. Use it to simulate biogeographic histories on a dated
tree, or to make ground truth for testing biogeographic inference (Lagrange, BioGeoBEARS).

| Model | Range dynamics | Reach for it when |
| --- | --- | --- |
| **DEC** | dispersal + local extinction along branches, cladogenetic split at nodes | you want ancestral ranges and range shifts on a dated tree |

## The models

### DEC (dispersal–extinction–cladogenesis)

The state space is every non-empty subset of `areas` with at most `max_range_size` areas. Along a
branch a range in area set `R` gains an absent area `a` (dispersal) at rate `Σ_{b∈R} dispersal[b,a]`
and, when `|R| ≥ 2`, drops an area `a` (local extinction) at rate `extinction[a]` — a range never
becomes empty. `dispersal` is a scalar or an `n×n` matrix; `extinction` is a scalar or a length-`n`
vector; both default to `0.1`. At each speciation `cladogenesis` splits the range: a single-area
range is inherited whole by both daughters (narrow sympatry), while a widespread range yields one
single-area daughter and another that is either the full range (subset sympatry) or its complement
(vicariance), drawn uniformly over the `2·|R|` outcomes (Ree & Smith 2008).

## Command line

`--model dec` on the `trait` command evolves a range down an existing tree; areas are given as a
count or as comma-separated labels.

```bash
# a dated tree to evolve the range along
zombi2 species --birth 1 --death 0.3 --tips 30 --age 5 --seed 1 -o run/

# DEC over 3 areas, starting the root in area 0
zombi2 trait -t run/species_tree.nwk --model dec --areas 3 \
    --dispersal 0.2 --extinction 0.1 --root-range 0 --seed 1 -o run/bio/
```

`--areas` takes a count (`3`) or labels (`A,B,C`); `--max-range-size N` caps how many areas a range
may span (default: all); `--root-range` gives the root range as comma-separated labels (default:
random).

## Python

Models live in `zombi2.traits` (and re-export at the top level, so `zombi2.DEC` also works):

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import DEC, simulate_biogeography

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=5.0, seed=1)
model = DEC(areas=["A", "B", "C"], dispersal=0.2, extinction=0.1)
result = simulate_biogeography(tree, model, root_state={"A"}, seed=1)

result.labeled_values()          # tip ranges, e.g. {leaf: ("A", "C")}
result.ancestral_states()        # internal-node ranges
```

`root_state` is an `Mk` state index or an iterable of area labels (e.g. `{"A"}`); omit it to follow
the model's `root` policy.

## Output

A discrete [`TraitResult`](../guide/traits.md): `labeled_values()` gives each tip's range as a tuple
of area labels, `ancestral_states()` the internal-node ranges, and `.history` the anagenetic
(dispersal/extinction) map along each branch. The CLI writes `traits.tsv` (one `node`/`trait` row
per node, ranges as area sets like `{0,2}`), `trait_tree.nwk`, and the `trait.log` manifest. Nodes
use the [standard node names](../contributing/conventions.md#naming).

## Validation

- **DEC.** The rate matrix has the exact analytic entries — gaining an area at the summed dispersal
  rate, losing one at the extinction rate, and no empty range
  (`test_biogeography.py::test_dispersal_and_extinction_rates`); pure anagenetic evolution down a
  branch matches the transition matrix `expm(Q·t)` over 20 000 replicates, tying the simulator back
  to that rate matrix (`test_biogeography.py::test_dec_anagenesis_matches_transition_matrix`); a
  widespread range splits into its `2·|R|` subset-sympatry / vicariance outcomes with equal
  probability and nothing else (`test_biogeography.py::test_dec_cladogenesis_probabilities`); and a
  fixed `seed` reproduces every node's range
  (`test_biogeography.py::test_biogeography_reproducible`).

## References

- Ree, R. H. & Smith, S. A. (2008). Maximum likelihood inference of geographic range evolution by
  dispersal, local extinction, and cladogenesis. *Systematic Biology* 57(1): 4–14.
