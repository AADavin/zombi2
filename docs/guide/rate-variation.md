# Rate variation (relaxed clock)

Species trees and gene trees in ZOMBI2 are **timetrees** — branch lengths are time. To get
the branch lengths you would infer from *sequence evolution*, overlay a substitution rate
that varies across the tree. `RateVariation` implements the discrete-bin, Markov-switching
model from the GTDB archaea study, turning a chronogram into a **phylogram**.

## The model

- A set of rate **bins** — multipliers, some fast (> 1), some slow (< 1).
- A continuous-time Markov process runs **along the phylogeny**, switching bins at a
  constant rate; the current bin is inherited by both descendants at each node.
- A branch may thus be split into several **segments** in different bins; its substitution
  length is `Σ (segment_duration × bin_rate)`.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)

rv = z.RateVariation(bins=[0.5, 1.0, 2.0],   # slow / clock-like / fast
                     switch_rate=1.0,        # how often the rate changes along the tree
                     weights=[0.5, 0.3, 0.2])  # stationary probabilities (default uniform)
scaled = rv.scale(tree, seed=1)

print(scaled.to_newick())          # the phylogram (substitution lengths)
scaled.branch_lengths[node]        # substitution length of one branch
scaled.segments[node]              # [(bin_index, duration), ...] pieces of that branch
```

- `switch_rate=0` gives a strict clock: one bin along the whole tree.
- Over a large tree, the average `substitution / time` ratio approaches the stationary mean
  rate `Σ wᵢ·binᵢ`.

## Works on gene trees too

`RateVariation` operates on any `Tree`. Gene trees come out of reconstruction as Newick, so
load them with `read_newick` first:

```python
genomes = z.simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=1)
_, extant = genomes.gene_trees()["1"]
gene_tree = z.read_newick(extant)          # Newick -> Tree
phylogram = z.RateVariation(bins=[0.5, 2.0], switch_rate=1.0).scale(gene_tree, seed=1)
```

The same rate process can be applied independently to the species tree and to each gene
tree, or you could drive them from a shared process — a natural building block for
simulating realistic, non-clocklike branch lengths.
