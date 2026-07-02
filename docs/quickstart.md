# Quickstart

## A full simulation in ten lines

```python
import zombi2 as z

# 1. species tree (backward), conditioned on 20 extant tips and a crown age of 5
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3),
                               n_tips=20, age=5.0, seed=1)

# 2. gene families (forward) along that tree
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_size=40, seed=42)

# results
print(tree.to_newick())                 # timed Newick
print(genomes.profiles.matrix.shape)    # (n_families, n_species) copy numbers
genomes.write("out/")                   # write every output to out/
```

Everything is seeded: the same `seed` gives byte-identical results.

## What you get from `write("out/")`

| File / folder | Contents |
|---|---|
| `species_tree.nwk` | the timed species tree (Newick) |
| `species_nodes.tsv` | node name, time, leaf/extant flags |
| `gene_family_events/<fid>_events.tsv` | per-family event log (O/D/T/S/L with lineage ids) |
| `gene_trees/<fid>_complete.nwk` / `_extant.nwk` | reconstructed gene trees (with / without losses) |
| `Transfers.tsv` | every transfer (donor, recipient, ids) |
| `Gene_family_summary.tsv` | per-family event counts and extant copies |
| `Profiles.tsv` / `Presence.tsv` | families × species copy-number / presence matrix |

## Command line

```bash
# species tree only
zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o out/

# species tree + gene families
zombi2 all --birth 1 --death 0.2 --tips 20 --age 5 \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

## Next steps

- [Species trees](guide/species-trees.md) — models and conditioning.
- [Gene families & rates](guide/gene-families.md) — uniform vs per-family sampled rates.
- [Transfers](guide/transfers.md) — replacement, distance, self-transfer.
- [Bounding growth](guide/growth.md) — caps and carrying capacity.
- [Gene trees & output](guide/gene-trees-and-output.md).
- [Ordered genomes](guide/ordered-genomes.md) — gene order, inversions, transpositions.
