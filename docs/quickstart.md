# Quickstart

## A full simulation in ten lines

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import simulate_genomes

# 1. species tree (backward), conditioned on 20 extant tips and a crown age of 5
tree = simulate_species_tree(BirthDeath(birth=1.0, death=0.3),
                             n_tips=20, age=5.0, seed=1)

# 2. gene families (forward) along that tree
genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, initial_families=40, seed=42)

# results
print(tree.to_newick())                 # timed Newick
print(genomes.profiles.matrix.shape)    # (n_families, n_species) copy numbers
genomes.write("out/")                   # write every output to out/
```

Everything is seeded: the same `seed` gives byte-identical results.

## What `genomes.write("out/")` writes

| File / folder | Contents |
|---|---|
| `species_tree.nwk` | the timed species tree (Newick) |
| `species_nodes.tsv` | node name, time, leaf/extant flags |
| `gene_family_events/<fid>_events.tsv` | per-family event log (O/D/T/S/L with lineage ids) |
| `gene_trees/<fid>_complete.nwk` / `_extant.nwk` | reconstructed gene trees (with / without losses) |
| `Transfers.tsv` | every transfer (donor, recipient, ids) |
| `Gene_family_summary.tsv` | per-family event counts and extant copies |
| `Events_trace.tsv` | compact one-row-per-event trace (O/D/T/L), replayable by `zombi2 sequence` |
| `Profiles.tsv` / `Presence.tsv` | families × species copy-number / presence matrix |

The Python `write()` saves **all** of these by default. The `zombi2 genomes` **command line**
writes a subset — `Profiles.tsv`, `Presence.tsv` and `gene_trees/` (its `--write` default is
`profiles trees`); pass `--write all` for the full set above.

## Command line

```bash
# species tree only
zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o out/

# gene families along a supplied tree
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 -o out/
```

The built-in model runs on Rust automatically. `--write` selects the outputs — the default is
`profiles trees`, `--write all` gives the full set above, and `--write profiles` gives just the
counts. See the full [command-line interface](cli.md) for every option.

## Next steps

- [Species trees](guide/species-trees.md) — models and conditioning.
- [Gene families & rates](guide/genomes.md) — shared vs per-family sampled rates.
- [Transfers](guide/genomes.md) — replacement, distance, self-transfer.
- [Bounding growth](guide/genomes.md) — caps and carrying capacity.
- [Gene trees & output](guide/genomes.md#gene-trees-output).
- [Ordered genomes](guide/genomes.md) — gene order, inversions, transpositions.
