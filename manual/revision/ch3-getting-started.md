# Getting started

This chapter installs ZOMBI2, runs a first simulation, and explains what the run leaves on disk.

## Installation

ZOMBI2 needs Python 3.10 or newer. It is written in Python and its only dependency is NumPy, so there is no compiler and no toolchain to set up.

Install the released version from PyPI:

```bash
pip install zombi2
```

Or install from a source checkout:

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2
pip install .
```

Two optional extras cover work on ZOMBI2 itself rather than simulation: `pip install "zombi2[dev]"` adds the test suite's requirements (pytest, scipy, ruff), and `pip install "zombi2[docs]"` adds the documentation build (MkDocs).

Installing puts a `zombi2` command on your path. Check it:

```bash
$ zombi2 --version
ZOMBI2 0.2.0
```

The library reports the same version:

```python
>>> import zombi2
>>> zombi2.__version__
'0.2.0'
```

From a source checkout you can also run the test suite:

```bash
pytest -q
```

`zombi2 -h` lists the commands, one per level. `zombi2 <command> -h` documents one command, its flags, and the rate modifiers that level accepts.

## Your first simulation

Every ZOMBI2 workflow starts with a species tree, because the other levels evolve along it. `zombi2 species` grows one. You give it a speciation rate, an extinction rate, a stopping condition, a seed, and an output directory:

```bash
$ zombi2 species --birth 1.0 --death 0.3 --n-extant 20 --seed 1 -o first_run
wrote first_run/ (20 extant + 9 extinct (29 tips, 57 nodes)) in 0.00031 s
```

The tree grew forward until 20 lineages were alive. Nine other lineages went extinct on the way. They are part of the true history, so ZOMBI2 keeps them: the run writes both the **complete** tree, which carries the extinct lineages, and the **extant** tree, which does not.

Now evolve gene families along that tree with `zombi2 genomes`. The four events are origination, duplication, transfer and loss, each with its own rate:

```bash
$ zombi2 genomes -t first_run/species_complete.nwk \
    --duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5 \
    --initial-families 20 --seed 42 -o first_run
wrote first_run/ (30 gene families across 20 extant genomes (unordered)) in 0.0117 s
```

`--initial-families 20` starts the crown lineage with 20 families; the rest arose by origination during the run. `-t` takes the **complete** tree, not the extant one, because genes evolve along extinct lineages too and a transfer can come from a donor that later dies.

Both commands take `--seed`. The same flags and the same seed give the same output, so a run is reproducible and a bug report is repeatable.

## The output folder

Both commands wrote into `first_run`:

```
$ ls first_run
genome_events.tsv
genome_species_tree.nwk
genomes.log
profiles.tsv
species.log
species_complete.nwk
species_events.tsv
species_extant.nwk
```

| File | What it is |
|---|---|
| `species_complete.nwk` | the whole tree that grew, extinct lineages included |
| `species_extant.nwk` | the survivors only |
| `species_events.tsv` | every speciation and extinction, with its time |
| `genome_events.tsv` | every origination, duplication, transfer and loss — the gene families' true history |
| `profiles.tsv` | each family's copy count in each extant species |
| `genome_species_tree.nwk` | the species tree the genome run used, labelled to match the `lineage` column of `genome_events.tsv` |
| `species.log`, `genomes.log` | one per command: the version, the timestamp, the command line, and every resolved parameter |

Trees are Newick with branch lengths in time. Tables are TSV with a header row. The two event logs record the true history of the run. `profiles.tsv` is the phyletic profile matrix, the copy-count table used in comparative genomics:

```
$ head -4 first_run/profiles.tsv | cut -f1-6
family	n19	n20	n21	n26	n27
0	1	1	2	1	3
1	2	1	1	1	0
2	0	1	2	2	2
```

Each command has a `--write` flag that chooses which outputs to keep, so `--write events` writes the event log alone. Other rates and resolutions add files of their own: `--resolution ordered`, for instance, also writes the gene order of every leaf genome. Appendix B lists every file every level can write.

## Reading the results back in Python

Both levels are also functions. Each returns a result object rather than a bare tree, because a run produces more than one thing:

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_unordered

sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
gen = simulate_genomes_unordered(sp, duplication=0.2, transfer=0.1, loss=0.25,
                                 origination=0.5, initial_families=20, seed=42)
```

Here the tree is handed to the genome level as an object. On the command line it is handed over as a Newick file, which does not carry the root's own branch length, so the same seeds give a different draw and different numbers from the run above.

`sp` is a `SpeciesResult`. It holds the two trees, the event log, and the seed:

```python
>>> sp.n_extant
20
>>> len(sp.complete_tree.extinct())
9
>>> len(sp.events)
37
>>> sp.events[0]
Event(time=0.2818670099229476, kind='speciation', node=0, children=(1, 2))
```

`gen` is a `GenomesResult`. It holds the genome at every node, the event log, the profiles, and the gene trees derived from the log:

```python
>>> gen.profiles.shape          # families × extant species
(25, 20)
>>> gen.profiles.matrix[:3, :5]
array([[1, 2, 0, 0, 1],
       [2, 1, 2, 0, 0],
       [2, 1, 1, 0, 1]])
>>> gen.family_counts(19)       # the genome at node 19, as family → copies
Counter({4: 5, 13: 4, 12: 3, 1: 2, 2: 2, 16: 2, 10: 2, 17: 2, 0: 1, 14: 1, 19: 1, 18: 1, 3: 1})
>>> gen.gene_trees[0].to_newick()[:64]
'(((g566:1.75259,(g691:1.32667,g692:1.32667)duplication_n20:0.425'
```

Every result writes the same way the commands do, so you can work in memory and only then materialise the files:

```python
sp.write("first_run")
gen.write("first_run")
```

To pick a run back up, read its tree with `read_newick`, which returns the tree and a map from node ids to the labels of the original file. The map is empty for a ZOMBI2 tree, whose labels already are the ids:

```python
>>> import pathlib
>>> from zombi2.species import read_newick
>>> tree, names = read_newick(pathlib.Path("first_run/species_complete.nwk").read_text())
>>> len(tree.nodes), len(tree.extant()), len(tree.extinct())
(57, 20, 9)
```

That tree goes straight into the next level, which is exactly what `zombi2 genomes -t` does:

```python
gen = simulate_genomes_unordered(tree, duplication=0.2, loss=0.25,
                                 origination=0.5, seed=7)
```

`read_newick` also accepts trees from elsewhere, so a genome or trait run can start from a published phylogeny instead of a simulated one. The remaining outputs are plain tab-separated tables with a header row; read them with whatever you already use.
