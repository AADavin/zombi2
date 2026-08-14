# Introduction

## Why simulate

Evolutionary biology infers the past from what survives into the present: a gene tree from an alignment, a rate of gene loss from a set of genomes, an ancestral body size from the sizes of living species. The true history is gone, so there is nothing to check the answer against.

Simulation is the way around this. You create a dataset whose history you already know: which lineages went extinct, which gene was transferred and when, what the sequence at each internal node was. Run a method on that dataset and you can measure how much of the history it recovers. This is how phylogenetic methods are tested, calibrated and compared.

## What ZOMBI2 is

ZOMBI2 simulates four levels of evolution: **Species**, **Genomes**, **Sequences** and **Traits**. You can run one level, several in sequence, or let one drive another. It is a Python library and a command-line tool over the same engine, taking the same parameters, so a run can be written either way.

## What it can do

Some questions ZOMBI2 is built to answer:

- **How well does a reconciliation method recover the truth?** Evolve gene families under duplication, transfer and loss, and you get every family's true gene tree with the event behind every node. Reconcile against the species tree and score what the method found against what happened.

- **Can a transfer be detected when the donor is gone?** Transfers can come from lineages that later go extinct, so a gene arrives in a survivor from a donor that leaves no other trace. The event log names that donor, so you can ask how often a detection method finds a transfer whose donor is no longer on the tree.

- **What does genome reduction look like in host-restricted bacteria?** Evolve a lifestyle trait, free-living or host-restricted, and let it drive the loss rate. Lineages that move inside a host shed genes faster, so you can measure how much of the genome-size pattern a method attributes to lifestyle rather than to shared ancestry.

- **Can Bayesian methods of inference (relaxed molecular clocks) recover the true ages of a tree?** Give the substitution rate a relaxed clock, so lineages evolve at different paces, and compare the dates a method infers from the alignment against the true node ages.

- **How accurate is ancestral sequence reconstruction?** A run records the sequence at every internal node, not just the tips, so a reconstruction can be compared residue by residue against the one that really sat there.

- **Is a trait correlation real, or an artefact of the tree?** Two traits evolving on the same tree look correlated at the tips whether or not either drives the other, because they share ancestry. Simulate with the correlation switched off, on the same tree, and you have the baseline any comparative method has to beat.

- **Does a trait actually drive diversification?** Let a trait set the speciation rate, so the tree and the trait grow together, and test whether a state-dependent method recovers the effect.

- **What signal survives in gene order?** At the ordered resolution, inversions, transpositions and translocations rearrange genes on chromosomes, and fissions and fusions change the karyotype itself, so synteny and rearrangement methods can be tested against the moves that were actually made.

## Installing it

ZOMBI2 needs Python 3.10 or newer and depends on NumPy and tqdm, plus the `tomli` backport on Python 3.10:

```bash
pip install zombi2
```

`zombi2 --version` confirms the install, and `zombi2 -h` lists the commands: one per level, plus `joint`, which grows a tree and the level driving it together (Chapter 10), and `tools`, analyses that read a finished run (Appendix C).

ZOMBI2 is pure Python over NumPy, with no compiled part to build, and the test suite runs on **Linux,
macOS and Windows** on every change, so the same command does the same thing on all three. Two
things differ on Windows and are worth knowing before they bite:

- **Paths in a rate expression.** A driver path goes inside the rate. On the command line ZOMBI2
  reads the backslashes as written, so `PerCopy(0.25).scaled_by('C:\Users\me\trait_events.tsv',
  {...})` works as pasted. In a Python script the same line is Python source, and Python reads `\U`
  as an escape, so it is a `SyntaxError` — worse, `C:\temp` is read quietly as `C:` followed by a
  tab. Write the path as a raw string there: `r'C:\Users\me\trait_events.tsv'`. Forward slashes work
  in all three places, and Windows accepts them.
- **Paths in a `--params` file.** TOML, not ZOMBI2, reads that file, and TOML's ordinary `"…"` string
  treats a backslash as an escape, so `C:\Users` fails there with a message about a hex value. Put a
  value containing a path in a TOML **literal** string instead, `'''…'''`, which is taken exactly as
  written:

  ```toml
  transfer-to = '''Recipients().weighted_by('C:\Users\me\trait_events.tsv', {'competent': 2.0})'''
  ```

## The examples gallery

Every figure in this book is drawn from a run, and the code behind a great many of them lives in the
[examples gallery](https://aadavin.github.io/zombi2/gallery.html) — a page of worked examples, each
with the exact commands that made it. The gallery numbers them by level: `Sp` species, `Ge` genomes,
`Sq` sequences, `Tr` traits, `Co` conditioning, `Jo` joining. This book cites them the way it cites a
figure, so a paragraph that describes what an inversion does to a chromosome ends **(Ge7)**.

## For the impatient

Two commands: build a species tree, then evolve gene families along it.

```bash
# 1. a species tree
zombi2 species out/

# 2. gene families along it, by duplication, transfer, loss and origination
zombi2 genomes out/
```

`out/` now holds a report of the whole run, and one directory per level:

```
out/run.zombi2                         every file the run wrote, and how to reproduce it
out/species/    species_complete.nwk   the tree, extinct lineages kept
                species_extant.nwk     the survivors only
                species_events.tsv     
                species_fates.tsv      which tips survived, died out, or went unsampled
                species_summary.json   what the run produced, in numbers
                species.log            what was run, and with which parameters
out/genomes/    genome_events.tsv      every duplication, transfer, loss and origination
                genomes.tsv            the genes each species ends up with, ancestors included
                profiles.tsv           gene-family copy numbers, extant species only
                initial_genome.tsv     the genome the run started from
                gene_trees/            one true gene tree per family
                genome_summary.json    what the run produced, in numbers
                genomes.log
```

Open `run.zombi2` first: it names every file the run wrote and carries the commands that rebuild it. Below it each level keeps to its own directory, run log included, and outputs that run to one file per gene family get a directory of their own inside it, so a few hundred families stay legible. `--flat` writes every file straight into `out/` instead, and writes no run report; `--quiet` turns off the progress bar.

Here is the same run from Python. The rates are written out because the defaults above are the command line's: in Python every genome rate starts at 0, so a call that names none of them still runs, and gives a history with no duplication, transfer or loss — just the families it started with, copied down the tree.

```python
from zombi2 import species, genomes

sp  = species.simulate_species_tree(birth=1.0, n_extant=20, seed=1)
gen = genomes.simulate_genomes_family(sp, duplication=0.2, transfer=0.1,
                                      loss=0.25, origination=0.5, seed=42)

sp.write("out/species/")
gen.write("out/genomes/")
```
