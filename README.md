<h1>
  ZOMBI2&nbsp;<img src="https://raw.githubusercontent.com/AADavin/zombi2/main/assets/logo.svg" alt="ZOMBI2 logo" height="45" align="absmiddle">
</h1>

**[🌐 Website](https://aadavin.github.io/zombi2/)** · [Gallery](https://aadavin.github.io/zombi2/gallery.html) · [Documentation](https://aadavin.github.io/zombi2/docs/) · [Manual (pdf)](https://aadavin.github.io/zombi2/zombi2-manual.pdf)

[![CI](https://github.com/AADavin/zombi2/actions/workflows/ci.yml/badge.svg?branch=main&event=push)](https://github.com/AADavin/zombi2/actions/workflows/ci.yml?query=branch%3Amain)
[![PyPI](https://img.shields.io/pypi/v/zombi2)](https://pypi.org/project/zombi2/)
[![Docs](https://img.shields.io/badge/docs-mkdocs--material-blue)](https://aadavin.github.io/zombi2/docs/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)

**Simulating the evolution of species, genomes, sequences and traits.**

ZOMBI2 simulates evolution at **four levels**: the **species** tree of lineages, the
**genomes** that evolve along it, the **sequences** inside each gene, and the **traits** a
lineage carries. Use it to generate benchmark datasets with known ground truth for phylogenetic and comparative methods.

---

## Install

```bash
pip install zombi2
```

---

## Quickstart

```bash
zombi2 species   out/ --birth 1 --death 0.3 --n-extant 20 --seed 1
zombi2 genomes   out/ --duplication 0.2 --transfer 0.1 --loss 0.25 --seed 1
zombi2 sequences out/ --model hky85 --length 1000 --divergence 0.2 --seed 1
```

From Python, each level is one function, and the result object carries the history:

```python
from zombi2 import species, genomes, sequences
from zombi2.sequences.substitution_models import hky85

sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
g  = genomes.simulate_genomes_family(sp, duplication=0.2, transfer=0.1, loss=0.25, seed=1)
s  = sequences.simulate_sequences(g, model=hky85(), length=1000, divergence=0.2, seed=1)

g.gene_trees                    # the true gene tree of every family

sp.write("out/")                # the trees, the event log and the fates
g.write("out/")                 # gene trees, the event log, profiles
s.write("out/")                 # alignments and phylograms (ancestral is opt-in)
```

---

## Levels

ZOMBI2 is organized around **four levels of evolution**. A genome, a sequence or a trait always
evolves along a species tree, so you run whichever you need, composed into one seeded,
[reproducible](https://aadavin.github.io/zombi2/docs/reproducibility/) run.

<p align="center">
  <img alt="One simulated dataset at all four levels: a species tree with its extinct lineages, the gene order of every surviving genome with homologues linked, the alignment behind one gene family, and two traits drifting together" src="https://raw.githubusercontent.com/AADavin/zombi2/main/assets/overview.png" width="900">
</p>

- **[Species trees](https://aadavin.github.io/zombi2/docs/guide/species-trees/)** — a
  birth–death process with rates that can shift in time, saturate with diversity or drift down
  the tree, plus mass extinctions, incomplete sampling and fossils.
- **[Genomes](https://aadavin.github.io/zombi2/docs/guide/genomes/)** — gene families under
  duplication, transfer, loss and origination, at three resolutions: gene families,
  [ordered](https://aadavin.github.io/zombi2/docs/guide/genomes-ordered/) chromosomes with
  rearrangements, and
  [nucleotide](https://aadavin.github.io/zombi2/docs/guide/genomes-nucleotide/) genomes, with real
  DNA along each chromosome.
- **[Sequences](https://aadavin.github.io/zombi2/docs/guide/sequences/)** — nucleotide (JC69,
  K80, HKY85, GTR) and protein substitution models run down each gene tree, with ancestral
  sequences at every node.
- **[Traits](https://aadavin.github.io/zombi2/docs/guide/traits/)** — continuous traits that
  diffuse, revert to an optimum or shift at speciation, and discrete traits switching between
  states.

<p align="center">
  <img alt="The four levels of evolution ZOMBI2 simulates: the species tree forks into genomes and traits, and sequences continue below genomes" src="https://raw.githubusercontent.com/AADavin/zombi2/main/manual/book/figures/fig-2-1-four-levels.svg" width="360">
</p>

## Conditioning

**[Conditioning](https://aadavin.github.io/zombi2/docs/guide/conditioning/)** is a run reading a value
that has already been grown. There are four parts: the **driver**, the value that is read; the
**target**, what the factor is attached to (a rate, an extent, or which lineage receives a transfer);
the **verb**, `scaled_by`, which joins them; and the **mapping** it carries, which says what each
value of the driver becomes.

<p align="center">
  <img alt="Conditioning: a habitat trait on the left, an arrow labelled drives running right to the gene loss rate and carrying a multiplier for each habitat state, and under the loss rate the expression you write on it, a per-copy loss rate of 0.25 scaled by habitat" src="https://raw.githubusercontent.com/AADavin/zombi2/main/manual/book/figures/conditioning.svg" width="560">
</p>

```bash
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1
zombi2 traits  out/ --kind discrete --states aquatic,terrestrial --switch 0.4 --seed 1
zombi2 genomes out/ --loss "PerCopy(0.25).scaled_by('out/traits/trait_events.tsv', {'aquatic': 4.0})" --seed 1
```

## Joining

**[Joining](https://aadavin.github.io/zombi2/docs/guide/joining/)** is what to reach for when neither
level can be grown first, because each drives the other: one run grows both. A trait
that speeds up speciation is the standard case: the tree shapes the trait's history and the trait
shapes the tree, so the tree is an *output* of the joint run rather than an input to it.

<p align="center">
  <img alt="Joining: body size drives the speciation rate, with the expression written on that rate shown beneath it, the rate creates the tree, and an arrow runs back from the tree to body size because the two grow at the same time" src="https://raw.githubusercontent.com/AADavin/zombi2/main/manual/book/figures/joining.svg" width="700">
</p>

```bash
zombi2 joint out/ --birth "PerLineage(1.0).scaled_by('trait', {'small': 1.0, 'large': 3.0})" \
    --states small,large --switch 0.3 --n-extant 100 --seed 1
```

---

## Performance

A species tree of a million leaves takes a few seconds. On the same gene-family task, ZOMBI2 runs
about **183× faster than the legacy ZOMBI v1** — both pure Python.

<p align="center">
  <img alt="ZOMBI2 performance overview: (a) species-tree simulation scaling to millions of tips; (b) genome simulation at the family, ordered and nucleotide resolutions; (c) ZOMBI2 about 183 times faster than the legacy ZOMBI v1 on one shared 1,000-tip tree" src="https://raw.githubusercontent.com/AADavin/zombi2/main/assets/performance-overview.svg" width="840">
</p>

## Gallery

The [gallery](https://aadavin.github.io/zombi2/gallery.html) is a page of worked examples, one figure
each, with the code that produced it: [species trees](https://aadavin.github.io/zombi2/gallery.html#species),
[genomes](https://aadavin.github.io/zombi2/gallery.html#genomes),
[sequences](https://aadavin.github.io/zombi2/gallery.html#sequences),
[traits](https://aadavin.github.io/zombi2/gallery.html#traits),
[conditioning](https://aadavin.github.io/zombi2/gallery.html#conditioning) and
[joining](https://aadavin.github.io/zombi2/gallery.html#joining). Every figure is drawn with
[Phylustrator](https://pypi.org/project/phylustrator/).

## Citation

A dedicated ZOMBI2 paper is in preparation. Until then, cite the original
[ZOMBI](https://github.com/AADavin/Zombi).

## License

ZOMBI2 is released under the [MIT License](LICENSE).
