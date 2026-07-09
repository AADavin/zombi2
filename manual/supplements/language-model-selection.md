# Language-model-guided selection

ZOMBI2 evolves protein-coding sequences under a *neutral* substitution model. Over a long branch that
is a problem: a neutral model has no notion of what a protein is *for*, so a gene left to drift long
enough wanders into sequence that no longer folds or functions — statistically a protein, biologically
noise. Real coding sequences are held in place by **selection**: most amino-acid changes are removed,
a few are tolerated, and which is which depends on the site and on the rest of the protein.

This supplement describes an **experimental** ZOMBI2 feature that overlays selection on sequence
evolution using a **protein language model** (a PLM such as ESM2) as an empirical, learned model of
what proteins look like. The language model plays the role of a fitness landscape: sequences it scores
as protein-like are favoured, sequences it scores as unlikely are suppressed. The result is a simulator
whose coding sequences drift *toward realistic proteins* instead of away from them, with a single knob —
`beta` — that turns selection from off (neutral) up to strongly purifying, and with **dN/dS emerging as
an output** rather than being imposed.

::: warning
Everything in this supplement lives in `zombi2.experimental`. It is shipped so you can use and iterate
on it, but it has **not** yet cleared ZOMBI2's core bar: the API may change, and the outputs are not yet
validated for publication. Install the optional dependencies with `pip install 'zombi2[selection]'`
(PyTorch, `fair-esm`, SciPy). This is a separate, self-contained document; it is **not** part of the
main manual.
:::

## Substitution = mutation times selection

The model rests on one decomposition. A substitution that reaches fixation is the product of two
things: a **mutation** happening, and that mutation being **fixed** rather than lost. ZOMBI2's ordinary
neutral model already supplies the first. The language model supplies the second.

Concretely, write the neutral (mutation) rate from one state $a$ to another $b$ as $\mu_{ab}$, and let
the language model assign each state a *preference* $\pi^{\text{pref}}$ — how protein-like that choice is
at that site. The selection-aware substitution rate is

$$Q_{ab} \;=\; \mu_{ab}\;\cdot\; h\!\left(F_b - F_a\right), \qquad F_x = \beta \,\ln \pi^{\text{pref}}_x ,$$

with the fixation factor $h(x) = x / (1 - e^{-x})$ (and $h(0)=1$). This is the classical
**Halpern–Bruno / Sella–Hirsh** mutation–selection form: $F$ is a scaled fitness read off the language
model, and $h$ is the probability that a mutation of fitness difference $F_b - F_a$ fixes in a
population. Its stationary distribution is the mutation process tilted by preference,
$\pi^{\text{target}} \propto \pi^{\text{mut}} \cdot (\pi^{\text{pref}})^{\beta}$, and — crucially — at
$\beta = 0$ the factor $h$ is identically $1$ and the kernel reduces **exactly** to the underlying
neutral model. Selection is a strict overlay: turn it off and you get plain ZOMBI2 back.

The one knob is `beta`. It is the population-scaled strength of selection (the $2N_e s$ axis): $0$ is
neutral drift, larger values pull ever harder toward the sequences the language model prefers. We treat
`beta` as the raw dial and *measure* its consequences (below), rather than calibrating it away.

::: note
The critic is pluggable. `ESM2Critic` is the first implementation, but anything that can turn a protein
into a per-site amino-acid preference satisfies the `Critic` interface — a smaller PLM, a hand-supplied
profile (`FixedProfileCritic`, useful for tests and for injecting a known constraint), or a future
model. Everything below is written against `Critic`, not against ESM2 specifically.
:::

## Two ways to read the critic: frozen and live

A language model scores a *whole* sequence, so *when* you ask it matters. There are two modes, and they
differ only in how often the critic is consulted as a family evolves down its gene tree.

![The two modes. **Frozen** (default): the critic is read once, on the root protein; each site's
preference is baked in and the sites then evolve independently — no epistasis, one language-model call
per gene, embarrassingly parallel. **Live**: the critic is re-read on the *current* sequence every
`refresh` substitutions per site along each lineage, so every site feels the others' current states
(within-gene epistasis) — many calls, and as `refresh` grows the live mode returns to frozen.](figures/selection_frozen_live.pdf){width=100%}

**Frozen** is the default and the workhorse. Reading the critic once, at the root, gives a per-site
amino-acid preference that is fixed for the whole family. Each site then evolves under its own
Halpern–Bruno process independently of the others — closed-form, cheap, and trivially parallel across
sites and families. Frozen has no epistasis by construction: a site's preferred residue does not depend
on what the neighbouring sites currently are.

**Live** captures within-gene epistasis at a cost. Every `refresh` substitutions per site along a
lineage, the critic is re-run on the sequence *as it currently stands* and the per-site preferences are
refreshed. A site that has drifted now re-scores its neighbours in their new context, so compensatory
and coupled changes can arise. This is a forward-simulation of the epistatic process, not an importance
sampler over a tilted distribution — a deliberate design choice, so the trajectory stays a faithful
sample of the model rather than collapsing onto the language model's single most "typical" protein.
Setting `refresh` to infinity recovers frozen exactly.

::: tip
Frozen recovers most of the realism gain for a tiny fraction of the cost, and it has no ensemble
degeneracy at high `beta`. Reach for live only when within-gene epistasis is the point of the
experiment.
:::

## Coding DNA, codons, and emergent dN/dS

Amino-acid-level selection is the right picture for a protein, but the observable that molecular
evolutionists care about — **dN/dS** — lives at the level of *codons*. The codon version of the model
evolves the coding DNA one codon at a time: **mutation acts on the nucleotides** (an ordinary
nucleotide model — JC69/K80/HKY85/GTR — supplies $\mu$), while **selection acts on the amino acid the
codon encodes** (the language model's preference for that residue).

That split gives the headline result for free. A **synonymous** mutation does not change the amino acid,
so its fitness difference is zero and its fixation factor is $h(0)=1$: synonymous sites evolve
neutrally, and $dS = 1$ by construction. A **non-synonymous** mutation changes the residue and is
scrutinised by the language model. The genome-wide ratio $\omega = dN/dS$ is therefore not a parameter
of the model — it is an **output**, an emergent consequence of how hard selection is pushing.

![Emergent dN/dS. **Left**: a codon substitution factorises into a nucleotide mutation times a
language-model fixation factor. A synonymous change (same amino acid) has $h=1$ and is neutral; a
non-synonymous change (new amino acid) is weighted by the critic's preference. **Right**: the model's
expected genome-wide $\omega$ as `beta` rises, computed from the codon model on a natural protein. It is
exactly $1$ at $\beta=0$ (neutral) and decreases monotonically toward $0$ under stronger purifying
selection — with realistic values ($\omega \approx 0.15$–$0.3$) reached around $\beta \approx 1$.](figures/selection_dnds.pdf){width=100%}

`CodonSelection.dnds(protein)` returns this expected $\omega$ analytically, so you can read the dN/dS a
given `beta` implies without simulating; `calibrate_beta(critic, protein, target_dnds)` inverts it, so
you can instead *ask* for a target $\omega$ and get the `beta` that delivers it.

::: note
There is a units subtlety worth stating. Both coding and non-coding evolution share one branch-length
scale — substitutions per **nucleotide** site — so the codon clock is normalised per nucleotide site
(three per codon). This makes a coding block at $\beta=0$ diverge at exactly the same neutral rate as an
intergenic block on the same tree; conserved sites then accrue proportionally fewer substitutions, which
is the visible face of $\omega < 1$.
:::

## Whole genomes: selection block by block

The flagship use is the nucleotide genome model: start from a **real annotated genome** at the root and
let it evolve, with selection acting on its genes. This raises a structural question — genes duplicate,
transfer, invert, and are lost, so "the tree of a gene" is not the species tree — which ZOMBI2 already
answers with its **block** decomposition.

![The block-based pipeline. **(1)** The root genome is partitioned by its GFF into coding genes and the
intergenic gaps between them. **(2)** The nucleotide model runs the full structural simulation and the
outcome is traced back into *blocks* — maximal never-cut intervals, each carrying its own gene tree.
Because a gene is never split by a breakpoint (ZOMBI2's "Design S"), a whole coding sequence is exactly
one block and evolves as one unit down one tree: a **gene block** evolves under language-model codon
selection, an **intergene block** drifts neutrally. **(3)** The evolved blocks are reassembled, in
genome order, into the DNA at every node — the root reproducing the input genome exactly.](figures/selection_genome_blocks.pdf){width=100%}

The important guarantees are that selection always runs **on a gene's own gene tree** (never the species
tree), so duplications, transfers and speciations within a gene family all inherit the ancestral
protein's constraint; and that mixing selected gene blocks with neutral intergene blocks reassembles
cleanly, because every block is evolved and stored in the same coordinate frame. Genes that cannot be
put under selection — a novel gene that originated mid-tree with no real coding sequence, a frame or
premature-stop problem, or a run with no root sequence — fall back to neutral evolution and are reported,
never silently mishandled.

## Using it

The library entry points live in `zombi2.experimental`:

```python
from zombi2.experimental import (
    ESM2Critic, PLMSelection,            # amino-acid selection over a gene tree
    CodonSelection, calibrate_beta,      # codon selection + emergent dN/dS
    read_cds_gff, simulate_nucleotide_selection,   # whole annotated genomes (block-based)
)

critic = ESM2Critic("esm2_t6_8M_UR50D")           # small by default; swap in a bigger model to taste
cds = read_cds_gff("genome.gff")                  # CDS with strand + reading frame
result, report = simulate_nucleotide_selection(
    species_tree, genome_str, cds, critic=critic, beta=1.0,
    inversion=0.01, duplication=0.01, loss=0.01, seed=1,
)
evolved_root = result.node_sequence(species_tree.root)   # == genome_str
```

The same pipeline is exposed on the command line as an experimental subcommand:

```
zombi2 experimental selection -t species_tree.nwk --gff genome.gff --genome-fasta genome.fna \
    --beta 1.0 --dup 0.01 --loss 0.01 --seed 1 -o out/
```

It writes the evolved genome at every node (`Genomes/<node>.fasta.gz`), the block architecture, the
extant per-gene alignments, and a selection report. Choose the critic model with `--esm-model` (the
small `esm2_t6_8M_UR50D` by default, or a large one such as `esm2_t33_650M_UR50D`), and set the
selection strength either directly with `--beta` or by target with `--target-dnds`, which calibrates a
single genome-wide `beta` from the root proteins.

::: warning
The critic is a neural network: the large ESM2 models want a GPU, and the feature targets small-to-medium
trees on a cluster rather than millions of tips. State this limitation plainly when you use it.
:::

## Status and limitations

- **Experimental.** In `zombi2.experimental`, behind `pip install 'zombi2[selection]'`; APIs may change
  and outputs are not yet validated.
- **Pseudogenization.** A gene that is pseudogenized on some lineage currently stays under selection on
  the whole gene-block tree; a per-lineage switch to neutral evolution after the loss-of-function edge is
  a planned refinement.
- **Codon epistasis.** The live (epistatic) mode is implemented at the amino-acid level; codon-level
  selection is frozen (per-site) for now.
- **Realism, measured.** With a real ESM2 critic the model reproduces the expected behaviour — $\omega=1$
  at $\beta=0$, decreasing to realistic purifying values as `beta` rises — and a Fréchet-style distance
  in the language model's own embedding space (`frechet_esm_distance`) quantifies how close simulated
  proteins land to real ones.
