# DNA substitution models

> **Tutorial:** see the [Sequences](../guide/sequences.md) guide.

Once `genomes` has produced gene trees and `sequence` has rescaled them into **substitutions per
site** (a phylogram), ZOMBI2 evolves an actual DNA alignment down each tree under a **continuous-time
Markov substitution model** over the four bases `ACGT`. Every model here is time-reversible and
normalised to one expected substitution per site per unit branch length, so the branch lengths mean
exactly what they say. The transition matrix over a branch of length `t` is `P(t) = exp(Q·t)`,
computed numpy-only by a reversible eigendecomposition. The models below are the classic nested
family — JC → K80 → HKY85 → GTR, each relaxing one assumption of the last — and any of them can carry
**+Gamma** among-site rate variation.

| Model | Parameters | Reach for it when |
| --- | --- | --- |
| **JC69** (JC) | none — equal rates, equal base freqs | the simplest baseline, no free parameters |
| **K80** (K2P) | `kappa` (transition/transversion ratio) | transitions and transversions differ but bases are even |
| **HKY85** (HKY) | `kappa` + unequal base freqs | a ti/tv bias *and* skewed base composition |
| **GTR** | 6 exchangeabilities + base freqs | the most general reversible DNA model |
| **+Gamma** (GammaRates) | `shape` (α), overlaid on any of the above | rates vary from site to site |

## The models

### Jukes–Cantor (JC69)

`jc69()` — the one-parameter-free baseline: all six exchangeabilities equal, all four base
frequencies `0.25`. Every substitution is equally likely, so a branch of length `t` gives the closed
form `P_ii = 1/4 + 3/4·e^{-4t/3}`. Reach for it as the null model or a sanity check.

### Kimura 2-parameter (K80)

`k80(kappa=2.0)` — separates **transitions** (A↔G, C↔T) from **transversions** by a single ratio
`kappa`, keeping base frequencies equal at `0.25`. `kappa 1` reduces to JC69. Use it when the ti/tv
bias matters but base composition is even.

### HKY85 (HKY)

`hky85(kappa=2.0, freqs=(0.25,0.25,0.25,0.25))` — K80's transition bias `kappa` plus **unequal
equilibrium base frequencies** `(A,C,G,T)`. Use it when the sequences have both a ti/tv bias and a
skewed base composition; the long-branch stationary distribution recovers `freqs`.

### GTR

`gtr(rates=(1,1,1,1,1,1), freqs=(0.25,0.25,0.25,0.25))` — the **general time-reversible** model: six
free exchangeabilities `[AC,AG,AT,CG,CT,GT]` and arbitrary base frequencies. It is the most general
reversible DNA model and the superclass of all the others (`rates` all `1`, `freqs` all `0.25` gives
JC69). Reach for it when you want no built-in symmetry assumptions.

### +Gamma among-site rate variation (GammaRates)

`GammaRates(shape, k=4)` — not a substitution model on its own but an **overlay** on any of the four
above. Sites are binned into `k` equal-probability discrete-Gamma categories (Yang 1994) with mean
rate 1; a small `shape` (α) makes rates highly heterogeneous across sites, a large `shape` makes them
nearly uniform. Use it whenever real substitution rates vary from one site to the next.

## Command line

Sequence simulation is its own step: run [`zombi2 sequence`](../cli.md#sequence) on a `genomes` run
that was written with `trace` in `--write`. `--subst-model` picks the model (`jc69`/`k80`/`hky85`/`gtr`
for DNA; the name auto-detects DNA vs protein), and the model-specific knobs are `--kappa`,
`--base-freqs`, `--gtr-rates`, and `--gamma-shape`. Omit `--subst-model` to only rescale the trees
without simulating sequences.

```bash
# a genomes run written with the event trace (prerequisite)
zombi2 species --birth 1 --death 0.3 --tips 8 --age 3 --seed 1 -o run/
zombi2 genomes -t run/species_tree.nwk --dup 0.3 --trans 0.1 --loss 0.3 \
    --orig 0.5 --write trace --seed 1 -o run/

# JC69 — no free parameters
zombi2 sequence --genomes run/ --subst-model jc69 --branch-speed 0.4 --seed 7 -o jc/

# K80 with a transition/transversion ratio
zombi2 sequence --genomes run/ --subst-model k80 --kappa 4 --seed 7 -o k80/

# HKY85: ti/tv bias plus unequal base frequencies (A C G T)
zombi2 sequence --genomes run/ --subst-model hky85 --kappa 4 \
    --base-freqs 0.4 0.1 0.1 0.4 --seed 7 -o hky/

# GTR: 6 exchangeabilities [AC AG AT CG CT GT] + base freqs, with +Gamma across sites
zombi2 sequence --genomes run/ --subst-model gtr \
    --gtr-rates 1 2.5 1 1 2.5 1 --base-freqs 0.3 0.2 0.2 0.3 \
    --gamma-shape 0.5 --branch-speed 0.4 --seed 7 -o gtr/
```

`--seq-length N` sets the alignment length (default 300); `--root-fasta` seeds each family's root
from a FASTA instead of a random draw. The `--branch-speed`/`--family-speed`/`--clock` knobs govern
the relaxed clock that turns the time tree into a phylogram — see the
[substitution branch-lengths](../cli.md#sequence) section.

## Python

The models live in `zombi2.sequences` as factory functions that build a `SubstitutionModel`, plus
`GammaRates` and the `evolve_on_tree` simulator (the same objects also re-export at the top level, so
`zombi2.hky85` works too):

```python
import numpy as np
from zombi2.sequences import jc69, k80, hky85, gtr, GammaRates, evolve_on_tree

# a minimal tree node has a .gid and .children (real trees come from a genomes run)
class Node:
    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)

a, b = Node("a"), Node("b")
root = Node("r", [a, b])
subst = {root: 0.0, a: 0.2, b: 0.2}          # branch lengths in substitutions/site

# HKY85 (kappa + unequal freqs) with +Gamma across-site rate variation
model = hky85(kappa=4.0, freqs=(0.4, 0.1, 0.1, 0.4))
seqs = evolve_on_tree(root, subst, model, np.random.default_rng(0),
                      length=1000, gamma=GammaRates(shape=0.5))
# seqs maps each node's gid -> its DNA string (internal nodes and tips)

# GTR with six explicit exchangeabilities [AC AG AT CG CT GT]
g = gtr(rates=(1, 2.5, 1, 1, 2.5, 1), freqs=(0.3, 0.2, 0.2, 0.3))
```

In practice you drive this through the `zombi2 sequence` command (or `SequenceEvolution`), which
supplies the rescaled gene trees; `evolve_on_tree` is the low-level engine underneath.

## Output

`zombi2 sequence --subst-model MODEL` writes, under `-o`:
`gene_trees/` with the substitution-unit phylograms (`<fam>_complete_subst.nwk` and
`<fam>_extant_subst.nwk` per family), `alignments/<fam>.fasta` with the simulated per-family DNA
alignment (one record per surviving gene copy), `branch_rates.tsv` and `gene_family_speeds.tsv`
recording the clock rate on each branch and the per-family speed multiplier, and `sequence.log` the
run manifest (every parameter, including the substitution model and its `--kappa`/`--base-freqs`/
`--gtr-rates`/`--gamma-shape` settings).

## Validation

- **JC69.** The numpy-only `exp(Qt)` matches the JC69 closed form `P_ii = 1/4 + 3/4·e^{-4t/3}` at
  several branch lengths (`test_sequence_sim.py::test_p_matrix_matches_jc_closed_form`), and a pair of
  tips `0.2` from the root recovers the true `0.4` JC-corrected distance to within `0.02`
  (`test_sequence_sim.py::test_jc_distance_recovered`).
- **K80.** The observed transition/transversion structure across an evolved branch matches the
  Kimura-1980 closed form implied by `kappa` — the per-site transition and transversion difference
  fractions, and their ratio, agree with the closed-form `p_ti`/`p_tv` derived from `kappa` (an
  oracle, not just matrix invariants)
  (`test_sequence_sim.py::test_k80_transition_transversion_matches_kappa`).
- **HKY85.** A single long-branch star run recovers the specified unequal base frequencies
  `(0.4,0.1,0.1,0.4)` to within `0.01`
  (`test_sequence_sim.py::test_stationary_frequencies_recovered`).
- **GTR.** A long GTR run, started from a uniform (25%-each) root that is *not* the target, recovers
  the asymmetric stationary base frequencies `pi=(0.1,0.2,0.3,0.4)` to within `0.01`
  (`test_sequence_sim.py::test_gtr_stationary_frequencies_recovered`).
- **+Gamma.** Overlaying the Gamma over-disperses the per-site divergence versus no-Gamma: the mean
  and variance of the per-site count of differing tips match the closed-form Binomial-vs-Gamma-mixture
  moments (law of total variance over the category mixture), and the +Gamma variance is strictly
  larger (`test_sequence_sim.py::test_gamma_overdisperses_per_site_divergence`). The discrete-Gamma
  category rates themselves average to exactly 1 with increasing, distinct categories
  (`test_sequence_sim.py::test_gamma_rates_numpy_only_mean_one`).

## References

- Jukes, T. H. & Cantor, C. R. (1969). Evolution of protein molecules. In *Mammalian Protein
  Metabolism*, 21–132.
- Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through
  comparative studies of nucleotide sequences. *Journal of Molecular Evolution* 16: 111–120.
- Hasegawa, M., Kishino, H. & Yano, T. (1985). Dating of the human–ape splitting by a molecular
  clock of mitochondrial DNA. *Journal of Molecular Evolution* 22: 160–174.
- Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences.
  *Lectures on Mathematics in the Life Sciences* 17: 57–86.
- Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates
  over sites: approximate methods. *Journal of Molecular Evolution* 39: 306–314.
