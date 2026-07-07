# Protein substitution models

> **Tutorial:** see the [Sequences](../guide/sequences.md) guide.

Once a `genomes` run's gene trees have been rescaled from time into substitutions per site by the
relaxed molecular clock, ZOMBI2 can evolve an **amino-acid alignment** along each tree. A protein substitution model is a 20-state **time-reversible** rate matrix `Q_ij = S_ij · π_j`,
built from a symmetric exchangeability matrix `S` and the stationary amino-acid frequencies `π`, and
normalised to one expected substitution per site. All the empirical models below are transcribed
byte-for-byte from the reference PAML data files (Ziheng Yang, `abacus-gene/paml`); ZOMBI2 differs
from a nucleotide model only in the alphabet (`K = 20`) — the same `evolve_on_tree` engine runs both.

| Model | Exchangeabilities `S` / frequencies `π` | Reach for it when |
| --- | --- | --- |
| **Poisson** | all equal / uniform (1/20) | you want a neutral, parameter-free baseline (F81-for-proteins) |
| **LG** | Le & Gascuel 2008 empirical | the modern default for most protein data |
| **WAG** | Whelan & Goldman 2001 empirical | a widely used, well-established general matrix |
| **JTT** | Jones, Taylor & Thornton 1992 empirical | matching an older analysis or JTT-based pipeline |
| **Dayhoff** | Dayhoff, Schwartz & Orcutt 1978 empirical | reproducing classic Dayhoff-PAM results |

## The models

### Poisson

Equal off-diagonal exchangeabilities and uniform frequencies (`π_i = 1/20`) — the protein analogue of
Felsenstein-81. No parameters and no amino-acid preference: every residue is equally likely and every
substitution equally rated. Use it as a neutral baseline or a sanity check, not as a realistic model
of protein evolution.

### LG

The Le & Gascuel (2008) exchangeabilities and frequencies, estimated across a large curated alignment
database with across-site rate variation accounted for during estimation. The current default choice
for most empirical protein work; reach for it unless a specific pipeline dictates otherwise.

### WAG

The Whelan & Goldman (2001) matrix, estimated by maximum likelihood from a broad set of globular
protein families. A long-established, general-purpose model; a reasonable choice when comparability
with the extensive WAG literature matters.

### JTT

The Jones, Taylor & Thornton (1992) matrix (PAML's `jones.dat`), derived from a large survey of
protein sequences. Use it to match an older analysis or a JTT-based reconstruction pipeline.

### Dayhoff

The original Dayhoff, Schwartz & Orcutt (1978) exchangeabilities (PAML's `dayhoff.dat`), the ancestor
of the PAM family. Mostly of historical and reproduction value; reach for it to recreate classic
Dayhoff-PAM results.

## Command line

Protein evolution runs through the `sequence` command: pass a prior `genomes` output directory and
select a protein model with `--subst-model` (DNA vs protein is auto-detected from the name). The
`genomes` run must have been done with `trace` in `--write`. Alignment length is set with
`--seq-length` (amino acids; default 300), and `--gamma-shape` adds discrete-Gamma across-site rate
heterogeneity. The empirical models take no further parameters — `S` and `π` are fixed. The lineage
clock knobs (e.g. `--branch-speed`, `--family-speed`) rescale the trees exactly as for DNA.

```bash
# 1) a genomes run recorded with the event trace
zombi2 species --birth 1 --death 0.3 --tips 8 --age 3 --seed 1 -o run/
zombi2 genomes -t run/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace --seed 1 -o run/

# 2) evolve LG protein alignments along the rescaled gene trees
zombi2 sequence --genomes run/ --subst-model lg --seq-length 200 \
    --branch-speed 0.4 --seed 7 -o seqs/

# WAG with across-site rate heterogeneity
zombi2 sequence --genomes run/ --subst-model wag --gamma-shape 0.5 --seed 7 -o seqs/
```

## Python

The model builders and the sequence engine live in `zombi2.sequences` (each name is re-exported at
the top level too, so `zombi2.lg` also works):

```python
import numpy as np
from zombi2.sequences import lg, poisson, make_model, evolve_on_tree, AMINO_ACIDS, PROTEIN_MODELS

model = lg()                       # or wag(), jtt(), dayhoff(), poisson()
model = make_model("wag")          # ...or build one by name
print(PROTEIN_MODELS)              # ('poisson', 'lg', 'wag', 'jtt', 'dayhoff')

# evolve amino-acid sequences over a gene tree (nodes carry .gid and .children;
# `subst` maps each node to the substitution length of the branch ending at it)
class Node:
    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)

tips = [Node(f"t{i}") for i in range(4)]
root = Node("r", tips)
seqs = evolve_on_tree(root, {t: 0.3 for t in tips}, model,
                      np.random.default_rng(0), length=60)
assert set("".join(seqs[t.gid] for t in tips)) <= set(AMINO_ACIDS)   # genuine 20-AA protein
```

## Output

`sequence` writes the rescaled trees to `gene_trees/` (`<family>_complete_subst.nwk` and
`<family>_extant_subst.nwk`, branch lengths in substitutions per site) and, with `--subst-model`, one
protein alignment per family to `alignments/<family>.fasta` (FASTA headers are the standard
[node names](../contributing/conventions.md#naming); each record is a 20-letter amino-acid sequence).
`branch_rates.tsv` and `gene_family_speeds.tsv` record the per-lineage clock rate and per-family speed
multiplier, and `sequence.log` is the run manifest.

## Validation

- **Empirical models (LG / WAG / JTT / Dayhoff).** The stored stationary frequencies match the
  published `π` vectors to 1e-4, catching any transcription error a mere reversibility check would
  miss (`test_sequence_sim.py::test_empirical_aa_frequencies_match_published`).
- **Poisson.** Uniform stationary frequencies and equal off-diagonal rates hold exactly by
  construction (`test_sequence_sim.py::test_poisson_is_exact`).
- **Stationarity (LG).** A star tree with long branches recovers the model's stationary amino-acid
  frequencies to within 0.02 (`test_sequence_sim.py::test_protein_stationary_recovered_on_long_branch`).

## References

- Le, S. Q. & Gascuel, O. (2008). An improved general amino acid replacement matrix.
  *Molecular Biology and Evolution* 25(7): 1307–1320.
- Whelan, S. & Goldman, N. (2001). A general empirical model of protein evolution derived from
  multiple protein families using a maximum-likelihood approach.
  *Molecular Biology and Evolution* 18(5): 691–699.
- Jones, D. T., Taylor, W. R. & Thornton, J. M. (1992). The rapid generation of mutation data
  matrices from protein sequences. *Computer Applications in the Biosciences* 8(3): 275–282.
- Dayhoff, M. O., Schwartz, R. M. & Orcutt, B. C. (1978). A model of evolutionary change in proteins.
  In *Atlas of Protein Sequence and Structure*, vol. 5, suppl. 3, pp. 345–352.
