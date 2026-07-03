# The gene-family coupling model (implementation note)

*Companion to `docs/non_independence.tex`, which derives the theory. This note records
**what was built** in `zombi2/coupling.py`, the design decisions behind it, and how to use
and validate it.*

## What it is

By default ZOMBI2 evolves each gene family independently, so the phylogenetic profile
(`ProfileMatrix`) correlates families **only through the shared species tree**. This module
makes families evolve **non-independently**: a prescribed pairwise coupling `J` (and fields
`h`) biases which families are present or absent *together*, so the simulated profiles carry
a **known ground-truth** coupling — the forward, generative counterpart of the inverse-Potts
/ DCA methods that infer functional linkage from profiles (Pellegrini 1999; Croce 2019;
Fukunaga & Iwasaki 2022).

## The model

Presence/absence of a fixed panel of `N` families inside one genome is an Ising vector
`σ ∈ {0,1}ᴺ`. Fields `hᵢ` and couplings `Jᵢⱼ` (symmetric, zero diagonal) define a **local
field** seen by family `i`:

$$f_i(\sigma) = h_i + \sum_{j} J_{ij}\,\sigma_j \qquad(\text{partners only}).$$

**Coupling enters through loss.** A *present* family is lost at rate

$$\text{loss}_i = \text{base\_loss}\cdot \exp(-\beta\, f_i).$$

So a present partner with `Jᵢⱼ > 0` raises `fᵢ` and **lowers** `i`'s loss (they protect each
other → co-occurrence); `Jᵢⱼ < 0` **raises** loss (they purge each other → avoidance);
`Jᵢⱼ = 0` is independence. `hᵢ` is the family's solo retention bias (a large positive `h` is
a near-universal "hub" gene); `β` is a global coupling strength.

**Gain is horizontal transfer.** A lost family returns only via the stock, **field-blind**
`TRANSFER` event — a donor that still has the family transfers a copy into a recipient. The
coupled loss then **selectively retains** it: kept where its partners are present (high
`fᵢ`), quickly purged where they are absent. That *differential retention of horizontally
acquired genes* is the mechanism that writes `J` into the profiles.

### Two design decisions (confirmed with Adrián)

1. **Regain = HGT** (not an explicit Glauber gain rate). A lost family reappears only from a
   donor that still carries it — mechanistically honest to how genes actually re-enter
   genomes. The design note's alternative — an explicit detailed-balance gain
   `gainᵢ = ρ·exp(+fᵢ/2)` — would make the profiles *exact* Boltzmann samples of `P(σ) ∝
   exp(Σ hᵢσᵢ + Σ Jᵢⱼσᵢσⱼ)`. We deliberately did **not** take that route.
2. **State = presence/absence** (Ising). The field reads only whether partners are present,
   matching the binary profiles the inference methods consume. Copy number is ignored beyond
   `>0` (transfers default to *replacement* so re-acquisition doesn't stack copies).

### Faithfulness caveat

Because the gain channel is field-blind, **detailed balance does not hold** and the process
has no exact Boltzmann stationary distribution. Consequently, couplings recovered from the
profiles track the injected `J` in **sign and rank**, *not* as a clean constant multiple —
this is the "coupling limited by donor availability" trade-off noted in
`non_independence.tex`. It is the right choice for a benchmark that mirrors real regain, but
users comparing recovered `Ĵ` to injected `J` should expect a monotone, not affine,
relationship.

## How it fits the architecture

`PottsRates` is an ordinary `RateModel`. The simulator already calls
`event_weights(genome, branch, time)` with the **whole genome**, so `PottsRates` reads the
presence vector, computes each present family's field, and returns coupled loss weights + a
transfer channel — with **no change to the simulator, sampler, genome, event types, or
output**. A custom rate model is automatically ineligible for the Rust fast path, so a
coupled run uses the pure-Python engine (the coupling breaks per-family independence, so the
fast path could not apply regardless). Cost is `O(N + nnz(J))` per event.

## Usage

```python
import numpy as np
from zombi2 import simulate_species_tree, BirthDeath
from zombi2.coupling import CouplingSpec, pathway_blocks, simulate_coupled

tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=100, age=6.0, seed=1)

# (a) pathway-block J: families 0–2 co-occur, 3–5 co-occur, blocks mutually exclusive
spec = pathway_blocks([3, 3], within=3.0, between=-1.0,
                      h=2.0, base_loss=1.0, transfer=0.2, beta=1.0)

# (b) or a dense matrix / sparse edges directly:
#   CouplingSpec.from_dense(J, h=..., base_loss=..., transfer=...)
#   CouplingSpec.from_edges(N, {(0, 1): 3.0, (0, 2): -2.0}, h=...)

res = simulate_coupled(tree, spec, seed=1)
profiles = res.profiles          # ProfileMatrix, N panel rows × species (all rows kept)
```

`simulate_coupled` seeds the root with the whole panel present (override with
`initial_presence`, a length-`N` 0/1 mask) and returns a `CoupledResult`
(`.profiles`, `.leaf_genomes`, `.event_log`, `.spec`). Transfers default to
`TransferModel(replacement=1.0)`; pass `transfers=` to customise (e.g. `distance_decay`).

### Public API (to be documented on the mkdocs site)

`zombi2.coupling`: `CouplingSpec` (`.from_dense`, `.from_edges`, `.dense_J`),
`pathway_blocks`, `PottsRates`, `simulate_coupled`, `CoupledResult`. Not yet re-exported
from the top-level `zombi2` namespace — import from `zombi2.coupling`.

## Validation (`tests/test_coupling.py`)

- **Formula** — the emitted loss equals `base_loss·exp(-β·fᵢ)` exactly; a positive partner
  lowers a family's loss, a negative partner raises it; non-panel families stay uncoupled.
- **Ground-truth recovery** — inject `J` and confirm the profiles recover it: `+J` →
  co-occurrence, `-J` → avoidance, `J = 0` → no structure, uncoupled families uncorrelated.
- **Phylogeny control** — recovery is checked on a **near-star** tree (all lineages split
  just below the root, then evolve independently), which isolates the injected coupling from
  the shared-ancestry confounding that inflates co-occurrence on a normal birth–death tree.
  On an ordinary tree even uncoupled families co-occur (clade-restricted loss), which is
  exactly why real inference (Fukunaga & Iwasaki 2022) corrects for the phylogeny — a caveat
  any downstream benchmark should respect.

## Extensions (not built)

- **Exact-Boltzmann mode** — add an explicit Glauber gain channel (via the origination seam
  with a family chosen in `target_params`) for detailed-balance profiles when an exact target
  is wanted instead of mechanistic regain.
- **Copy-number (Potts, q>2)** — let copy number feed the field and modulate duplication.
- **Branch-varying `J`** — regime shifts across the tree (compose with a `BranchRates`-style
  wrapper); constant `J` matches the inference assumptions and is the current default.
