# Gene-family coupling

By default every gene family in ZOMBI2 evolves independently, so the phylogenetic profile
correlates families *only* through the shared species tree. Real genomes correlate through
**function**: families in the same pathway or complex tend to be present or absent together
[@pellegrini1999profiles]. ZOMBI2 can inject that structure directly — a prescribed pairwise
coupling $J$ (and fields $h$) that makes families gain and lose non-independently — so the simulated
profiles carry a *known ground-truth* coupling. It is the forward, generative counterpart of the
inverse-Potts and direct-coupling-analysis methods that infer functional linkage from profiles
[@croce2019multiscale; @fukunaga2022ipm], and the benchmark those methods can be tested against.

## The model

Presence or absence of a fixed panel of `N` families inside one genome is an Ising vector
$\sigma \in \{0,1\}^N$. Fields $h_i$ and couplings $J_{ij}$ (symmetric, zero diagonal) define the
**local field** family $i$ feels — its intrinsic bias plus a contribution from every *present*
partner:

$$f_i = h_i + \sum_j J_{ij}\,\sigma_j$$

where the sum runs over present partners.

![A genome is a present/absent vector $\sigma \in \{0,1\}^N$ over a fixed panel of gene families; colour marks which pathway module each present family belongs to.](figures/potts_genome.pdf)

**Coupling enters through loss.** A present family is lost at a rate that is modulated by its local
field:

```
loss_i = base_loss * exp(-beta * f_i)
```

![Coupling enters through the loss rate: a present family is lost at rate $\mathrm{base\_loss}\cdot e^{-\beta f}$, so a high local field (partners present) protects it while a low field lets it go quickly.](figures/potts_lossrate.pdf)

So a present partner with $J_{ij} > 0$ raises $f_i$ and **lowers** family $i$'s loss: the two
families protect each other and co-occur. A partner with $J_{ij} < 0$ **raises** the loss: the
families purge each other and avoid one another. $J_{ij} = 0$ is independence. The field $h_i$ is the
solo retention bias — a large positive $h$ is a near-universal "hub" gene — and $\beta$ is a global
coupling strength.

**Gain is horizontal transfer.** A lost family returns only via the stock, **field-blind** transfer
event: a donor that still carries the family passes a copy to a recipient, with no regard for the
recipient's local field. The coupled loss then **selectively retains** the arrival — kept where its
partners are present and $f_i$ is high, quickly purged where they are absent. That *differential
retention of horizontally acquired genes* is the mechanism that writes $J$ into the profiles.

Two deliberate modelling choices underlie this. Gain is horizontal transfer rather than an explicit
detailed-balance rate, which is mechanistically honest to how genes actually re-enter genomes; and
the state is presence/absence, so copy number is ignored beyond being greater than zero — transfers
default to *replacement* so that re-acquisition does not stack copies.

## Building a coupling

A `CouplingSpec` holds the panel size, the couplings $J$, the fields $h$, and the base rates. Build
it from pathway blocks, from a dense matrix, or from a sparse edge list:

```python
import zombi2 as z

# (a) pathway blocks: families 0-2 co-occur, 3-5 co-occur,
#     the two blocks mutually exclusive
spec = z.pathway_blocks([3, 3], within=3.0, between=-1.0,
                        h=2.0, base_loss=1.0, transfer=0.2, beta=1.0)

# (b) a dense N x N coupling matrix (diagonal ignored)
J = [[0, 3, 0], [3, 0, 0], [0, 0, 0]]
spec = z.CouplingSpec.from_dense(J, h=2.0, base_loss=1.0, transfer=0.2)

# (c) a sparse edge list {(i, j): J_ij} (symmetrised)
spec = z.CouplingSpec.from_edges(6, {(0, 1): 3.0, (0, 2): 3.0,
                                     (0, 3): -2.0}, h=2.0)
```

![The coupling graph: families are nodes and couplings are edges. Positive couplings (solid) tie families into co-occurring modules; a negative coupling (dashed) makes two modules mutually exclusive.](figures/potts_coupling.pdf)

The parameters have direct interpretations:

- **`within`** and **`between`** (in `pathway_blocks`) set the coupling *inside* a block — positive
  values make pathway members co-occur — and *across* blocks — negative values make blocks mutually
  exclusive "rival" pathways, while `0` leaves the blocks independent.
- **`h`** is a scalar applied to every family, or a length-`N` vector of per-family biases.
- **`base_loss`** is the loss at $f_i = 0$, **`transfer`** is the per-copy horizontal gain rate, and
  **`beta`** is the global coupling strength.
- **`origination`** is an optional background rate of brand-new *uncoupled* families; leaving it at
  `0` keeps the panel closed.

Panel families are named `F0` through `F{N-1}`, and `spec.dense_J()` materialises the coupling matrix
for inspection.

## Running a coupled simulation

Simulate a species tree, then evolve the coupled panel over it:

```python
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=100, age=6.0, seed=1)
res = z.simulate_coupled(tree, spec, seed=1)

res.profiles          # ProfileMatrix: N panel rows x extant species
res.profiles.presence()
```

`simulate_coupled` seeds the root with the whole panel present — override this with
`initial_presence`, a length-`N` 0/1 mask — and returns a `CoupledResult` exposing `.profiles`,
`.leaf_genomes`, `.event_log`, and `.spec`. Every panel row is kept in the profile matrix, including
families absent from all extant species. Transfers default to `TransferModel(replacement=1.0)`; pass
`transfers=` to customise them, for example with distance-weighted recipients.

![Blind gain, selective retention: horizontal transfer re-acquires a family regardless of context, and the coupled loss then keeps it where its partners are present but purges it where they are absent — the differential retention that writes the coupling $J$ into the profiles.](figures/potts_retention.pdf)

Because `PottsRates` is a custom rate model, a coupled simulation runs on the pure-Python engine: the
coupling breaks the per-family independence that the Rust fast path assumes. The cost is
$O(N + \mathrm{nnz}(J))$ per event, which is comfortable at benchmark scale.

## Two caveats worth knowing

::: note
**Recovered couplings track the injected $J$ in sign and rank, not as a clean multiple.** Because the
gain channel is field-blind, detailed balance does not hold and the process has no exact Boltzmann
stationary distribution. Couplings recovered from the profiles follow the injected $J$
*monotonically* — positive to co-occurrence, negative to avoidance — but not as an affine constant.
This is the price of keeping regain mechanistic, since a lost family returns only from a donor that
still carries it. An exact-Boltzmann mode, using an explicit Glauber gain, is a documented extension
rather than the default.
:::

::: warning
**Control for the phylogeny.** On an ordinary birth–death tree even *uncoupled* families co-occur,
because loss is clade-restricted and families share ancestry. The coupling model's own tests isolate
the injected $J$ on a **near-star** tree — all lineages split just below the root and then evolve
independently. This is exactly why real inference corrects for the tree [@fukunaga2022ipm], and any
downstream benchmark should respect the same caveat.
:::

## Coupling a trait instead of another family

The same retention mechanism drives trait-linked gene families (the coevolution chapter). There a
family's loss is modulated by a **trait** value, `loss_i = base_loss * exp(-effect * w_i * s)`,
rather than by its partners' presence — the coupling field with the trait standing in for the coupled
state. Use this chapter's model when you want family-to-family structure, and the trait-linked model
when you want gene content to track a phenotype.
