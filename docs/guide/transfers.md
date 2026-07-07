# Transfers

Rates decide *how often* a transfer fires (the [rate model](gene-families.md)); a
**`TransferModel`** decides *what a transfer does*.

```python
from zombi2.genomes import simulate_genomes, TransferModel

genomes = simulate_genomes(
    tree, transfer=0.3, ...,
    transfers=TransferModel(
        replacement=0.2,      # additive vs replacement
        distance_decay=2.0,   # recipient choice by phylogenetic distance
        allow_self=False,     # self-transfer = duplication
    ),
)
```

The default `TransferModel()` is additive, uniform-recipient, no self-transfer.

## Additive vs replacement

- **Additive** (`replacement=0`): the recipient gains a copy (net +1).
- **Replacement** (`replacement=p`): with probability `p`, the transfer also removes one
  pre-existing copy of that family in the recipient — a net-zero swap. It is only possible
  when the recipient already has the family (otherwise the transfer is additive).

`replacement=1` makes every possible transfer a replacement. Replacement transfers appear
in the log as a transfer **plus** a compensating loss.

## Distance-dependent recipient choice

By default the recipient is drawn **uniformly** among lineages alive at the transfer time.
Set `distance_decay=λ` to favour phylogenetically close recipients: candidate `r` is
weighted by `exp(-λ · d)`, where `d = 2·(t − t_MRCA)` is the patristic distance between
donor and candidate at the transfer time `t`. Larger `λ` = more local transfers; distant
transfers are damped but never forbidden.

!!! note "Performance"
    Distance weighting computes, per transfer, the MRCA time of the donor with every
    co-existing lineage (`O(alive · depth)`). It is the one hot spot; for very large trees
    it can be swapped for a sparse-table LCA later.

## Self-transfer

With `allow_self=True` the donor lineage is an eligible recipient. A self-transfer creates
a second copy in the same genome — mechanically a **duplication**. This lets you drop
explicit duplications and run a transfer/loss-only model:

```python
from zombi2.genomes import simulate_genomes, TransferModel

simulate_genomes(tree, transfer=1.0, duplication=0.0,
                 transfers=TransferModel(allow_self=True),
                 max_family_size=0.5, seed=1)   # cap it — self-transfers grow like D
```

!!! warning
    Self-transfers grow families exactly as duplications do, so pair them with a growth
    cap (see [Bounding growth](growth.md)).
