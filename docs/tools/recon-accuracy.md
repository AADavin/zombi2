# Reconciliation accuracy

**`recon-accuracy`** is the mirror of [ALElite](reconciliation-likelihood.md): ALElite asks
*how probable* a reconciliation is; `recon-accuracy` asks *how close an inferred reconciliation
is to the truth*. Given the true reconciliation a simulation emits and an inferred reconciliation
for the **same gene tree**, it reports — node by node — how much of the event history was
recovered. It is the scoring half of a *simulate → infer → score* benchmark of reconciliation
methods (ALE, GeneRax, ecceTERA, …).

## What it computes

For every internal node of the shared extant gene tree, `recon-accuracy` compares the true and
inferred annotation and reports:

| Quantity | Measures |
| --- | --- |
| **Event accuracy** | fraction of nodes whose inferred event (S / D / T, plus G pseudogenization, C conversion) is correct |
| **Per-class precision / recall** | for each event class, the standard precision, recall, and F1 — the usual way to score D/T/L recovery against a simulated truth |
| **Mapping accuracy** | fraction of nodes mapped to the correct species branch (the LCA / MRCA mapping) |
| **Joint accuracy** | fraction of nodes with **both** the event *and* the mapping correct |
| **Transfer recovery** | for each true transfer: was it detected (called a transfer), and were its **donor** and **recipient** branches recovered? |

## Scope: fixed topology, node-by-node

Both reconciliations must annotate the **same** extant gene tree — identical tip labels,
identical branching — as when an inference method is run on the fixed simulated gene tree. The
tool aligns the two trees **structurally**, matching children by the set of tips beneath them, so
it is insensitive to Newick child order and correctly handles the unary pseudogenization nodes
ZOMBI2 keeps in the extant tree. A difference in topology or tip labels is an error, not a
partial score.

!!! note "Losses"
    A loss leaves no node in the extant gene tree, so losses fall outside node-by-node scoring.
    (A per-species-branch loss-count comparison would be a separate, topology-agnostic tool.)

Inputs are ZOMBI2 **annotated reconciled Newicks** — a `Reconciliation` (its `.extant` string is
used) or the string itself, with internal labels `"<species branch>|<EVENT>"`
(`"<donor>|T>recipient"` for transfers) and tips `"<species>|<gid>"`, exactly as written by
`zombi2 genomes --write reconciliations` / `reconcile()`. Converting a third-party method's
output into that format is glue the caller supplies; this tool **scores, it does not run
inference**.

## When to use

- **Benchmark a reconciliation method** against a ZOMBI2 truth: how well does it recover events,
  species mappings, and transfers under known conditions?
- **Measure the cost of a modelling choice** (e.g. undated vs dated inference) as event / mapping
  accuracy against the generating scenario.
- **Stress-test transfer inference** specifically — the donor/recipient recovery isolates the
  hardest part of DTL reconciliation.

## Function

`reconciliation_accuracy(truth, inferred) -> ReconAccuracy`

`truth` and `inferred` are each a `Reconciliation` or an annotated reconciled-Newick string.
The result is a namedtuple:

| Field | Meaning |
| --- | --- |
| `n_nodes` | internal gene-tree nodes compared |
| `event_accuracy`, `mapping_accuracy`, `joint_accuracy` | fractions in `[0, 1]` |
| `per_event` | `{event_char: EventPR(precision, recall, f1, tp, support_true, support_pred)}` |
| `transfer` | `TransferRecovery(n_true, detected, donor_correct, recipient_correct, both_correct)` (node counts) |

A gene tree with no internal nodes (a lone surviving copy) yields zero comparisons and all-zero
accuracies.

## Python usage

```python
from zombi2.tools import reconciliation_accuracy

# truth: a ZOMBI2 Reconciliation; inferred: an inference method's output for the SAME gene tree,
# converted to ZOMBI2's annotated reconciled-Newick format.
a = reconciliation_accuracy(truth_recon, inferred_recon)

print(a.event_accuracy, a.mapping_accuracy)
print(a.per_event["T"].recall)            # transfer detection recall
print(a.transfer)                         # TransferRecovery(n_true, detected, donor, recipient, both)
```

A reconciliation scored against itself is perfect by construction — a useful sanity check:

```python
a = reconciliation_accuracy(recon, recon)
assert a.event_accuracy == 1.0 and a.transfer.n_true == a.transfer.both_correct
```

## Command line

```bash
# 1. simulate a truth — --write reconciliations writes reconciled_extant.nwk (tips <species>|<gid>)
zombi2 genomes -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15 \
    --initial-families 50 --seed 7 --write reconciliations -o truth/

# 2. run your method on the same extant gene trees, then score it against the truth
zombi2 tools recon-accuracy -t truth/reconciled_extant.nwk -i inferred_recon.nwk
```

`-t` and `-i` are files of annotated reconciled Newicks, **one family per line**, paired line by
line (same count). `-o DIR` writes `DIR/reconciliation_accuracy.tsv`; otherwise the table prints
to stdout. Blank/`#`-comment lines are skipped.

## Output

One row per family, then a pooled (micro-averaged over all nodes) summary line:

```
family  n_nodes  event_acc  mapping_acc  joint_acc  transfers  transfers_recovered
1       9        0.888889   1.000000     0.888889   0          0
2       4        1.000000   1.000000     1.000000   0          0
...
# pooled over 22 family(ies), 212 node(s): event_acc=0.9623 mapping_acc=1.0000 joint_acc=0.9623 | transfers: 20/20 detected, 15/20 donor+recipient recovered
```

The pooled line is prefixed with `#` so TSV parsers skip it. `transfers` is the number of true
transfers in a family and `transfers_recovered` the number with donor **and** recipient correct.

## Validation

Known-answer cases on hand-built reconciliations plus round-trips on real simulations
(`tests/test_recon_accuracy.py`):

- **Perfect against itself.** Every extant family of a real simulation scores `event_accuracy =
  mapping_accuracy = joint_accuracy = 1.0` with all transfers fully recovered —
  `::test_real_reconciliations_perfect_against_themselves`,
  `::test_identical_reconciliation_is_perfect`.
- **Order-insensitive.** Swapping Newick child order does not change the score —
  `::test_child_order_does_not_matter`.
- **One wrong event.** Mislabelling one speciation as a duplication drops event accuracy to `2/3`
  and shows exactly in the per-class precision/recall — `::test_one_wrong_event_type`.
- **Mapping vs event.** A right event on the wrong species branch lowers `mapping_accuracy` but
  not `event_accuracy` — `::test_wrong_species_mapping_only`.
- **Transfer recovery.** A missed transfer (`detected = 0`), a wrong recipient (`(1,1,1,0,0)`),
  and a wrong donor (`(1,1,0,1,0)`) are each accounted correctly —
  `::test_transfer_missed`, `::test_transfer_wrong_recipient`, `::test_transfer_wrong_donor`.
- **Contract.** A topology or tip-label mismatch raises a clear error —
  `::test_topology_mismatch_raises`, `::test_tip_label_mismatch_raises`.

## References

- Szöllősi, G. J., Rosikiewicz, W., Boussau, B., Tannier, E. & Daubin, V. (2013). Efficient
  exploration of the space of reconciled gene trees. *Systematic Biology* 62(6): 901–912.
- Morel, B., Kozlov, A. M., Stamatakis, A. & Szöllősi, G. J. (2020). GeneRax: a tool for
  species-tree-aware maximum likelihood-based gene family tree inference under gene duplication,
  transfer, and loss. *Molecular Biology and Evolution* 37(9): 2763–2774.
- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a
  phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages.
  *Bioinformatics* 36(4): 1286–1288.
