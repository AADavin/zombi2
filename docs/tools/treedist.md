# Tree distances

**`treedist`** measures how far apart two trees are — the standard currency for benchmarking a
phylogenetic method against a known ground truth. ZOMBI2 emits a *true* tree (a species tree, or
a reconciled gene tree); an inference method emits its *estimate*; `treedist` reports the
distance between them, one number with a right answer.

Every metric is exact and defined over the **shared leaf-label set** of the two trees. A
mismatch in the leaf sets is an error, not a silent partial score.

## What it computes

Four classical distances between a reference tree and one or more comparison trees:

| Metric | Measures | Notes |
| --- | --- | --- |
| **Robinson–Foulds** | symmetric difference of the trees' clade (rooted) or bipartition (unrooted) sets | Fast, exact, the most widely reported topology distance — but coarse: one misplaced leaf can saturate it. |
| **Branch score** | the L2 (Kuhner–Felsenstein) or L1 norm of per-clade branch-**length** differences | Uses ZOMBI2's time-derived branch lengths, so it reads a dated tree directly. Matched clades contribute `(ℓ₁ − ℓ₂)²`; unmatched clades their full length against zero. |
| **Quartet distance** | fraction of four-leaf subsets whose induced unrooted topology differs | Finer-grained than RF (degrades gracefully). Exact but `O(n⁴)`, so guarded by `max_leaves`. |
| **Matching distance** | minimum-cost matching of one tree's splits/clusters to the other's | Does **not** saturate: a slightly displaced clade costs a little, not a whole unit. Solved as the assignment problem via **optional** SciPy. |

### Rooted vs unrooted

ZOMBI2 trees are rooted (`node.time` runs forward from a root at 0). Robinson–Foulds and the
matching distance therefore default to the **rooted** form (comparing clades / clusters — the
honest distance when the root is known, as it is for a simulated tree); pass `rooted=False` for
the unrooted form (bipartitions / splits) when comparing against a method that does not infer a
root. The quartet distance is intrinsically unrooted (a quartet has no root).

### The matching distance and SciPy

The matching distance is the only metric with an optional dependency: it is a minimum-weight
bipartite matching (the assignment problem), solved with
[`scipy.optimize.linear_sum_assignment`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html).
SciPy is **not** a core ZOMBI2 dependency — install it with `pip install scipy` or the
`zombi2[dev]` extra. When SciPy is absent the metric is *skipped* (its fields become `None` /
blank), never fabricated; the other three metrics need only NumPy.

## When to use

- **Benchmark a tree-inference method** — simulate a species or gene tree, run the method, and
  score the estimate against the truth.
- **Quantify the effect of a setting** on reconstructed topology (e.g. how a rate or a model
  choice moves the estimate away from the generating tree).
- **Compare a batch of replicates** (bootstrap or repeated inferences) against one reference in a
  single call.

## Functions

All functions accept a `zombi2.tree.Tree` **or** a Newick string for each tree.

| Function | Returns |
| --- | --- |
| `robinson_foulds(t1, t2, *, rooted=True)` | `RFResult(rf, max_rf, normalized)` |
| `branch_score(t1, t2, *, order=2)` | `float` (L2 for `order=2`, L1 for `order=1`) |
| `quartet_distance(t1, t2, *, max_leaves=100)` | `QuartetResult(differing, total, normalized)` |
| `matching_distance(t1, t2, *, rooted=True, max_leaves=2000)` | `MatchingResult(distance, max_distance, normalized)` |
| `compare_trees(t1, t2, *, quartet=True, matching=True, max_leaves=100, matching_max_leaves=2000, branch_score_order=2)` | `TreeComparison` (all metrics; skipped ones are `None`) |

`normalized` divides each raw distance by its maximum (for RF, `|splits₁| + |splits₂|`; for
quartet, `C(n, 4)`; for matching, the cost of leaving every split unmatched), giving a value in
`[0, 1]`.

## Python usage

```python
from zombi2.tree import read_newick
from zombi2.tools import compare_trees, robinson_foulds

truth = read_newick(open("true_tree.nwk").read())
est   = read_newick(open("inferred_tree.nwk").read())

# every metric in one call
c = compare_trees(truth, est)
print(c.rf, c.rf_normalized, c.branch_score, c.quartet, c.matching)

# or an individual metric (strings work too)
r = robinson_foulds("((A,B),(C,D));", "((A,C),(B,D));")
print(r.rf, r.normalized)          # 4 1.0  (rooted clusters, fully disjoint)
```

## Command line

```bash
# distances between a true species tree and an inferred one
zombi2 tools treedist -r true_tree.nwk -e inferred_tree.nwk

# score many replicate/bootstrap trees against one reference, saved to out/
zombi2 tools treedist -r true_tree.nwk -e replicates.nwk -o out/
```

`-r` is a Newick file with exactly one **reference** tree; `-e` is a Newick file of one or more
**comparison** trees (one per line), each compared to the reference and each sharing its leaf
labels. Options:

| Flag | Meaning |
| --- | --- |
| `-o DIR` | write `DIR/Tree_distances.tsv` (default: print the table to stdout) |
| `--no-quartet` | skip the quartet distance (it is `O(n⁴)`) |
| `--max-leaves N` | quartet guard: skip it above `N` leaves (default 100); raise to force it |
| `--branch-order {1,2}` | branch-score norm: `2` = L2 / Kuhner–Felsenstein (default), `1` = L1 |

## Output

One row per comparison tree, with columns:

```
tree  n_leaves  rf  rf_norm  rf_unrooted  branch_score  quartet  quartet_norm  matching  matching_norm
```

A blank cell means the metric was skipped: `quartet` when the tree exceeds `--max-leaves`,
`matching` when it exceeds the matching guard or SciPy is not installed. Example (a `((A,C),(B,D))`
estimate against an `((A,B),(C,D))` reference):

```
tree  n_leaves  rf  rf_norm   rf_unrooted  branch_score  quartet  quartet_norm  matching  matching_norm
1     4         4   1.000000  2            0.000000      1        1.000000      4         0.500000
```

## Validation

Each metric is checked against **hand-computed** known answers (`tests/test_treedist.py`), not
merely a reference binary:

- **Incompatible resolutions.** `((A,B),(C,D))` vs `((A,C),(B,D))`: rooted RF `= 4`, unrooted RF
  `= 2`, quartet `= 1/1`, matching `= 4` (normalized `0.5`) —
  `::test_rf_hand_computed_incompatible_quartet`, `::test_quartet_single_incompatible_quartet`,
  `::test_matching_hand_computed_incompatible_quartet`.
- **Root sensitivity.** Two trees with the same unrooted topology but different roots have
  unrooted RF `= 0` yet rooted RF `> 0` — `::test_rooted_differs_but_unrooted_agrees_under_reroot`.
- **Quartet topologies.** A caterpillar `(((A,B),C),D)` induces the split `AB|CD` (agreeing with
  `((A,B),(C,D))`), and an unresolved star differs from any resolved tree —
  `::test_quartet_caterpillar_topology`, `::test_quartet_unresolved_star_differs_from_resolved`.
- **Branch score.** A tree with two branch lengths changed gives L1 `= 2` and L2 `= √2` —
  `::test_branch_score_hand_computed`.
- **Matching padding.** A resolved tree vs a star leaves every cluster unmatched, so the matching
  distance equals its maximum (normalized `1.0`) — `::test_matching_padding_against_star`.

Identity (`d(t, t) = 0` for every metric) and symmetry are checked on both hand-built and
simulated trees; mismatched leaf sets and over-`max_leaves` inputs raise clear errors
(`::test_identical_trees_are_zero`, `::test_mismatched_leaf_sets_raise`,
`::test_quartet_guard_on_large_tree`, `::test_matching_guard_on_large_tree`).

## References

- Robinson, D. F. & Foulds, L. R. (1981). Comparison of phylogenetic trees. *Mathematical
  Biosciences* 53(1–2): 131–147. (Robinson–Foulds distance.)
- Kuhner, M. K. & Felsenstein, J. (1994). A simulation comparison of phylogeny algorithms under
  equal and unequal evolutionary rates. *Molecular Biology and Evolution* 11(3): 459–468. (The
  branch-score distance.)
- Estabrook, G. F., McMorris, F. R. & Meacham, C. A. (1985). Comparison of undirected
  phylogenetic trees based on subtrees of four evolutionary units. *Systematic Zoology* 34(2):
  193–200. (Quartet distance.)
- Bogdanowicz, D. & Giaro, K. (2012). Matching split distance for unrooted binary phylogenetic
  trees. *IEEE/ACM Transactions on Computational Biology and Bioinformatics* 9(1): 150–160.
- Lin, Y., Rajan, V. & Moret, B. M. E. (2012). A metric for phylogenetic trees based on matching.
  *IEEE/ACM Transactions on Computational Biology and Bioinformatics* 9(4): 1014–1022. (The
  matching cluster distance for rooted trees.)
