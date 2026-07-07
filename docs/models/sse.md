# State-dependent diversification (SSE)

State-dependent speciation and extinction (SSE) models let a **trait drive the shape of the tree**:
a lineage's character state sets its speciation and extinction rates, so the tree and the trait
must be grown **together** (Maddison, Midford & Otto 2007). This is the `traits:species` edge of the
[coevolution model](../coevolution_models.md) — an arrow *into* species, so the run is
**forward-only** and **produces** the tree rather than taking one. The models below differ only in
the kind of trait doing the driving: binary, k-state, hidden, continuous, or a trait that also jumps
*at* each speciation.

| Model | Trait | Reach for it when |
| --- | --- | --- |
| **BiSSE** | binary (0/1), each state its own λ, μ | one two-state character shapes a radiation |
| **MuSSE** | k discrete states | a multi-state character drives diversification |
| **HiSSE** | binary observed + hidden classes | rate heterogeneity might be off the observed trait |
| **QuaSSE** | one continuous trait | a quantitative character sets a bounded λ(x) |
| **ClaSSE** | discrete/continuous, jumping at nodes | change is concentrated *at* speciation (cladogenetic) |

## The models

### BiSSE

Two states (`0`, `1`), each with its own speciation rate (`lambda0`/`lambda1`) and extinction rate
(`mu0`/`mu1`), plus asymmetric anagenetic transitions (`q01`, `q10`). The classic binary
state-dependent model (Maddison, Midford & Otto 2007); the fast-speciating state comes to dominate
the standing tips. The default `--sse-model`.

### MuSSE

The k-state generalisation: length-`k` `birth` and `death` rate vectors and a `k × k` anagenetic
`Q` matrix (off-diagonals ≥ 0; the diagonal is recomputed so rows sum to zero, exactly as in
[`Mk`](../guide/traits.md)). Use it when a multi-state character — not just a binary one — drives
diversification (FitzJohn 2012, *diversitree*). BiSSE is the `k = 2` special case.

### HiSSE

Extends BiSSE with unobserved **hidden classes**: each observed state comes in `H` variants, each its
own diversification regime, with switch rates between classes (Beaulieu & O'Meara 2016). It is the
honest null for SSE inference — rate heterogeneity that lives on a *hidden* class is not falsely
pinned on the observed character. Build it from one `BiSSE` per hidden class plus a `hidden_transition`
matrix (or a scalar for a symmetric rate); the tips report the **observed** state, with the
`(observed, hidden)` pair still available per node. **Python-only** (no `--sse-model hisse`).

### QuaSSE

A **continuous** trait diffuses (Brownian motion, `sigma2`, optional `drift`) along every lineage and
the rates are functions of its current value (FitzJohn 2010). The rate functions must be **bounded** —
an unbounded λ(x) under a diffusing x has no valid thinning bound — so you pass a `rate_bound` on
λ(x) + μ(x); `QuaSSE.sigmoid(low, high, center, slope)` builds a convenient bounded speciation curve.
On the CLI the trait is a sigmoidal speciation (`--spec-low/high/center/slope`) plus a constant
extinction (`--qmu`).

### ClaSSE

Not a separate class but the **both-arrows** combination: a discrete or continuous SSE model *plus* a
[`Cladogenesis`](../guide/traits.md) kernel that jumps each daughter's state **at** speciation
(Goldberg & Igić 2012). The trait both shapes the tree (`traits:species`) *and* is kicked by its
branching (`species:traits`), so change is concentrated at nodes rather than spread along branches.
`shift` is the per-daughter state-hop probability (discrete); `jump_sigma2` is the Gaussian jump
variance (continuous, `quasse`). With `Q = 0` the cladogenesis kernel supplies all the state dynamics.

## Command line

`--couple traits:species` selects SSE; `--sse-model` picks `bisse` (default), `musse`, or `quasse`.
The run **grows** the tree, so it takes no `-t` — give exactly one stopping condition (`--tips` or
`--age`). Add `--couple species:traits` (with `--clado-shift`/`--clado-jump`) to turn it into ClaSSE.
`--root-state` sets the root state index for `bisse`/`musse` (default: the character's stationary
distribution).

```bash
# BiSSE: state 1 speciates 3x faster, so it comes to dominate the standing tips
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/bisse

# MuSSE: k-state — birth/death vectors + a k x k transition-rate matrix file
printf "0 0.2 0.2\n0.2 0 0.2\n0.2 0.2 0\n" > q3.txt
zombi2 coevolve --couple traits:species --sse-model musse \
    --birth 1 1 3 --death 0.2 0.2 0.2 --q-matrix q3.txt --tips 200 --seed 1 -o out/musse

# QuaSSE: continuous trait — sigmoidal speciation + constant extinction + Brownian diffusion
zombi2 coevolve --couple traits:species --sse-model quasse \
    --spec-low 0.4 --spec-high 3 --spec-center 0 --spec-slope 3 \
    --qmu 0.2 --diffusion 0.5 --root-value -1.5 --tips 200 --seed 1 -o out/quasse

# ClaSSE: both arrows — BiSSE rates + a cladogenetic state hop at each speciation
zombi2 coevolve --couple traits:species --couple species:traits \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.05 --q10 0.05 \
    --clado-shift 0.3 --tips 200 --seed 3 -o out/classe
```

Run the CLI as `python -m zombi2 coevolve ...` (not a bare `zombi2`) if the entry point is not on
your PATH. HiSSE is not exposed on `--sse-model`; use the Python API for it.

## Python

The models live in `zombi2.coevolve`; the cladogenetic kernel in `zombi2.traits` (each also
re-exports at the top level, so `zombi2.BiSSE` / `zombi2.Cladogenesis` work too):

```python
from zombi2.coevolve import BiSSE, MuSSE, HiSSE, QuaSSE, simulate_sse
from zombi2.traits import Cladogenesis

# BiSSE: state 1 speciates 3x faster -> tips biased toward state 1
res = simulate_sse(BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2, q01=0.1, q10=0.1),
                   n_tips=200, seed=1)
res.tree                    # complete tree (extinct lineages kept; z.prune() for the reconstructed one)
res.labeled_values()        # the observed trait at the extant tips

# MuSSE: k-state — birth/death vectors + a k x k Q
Q = [[0, 0.2, 0.2], [0.2, 0, 0.2], [0.2, 0.2, 0]]
musse = simulate_sse(MuSSE(birth=[1, 1, 3], death=[0.2, 0.2, 0.2], Q=Q), n_tips=200, seed=1)

# HiSSE: hidden classes drive the tree while the observed character stays neutral
fast, slow = BiSSE(2.5, 2.5, 0.2, 0.2, 0.3, 0.3), BiSSE(0.4, 0.4, 0.2, 0.2, 0.3, 0.3)
hisse = simulate_sse(HiSSE([fast, slow], hidden_transition=0.15), age=1.5, seed=0)

# QuaSSE: a bounded sigmoidal speciation on a Brownian continuous trait
spec = QuaSSE.sigmoid(low=0.4, high=3.0, center=0.0, slope=3.0)
quasse = simulate_sse(QuaSSE(spec, lambda x: 0.2, sigma2=0.5, rate_bound=3.2, x0=-1.5),
                      age=2.5, seed=0)

# ClaSSE: BiSSE + a cladogenetic jump at each speciation (both arrows)
classe = simulate_sse(BiSSE(1, 3, 0.2, 0.2, 0.05, 0.05),
                      cladogenesis=Cladogenesis(shift=0.3), n_tips=200, seed=3)
```

Provide exactly one stopping condition: `age` (fixed crown age, random tip count) or `n_tips` (grow
until this many extant tips first coexist, random age). The run is conditioned on at least two extant
survivors. `simulate_sse` returns a `TraitResult`: `.tree` is the complete tree (extinct leaves carry
`is_extant=False`), `.values` are the extant tips' states, and `.history` is the realized character
map (discrete models; `None` for QuaSSE).

## Output

The three shared coevolve files, written to the output directory: `species_tree.nwk` (the complete
tree the trait's rates shaped, extinct lineages kept), `traits.tsv` (every node — tips *and*
ancestral states), and `trait_tree.nwk` (the trait annotated on every node). Prune to the
reconstructed, survivors-only tree with `zombi2.prune(result.tree)` for downstream analysis.

## Validation

- **BiSSE** — a state that speciates 3× faster strongly biases the standing tips toward it
  (`test_sse.py::test_sse_faster_speciation_biases_tips`).
- **MuSSE** — with three states sharing extinction and a symmetric transition matrix, the
  fastest-speciating state (3× the others) is over-represented among the standing tips, far above the
  1/3 state-independent baseline
  (`test_sse.py::test_musse_fastest_speciation_state_over_represented_in_tips`).
- **HiSSE** — the fast *hidden* class dominates the tips while the *observed* character stays neutral
  (`test_sse.py::test_hisse_hidden_drives_diversification_not_observed`).
- **QuaSSE** — when speciation rises with the trait, surviving tips are biased to high values versus a
  constant-rate null (`test_sse.py::test_quasse_x_dependent_speciation_biases_the_trait`).
- **ClaSSE** — with anagenetic diffusion switched off, each parent→child step is exactly one
  cladogenetic jump, and those jumps are distributed `Normal(0, jump_sigma2)`: the empirical mean is
  ~0 and the empirical variance matches the `jump_sigma2` parameter to several sigma
  (`test_sse.py::test_classe_continuous_jumps_are_normal_zero_jump_sigma2`).

## References

- Maddison, W. P., Midford, P. E. & Otto, S. P. (2007). Estimating a binary character's effect on
  speciation and extinction. *Systematic Biology* 56(5): 701–710. (BiSSE)
- FitzJohn, R. G. (2010). Quantitative traits and diversification. *Systematic Biology* 59(6):
  619–633. (QuaSSE)
- FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R.
  *Methods in Ecology and Evolution* 3(6): 1084–1092. (MuSSE)
- Beaulieu, J. M. & O'Meara, B. C. (2016). Detecting hidden diversification shifts in models of
  trait-dependent speciation and extinction. *Systematic Biology* 65(4): 583–601. (HiSSE)
- Goldberg, E. E. & Igić, B. (2012). Tempo and mode in plant breeding system evolution. *Evolution*
  66(12): 3701–3709. (ClaSSE — cladogenetic state change + SSE)
</content>
</invoke>
