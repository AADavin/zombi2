```{=latex}
\part{Species trees}
```

# Species trees (basic models)

Every ZOMBI2 simulation has a species tree: the branching history of the species onto
which gene families, traits, and sequences are later layered. This chapter covers the
constant-rate birth–death process and its pure-birth (Yule) special case, how to condition a
simulation on the number of tips and on the tree's age, the two ways ZOMBI2 can generate a tree
— backward (the reconstructed process) and forward (the complete process, keeping extinct
lineages) — and the `Tree` object the simulator returns.

## The birth–death model

Under the constant-rate birth–death process, each lineage independently speciates at rate
$\lambda$ (birth) and goes extinct at rate $\mu$ (death). A pure-birth process, in which lineages
never die ($\mu = 0$), is the **Yule** model. Because extinction removes lineages, a tree grown
forward in time contains two kinds of tips: *extant* leaves that survive to the present, and
*extinct* leaves that died along the way. The tree obtained by pruning away the extinct lineages
— keeping only the branches ancestral to the survivors — is the **reconstructed tree**
[@nee1994reconstructed]. It is what one recovers from present-day data, and it is what ZOMBI2
produces by default.

The reconstructed birth–death tree is a *coalescent point process*: conditional on the number of
extant tips and the tree age, its internal-node ages are independent and identically distributed
draws from a known distribution [@hartmann2010sampling]. ZOMBI2 exploits this directly. It draws
the internal-node ages i.i.d. from the reconstructed-process CDF and assembles a ranked tree by
uniform coalescence [@hartmann2010sampling]. No lineages are simulated and discarded, so
producing a tree with a prescribed number of tips is exact and fast, however small the survival
probability.

The two models are constructed by naming their rates:

```python
from zombi2.species import BirthDeath, Yule, simulate_species_tree

BirthDeath(birth=1.0, death=0.3)   # speciation lambda, extinction mu
Yule(birth=1.0)                    # pure birth == BirthDeath(birth, death=0)
```

## Simulating a tree

`simulate_species_tree` takes a model and the conditioning that pins down its size and timescale:

```python
tree = simulate_species_tree(
    BirthDeath(1.0, 0.3),
    n_tips=20,          # condition on the number of extant species (>= 2)
    age=5.0,            # tree age
    age_type="crown",   # "crown": age of the root; "stem": time of origin
    seed=1,             # or rng=<numpy Generator>
)
```

The two knobs fix the two quantities the reconstructed process needs:

- **`n_tips`** — the tree has exactly this many extant leaves. The conditioning is exact, not the
  result of rejection sampling, so it is cheap even when extinction is high.
- **`age`** and **`age_type`** — with `age_type="crown"` the root sits at time 0 and every extant
  leaf sits at `age`; the tree is ultrametric. With `age_type="stem"` the age is instead the
  *origin* time, and a stem branch precedes the crown (the first speciation)
  [@stadler2009incomplete].

![What `age` measures. With `age_type="crown"` (left) the age is the depth from the crown — the root of the reconstructed tree — to the present. With `age_type="stem"` (right) it is measured from the origin instead, so a stem branch precedes the crown and the crown subtree is correspondingly shorter.](figures/age_crown.pdf){width=100%}

## Backward versus forward simulation

By default ZOMBI2 simulates the *reconstructed* tree backward in time, sampling node ages from
the coalescent point process as described above. This is the right object when you want a tree
that looks like an inferred phylogeny of extant species: it contains no extinct lineages, because
none survive to be observed.

Sometimes, though, the extinct lineages matter. Gene families evolving on the tree can be
transferred *from* a since-extinct donor, so the dead branches leave a genomic signature even
though they contribute no tips. For this ZOMBI2 can instead grow the **complete** tree forward in
time, retaining every lineage — extant and extinct alike:

```python
# grow the full tree forward, keeping extinct lineages
tree = simulate_species_tree(
    BirthDeath(1.0, 0.3),
    age=5.0,
    direction="forward",
    seed=1,
)
```

A forward run realizes the process by an exact event-by-event (Gillespie) loop — the engine
described in full in Appendix A: starting from the root it draws waiting times to the next birth
or death and applies them until it reaches `age`.
The result keeps the extinct leaves — conventionally named `e*` — so a lineage that died before
the present is still a tip of the complete tree. The run is conditioned to leave at least two
survivors; a realization in which every lineage dies is rejected and redrawn.

![A complete forward tree that keeps its extinct lineages, drawn dashed and named `e1, e2, …`.](figures/species_tree_extinct.pdf){width=100%}

The complete-versus-reconstructed distinction is the same one drawn analytically for these
processes: the reconstructed tree is the complete tree with its extinct subtrees pruned away
[@nee1994reconstructed; @lambert2013birthdeath].

## The `Tree` object

Both modes return the same `Tree` object. Its nodes carry times, and it exposes the traversals
the rest of ZOMBI2 needs:

```python
tree.to_newick()          # timed Newick (branch lengths from node times)
tree.leaves()             # extant leaves
tree.internal_nodes()
tree.branches_alive_at(t) # lineages crossing time t (used by the gene-family loop)
tree.total_age
```

Node times increase forward from the root at time 0 to the extant leaves at `total_age`. A branch
is identified by its child node and spans the half-open interval `(parent.time, node.time]`. In a
backward (reconstructed) tree every leaf sits at `total_age` and the tree is ultrametric; in a
forward tree the extinct leaves sit at earlier times, wherever their lineages died.

The `branches_alive_at(t)` method returns the set of lineages crossing a given time, which is
exactly what the forward gene-family simulator iterates over as it walks the species tree. Because
the timed Newick encodes node times as branch lengths, the tree round-trips through standard
phylogenetics tooling.

Beyond the constant-rate birth–death and Yule models, ZOMBI2 offers a family of richer
diversification models — episodic (skyline) rates, instantaneous mass extinctions, per-lineage
ClaDS rates, diversity-dependent diversification, and scheduled clade-specific shifts. Those are
the subject of the next chapter.
