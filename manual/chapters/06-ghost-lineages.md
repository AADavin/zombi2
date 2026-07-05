# Ghost lineages

A reconstructed species tree contains only the sampled, extant lineages — that is the whole point
of the backward, conditioned birth–death process behind `simulate_species_tree`. But the real
diversification process also produced lineages that went extinct before the present, or that are
alive today yet were never sampled. These *dead* lineages are invisible in the reconstructed tree,
and dropping them is not harmless. A gene can be laterally transferred from a lineage that later
vanished — a transfer from the dead [@szollosi2013lgtdead]. If the simulated tree carries only
survivors, every transfer is forced to run between sampled, co-existing lineages, which biases any
benchmark of transfer inference.

`add_ghost_lineages` puts the dead lineages back. It *un-prunes* the reconstructed tree: it grafts
the extinct and unsampled ("ghost") lineages back on, so the tree again reflects the full process
rather than only its survivors. The result is an ordinary `Tree`; the extra tips are just leaves
marked `is_extant=False`.

![The reconstructed tree of survivors (solid) with the extinct and unsampled ghost lineages grafted back on (dashed), recovering the full diversification history.](figures/model_ghosts.pdf){width=100%}

## Un-pruning versus simulating forward

There are two ways to obtain a tree that includes dead lineages, and they serve different starting
points.

::: note
Un-pruning is for when you already have — or want — a reconstructed tree conditioned on $N$ extant
tips, and need the dead lineages added back. If instead you want the complete tree from the start,
`simulate_species_tree(..., direction="forward")` grows it forward with extinct lineages included
natively, and no un-pruning is needed. The two are distributionally equivalent for the
reconstructed part; pick by which tree you have in hand.
:::

## Basic use

Pass the *same model* you used to build the tree, so the ghost process runs at the matching rates.
The function mutates the tree in place and also returns it.

```python
model = z.BirthDeath(birth=1.0, death=0.5)
tree = z.simulate_species_tree(model, n_tips=50, age=5.0, seed=1)

z.add_ghost_lineages(tree, model, seed=7)   # grafts ghosts in place
```

Each grafted ghost attaches along an edge and roots a birth–death subtree that leaves no sampled
descendant. Every new node is marked `is_extant=False`; the dead leaves are named `e*` and the new
internal nodes `i*`, following ZOMBI2's convention (extant leaves stay `n*`). The sampled leaves
are left untouched, so pruning back to the extant tips recovers the original reconstructed tree
exactly.

```python
extant = [n for n in tree.leaves() if n.is_extant]      # sampled tips (n*)
ghosts = [n for n in tree.leaves() if not n.is_extant]  # dead tips (e*)
print(len(extant), "extant +", len(ghosts), "ghost leaves")
```

Ghosts flow straight into the forward gene simulation with no extra step: they are just branches,
so speciation gives each a genome and transfers can now draw ghost donors and recipients
automatically. Because gene trees still prune to sampled tips, a gene from a ghost donor surfaces
as a transfer from a lineage absent from the species tree — exactly the transfer-from-the-dead
case.

On the command line, add ghosts with the `--ghosts` flag, which un-prunes the tree before the gene
layer runs (`--ghost-method` picks the sampler):

```bash
zombi2 species --birth 1.0 --death 0.5 --tips 50 --age 5.0 \
    --ghosts --seed 7 -o run/
```

## When do ghosts appear?

Ghosts appear only where lineages could actually have been lost — the process must allow extinction
or incomplete sampling.

| Model | Ghosts? |
| --- | --- |
| `Yule` / `BirthDeath(death=0)`, full sampling | none — nothing goes extinct |
| `BirthDeath(death>0)` | yes — extinct lineages grafted back |
| `EpisodicBirthDeath(..., sampling_fraction<1)` | yes — extinct *and* unsampled lineages |

Extinct lineages need $\mu > 0$; unsampled-extant lineages need $\rho < 1$. Both reduce to "leaves
no sampled descendant," so a tree with no possible loss returns from `add_ghost_lineages`
unchanged. The supported models are `BirthDeath` / `Yule` and `EpisodicBirthDeath`, the latter with
time-varying rates and incomplete sampling.

## The conditional law, briefly

Un-pruning is exact because it uses the known conditional law of the complete birth–death tree
given its reconstructed version [@nee1994reconstructed; @stadler2009incomplete;
@lambert2013birthdeath]. Along each edge of the reconstructed tree, dead lineages attach as an
inhomogeneous Poisson process with attachment intensity proportional to $\lambda(t)\,E(t)$, where
$E(t)$ is the probability that a lineage present at time $t$ leaves no sampled descendant at the
present.

The intuition is short. A surviving lineage still speciates at rate $\lambda$; each side-branch
independently leaves a sampled descendant (probability $1 - E$) or nothing (probability $E$). The
reconstructed tree kept only the branchings where both sides survived — those are its nodes. So
*between* its nodes, every branching had a dying sibling, and those attach at rate
$\lambda\,E(t)$. The quantity $E(t)$ is the same survival probability ZOMBI2 already solves for the
birth–death process, with the boundary condition $E(0) = 1-\rho$ — which is why the extinct and
unsampled cases fall out of one equation.

Each attachment then roots a subtree grown forward and conditioned on leaving no sampled
descendant. Two equivalent samplers are available through `method=`:

- `"rejection"` (default) — grow an ordinary birth–death subtree from the attachment time and
  reject any that leaves a survivor. Simple and provably correct; ghosts are small because they are
  conditioned to die out, so the work stays cheap.
- `"htransform"` — rejection-free, via Doob's h-transform. The conditioned process is itself a
  birth–death with per-lineage birth $\lambda E$ and death $\mu/E$, so each subtree is drawn in a
  single pass. It is statistically equivalent to rejection and faster when heavy extinction or
  incomplete sampling would otherwise force many retries.

```python
z.add_ghost_lineages(tree, model, method="htransform", seed=7)
```

::: tip
Pruning the ghosts back off returns the original reconstructed tree unchanged, since ghosts have no
sampled descendants. This is an exact augmentation, not an approximation: the sampled-tree
statistics — the set and times of the extant tips — are untouched.
:::
