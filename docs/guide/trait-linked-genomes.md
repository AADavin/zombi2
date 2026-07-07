# Trait-linked gene families

Gene families and phenotypic traits do not evolve independently. A lineage that becomes
aerobic retains and acquires oxygen-using gene families; one that reverts to anaeroby sheds
them. ZOMBI2 can simulate that link directly: evolve a trait down the tree, then evolve a
panel of gene families whose **loss and gain depend on the local trait value**. The resulting
phylogenetic profile carries a *known*, trait-linked signal — the forward generator behind
studies that read gene content as a record of a trait's history (e.g. timing the bacterial
tree from the Great Oxidation Event, Davin 2025).

```python
from zombi2.species import simulate_species_tree, BirthDeath
from zombi2.traits import Mk
from zombi2.coevolve import TraitGeneCoupling, simulate_trait_linked_genomes

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=60, age=6.0, seed=1)

# a binary aerobic(1)/anaerobic(0) trait, then genes conditioned on it
coupling = TraitGeneCoupling.build(n_families=40, responsive=0.3, weight=1.0,
                                   effect_loss=3.0, base_loss=0.5, transfer=1.0, seed=1)
res = simulate_trait_linked_genomes(tree, Mk.equal_rates(2, 0.4), coupling, seed=2)

res.profiles.presence()        # panel families × extant species (0/1) — the trait-linked data
res.trait.labeled_values()     # the trait at the tips, from the same run
```

This is the *genes-conditioned-on-a-trait* direction. It reuses the trait models of
[Trait evolution](traits.md) and the coupled rate machinery of the
[gene-family coupling model](../coupling_model.md); only the family-side rate model is new, so
the whole output pipeline (profiles, gene trees, reconciliations) applies unchanged.

## The model

A fixed **panel** of `N` gene families is seeded present at the root. Each family carries a
coupling **weight** `w_i` (`0` = inert). With the trait value on a branch at time `t` written
`s`, a *present* family is lost at rate

```
loss_i = base_loss · exp(-effect_loss · w_i · s)
```

so where the trait favours a responsive family (`w_i·s` large and positive) it is retained,
and where it does not (`w_i·s` negative) it is purged faster than the baseline. Inert families
(`w_i = 0`) always lose at `base_loss`.

**Gain is horizontal transfer** — a field-blind influx, exactly as in the
[coupling model](../coupling_model.md): a family flows into a lineage at a constant rate, and
the trait-modulated *loss* then selectively retains it — kept where the trait favours it,
purged where it does not. So the **net** gene content of a lineage tracks its trait even though
the influx itself is trait-blind. That differential retention is what writes the trait↔gene
association into the profiles.

!!! note "Why retention, not a gain switch"
    Coupling through *loss* is the mechanism that produces a clean, datable signal (it is the
    validated device of the Potts coupling model). `effect_gain` optionally scales a lineage's
    transfer (HGT) *activity* by `exp(effect_gain · s)`, but it is a donor-side effect and is
    **off by default** — the retention channel already makes net gene content track the trait.

## Choosing the responsive families

`TraitGeneCoupling.build(n_families, responsive, ...)` populates the weight vector. The
`responsive` selector is the flexible part:

```python
from zombi2.coevolve import TraitGeneCoupling

TraitGeneCoupling.build(50, 8)                       # 8 families, chosen at random
TraitGeneCoupling.build(50, 0.3)                     # a random 30% of the panel
TraitGeneCoupling.build(50, ["F3", "F7", 12])        # exactly these families (id or index)
TraitGeneCoupling.build(50, 10, signed=True)         # half favoured by a high trait value,
                                                     # half by a low one
```

`weight` sets each responsive family's magnitude; `signed=True` randomises its sign so some
families co-occur with a high trait value and others with a low one. `effect_loss` is the
overall coupling strength (`0` recovers plain, uncoupled gene-family evolution). The remaining
rate parameters — `base_loss`, `transfer`, `duplication`, `origination` — are the panel's base
DTL rates.

## The trait as a covariate in time

The trait value varies *along* each branch, and the simulation follows it exactly:

- **Discrete traits** (`Mk`, threshold, …) contribute their exact *stochastic character map* —
  the per-branch `(state, duration)` segments — so a mid-branch state change is honoured to the
  instant it happens (it becomes a rate-refresh point in the Gillespie loop).
- **Continuous traits** (`BrownianMotion`, `OrnsteinUhlenbeck`, …) are sub-segmented into
  `trait_steps` pieces per branch (default 16), with the value interpolated between the node
  endpoints and held constant across each piece.

```python
from zombi2.traits import BrownianMotion
from zombi2.coevolve import simulate_trait_linked_genomes

simulate_trait_linked_genomes(tree, BrownianMotion(0.6), coupling, trait_steps=24, seed=1)
```

For a binary trait it is usually best to **center** the two states around zero
(`state_values=[-1.0, 1.0]`), so the trait pushes a responsive family's retention *up* in one
state and *down* in the other — a symmetric, two-sided coupling — rather than only lowering
loss in the "on" state:

```python
from zombi2.coevolve import TraitGeneCoupling

coupling = TraitGeneCoupling.build(40, 0.3, weight=1.0, effect_loss=3.0,
                                   base_loss=0.5, transfer=1.0,
                                   state_values=[-1.0, 1.0], seed=1)
```

## Reusing an already-simulated trait

`simulate_trait_linked_genomes` accepts either a trait **model** (evolved for you) or an
already-simulated `TraitResult`, so you can inspect or reuse the exact trait the genes were
conditioned on:

```python
from zombi2.traits import simulate_traits, Mk
from zombi2.coevolve import simulate_trait_linked_genomes

trait = simulate_traits(tree, Mk.equal_rates(2, 0.4), seed=2)
res = simulate_trait_linked_genomes(tree, trait, coupling, seed=3)
assert res.trait is trait
```

## The result

`simulate_trait_linked_genomes` returns a `TraitLinkedResult`:

| Access | Meaning |
| --- | --- |
| `res.profiles` | the `N × extant-species` panel `ProfileMatrix` (every panel row kept, even all-absent ones) |
| `res.trait` | the `TraitResult` the genes were conditioned on |
| `res.leaf_genomes` / `res.event_log` | the raw gene-family state and event log |
| `res.coupling` | the `TraitGeneCoupling` used (weights + effect sizes) |
| `res.genomes()` | promote to a standard [`Genomes`](../reference/api.md#simulation-driver) for gene trees, reconciliations and `write()` |

## From the command line

`zombi2 coevolve --couple traits:genes` runs the whole thing on a species tree you provide. It
simulates the trait (`--trait-model`, reusing every `zombi2 trait` model), builds the coupling
(`--panel`, `--responsive`, `--weight`, `--effect-loss`), and writes the gene-family output
plus the trait and a coupling manifest:

```bash
T=out/species_tree.nwk

# a binary aerobic/anaerobic trait; 30% of a 40-family panel respond to it
zombi2 coevolve --couple traits:genes -t $T \
    --trait-model mk --states 2 --rate 0.3 --trait-center \
    --panel 40 --responsive 0.3 --weight 1 --effect-loss 3 \
    --loss 0.4 --trans 1.0 --write all --seed 7 -o out/
```

Besides the usual gene-family files (chosen with `--write`, exactly as in
[`genomes`](../cli.md#choosing-the-output-and-the-rust-engine)), it always writes:

- **`traits.tsv`** / **`trait_tree.nwk`** — the trait at every node (as `zombi2 trait` writes).
- **`coupling.tsv`** — the per-family coupling weights and the effect sizes, so the exact
  trait↔gene linkage that generated the profiles is on record for downstream inference.

Useful options:

- `--responsive` takes a count (`8`), a fraction (`0.3`), an id/index list (`F3,F7,12`), or
  `@file` of ids; `--signed` randomises the weight signs.
- `--trait-center` centers a discrete trait's states (recommended for binary characters).
- `--trait-steps K` sets the within-branch resolution for a continuous trait.
- `--trait-file traits.tsv` reuses a precomputed trait instead of simulating one — a
  `node`/`value` table over **every** node (numeric values; encode discrete states as numbers),
  as `zombi2 trait` writes with its all-nodes output.
- `--effect-gain` turns on the optional donor-side HGT-activity coupling.

See the [CLI reference](../cli.md#coevolve-couple-traitsgenes-trait-conditioned-gene-families)
for the full option table.

## What it recovers

Inject a strong coupling and the trait shows up in the profiles: responsive families are
present where the trait favours them and absent where it does not, while inert families do not
distinguish the states. Concretely, with a two-clade trait (half the tips aerobic, half
anaerobic) and `effect_loss = 3`, responsive families sit at ~0.7 prevalence in the aerobic
tips and ~0.1 in the anaerobic ones, whereas inert families are indistinguishable between the
two — the signal is entirely in the responsive panel, which is exactly what a downstream
inference should be able to pick out.

Keep `base_loss` moderate relative to `transfer` so the *inert* families persist as a control:
with an over-large `base_loss` an unprotected family, having only the field-blind influx to
regain it, is lost tree-wide and the inert rows go all-zero.

## Roadmap

`coevolve --couple traits:genes` is the **`traits:genes`** edge of a broader coupled-simulation
design — see [Coevolution (coupled models)](../coevolution_models.md), which generalises it to a
directed graph over species, traits and gene families (`zombi2 coevolve --couple driver:target`).
It was formerly the standalone `coevolve-genetrait` command, now folded into `coevolve`.
Planned next:

- an **environmental clock** — a trait (and its coupled families) gated by a dated event such
  as the GOE, which is what turns the coupled dynamics into a *time* signal;
- a **recipient-side gain** channel (trait-dependent acquisition, not only retention);
- the into-species edges (`traits:species` = SSE, `genes:species`) that couple traits and gene
  families *back* to the diversification process, up to the fully joint `--all` model.
