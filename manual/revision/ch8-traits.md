# Trait evolution

The trait level is where ZOMBI2 leaves genes behind and evolves a **phenotype**: a body size, a habitat, the presence or absence of a structure. A trait rides along a tree the way everything else in the book does, but it is a different sort of thing from the levels before it, and this chapter opens by saying so plainly, because naming the difference is what makes the rest simple. There are two entry points, one for each kind of trait:

```python
from zombi2 import traits
result = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)   # a real value
result = traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                                  switch=0.1, seed=1)                     # a discrete state
```

*[`simulate_continuous` and `simulate_discrete` are built (`zombi2/traits`), following `docs/design/trait-api.md`; the thirteen-class model zoo they replace (`BrownianMotion`, `OrnsteinUhlenbeck`, `EarlyBurst`, `Mk`, `ThresholdModel`, `CorrelatedBinary`, and the rest) is gone. Brownian motion, Ornstein–Uhlenbeck, early burst, variable-rates BM, diversity-dependent σ², the Mk model, the threshold model, and correlated traits (continuous and discrete) all run today. Not yet wired: per-trait modifiers under `correlation=`, multivariate OU, and the deferred cases (regime shifts, cladogenetic jumps, hidden rate classes, DEC → experimental); each divergence is flagged where it arises.]*

## A trait is a different kind of object, and that is fine

The other three levels are **genealogies**. Lineages, gene copies, and sites are things that are *born and lost* along a tree, and you count them as events: a speciation here, a duplication there, a substitution at this site. A trait is not born or lost. It is a **value that rides the tree**, and you observe the value itself, a number or a state at each tip, not a tally of events. So the trait level has no "rate of events" in the way the others do, and there is no "per what?" to answer: a trait is one value per lineage, full stop.

That is a real seam, and it is worth naming rather than papering over. What keeps traits inside the same book is that the *ways* a value evolves are exactly the ways a rate varies one level up. The knobs that bent a substitution rate across the tree in the previous chapter are the knobs that bend a trait here, so once you have read the clocks section you already know most of this one.

## Two entry points, split by state space

The one thing that genuinely differs between traits is the state space: a real-valued measurement behaves nothing like a finite set of states, and the two need genuinely different arguments. So the *kind is the function*, not a flag you pass, the same principle that gave the genome level three functions rather than one with a `--resolution` switch:

```python
traits.simulate_continuous(tree, …)   # real-valued: BM / OU / early burst
traits.simulate_discrete(tree, …)     # finite states: Mk / threshold
```

A trait rides any tree, so `tree` is a species tree (or a gene tree, or any tree you hand it). And because a trait is one value per branch of one tree, there is no per-family multiplicity to manage, no "which of the many trees does this ride" question that the sequence level had to answer. Traits are, in that sense, the *cleanest* level in the book.

## Continuous traits: Brownian motion is the native process

A continuous trait does **Brownian motion** natively. That is not a modifier bolted onto a simpler default; it is what a diffusing value *is*. You give it a starting value and a diffusion rate, and it wanders down every branch, its variance growing in proportion to elapsed time:

```python
# BM — a body size diffusing from 0 at variance-rate σ² = 1.0
traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)
```

Here `rate` is the Brownian variance-rate σ², the trait level's reading of "how fast", and like every rate in the book it accepts modifiers. The two classical variants on Brownian motion are then nothing but the knobs you already met on the molecular clock, one level over:

```python
# OU — the same diffusion, pulled toward an optimum value
traits.simulate_continuous(tree, start=0.0, rate=1.0,
                           reverts_to=2.0, pull=0.5, seed=1)

# early burst — the diffusion rate itself decays through time
traits.simulate_continuous(tree, start=0.0,
                           rate=1.0 * mod.OnTime({0: 1.0, 5: 0.2}), seed=1)
```

The **Ornstein–Uhlenbeck** process is Brownian motion with a rubber band: `reverts_to` is the optimum it is pulled back toward, and `pull` is how hard. These are the same-named `reverts_to` and `pull` that turn the autocorrelated clock into a CIR clock in the previous chapter, a rate that drifts but is held near a mean. **Early burst** (or ACDC) is a diffusion rate that decays as the tree ages, so most of the divergence happens near the root: it is written with the very same `mod.OnTime` modifier that gives the species tree its skyline. Adaptive radiations are the usual motivation, and the shape is one modifier deep.

The point worth pressing is that the unification is at the level of the *knobs*, not a shared wrapper class. `reverts_to`/`pull` and `mod.OnTime` are literally the same knobs used at the species and sequence levels, reused here because a value and a rate answer the same questions about what they remember and how they change. The seam from the opening does show through in the spelling: because a trait *is* a value rather than a rate riding on something else, `reverts_to` and `pull` are direct arguments, while a rate-shaping knob like `OnTime` multiplies the `rate` exactly as it would anywhere. Same knobs, landing in the two natural places.

There is a deeper echo here, and it is exact. The autocorrelated molecular clock of the last chapter was described as a rate doing Brownian motion down the tree. A continuous trait is that same object with the disguise removed: the value doing Brownian motion is now the whole thing you simulate, not a modifier on a substitution rate. The clock was a trait all along.

## Discrete traits: a state switching along the tree

A discrete trait takes a finite set of states and switches between them along the branches, a continuous-time Markov chain, which the field calls the **Mk model**:

```python
# Mk — habitat flips between two states at rate 0.1
traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                         switch=0.1, start="marine", seed=1)
```

This is the literal twin of the discrete molecular clock: `simulate_discrete(states=…, switch=…)` reads the same as `mod.Markov(states=…, switch=…)`, a set of states plus a switching rate. When the flips are not symmetric, replace the single rate with a small matrix of directed rates:

```python
# asymmetric — gains are commoner than losses
traits.simulate_discrete(tree, states=["absent", "present"],
                         switch={"absent->present": 0.2, "present->absent": 0.05},
                         seed=1)
```

A **threshold** trait is the third case, and it is a bridge back to the continuous world. An observed discrete state can be driven by an underlying continuous **liability** that itself does Brownian motion; the state you see is which side of a threshold the liability currently sits on. It is Wright's threshold character, later Felsenstein's, and it earns its place here because it is what makes correlated discrete traits fall out for free below:

```python
# threshold — a discrete state read off an underlying continuous liability
traits.simulate_discrete(tree, states=["absent", "present"],
                         liability=1.0, threshold=0.0, seed=1)
```

## Correlated traits: the joint rule, inside a level

Two traits that evolve independently are two separate calls, in either order. Two traits that drift *together* cannot be simulated one before the other, because each is entangled with the other as it unfolds. That is exactly the **joint** case of Chapter 2, `P(size, limb)` rather than `P(size) · P(limb)`, and the same rule applies: joint means one command. The novelty is only that here the joint pair sits *within* a single level rather than across two, so there is nothing new to teach, just the Chapter 2 rule applied one scale down.

Correlation is specified as **per-trait rates plus a correlation overlay**, not as a full covariance matrix. Each trait keeps its own rate and its own modifiers, so the whole grammar survives intact; the coupling is a separate, dimensionless number in `[−1, 1]` that reads the way people actually think about it:

```python
traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8 * mod.OnTime({0: 1, 5: 0.3})},   # each keeps its modifiers
    correlation={("size", "limb"): 0.6},                          # the overlay, ∈ [−1, 1]
    seed=1)
```

*[This slice wires correlated **Brownian motion** with bare per-trait rates and the `correlation=` overlay (the runnable form is in *Usage from Python* below). Per-trait modifiers under `correlation=`, and multivariate OU, are later slices — the `mod.OnTime` on `limb` above shows the intended surface, not yet-runnable code.]*

The same `correlation=` overlay handles *discrete* correlation with no extra machinery, through the threshold model: give each discrete trait a liability, correlate the liabilities, and put the thresholds on top. Correlated presence/absence characters, the setting Pagel's method was built for, are then one call:

```python
traits.simulate_discrete(tree, states=["absent", "present"],
    liability={"wings": 1.0, "flight": 1.0},
    correlation={("wings", "flight"): 0.7}, threshold=0.0, seed=1)
```

One overlay, working on continuous traits directly and on discrete traits through their liabilities. (A full covariance matrix may still be accepted as an alternative input for readers who arrive from comparative methods already thinking in Σ, but the per-trait-plus-`correlation=` form is the surface the chapter teaches.)

## The literature → command bridge

Trait models arrive under a thicket of names, and a reader who wants "an OU model" or "a threshold model" should be able to find it. The names live here, in one table, and organise nothing else in the chapter.

| Literature name | What it does | ZOMBI2 |
|---|---|---|
| Brownian motion (BM) | a value diffusing | `simulate_continuous(rate=…)` |
| Ornstein–Uhlenbeck (OU) | diffusion pulled to an optimum | `simulate_continuous(rate=…, reverts_to=…, pull=…)` |
| Early burst (EB / ACDC) | diffusion rate decays through time | `simulate_continuous(rate=1.0 * mod.OnTime({…}))` |
| Multivariate BM / OU | traits evolving together | one `simulate_continuous(rate={…}, correlation={…})` call |
| Mk (k-state Markov) | a discrete state switching | `simulate_discrete(states=…, switch=…)` |
| Threshold / liability (Wright–Felsenstein) | discrete driven by continuous liability | `simulate_discrete(liability=…, threshold=…)` |
| Correlated binary / Pagel | discrete traits evolving together | `simulate_discrete(liability={…}, correlation={…})` |
| DEC biogeography | a range = a set of areas | → **experimental** (purged from the trait level for now) |
| BiSSE / MuSSE / QuaSSE / HiSSE | a trait drives speciation | **not a trait model** — trait ↔ species *joint*, Part III |

Two rows point off the chapter and are worth a sentence each. **DEC**, the dispersal–extinction–cladogenesis model of geographic ranges, is a genuine discrete model but a heavy and specialised one; it moves to `zombi2.experimental`, recoverable, rather than cluttering the trait level. The **SSE** family (BiSSE and its relatives) is the more important cut: it is *not* a trait model at all. A trait that drives speciation reaches back into the tree it rides on, so the tree becomes an output grown together with the trait. That is the joint case where the coupling crosses into the species level, and it lives in Part III (Chapter 10), not here.

## Still to design

A couple of corners of the trait level are now settled; the rest are agreed in concept but not yet nailed down, and the chapter should not pretend otherwise.

- **Decided: the OU trait is *not* a `FromParent` modifier; the mechanism differs, only the names are shared.** OU still needs both an optimum (`reverts_to`) and a pull strength (`pull`), but these stay **function arguments** that revert the trait's *value* continuously along a branch (`θ + (x−θ)·e^{−α·dt}`), which cannot be written as a multiplier on σ². The same-named `reverts_to`/`pull` reappear one level over as the sequences-side **CIR clock**, `mod.FromParent(spread=, reverts_to=, pull=)`: a mean-reverting *rate*, where plain `spread` is pure drift (variable-rates BM here, ClaDS at the species level, the autocorrelated clock) and adding `reverts_to` + `pull` makes that rate mean-reverting. Shared vocabulary, different mechanism (a value versus a rate); the names `reverts_to`/`pull` stand.
- **Regime shifts (multi-optimum OU).** An OU optimum that jumps on painted branches, so different clades pull toward different values, is the natural home of adaptive-regime studies. It is an advanced case, probably a `regimes=` argument on `simulate_continuous`; deferred, and named here so it is not mistaken for missing by oversight.
- **Traits that jump at speciation (cladogenesis).** A trait can change *at a split* rather than along a branch. By the spec's own reasoning this is the trait reading the tree it already lives on, an option of the trait's own model rather than a coupling to the species level, so it belongs in this chapter; the open question is only its spelling (likely `at_speciation=`). It is not built.
- **Hidden rate classes under an Mk trait.** The discrete twin of the clock's hidden categories, letting the switching rate itself vary invisibly across the tree. Likely a hidden-state option on `simulate_discrete`; deferred.
- **Decided: `switch=` accepts both forms.** The `"a->b"` string-keyed dict written above (readable, few states) and a numeric matrix paired with `states=[...]` (many states) are both valid, and `switch=0.1` stays the symmetric shortcut.

## The objects

*[`TraitsResult` is built (`zombi2/traits`). It does not yet share a common result spine with the other levels — the accessors below are its own; the cross-level spine of `docs/design/result-api.md` is still to come.]*

A run returns a **`TraitsResult`** bundle: the tree it ran on (`.complete_tree`), the `.seed`, the `.kind` (`"continuous"` or `"discrete"`), and `.write(directory, outputs=[...])` to materialise the chosen outputs to disk — plus a payload that follows directly from what a trait is. Because a trait is a value at *every* node, the result records the value at every node, so the ancestral states are not a separate reconstruction step but a byproduct of the simulation: they are exact, drawn from the same process that produced the tips, not inferred after the fact. This is the trait seam again: the value at every node (`.node_values`) *is* the source of truth here, held directly rather than replayed from an event log — for a continuous trait there is no event log to replay.

- `.values` — the observable vector: the trait's value at each **extant tip**. This is the comparative-data matrix a method would be handed.
- `.node_values` — the value at **every** node (extant, extinct, and internal alike), the true ancestors at each split, from the same process that produced the tips.
- `.history` — for a **discrete** trait, the per-branch stochastic character map, the ordered list of `(state, duration)` segments each branch passed through; `.events` reads off the individual transitions with their times. Both are empty/`None` for continuous traits, which have no jumps to record.

For discrete traits the stored values are the state labels you gave (not integer indices), so `.values` and `.node_values` already read back in your own vocabulary. Writing the result as an annotated tree — a **trait tree**, each node carrying its value — is a forthcoming output.

## Usage from Python

```python
from zombi2 import species, traits
from zombi2.rates import modifiers as mod

# a species tree from the previous chapters, then a trait riding along it
tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=1).extant_tree

# continuous: body size under Brownian motion
size = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)
size.values                 # {extant tip: float}

# continuous: an Ornstein–Uhlenbeck trait pulled toward an optimum of 2
temp = traits.simulate_continuous(tree, start=0.0, rate=1.0,
                                  reverts_to=2.0, pull=0.5, seed=1)

# discrete: habitat flipping between two states
habitat = traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                                   switch=0.1, start="marine", seed=1)
habitat.values              # {extant tip: "marine" | "terrestrial"}
habitat.events              # the realized flips, in time order

# two continuous traits that drift together — one joint call
bodyplan = traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8},
    correlation={("size", "limb"): 0.6}, seed=1)
```

## Usage from the CLI

*[Draft — the trait CLI is not built in the clean core yet; the commands below show the intended surface, not a shipped interface.]*

```bash
# a continuous (Brownian) trait along a species tree
zombi2 traits --continuous --tree species_tree.nwk \
    --start 0.0 --rate 1.0 --seed 1 -o my_trait

# a discrete two-state trait
zombi2 traits --discrete --tree species_tree.nwk \
    --states marine,terrestrial --switch 0.1 --seed 1 -o my_habitat
```

## Outputs

A run writes the **trait values** at the extant tips (`trait_values.tsv`, the observable comparative-data vector). For a discrete trait it also writes the **change history** (`trait_changes.tsv`, the realized transitions along each branch, read off the stochastic character map), the ground truth against which an ancestral-state or stochastic-mapping method would be scored. The exact ancestral values at every internal node are kept in `.node_values`, and an annotated **trait tree** — the input tree carrying the value at each node — is a forthcoming output. The full list of files lives in Appendix B.
