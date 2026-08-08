```{=latex}
\appendix
```

# Rates in detail, and the Gillespie algorithm

Chapter 2 introduced the shape every rate takes, `effective rate = scope(base) × modifiers`. This
appendix is the full reference: how a rate's units work, the default scope at each level, the catalogue
of modifiers and which levels accept them, and the Gillespie algorithm that turns rates into events.

## How a rate is counted: the scope

A rate always has units of time⁻¹, on the scale imposed by the species tree. In a phylogenetic context,
though, a single global rate rarely makes sense for most events. A substitution happens at a **site**, so
a mutation rate is counted per site (mutations × time⁻¹ × per site): each site is an independent chance to
mutate. A speciation happens to a **lineage**, so the speciation rate is counted per lineage (speciations
× time⁻¹ × per lineage): each branch alive is an independent chance for the tree to split. And a gene is
lost one gene copy at a time, so gene loss is counted per copy (loss × time⁻¹ × per gene-copy). The unit a
rate is counted in, whether per lineage, per copy or per site, is what we call its **scope**.

By default, this is the scope ZOMBI2 uses at each level:

| Level | Counted per | The rates it applies to |
|---|---|---|
| Species | lineage | `birth`, `death` |
| Genomes, family and ordered | copy | `duplication`, `transfer`, `loss` |
| Genomes, family and ordered | lineage | `origination` |
| Genomes, ordered | copy | `inversion`, `transposition`, `translocation` |
| Genomes, ordered | chromosome | `fission`, `fusion`, `chromosome_loss` |
| Genomes, ordered | lineage | `chromosome_origination` |
| Genomes, nucleotide | lineage | `duplication`, `transfer`, `loss`, `origination`, `inversion`, `transposition`, `translocation`, `chromosome_origination` |
| Genomes, nucleotide | chromosome | `fission`, `fusion`, `chromosome_loss` |
| Sequences | site | `substitution` (times a clock) |
| Traits | lineage | `rate` (continuous), `switch` (discrete) |

A bare number takes the default. Writing the default wrapper explicitly is accepted at every level
but one: a discrete trait's `switch` takes a number, a `{'from->to': rate}` dict or a matrix, and a
bare `PerLineage(0.5)` is none of those.

Three levels accept more than their default, and each extra scope is a different model rather than a
different spelling:

| Rate | Also accepts | What that means |
|---|---|---|
| Species `birth`, `death` | `Global` | one budget for the whole tree: linear growth, not exponential |
| Genomes, family and ordered: `duplication`, `transfer`, `loss` | `PerLineage` | a fixed budget per lineage — the rate is the same however much the genome holds |
| Genomes, ordered: `inversion`, `transposition`, `translocation` | `PerLineage` | the same, for the rearrangements |

The difference is large and worth stating plainly. `loss = PerCopy(0.25)` puts every copy
independently at risk, so a genome of a thousand genes loses ten times as often as one of a hundred.
`loss = PerLineage(0.25)` is a deletion budget: the lineage loses at 0.25 whatever it holds, and the
genome's size never enters. The same number, a hundredfold different model.

`origination` stays per lineage — it is the rate at which *new* families arrive, so per copy it would
be zero in an empty genome — and the chromosome tier stays per chromosome. The nucleotide resolution
already counts its gene events per lineage, which is why it has no second scope to offer: there the
rate says how often a lineage does the thing and the extent says how much DNA it touches.

```python
from zombi2 import species
from zombi2.rates import scope

# a death rate applied to the whole tree at once, not once per lineage
species.simulate_species_tree(birth=1.0, death=scope.Global(0.3), total_time=8.0, seed=2)
```

`Global(x)` detaches the rate from the amount of material present, so a `Global` death rate does not
grow as the tree does. The wrappers are `Global`, `PerLineage`, `PerCopy`, `PerSite` and
`PerChromosome`. Everywhere else — and for every scope but `Global` at the species level — a rate
handed a scope other than the one in the table above is refused, and the error names the scope that
level takes — except on a discrete trait's `switch`, where a bare wrapper is turned away by the shape
check instead and the message is about the shapes `switch` accepts. No other scope is implemented
anywhere.

## Bending a rate: modifiers

A **modifier** alters a rate in context. You might give a gene family a constant loss rate across the
tree except in one clade known to shed genes, a symbiotic bacterium, by multiplying the rate there.
Or let families evolve at different speeds: an antimicrobial-resistance family prone to transfer, a
ribosomal-protein family the opposite.

The modifiers are:

| Modifier | What it does to the rate |
|---|---|
| `OnTime` | Follows a **time schedule**: one factor up to a breakpoint, another after it. |
| `OnTotalDiversity` | **Slows as the tree fills up**: the factor falls from 1 toward 0 as the number of lineages approaches a carrying capacity, and stays there. |
| `Inherited(per='lineage')` | Is **inherited from the parent lineage and nudged at each split**, so the rate drifts gradually down the tree and close relatives keep similar rates. |
| `Drawn(per='lineage')` | Is an **independent draw for each lineage**, with no memory of its parent, so nearby branches are no more alike than distant ones. |
| `Drawn(per='family')` | Is an **independent draw for each gene family**, so one family is prone to transfer and another is not, whatever lineage either sits in. |
| `ScaledBy` | **Reads an evolved value**: the factor is looked up from a driver's state, either another level or another run of the same one, which is how one thing conditions another (Chapter 9). |
| `SetBy` | **Replaces the base** instead of multiplying it: the mapping gives the rate itself, in the rate's own units, so nothing is written in front of it. It reads a driver just as `ScaledBy` does, the scope still applies, and a rate carries one. |

The first two are **deterministic**: `OnTime` and `OnTotalDiversity` are fixed functions of the state of
the world, so every lineage that meets the same time, or the same diversity, gets the same factor. The
next two are **random and vary from lineage to lineage**, and they differ in *memory*: `Inherited(per='lineage')` is
passed down and drifts, so the rate is autocorrelated along the tree, a slowly wandering clock or a
clade that inherits a fast tempo, whereas `Drawn(per='lineage')` is drawn afresh on every branch, so the variation
is scattered, an uncorrelated ("relaxed") clock. `Drawn(per='family')` is the same independent draw made over gene
families rather than lineages, so it varies what a family does wherever it sits. All three random
modifiers are **mean-corrected**, meaning their factors average to 1, so a lineage is no likelier to be
sped up than slowed down.

Be careful about what that does and does not buy you. It fixes the *factor*, not the *tree*. A birth
rate that drifts is multiplicative, and a branching process is convex in its rate: the fast lineages
branch, and their descendants inherit the fast tempo, so they come to dominate the tree while the slow
ones contribute almost nothing. Standing diversity therefore **rises** as you turn `spread` up, even
though every individual factor averages to 1, and at moderate spread a run can grow explosively enough
to hit the `max_lineages` guard. Mean-correcting keeps the rate honest per lineage; it does not hold
`E[N(t)]` fixed, and nothing could. A driven factor is neither random nor corrected: it is whatever
the driver's state says it is. A driver need not come from another level: `Clade({"fast": ["n12",
"n27"]})` names a subtree by its tips, or by a node id, and reads membership off the tree the run is
already walking, with every unnamed lineage in `"rest"`. A driver is not confined to a rate: the same
factor multiplies an **extent** at the ordered and nucleotide resolutions. An extent takes only
`OnTime` and `ScaledBy` — a per-family draw has no one family to read (a run covers several), and
`SetBy` has no base to replace. Drivers also go on `transfer_to`, written `Weights` rather than
`ScaledBy` because there the numbers are normalised **weights** over the candidate recipients,
compared against each other, with no base in front of them (Chapter 9).

**How wide the variation is, and what shape it takes.** A drawn or inherited value takes `spread=σ`,
which means a lognormal of that log-scale, and the number it means is different in the two: under
`Drawn` it is the spread of the values themselves, and under `Inherited` it is the step taken at each
split, which accumulates down the tree. `Drawn` also takes `dist=` instead — any built-in
distribution, so `Drawn(per='family', dist=Gamma(shape=4.0, scale=0.25))` gives the draws a gamma's
shape. Whatever the distribution, the draw is divided by its own mean, so what a distribution
contributes is its shape and not where it sits; give a spread or a dist, never both.

`Inherited` takes one more, `bins=`. With it the value moves in **steps** rather than continuously:
it takes one of `bins` values on a log-spaced ladder, and a daughter moves to a neighbouring rung or
stays. `Inherited(per='lineage', spread=0.45, bins=6)` is the discrete-bin clock, where a clade's
rate is one of a handful rather than anything at all.

Modifiers **stack by multiplication**, so they combine: `1.0 * mod.OnTime({0: 1, 5: 0.3}) *
mod.Inherited(per='lineage', spread=0.3)` is a rate that both follows a schedule and drifts between lineages.

### Which level accepts which

A modifier only makes sense where the level can act on it, and a level **rejects** one it does not
accept rather than silently ignoring it. `SetBy` is listed separately from `ScaledBy` even though both
read a driver, because replacing a base is a capability an engine has or has not: the three levels
below can do it and the rest refuse. This is what each accepts today:

| Level | The modifiers it accepts |
|---|---|
| Species | `OnTime` · `OnTotalDiversity` · `inherited per lineage` · `drawn per lineage` |
| Genomes, family and ordered | `OnTime` · `ScaledBy` · `SetBy` · `drawn per family` |
| Genomes, nucleotide | `OnTime` · `ScaledBy` |
| Sequences | `drawn per lineage` · `inherited per lineage` · `ScaledBy` |
| Traits, continuous `rate` | `OnTime` · `inherited per lineage` · `OnTotalDiversity` · `ScaledBy` · `SetBy` |
| Traits, discrete `switch` | `ScaledBy` |
| Joint, `birth` / `death` | `OnTime` · `OnTotalDiversity` · `ScaledBy` |

`zombi2 <command> -h` prints the same list for that command, read from the engine itself, so the two
cannot disagree, and a test asserts this table against them.

A modifier missing from a row is one that level does not read **yet**. It is not a claim that the
combination would be meaningless; each engine gains a modifier when its own code learns to read it, and
some have not got there. You never have to guess which: give a level a modifier it does not accept and
the error names the ones it does, so this table can always be read back off the tool itself.

### Writing your own

The modifiers above are the ones ZOMBI2 ships. If none of them says what your rate depends on,
you can write your own. It is a small class, and it needs two things.

**A `factor()` method, returning the multiplier.** The engine calls it as the run goes and hands it
the context that engine has. Every engine passes `time` and `lineages`; the rest differ:

| Engine | Context passed to `factor` |
|---|---|
| `species` | `time`, `lineages`, `diversity` |
| `genomes.family` | `time`, `lineages`, `copies`, `drivers` |
| `genomes.ordered` | `time`, `lineages`, `copies`, `chromosomes`, `drivers` |
| `genomes.nucleotide` | `time`, `lineages`, `copies`, `chromosomes`, `drivers` |
| `traits.continuous` | `time`, `lineages`, `diversity`, `drivers` |
| `traits.discrete` | `time`, `lineages`, `drivers` |
| `joint` | `time`, `lineages`, `diversity`, `drivers` |

Take `**_` and default every key you read. A key your engine does not supply never arrives, and the
genome engines pass `drivers` only on a rate that is itself driven, so it can be missing even on an
engine that lists it. One key is there and is always 0: `copies` at the nucleotide resolution, where
gene events are counted per lineage rather than per copy, so there is no copy count to pass. What you
return multiplies the rate, so 1.0 leaves it alone, 2.0 doubles it, 0.0 switches it off.

**An `implemented_for` list, naming the engines that read it.** A modifier an engine does not read
would silently return 1.0 and hand you a run that is not the model you asked for, so a level takes
only the modifiers that name it and refuses the rest. The engine names are `species`,
`genomes.family`, `genomes.ordered`, `genomes.nucleotide`, `traits.continuous`, `traits.discrete`
and `joint`.

Here is a complete one. `OnTotalDiversity` makes speciation *slow down* as lineages pile up. This
makes extinction *speed up* instead, which is the other half of the same idea and something no
shipped modifier can say:

```python
from zombi2 import species
from zombi2.rates.modifiers import Modifier

class OnCrowding(Modifier):
    """Extinction rises as lineages accumulate: at `crowd` standing lineages it has doubled."""
    implemented_for = ("species",)

    def __init__(self, crowd):
        self.crowd = float(crowd)

    def factor(self, *, diversity=0.0, **_):
        return 1.0 + diversity / self.crowd
```

That is the whole modifier. Use it the way you use a shipped one, by multiplying it onto a rate:

```python
plain   = species.simulate_species_tree(birth=1.0, death=0.2, total_time=6, seed=3)
crowded = species.simulate_species_tree(birth=1.0, death=0.2 * OnCrowding(crowd=50),
                                        total_time=6, seed=3)

print(len(plain.complete_tree.extant_leaves()))     # 290 lineages survive
print(len(crowded.complete_tree.extant_leaves()))   # 142 when extinction rises with crowding
```

Two limits. A modifier of your own is Python-only: `--death` and a `--params` file know the names
ZOMBI2 ships and cannot build a class you wrote. And if your `factor` reads `time`, making it a rate
that changes continuously rather than only when an event fires, give the class a `next_change(time)`
method returning the next moment the rate changes, so the engine stops and re-evaluates there instead
of holding it at whatever it was. That is the horizon stepping described at the end of this appendix.
`OnCrowding` needs none, because diversity only changes when something is born or dies, and the engine
already re-evaluates then.

## The Gillespie algorithm

Almost every simulation here comes from one small engine, run over and over: the same loop grows a
birth–death tree, duplicates and loses genes, and switches a discrete trait between states, given a
different list of events each time.

That engine is the **Gillespie algorithm** [@gillespie1976; @gillespie1977], an exact, event-by-event
recipe for a continuous-time process defined by rates. This section builds it from scratch, covering what a
rate is, why waiting times are exponential, how competing events race, and how those assemble into the
loop, assuming no prior exposure to continuous-time Markov chains. **Rates in, a timed history out.**

### From a rate to a waiting time

A rate says *how fast* something tends to happen: the expected number of events per unit of time. If a
gene copy is lost at rate $\mu = 0.25$, then, left alone, it is lost on average once every $1/0.25 = 4$
time units.

More precisely, a rate $\lambda$ is defined by what happens over a very short slice of time $\Delta t$.
The chance that one event fires during that slice is proportional to its length,

$$P(\text{an event in the next } \Delta t) \approx \lambda\,\Delta t,$$

and the shorter the slice, the better the approximation. Notice this says nothing about a clock ticking
down to the next event: in any instant the chance of firing is the same regardless of how long we have
already been waiting. A rate has no memory. That single fact is what makes the whole algorithm work.

It also means a rate is a statement about probability, not a fixed schedule. A loss rate of one per unit
of time does not deliver exactly one loss in every unit: over any given unit the number of events that
fire is random. That count follows the **Poisson distribution**, with mean $\lambda T$ over an interval of
length $T$ (Figure A.1).

![The count of events in a fixed window is random, not fixed. With rate $\lambda$, the number of events in one unit of time is Poisson-distributed with mean $\lambda$ (dashed line). At a low rate (left) most windows see zero or one event and a few see more; at a higher rate (right) the count spreads out around the mean. The rate fixes only the average.](figures/gillespie_poisson.pdf){width=100%}

Now fix a single event with a constant rate $\lambda$ and ask: starting now, how long until it fires?
Call that waiting time $W$. Because the chance of firing in each little slice is $\lambda\,\Delta t$ and
slices are independent, the chance of surviving without an event up to time $t$ decays to an exponential:

$$P(W > t) = e^{-\lambda t}.$$

$W$ follows an **exponential distribution** with rate $\lambda$, whose mean is $1/\lambda$. Short waits
are the most common. Because the exponential is memoryless, we never have to simulate the empty time
*between* events tick by tick: we can draw the waiting time in one shot and jump to the next event.
Drawing it is a single line: given a uniform random number $u$ on $(0, 1)$,

$$\Delta t = -\frac{\ln u}{\lambda},$$

which is what `rng.exponential(1 / lambda)` returns, its argument being the mean. Every waiting time in
ZOMBI2 is drawn this way.

### When several things can happen: the race

A real simulation never has just one possible event. A genome with many gene families can duplicate any
of them, transfer any of them, or lose any of them; a species tree with many lineages can speciate or go
extinct on any branch. At a given moment there is a whole menu of possible events, event $i$ with its own
rate $r_i$.

Treat every possible event as an independent alarm clock, each set to go off after its own exponential
waiting time. They all start together and race; the first alarm to ring is the event that happens. Two
facts govern that race, and together they *are* the Gillespie algorithm.

**When does the first event fire?** The minimum of independent exponential waiting times is itself
exponential, with a rate equal to the sum of the individual rates. So with a **total rate**

$$R = \sum_i r_i,$$

the time to the next event, whichever it turns out to be, is a single exponential draw with rate $R$.
More possible events, or faster ones, means a larger $R$ and therefore shorter waits. This is why we need
only one waiting-time draw per step, however long the menu.

**Which event fires?** The winner is event $i$ with probability equal to its share of the total rate,

$$P(\text{event } i \text{ fires}) = \frac{r_i}{R},$$

and which event wins is independent of when it happens. So the two are decided separately: draw the time
from the total rate, then pick the event on a weighted roulette wheel, each slice sized to a rate.

![The two draws that make up one Gillespie step. **(1)** Each possible event has a rate; here duplication, transfer and loss have rates 3, 2 and 1, summing to a total rate $R = 6$. **(2)** The waiting time to the *next* event is a single exponential draw with rate $R$; larger total rates give shorter waits, with mean $1/R$. **(3)** Which event fires is a second, independent draw: event $i$ wins with probability $r_i/R$: the rates laid end to end as a roulette wheel, here landing on transfer. The step then advances the clock by $\Delta t$, applies the chosen event to the state, and repeats.](figures/gillespie_step.pdf){width=100%}

In code the roulette wheel is a running sum: lay the rates end to end, draw a point uniformly along their
combined length $R$, and see which segment it lands in.

### The loop

Assembling the pieces gives the loop below. Starting from an initial state at time $t = 0$, repeat: read
off the current rates and their total $R$; draw a waiting time and advance the clock; stop if the clock
has run past the target time, or the process has died out; otherwise pick one event in proportion to its
rate, apply it to the state, record it, and go round again.

![The Gillespie loop. Each pass computes the current total rate, draws one exponential waiting time, and, unless the clock has passed the target age, fires a single event chosen in proportion to its rate, updates the state, and repeats. The output is a list of events with the exact times at which they occurred: a timed history.](figures/gillespie_loop.pdf){width=68%}

As pseudocode, the whole engine is short:

<!-- doc-test: skip - pseudocode, deliberately not the real engine -->
```python
t = 0.0
state = initial_state
history = []
while t < total_time:
    rates = event_rates(state)      # every possible event's rate, given the state
    R = sum(rates)                  # the total rate
    if R == 0:                      # nothing can happen; the process is frozen
        break
    t += rng.exponential(1 / R)     # WHEN: draw the waiting time, advance the clock
    if t >= total_time:
        break                       # the next event would fall past the horizon
    i = choose(rates, R, rng)       # WHAT: pick an event with probability r_i / R
    state = apply(state, i)         # update the state
    history.append((t, i))          # record the timed event
```

The result is not a snapshot but a **timed history**: the exact sequence of events and the exact
real-valued times at which they happened. That is what a phylogenetic simulator needs: a species tree
*is* the history of its speciation and extinction events, and a gene family *is* the history of its
duplications, transfers and losses. Because the times are drawn from continuous exponentials rather than
stepped through a fixed grid, the histories are exact: no time-step to tune, and no discretisation error.

::: note
The `rng` is a seeded `numpy` random generator. Because the waiting times and event choices are its only
source of randomness, the same seed reproduces the same history exactly. The clean core runs this loop in
plain Python.
:::

### When the rate changes with the clock

The loop above assumes the rates hold still *between* events, so that a single $\text{Exponential}(R)$
draw lands exactly on the next event. That holds whenever the rates depend only on the current state,
which changes only when an event fires. But some rates move with the clock itself, even while nothing is
happening: an `OnTime` schedule steps at fixed breakpoints, a scheduled mass extinction arrives at a set
time, and a driven rate follows a driver that changes state on its own timetable. Now $R$ is a moving
target, and a draw at today's $R$ would be wrong.

ZOMBI2 keeps the draw exact by never letting it cross a change. Every rate can report the next time it
changes on its own, and the engine takes the earliest such time, together with the next scheduled pulse
and the end of the run, as a **horizon**:

1. Compute $R$ from the rates as they stand, and the horizon.
2. Draw $\Delta t \sim \text{Exponential}(R)$.
3. If the event lands **before** the horizon, fire it: the rates really were constant over that stretch,
   so the draw is exact.
4. If it lands **after**, discard it, advance the clock to the horizon, and start again with the rates as
   they are there.

Discarding is sound for the same reason the loop works at all: the exponential is memoryless, so a
partial wait carries no information into the next stretch. The rates are piecewise constant, each piece
gets its own exact draw, and no integral or rejection step is ever needed.

### It's all Gillespie

Each level supplies its own events and rates, and the same loop realises all of them:

| Level | The events | Their rates |
|---|---|---|
| Species | speciation, extinction | `birth`, `death`, per lineage |
| Genomes | duplication, transfer, loss, origination | one rate per event |
| Traits (discrete) | switches between character states | the entries of the switch matrix |

![One engine, many events. Each level supplies its own events and rates, but all are realised by the identical loop on the right: total rate, exponential waiting time, an event chosen in proportion to its rate, apply, repeat.](figures/gillespie_everywhere.pdf){width=100%}

Swapping levels swaps the list of events and how their rates are computed; the timing machinery of total rate,
exponential wait and proportional choice never changes.

### …except when it isn't

ZOMBI2 steps outside the event-by-event loop in two places, for the same reason both times: when only
the endpoints are needed, not the whole timed history, an exact shortcut beats simulating events you
would throw away.

The first is **sequence substitution along a branch**. Once a gene tree and its branch lengths are
settled, evolving a sequence down a branch does not require the individual substitution events, only the
state at each end. The probability of ending in each state after a branch of length $t$ is given exactly
by the matrix exponential $P(t) = e^{Qt}$, so ZOMBI2 draws each site's descendant state straight from
$P(t)$ in one step (Chapter 7). Running Gillespie here would generate,
and then discard, thousands of intermediate substitutions.

The second is a **continuous trait**. Brownian motion has no events to fire, since it moves at every instant,
so there is nothing for the loop to enumerate. ZOMBI2 walks the tree branch by branch instead, in
preorder, and draws one normal per branch: mean the parent's end value, variance the exact integral of
$\sigma^2$ along the branch (Chapter 8). With a constant $\sigma^2$ that integral is just
$\sigma^2 \times$ branch length, and the values that come out have the Brownian structure exactly —
variance $\sigma^2 \times$ root-to-tip depth at a tip, covariance $\sigma^2 \times$ shared path length
between two — with no event drawn anywhere. A $\sigma^2$ that varies along the branch changes only the
integral, not the walk.

The rule of thumb is the same each time. Reach for Gillespie when you need the whole history, every
branching, gain and loss at its exact time; reach for a shortcut when the endpoints are all you need. For
the trees and genomes that are ZOMBI2's real subject the history *is* the result, so the loop is the norm
and these two shortcuts are the exceptions.
