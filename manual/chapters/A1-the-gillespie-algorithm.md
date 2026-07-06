```{=latex}
\appendix
```

# The Gillespie algorithm

Chapter 3 produced a species tree and evolved gene families along it with a
few lines of code. Behind those lines is a single small engine, run over and over. It is
the same engine that grows a birth–death tree, that duplicates and loses genes, that flips
a discrete trait between states, and that grafts on ghost lineages. Learn it once and you
understand how *every* stochastic process in ZOMBI2 is realised.

That engine is the **Gillespie algorithm** [@gillespie1976; @gillespie1977] — an exact,
event-by-event recipe for simulating a continuous-time process defined by *rates*. This
appendix builds it up from scratch: what a rate is, why waiting times are exponential, how
competing events race to fire, and how those pieces assemble into the loop that drives the
rest of the manual. No prior exposure to continuous-time Markov chains is assumed. The one
running idea is simple: **rates in, a timed history out.**

## Rates

Everything starts with a rate. A rate is a number that says *how fast* something tends to
happen — the expected number of events per unit of time. If a gene copy is lost at rate
$\mu = 0.25$, then, left alone, it is lost on average once every $1/0.25 = 4$ time units.
Double the rate and events come twice as often; halve it and they come half as often.

More precisely, a rate $\lambda$ is defined by what happens over a *very short* slice of
time $\Delta t$. The chance that one event fires during that slice is proportional to its
length,

$$P(\text{an event in the next } \Delta t) \approx \lambda\,\Delta t,$$

and the shorter the slice, the better the approximation. This is the only property of a
rate we will need. Notice it says nothing about a *clock ticking down* to the next event:
in any instant the chance of firing is the same, $\lambda\,\Delta t$, regardless of how
long we have already been waiting. A rate has no memory. That single fact is what makes the
whole algorithm work.

Rates carry units of *one over time*, so their reciprocal $1/\lambda$ is a time — the mean
waiting time we are about to meet. In ZOMBI2 the rates are the numbers you supply on the
command line or to the Python models: the speciation and extinction rates $\lambda$ and
$\mu$ of a species tree, the duplication, transfer, loss and origination rates of a gene
family, the entries of a trait's transition matrix. They may be constant, or they may
depend on the current state (how many lineages are alive, how many copies a family has);
either way, at any moment the process has a well-defined set of rates, and that is enough.

## From a rate to a waiting time

Fix a single event with a constant rate $\lambda$ and ask: starting now, how long until it
fires? Call that waiting time $W$. Because the chance of firing in each little slice is
$\lambda\,\Delta t$ and slices are independent, the chance of having survived *without*
an event up to time $t$ decays geometrically, slice after slice, to an exponential:

$$P(W > t) = e^{-\lambda t}.$$

$W$ follows an **exponential distribution** with rate $\lambda$. Its density is largest at
$t = 0$ — short waits are the most common — and its mean is

$$\mathbb{E}[W] = \frac{1}{\lambda}.$$

The exponential is the memoryless distribution, and it inherits that property directly from
the rate: having already waited a while tells you nothing about how much longer you will
wait, $P(W > s + t \mid W > s) = P(W > t)$. This is why we never have to simulate the empty
time *between* events tick by tick. We can draw the waiting time in one shot and jump
straight to the next event.

Drawing it is a single line. Given a uniform random number $u$ on $(0, 1)$, inverting the
formula above gives

$$\Delta t = -\frac{\ln u}{\lambda},$$

which is exactly what a call to `numpy`'s `rng.exponential(1/lambda)` returns (its argument
is the mean, $1/\lambda$). Every waiting time in ZOMBI2 is drawn this way.

## When several things can happen: the race

A real simulation never has just one possible event. A genome with many gene families can
duplicate any of them, transfer any of them, or lose any of them; a species tree with many
lineages can speciate or go extinct on any branch. At a given moment there is a whole menu
of possible events, event $i$ with its own rate $r_i$. What happens next?

Picture every possible event as an independent alarm clock, each set to go off after its
own exponential waiting time. They all start together and *race*; the first alarm to ring is
the event that happens. Two beautifully simple facts govern that race, and together they
*are* the Gillespie algorithm.

**When does the first event fire?** The minimum of independent exponential waiting times is
itself exponential, with a rate equal to the *sum* of the individual rates. So if we define
the **total rate**

$$R = \sum_i r_i,$$

then the time to the next event — whichever it turns out to be — is a single exponential
draw with rate $R$:

$$\Delta t \sim \text{Exponential}(R).$$

More possible events, or faster ones, means a larger $R$ and therefore shorter waits. This
is why we only ever need *one* waiting-time draw per step, no matter how many events are on
the menu.

**Which event fires?** The winner of the race is event $i$ with probability equal to its
share of the total rate,

$$P(\text{event } i \text{ fires}) = \frac{r_i}{R},$$

and — crucially — *which* event wins is independent of *when* it happens. So we can decide
the two separately: first draw the time from the total rate, then spin a weighted roulette
wheel to pick the event, each slice of the wheel sized to a rate. The figure below shows the
two draws for a small menu of three events.

![The two draws that make up one Gillespie step. **(1)** Each possible event has a rate; here duplication, transfer and loss have rates 3, 2 and 1, summing to a total rate $R = 6$. **(2)** The waiting time to the *next* event is a single exponential draw with rate $R$; larger total rates give shorter waits, with mean $1/R$. **(3)** Which event fires is a second, independent draw: event $i$ wins with probability $r_i/R$ — the rates laid end to end as a roulette wheel, here landing on transfer.](figures/gillespie_step.pdf){width=100%}

In code the roulette wheel is a cumulative sum: lay the rates end to end, draw a point
uniformly along their combined length $R$, and see which segment it lands in. ZOMBI2 does
exactly this — a `numpy.cumsum` of the rates followed by a binary search for the drawn
point — inside a small shared component (`zombi2/_sampling.py`) that every simulator calls.

## The algorithm

We now have every piece. Assembling them gives the loop shown below. Starting from an
initial state at time $t = 0$, repeat: read off the current rates and their total $R$;
draw a waiting time and advance the clock; stop if the clock has run past the target age
(or the process has died out); otherwise pick one event in proportion to its rate, apply it
to the state, record it, and go round again.

![The Gillespie loop. Each pass computes the current total rate, draws one exponential waiting time, and — unless the clock has passed the target age — fires a single event chosen in proportion to its rate, updates the state, and repeats. The output is a list of events with the exact times at which they occurred: a timed history.](figures/gillespie_loop.pdf){width=68%}

As pseudocode, the whole engine is short:

```python
t = 0.0
state = initial_state
history = []
while t < age:
    rates = event_rates(state)      # every possible event's rate, given the state
    R = sum(rates)                  # the total rate
    if R == 0:                      # nothing can happen; the process is frozen
        break
    t += rng.exponential(1 / R)     # WHEN: draw the waiting time, advance the clock
    if t >= age:
        break                       # the next event would fall past the horizon
    i = choose(rates, R, rng)       # WHAT: pick an event with probability r_i / R
    state = apply(state, i)         # update the state
    history.append((t, i))          # record the timed event
```

The result is not a snapshot but a *timed history*: the exact sequence of events and the
exact real-valued times at which they happened. That is precisely what a phylogenetic
simulator needs — a species tree *is* the history of its speciation and extinction events,
a gene family *is* the history of its duplications, transfers and losses. Because the times
are drawn from continuous exponentials rather than stepped through a fixed grid, the
histories are *exact*: no time-step to tune, and no discretisation error.

::: note
The `rng` is a seeded `numpy` random generator. Because the waiting times and event choices
are its only source of randomness, the same seed reproduces the same history exactly — the
reproducibility you saw in Chapter 3. ZOMBI2's built-in models may run this loop
in a fast Rust engine rather than in Python, but the algorithm, and the results, are
identical.
:::

## A worked example: birth and death

Make it concrete with the simplest possible menu of events — the birth–death process of
Chapter 4, which grows a species tree. Each lineage alive at the moment
speciates at rate $\lambda$ and goes extinct at rate $\mu$. Suppose the tree currently has
$n$ lineages. Every lineage contributes one possible birth and one possible death, so the
menu has $2n$ events and the total rate is

$$R = n\lambda + n\mu = n(\lambda + \mu).$$

One pass of the loop then reads:

1. **When.** Draw $\Delta t \sim \text{Exponential}\big(n(\lambda+\mu)\big)$ and advance the
   clock. With more lineages alive, $R$ is larger and the next event comes sooner — which is
   why a growing clade branches faster and faster.
2. **What.** The event is a birth with probability
   $n\lambda / \big(n(\lambda+\mu)\big) = \lambda/(\lambda+\mu)$ and a death otherwise. Since
   the lineages are interchangeable, pick the one it happens to uniformly at random.
3. **Apply.** A birth replaces the chosen lineage with two daughters (the tree branches); a
   death marks it extinct. Either way $n$ changes, so the rates are recomputed and the next
   pass uses the new total.

Run this until the clock reaches the target age and you have grown a *complete* tree, extinct
lineages and all — exactly ZOMBI2's forward birth–death mode. This is the general engine at
work. Some processes additionally admit exact shortcuts: ZOMBI2's *default* species tree is
the reconstructed tree (extinct lineages already pruned away), which a special result lets it
sample directly without the event-by-event loop. The Gillespie loop is what you reach for
when no such shortcut exists — which is almost always.

## When the rates change within a step

The algorithm above assumes the rates hold constant *between* events, so that a single
exponential draw with rate $R$ is exact. That is true whenever the rates depend only on the
current state, which covers most of ZOMBI2. But some models have rates that drift with
*time itself*, even while no event fires: episodic (skyline) birth–death rates that step at
fixed epochs, diversity-dependent rates that ease off as a clade fills up, ghost-lineage
grafting along a time-varying hazard. Then $R$ is a moving target and the plain draw would
be wrong.

The fix is **thinning** [@lewis1979thinning], also called the rejection method. Pick a
ceiling $\bar R$ that is at least as large as the true total rate over the interval. Propose
candidate events at the constant rate $\bar R$ — the easy, exact draw again — but *accept*
each candidate only with probability $R(t)/\bar R$, the true rate at that instant over the
ceiling; reject it otherwise and simply carry on. The accepted events occur with exactly the
right time-varying intensity, and no integral of the rate is ever needed. ZOMBI2 uses
thinning for its time-varying species-tree models, for state-dependent diversification
(the SSE models), and for grafting ghost lineages; the exact direct method of the previous
sections drives everything whose rates are constant between events — gene families, discrete
traits, and the per-lineage diversification models.

::: tip
Thinning is worth recognising because a poorly chosen ceiling $\bar R$ makes it slow: if the
bound sits far above the true rate, most candidates are rejected and the simulation spins.
ZOMBI2 chooses tight bounds internally, but it is the reason a model with wildly varying
rates can run slower than a constant-rate one with the same number of events.
:::

## It's all Gillespie

Step back and the unity is the point. Every level of a ZOMBI2 simulation is the *same loop*
fed a different bag of events (the figure below):

| Level | The events | Their rates |
|---|---|---|
| Species tree | speciation, extinction | $\lambda$, $\mu$ per lineage |
| Gene families | duplication, transfer, loss, origination | the four DTLO rates |
| Discrete traits | changes between character states | the entries of the $Q$ matrix |

![One engine, many events. Each level supplies its own events and rates, but all are realised by the identical loop on the right: total rate, exponential waiting time, an event chosen in proportion to its rate, apply, repeat. The lone exception is sequence substitution along an already-fixed branch, which needs only the branch's endpoints and so is computed in a single matrix step rather than event by event.](figures/gillespie_everywhere.pdf){width=100%}

Swapping levels means swapping the list of events and how their rates are computed; the
timing machinery — total rate, exponential wait, proportional choice — never changes. That
is why ZOMBI2 factors it into one shared component that every simulator reuses, and why
understanding this one appendix carries you through all the chapters that use it.

There is a single, telling exception. Once a gene tree's branches are fixed, evolving a
DNA sequence *along* a branch does not need the individual substitution events — only the
sequence at the two ends. For that, ZOMBI2 skips the event-by-event loop and jumps straight
to the answer with a single matrix operation, $P(t) = e^{Qt}$ (Chapter 15). The contrast
is instructive: you reach for Gillespie precisely when you need the *whole history* — every
branching, every gain and loss, at its exact time — and not merely a before-and-after
snapshot. For the trees and genomes that are ZOMBI2's subject, the history *is* the result.

This engine is everywhere in the manual: it grows the species trees of Chapter 4, races the
duplications, transfers and losses of the genome chapters, and flips the discrete traits of
Chapter 13 — the same loop throughout, differing only in the events on its menu.
