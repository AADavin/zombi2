# Proposal — the parameter grammar

**Status: ratified and built.** This is the wording of [SPEC](SPEC.md) §5–§7. It expresses the same
models ZOMBI2 expressed before it, plus several that were unsayable; what changed is what you write.

Some of what it proposed is **not** built. Every such place is marked on the line, and §13's worked
examples all construct and run as written, so a line carrying no mark is a line you can type. Two of
the gaps are structural rather than local: the package is still `zombi2/params/` rather than
`zombi2/params/` (§11, §12), and the forgiveness rule's refusal of a bare number where two scopes are
legal (§3) is not enforced — a bare number still takes the level's default scope.

---

## 1. Why

Three defects, all the same defect: one word doing several jobs.

- **"Modifier" names two things.** The `Modifier` docstring says a modifier returns "a dimensionless,
  non-negative multiplier". `SetBy` does not — `Rate.effective` special-cases it and takes its value
  *as* the base, units and all. `Weights` does not either — its numbers are compared with each other
  to pick a recipient, and `transfer_to` is not a `Rate` at all.

- **"per" names three things.** It is the scope (per what a rate is *counted* — `PerCopy`), the unit a
  random value is *drawn and held for* (`Drawn(per='family')`), and the unit a driver is *looked up on*.
  The first is about how many chances exist; the second and third about who differs from whom. SPEC §5
  legislates against exactly this — *"'per' is the scope word, never a modifier's"* — and then the
  modifier's own keyword is `per`.

- **The same mechanism appears twice under different names.** `OnTime({0: 1.0, 3: 0.3})` and
  `ScaledBy(Time(), {0: 1.0, 3: 0.3})` build the same object. The tour lists them in a six-row table
  that reads as six unrelated mechanisms.

And a fourth, smaller: `spread` means the σ of the lognormal a **value** is drawn from in `Drawn`, and
the σ of the per-split **step** of a geometric Brownian motion in `Inherited`. Nothing on the page says
which distribution either is.

---

## 2. The shape

```
    scope( base ) . verb( driver , mapping )
    └── what you set ──┘  └── what it depends on ──┘
```

A **parameter** is what you set. A **driver** is what it reads. A **verb** says what the reading does
to it. A **mapping** turns the driver's value into a number.

---

## 3. Parameters, and their entry points

| Parameter | What it is | Written from | Legal verbs |
|---|---|---|---|
| **rate** | how often an event fires | a **scope** | `scaled_by`, `set_by` |
| **extent** | how much an event takes once started | `Extent(...)` | `scaled_by` |
| **choice** | which candidate receives it (`transfer_to` is the only one) | `Recipients()` | `weighted_by` |
| **value** | a number that is not a rate — a trait's optimum, its pull | `Value(...)` *(not built)* | `scaled_by`, `set_by` |

A rate is always written from its scope, so *per what?* is answered on the page rather than by a
default nobody types:

```python
birth        = PerLineage(0.5)
loss         = PerCopy(0.25)
substitution = PerSite(0.01)
birth        = Global(0.5)
```

The scopes are `Global`, `PerLineage`, `PerCopy`, `PerChromosome`, `PerSite`.

**The forgiveness rule.** Write the scope wherever there is a genuine choice; a bare number is accepted
where only one scope is legal. Origination happens per lineage and nothing else, so `origination = 0.1`
is fine. Loss can be per copy *or* per lineage, and the two differ by a factor of the genome's size, so
`loss = 0.25` is refused and asks which you meant. The rule teaches the distinction by where it bites.
It presumes scope overrides exist at Genomes, and since §18 shipped they do — `loss = PerLineage(0.25)`
and `loss = PerCopy(0.25)` both run. What is still missing is the refusal: `loss = 0.25` is accepted
and takes the level's default.

The other three entry points work the same way — written only when something is chained onto them:

```python
loss_extent = Gamma(2.0, 250.0)
loss_extent = Extent(Gamma(2.0, 250.0)).scaled_by("habitat.tsv", {'aquatic': 2.0})
transfer_to = Recipients().weighted_by("competence.tsv", {'competent': 3.0, 'normal': 1.0})
reverts_to  = Value().set_by("habitat.tsv", {'aquatic': 5.0, 'terrestrial': -5.0})   # not built
```

**`Value` is not built**, so the fourth row of the table above has no entry point behind it: today
`reverts_to` takes a plain number and refuses anything else, and a trait's optimum therefore cannot be
driven at all. §7 says what stands in the way.

**Why an extent cannot be `set_by`:** an extent's base is a *distribution* over sizes, so a `set_by`
supplying one scalar would silently discard it. `extent.py` refuses this today and the proposal keeps
the refusal.

---

## 4. The verbs

Three, and only three, because the engine can only do three things with a number.

| Verb | What it does | On |
|---|---|---|
| `.scaled_by(driver, mapping)` | multiplies the base | rate, extent, value |
| `.set_by(driver, mapping)` | replaces the base, in the parameter's own units | rate, value |
| `.weighted_by(driver, mapping)` | compares the candidates of a choice | choice |

A `set_by` takes **no base in front**; the scope still stands, because replacing *how fast* says
nothing about *per what*:

```python
loss = PerCopy().set_by("habitat.tsv", {'aquatic': 1.0, 'terrestrial': 0.25})
```

Adding a fourth verb is not a naming decision — it is a change to how parameters compose (§8, §10).

### 4.1 The two shortcuts

`Random` and `Time` are the only drivers that are built in and take no naming, and they are by far the
most written. Each has one verb of its own, and it is the **only** way to write that driver:

```python
substitution = PerSite(0.01).varying_among('lineages', LogNormal(0.0, 0.3))
birth        = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})
```

`scaled_by` **refuses** a `Random` or a `Time` and names the shortcut instead — the same
mismatched-verb refusal that already exists for `ScaledBy` on `transfer_to`. So there is exactly one
spelling for each, not two.

`varying_among` also takes a named `Random` object, which is how two rates share one draw (§6).

`Time` is the one driver that can also *replace* rather than scale, and that half keeps the general
verb, because the numbers mean something different:

```python
birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})     # multipliers: 30% of what it was
birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})    # the rates themselves
```

> **There are exactly two shortcuts and there will not be a third.** A new driver — `Age`,
> `GeneContent`, `Distance` — uses `scaled_by` and gets no verb of its own. This is a rule, not a
> pattern to extend by analogy; the ten-verb version of this grammar is what it exists to avoid.

---

## 5. Drivers

What a parameter can read. Each has a unit, which is when it is evaluated or drawn, and each is either
**standing** (settled once for that unit, constant thereafter) or **changing** (moves during a unit's
life, so the Gillespie must step at its breakpoints).

| Driver | Unit | Kind | Reached by | |
|---|---|---|---|---|
| `Time()` | run | changing | `changing_at` / `set_by` | the run's clock |
| `Random(unit, law)` | lineages, families, … | per its law | `varying_among` | a value drawn for that unit (§6) |
| `TotalDiversity(cap=…)` | run | changing | `scaled_by` | lineages standing right now |
| `Age()` | lineage | changing | `scaled_by` | time since this lineage was born *(not built)* |
| `Clade({...})` | lineage | standing | `scaled_by` | which named group a lineage belongs to |
| `Distance()` | pair | standing | *(a whole `transfer_to` rule, not a driver — see below)* | patristic distance donor→recipient, in tree depths |
| a level's output | lineage | changing | `scaled_by` / `set_by` | a trait, gene content, a sequence (§9) |

`Clade` is refused at Species, and not because nobody has built it: at Species the tree is the run's
output, so while birth and death are being read there is no tree to name a clade on. Every level that
takes an already-grown tree — Genomes, Sequences, Traits — reads it.

**`Distance` is not a driver today, and neither is `Between`.** `Distance(decay=…)` is a whole
`transfer_to` rule with its own `exp(-decay × d / depth)` built in rather than a mapping to write, so
`weighted_by(Distance(), mapping)` is refused and names it (§7); and `Between({...})` is the pair
*mapping* `Clades` takes, not a wrapper that turns a per-lineage driver into a pair one (§12).

There are **two** clade classes, and they stayed two: `Clade` in `rates/clade.py` is a per-lineage
value any rate can read, and `Clades` in `rates/choice.py` is the `transfer_to` rule that weights the
**pair** (donor's clade, recipient's clade). The proposal wanted one, on the grounds that only the
verb differs; a pair rule and a per-lineage value turned out to differ in what they *return*, not only
in who reads them. They share `resolve_groups`, so a clade means the same thing whichever way it is
read, which is what the unification was for.

---

## 6. `Random` and its laws

`Random(unit, law)` is a value drawn for each unit of that kind. The **law** says what happens to it
afterwards, which is a separate question from what it starts as. Units are **plural**: `'lineages'`,
`'families'`, `'copies'`, `'sites'`, `'chromosomes'`.

| Law | | Built |
|---|---|---|
| a bare distribution | drawn once, held for that unit's life | yes (`Drawn`) |
| `Drift(dist)` | the parent's value times a mean-corrected draw at each split | yes (`Inherited`) |
| `Drift(dist, bins=n)` | the same, discretised onto a ladder — the rate-category clock | yes |
| `Reverting(dist, pull=k)` | drift pulled back toward the average instead of wandering | no |
| `Markov(dist, rate=μ)` | inherited, but re-drawn at rate μ — the random local clock | no |
| `WhiteNoise(sigma)` | drawn per branch with variance ∝ 1/duration | no |

```python
substitution = PerSite(0.01).varying_among('lineages', LogNormal(0.0, 0.3))                    # drawn, held
substitution = PerSite(0.01).varying_among('lineages', Drift(LogNormal(0.0, 0.3)))             # Brownian
substitution = PerSite(0.01).varying_among('lineages', Drift(LogNormal(0.0, 0.3), bins=8))     # rate categories

# the three laws the table marks "no" — these raise NameError today:
substitution = PerSite(0.01).varying_among('lineages', Reverting(LogNormal(0.0, 0.3), pull=0.3))
substitution = PerSite(0.01).varying_among('lineages', Markov(LogNormal(0.0, 0.5), rate=0.1))  # random local clock
substitution = PerSite(0.01).varying_among('lineages', WhiteNoise(0.3))                        # variance ∝ 1/duration
```

**A bare distribution is deliberate, not an oversight.** It follows the convention the grammar already
uses everywhere — a bare dict is a `Table`, a bare function is a `Curve` — where the plain case is
written plainly and anything else is named. Each law **owns and documents its own argument**: `Drift`'s
distribution is the per-split step, `Markov`'s is the value drawn at each reset, and a bare one is the
value itself. Nobody has to infer the role from the slot. `pull` is the word this library already uses
for reversion, on the OU trait, so `Reverting(…, pull=0.3)` is one word for one concept.

Every form is normalised to mean 1, because a multiplier whose average is not 1 changes what the base
means. A consequence worth stating plainly: **only the shape of the distribution is read**, so
`LogNormal(0.0, 0.3)` and `LogNormal(5.0, 0.3)` are the same law, and `Exponential`'s single parameter
is normalised away entirely.

A `Random` can be **named and reused**, which is the only way to say that two rates share one draw:

```python
family_speed = Random('families', LogNormal(0.0, 0.5))
duplication  = PerCopy(0.20).varying_among(family_speed)
loss         = PerCopy(0.10).varying_among(family_speed)     # exactly half, in every family
```

Two separately built `Random` objects are two draws even with identical arguments. The question is
whether you wrote one thing or two.

Not available, by decision rather than omission: **two memory structures on one axis.** A value is
either drawn afresh or carried and perturbed, never both — SPEC §5, enforced by `check_one_memory` at
every level. `Markov` would be a third law, not a composition of the first two.

---

## 7. Mappings

How a driver's value becomes a number. The mapping follows from the **driver**, never from the target.

| | |
|---|---|
| a `dict` → `Table` | a driver with named states; `default=` sets what an unlisted state gets (1.0 unless said) |
| a function → `Curve` | a driver whose value is a number |
| `Scalar(strength)` | the log-link `exp(strength · value)` — a `Curve` with the exponential chosen |
| a `dict` keyed by pairs → `Between` | the pair (donor's group, recipient's group), for a `Clades` rule |

A mapping is omitted when the driver's value is already the number wanted — a `Random` is normalised to
a dimensionless multiplier by construction.

**A mapping holds numbers; the verb and the target decide what is legal.** Today a `Table` refuses
negative numbers when it is *built*, before it knows what it will be attached to — which is why an OU
optimum cannot be driven (`Table factor for 'terrestrial' must be a finite non-negative number, got
-5.0`). That constraint belongs to the attachment, not to the table:

| Attached with | Legal numbers | Because |
|---|---|---|
| `scaled_by` | ≥ 0, dimensionless | they are multipliers |
| `weighted_by` | ≥ 0 | they are weights, normalised across candidates |
| `set_by` a rate or extent | ≥ 0, the parameter's units | a rate cannot be negative |
| `set_by` a value | anything the parameter allows | an optimum can be negative |

Nothing would become later or vaguer: attaching happens on the same line the mapping is written, so
the error would still point at the right place — and it could say *why*, which "a table cannot be
negative" never could, because it was not true.

**Not built.** The check has not moved, so the four rows above are the design and not the code, and an
OU optimum still cannot be driven. Moving it is cheap on its own; what it is worth waiting for is
`Value`, since a driven optimum is the only parameter that wants a negative number and there is
nothing else for the loosened rule to admit.

Separating the driver from the mapping is also what would make a kernel writable. `Distance(decay=1.0)`
bundles both, so `exp(-decay · d / depth)` is the only available shape:

```python
transfer_to = Distance(decay=1.0)          # the whole rule, its shape included — what ships

# not built: Distance is refused as a driver, and the refusal names the line above
transfer_to = Recipients().weighted_by(Distance(), lambda d: math.exp(-d))       # the same fall, written out
transfer_to = Recipients().weighted_by(Distance(), lambda d: 1 / (1 + d ** 2))   # a power law
transfer_to = Recipients().weighted_by(Distance(), {0.0: 1.0, 0.5: 0.2, 1.0: 0.0})
```

What `Recipients().weighted_by(...)` does take is a driver another level grew, which is the rest of
the grammar and is built:

```python
transfer_to = Recipients().weighted_by("competence.tsv", {'competent': 3.0, 'normal': 1.0})
```

---

## 8. Composition

- **Verbs chain, and their factors multiply.** `PerCopy(0.25).scaled_by(a, m).varying_among('families', d)`
  is one rate reading two drivers. Which chains a given level *reads* is that level's declaration, and
  this one the genome engines refuse: a driver weights lineages and a per-family draw weights copies by
  their family, so combining them means weighting by the product, which no resolution implements yet.
- **One base.** A parameter carries at most one `set_by`, written first, because everything to its left
  is a base it would discard.
- **One memory structure per axis** (§6).
- **Weights on a choice multiply and are then normalised across the candidates**, so chaining two is
  meaningful: prefer close relatives *and* run a highway between two distant clades. `Choice` chains
  them, but no genome resolution reads more than one, so a second is refused at the engine rather than
  dropped — *"transfer_to carries 2 weightings, and the genome engines read one."*

```
effective rate  =  scope(base)  ×  every scaled_by / varying_among / changing_at
extent          =  base × the same                       (no scope: already absolute)
choice          =  the product of the weights, normalised
```

with one line beside it: a `set_by` replaces `scope(base)` rather than multiplying it.

**A modifier is a connection whose effect is scaling. That is the whole of it.** With the word narrowed
that far the central equation is true as written, and `set_by` and `weighted_by` sit outside it as the
other two verbs rather than as exceptions to it.

---

## 9. Conditioned and joint

The same driver, two supplies, and that is the *only* difference:

```python
loss  = PerCopy(0.25).scaled_by("habitat.tsv", {'aquatic': 4.0, 'terrestrial': 1.0})   # conditioned
birth = PerLineage(1.0).scaled_by("trait",     {'aquatic': 4.0, 'terrestrial': 1.0})   # joint
```

**Conditioned** — the driver was grown first and written to a file, so two ordinary runs in order do
the whole job and nothing comes back. **Joint** — driver and target advance in one engine and each
feels the other, so the driver is named as a live quantity.

The two lines are different parameters because the pairings are: at Genomes a driver is a file, and
the live names — `"trait"`, `"genomes:count"`, `"genomes:<family>"` — are read by `birth` and `death`
in a joint run and nowhere else. Which pairs are legal at all is SPEC §3.

A model never references another model object. A driver names a **quantity**, and the joint call binds
it to whatever is running:

```python
habitat = DiscreteTrait(states=("aquatic", "terrestrial"), switch=0.4, start="aquatic")
run = simulate_joint(birth=birth, death=PerLineage(0.1), trait=habitat, n_extant=200, seed=1)
```

so the object graph has no cycle even when the model does, each model validates alone, and construction
order stops mattering. A named quantity a pairing cannot supply fails at the call and says what that
pairing does supply: naming `"genomes:count"` beside `trait=` answers *"with trait=, drive from the
live trait — scaled_by("trait", ...); got driver(s) ['genomes:count']"*.

---

## 10. What cannot be expressed

Stated so nobody has to discover it.

- **Adding.** Composition is multiplicative and a verb never sees the base, so `+2` cannot be written —
  it would be `× (1 + 2/base)` and nothing knows `base`. The only route is `set_by` with the base folded
  into your own function, at which point there is no base. A real `.added_to` is a fourth verb and a
  change to §8, not a mapping.
- **A smooth function of time.** `Time()` takes a schedule, because a rate that changes in steps can be
  stepped to exactly; a continuously varying rate needs the engine to integrate rather than sample at a
  point. This is also what blocks a smooth Weibull hazard.
- **A level reading an aggregate of itself** — every family's rate reading the genome's total size, a
  pathway's completeness driving origination. SPEC §3 classifies this as a joint model needing an
  engine, not a knob.
- **A user function in the written form.** `Curve(<lambda>)` does not reparse, so a run's record of its
  own parameters is faithful only when every value is writable. True today, unchanged here, and it
  limits what "the run records what it used" can promise.

---

## 11. Modules

**Not built.** The grammar shipped inside the existing `zombi2/params/`, whose files kept their names:
`Random`, `Drift` and `TotalDiversity` live in `modifiers.py` beside `Drawn` and `Inherited`, and the
verbs in `verbs.py` with the methods on `rate.py` calling them. The tree below is the reorganisation
this proposed, left as written.

```
zombi2/params/
    parameter.py      Rate, Extent, Choice, Value — the chainable objects and the verbs
    scope.py          Global, PerLineage, PerCopy, PerChromosome, PerSite
    driver.py         Time, Random, TotalDiversity, Age, Clade, Distance, Between, Trait, GeneContent
    law.py            Drift, Reverting, Markov, WhiteNoise — what a Random value does over the tree
    mapping.py        Table, Curve, Scalar
    distributions.py  Fixed, Exponential, Gamma, LogNormal, Uniform, Geometric
    conditioned.py    reading a driver's trajectory back off a file
    parse.py          the written form, both directions
```

Eight files, one noun each, from twelve. `modifiers.py` splits three ways — the verbs become methods on
`parameter.py`, `OnTime`/`OnTotalDiversity` become `Time`/`TotalDiversity` in `driver.py`, and
`Drawn`/`Inherited` become `Random` plus two entries in `law.py`. `verbs.py`, `values.py` and
`clade.py` vanish into their nouns. Today's `driver.py` is a file reader wearing the concept's name; it
becomes `conditioned.py`.

The verbs living as methods in the core module is deliberate: a driver, a law, a mapping or a
distribution is added by dropping in a class, and a **verb** cannot be added without editing how
parameters compose.

If `parameter.py` grows past about 700 lines it splits into `parameter.py` (what you write) and
`evaluate.py` (what the engine calls) — not into another modifier grab-bag.

---

## 12. Migration

| Today | Write instead |
|---|---|
| `0.25 * OnTime({0: 1.0, 3: 0.3})` | `PerCopy(0.25).changing_at({0: 1.0, 3: 0.3})` |
| `1.0 * OnTotalDiversity(cap=100)` | `PerLineage(1.0).scaled_by(TotalDiversity(cap=100))` |
| `0.25 * Drawn(per='family', spread=σ)` | `PerCopy(0.25).varying_among('families', LogNormal(0.0, σ))` |
| `1.0 * Inherited(per='lineage', spread=σ)` | `PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, σ)))` |
| `0.25 * ScaledBy(driver, m)` | `PerCopy(0.25).scaled_by(driver, m)` |
| `SetBy(driver, m)` | `PerCopy().set_by(driver, m)` |
| `Weights(driver, m)` | `Recipients().weighted_by(driver, m)` |
| `spread=` | a written distribution |
| `per='family'` | `'families'` — the units go plural |
| `zombi2.params` | `zombi2.params` — the rename to `zombi2.params` (§11) is not built |

**`Distance` and `Clades` are not on that list.** The proposal wanted both dissolved into
`weighted_by` — `Distance(decay=k)` into a driver plus a kernel, `Clades({...}, Between({...}))` into
`weighted_by(Between(Clade({...})), {...})` — and neither is built (§5, §7). `Distance(decay=k)` and
`Clades(groups, Between(pairs))` are the shipped spellings, and for a clade-to-clade transfer weight
`Clades` is the only spelling there is. Removing them would take a model away.

`parse.py` keeps every retired spelling in its refusal table, as it does today for `DrivenBy`,
`ByFamily`, `ByLineage` and `FromParent`: the error names the replacement rather than guessing one. The
singular-to-plural switch is a small safety win — `'family'` and `'families'` are different strings, so
old code fails loudly instead of quietly meaning something new.

---

## 13. Worked examples

Every line here constructs, and every line is accepted by the engine that reads it — the parameter
names are the keyword arguments of `simulate_species_tree`, `simulate_genomes_*`,
`simulate_sequences`, `simulate_continuous` / `simulate_discrete` and `simulate_joint`.

```python
# ── Species ────────────────────────────────────────────────────────────────
birth = PerLineage(0.5)
death = PerLineage(0.1)
birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})                       # skyline
birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})                      # the rates themselves
birth = PerLineage(1.0).scaled_by(TotalDiversity(cap=100))                  # diversity-dependent
birth = PerLineage(0.5).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))   # ClaDS
birth = Global(0.5)

# ── Genomes ────────────────────────────────────────────────────────────────
duplication = PerCopy(0.10)
loss        = PerCopy(0.20).varying_among('families', LogNormal(0.0, 0.5))
loss        = PerLineage(0.25)                                              # the deletion budget
loss        = PerCopy(0.25).scaled_by(Clade({'burrowers': ['n12']}), {'burrowers': 2.0})
loss        = PerCopy().set_by("habitat.tsv", {'aquatic': 1.0, 'terrestrial': 0.25})
loss_extent = Extent(Gamma(2.0, 250.0)).scaled_by("habitat.tsv", {'aquatic': 2.0})

# duplication always exactly twice loss, family by family
family_speed = Random('families', LogNormal(0.0, 0.5))
duplication  = PerCopy(0.20).varying_among(family_speed)
loss         = PerCopy(0.10).varying_among(family_speed)

# who receives a transfer: close relatives; a highway between two clades; or a driver
transfer_to = Distance(decay=1.0)
transfer_to = Clades({'A': ['n12', 'n27'], 'B': 'n40'},
                     Between({('A', 'B'): 100.0, ('B', 'A'): 100.0}, default=1.0))
transfer_to = Recipients().weighted_by("competence.tsv", {'competent': 3.0, 'normal': 1.0})

# ── Sequences ──────────────────────────────────────────────────────────────
substitution = PerSite(0.01)
substitution = PerSite(0.01).varying_among('lineages', LogNormal(0.0, 0.3))          # relaxed clock
substitution = PerSite(0.01).varying_among('lineages', Drift(LogNormal(0.0, 0.3)))   # autocorrelated
substitution = PerSite(0.01).scaled_by(Clade({'rodents': ['n12']}), {'rodents': 3.0})  # a local clock
substitution = PerSite(0.01).scaled_by("habitat.tsv", {'aquatic': 2.5}) \
                            .varying_among('lineages', LogNormal(0.0, 0.3))

# ── Traits ─────────────────────────────────────────────────────────────────
rate       = PerLineage(1.0).scaled_by("habitat.tsv", {'aquatic': 3.0, 'terrestrial': 1.0})  # a diffusion
switch     = PerLineage(0.4).scaled_by("habitat.tsv", {'aquatic': 2.0, 'terrestrial': 1.0})  # a discrete trait
reverts_to = 5.0     # an OU optimum is a plain number — there is no `Value` to drive it (§3, §7)

# ── Joint ──────────────────────────────────────────────────────────────────
birth = PerLineage(0.4).scaled_by("trait", {'aquatic': 2.0, 'terrestrial': 1.0})
death = PerLineage(0.1).scaled_by("trait", {'aquatic': 0.5, 'terrestrial': 1.0})
birth = PerLineage(0.3).scaled_by("genomes:count", lambda n: 1 + n / 1000)
```

Four forms are deliberately absent. Each of the first three constructs and is then refused by the
engine that would read it, which is why constructing is not the test this section is held to:

- `birth = PerLineage(0.5).scaled_by(Clade({...}), {...})` — the species engine reads no driver at
  all, and a clade in particular cannot be read while the tree that defines it is the thing being
  grown (§5). The ClaDS line above is the species rate that does vary between lineages.
- `loss = PerCopy(0.25).scaled_by(Clade(...), {...}).varying_among('families', ...)` — the genome
  engines take a per-family draw or a driver, not both: *"combining them means weighting by the
  product. Use one or the other for now."* The two halves are shown on separate lines instead.
- `transfer_to = Recipients().weighted_by(...)` chained twice — no genome resolution reads more than
  one weighting (§8). The three transfer lines above are three separate rules, not a chain.
- `size = BM(...)` / `size = OU(...)` — there are no `BM` and `OU` classes to construct in the first
  place. A continuous trait is one `simulate_continuous` call, Brownian by default and
  Ornstein–Uhlenbeck once `reverts_to` and `pull` are given.

---

## 14. Settled, and still open

**Settled, and built.** Named laws with a bare distribution for the plain case (§6). Exactly two
shortcuts, `varying_among` and `changing_at`, each the only spelling for its driver, with `scaled_by`
refusing those two and no third shortcut ever (§4.1). `among`, with plural units (§6). Scope overrides
at Genomes (§18) — `loss` takes `PerCopy` or `PerLineage`, which is what the forgiveness rule in §3
needed to have anything to teach.

**Settled, and not built.** The legality check moving from the mapping to the attachment (§7). It is
worth doing with `Value`, not before: a driven trait optimum is the only parameter that wants a
negative number, so until `Value` exists the loosened rule would admit nothing.

**Open.**

- **`Value`, the fourth parameter kind.** Not built at all — a trait's optimum takes a plain number
  and cannot be driven (§3, §7). The name was the open question when this was written; the entry
  point itself is the open question now.
- **What actually ships.** `Age`, `Reverting`, `Markov` and `WhiteNoise` are written here because the
  shape of the law and driver menus should be settled with them in view. Nothing here commits to
  building them, and none of the four is built.
- **`Distance` and `Between` as drivers.** §5 lists them in the driver menu and §7 uses a `Distance`
  kernel to argue that separating the driver from the mapping is what makes a kernel writable. Neither
  is built: `Distance(decay=…)` is a whole `transfer_to` rule and `Between({...})` is the pair mapping
  `Clades` takes. Dissolving them into `weighted_by` would have to keep both models writable.

## 15. Cost

Every rate in the manual, the docs site, the gallery, `analyses/`, the tests and the CLI help is
written in the old vocabulary. This is wide, shallow, and should be done in one change rather than
spread over several, so that no chapter is ever half-migrated. It is cheapest now — there are almost no
users to break — and it gets dearer the moment the paper lands.

---
---

# Part II — implementation decisions

Sections 1–15 settle **what you write**. These four settle **what the objects are**, which is what a
migration actually needs and what a fleet of agents would otherwise each answer differently.

> **Ratified by Adrián, 2026-08-09.** §16 and §17 as written: `*` is removed rather than deprecated,
> and `Random` is a new spelling of the object the engines already carry. §18 and §19 are **built and
> merged** — the scope table ships and the generative round-trip test guards rates, extents and
> choices. §20 step 3, the migration, has since landed; what remains is step 4, the prose pass.
>
> The two answers matter together: verbs-only means every call site must change, and same-object-underneath
> means none of the five engines does. That is what makes this wide and shallow rather than deep.

---

## 16. What `PerCopy(0.25)` returns, and how it composes

**Decision: the verbs live on `Rate`, and a scope constructor returns a `Rate` directly.**

`PerCopy(0.25)` evaluates to `Rate(base=0.25, scope=PerCopy, modifiers=())`. There is no intermediate
`Scope` object to be coerced later, so `as_rate` keeps exactly two cases: a bare number, which gets the
level's default scope, and an already-built `Rate`.

**A scope stops carrying the base.** Today it carries one, and `with_default_scope` does
`default(self.base)` — so the same number lives in `Rate.base` and in `Rate.scope.base`, and nothing
enforces that they agree. A scope becomes a **unit marker**: `PerCopy.unit == "copies"`, no fields.
`Scope.total_of(base, **counts)` stays as it is and is the only thing the engine calls.

**Every verb returns a new frozen `Rate`.** `a = PerCopy(0.25)` then `b = a.scaled_by(...)` leaves `a`
untouched. That is worth stating because it is the one thing the mutable-model sketch could not
promise: a parameter sweep built by chaining cannot alias.

**`Rate.__repr__` becomes the written form.** Today it is the dataclass repr, so a log or an error
prints `Rate(base=0.25, scope=PerCopy(base=0.25), modifiers=(Drawn(...),))` rather than what was typed.
One renderer instead of two also means `written_form` cannot drift from `__repr__` (§19).

**`__mul__` is removed.** One way to compose, and the `*`-specific refusals in `Rate.__mul__` and
`Scope.__mul__` become verb-level checks instead — `set_by` after anything else is refused because
`Rate` already carries a base, which is the same rule stated once rather than at two operand positions.

*Cost:* `RateCompositionError` and its two carefully written messages are rewritten as verb errors, and
`as_rate`'s `Scope` and `Modifier` branches go. About 200 lines move; none of it is subtle.

---

## 17. Is `Random` a driver or sugar?

**Decision: a driver in the written form, the same carried modifier internally.**

This is the load-bearing decision of the whole migration, and it is what makes it a re-spelling rather
than a rewrite.

`.varying_among('families', law)` builds an object whose `reads` is `(DRAWN, 'families')` for a bare
distribution and `(INHERITED, 'families')` for `Drift` — exactly what `Drawn` and `Inherited` produce
today. Therefore **nothing in any engine moves**: `Rate.carried_modifiers`, `values_at_birth`,
`values_at_split`, `check_one_memory` and every level's `IMPLEMENTED_MODIFIERS` keep working against
`reads` without knowing the class names changed.

Shared draws are unchanged for the same reason. `values_at_birth` caches on `id(m)`, and
`family_speed = Random('families', LogNormal(0.0, 0.5))` **is** that object, so

```python
duplication = PerCopy(0.20).varying_among(family_speed)
loss        = PerCopy(0.10).varying_among(family_speed)
```

shares one draw by the mechanism that already exists, with no new concept.

**Units go plural, everywhere at once.** `UNITS` becomes `('run', 'lineages', 'chromosomes',
'families', 'copies', 'sites')` and `reads[1]` with it, so there is one spelling rather than a plural
one on the page and a singular one underneath. That touches about eight `carried_modifiers(unit=…)`
call sites; mechanical.

**Correction to §6.** Adding a law is *not* always "drop in a class". `Drift` and `Reverting` are, and
`WhiteNoise` is at Sequences, because all three are settled when the unit is born or when its branch is
known. `Markov` is not: it re-draws **mid-branch**, so it needs a third kind alongside `DRAWN` and
`INHERITED`, an entry in `CARRIED_KINDS`, and a `next_change` for carried values so the Gillespie stops
at each reset. That machinery exists for conditioned drivers and would have to be pointed at carried
ones. `Markov` is an engine change wearing a law's clothes, and §6 should say so.

---

## 18. The scope legality table

The forgiveness rule needs this and it does not exist. Read from today's `as_rate(default_scope=…)`
calls, with the second column being what the proposal admits.

| Level | Parameter | Legal scopes | Default | Bare number? |
|---|---|---|---|---|
| Species | `birth`, `death` | `PerLineage`, `Global` | `PerLineage` | **no** |
| Genomes (family, ordered) | `duplication`, `transfer`, `loss` | `PerCopy`, `PerLineage` | `PerCopy` | **no** |
| Genomes (ordered) | `inversion`, `transposition`, `translocation` | `PerCopy`, `PerLineage` | `PerCopy` | **no** |
| Genomes (all) | `origination` | `PerLineage` | — | yes |
| Genomes (ordered) | `fission`, `fusion`, `chromosome_loss` | `PerChromosome` | — | yes |
| Genomes (ordered) | `chromosome_origination` | `PerLineage` | — | yes |
| Genomes (nucleotide) | every gene event | `PerLineage` | — | yes |
| Sequences | `substitution` | `PerSite` | — | yes |
| Traits | `rate`, `switch` | `PerLineage` | — | yes |

**Built, as of the change that added it.** The nucleotide resolution already counted its gene events
per lineage by design — the rate says how often a lineage acts and the extent says how much DNA it
touches — so it has no second scope to offer rather than a missing one. The chromosome tier keeps
`PerChromosome` alone: a per-lineage fission is coherent but nobody has asked for it, and the rule
here is to admit a scope because someone wants the model, not because the class exists.

**The rule:** a bare number is accepted exactly where the legal set has one member. Everywhere else the
scope must be written. Nothing arbitrary is being chosen — the rule falls out of the table, and the
table falls out of which scopes are meaningful.

**`Global` stays at Species only.** One shared budget for the whole system is a statement about the
tree; at Genomes it would mean one loss budget shared across every lineage at once, which nobody has
asked for and which would make a genome's rate depend on how many other lineages exist. Excluded until
someone wants it, rather than admitted because the class exists.

*Cost:* every non-default cell is engine work. Nine parameters gain a second legal scope, and each one
needs the two-step Gillespie described in task #7 — the total computed from the scope, then the victim
chosen within it — plus a test that the `PerCopy` path is numerically unchanged. This is the largest
single piece of work in the proposal and it is a **precondition**, not a follow-up: without it the
forgiveness rule refuses `loss = 0.25` while offering nothing to write instead.

---

## 19. The written-form round trip

**What exists.** `test_written_form_round_trips` renders, reparses, and compares *objects* — which is
the right comparison — but over six hand-written strings. Under a migration where every `__repr__`
changes, enumeration is not coverage.

**Object equality is not semantics here.** `Rate` is a frozen dataclass, so `==` compares fields, and
some fields are deliberately outside `__eq__` — a `Driven`'s verb, for one. So a rendering can compare
equal while behaving differently.

**Decision: a generative round-trip test, checked numerically.** Enumerate the legal space
programmatically — scope × verb × driver × mapping × law, drawn from the §18 table and the §5–§7 menus
— and for each expression: render, reparse, compare the objects, **and** compare `effective()` over a
fixed grid of contexts (a few times, diversities, counts, driver states). Numeric agreement is the
check that matters; object equality is the cheap pre-filter.

**The parser is the highest-risk file.** Expressions move from `BinOp` (`0.25 * Drawn(...)`) to a call
on an attribute (`PerCopy(0.25).varying_among(...)`), so the ast whitelist changes shape: attribute
names must be checked against the verb set, not just call names against a class list. Nothing else in
the migration has that property.

**Keep `*` parseable for one release, as a reader only.** An old `--params` file says
`0.25 * Drawn(per='family', dist=LogNormal(0.0, 0.5))`. If `BinOp` stops parsing entirely, that fails with a syntax
error instead of the retired-name message that names the replacement. The parser should still
*recognise* a product far enough to reach the name and refuse it properly.

*Cost:* the generative test is perhaps 150 lines and pays for itself the first time it catches a
renderer that parses back to something valid but different — which is the failure mode no hand-written
list finds.

---

## 20. Sequencing

1. ~~**Scope overrides at Genomes** (§18), alone, with the numeric-equivalence check.~~ **Done** (#326).
   14,236 events across 15 configurations are byte-identical on the per-copy path.
2. ~~**The generative round-trip test** (§19) against *today's* grammar, so it is known-good before it
   is the thing guarding the migration.~~ **Done** (#326, extended in #327). It has since caught two
   real defects — a `SetBy` that dropped its `step`, and a choice written in a form the CLI refuses.
3. ~~**The migration.** The core is **not** a fan-out job: one pass does §16 and §17 serially — `Rate`'s
   verbs, the scope constructors, `Random` — because every other worker needs the exact object shapes
   to exist first. Only then do the call sites fan out by area: each engine, the CLI, the manual, the
   docs, the gallery, `analyses/`, the tests. Code and prose in the same change, same sentences with
   new spellings, so the diff is reviewable.~~ **Done**, on `feat/parameter-grammar`; the engines are
   byte-identical to `main` across 39 event streams. The adversarial review after it is what found the
   gaps this document now marks.
4. **The prose pass** — simplifying the chapters — separately, afterwards. Rewriting the words and
   renaming everything in one diff produces something nobody can check, including me.
