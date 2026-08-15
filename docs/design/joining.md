# Joining — a design note

**Status: steps 1 to 9 of §14 are built; the rest is designed.** This note records what joining is
and which pairs can be joined. It also says where each joint model lives in the API, and the order
the work is done in. It is subordinate to
[`SPEC.md`](SPEC.md): where the two disagree, SPEC wins.

---

## 1. What joining is

A joint run simulates two things at once, because neither can be finished before the other starts.

That is the whole definition. The test is one question: can the driver be grown on its own and handed
over? If it can, condition. If it cannot, join.

The species tree is an **output** when the Species level is one of the things being simulated. That is
one case of joining, not what joining means. A trait and a gene family can drive each other on a tree
you hand the run, and nothing about the tree changes.

This wording replaces the older one in several places. See §13.

---

## 2. Which pairs can be joined

| Pair | Joinable | Built |
|---|---|---|
| Species – Traits | yes | yes |
| Species – Genomes | yes | yes, family resolution |
| Genomes – Traits | yes | no |
| Traits – Sequences | yes | no |
| Genomes – Sequences | yes | no, out of scope here |
| Species – Sequences | no | — |
| a level with itself | yes | yes |

**Traits – Sequences was the one correction to SPEC §3.** SPEC used to say the pair cannot be
joined, and Chapter 10 gave the reason as "a sequence never feeds back into a trait". A sequence does
feed back into a trait: `composition()` drives a trait's switch rate, and that is Co18 in the
gallery. Both directions already run as conditioning, so a cycle between them is possible. SPEC's own
generating rule allows it, because the two sit on separate branches. The table was what was wrong,
not the rule, and step 1 fixed it. Figure 10.2 gained an arrow for that pair.

---

## 3. Where a joint model lives

One rule.

**One level reading itself stays on that level's own function**, with an explicit `joint=True`.
**More than one level goes to `joint.simulate`.**

This matches conditioning, which already folds into the driven level's own run. A self-join produces
one level's result, so there is nothing to reconcile between two functions.

`joint=True` is checked both ways. Asking for it when no rate reads a live driver is an error.
A rate reading a live driver without it is an error. The flag is also where a joint run's own
arguments live, which are the runaway ceiling and `step`.

`joint.simulate` replaces `simulate_joint`. Participants are process specs, and you give what you are
not simulating.

```python
# the tree is an output, because Species is a participant
joint.simulate(species.birth_death(birth=faster_if_large, death=0.2, n_extant=100),
               traits.discrete(name="size", states=["small", "large"], switch=0.1), seed=1)

# the tree is an input, because it is given
joint.simulate(genomes.genome(loss=faster_in_caves, duplication=0.1),
               traits.discrete(name="habitat", states=["cave", "surface"], switch=slower_with_toxin),
               tree=tree, seed=1)
```

`params/retired.py` answers the old spelling with the sentence naming the new one.

---

## 4. How a driver is named

A live driver reads `"<level>:<handle>"`. Each participant declares the handle it offers.

| Driver | Reads |
|---|---|
| `"traits:<name>"` | that trait's current state or value |
| `"genomes:<family>"` | whether that family is present |
| `"genomes:count"` | the lineage's gene count |
| `"sequences:<name>"` | that gene's declared statistic |

`"trait"` becomes `"traits:<name>"`. The old spelling is singular where SPEC §7 says level names are
plural, and it cannot name one of two traits.

`TotalDiversity` also reads a live aggregate of the level being simulated, and keeps its own spelling.
Two spellings for one idea is accepted here rather than unified.

---

## 5. Per-family rates

This is a prerequisite, and it is useful on its own.

Today a genome rate applies to every family in the run. `varying_among('families', ...)` gives each
family its own rate, but the number is drawn rather than written. There is no way to name a family
and give it a rate.

```python
simulate_genomes_family(
    tree, initial_families=100, duplication=0.2, loss=0.25, seed=1,
    families=[
        family("IS1", transfer=PerCopy(1.5), loss=0.02),
        family("toxin", loss=0.4, origin=("n5", 0.4)),
        family("nuoA", module="aerobic"),
    ])
```

A family uses the run's rate for anything it does not write.

**Three arguments collapse into one.** `family_names=` becomes a list with no rates. `origins=`
becomes `origin=` on the family, which also removes the arithmetic a caller does today to work out a
placed family's id. `modules=` becomes `module=` on each member, and `resolve_modules` goes with it,
since that function exists only to cross-check two lists that are now one.

**A written family rate is the rate.** The run's `varying_among('families', ...)` draw does not
multiply it. This is SPEC §5's argument for `set_by`, applied to a family.

**Do not store a family's rate as a multiplier on the run's rate**. That encoding works for bare
numbers only. It breaks as soon as a family's rate carries its own verbs, which is where this is
heading.

**One rename.** `genomes.family(...)` currently means a whole-genome process spec. Everywhere else in
ZOMBI2 a family is a gene family. That spec becomes `genomes.genome(...)`, and `family` takes the
meaning above.

---

## 6. The three loops

### Genomes with itself

One family's presence drives the rates of another family in the same genome.

```python
simulate_genomes_family(tree, initial_families=25, seed=1, joint=True, families=[
    family("IS1", loss=0.05),
    family("cargo", transfer=PerCopy(0.02).scaled_by("genomes:IS1",
                                                     {"present": 25.0, "absent": 1.0}))])
```

The gallery has this today as two genome runs, one per family, and the code says so. The element
finishes before the second run starts. The cost is that one organism has two genomes, and the element
can never be lost as a consequence of the transfer it caused.

The engine change is small. The serial loop already recomputes every rate at the top of every pass. It
already has a per-lineage path for driven rates, and it already holds every lineage's gene content.
The driver's value is read off the live genome instead of a replayed trajectory.

**It is exact.** A live driver changes only when a genome event fires, and an event ends the current
step. So every rate is constant between events, which is what the waiting-time draw assumes.

**SPEC §3 is stronger than it needs to be here**. It says a rate driven by an aggregate of its own
level needs a new engine. The serial engine never relied on families being independent. The accurate
statement is that this costs the parallel and streaming engines, which do.

### Traits with itself

Two traits, each reading the other's current value.

```python
simulate_traits(tree, seed=1, joint=True, traits=[
    traits.discrete(name="habitat", states=["cave", "surface"], switch=faster_if_big),
    traits.continuous(name="size", rate=faster_in_caves)])
```

`simulate_traits` is new. It is needed whatever we decide, because `simulate_discrete` grows one
trait and `correlation=` only makes several parts of that one trait. A cyclic pair may also be one
discrete trait and one continuous one, so it cannot live on either existing function.

`traits.continuous(...)` is a real process spec as of step 7, where it drives speciation.

Two discrete traits reading each other are exact. The pair is a Markov chain on the product of the
two state spaces, so the existing branch walk runs it unchanged. That pair is what step 6 built. A
pair with a continuous trait in it needs the slicing of §8 and is not built.

### Sequences with itself

Two genes, each one's substitution rate reading the other's composition.

```python
joint_run = simulate_sequences(g, genomes=g, step=0.05, seed=1, joint=True, genes=[
    sequences.gene(name="chaperone", model=lg(), length=300,
                   offers=composition("KR", absent=0.08), substitution=tracks_client),
    sequences.gene(name="client", model=lg(), length=300,
                   offers=composition("KR", absent=0.08), substitution=tracks_chaperone)])
```

Both gene trees already exist, so no topology is grown. The engine walks species time in slices. See
§8 for the method and §9 for the driver it reads.

---

## 7. The cross-level joins

**Genomes with Traits.** Gene content drives a trait's switch rate while the trait drives the genome's
rates, on a tree you hand the run. Both halves exist in different shapes. The genome engine is one
race over the live set; the discrete trait engine walks branch by branch. The merge moves the trait
into the genome's race, which `_grow_joint` already does against speciation.

Transfer is available here, because on a fixed tree the set of lineages alive at an instant is known.
The tree-growing joint refuses transfer for exactly that reason.

**Traits with Sequences.** The new arrow. A trait drives a gene's substitution rate while that gene's
composition drives the trait's switch rate.

This is the first cell to be built after Genomes with Traits, and it reuses that engine's shape plus
the slicing in §8.

---

## 8. Slicing, and continuous drivers

A Gillespie needs the rate to hold still between events. A continuously changing driver breaks that.

**We slice**. Cut time into steps of `step` and hold the driver fixed inside each step. Each level
then advances across the step by whichever method it already uses.

This is one mechanism used in three places:

- a continuous trait driving speciation — **built** (step 7), and the first user of the mechanism;
- two genes' compositions driving each other;
- a trait and a gene's composition driving each other.

**Thinning is the alternative and we are not using it**. Thinning draws candidate events at a ceiling
rate and rejects each with the right probability, which is exact. It needs a ceiling, and a
continuously diffusing trait has none unless the response curve is capped with `bound=`. Capping the
curve changes the model. Slicing changes no model and approximates the simulation instead. Thinning
can be added later for anyone willing to cap their curve.

**`step` belongs to the driver, not to the run**. `Driven` already takes it and keys its trajectory by
`(driver, step)`, so the same driver read at two resolutions is two trajectories. Two participants can
need different values, because a steep response curve needs a finer step than a flat one.

**The error is first-order.** A conditioned continuous driver interpolates between two known values,
because the driver is finished and both ends are on disk. A joint run has no future to read, so it
carries the last value forward. For the same `step` the joint version is the cruder of the two.

**The rule for choosing `step`.** The driver must not move much within one slice. For a composition
driver that means `step × substitution rate` well under 1. With rate 1.0 and `step = 0.05`, a slice is
0.05 substitutions per site. On 300 sites that is 15 substitutions, and a composition change of up to
5%. Under a curve as steep as `20^(x − 0.1)` that is a 16% error on the rate. A steep mapping needs a
finer step.

**The check.** Halve `step`, rerun the same seed, and see whether the answer moves. The manual should
say this.

---

## 9. Composition for one family

This is conditioning, not joining, and it is in scope because the sequence loops need it.

`Composition` was pooled over every family in the run, and one family's was refused. The stated
reason is that a per-family composition is undefined wherever that family is absent. A driver has to
answer for every branch the target walks. `gc()` carries a `family` argument whose only job is to
refuse one, and it still does — it now names the way instead.

**The fix is a declared default.** `composition("KR", absent=0.08)` says what a branch reads where the
family is not there.

**One family per run, through `families=`.** `simulate_sequences` evolves every gene tree in the
genome run it is handed. Add `families=` and a one-family run's pooled composition is that family's
composition, which is what `gc()`'s refusal already tells the user to do. No `family=` argument is
then needed on the driver, and the sequences result does not have to learn the genome's name-to-id
map.

`families=` is also the piece that makes gene A driving gene B runnable at all. Without it the second
run re-evolves A independently, so the driver reads one history of A while the output on disk holds
another.

```python
chaperone = simulate_sequences(g, families=["chaperone"], model=lg(), length=300, seed=3)
basic = chaperone.composition("KR", absent=0.08)

client = simulate_sequences(
    g, families=["client"], model=lg(), length=300, seed=4,
    substitution=PerSite(1.0).scaled_by(basic, Curve(lambda x: 20.0 ** (x - 0.1))))
```

**A gate.** If the run was restricted with `families=`, the family is missing on some branch, and no
`absent=` was written, raise. The alternative is that the driver quietly inherits the parent's value,
which is a different model from the one the user asked for.

**Python only.** `load_driver` reads trait event logs and trait value tables. A composition is an
in-memory object. Reaching the command line would mean writing a per-node composition table and
teaching `load_driver` its header, and that is not in this project. Chapter 9's target table already
marks things Python-only, so there is a precedent for saying so.

---

## 10. The substitution log

`record=` is off by default. The sequence level records the letters at every node and no events. It is
the only level that does not record its own history.

**Size is why it stays off.** Take 300 sites and a tree whose branches total 30 time units at rate
1.0. That is 9000 substitutions for one family. A hundred families is close to a million rows, where a
genome log for the same run holds a few thousand.

**Indels get rows too.** Otherwise the log says a site changed and never says when that site arrived
or left.

**Coordinates.** A site is named by where it came from and its offset there. A position in a lineage's
own sequence shifts with every indel above it and cannot go in a log row. The alignment column index
is stable only once the run has finished, so the writer adds it rather than the engine.

The nucleotide level already uses this frame for gene events, against `(source, start, end)`, and it
already learned the rule the hard way. `nucleotide-engine.md` records it: a provenance scheme has to
record the cut at the moment it is made, because the logs cannot rebuild it afterwards.

So an insertion row carries the ids it minted **and the id they follow**. A deletion row carries the
ids dropped. With those two, a reader rebuilds any lineage's column list at any time, which is what
turns an id back into a position.

**Strand.** A nucleotide block can be reverse-complemented, so the row carries a strand column. The
sequence level always sets it to +1.

**Which method.** Build the plain forward Gillespie. It is needed anyway for anything where sites or
genes stop being independent. There is a third method. It draws the two ends with `P(t)` as now, then
fills in the path between them. That is exact and keeps the speed. Add it later if `record=` measures
too slow, and not before.

---

## 11. What must refuse

- `parallel=True` and `stream_to=` refuse a joint genome run. Those engines evolve one family per
  process, and a rate reading the whole genome breaks that. Nothing else does.
- A joint genome run needs a ceiling on how far it can grow. A count driving duplication feeds itself.
  `max_family_size` caps one family and does not cap the number of families. The guard raises, as
  `max_lineages` does, because a run cut off at a size is no longer a sample from the process asked
  for.
- `joint=True` with no live driver raises. A live driver with no `joint=True` raises.
- An unbounded response curve is fine, because we slice. Nothing refuses it.

---

## 12. Out of scope

These are real and not being built here.

- Position at the ordered resolution driving a rate.
- Three-way joins.
- Contact maps between two proteins, which need an event-driven engine and interleaved indels.
- Genomes joined with Sequences, which needs every live gene copy to carry a sequence inside the
  genome race.
- Gene content driving speciation at the ordered and nucleotide resolutions. The family resolution is
  enough.

---

## 13. What changes outside the code

**SPEC.** §2 and §4 take the wording in §1. §3's Traits – Sequences row changes to yes, with the
reason. §3's sentence about a rate driven by an aggregate of its own level is softened to name the
cost rather than demand a new engine.

**The manual.** Chapter 2's Joining section and Chapter 10's opening take the wording. Chapter 10's
two figure captions follow. Figure 10.2 gains the Traits – Sequences arrow. Chapter 9's driver table
gains per-family composition. Appendix A's table of what each level accepts gains the new rows.
Appendix B gains the substitution log and the frame it is written in.

Write all of it straight. State what joining is. Do not write against the older wording.

**The gallery.** Each part below lands with its figure, and the figure is shown before it is merged.
New cards go in `gallery/joining.py` and `gallery/crosslevel.py`.

---

## 14. Order of work

Part by part, not all at once. Steps 1 to 5 are **built** (PR #375), then 6 to 9.

1. ~~**The wording.**~~ Done. SPEC, the manual, the two docstrings, Figure 10.2's new arrow.
2. ~~**Per-family rates.**~~ Done. §5, with the `genomes.genome` rename. No figure: per-family rates
   on their own have no driver, no target and no connection, so there was nothing to draw an arrow
   between. The mobile element became step 3's card instead, where it has a loop.
3. ~~**Genomes with itself.**~~ Done. §6, `joint=True` on `simulate_genomes_family`. Figure: the
   element that makes its own genome transfer-prone, and so spreads.
4. ~~**The API shape.**~~ Done. `joint.simulate`, `species.birth_death`, the `"<level>:<handle>"`
   drivers, and the retired spellings — `simulate_joint`, `family_names=`, `origins=`, `modules=`,
   each answering with its replacement. The two built cells moved over.
5. ~~**Genomes with Traits.**~~ Done. §7, on a tree the run is handed. Transfer works there, and is
   still refused where the tree is being simulated. Figure: a habitat and a genome each other's
   driver. The **command line stops here**. These two models are Python only. One needs a rate written
   for a single named family. The other needs two levels' flags plus a tree.
6. ~~**`simulate_traits` and the trait loop.**~~ Done. §6, `joint=True` on `traits.simulate_traits`,
   with a `name=` on each trait. The pair is one Markov chain over the pairs of their states, so the
   run is exact. Figure: body size and the cave, each reading the other.
7. ~~**Slicing, and a continuous trait driving speciation.**~~ Done. §8. `traits.continuous(...)`
   is a real process spec now, and `step=` on the driven rate is the slice the diffusion is held
   fixed across. Required rather than defaulted: a step is the size of the approximation, and any
   number the code invented would be a claim about a timescale only the model knows. Figure: body
   size diffusing, and the big lineages radiating.
8. ~~**Composition for one family, and `families=`.**~~ Done. §9. `families=` on
   `simulate_sequences` and `absent=` on `composition()` / `gc()`. Conditioning, not
   joining, and here because step 9 needs it. Figure: a chaperone's GC setting a
   client's substitution rate, both families out of one genome run.
9. ~~**The sequence loop.**~~ Done. §6, on the slicing from step 7 and the per-family
   composition from step 8. `joint=True` and `genes=[sequences.gene(...)]` on
   `simulate_sequences`; the walk is by time rather than by family. One addition to the
   sketch: `gene(start=...)`, a second model to found from. Without it a gene sits at
   its own equilibrium, its composition never moves, and the loop drives nothing a
   figure could show. Figure: two genes ameliorating together.
10. **Traits with Sequences.** §7, the new arrow.
11. **`record=`.** §10, last, because nothing above depends on it.
