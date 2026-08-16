# Genomes with Sequences — a design note

**Status: built.** This note covers the last empty cell of Figure 10.2. It is
subordinate to [`SPEC.md`](SPEC.md) and follows [`joining.md`](joining.md), whose §12 put this pair
out of scope. That decision is reversed here, for the narrower version below.

---

## 1. What the pair is

The genome decides which sequences exist. It says when a gene is copied, when a copy moves to
another lineage, and when a copy dies. Each of those makes a sequence, moves one, or ends one.

The sequences decide how fast the genome changes. A genome rate reads something about the
sequence — its composition — and speeds up or slows down.

Here is the loop in one lineage. The genome duplicates `hisA`, so the lineage holds two copies. Each
copy evolves its own sequence. Together they set the lineage's composition. That composition sets the
loss rate. A high loss rate removes a copy. Removing it changes the composition again.

Neither half can go first. To evolve the sequences you need to know which copies exist and when. To
know which copies exist you need the loss rate, which reads the sequences.

**This pair can be joined and never conditioned.** A sequence is grown along the gene trees the
genome level produces, so a genome reading a finished sequence run would be conditioned on its own
output. `Composition.refuses` says so by name, and it stays saying it.

---

## 2. Scope: the family resolution only

At the family resolution a gene copy is an atom. It exists or it does not, and its sequence is one
string. A duplication copies that string, a transfer moves it, a loss ends it. The race never has to
cut a sequence.

The **nucleotide** resolution is where it gets hard. An event there takes a run of base pairs and can
fall inside a gene, so a copy's sequence stops being one string. It becomes blocks with coordinates,
and an inversion reverse-complements them. Carrying that through a live race is different work.

The **ordered** resolution sits closer to family than to nucleotide, because a gene is never split.
It is still refused here, to keep the first version to one engine and one set of tests.

**Only declared families carry a sequence.** The genome spec names them with
`families=[family("hisA")]`, and only those live copies hold a states array. Every other family races
exactly as it does today. That is what keeps the cost proportional to what the model reads.

---

## 3. The API

```python
joint.simulate(
    genomes.genome(duplication=0.1, origination=0.2, families=[family("hisA")],
                   loss=PerCopy(0.2).scaled_by("sequences:hisA",
                                               Curve(lambda gc: 4.0 ** ((0.5 - gc) / 0.3)),
                                               step=0.05)),
    sequences.gene(name="hisA", model=jc69(), length=300,
                   offers=sequences.composition("GC", absent=0.5)),
    tree=ct, seed=1)
```

Both levels are participants, so both come out of the run. The **tree** is the one thing handed over,
which is the same rule every joint model follows: give what you are not simulating. The gene trees
are not handed over either — the genome participant produces them.

`sequences.gene(...)` is the same process spec §9 of the joining note built for the sequence loop.
`offers=` says what the gene publishes; the genome rate reads it as `"sequences:<name>"`; `step=`
rides on the connection, as it does for every sliced model.

What comes back is a `JointResult` carrying `.genome` and `.sequences`.

---

## 4. What a genome event does to a sequence

Every event does one of three things: clone a sequence, end one, or draw a fresh one.

| Event | What happens to the sequence |
|---|---|
| duplication | advance the picked copy to the instant, clone it, both carry on |
| transfer | the same, and the clone lands in the recipient lineage, on its clock |
| loss | the copy's sequence stops, and is kept at that gene-tree node |
| origination | a new family, so a founding sequence drawn from the model's frequencies |
| speciation | every copy goes to both daughters, each starting from the parent's |

None of that is new. It is what a gene tree already means, and what the sequence level already
replays afterwards. The difference is only that it happens while the race runs.

**The walk has to be interruptible.** A duplication halfway through a slice must advance the parent's
sequence to that instant *before* cloning. Otherwise the two copies each redraw the stretch they
actually shared. `sequences/_loop.py` already walks that way.

---

## 5. What a sequence does to a genome rate

The statistic is read **per lineage**, pooled over the copies that lineage holds. Every driver in
ZOMBI2 is read per lineage, and keeping that means the rate grammar needs no change.

Every family-resolution target takes it: `duplication`, `transfer`, `loss`, `origination`, and
`transfer_to`. A lineage that carries none of the declared family reads the `absent` the gene
declared, which is step 8's rule and the same reason for it.

**Per-copy is version two, and is not built here.** A copy's own composition driving its own loss is
the model people actually want — duplicate retention. It breaks the per-lineage rule, so it needs a
decision about the grammar first, and it is better made against a working engine than before one.

---

## 6. The engine

A new module, `zombi2/joint/_genomes_sequences.py`, exposing `grow` and returning `Grown`, like the
other five.

It is closest to `_genomes_traits`: the genome races on a tree the run is handed. The difference is
that a trait is *in* the race, as a fifth event class, while a sequence **walks**. So the loop is:

1. freeze each lineage's composition, for the whole slice;
2. race the genome across the slice with those values;
3. at every event, advance the picked copy's sequence to that instant, then clone or end it;
4. at the boundary, advance every live copy's sequence to it, and go to 1.

Steps 1 and 4 are the slicing every other sliced engine does, and the `step` contract comes from
`_runtime/slicing.py`. Step 3 is the interruptible walk from `sequences/_loop.py`.

**Why it slices.** A composition moves with every substitution, so a genome rate reading it is never
constant. There is no interval where the rate holds still on its own, and so nothing exact to draw
against.

---

## 7. What it costs

**Not memory.** Because only declared families carry a sequence, the live array is one family's
copies across the living lineages — tens of kilobytes, not megabytes. That corrects an earlier guess
in conversation, and it changes what to measure.

**Time.** Every live copy's sequence is advanced at every slice boundary. The work is
`slices × live copies × sites`.

Measured: sixty extant species, thirty families, one of them declared, three hundred sites, `step`
0.05 over a tree of height 3.6 — so 72 slices. The joint run takes 0.13 s; the same model as two runs
in order takes 0.01 s. Ten times slower, and still a tenth of a second. The joint run is doing
something the two runs cannot, so the comparison is a ceiling on the cost rather than a like-for-like
one.

That is cheap enough to ship. It stays cheap because only the **declared** family carries a sequence;
a run that declared thirty of them would pay thirty times this. The slice boundary already uses the
matrix jump rather than the site walk, so there is no easy factor left — only a `record=True` run
walks site by site.

---

## 8. What refuses

- an ordered or nucleotide genome run, by name, with the reason in §2;
- a gene naming a family the genome spec did not declare;
- a gene something reads that `offers=` nothing;
- a connection with no `step=`, or two connections that disagree on it;
- a genome rate reading a **finished** sequence run, which is conditioning a run on its own output —
  the refusal that already exists;
- `parallel=` and `stream_to=` at the genome level, which evolve one family per process.

---

## 9. Out of scope

- per-copy reads (§5);
- the ordered and nucleotide resolutions (§2);
- a sequence driving the *species* tree, which is not a joinable pair at all;
- three-way joins.

---

## 10. What changes outside the code

**The manual.** Figure 10.2's last dashed arrow becomes solid, and its caption and the two paragraphs
about it in Chapter 10 follow. Chapter 10 gains a section for the pair. Appendix B gains nothing: both
levels write what they already write.

**SPEC.** §3's Genomes – Sequences row says the pair can be joined and never conditioned, which is
already right. §4's list of models that must be joined gains a line.

**The gallery.** One card, in `gallery/joining.py`.

---

## 11. Order of work

1. **The engine and the front door.** `_genomes_sequences.grow`, `joint.simulate` routing, the
   refusals in §8, and the tests.
2. **Measure the time** on a small case, and record the number here.
3. **The figure and the chapter.** Figure 10.2, Chapter 10's section, the gallery card, the
   CHANGELOG.
