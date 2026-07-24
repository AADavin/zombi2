# Weighting transfers between clades — a design study

*Branch `weighting_transfers`. A working proposal, not a permanent design doc: on landing it
folds into `SPEC §5` (the choice slot) and manual Ch4/Ch9, and this file is deleted. Adrián has the
last word on every decision flagged **OPEN** below.*

## What a user wants and cannot get today

> "Simulate high rates of transfer **between two clades** I name — and not with the rest of the tree."

Trying this today (branch experiment, see the verdict that opened this task):

- **No `clade` primitive.** Transfer targeting reads a per-lineage *trait*, so the user must smuggle
  clade membership in as a fake trait: hand-build a `trait_events.tsv` with one row per node (97 rows
  for a 30-tip tree), in an exact format, joined by node id. There is no helper that paints a subtree.
- **The recipient weight cannot see the donor.** `transfer_to` weights each candidate recipient by a
  function of *that candidate only* (`recipient_index` in `zombi2/genomes/_transfer.py`). So you can
  confine transfer to `A ∪ B` (0 % leaks to the rest — that part works), and you can express a strict
  *directional* `A → B`, but you **cannot** express the symmetric "between A and B, not within":
  letting both clades donate and receive necessarily re-admits `A → A` and `B → B`. In the experiment,
  41 % of transfers landed within a clade, with no knob to forbid it.

Two gaps, then: **(1) no way to name a clade as a targeting group**, and **(2) the choice is
donor-blind.** The design closes both.

## The principle that decides the shape

Ch9 draws a line the rest of the codebase holds to:

> *"Ask whether the rate reads a value that some **other level** produced; if it only reads the tree,
> it is not [a coupling]."*

Clade membership is a pure function of the species tree. So **clade-targeted transfer is the genome
level reading the tree it already lives on — exactly like `"distance"`.** `"distance"` weights
recipients by patristic distance, computed from the tree, with no driver file and no coupling. A clade
rule is its sibling: another *topological* `transfer_to` rule, **not** a `DrivenBy`.

This is the crux. It keeps `DrivenBy` meaning what it means ("reads another level"), and it gives the
user a self-contained rule with no fake trait file. It also matches the model-admission discipline in
memory: this is a **within-level refinement of the existing choice slot**, not a new engine — the
transfer already couples two contemporaneous lineages and `recipient_index` already receives the
donor, so donor-conditioned weighting adds no new cross-lineage dependency the engine did not have.

## The design

Two additions, one shared kernel.

### 1. A shared kernel: `Between` — a weight over ordered group pairs

The 2-D analogue of `Table`. `Table` maps *one* state to a factor; `Between` maps an *(from-group,
to-group)* pair to a weight.

```python
Between({("A", "B"): 1.0, ("B", "A"): 1.0})          # A↔B only; unlisted pairs → default
Between({("A", "B"): 3.0}, default=1.0)              # A→B three times baseline, else baseline
```

- Keys are ordered pairs `(donor_group, recipient_group)`; matched by string form, like `Table`.
- `default` (**OPEN**: I propose `0.0`) is the weight for any pair not named. `0.0` makes the
  two-clade idiom a single line — name the flows you want, everything else cannot receive — and reuses
  the established "weight 0 = cannot receive" rule. A run that wants a baseline sets `default=1.0` and
  up-weights the pairs it cares about.
- Lives in `zombi2/rates/mapping.py` beside `Table`/`Curve`/`Scalar`. It is **not** a `Mapping`
  (a `Mapping.multiplier` takes one value); it is a sibling *kernel* whose `.weight(g_from, g_to)`
  takes two. The choice-slot caller detects it and passes both group labels.

### 2. `Clades` — a topological `transfer_to` rule (the user's answer)

A new recipient rule beside `Distance`, in `zombi2/genomes/_transfer.py`:

```python
from zombi2.genomes import Clades

g = simulate_genomes_unordered(
    tree, transfer=2.0, initial_families=40, seed=11,
    transfer_to = Clades(
        {"A": ["n12", "n27"], "B": ["n40", "n55"]},      # each group = MRCA of these tips
        Between({("A", "B"): 1.0, ("B", "A"): 1.0}),     # between-only, default 0
    ),
)
```

- **Groups** are named by the **MRCA of a set of tips** (robust, and how a user thinks) or directly by
  an internal node id. The engine resolves each to its subtree in the **complete** tree (internal and
  extinct nodes included — transfers happen among all contemporaneous lineages). Groups must be
  **disjoint**; a lineage in none of them is in an implicit group `"rest"`, usable as a kernel key.
- Membership is precomputed **once** per run (like `depth = mean_root_to_tip(tree)` already is), as a
  plain `{node_id: group_label}` dict — no `DriverTrajectory`, no file, because membership is constant
  along every branch.
- `recipient_index` gains one branch: read `g_d = group[donor]`, weight candidate `k` by
  `kernel.weight(g_d, group[alive[k]])`, normalise over candidates, pick. This is the same
  normalise-and-pick that `"uniform"`, `Distance` and `DrivenBy` already use, so it inherits the
  weight-0 → no-op thinning verbatim (`recipient_index` returning `None` when every candidate weighs 0,
  already handled in `_do_transfer`).

That single expression produces exactly what the user could not get before: `A ↔ B` only, `0 %` within
a clade, `0 %` to the rest.

### 3. `DrivenBy(trait, Between(...))` — the donor-conditioned *coupling* (natural follow-on)

When the groups are a genuine **trait** (habitat, competence — another level), donor-conditioned
weighting *is* a coupling, and it reuses `DrivenBy` with the same kernel:

```python
transfer_to = mod.DrivenBy(habitat, Between({("marine", "marine"): 3.0,
                                             ("terrestrial", "terrestrial"): 3.0}, default=1.0))
```

Same code path: `recipient_index` already holds the driver trajectory `to_traj`; it reads
`to_traj.value(donor, t)` for the donor group and `to_traj.value(alive[k], t)` per candidate, and
applies the kernel. Today's 1-D `DrivenBy(trait, Table)` is the special case where the kernel ignores
the donor. This is a small addition once `Between` exists, and it is the honest home for the "assortative
transfer within an ecological guild" model.

### Why this is coherent

| reads… | topological rule (no file) | coupling (`DrivenBy`) |
|---|---|---|
| the **tree** | `"distance"`, **`Clades`** | — |
| **another level** (a trait) | — | `DrivenBy(trait, Table)` → **`DrivenBy(trait, Between)`** |

`Clades` sits under "reads the tree", where `"distance"` already lives. `Between` is the one new kernel
both columns share. Nothing new in the rate grammar; nothing that changes *how fast* or *how many* —
`Between` only redistributes *who receives*, so `SPEC §5`'s choice-slot invariant holds unchanged.

## Semantics and edge cases

- **The choice slot still only redistributes.** Weights normalise over the live candidates; the
  transfer *rate* is untouched. `SPEC §5` holds.
- **"High rate" is composed, not baked in.** `Clades` decides *who receives*; the *amount* is the
  `transfer` rate. With `default=0`, a donor copy drawn from `rest` finds no eligible recipient and the
  event is a no-op (Poisson thinning — correct, the kept process is exactly the intended one). So the
  realised `A↔B` rate is `transfer × (copies in A ∪ B)`, lower than `transfer × all copies`. To make it
  genuinely high, **raise the base `transfer`**. Deliberately *not* in scope: a topological *rate*
  covariate ("only A and B donate") — that would be a rate reading the tree, i.e. a new modifier family
  (`OnClade`?), a separate and larger question. The kernel alone is correct; base-rate is the volume knob.
- **`self_transfer`.** If on, the donor is among its own candidates; `kernel.weight(g_d, g_d)` applies,
  so the two-clade `default=0` idiom excludes self-transfer too, consistently.
- **One rule per slot.** `Clades` cannot combine with `"distance"` or `DrivenBy` in the same
  `transfer_to`, exactly as those cannot combine today. The existing "takes one recipient rule" error
  covers it.
- **Empty / mismatched kernel.** A `Between` naming only pairs that never co-occur leaves every
  transfer a no-op — the choice-slot analogue of `check_mapping_fires`. Add a `check_kernel_fires`:
  raise if no named pair's two groups both exist.
- **Determinism.** Membership is a deterministic precompute; the pick uses the existing
  `_weighted_index`. Byte-identical to today when `transfer_to` is left at `"uniform"`.

## What changes, file by file

| File | Change |
|---|---|
| `zombi2/rates/mapping.py` | add `Between` (kernel: `{(from,to): weight}`, `default`, `.weight(a,b)`, string-keyed, validated like `Table`); export it. **Not** a `Mapping`. |
| `zombi2/genomes/_transfer.py` | add `@dataclass Clades` (groups → MRCA subtrees, a `Between`); a `resolve_groups(tree, groups) -> {node_id: label}` precompute; extend `recipient_index` with the `Clades` branch and the donor-conditioned `Between` branch for `DrivenBy`. |
| `zombi2/genomes/__init__.py` | `transfer_to` validation: accept `Clades`; thread the precomputed group map into `_do_transfer`/`recipient_index`; when `transfer_to` is `DrivenBy` **with a `Between` mapping**, pass the donor value. Add `check_kernel_fires`. Export `Clades`, `Between`. |
| `zombi2/rates/parse.py` | add `Between` to `_NAMES` so `DrivenBy(file, Between({...}))` is spellable on the CLI/TOML. (`Clades` with tip-lists stays Python-first — see CLI note.) |
| `manual/book/ch4.md` | `transfer_to` gains a third bullet: `Clades(...)` — topological, between named clades. A worked "two clades" example. |
| `manual/book/ch9.md` | the donor-conditioned *coupling* row: `DrivenBy(trait, Between(...))`; state plainly that clade targeting is **not** here (it reads the tree — Ch4). |
| `docs/design/SPEC.md §5` | one line: the choice slot's weight may be a **kernel** over (donor, recipient) groups, not only a per-recipient factor. |
| `manual/book/appendix-b.md` | no new output files; the `genome_events.tsv` transfer rows already carry `donor`/`recipient`, which is what an audit of "did it stay between the clades" reads. |
| `tests/` | new: `Clades` two-clade → 100 % between / 0 % within / 0 % rest; directional; `default` semantics; disjointness + unknown-tip + empty-kernel errors; `DrivenBy(trait, Between)` donor-conditioning; determinism vs `"uniform"` baseline. |

## CLI and `--params`

The **driven** form is CLI-ready once `Between` is whitelisted:
`--transfer-to "DrivenBy('habitat.tsv', Between({('marine','marine'): 3.0}))"`.

The **topological `Clades`** form takes tip-lists, which do not sit naturally on a flag — and `Distance`
already sets the precedent that a parameterised topological rule is Python-first (the CLI exposes only
the `"distance"` string, not `Distance(decay=)`). **Recommendation (matches the "experimental =
Python-first, promise-what-you-keep" rule in memory): ship `Clades` in the Python API first.** A later
CLI surface, if wanted, is a small **groups file** (`group<TAB>tip` rows) plus
`--transfer-to "Clades('groups.tsv', Between({...}))"` — additive, no new subcommand. **OPEN** whether
to do that now or defer.

## Recommended build order

1. **`Between` kernel** (`mapping.py`) + its `check_kernel_fires`. Small, self-contained, tested alone.
2. **`Clades` topological rule** (`_transfer.py` + `__init__.py` wiring). This *is* the user's ask:
   self-contained, no file, conceptually correct. Ship + document (Ch4).
3. **`DrivenBy(trait, Between)`** donor-conditioned coupling + Ch9 + CLI whitelist. Natural follow-on,
   reuses everything from 1–2.

Steps 1–2 fully answer the opening request. Step 3 generalises it to trait-painted groups.

## Open questions for Adrián

1. **`Between` default weight** — `0.0` (name-the-flows-you-want, the two-clade one-liner) or `1.0`
   (baseline, up-weight pairs)? I lean `0.0`.
2. **Group naming** — MRCA-of-tips, a bare internal node id, or both? I lean **both** (tips for humans,
   id for scripts).
3. **Name** — `Clades` (true monophyletic groups; the engine can verify) vs `Groups` (any partition,
   incl. paraphyletic sets)? I lean `Clades` for honesty, with `Groups` reserved if paraphyletic
   targeting is ever wanted.
4. **CLI now or later** — ship `Clades` Python-first and add a groups-file CLI surface later, or wire
   the CLI in the same pass?
5. **Scope** — build all three steps, or land `Between` + `Clades` (steps 1–2) first and treat step 3
   as a separate change?

No code beyond this study has been written yet — awaiting your call on the five points above.
