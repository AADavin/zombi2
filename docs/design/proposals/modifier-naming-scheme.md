# Proposal — the modifier families (`On` / `By` / `From`)

**Status: ADOPTED (2026-07-19).** The scheme is in `SPEC.md` §5 (the *modifier families* block), and
the rename is applied across the code and docs: `Time → OnTime`, `Diversity → OnTotalDiversity`,
`Inherited → FromParent`, `ByBranch → ByLineage`; `Speed` retired (per-family = `ByFamily`). This
document stays the **full analysis** SPEC §5 links to. (`sequence-api.md` already carried `ByLineage`
from PR #183; `genome-api.md` already retired `Speed`.) Propagation is a cross-file
rename (`zombi2/rates/modifiers.py` + every importer + `SPEC.md` + the level design docs) that touches
**shipped species code** and **multiple open PRs** (the sequences build uses `Inherited`), so it must
be *sequenced*, not landed piecemeal — see **Propagation** at the end.

Designed with Adrián on 2026-07-19, out of the trait-level modifier work (the trait `rate` took
`Time` / `Inherited` / `Diversity`). The question was whether a consistent naming scheme would read
more clearly across the four levels. It does — with honest limits, recorded below.

It also fits the intended **introduction arc for ZOMBI2 — levels → rates → modifiers → joining /
conditioning** (`SPEC.md` §1 → §5 → this §5 addendum → §2): once the reader has the four levels and
the `scope(base) × modifiers` rate, the preposition families are the natural next beat, and joining /
conditioning close the exposition.

---

## The proposed §5 addendum (the text that would go into `SPEC.md`)

> ### Modifiers name their family with a preposition
>
> A modifier is a dimensionless multiplier on **one** rate (§5: *how fast*, never *how many*; never
> starts with "per"). Its **name begins with the preposition that fixes its family**:
>
> | Preposition | Family | The factor is… | Members |
> |---|---|---|---|
> | **`On`** | covariate | a deterministic function of a measured quantity | `OnTime`, `OnTotalDiversity` |
> | **`By`** | independent | an i.i.d. draw, one per unit — **no memory** (uncorrelated) | `ByLineage`, `ByFamily` |
> | **`From`** | inherited | inherited along a genealogical edge — **continuous memory** (autocorrelated) | `FromParent` |
>
> The preposition *is* the mechanism: `On` responds to something measurable, `By` draws fresh per
> unit, `From` inherits from an ancestor. So the uncorrelated / autocorrelated relaxed-clock split is
> just `ByLineage` vs `FromParent`; and ClaDS, the autocorrelated clock, and variable-rates BM are one
> modifier — `FromParent` — at three levels. `ByLineage` names the unit that draws independently (one
> value per lineage — a whole species runs hot or cold), matching the scope `PerLineage`.
>
> Four rules:
>
> - **Fully-qualify an `On` covariate.** The preposition does not fix the covariate's *scope*, so the
>   name must: `OnTotalDiversity` (whole-tree standing diversity), not a bare `OnDiversity` that could
>   equally mean clade-local or (under an SSE model) state-specific. `By` / `From` name a concrete
>   unit / source and need no qualifier.
> - **One memory-structure per axis.** On the tree axis a rate is `By…` (none), `From…` (continuous),
>   or `Markov` (discrete) — never two at once. Families on *orthogonal* axes compose freely
>   (`OnTime * ByLineage * ByFamily`).
> - **Named exceptions.** Not every mechanism has a preposition. **Discrete memory** — a rate
>   switching between categories via a CTMC — is **`Markov`** (a mechanism name), the third memory
>   structure beside `By` (none) and `From` (continuous). Name any future prepositionless mechanism
>   the same way.
> - **Outside the grammar.** A modifier multiplies *one* rate. A construct that evolves a *value*
>   instead of a rate (the OU trait's `reverts_to` / `pull`) is **not** a modifier — it is an
>   argument, and takes no preposition. (`From…` may carry `reverts_to` / `pull` to mean
>   *rate*-reversion — the CIR clock — sharing only its knob-*names* with the OU trait, not its
>   mechanism.)
>
> **The rename:** `Time → OnTime`, `Diversity → OnTotalDiversity`, `Inherited → FromParent`,
> `ByBranch → ByLineage`; `ByFamily`, `Markov` stand. `Speed` is **dropped** — per-family
> heterogeneity is `ByFamily`. (`ByBranch` is then free for its literal future meaning — a
> per-*gene-tree*-branch clock, deferred as exotic in `sequence-api.md`.)

---

## Supporting material (context for the review — would *not* go into SPEC)

### Level map — the scheme across all four levels

`*` = hypothetical / future, included to probe the grammar's reach.

| Modifier | Family | Species | Genomes | Sequences | Traits |
|---|---|---|---|---|---|
| `OnTime` | covariate | episodic diversification | time-varying D/T/L | time-varying substitution | early burst |
| `OnTotalDiversity` | covariate | diversity-dependent diversification | transfer ∝ contemporaries | (rare) | ecological limits |
| `OnAge`\* | covariate | age-dependent speciation | — | — | (rare) |
| `ByLineage` | i.i.d. | uncorrelated rate heterogeneity | per-lineage D/T/L | **uncorrelated clock (UCLN)** | uncorrelated variable rates |
| `ByFamily` | i.i.d. | — *(no families)* | per-family D/T/L | per-gene rate | — *(no families)* |
| `ByBranch`\* | i.i.d. | — | — | per-*gene-tree*-branch clock *(exotic, deferred)* | — |
| `FromParent` | inherited | **ClaDS** | family-rate drift | **autocorrelated clock** | variable-rates BM |
| `FromDonor`\* | inherited | — *(no HGT)* | transferred gene keeps donor's rate | donor-inherited clock | — |
| `Markov` | *(exception)* | clade-shift categories | — | discrete-category clock | rate-shift categories |

The `From` family is not a singleton: **`FromDonor`** (a horizontally-transferred gene inheriting its
*donor's* rate, not its parent's) is a real second member at the genome / sequence levels — the
preposition generalises to any genealogical edge, vertical or horizontal.

**`ByBranch → ByLineage`** (raised by Adrián, also under discussion on the sequences PR #183): the
clock is *"a property of a lineage — one value per species branch, shared by all its genes"*
(`sequence-api.md`). So the unit is the **lineage**, not an unqualified "branch" (species-tree branch?
gene-tree branch?), and `ByLineage` both matches the scope `PerLineage` and lets `sequence-api.md`
drop its disambiguation paragraph. It also *reserves* `ByBranch` for the literal per-gene-tree-branch
clock that doc defers as exotic — a real future modifier, distinct from `ByLineage`.

### Limitations the scheme carries (acknowledged, not hidden)

1. **Not a closed taxonomy.** `Markov` (discrete memory) has no honest preposition — the real scheme
   is "`On` / `By` / `From` **+ mechanism names** for the rest." Three clean families, one named
   exception; state it, don't oversell it.
2. **`On` underspecifies scope.** `OnDiversity` could mean total / clade-local / state-specific;
   `OnAge` lineage-age vs tree-age. The precision lands on the covariate word, and several covariates
   have more than one sensible reading — hence the *fully-qualify* rule.
3. **Blind to cross-rate coupling — now moot.** `By` / `From` are the *tree-memory* axis; whether a
   lineage's *several* rates move together vs independently is an orthogonal axis the prepositions
   cannot express. This was the reason `Speed` sat outside the grammar — but `Speed` is now **dropped**
   (per-family heterogeneity is `ByFamily`), so the design carries no cross-rate construct and the gap
   is only theoretical.
4. **Composition is not policed.** The grammar happily multiplies `ByLineage * FromParent` — two
   mutually-exclusive memory-structures on one axis, conceptually incoherent. The scheme makes
   families legible; an engine rule (or validation) still has to enforce "one memory-structure per
   axis."
5. **Value-level seam.** OU (`reverts_to` / `pull`) evolves the *value* and is a function argument;
   `FromParent(reverts_to, pull)` = CIR reverts a *rate*. Same knob-names, two homes — the name alone
   won't tell a reader which they are touching.
6. **Minor:** `From` shadows Python's `from` (as `Global` shadows `global` — livable, capitalised);
   and membership is unbalanced (`On` ~3, `By` ~2, `From` ~1–2), so it reads as a *system*, not a
   *symmetry*.

### Verdict

A strong **spine**, not a total grammar: it names the three smooth families durably across every
level, respects §5's existing constraints, and even anticipates horizontal inheritance. Adopt
`On` / `By` / `From` as the naming spine, keep `Markov` as the one explicitly-named exception
(discrete memory), and require `On` covariates to be fully qualified.

---

## Propagation (how the rename lands safely)

**What it touches** — one mechanical rename, wide surface:

- `zombi2/rates/modifiers.py` — rename the classes (`Time → OnTime`, `Diversity → OnTotalDiversity`,
  `Inherited → FromParent`, `ByBranch → ByLineage`) and their `__all__`.
- **Every importer + its tests** — `zombi2/species` (shipped), `zombi2/genomes`, `zombi2/traits`,
  `zombi2/sequences` (in flight, uses `Inherited` + `ByBranch`), and `tests/test_modifiers.py` + each
  level's tests.
- **The words** — `SPEC.md` (§5 addendum + the §6 vocabulary row), `MAP.md` (the `modifiers.py`
  line), and the level design docs: `sequence-api.md` (names `ByBranch` → `ByLineage`, and can drop
  its disambiguation paragraph), `genome-api.md` (also drops `Speed`), `trait-api.md`. Any manual
  chapter that names a modifier.

**The conflict problem.** Every open PR that imports these names (this traits PR; the sequences PR,
which uses `Inherited`; and any other in-flight level work) will collide with a rename to `main`. So
it cannot be dropped in piecemeal.

**Recommended sequence** (least churn):

1. **Ratify the words** — move this addendum into `SPEC.md` §5 and update the §6 vocabulary row.
   (Words only; no code — low conflict with in-flight PRs.)
2. **Land the in-flight level PRs on their current names** (`Time`/`Diversity`/`Inherited`) — traits,
   sequences, and anything else close to merge.
3. **One atomic rename PR to `main`** once (2) has landed: `rates/modifiers.py` + all importers +
   tests + `MAP.md` + the level docs, in a single commit. With the level PRs already in, there are no
   open-PR conflicts to fight. (I can drive this — it is mechanical and I know the tree — once you say
   the coast is clear.)
4. **New work uses the new names** from that point.

**The one transition question** — new code opened between (1) and (3) still sees the old class names,
so it either uses the old names (and joins the rename surface) or we ship a **short-lived alias**
(`Time = OnTime`, …) so new work can adopt the new names immediately. Aliases cut against the clean
core's "one canonical path per name," but there is precedent (`SPEC.md` §12 keeps a deprecated alias
for one release for the `coevolve` rename). Adrián's call: strict (no alias, everyone renames at
step 3) vs. bridged (one-release aliases to smooth the window).
