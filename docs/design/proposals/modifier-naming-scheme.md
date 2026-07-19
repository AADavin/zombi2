# Proposal — the modifier families (`On` / `By` / `From`)

**Status: DRAFT for Adrián's review. NOT yet in `SPEC.md`.** Adopting it is a `SPEC.md` §5 edit
(words first) **plus** a coordinated rename in `zombi2/rates/modifiers.py` and every level that uses
these names — species is **shipped**, and the **sequences build is actively using `Inherited`**. So:
ratify the words here, then propagate serially, coordinating with the sequences instance. Nothing on
this branch renames the shared modifiers yet.

Designed with Adrián on 2026-07-19, out of the trait-level modifier work (the trait `rate` now takes
`Time` / `Inherited` / `Diversity`). The question was whether a consistent naming scheme would read
more clearly across the four levels. It would — with honest limits, recorded below.

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
> | **`By`** | independent | an i.i.d. draw, one per unit — **no memory** (uncorrelated) | `ByBranch`, `ByFamily` |
> | **`From`** | inherited | inherited along a genealogical edge — **continuous memory** (autocorrelated) | `FromParent` |
>
> The preposition *is* the mechanism: `On` responds to something measurable, `By` draws fresh per
> unit, `From` inherits from an ancestor. So the uncorrelated / autocorrelated relaxed-clock split is
> just `ByBranch` vs `FromParent`; and ClaDS, the autocorrelated clock, and variable-rates BM are one
> modifier — `FromParent` — at three levels.
>
> Four rules:
>
> - **Fully-qualify an `On` covariate.** The preposition does not fix the covariate's *scope*, so the
>   name must: `OnTotalDiversity` (whole-tree standing diversity), not a bare `OnDiversity` that could
>   equally mean clade-local or (under an SSE model) state-specific. `By` / `From` name a concrete
>   unit / source and need no qualifier.
> - **One memory-structure per axis.** On the tree axis a rate is `By…` (none), `From…` (continuous),
>   or `Markov` (discrete) — never two at once. Families on *orthogonal* axes compose freely
>   (`OnTime * ByBranch * ByFamily`).
> - **Named exceptions.** Not every mechanism has a preposition. **Discrete memory** — a rate
>   switching between categories via a CTMC — is **`Markov`** (a mechanism name), the third memory
>   structure beside `By` (none) and `From` (continuous). Name any future prepositionless mechanism
>   the same way.
> - **Outside the grammar.** A modifier multiplies *one* rate. A construct that couples *several* of a
>   lineage's rates (a family-wide `Speed`) or that evolves a *value* instead of a rate (the OU trait's
>   `reverts_to` / `pull`) is **not** a modifier — it is an argument, and takes no preposition.
>   (`From…` may carry `reverts_to` / `pull` to mean *rate*-reversion — the CIR clock — sharing only
>   its knob-*names* with the OU trait, not its mechanism.)
>
> **The rename:** `Time → OnTime`, `Diversity → OnTotalDiversity`, `Inherited → FromParent`;
> `ByBranch`, `ByFamily`, `Markov`, `Speed` stand.

---

## Supporting material (context for the review — would *not* go into SPEC)

### Level map — the scheme across all four levels

`*` = hypothetical / future, included to probe the grammar's reach.

| Modifier | Family | Species | Genomes | Sequences | Traits |
|---|---|---|---|---|---|
| `OnTime` | covariate | episodic diversification | time-varying D/T/L | time-varying substitution | early burst |
| `OnTotalDiversity` | covariate | diversity-dependent diversification | transfer ∝ contemporaries | (rare) | ecological limits |
| `OnAge`\* | covariate | age-dependent speciation | — | — | (rare) |
| `ByBranch` | i.i.d. | uncorrelated rate heterogeneity | per-branch | **uncorrelated clock (UCLN)** | uncorrelated variable rates |
| `ByFamily` | i.i.d. | — *(no families)* | per-family D/T/L | per-gene rate | — *(no families)* |
| `FromParent` | inherited | **ClaDS** | family-rate drift | **autocorrelated clock** | variable-rates BM |
| `FromDonor`\* | inherited | — *(no HGT)* | transferred gene keeps donor's rate | donor-inherited clock | — |
| `Markov` | *(exception)* | clade-shift categories | — | discrete-category clock | rate-shift categories |
| `Speed` | *(argument)* | — | family-wide speed | — | — |

The `From` family is not a singleton: **`FromDonor`** (a horizontally-transferred gene inheriting its
*donor's* rate, not its parent's) is a real second member at the genome / sequence levels — the
preposition generalises to any genealogical edge, vertical or horizontal.

### Limitations the scheme carries (acknowledged, not hidden)

1. **Not a closed taxonomy.** `Markov` (discrete memory) has no honest preposition — the real scheme
   is "`On` / `By` / `From` **+ mechanism names** for the rest." Three clean families, one named
   exception; state it, don't oversell it.
2. **`On` underspecifies scope.** `OnDiversity` could mean total / clade-local / state-specific;
   `OnAge` lineage-age vs tree-age. The precision lands on the covariate word, and several covariates
   have more than one sensible reading — hence the *fully-qualify* rule.
3. **Blind to cross-rate coupling.** `By` / `From` are the *tree-memory* axis. Whether a lineage's
   *several* rates move together (`Speed`) vs independently (`ByFamily`) is an orthogonal axis the
   prepositions cannot express — which is why `Speed` sits outside the `scope × modifiers` grammar as
   an argument.
4. **Composition is not policed.** The grammar happily multiplies `ByBranch * FromParent` — two
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
`On` / `By` / `From` as the naming spine, keep `Markov` (and cross-rate `Speed`) as explicitly-named
exceptions, and require `On` covariates to be fully qualified.
