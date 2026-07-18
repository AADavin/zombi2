> ⚠️ **SUPERSEDED by [SPEC.md](SPEC.md).** This document encodes the old lexicon (e.g. coevolve / diamond / opportunity / propensity) and is kept only for history. Do NOT treat it as current — SPEC.md is authoritative.

# Coevolution → 4-Level Diamond Migration Guidelines

*Design/engineering guideline for re-packaging the `zombi2.coevolve` subsystem onto one unified coupling grammar over the S/T/G/Σ diamond. Target home for this document: `docs/design/coevolve-grammar.md` (lowercase-kebab, per D2).*

> **Implementation status** (branch `claude/coevolve-grammar-migration`) — the *disjoint* slices (independent of the in-flight `Rates(per=…)` rate rename) are landing:
> - ✅ **Grammar core** — `zombi2/coevolve/grammar.py`: `Response` (`Scalar`/`Table`/`Curve`), `Driver`, `TargetVariable`, `Coupling`, `CouplingGraph` (layer/fuse via SCC), the null legality layer (`neutral`/`cid`/`timing`), `solve_plan`, the `DriverSignal` protocol.
> - ✅ **Genome rate bridge** — `rate_bridge.py`: `CouplingModifier` (a `Modifier`); verified `null == base` byte-identical.
> - ✅ **P1 trait–gene pair reframed** — `genomes:traits` (`walk_optimum_coupled_trait`, **byte-identical**) and `traits:genomes` (`TraitGeneRates` → a `ModifiedRates` over `PerCopyRates` + `CouplingModifier`, loss values byte-exact, per-family gain = statistical identity per §6).
> - ✅ **Sequence tier** — `sequence_bridge.py`: `DriverClock` (`T→Σ` substitution-speed) · `OmegaSelector` (`T→Σ` selection, ω/dN-dS) · `GeneEventOmega` (`G→Σ` selection, post-event relaxed selection). ω rides a small `model_for` hook added to `evolve_on_tree` (backward-compatible) with an ω-class cache.
> - ✅ **P2 single-direction edges** — `species_bridge.py`: `traits:species` (SSE, `Table`→per-state birth/death) and `genes:species` (key innovation, `Scalar` exp-link) grow the tree (the **fuse/into-species** solver); `species:traits` and `species:genomes` (cladogenetic **overlays**, event-driver→jump) reuse the `Cladogenesis`/`CladogeneticGenome` engines. Rebased onto current main (rate rename PR #151), using canonical `Rates(per="copy")`.
> - ⚪ **Remaining** — the joints (ClaSSE / co-diversification / trait-gene feedback) already exist as composable functions; a formal grammar `Jump` response for event drivers; the `G↔Σ` showcase (deferred); and CLI/public-API wiring of all the new edges.
> - ⚪ **Remaining follow-ups** — `G→Σ` substitution-speed (per-node `_annotate` factor), the `G↔Σ` concerted-evolution showcase (research-grade, deferred), and CLI/public-API wiring of the new edges.

---

## 1. Purpose & scope

### 1.1 Goal

Replace the current ~3 per-edge "dialects" in `zombi2/coevolve/` (each edge class carries its own bespoke constructor vocabulary — `effect_loss`, `driver_speciation`, `theta_present/absent`, `theta_low/high`, `cladogenetic_loss/gain`, `birth/death`, …) with **one grammar**: every coupling is the sentence

> **driver → target-variable : response**

and two rules read the sentence (a solve rule and a topology rule). At the same time, extend the current three-node **triangle** {species, traits, genomes} to the four-node **diamond** {species (S), traits (T), genomes (G), sequences (Σ)}, formalising the sequence tier that already exists in the codebase (`zombi2/sequences/`, `zombi2/experimental/`) but is not yet expressed as coupling edges.

### 1.2 Guiding principles (weave these through every PR)

- **Generality lives in the ONTOLOGY, not the ENGINE.** A general language does not require a general engine. We keep the existing solvers — the forward exact-Gillespie (`gene_diversification._simulate_once`, `sse._simulate_sse`/`_simulate_quasse`) and the small-time-step co-integrator (`trait_gene_feedback._evolve_branch`) — and only replace the *front door*.
- **Adding a level adds one noun + its target-variables.** It does not change the sentence or the two rules. Σ is the proof: adding it is near-zero engine code.
- **Clarity dies from DIALECT sprawl, not level count.** The migration's entire point is to collapse the dialects into one grammar, not to add models.
- **This is RE-PACKAGING, not a science rewrite.** Every mechanism/matrix/RNG-draw-order is preserved. Neutral runs must stay byte-identical; coupled runs must stay statistically identical under fixed seeds.

### 1.3 Explicit non-goals

- **Populations (a 5th level) are OUT OF SCOPE.** See the single deferred paragraph in §5.5. No populations code, flag, or noun may enter this migration.
- **No new science.** No new mechanism ships except the one G↔Σ concerted-evolution *showcase* (§5.3), which is itself a reframe over **core** codon ω models + core `ConversionModel`.
- **SSE class names are NOT renamed.** `BiSSE/MuSSE/HiSSE/QuaSSE/CID` stay (ratified C8; field-standard, documented as the `traits:species` edge).
- **Gene-stem UNIT class names are NOT genome-ified.** `GeneDiversification`, `CladogeneticGenome`, `TraitGeneCoupling` keep the `Gene` stem (ratified C6 — the *node token* is `genomes`, the *unit* is `gene family`).

---

## 2. The target design ("the grammar")

### 2.1 The sentence

Every edge is exactly:

```
driver  →  target-variable  :  response
```

- **driver** — whose state pushes.
  - **STATE driver**: a value along the tree (a trait value `s`, a gene count / presence). Reaches the engine as a `TraitTrajectory`-style piecewise-constant `value(lineage, t)` signal.
  - **EVENT driver**: a set of instants (speciation). Reaches the engine as a per-split hook.
- **target-variable** — the variable of the target that gets bent. Each level exposes a **fixed, closed list**, read off its process definition (§2.2). Two kinds:
  - a **RATE** (continuous flow: λ, μ, gene loss/gain, OU pull/diffusion, substitution rate, dN/dS ω) → the effect is a **MODULATION** (a multiplier on the rate);
  - a **STATE / level** (trait value, gene presence/copy number, residues) → the effect is a **JUMP**.
- **response** — how the driver value maps to effect size:
  - **SCALAR** (default): `rate = base · exp(strength · driver)` — the exp-link form.
  - **TABLE** (per discrete state): recovers MuSSE / free per-state rates / `theta_present/absent`.
  - **CURVE** (nonlinear continuous, bounded): recovers QuaSSE.
  - A **latent/hidden driver** recovers HiSSE/CID.
  - Interactions/epistasis are opt-in/advanced.

> **What "JUMP" spans.** JUMP is the state-effect *kind*, and it covers two shapes: (a) an **instantaneous at-split state change** (the cladogenetic edges: a daughter drops a gene family, hops a discrete trait, or takes a Gaussian jump), and (b) a **continuous conditioning of the target state's trajectory** — most importantly the OU **optimum** for `genomes:traits` and the T↔G joint, where the driver's current value sets where the trait is pulled, not an instantaneous displacement. Both are "state effects" (contrast with rate MODULATION); do not read JUMP as instantaneous-only.

### 2.2 The closed per-level target-variable menu

Each level exposes exactly these target variables (read off its process). A coupling may only bend a variable on this list.

| Level | Target-variable | Kind | Effect | Backing code |
|---|---|---|---|---|
| **species (S)** | `speciation` (λ) | rate | modulation | `sse._simulate_sse`, `gene_diversification` |
| | `extinction` (μ) | rate | modulation | same |
| | *(cladogenetic split)* — the EVENT-driver source, not a bendable target | event | — | `Cladogenesis`, cladogenetic bursts |
| **traits (T)** | `optimum` (OU θ) | state | jump (optimum-shift) | `GeneConditionedTrait`, `TraitGeneFeedback` |
| | `pull` (OU α) / `diffusion` (σ²) | rate | modulation | `traits.models` |
| | `value` (state itself, at a split) | state | jump (at-split) | `Cladogenesis` |
| **genomes (G)** | `loss` | rate | modulation | `genomes.rates` LOSS |
| | `gain`/`origination` | rate | modulation | ORIGINATION |
| | `duplication` | rate | modulation | DUPLICATION |
| | `transfer` | rate | modulation | TRANSFER |
| | `presence`/`copy-number` (at a split) | state | jump (at-split) | `CladogeneticGenome` |
| **sequences (Σ)** | `substitution-speed` (R_b · s_g) | rate | modulation | `sequences/clocks.py`, `sequences/evolution.py` |
| | `selection` (ω = dN/dS) | rate | modulation | `sequences/codon_models.py` (`gy94`/`mg94`, M-series) |
| | `residues` (state, at an event) | state | jump | (deferred — see §5.4) |

> **Naming hazard baked into this table:** the two core Σ rates are genuinely different — `substitution-speed` (overall subs/site, from `clocks.py`) vs `selection` (ω = dN/dS, from `codon_models.py`). A T→Σ or G→Σ edge MUST name which one it drives. Never let them collapse into an overloaded "rate" (the exact trap the rate-vocabulary consolidation exists to prevent).

### 2.3 NULL = response 0

The null is not a separate mechanism per edge — it is the grammar's **response set to zero**, matched by shared ancestry for free, uniform across ALL edges. This **replaces** the per-edge `null(kind=...)` machinery and the "cid null is a workflow not a transform" wart. The three archetypes become response-layer operations:

- **neutral** = response strength 0 (unifies today's four spellings: `effect_loss=0`, `driver_speciation=0`, `λ0=λ1`, `theta_present=theta_absent`).
- **cid** = response driven by a HIDDEN driver of the same type (matched variance, no observed-state signal). **DECIDED (author-approved):** `cid` means this ONE clean thing everywhere — a matched hidden-driver *transform*, uniform across edges, legal only for **state-driver** edges. The old data-hiding "ground-truth-withholding" benchmark is retained but moves to its own named surface (§4.4), NOT under the `cid` knob. See §4.4 for how this removes the "workflow not a transform" wart.
- **timing** = a cladogenetic (at-split) response re-expressed as an equivalent anagenetic (along-lineage) rate. Legal only for **event-driver** edges (`species:traits`, `species:genomes`).

**Legality is a property of the driver archetype, not per-edge if/error ladders** (see §2.6, §4.4, §6).

### 2.4 Rule 1 — the SOLVE rule (directional vs bidirectional)

Collect edges into a graph:

- **Directional (acyclic)** → run in order (**LAYER**): simulate the driver history first, feed it as a pre-built trajectory into the target.
- **Bidirectional (a cycle)** → run together (**FUSE**). A cycle **includes any arrow into the tree** (everything rides the tree), so *any into-species coupling is bidirectional in disguise* and grows the tree.

The engine picks the concrete solver and **hides it**:
- cycle of rates in a birth-death → exact forward-Gillespie (`gene_diversification._simulate_once`, `sse` engines);
- a continuous OU inside the loop → small-time-step co-integration (`trait_gene_feedback._evolve_branch`, `steps` sub-pieces).

User-facing axis name: **DIRECTIONAL vs BIDIRECTIONAL**. "An arrow into species" ⇒ tree is an OUTPUT (grown, take `--age`/`--tips`); otherwise tree is an INPUT (overlay, `-t`).

### 2.5 Rule 2 — the TOPOLOGY rule (adjacent tier only)

A coupling connects levels **within one tier of each other**; tier-skips are forbidden. Tiers: **S (substrate) → {T, G} (characters) → Σ (sequences).**

- The **forbidden diagonal is species ↔ sequence (S–Σ)**: a sequence rides *gene* trees, not the species tree. To recapitulate the species tree you simulate one gene family with no events. This is *downstream-only*, **not a TODO** — the docs must say so explicitly.
- The **T–G "diagonal" IS allowed** (same tier {T,G}).

### 2.6 The diamond

```
              species (S)          ← substrate / timeline
             /          \
        traits (T) —— genomes (G)  ← characters (ride species tree)
             \          /
             sequences (Σ)         ← ride GENE trees
```

**Five pairs, each = 3 models (forward, reverse, joint):**

| Pair | Status | Literature aliases |
|---|---|---|
| S–T | BUILT | SSE / cladogenetic / ClaSSE |
| S–G | BUILT | key innovation / punctuational / co-diversification |
| T–G | BUILT | trait-linked / gene-conditioned / feedback |
| T–Σ | NEW (reframe) | trait-driven selection / dN-dS |
| G–Σ | NEW (reframe) + 1 showcase | post-duplication relaxation / concerted evolution |
| **S–Σ** | **FORBIDDEN** | — |

### 2.7 Canonical-name-plus-alias policy

- Canonical names are **STRUCTURAL** ("driver shapes target", one verb). The CLI already uses structural tokens (`traits:species`, `genes:species`, …) — **promote these to primary everywhere** (docs, figures, class names) and **demote borrowed names to subtitles/citations**.
- Literature names (SSE, key innovation, cladogenetic, punctuational, ClaSSE, co-diversification, trait-linked, gene-conditioned, feedback) are kept as **searchable aliases** for discoverability.
- Node token is **`genomes`** everywhere (ratified D1/C6); `genes` is a warned, SUPPRESS-ed alias.

---

## 3. The reframing table

Each current edge/class → its grammar expression. **No math changes; only the front door.** Params split into three grammar slots — DRIVER dynamics (the driver's own sub-model), TARGET base/intrinsics (the target's own dynamics), RESPONSE (the coupling coefficient). "Disappears" means folded into a generic slot, kept as a deprecated kwarg alias.

**Result wrappers all PRESERVE their output surfaces.** `GeneDiversificationResult`, `CladogeneticGenomeResult`, `GeneConditionedTraitResult`, `TraitGeneFeedbackResult`, and the renamed `TraitGeneResult` keep their existing adapters (`to_tsv`/`to_newick`/`profile_matrix`/`tip_prevalence`/`panel_occupancy`/`trait_gene_correlation`) and reuse the standard `Genomes`/`ProfileMatrix` writers. The grammar keeps the same output surface; only names change where §6 dictates.

### 3.1 Top triangle (BUILT today)

#### S→T — `species:traits` (`Cladogenesis`, in `zombi2/traits/models.py`)

| Grammar slot | From | Notes |
|---|---|---|
| driver | speciation **EVENT** | — |
| target-variable | trait **STATE** (jump) | |
| response | JUMP kernel: `jump_sigma2` (continuous scale) / `shift` (discrete flip prob) | the "event → state jump" response form |

Reframe hazard: `Cladogenesis` lives in `zombi2.traits`, is **NOT** re-exported by `coevolve/__init__.py`, has **NO** `null()` and **NO** dedicated `simulate_*` (only reachable via `simulate_sse(cladogenesis=...)` / `simulate_traits`). Giving it a first-class grammar surface is *net-new plumbing*, not a pure reframe.

#### T→S — `traits:species` (`sse.py`: `BiSSE`/`MuSSE`/`HiSSE`/`QuaSSE`)

| Class | Current params | Grammar slot mapping |
|---|---|---|
| `MuSSE` | `birth` (len-k), `death` (len-k), `Q` (k×k), `states` | driver=trait STATE (index `i`); target=diversification RATE; response=**free per-state TABLE** (`birth`/`death` → per-state rate table; `Q` → driver's own transition dynamics) |
| `BiSSE` | `lambda0/1`, `mu0/1`, `q01`, `q10` | 2-row Table; `q01/q10` = driver 2-state transition |
| `HiSSE` | `classes` (BiSSE per hidden class), `hidden_transition` | latent-driver Table; `discretize()` = observed-vs-latent projection — **PRESERVE**, but see the hidden-driver-unification note below |
| `QuaSSE` | `speciation`/`extinction` callables, `sigma2`, `rate_bound`, `x0`, `drift` | driver=continuous trait; response=**CURVE** (arbitrary bounded callable); `sigmoid()` = convenience builder; `rate_bound` carried for thinning |
| `CID`, `_CharacterIndependentMuSSE` | class rates + hidden transition | **PRESERVE** — the honest-null construct the grammar's `cid` layer produces |

Key contrast the grammar's response spec **must** support: SSE hard-codes a **free per-state Table**, while `traits:genomes`/`genes:species` hard-code an **exp-link Scalar**. Forcing SSE into exp-link would silently narrow expressiveness. `birth/death` → rate Table; `Q` → driver dynamics.

**HiSSE-the-model vs cid-null-the-construct (resolution).** `HiSSE`, `CID`, and `--null cid --hidden` are all instances of ONE underlying primitive: a **hidden-class driver variable** (a latent dimension of the driver whose value modulates the target). The grammar unifies them at the *primitive* level — the cid null (§4.4) inserts a matched hidden driver, and `HiSSE` is the same primitive with user-supplied per-class rates. **Both user-facing surfaces persist deliberately:** `HiSSE` stays a named, citable model (Beaulieu & O'Meara 2016) for discoverability and is *not* renamed (C8); `cid` stays the null *operation*. They are documented as two doors onto the same hidden-driver primitive, not two mechanisms.

Joint S↔T (ClaSSE) exists **only** as the `cladogenesis=` keyword on `simulate_sse` — no `simulate_classe`. Grammar makes it an explicit two-edge composition (T→S rates + S→T jump).

#### G→S — `genes:species` (`gene_diversification.py`: `GeneDiversification`)

| Current param | Grammar slot |
|---|---|
| `n_drivers` | number of driver variables |
| `driver_speciation` (βλ, scalar or len-K), `driver_extinction` (βμ) | **RESPONSE** — Scalar exp-link coefficients on λ/μ (log-additive over present drivers) |
| `lambda0`, `mu0` | target base rates (shared with `traits:species` — see reuse) |
| `loss`, `origination`, `transfer` | **DRIVER** gene dynamics (transfer is frequency-dependent `t·carriers/(n−1)`) |
| `root_drivers` | driver initial condition |
| `cladogenetic_loss`, `cladogenetic_gain` | belong to the **reverse** `species:genomes` half of the joint |

**UNIFY the coefficient naming**: `driver_speciation`/`driver_extinction` (β) vs `effect_loss`/`effect_gain` vs `birth`/`death` all mean "coefficient on a rate."

#### S→G — `species:genomes` (`cladogenetic_genome.py`: `CladogeneticGenome`)

| Current param | Grammar slot |
|---|---|
| `initial_families` | target size / IC |
| `loss`, `origination` | target's own anagenetic dynamics |
| `cladogenetic_loss` (per-family drop prob), `cladogenetic_gain` (Poisson mean of new families) | **RESPONSE** — the at-split burst (event → state jump), same family as `species:traits` |

UNIQUE: has a **`timing` null** (analytic branch-spread; requires `tree=`). Grammar's null layer must support `timing` **only for event-driver edges**; `cid` is meaningless here (currently raises `ValueError`).

#### T→G — `traits:genomes` (`trait_coupling.py`: `TraitGeneCoupling` + `TraitGeneRates`)

| Current param | Grammar slot |
|---|---|
| `n_families`, `weights`, `state_values`, `prefix` | driver-value map + target panel (PRESERVE) |
| `base_loss`, `transfer`, `duplication`, `origination` | target base rates on LOSS/TRANSFER/DUPLICATION/ORIGINATION |
| `effect_loss`, `effect_gain` | **RESPONSE** — Scalar exp-link coefficients: `loss = base_loss·cn·exp(−effect_loss·w_i·s)`, `gain = transfer·n·exp(effect_gain·s)` |

`TraitTrajectory` (piecewise-constant `value(branch,t)` + `refresh_times`) is **PRESERVE-verbatim** — the generic STATE-driver signal for *any* overlay edge. `TraitGeneRates` is **REWRITE** — the monolithic special case the grammar subsumes (formerly `TraitLinkedRates`, now a deprecated alias); it must be regenerated from a declarative response spec that emits the same `EventWeight` list (byte-identity hazard: it never uses the aggregate `family=None` fast path; port must preserve weight values, `_clamp`/`_MAX_EXPONENT=40.0`, and weight ordering).

**Output-surface hazard (do not lose):** `traits:genomes` uses the genome-style `--write` surface — `Genomes.WRITE_PARTS` choices, `--sparse`, `--annotate-species`, default `[profiles, trees]` — *unlike* every other edge, which writes a fixed file set. The unified grammar output must PRESERVE this `--write` flexibility for this edge; do not flatten it into a fixed file list.

**Ratified rename (C8/PR6, ALREADY APPLIED on current main):** the primary names are `TraitGeneRates` / `TraitGeneResult` / `simulate_trait_conditioned_genomes` (`trait_coupling.py:313/367/418`). The old spellings resolve as deprecation-warned aliases: `coevolve/__init__.py` `__getattr__` reads them off `_DEPRECATED_ALIASES` (`TraitLinkedRates`/`TraitLinkedResult`), plus module-level aliases in `trait_coupling.py` (`TraitLinkedRates = TraitGeneRates`, etc.). Nothing left to apply here — the rename is done.

#### G→T — `genomes:traits` (`gene_conditioned_trait.py`: `GeneConditionedTrait`)

| Current param | Grammar slot |
|---|---|
| `gene_gain`, `gene_loss`, `root_gene` | **DRIVER** (modifier gene) own dynamics (2-state Mk) |
| `theta_absent`, `theta_present` | **RESPONSE** — free 2-entry per-state **TABLE** (contrast the exp-link of `traits:genomes`); sets the OU optimum (an optimum-shift JUMP, §2.1) |
| `alpha`, `sigma2`, `x0` | target trait's intrinsic OU params |

### 3.2 Three joints

#### T↔G joint — `simulate_trait_gene_feedback` (`trait_gene_feedback.py`: `TraitGeneFeedback`)

Bidirectional; tree is INPUT (neither arrow into S). Composes a `traits:genomes` response and a `genomes:traits` response on the SAME objects, co-integrated in `steps` sub-pieces. `effect_loss=0` recovers `genomes:traits`; `theta_high=theta_low` recovers `traits:genomes`.

**NAMING HAZARD (must reconcile):** `theta_low/theta_high` here (a 2-endpoint *interpolation* `optimum(m)=theta_low+(theta_high−theta_low)·m/N`) vs `theta_absent/theta_present` in `gene_conditioned_trait` (a *lookup*); and `gain` here vs `transfer`/`origination` elsewhere. Unify under the Table-response vocabulary. **PRESERVE** `trait_gene_correlation()` — the headline diagnostic.

#### S↔G joint — `simulate_co_diversification` (`gene_diversification.py`)

Bidirectional; tree GROWN. `genes:species` rates (forward) + `species:genomes` cladogenetic burst (`_burst_drivers` reshuffles the same fixed driver panel at each split). Requires `is_co_diversification` (a `cladogenetic_*>0`) else `ValueError`. Fold the `is_co_diversification` guard into grammar validation. **Byte-identity oracle:** both bursts 0 ⇒ identical newick+tsv to plain `genes:species` (`test_codiv_reduces_to_genes_species_when_no_burst`).

#### S↔T joint — ClaSSE (`simulate_sse(cladogenesis=...)`)

Bidirectional; tree GROWN. T→S state-dependent rates + S→T at-split jump. Grammar makes the two-edge composition explicit. `root_state` defaults to the stationary distribution.

### 3.3 Reframing summary — the coefficient/response rename map

| Old name(s) | Concept | Grammar slot |
|---|---|---|
| `effect_loss`, `effect_gain`, `driver_speciation`, `driver_extinction` | log-linear coefficient on a rate | Scalar RESPONSE coefficient |
| `birth`, `death` (per-state vectors) | free per-state rate | Table RESPONSE |
| `theta_absent`/`theta_present`, `theta_low`/`theta_high` | per-driver-state target value | Table RESPONSE |
| `cladogenetic_loss`/`gain`, `clado_shift`/`jump`, `jump_sigma2`/`shift` | at-split jump magnitude | JUMP RESPONSE (event driver) |
| `speciation`/`extinction` callables + `sigmoid` | continuous driver→rate map | Curve RESPONSE |
| `loss`/`origination`/`transfer`/`duplication`, `gene_gain`/`gene_loss` | the driver's or target's OWN dynamics | DRIVER dynamics / TARGET base — **not** the coupling |
| `alpha`/`sigma2`/`x0`, `q01`/`q10`, `Q` | intrinsic process params / driver transition | DRIVER dynamics / TARGET intrinsics |

---

## 4. New core to build

### 4.1 The four objects

Build a small declarative core (proposed module `zombi2/coevolve/grammar.py`), each object thin over existing infrastructure:

- **`Driver`** — wraps a driver sub-model + how it enters the engine.
  - STATE driver → produces a `TraitTrajectory` (reuse `trait_coupling.TraitTrajectory` **verbatim**; conceptually a `DriverTrajectory` — a driver need not be a trait). `value(lineage,t)` + `refresh_times(t0,t1)` already emit the exact discrete stochastic-map change points.
  - EVENT driver → a per-split hook (speciation instants).
  - Driver's own dynamics = a call into the existing sub-simulators (`traits.models.simulate_traits`/`Mk`, the gene sub-model, `species.forward`).
  - **One driver-variable-type vocabulary.** Today's `--sse-model` and `--trait-model` are two spellings of the same idea — *which stochastic process the driving variable follows*. Unify them into a single driver-variable-type selector (binary / k-state / continuous / Mk / OU / …), routing `_build_sse_model`, `_build_anagenetic_trait`, and `_build_trait_model` through it. `--sse-model` and `--trait-model` become aliases onto this one vocabulary.
- **`TargetVariable`** — a `(level, variable)` pair from the closed menu (§2.2). Resolves to the concrete rate/state seam: gene events via `genomes.rates._resolve_events` (currently D/T/L-only — widen or gate for origination); diversification via the `sse`/`gene_diversification` engines; OU optimum via the trait integrator; Σ via `sequences` (see §5). When `target ≠ diversification` on a given tree, λ/μ are automatically inert — this is exactly the `_build_anagenetic_trait` "strip the diversification, keep only transition/diffusion structure" logic, which becomes automatic in the grammar rather than a special-case builder.
- **`Response`** — one of `Scalar(strength)` / `Table(per_state)` / `Curve(callable, rate_bound)`; carries the `_clamp`/`_MAX_EXPONENT` guard for exp-links and the `rate_bound` for thinning. `response.strength = 0` ⇒ null (neutral). A hidden-driver `Response` ⇒ `cid`. **Build `Table` and `Curve` first-class, not as exp-link special cases** (see §10.1 risk 1).
- **`Coupling(driver, target_variable, response)`** — the compiled sentence. For a rate target it compiles to a **`Modifier`** (see §4.3); for a state/jump target it compiles to an at-split or trajectory-conditioning hook.

### 4.2 The graph engine

- Collect `Coupling`s into a directed graph over the diamond nodes.
- **Validate** with the two rules: TOPOLOGY (reject S–Σ and tier-skips; allow T–G); SOLVE (detect cycles; an arrow into S is a cycle).
- **DIRECTIONAL (acyclic)** → topological order; simulate each driver history, hand pre-built trajectories forward. Tree is INPUT.
- **BIDIRECTIONAL (cycle)** → FUSE; dispatch to the concrete solver and hide it:
  - rates-only cycle in birth-death → `gene_diversification._simulate_once` / `sse` engines (tree GROWN);
  - continuous OU in the loop → `trait_gene_feedback._evolve_branch` co-integrator.
- **Tree input/output dispatch** is a pure predicate: *is `species` a target-variable in any edge?* Yes ⇒ grow (`--age`/`--tips`); No ⇒ overlay (`-t`). This single predicate replaces the ~14 duplicated `parser.error` guards in `cli.py` (L913-918, 1300-1303, 1336-1338, 1380-1382, 1507-1512, 1542-1546, 1585-1589, 1629-1634, 1668-1673, 1721-1726, 1770-1775).
- **Per-edge output-file contracts are preserved.** The engine dispatches each edge/joint to its existing writer set; downstream tooling depends on the exact filenames (recased per D2, §7). The `traits:genomes` edge keeps its genome-style `--write`/`--sparse`/`--annotate-species` flexibility (§3.1); do not homogenise it into a fixed file list.

### 4.3 Build the rate-target coupling as a `Modifier` (the north-star unification)

The `Modifier` ABC in `zombi2/genomes/rates.py` is *already* "a context-keyed multiplier on an event's rate," composed as `∏factor` by `ModifiedRates`. Today it keys only on a lineage's **own** identity (`branch`, `family`, `time`). **The one genuinely-missing generalisation** is a `CouplingModifier` whose `factor(event, family, branch, time)` reads `driver_trajectory.value(branch, time)` and applies the compiled `Response`. That is a ~15-line `Modifier` subclass — no new composer.

`CouplingModifier` itself does **not** yet exist (grep confirms — genuinely net-new). Everything it composes onto is already present under the exact cited names: the `Modifier` ABC, `ModifiedRates`, `ModifiedRates._per_family`, `Modifier.bind`, and `Modifier.refresh_times` all ship in `zombi2/genomes/rates.py` today. So the only new code is the subclass — no new composer.

Concretely, `traits:genomes` collapses from the monolithic `TraitGeneRates` (still a monolithic `RateModel` subclass today; the reframe onto `ModifiedRates` + `CouplingModifier` below remains this migration's job) to:

```
ModifiedRates(
    PerCopyRates(loss=base_loss, transfer=...),          # target base = base rate × opportunity  (Part 3: → Rates(per="copy"))
    [CouplingModifier(trajectory, response=Scalar(effect_loss))],  # the coupling = a modifier
)
```

Reuse points (cite in code comments):
- `EventWeight(event, family, rate)` — the stable multiplication atom the response multiplies. **PRESERVE.**
- `RateModel.refresh_times` — the exact-stochastic-map seam; a `CouplingModifier` forwards its trajectory's change points here. Already wired end-to-end into `GenomeSimulator._evolve_interval`/`_apply_breaks` — **no simulator change needed** as long as the trajectory exposes `refresh_times`.
- `ModifiedRates._per_family` — splits an aggregate `family=None` weight into per-family weights so a family-keyed response resolves targets.
- `LineageModifier`/`FamilyModifier` — reference implementations of "response as a Modifier"; reuse `FamilyModifier`'s lazy per-key sample+cache pattern for a driver-keyed coupling.
- `state_values` + `_clamp`/`_MAX_EXPONENT` — reusable driver-scalarisation and safe exp guard.
- `bind(rng, tree)` / `time_dependent` — threaded through both `RateModel` and `Modifier`, so a coupling inherits RNG, tree access, and fast/slow-path routing for free.

**Boundaries to encode, not paper over:**
- Joint/feedback edges are **inherently bidirectional** — driver and target co-evolve, so they **cannot** use a pre-built trajectory. They stay on the co-integrator (`_evolve_branch`), outside the thin Modifier layer.
- `_resolve_events` admits only duplication/transfer/loss as scalable; origination mints a new family and is skipped by `_per_family`. Couplings on origination/gain need explicit widening.
- **Recipient/gain-channel seam EXISTS now.** `PairModifier` (donor→recipient transfer highways) is BUILT — `zombi2/genomes/transfers.py:92`, attached to a `TransferModel` via `pair=`, supporting both exact `pairs` `{(donor,recipient): factor}` and clade→clade `blocks` (highways). Emission-seam couplings (loss/retention, duplication, donor-side transfer) are ready via `CouplingModifier`; a driver→recipient coupling builds ON `PairModifier` (key into `TransferModel.pair`), it does **not** wait on a new seam. The one nuance to encode: `PairModifier.factor(donor, recipient)` has a different signature from the emission `Modifier.factor(event, family, branch, time)` — a driver→recipient coupling must adapt to the `(donor, recipient)` signature rather than reuse the emission one.
- Any `ModifiedRates`/coupling routes to the **Python engine** (not the Rust fast path). Acceptable but explicit performance boundary.

**Relationship to the opportunity axis (Part 3).** The rate model is `base × opportunity × modifiers` ([opportunity-knob.md](opportunity-knob.md), [rate-vocabulary.md](rate-vocabulary.md)). A coevolution coupling lives **entirely on the modifier axis** — a `CouplingModifier` is one more dimensionless factor in the `∏factor` product — so it is **orthogonal to the opportunity (`per=`) knob** and never sets it: a lineage can have its `loss` driven by a trait *and* its clock set to `per="copy"|"lineage"|"shared"` independently. Two consequences for this migration:
- **Shared declarative refactor.** Part 3 step 4 ("the simulator owns the count; a model declares `(event, base_rate, opportunity)` instead of baking `r × n` into `event_weights`") is the *same* architectural move as the grammar's declarative `Coupling`/`Response` layer. Co-design them — once rates are declarative data, a `CouplingModifier` reads the count-free base rate and the engine applies opportunity and modifiers uniformly. Do **not** build two parallel declarative layers.
- **Preset → knob naming (incoming).** The `PerCopyRates`/`PerLineageRates` in the sketch are today's classes; Part 3 folds them into `Rates(per="copy"|"lineage"|"shared")` presets (classes retained as deprecated aliases). When that code lands, restate the target base as `Rates(per="copy", loss=…)`; the sketch stays valid meanwhile via the aliases.

### 4.4 One null layer (and the cid-transform-vs-benchmark-workflow split)

Lift the currently-uniform-but-scattered `null(kind=...)` protocol into ONE response-layer keyed on **driver archetype**:

| Driver archetype | `neutral` | `cid` | `timing` |
|---|---|---|---|
| STATE driver (traits, genes) | ✔ response=0 | ✔ hidden driver (transform) | ✘ (`ValueError` "no timing null") |
| EVENT driver (speciation) | ✔ | ✘ (`ValueError` "no cid null") | ✔ (analytic branch-spread; needs `tree=`) |
| continuous STATE (QuaSSE) | ✔ (constant λ) | ✘ (`TypeError` discrete-char) | ✘ |

This **encodes** the legality matrix rather than flattening it — preventing the grammar from resurrecting the very `TypeError('workflow not a transform')` guards these edges deliberately raise.

**Resolving the two-mechanism `cid` hazard.** Today `--null cid` means two genuinely different things, and the migration must NOT collapse them naively:

1. For `traits:species`, cid is a native **transform**: `model.null('cid', n_hidden=...)` inserts a hidden class carrying the rate spread (→ `CID`/`_CharacterIndependentMuSSE`).
2. For the gene/trait edges (`genes:species`, `genes:traits`, `traits:genomes`), cid is a **benchmarking workflow** implemented by the three near-duplicate `_run_*_cid_null` runners: simulate the coupled model, hand back a *neutral observed channel* (a decoupled genome/trait), and *withhold the true driver* as `*_ground_truth.tsv`, plus a `null_manifest.tsv`. This is the "cid null is a workflow not a transform" wart.

The grammar splits these cleanly. **This split is DECIDED (author-approved), not a proposal:**

- **`cid` (the null) IS the uniform hidden-driver transform** for *all* state-driver edges — a matched latent driver of the same type, matched variance, no observed-state signal. For the gene/trait edges this is a genuine transform they did not previously have (net-new but conceptually clean), replacing the mislabeled workflow. This is precisely how the wart is removed: `cid` means one clean thing — a matched hidden-driver transform — everywhere.
- **The ground-truth-withholding benchmark is retained but MOVED to its own distinct named surface** — provisional name `--benchmark neutral-channel` (final token TBD, a naming detail only) — *not* under the `cid` null knob. The three `_run_*_cid_null` runners fold into ONE generic benchmarking-workflow implementation. **Its output contract is preserved byte-for-byte:** the `*_ground_truth.tsv` files (`drivers_ground_truth.tsv` / `modifier_ground_truth.tsv` / `trait_ground_truth.tsv`, recased per D2), the neutral observed channel, the `null_manifest.tsv`, and the stdout "overlay the neutral genome with…" hint all remain — only their trigger moves off `--null cid` onto the named benchmark surface.

`neutral` unifies the four no-effect spellings (§2.3). `timing` re-expresses an at-split response as a matched anagenetic rate (analytic branch-spread; event-driver edges only; needs `tree=`).

**Default-null quirks to reconcile explicitly (do not silently flip):** `BiSSE.null` currently defaults to `kind='cid'` while every other edge defaults to `neutral`; `HiSSE.null` allows only `neutral`; `CladogeneticGenome.null` requires `tree=` for `timing`. Decide per-quirk: preserve the default (safest) or normalise as a *documented* behaviour change. Whichever, encode it in the legality/default table above, not in scattered `if/error` ladders. This is an open question for the author (§10.2).

---

## 5. The sequence tier (the diamond's bottom)

**Present this tier as UNIFICATION, not new inventory.** The mechanisms already exist; the grammar *locates* them as edges.

### 5.1 T→Σ — trait sets selection strength / substitution rate — **BUILD NOW (reframe)**

Two sub-variants, one per **core** Σ target-variable:

- **T→Σ (selection):** a trait drives per-lineage **ω = dN/dS**. Existing core surface: `sequences/codon_models.py` `gy94`/`mg94` (`--omega`), M-series site models (`--omega-model`), `.expected_dnds()`. Today ω is a fixed scalar baked into the 61×61 Q at construction. Reframe: build a codon matrix per lineage keyed on the trait value (relaxed vs purifying). `evolve_on_tree` already takes one model per family; the concrete change is **model-per-branch** with an ω-class cache to stay tractable (per-lineage eigendecomposition is a real cost).
- **T→Σ (substitution-speed):** a trait sets/scales the per-lineage rate **R_b**. Existing surface: `sequences/clocks.py` `Clock.lineage_segments` → `{branch_name: [(rate,t0,t1)]}`, integrated by `SequenceEvolution` (`sequences/evolution.py`). Reframe: a trait-conditioned `Clock` subclass fills the same segment contract; the integrator, pruning, and Newick output are unchanged.

These are the only two core Σ rate target-variables, so a core T→Σ edge drives one of exactly these two. There is no third substrate.

### 5.2 G→Σ — gene event drives sequence — **BUILD NOW (reframe)**

A gene event (duplication/transfer) bumps a family's sequence rate/selection: post-duplication relaxed selection; copy-number → ω. Existing attach points (no new event plumbing):
- `SequenceEvolution._annotate` walks the reconciled gene tree; each `_Node` already carries `node.kind` (EventType DUPLICATION/TRANSFER/LOSS/CONVERSION). A per-event rate/ω multiplier keyed on `node.kind` lands **exactly here**.
- `SequenceEvolution.family_factors` (the sequence analogue of `FamilyModifier`) is a ready multiplicative override slot: `s_g = family_factors.get(fam,1.0) · family_speed.sample()`. A G→Σ coupling writes into this slot, composing without touching the clock. `--family-speeds` already proves the file-driven external-override path end-to-end.
- Gene events are already in scope: `zombi2 sequence` replays `Events_trace.tsv`. The core G→Σ path is exactly the reframe above: `codon_models.py` codon ω + `SequenceEvolution._annotate` (a per-event multiplier keyed on `node.kind`) + `family_factors`, all over the `zombi2 sequence` replay of `Events_trace.tsv`. No new event plumbing.

### 5.3 G↔Σ — concerted evolution — **THE ONE NEW SHOWCASE BUILD**

The single genuinely-new model. Dosage relaxes selection ↔ sequence decay drives gene loss. Ties to the ADH1 / gene-conversion work. **Core substrate (D3-compatible):** the core codon ω models (`codon_models.py` `gy94`/`mg94`, `.expected_dnds()`) evolving coding DNA down gene trees, plus core `ConversionModel` (whole-copy homogenisation, a G event in `Events_trace`). The G-loss half rides the existing gene-content layer; the Σ half reads out as shifted/emergent ω.

- **Constraint from the sequence-coupling lab:** identity is NOT a sufficient statistic under homology-biased conversion — the showcase must **simulate** the homogenised sequences, not infer them from event counts.
- It ships as a **default-shipped core edge** (notebooks are not shipped — author decision 2026-07-18) because it is the one genuinely-new model — it *composes* existing CORE pieces (core codon ω from `codon_models.py` + core `ConversionModel`) into a new coupled behaviour — not because of any experimental dependency.

### 5.4 Deferred edges — with principled reasons

- **Σ→G (decayed sequence → gene loss):** only meaningful *inside* the G↔Σ loop; **not a standalone edge**.
- **Σ→T (a residue flips a discrete phenotype):** **DEFER** — niche; needs a genotype→phenotype map.
- **T↔Σ (molecular phenotype ↔ sequence):** **DEFER** — this *is* the protein-fitness-landscape research program. The framework **locates** it but does not implement it. (A reviewer asking "why not T↔Σ?" is answered by exactly this sentence.)

### 5.5 Forbidden & out-of-scope

- **S–Σ is FORBIDDEN** (topology rule). Sequences ride gene trees, not the species tree. Molecular clocks (`Clock`/`--clock`) are an **S–G–Σ path** — the clock is the *metric on the gene tree*, **not** a direct S↔Σ edge. Codon ω is a T→Σ / G→Σ edge. Present the clock this way to close the S–Σ question.
- **Populations (5th level) — DEFERRED, single paragraph:** The grammar *could* admit populations as the level that makes the substrate stochastic (coalescent/ILS; cf. `experimental ils`, `MultispeciesCoalescent`). It is a separate, massive, later decision and **must not touch this migration**. No populations noun, flag, or code enters here.

### 5.6 Sequence-tier hazards (encode explicitly)

- **RNG/byte-identity:** `SequenceEvolution` draws family speeds in `sorted(families, key=_natkey)` order; the clock is drawn once up front. Any per-event hook that consumes/reorders the RNG stream breaks reproducibility. A coupling must be **deterministic given already-drawn quantities**, or draw from a **separate appended substream**; neutral runs stay byte-identical.
- **ω is baked into Q at construction** — per-lineage/per-event ω ⇒ many matrices ⇒ discretise into ω-classes / cache.
- **`CodonSiteModel` enforces a shared stationary distribution** (`__post_init__` raises otherwise) — a coupling may vary only the non-synonymous weighting, not frequencies.
- **Precondition:** a gene-event-driven Σ edge requires `Events_trace.tsv` (genomes run with `trace` in `--write`) — must fail loudly (matching the existing `FileNotFoundError` guard in `_run_sequence`), never silently no-op.
- **Composition order** for `family_factors × random draw × clock` must be fixed explicitly when a G→Σ coupling also multiplies `s_g`.
- **No-event contrast baseline (core terms).** If a G→Σ showcase wants a "no-event" control, express it in CORE terms: evolve a family's coding DNA down its gene tree with a fixed core codon ω model (`codon_models.py`) and **no** DTL-conditioned multiplier (`family_factors=1`, no `node.kind` hook). That fixed-ω, no-gene-event run is the contrast the coupled G→Σ / G↔Σ variants are measured against (does adding the gene-event coupling shift emergent dN/dS vs this control). No experimental surface involved.

---

## 6. Backward compatibility & deprecation

**Precedent to copy verbatim:** the rate-vocabulary consolidation kept aliases (`SharedRates`→`PerCopyRates`, `PerGenomeRates`→`PerLineageRates`, `BranchModifier`→`LineageModifier`, `--rate-per {copy,genome}`→`{copy,lineage}`) with zero behaviour/RNG change. `docs/design/naming-consolidation.md` P2/C1 specifies the exact mechanism:

- **Python names** — PEP-562 module `__getattr__` map: canonical name in `__all__`; old name **dropped** from `__all__`, resolves for one minor version, emits `DeprecationWarning` naming the canonical. Reuse the existing C1 `__getattr__` map in `coevolve/__init__.py`.
- **CLI flags/value tokens** — kept accepted, `help=argparse.SUPPRESS`, one stderr line via the existing `_deprecated_flag(old,new)` helper.
- **Output filenames** — the one surface with no transparent alias (see §7 for the D2 lowercase recasing).
- **Window:** aliases land 0.3.0, removed 0.4.0. Slot new aliases into the **same** window and the **same** `__getattr__` map — no open-ended parallel alias set.

**Mandatory cross-edge flag-reuse audit (every flag rename).** Flag reuse across edges is pervasive and *asymmetric*; renaming a flag under one edge's grammar slot can silently change another edge's behaviour. Before renaming any flag, audit and document every consumer:
- `--lambda0`/`--mu0` — `traits:species` (SSE base) **and** `genes:species` base rates.
- `--diffusion`/`--root-value` — QuaSSE **and** `species:traits`-alone anagenetic trait.
- `--genome-size` — `species:genomes` **and** two cid-null neutral overlays.
- `--theta-*`/`--trait-*` — `genes:traits`, the T↔G feedback joint, **and** the genes:traits cid null.
- `--panel`/`--effect-loss`/`--loss`/`--trans` — `traits:genomes` **and** the feedback joint.

Renaming any of these requires touching every listed consumer in lockstep (canonical + alias) and a reproducibility check on each affected edge. This audit is a required checklist item on P1–P4 (§9).

**Aliases this migration introduces:**
- `genes`→`genomes` node token (D1/C6) — the only node-token alias this migration still needs to introduce. It is currently **inverted** in the CLI: `_COEVOLVE_EDGES` holds the `genes:*` tokens internally and input `genomes:*` is normalised *down* to `genes` (see §9 P0, item 6). The flip makes `genomes` canonical internally and `genes` the warned input alias.
- `--sse-model`/`--trait-model` → the unified driver-variable-type selector — kept as warned flag aliases.
- Any grammar-level renames of edge dialect kwargs — kept as deprecated constructor kwargs.

The `TraitLinkedRates`→`TraitGeneRates`, `TraitLinkedResult`→`TraitGeneResult`, `simulate_trait_linked_genomes`→`simulate_trait_conditioned_genomes` (C8/PR6) Python aliases already exist on current main (`coevolve/__init__.py` `_DEPRECATED_ALIASES` + `__getattr__`, and module-level aliases in `trait_coupling.py`) — this migration does **not** introduce them.

**Do NOT rename:** SSE class names (`BiSSE/MuSSE/HiSSE/QuaSSE/CID` — C8 exempt); Gene-stem unit classes (`GeneDiversification`, `CladogeneticGenome`, `TraitGeneCoupling` — C6). Do not blanket-sweep `genes`→`genomes` (would corrupt these class names, the `Genomes` result class, and the `Genomes` ancestral-genomes output dir — do surgical replaces only).

**Identity policy — DECISION (Adrián, 2026-07-15): clarity/uniformity over byte-identity.**

The machinery-swap reframes (a bespoke per-edge rate class → the shared `PerCopyRates` / `Rates(per=…)` base + a `CouplingModifier`) change the *order* the base emits events, and the Gillespie's next-event pick is order-sensitive — so a fixed-seed reframed run **diverges** from the current fixed-seed run. Matching the old class's emission order would restore byte-identity, but that nudge is itself edge-specific special-casing — exactly the one-off behaviour the migration exists to delete. So we do **not** chase byte-identity for these reframes.

- **Gate = STATISTICAL identity for the machinery-swap reframes.** The proof is that the *model* is unchanged: exact per-event **formulas** (e.g. `loss = base·cn·exp(-strength·w·s)`), **total rates**, and the **`null == base`** property are checked exactly (the `CouplingModifier` tests already do this), and the analytic + many-σ oracles (ClaSSE jump ~Normal(0,`jump_sigma2`), honest-null suite) stay green. No emission-order matching is added.
- **Consequence to accept:** a fixed-seed reframed run differs event-for-event from the pre-reframe run — same distribution, different single trajectory. Run-twice-same-seed reproducibility (`test_*_reproducible`) still holds (internal consistency); any test that pins *golden* output values for a reframed edge gets its expected values **regenerated** as a documented, model-preserving trajectory change (a re-seed, not a science change).
- **Byte-identity is still required where it is *natural*** — the true algebraic reductions that neither swap machinery nor reorder: the zero-burst `co_diversification`→`genes:species` and default-no-cladogenesis ClaSSE reductions (`test_codiv_reduces_to_genes_species_when_no_burst`, `test_classe_default_no_cladogenesis_matches_phase1`) must still diff to zero.
- **RNG draw order** is still replicated exactly *within* the co-integrator (fused) paths, where nothing is reordered; the accepted divergence is only the emission-order effect of the base-class swap in the *layered* rate paths.

---

## 7. Rewrite / redo list

Cite exact paths. Guide and manual are near-duplicates — **keep in lockstep**.

### 7.1 Docs (prose)

| Path | Change |
|---|---|
| `docs/guide/coevolution.md` | REWRITE: lead with the single grammar rule; demote the nine named models to a lookup menu keyed on `driver:target` (keep citations). Add Σ → the diamond; state the **forbidden S–Σ diagonal** as downstream-only (not a TODO). Replace the `coevolve_modes.svg` triangle ref with the diamond figure. Change all node tokens to `genomes`. Document the cid transform vs the ground-truth-withholding benchmark workflow split (§4.4). |
| `docs/guide/coevolution_nulls.md` | Reframe ① neutral as the single `response=0` knob (unify the four spellings). Keep cid/timing as the two matched-variance menu options; clarify cid is now a uniform hidden-driver transform and the ground-truth benchmark moved off `--null cid` onto its own surface. Extend the edge×archetype table with Σ; mark S–Σ N/A. Preserve the `rates.md` cross-link and the `fig_sse_hisse` "honest null" ref. |
| `docs/cli.md` (§222-245, L24) | Update "six edges among three levels" → grammar framing; mention Σ + forbidden S–Σ; `--null` prose gains `response=0` + the cid/benchmark split; note `--sse-model`/`--trait-model` unify. Keep short. |
| `docs/index.md` | Light reframe of the coupling bullet; align `four_levels.svg` caption to the S/T/G/Σ diamond. |
| `docs/comparison.md` (L36) | Reword to unified-grammar phrasing. Minor. |
| `docs/validation.md` (L38-40) | PRESERVE claims (behaviour); at most reword "SSE" → "the `traits:species` entry." |
| `docs/img/four_levels.svg` (+ `_dark`, `levels.svg`, `event_levels.svg`) | Hand-authored (no script) — edit raw SVG to introduce Σ as the fourth node + forbidden S–Σ diagonal; keep dark variant in sync. Manual rebuilds via `make figures`. |

### 7.2 Manual (book)

| Path | Change |
|---|---|
| `manual/chapters/13-coevolution.md` | Same rewrite as the guide, in prose. Factorisation math gains the Σ term + diamond; six-edge table → grammar + menu; null section → `response=0` + cid/benchmark split. Update all embedded `.pdf` refs. **Sync with `docs/guide/coevolution.md`.** |
| `manual/chapters/01-introduction.md` | Name sequences as the fourth node whose S coupling is forbidden (downstream-only); set up the grammar. |
| `manual/chapters/12-trait-evolution.md` | PRESERVE mostly; ensure "related nulls elsewhere" cross-refs (`HiddenStateMk`, `CorrelatedBinary.independent`) still resolve. Peripheral. |

### 7.3 Figures (redraw triangle → diamond)

| Script | Change |
|---|---|
| `figures/scripts/fig_coevolve_modes.py` | **REDRAW as the diamond** (headline change): add Σ node (fourth hue), draw S–Σ as the dashed/greyed forbidden diagonal, keep S/T/G arrows + joint double-arrows. Reuse `node()/arrow()/biarrow()/_head()` + teal/sage/terracotta palette. Cascades to every embedding. |
| `figures/scripts/fig_coevolve_nulls.py` | Relabel NEUTRAL panel → `response=0`. Geometry unchanged. |
| `figures/scripts/fig_coevolve_null_timing.py` | Caption only → `timing` menu option. |
| `fig_sse.py`, `fig_sse_quasse.py`, `fig_sse_cladogenetic.py`, `fig_sse_hisse.py`, `fig_key_innovation.py`, `fig_gene_conditioned_trait.py`, `fig_trait_linked_genes.py`, `fig_punctuational_genome.py` | PRESERVE geometry; recaption as **instances** of the grammar (not standalone "models"). `fig_sse_hisse` stays the canonical `cid` figure (the hidden-driver primitive); `fig_sse_cladogenetic`/`fig_punctuational_genome` tie to the `timing` null. |

Pipeline: `fig_*.py` → `docs/img/*.svg` (+`_bw`/`_dark`) → `make figures` → `figures/*.pdf`. Expect two-step rebuild and non-trivial byte churn; keep bw/dark variants from going stale. Uses `zombi_style.py` + `model_common.py` — Σ node inherits the visual system.

### 7.4 Infra & external

- `mkdocs.yml` (L86-87) — update only if page titles change; else preserve.
- **External manuscript** `CLAUDE/ZOMBI2_workspaces/ZOMBI2_PAPER` (off this worktree/cwd) — its coevolution section + Table 1 capability matrix describe the six-edge triangle and must be rewritten to the grammar+diamond. **Cannot be edited from here — treat as a separate deliverable; easy to forget.**

---

## 8. Test migration plan

The suite is the safety net. Three buckets.

### 8.1 PRESERVE (behaviour oracles — must stay green unchanged)

Run before/after each PR on the **same seeds**:
- Every `test_*_reproducible` (byte-identical `to_newick`/`to_tsv`/matrix): `test_gene_div_reproducible`, `test_cladogenetic_genome_reproducible`, `test_tgf_reproducible`, `test_gene_conditioned_trait_reproducible`, `test_codiv_reproducible`, `test_sse_reproducible`, `test_classe_reproducible`, `test_cid_null_is_reproducible`, `test_driver_shape_and_reproducible`.
- Exact-reduction tripwires: `test_codiv_reduces_to_genes_species_when_no_burst`, `test_classe_default_no_cladogenesis_matches_phase1`.
- The honest-null statistical suite in `tests/test_coevolve_nulls.py` (neutral removes tip bias, CID hidden-drives-diversification, timing analytic `loss=base+clado·nb/L`).
- `tests/test_sse.py` analytic oracles (Yule mean `2·e^{λ·age}`, πQ=0, HiSSE collapse, QuaSSE bounds, ~8σ ClaSSE jumps) — slow to re-baseline; do not perturb the RNG stream.
- Exact-formula units: `test_loss_modulation_formula`, `test_effect_gain_scales_transfer`, `test_tgf_optimum_is_monotone_in_panel`.

### 8.2 UPDATE + add alias tests (API-surface locks)

- `tests/test_api_namespaces.py` — the authoritative rename ledger. Edit `NAMESPACES['coevolve']`/`coevolve.__all__` for any rename; **bump the hardcoded `len(z.__all__)==152`** (`test_api_namespaces.py:147`, current comment reads "151 + SharedBirthDeath"). Add a companion test asserting old names still import (aliases) with a `DeprecationWarning`.
- `tests/test_cli.py` (`test_coevolve_*`, ~L894-1082 — the coevolve section banner is ~L894, the last `test_coevolve_*` def ~L1082) — every `--couple`/`--sse-model`/flag token and asserted filename. Add deprecation-alias tests for renamed flags/tokens. **Preserve** seed-reproducibility asserts. Preserve the default: `zombi2 coevolve --tips 50 -o out` (no `--couple`) = BiSSE. Preserve `--couple`'s dual parsing (`--couple a --couple b` AND `--couple a b`).
- Per-edge API asserts in `test_gene_diversification.py`, `test_cladogenetic_genome.py`, `test_gene_conditioned_trait.py`, `test_trait_gene_feedback.py`, `test_co_diversification.py`, `test_trait_coupling.py`, `test_sse.py` — class/param names, `to_tsv` headers (`'node\tD0\tD1\tD2'`, `'node\tmodifier\ttrait'`), `repr` (`"CID(classes=2)"`).
- `test_trait_coupling.py:246` imports `SharedRates` from `zombi2` — a live alias oracle; keep it, it doubles as the pattern to follow.

### 8.3 Guard traps

- **Ground-truth-file trigger move (§4.4):** output-filename recasing is already DONE on main (D2/C7) — every coevolve output and its test asserts are already lowercase_snake (no capitalised output filenames remain). The one thing still to do here: the ground-truth files (`drivers_ground_truth.tsv`, `modifier_ground_truth.tsv`, `trait_ground_truth.tsv`) and their `null_manifest.tsv` *move trigger* off `--null cid` onto the named benchmark surface (§4.4) — update the invoking flag in the asserting tests (`test_coevolve_nulls.py`, `test_cli.py`), not the filenames.
- **Manifest text** (`'null\tcid'`, `'hidden_classes\t2'`, `'null\ttiming'`) — public surface; update deliberately.
- **Error-message substrings** matched via `pytest.raises(match=...)`: `'still evolve'`, `'exactly four'`, `"no 'timing' null"`, `'unknown null kind'`, `'workflow'`, `"no 'cid' null"`, `'already a hidden-state'`, `'n_hidden >= 2'`, `'needs the tree'`, `'univariate'`, `'out of range'` — reworded grammar messages break these even when behaviour is identical. Update in lockstep.
- **Cross-module alias coupling:** `test_default_rate_model_has_no_refresh_times` imports `SharedRates` — if an unrelated subsystem drops that alias, this coevolve test reddens. Watch.

### 8.4 ADD (new grammar tests)

- **`response=0` ≡ null**: for each edge, `Coupling(..., Response.zero())` produces output identical to the legacy `.null('neutral')`.
- **cid transform ≡ native**: for `traits:species`, the grammar's `cid` matches legacy `model.null('cid', n_hidden=...)`; for gene/trait edges, assert the grammar's `cid` now yields a genuine hidden-driver transform (matched variance, observed char ~uninformative).
- **benchmark workflow preserved**: the named neutral-channel benchmark reproduces the old `_run_*_cid_null` outputs (neutral observed channel + `*_ground_truth.tsv` + `null_manifest.tsv`) byte-for-byte (modulo D2 recasing).
- **Solve-rule equivalence**: a directional edge run through the grammar (LAYER) matches the legacy overlay entrypoint byte-for-byte; a bidirectional edge (FUSE) matches the legacy joint entrypoint.
- **Topology-rule validation**: constructing an S–Σ coupling raises a clear grammar error; a tier-skip raises; a T–G edge is accepted.
- **Null-legality matrix** (§4.4): `cid` on an event-driver edge and `timing` on a state-driver edge raise the encoded errors.

---

## 9. Phased rollout (PR plan)

Smallest-first, low-risk-first. **Hard ordering: nothing starts before the grammar core (P1) lands.** Every PR keeps the full suite green (behaviour oracles unchanged), adds alias tests for any rename, and completes the §6 cross-edge flag-reuse audit for any flag it touches.

**P0 — Prerequisite consolidation (largely DONE).** The C-series (C2–C9) is merged on current main, so most of this prerequisite has already landed: the `TraitLinked*`→`TraitGene*` aliases exist (`coevolve/__init__.py` `_DEPRECATED_ALIASES` + `__getattr__`, and module-level aliases in `trait_coupling.py`), and D2/C7 output-filename lowercasing is complete. The **only** remaining prerequisite is to flip the CLI edge token so `genomes` is the canonical INTERNAL spelling: today `cli.py` `_COEVOLVE_EDGES` (~L990-993) still holds the `genes:*` tokens, and `_normalise_couple_edge` (~L1448-1454) normalises input `genomes:*` **down** to `genes` (the inverse of D1/C6). Flip both so the set holds `genomes:*` and `genes` becomes the warned input alias normalised *up* to `genomes`. No new grammar. Gate: suite green + alias tests.

**P1 — Grammar core + the T–G pair (proof).** Build `Driver`/`TargetVariable`/`Response`/`Coupling` + `CouplingModifier` (§4), with `Table`/`Curve` first-class. Re-express `traits:genomes` (`TraitGeneRates`) and `genomes:traits` (`GeneConditionedTrait`) on the grammar — the most "ad-hoc"-reading edges (exp-link Scalar vs 2-entry Table) prove the response spec spans forms. Preserve `traits:genomes`'s genome-style `--write` surface. Legacy entrypoints become thin wrappers. Gate (**statistical identity**, per §6 — *not* before/after byte-identity): exact per-event formula + total-rate + `null == base` checks (the `CouplingModifier` tests) plus run-twice reproducibility on `test_trait_coupling.py` / `test_gene_conditioned_trait.py`; regenerate any golden-value expectations as a documented model-preserving re-seed; new `response=0` tests; flag audit for `--theta-*`/`--panel`/`--effect-loss`/`--loss`/`--trans`.

**P2 — Into-species edges onto the grammar.** Bring `traits:species` (SSE **free-table** → Table response), `genes:species` (Scalar), `species:traits` (`Cladogenesis` — give it its first-class simulate/null surface), `species:genomes` onto the grammar via the FUSE path (reuse the existing forward-Gillespie engines unchanged). Unify `--sse-model`/`--trait-model` into the driver-variable-type selector. Merge `HiSSE`/`cid`/`--hidden` onto the one hidden-driver primitive (both surfaces retained, §3.1). Reconcile the `BiSSE.null`/`HiSSE.null` default quirks. Gate: `test_sse.py` many-σ oracles + reproducibility unchanged; flag audit for `--lambda0`/`--mu0`/`--diffusion`/`--root-value`/`--genome-size`.

**P3 — CLI graph engine + unified null + benchmark split.** Replace `_run_coevolve_mode`'s `if eset==...` waterfall + the 14 tree-input/output `parser.error` guards + the ~10 null-legality strings with the table-driven dispatch (§4.2) and the one null layer (§4.4). Make `cid` the uniform hidden-driver transform; fold the three `_run_*_cid_null` runners into one **named benchmark workflow** (not the `cid` knob), preserving `*_ground_truth.tsv`/`null_manifest.tsv`. Keep `--couple` (dual parsing), `--null`/`--hidden`, the default BiSSE, the per-edge output-file contracts, and the `traits:genomes` `--write` flexibility. Gate: `test_cli.py` coevolve suite + manifests + benchmark-preserved test.

**P4 — Joints as two-edge compositions.** Express T↔G (`trait_gene_feedback` co-integrator), S↔G (`co_diversification`), S↔T (ClaSSE) as explicit two-edge compositions; keep the co-integrator/Gillespie solvers. Gate: joint reproducibility + zero-burst/zero-cladogenesis reductions.

**P5 — Sequence tier.** T→Σ and G→Σ reframes on **core** codon ω / clock (§5.1-5.2, near-zero code: `Clock` subclass, `family_factors`/`_annotate` hooks, ω-per-lineage cache) + the one G↔Σ concerted-evolution showcase (§5.3, entirely on core substrate — core codon ω + core `ConversionModel`, no experimental dependency). A no-event contrast baseline, if needed, is a fixed-ω core codon run with no gene-event multiplier (`family_factors=1`, no `node.kind` hook) — no experimental surface. Gate: neutral sequence runs byte-identical; new Σ edge tests; `Events_trace` precondition fails loudly.

**P6 — Docs / figures / manuscript.** Rewrite guide + manual ch13 in lockstep; redraw `fig_coevolve_modes` → diamond + edit `four_levels.svg`; recaption per-edge figures; update `cli.md`/`index.md`/`comparison.md`; flag the external `ZOMBI2_PAPER` as a separate deliverable. Gate: `make figures` clean; mkdocs builds.

Dependencies: P1 blocks all; P3 depends on P1+P2; P4 depends on P2; P5 depends on P1 (and the Σ modules, independent of P2-P4); P6 last (documents the shipped surface).

---

## 10. Risks, open questions & consistency with prior ratified decisions

### 10.1 Top risks

1. **Response-form heterogeneity is the crux.** Scalar exp-link (`traits:genomes`, `genes:species`), free per-state Table (SSE), free Curve (QuaSSE), free 2-entry value Table (`genes:traits`), JUMP kernel (event edges). The one `Response` spec must express all five or silently narrow expressiveness (e.g. forcing SSE into exp-link loses free per-state rates). Build Table/Curve first-class, not as exp-link special cases.
2. **RNG draw order — scoped by the identity policy (§6).** The feedback co-integrator (OU step then n uniforms per step), `gene_diversification` (lineage-then-event-then-target weighting), SSE thinning, and `SequenceEvolution`'s `sorted(_natkey)` family draws each have a specific stream. For the **fused/co-integrator** paths and the **natural reductions**, the stream must be replicated exactly (byte-identity). For the **layered rate-path reframes**, the base-class swap reorders emissions and a fixed-seed run diverges — this is **accepted** (statistical identity, §6), not a bug; it is gated by run-twice reproducibility + the formula/total-rate/`null==base` checks, not a before/after diff.
3. **Driver-dynamics vs coupling-params entanglement.** Every current constructor mixes the driver's own dynamics, the target's intrinsics, and the coupling coefficient in one flat signature (e.g. `GeneDiversification` mixes `loss/origination/transfer` + `driver_speciation/extinction` + `cladogenetic_loss/gain`). The clean grammar split risks breaking byte-identity if wiring order (burst-then-anagenesis, trait-first-then-genes) is not preserved exactly.
4. **Null legality flattening.** Auto-exposing `null(kind='cid')` uniformly resurrects the very `TypeError('workflow not a transform')` guards. Encode the driver-archetype→allowed-kinds matrix; do not flatten. Keep the cid *transform* and the ground-truth *benchmark* distinct (§4.4).
5. **Cross-edge flag-reuse.** Renaming a flag under one edge silently changes another (asymmetric reuse of `--lambda0`/`--mu0`, `--diffusion`/`--root-value`, `--genome-size`, `--theta-*`/`--trait-*`, `--panel`/`--effect-loss`). Every rename needs the §6 audit.
6. **Surgical vs blanket renames.** `genes`→`genomes` must not touch Gene-stem unit classes, the `Genomes` result class, or the `Genomes` ancestral-genomes output dir.
7. **Recipient-seam couplings — signature-adaptation, not seam-building.** The recipient seam is BUILT: `PairModifier` (`zombi2/genomes/transfers.py:92`, attached via `TransferModel.pair=`, exact `pairs` + clade→clade `blocks`/highways) already exists. A driver→recipient coupling is now a matter of keying into `TransferModel.pair`, not building the seam. The only residual risk is that `PairModifier.factor(donor, recipient)` has a different signature from the emission `Modifier.factor(event, family, branch, time)`, so a recipient-side coupling must adapt to `(donor, recipient)` rather than reuse the emission form.
8. **External manuscript** drifts silently (off-worktree).

### 10.2 Open questions (for the author)

- **`BiSSE.null` default = `cid`** while all others = `neutral`, and `HiSSE.null` allows only `neutral`. Preserve these quirks or normalise (documented behaviour change)?
- **Benchmark surface name (minor).** The cid-transform-vs-benchmark split is DECIDED (§4.4): `cid` is the uniform hidden-driver transform, and the ground-truth-withholding benchmark moves to its own named surface. The only thing left open is the final name for that surface (provisional: `--benchmark neutral-channel`) — a naming detail, not a design question.
- **`--all` joint run** (currently "Not yet implemented"): **STAYS DEFERRED (author decision)** — out of scope for this migration. The graph engine may make it fall out later, but it is not a deliverable here.
- **ω-class discretisation** granularity for per-lineage codon models — a tunable, or a fixed default? **Non-blocking implementation detail** (does not gate any PR).
- **G↔Σ showcase placement** — **default-shipped core edge** (notebooks are not shipped — author decision 2026-07-18). It is pure-core; there is no opt-in selection extra to gate the showcase behind.

### 10.3 Consistency with ratified decisions (honor exactly)

- **D1/C6** — node token `genomes` everywhere; `genes` = warned alias. Grammar tokens read `traits:genomes`/`species:genomes`/`genomes:traits`/`genomes:species`. **Gap remains:** `cli.py` `_COEVOLVE_EDGES` still holds the `genes:*` tokens internally and normalises input `genomes`→`genes` (`cli.py:~1448-1454`) — `genes` is still canonical *internally*, so the genomes-everywhere end-state is not yet reached in the CLI edge set. This is the one remaining P0 item (§9).
- **C6/C8** — do NOT genome-ify Gene-stem unit classes; use `TraitGene*` (not `TraitLinked*`); SSE names left alone, documented as `traits:species`; the guide's edge→class table is canonical.
- **D4/C5** — command stays the verb `coevolve`; coevolve levels are {species, genomes, traits} (Σ is downstream, does not drive into S).
- **P1 vocabulary** — `level` = the four domains only; `model` = a process; `lineage` = the tree-entity (never "branch"/"per-genome"); `rate` = per-time; `modifier` = dimensionless multiplier.
- **rate-vocabulary.md** — use `lineage`, `PerLineageRates`/`LineageRates`/`LineageModifier`; per-lineage bursts, not per-branch. A driver's influence is a **modifier** (multiplier), distinct from the opportunity count (`rate = base × opportunities × modifiers`).
- **opportunity-knob.md (opportunity, Part 3)** — the three axes are `base × opportunity × modifiers`; a coupling is a **modifier**, orthogonal to the opportunity (`per=`) knob (`site ⊂ copy ⊂ lineage`, or `shared`). Keep them separable — a coupling never sets `per=`. Track the incoming rename direction: `PerCopyRates`/`PerLineageRates`/`SharedBirthDeath` → `Rates(per=…)`/`BirthDeath(per=…)` presets; target-base references in this doc should move to the `per=` form once the code lands (classes stay as deprecated aliases). The grammar's declarative `Coupling`/`Response` layer and Part 3's declarative rate interface (step 4) are the **same** refactor — build once.
- **rate-modifiers.md** — express a driver→target coupling as a `Modifier` on the target's rate (emission or recipient seam); cite this doc, don't invent a parallel mechanism; respect the byte-identity/RNG-order constraint.
- **D2/C7** — lowercase_snake for every generated table/tree/dir; this design doc is `docs/design/coevolve-grammar.md` (lowercase-kebab).
- **D3/C9** — the experimental ESM `selection` family is PURGED (done); "selection" means core codon ω (`codon_models.py`) everywhere; `experimental ils` is the lone remaining experimental entry.
- **P2 deprecation window** — new aliases (`genes`, `TraitLinked*`, `--sse-model`/`--trait-model`, renamed flags/kwargs) land 0.3.0, removed 0.4.0, via the existing C1 `__getattr__` map + `_deprecated_flag` helper — no open-ended parallel alias set.