> ⚠️ **Largely superseded by [SPEC.md](SPEC.md).** Parts of this document are stale; SPEC.md is authoritative for all conventions, vocabulary, and the coupling model. Read SPEC.md first.

# ZOMBI2 model architecture: a shared model-plugin pattern

**Status:** Draft / proposal — for review
**Date:** 2026-07-13
**Author:** Claude (with Adrián A. Davín)
**Scope:** Internal reorganization of how the four simulation levels — **species**, **genomes**, **traits**, **sequences** — declare, select, and dispatch their model families. No user-facing Python API or CLI flag changes.

> This is a *design/dev* document, deliberately not wired into the published mkdocs site. It records a plan to be reviewed before any code moves. Every `file:line` reference below was verified against the current `main` in a code audit on 2026-07-13.

---

## 1. TL;DR

ZOMBI2 has grown a large, healthy family of models at every level. The problem is **not** that there are too many models — a genome-evolution simulator *should* have many. The problem is that the four levels each independently invented a *different* convention for "a family of interchangeable models + a driver that runs them," ranging from genuinely clean (sequences) to entirely ad hoc (species). Knowledge doesn't transfer between levels, and the weakest level (species) pays a real, error-prone dispatch tax on every new model.

This document proposes:

1. **One tiny shared registry primitive** (`zombi2/_registry.py`, ~150 lines) + **per-level capability descriptors**, applied uniformly across all four levels. It replaces scattered `isinstance` / `if-elif` / `getattr`-with-default dispatch. It is *additive* and *internal*: no shipped public API or CLI flag changes.
2. **A `SpeciesModel` capability contract** that collapses the ~30 species-side dispatch sites into declarative data and turns "silently routed to the wrong engine" into a loud error. This is the priority fix — it's the worst offender *and* it blocks the roadmap of adding more diversification models.
3. **Unification of the four forward growth loops** (`_grow`, `_grow_gillespie`, `_simulate_sse`, `_simulate_quasse`) into one parameterized engine — the largest code-reuse win, done carefully behind the existing analytic tests.
4. **A Python↔Rust event-parity test** that makes the second (Rust) engine's opcode table unable to drift silently from the Python `EventType` enum.

**A note on scope — this is not one commitment.** An adversarial review of this doc rightly flagged that the heavier pieces risk the same over-abstraction that sank the coupling subsystem. So the plan is split into two tiers (see §15):

- **Recommended lean slice** (value clearly exceeds cost): (1) the Rust parity **+ default-deny eligibility** test (§14); (2) merging forward loops 1+2 by extracting shared helpers (§12, Step 1); (3) a **per-class species capability attribute** that kills the silent engine misroute (§8). These need no shared framework.
- **Optional extensions** (adopt only if §17's open questions resolve in favor): the shared cross-level `_registry.py`, `build_from_args` foreign-flag rejection, the auto-generated capability matrix, and the full 4-loop `Sampler`/`View`/`Action` engine.

The capability matrix (Table 1 auto-generation) is a *nice-to-have*, not a driver — and as written it produces one table **per level**, not the single combined Table 1 (see §5.2, §17 Q2). Don't let a paper table dictate that every model carries a `summary`-string descriptor.

---

## 2. The diagnosis: inconsistency, not mess

Every level is fundamentally the same shape: *a set of interchangeable models + a driver that runs a model on the tree + a way to pick a model by name from the CLI*. Each level solved this differently:

| Level | How it organizes its models | Grade |
|---|---|---|
| **Sequences** (`zombi2/sequences/models.py`) | Real registry: `_NT_MODELS`/`_AA_MODELS` dicts merged into `_MODELS`, a single `make_model(name, …)` front door, one normalization kernel `_reversible_model`, an `is_protein_model` predicate | **Clean — the reference** |
| **Traits** (`zombi2/traits/models.py`) | Duck-typed protocol: `kind` + `root_value` + `evolve`/`evolve_branch`; one driver `simulate_traits` | **Clean engine**, but CLI is a flat if/elif and ~half the models are unexposed |
| **Genomes** (`zombi2/genomes/…`) | Partial registry: `EventType` enum + `supported_events()` capability set ∩ rate-model `event_weights` — but `_fire` still hardcodes some events, and a **second engine in Rust** duplicates the event dispatch | **Mixed** |
| **Species** (`zombi2/species/…`) | *Nothing.* Standalone classes + scattered `isinstance`/`getattr`; **four** parallel forward growth loops (two here, two in `coevolve/sse.py` — see §9) | **Ad hoc — the offender** |

The missing thing is not code sharing across levels (those processes are genuinely different — see §9). It is a **shared *pattern***: a contributor who learns one level can't navigate the others, and the species level, having no abstraction at all, makes every addition expensive and fragile.

### The two "taxes" this creates

- **Tax #1 — the species dispatch tax (a fixable smell).** Adding one species model means editing ~6 scattered `isinstance`/`getattr` sites with nothing forcing you to find them all. Miss `forward.py:439` and a constant-between-events model silently routes into the *thinning* loop instead of the *Gillespie* loop — a plausible-looking but wrong tree, no exception.
- **Tax #2 — the Rust parity tax (a permanent structural cost).** The default genome model runs only on the compiled Rust engine (no Python fallback). Every genome-level event is implemented twice, in two languages, with the opcode↔`EventType` mapping hardcoded in **three** independent places bound only by index convention. Implement a new event in Python only, and it's a silent no-op on the fast path most users actually run.

Tax #1 is worth paying down once (pure refactor). Tax #2 can't be removed (it's the price of the ~581× speedup) — but it can be made *undriftable* with a test.

---

## 3. Goals and non-goals

**Goals**

- Replace scattered type-switch dispatch with **declarative, co-located** model metadata.
- Make model **capabilities** explicit (backward-capable? which growth engine? forward-only features?) so illegal combinations fail loudly and early.
- Make it cheap and safe to **add a model** at any level: write the class, declare its capabilities, done — CLI exposure and validation follow automatically. *Caveat:* this holds fully only for a model that reuses an existing simulation path; one needing a genuinely new event menu (e.g. protracted speciation's incipient→good transition, which is neither birth nor death) also needs engine work (see §10 / Part C).
- Reduce the four forward growth loops to shared machinery (fully or partially — see §10).

**Non-goals (explicit anti-scope)**

- **No public API change.** `z.BirthDeath(...)`, `z.make_model(...)`, `z.simulate_species_tree(...)`, `z.simulate_traits(...)` are untouched. The registry is a new *internal* seam. (v0.2.0 is on PyPI; backward compatibility is required.)
- **No CLI flag change.** Every existing flag keeps its name, default, and meaning.
- **No grand unified cross-level engine.** Traits (continuous diffusion), genomes (discrete events), and sequences (per-site CTMC) are different processes; forcing them through one superclass is exactly the coupling that was purged once already.
- **No `_fire` *refactor*.** Its `TRANSFER`/`CONVERSION` branches carry cross-genome orchestration (recipient selection, establishment probability, paired log records) that is not pure per-genome state — don't collapse them into a per-genome handler protocol. **But this is *not* a claim that new events never touch `_fire`:** a single-target event (a new rearrangement reducing to one `genome.apply(event, selection)`) needs no `_fire` change, whereas a whole-genome/orchestrated event (WGD, fission/fusion) *does* require its own `_fire` branch — `_fire`'s generic tail draws exactly one `selection` (`genome_sim.py:414-420`), so routing WGD through it would duplicate a single gene, not the genome. The `EventType` docstring's blanket "no call site changes" (`events.py:24-27`) is only true for the single-target subset and should be corrected.
- **No metaclasses, plugin auto-discovery, or import-time magic.** The whole point is *lean*.

---

## 4. Design principles

1. **Minimum viable structure.** The shared core is a handful of dataclasses and free functions. If a construct isn't needed by code that exists today, it isn't added.
2. **Structural typing over inheritance.** `typing.Protocol` (`runtime_checkable` only for optional asserts) matches the existing duck-typing; no model class is forced to inherit anything.
3. **Additive, behavior-preserving migration.** Every step keeps the suite green. Old `parser.error` guards can stay as belt-and-suspenders until the registry demonstrably subsumes them.
4. **Capabilities are per-level data.** The shared core is capability-agnostic (`caps: Any`); each level defines its own small `Capabilities` dataclass. The core never learns level-specific detail — the failure mode that sank the coupling subsystem.

---

## Part A — The shared model-plugin pattern

## 5. What a registry entry carries

The dispatch audit shows the current code answers **three distinct questions in three different places**. The registry co-locates them per model:

1. **"Build the object."** A `build(**params) -> model` callable — usually the class itself; a small closure when name→class is not 1:1 (see §7).
2. **"What parameters does *this* model take, and how do they appear on the CLI?"** A tuple of `Param` specs. This is the hardest current gap: today every model's flags are flat, always-present, with hardcoded defaults, and each builder both hand-picks its own subset *and* asserts the others are unset via `parser.error` guards (`cli.py:450-455, 475-479`; note `:461` is a *required*-param check, a different shape).
3. **"What can this model do?"** A `Capabilities` descriptor (a per-level frozen dataclass) that replaces the `isinstance` engine split, the `supports_backward`/`supports_ghosts` probes, and the forward-only-feature rejections.

### 5.1 The core (`zombi2/_registry.py`)

```python
@dataclass(frozen=True)
class Param:
    """One tunable parameter and how it appears on the CLI. `flags`/`kwargs`
    forward verbatim to argparse.add_argument, so a model reuses the EXACT
    current flag spec. Shared params (e.g. --sampling-fraction) are declared
    with the same dest by several entries and de-duplicated on add."""
    dest: str
    flags: tuple[str, ...]
    kwargs: Mapping[str, Any] = field(default_factory=dict)
    group: str = "model parameters"
    def add_to(self, container): container.add_argument(*self.flags, dest=self.dest, **dict(self.kwargs))
    @property
    def default(self): return self.kwargs.get("default")

@dataclass(frozen=True)
class ModelEntry:
    name: str
    build: Callable[..., Any]        # (**{param.dest: value}) -> model object
    params: tuple[Param, ...] = ()
    caps: Any = None                 # a per-level Capabilities dataclass
    help: str = ""
    def read_params(self, args): return {p.dest: getattr(args, p.dest) for p in self.params}

class Registry:
    def __init__(self, level): self.level, self._by_name = level, {}
    def register(self, entry):
        if entry.name in self._by_name: raise ValueError(f"{self.level}: duplicate {entry.name!r}")
        self._by_name[entry.name] = entry; return entry
    def get(self, name):
        try: return self._by_name[name]
        except KeyError:
            raise ValueError(f"unknown {self.level} model {name!r} "
                             f"(choose from {sorted(self._by_name)})") from None
    def entries(self): return list(self._by_name.values())
    def names(self): return tuple(self._by_name)
```

Two free helpers replace the hand-written argparse blocks and the per-builder guards:

```python
def add_model_args(parser, registry):
    """Add every entry's params, de-duplicated by dest, grouped as declared.
    Replaces _add_species_args / _add_trait_args' monolithic flag lists."""
    seen, groups = {}, {}
    for entry in registry.entries():
        for p in entry.params:
            if p.dest in seen: continue
            seen[p.dest] = p
            groups.setdefault(p.group, parser.add_argument_group(p.group))
            p.add_to(groups[p.group])

def build_from_args(registry, name, args, *, on_foreign="warn"):
    """Resolve name -> entry, reject params belonging to OTHER models that were
    changed from their default, then build with just this model's params.
    Generic replacement for the scattered parser.error guards AND for
    make_model's silent param-swallowing."""
    entry = registry.get(name)
    mine = {p.dest for p in entry.params}
    for other in registry.entries():
        if other.name == name: continue
        for p in other.params:
            if p.dest in mine: continue
            if getattr(args, p.dest, p.default) != p.default:
                msg = f"{p.flags[0]} does not apply to {registry.level} model {name!r}"
                if on_foreign == "error": raise ValueError(msg)
                if on_foreign == "warn":  warnings.warn(msg, stacklevel=2)
    return entry.build(**entry.read_params(args))
```

> **`build_from_args` is an *optional* extension, not part of the lean slice.** The existing per-builder `parser.error` guards produce *better, teachable* messages (e.g. "`--fossilization`/`--removal` require `--mode forward` — the backward reconstructed sampler assumes complete sampling") than the generic "X does not apply to Y", and the foreign-flag check has a limitation the risk table concedes: it treats "set" as `value != default`, so passing the default reads as unset. **Keep the specific guards**; adopt `build_from_args` only where hand-written guards become unmanageable. The `Param` / `add_model_args` half (per-model flag *declaration*) is the more clearly worthwhile part.

### 5.2 The capability matrix (a nice-to-have, per level)

```python
def capability_matrix(registry, *, fields=None):
    """Reflect over each entry's Capabilities dataclass -> one row per model.
    Rendered to Markdown, this is a per-level capability table."""
    rows = []
    for e in registry.entries():
        row = {"model": e.name}
        if dataclasses.is_dataclass(e.caps):
            for f in dataclasses.fields(e.caps):
                if fields is None or f.name in fields:
                    row[f.name] = getattr(e.caps, f.name)
        rows.append(row)
    return rows
```

**Why this is a nice-to-have, not a driver:** `capability_matrix` takes *one* registry and reflects over *one* level's caps fields, so the four levels yield four tables with **disjoint columns** — not the single combined Table 1 the paper shows (see §17 Q2). Producing one table needs either a single cross-level `Capabilities` god-struct or an explicit projection of each level's fields onto a shared column vocabulary; and `fields=None` dumps non-tabular fields (`summary` prose, tuples), so a real render needs a per-level column whitelist. Worth doing *only* if Table 1 actually drifts in practice — don't let it force every model to carry a `summary` descriptor.

---

## 6. Reconciling with the sequences registry (the reference — and its wart)

`zombi2/sequences/models.py` is already the closest thing to the target and is the model to generalize. The **reusable** parts become the generic `Registry`:

- merged dict `_MODELS = {**_NT_MODELS, **_AA_MODELS}` (`models.py:194`);
- a single public constructor `make_model(name, *, kappa, freqs, rates)` (`models.py:207`);
- a routing predicate `is_protein_model` (`models.py:202`);
- one normalization kernel `_reversible_model` (`models.py:88`) every model funnels through.

But it also exposes the exact **wart** this proposal must fix generically. `make_model` dispatches the 4 DNA models through a hand-written `if name == …` ladder (`models.py:214-221`) — *not* the dict — because **each DNA model consumes a different parameter set** (`kappa`; `kappa+freqs`; `rates+freqs`). So the DNA registry entries are decorative for dispatch, and `make_model` **silently ignores** inapplicable params (`--kappa` on `lg`, `--base-freqs` on `k80`).

The fix is exactly the recommendation the audit surfaced: **separate "which class" (the dict) from "with what params" (a per-model `Param` spec)**, instead of baking a fixed keyword vocabulary into one constructor signature. `make_model` becomes a thin wrapper over the registry — same public signature, but now it *knows* which params each model consumes and can warn on misconfiguration:

```python
def make_model(name, *, kappa=2.0, freqs=None, rates=None):     # SAME signature
    entry = SEQ.get(name.lower())
    supplied = {"kappa": kappa != 2.0, "freqs": freqs is not None, "rates": rates is not None}
    extra = [k for k, on in supplied.items() if on and k not in entry.caps.consumes]
    if extra: warnings.warn(f"{name} ignores {extra}")          # was: silent
    return entry.build(kappa=kappa, freqs=freqs, rates=rates)

def is_protein_model(name):                                     # SAME behavior, single-sourced
    return SEQ.get(name.lower()).caps.alphabet != BASES
```

> Note the single-character-alphabet assumption (`SubstitutionModel.alphabet` is a string indexed one char per state) is orthogonal to the registry and is *not* addressed here — it's the thing to fix when/if codon models land, and belongs in its own note.

---

## Part B — Species: the capability contract (priority fix)

## 7. The species dispatch census

The species layer branches on 6 model classes (`BirthDeath`/`Yule`, `EpisodicBirthDeath`, `ClaDS`, `DiversityDependent`, `CladeShiftBirthDeath`) across ~30 sites. The load-bearing ones:

| Site | What it decides |
|---|---|
| `species/sim.py:94` | `isinstance(model, (ClaDS, DiversityDependent, CladeShiftBirthDeath))` → reject forward-only from backward |
| `species/sim.py:102-105` | forward-only feature rejections (`fossilization`, `removal`, `mass_extinctions`) — via **attribute-type** check + `getattr` |
| `species/sim.py:110` | `isinstance(model, BirthDeath) and sampling_fraction < 1` → reject incomplete sampling in constant-rate backward |
| `species/forward.py:57-81` | `_ForwardRates.__init__`: `if Episodic / elif BirthDeath / else NotImplementedError` |
| `species/forward.py:439` | `heterogeneous = isinstance(...)` → **thinning vs Gillespie engine** |
| `species/forward.py:440, 446, 452` | reject `n_tips` mode for Episodic / CladeShift; catch-all `NotImplementedError` for unknown types |
| `species/forward.py:467-473` | reject `n_tips` mode when `sampling_fraction < 1` |
| `species/forward.py:474, 479, 484` | per-model n_tips/age guards (`DiversityDependent.K`, `CladeShift.clade_shifts`, `mass_extinctions`) |
| `species/forward.py:500-510` | build `_ClaDSView`/`_DDView`/`_ShiftView` by `isinstance` |
| `species/ghosts.py:38, 57, 65` | episodic vs constant ghost rates; `else NotImplementedError` |
| `cli.py:410-483` | multi-key model selection (`--diversification`, `--clade-shift`, mode) + the **arity heuristic** |

**Two dispatches are NOT type-switches** — and are the trap for a naive name→class registry:

1. **Arity heuristic** (`cli.py:428`): `episodic = args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1`. `BirthDeath` vs `EpisodicBirthDeath` is inferred from *how many values* `--birth/--death` got, not an explicit choice.
2. **Attribute-type check** (`sim.py:103`): `sum(foss) if isinstance(foss, list) else foss` — branches on whether `fossilization` is a scalar (`BirthDeath`) or a list (`EpisodicBirthDeath`).

Both are preserved for *construction* (see §8) — with one subtlety the review caught: collapsing `BirthDeath` and `EpisodicBirthDeath` behind a single capability descriptor is **not** enough, because the two have **different downstream capabilities** (backward `BirthDeath` rejects ρ<1 at `sim.py:110` while Episodic *allows* it; Episodic forbids `n_tips` mode at `forward.py:440` while `BirthDeath` allows it). §8 handles this by giving them **separate entries / caps**, with the arity heuristic living in `resolve_species_name`.

## 8. The species capability contract

**Recommended (lean) form: a capability descriptor on each model class — no shared framework.** The bug to kill is the silent misroute at `forward.py:439`, where a forgotten `isinstance` branch routes a model into the wrong engine and yields a plausible-but-wrong tree. The minimal fix that forces every *new* model to answer the routing questions **at definition time, co-located with the class**, is a small descriptor on the class itself:

```python
class GrowthEngine(enum.Enum):
    THINNING = "thinning"; GILLESPIE = "gillespie"

@dataclass(frozen=True)
class SpeciesCaps:
    growth: GrowthEngine                        # replaces the forward.py:439 isinstance
    supports_backward: bool = False             # replaces sim.py:94
    supports_ghosts: bool = False               # replaces ghosts.py:65
    supports_n_tips: bool = True                # Episodic/CladeShift = False (forward.py:440,446)
    incomplete_sampling_backward: bool = False  # Episodic = True, BirthDeath = False (sim.py:110)
    make_view: Callable | None = None           # replaces the forward.py:500-510 view ladder
    forward_only_features: tuple["Feature", ...] = ()   # replaces sim.py:102-105

class ClaDS:                                     # each class declares its OWN caps (locality!)
    caps = SpeciesCaps(GrowthEngine.GILLESPIE,
                       make_view=lambda m, present: ClaDSView(m, present),
                       forward_only_features=(MASS_EXTINCTIONS,))
```

`species_caps(model)` walks `type(model).__mro__` (so `Yule(BirthDeath)` inherits). The dispatch sites collapse to reads, and an *unregistered* `growth` is a **loud error**, not a wrong tree:

```python
# forward.py:439   view, grow = engine_for(species_caps(model).growth)   # raises on any unknown GrowthEngine
# forward.py:500   view = species_caps(model).make_view(model, present)  # replaces the isinstance view ladder
# sim.py:94        if not species_caps(model).supports_backward: raise ValueError(...)
# ghosts.py:65     if not species_caps(model).supports_ghosts: raise NotImplementedError(...)
```

**Per-feature `is_active`, not a global `_NEUTRAL` + branchy function.** The instance-level "is fossilization *active*?" check must not re-scatter the `sim.py:103` scalar-vs-list special-case (an earlier draft's `_NEUTRAL` dict + `if feat == …` branches did exactly that — three edits to add one feature). Co-locate it with each feature so a new forward-only feature is *one* declaration:

```python
@dataclass(frozen=True)
class Feature:
    name: str
    is_active: Callable[[Any], bool]         # knows this feature's neutral value AND its shape

FOSSILIZATION    = Feature("fossilization",    lambda v: (sum(v) if isinstance(v, list) else v) > 0)  # scalar OR Episodic list
MASS_EXTINCTIONS = Feature("mass_extinctions", lambda v: bool(v))
REMOVAL          = Feature("removal",          lambda v: v != 1.0)

def active_forward_features(model, caps):
    return [f.name for f in caps.forward_only_features
            if f.is_active(getattr(model, f.name, None))]
```

> `sampling_fraction` is deliberately **not** a `Feature`: it is the separate per-class check at `sim.py:110` (constant-rate `BirthDeath` only — hence the `incomplete_sampling_backward` cap above), and Episodic even stores the extant fraction under a *different* attribute (`rho`, `model.py:186`). Folding it into `active_forward_features` would misfire.

**`BirthDeath` and `EpisodicBirthDeath` need *two* entries, not one.** They share the CLI arity heuristic but differ in the downstream capabilities just listed (ρ<1 backward; `n_tips` mode), which a single shared caps cannot express. So each is its own key with its own `SpeciesCaps`, and the arity heuristic lives in the CLI name resolver that picks *which* to build:

```python
def resolve_species_name(args):
    if args.clade_shift: return "clade-shift"
    if args.diversification != "constant": return args.diversification    # clads | diversity-dependent
    return "episodic" if (args.shifts is not None or len(args.birth) > 1
                          or len(args.death) > 1) else "constant"           # the arity heuristic
```

The mode check then reads the selected model's caps — replacing the per-builder `parser.error` guards at `cli.py:407-472`:

```python
def _build_species_model(args, parser):
    model = build_species(resolve_species_name(args), args)   # per-class build; or the Part-A registry
    caps = species_caps(model)
    if args.model == "backward":
        if not caps.supports_backward:
            parser.error(f"{type(model).__name__} is forward-only; add --mode forward")
        if (bad := active_forward_features(model, caps)):
            parser.error(f"{', '.join(bad)} are forward-only; use --mode forward")
        if not caps.incomplete_sampling_backward and getattr(model, "sampling_fraction", 1.0) < 1.0:
            parser.error("incomplete sampling (ρ<1) is forward-only for constant-rate BD; use Episodic or --mode forward")
    return model
```

**The heavier (optional) form** wraps all of this in the shared `Registry`/`ModelEntry` with `Param` specs and `build_from_args` (Part A). That buys uniform per-model CLI-flag declaration and a capability-matrix row — but it holds the caps on the *entry* rather than the class, which is *worse* locality than the lean form. Adopt it only if the `Param`/`add_model_args` machinery is independently wanted.

**Net (lean form):** ~30 dispatch sites collapse to `species_caps(model)` + `resolve_species_name` + `active_forward_features`, with **no shared framework**. Adding a species model that reuses an existing engine becomes: write the class, add its `caps` descriptor, done — a missing or unknown capability is a loud error, not a silently wrong tree. (A model needing a genuinely *new* event menu still needs Part C — see §10.)

---

## Part C — Forward-engine unification

## 9. The four loops share one skeleton

There are **four** forward growth loops duplicating the same lineage-lifecycle logic:

| Loop | File | Sampler | Models | Extra menu |
|---|---|---|---|---|
| `_grow` | `species/forward.py:88` | thinning (bound) | BirthDeath/Yule, Episodic | mass-ext, **serial-sampling ψ + sampled ancestors** |
| `_grow_gillespie` | `species/forward.py:300` | **exact Gillespie** | ClaDS, DiversityDependent, CladeShift | mass-ext, **clade shifts** |
| `_simulate_sse` | `coevolve/sse.py:481` | thinning (bound incl. anagenetic out-rate) | BiSSE/MuSSE/HiSSE/CID | **anagenetic Q transition** (lineage continues), cladogenesis |
| `_simulate_quasse` | `coevolve/sse.py:570` | thinning (bound over reachable x) | QuaSSE | **continuous diffusion**, cladogenesis |

All four share the identical skeleton: crown of 2 at `t=0` → `while` with identical guards (`n==0` reject / `n_tips` reached break / `max_lineages` raise) → draw exponential waiting time → boundary/scheduled handling → pick a lineage → resolve an event → speciation appends 2 daughters & swap-pops parent, extinction marks dead → finalize survivors at `end` under ρ, condition on ≥2 → wrapped by a `max_attempts` retry that builds `Tree` + `_name`.

They differ in **12 enumerated points**, the load-bearing ones being:

- **Waiting-time / acceptance:** thinning (`dt=Exp(1/(n·bound))`, accept w.p. `total/bound`) vs exact Gillespie (`dt=Exp(1/total_rate)`, always accept). **The deepest fork.**
- **Lineage selection:** thinning picks *uniformly*; Gillespie picks *rate-weighted*. **Selection-mode is rigidly coupled to sampler-mode** — mixing them silently biases which lineage branches.
- **Event menu:** `{speciate, extinct}` + optionally `{serial-sample (ψ), anagenetic transition (continue), continuous diffusion}`.
- **Daughter-state rule:** none / `view.split()` (ClaDS multiplicative jump) / `born()` cladogenetic kick. ClaDS's `split` and SSE's `born()` are the *same* abstraction via different mechanisms.
- **Present sampling ρ:** loops 1&2 use `_finalize_present` (ρ<1 → unsampled ghosts); loops 3&4 hardcode ρ=1.

## 10. Proposed unified engine: two tiers behind one skeleton

**Recommendation: two tiers, not one flat god-function.** Split responsibilities on the two axes that actually vary, so the shared core never grows a model-specific branch:

- **`Sampler`** owns *when the next event happens and which lineage it hits* — the deepest fork. Thinning and exact are **rigidly coupled** knobs (thinning ⟺ uniform pick + accept-by-`total_i/bound`; exact ⟺ rate-weighted pick + always-fire), so they live together as two strategy *objects*, `ThinningSampler(bound)` and `ExactSampler()`, behind one method `sampler.step(...) -> Outcome`. Making the fork a strategy object rather than an `if mode == …` inside the hot loop is what prevents the silent uniform-under-exact / weighted-under-thinning bias.
- **`View`** owns *the rate source, per-lineage state, the event menu, daughter-state rules, scheduled tree-wide events, present-sampling ρ, and any output payload*. One per model family.
- **Core** — `forward_grow(view, sampler, age, n_tips, rng, max_lineages)` — owns *only* the shared skeleton: init/guards/finalize/retry, schedule preemption + memoryless redraw, boundaries, and applying whatever `Action` the View returns. No model-specific branch.

**The event vocabulary.** Every event any loop can fire is one of four `Action`s the View returns from `resolve_event` — this is what unifies the disparate menus:

| Action | Absorbs | Note |
|---|---|---|
| `Speciate(state_a, state_b)` | plain BD `(None, None)`, ClaDS `split` (two jumps), DD/Shift `(s, s)`, ClaSSE `born(i)` | daughter states drawn **twice, independently** |
| `Extinct` | all loops | |
| `Sample(continue_lineage)` | serial ψ | the **only** Action making degree-2 nodes; only `_ForwardRates` returns it |
| `Anagenesis(new_state)` | discrete-SSE Q transition | lineage **continues**, records a segment; no branch |

Because event *selection* lives inside `resolve_event`, the core never sees the menu — which keeps ρ/removal/ψ out of the SSE views (they simply never construct `Sample`) and the anagenesis stochastic map out of the species views. Output payloads (`node_values`, `history`) are owned by the View and populated through no-op-by-default hooks (`on_anagenesis`, `on_close`, `finalize`); `forward_grow` always returns `(root, end)`, and `simulate_sse` reads `view.node_values`/`view.history` afterward.

`Sampler.step` returns a tagged `Outcome` — `FIRE(idx, rates)` / `THINNED` / `SCHEDULE` / `BREAK(end)` — so the core loop is byte-identical for both samplers; only the object differs. `ThinningSampler` additionally owns the continuous `on_interval(dt)` diffusion call (on **every** drawn dt, including thinned candidates) and the retained QuaSSE bound-violation guard; `ExactSampler` owns the `total_rate<=0` degenerate branch (apply next scheduled event / coast to present in age mode / `RuntimeError` in n_tips mode).

Each existing loop becomes a `(View, Sampler)` pair: BirthDeath/Episodic → `(_ForwardRates, Thinning)`; ClaDS / DiversityDependent / CladeShift → `(…View, Exact)`; BiSSE/MuSSE/HiSSE → `(MuSSEView, Thinning)` with bound incl. the anagenetic out-rate; QuaSSE → `(QuaSSEView, Thinning)` where state is a float and `on_interval` does the diffusion.

### Where unification is *not* worth forcing

Merge **loops 1+2 (species) and loop 3 (discrete SSE) unconditionally** — that collapse is clean and high-value. Bring **loop 4 (QuaSSE)** in *only* via the single `on_interval` hook, behind a **gate**: if in review that hook starts accreting continuous-only conditionals in the core, stop and keep `_simulate_quasse` as a thin **sibling** that still reuses the shared `Lineage`, `finalize`, retry, and `Action` machinery but keeps its own diffusion placement. Don't let the core become a god-function to save ~40 lines. Everything else (ψ, anagenesis, ClaDS split, scheduling) merges without leakage.

## 11. The correctness traps (must-preserve invariants)

The audit flagged the exact places a careless unification breaks — each has (or needs) a guarding test:

- **`born()`/`split` is drawn *independently per daughter* (two draws).** Collapsing to one shared draw halves the realized jump variance and fails `test_classe_continuous_jumps_are_normal_zero_jump_sigma2` (pins variance == `jump_sigma2`). The daughter uses the parent's state *at the speciation instant*; ordering of the `born()` rng draws vs daughter-node creation must be preserved or every `*_reproducible` byte-identity test breaks.
- **Thinning ⇒ uniform pick; exact ⇒ weighted pick — the two knobs move together.** A uniform pick under exact, or weighted under thinning, silently biases tree shape and tip-fractions. Only the *statistical* tests (Yule mean count, λ/μ tip-fraction bias) catch it, and only in aggregate.
- **QuaSSE diffuses on *every drawn dt*, including thinned events**, and on the final partial interval *before* the age-break. A hook that only advanced state on *accepted* events would be wrong. Continuous "state evolves over dt" is fundamentally different from discrete "state-change events" (a competing Poisson event subject to thinning) — the interface must model them separately (`on_step` vs `resolve_event`).
- **The QuaSSE runtime bound-violation guard (`sse.py:620`) has NO test.** A refactor could drop it silently. **Add a test before touching it.**
- **Serial-sampling ψ produces the only degree-2 nodes** (sampled ancestors) and must *not* leak ρ/removal semantics into the SSE loops (which assume ρ=1). Guarded by the FBD + sampled-ancestor test groups.
- **Scheduled events fire when `time ≤ t+dt`, discarding the drawn `dt` and redrawing** — valid only because the exponential is memoryless (true for both samplers). Reordering "check schedule / draw dt / test age boundary" changes correctness *and* reproducibility.
- **ρ default must be 1.0 for SSE views** — a stray ρ<1 leaking into the SSE path silently drops tips from `res.values`.

**Reproducibility tax:** every model has a `*_reproducible` test asserting byte-identical `to_newick` for a fixed seed. Any unification will reorder rng draws and require regenerating these golden files. **Plan validation around the durable oracles** (Yule reductions, λ/μ bias, stationary distribution, `history` tiling) — not the byte-identity tests.

## 12. Test-safe migration order

The guarding tests exist and are extensive (`tests/test_species_forward.py`: Yule / age / n_tips / extinction / episodic / FBD / sampled-ancestor / mass-extinction / ClaDS / DiversityDependent / clade-shift groups; `tests/test_sse.py`: Yule reductions, λ/μ tip-fraction bias, MuSSE stationary, HiSSE, QuaSSE, ClaSSE/cladogenesis incl. the `jump_sigma2`-variance test). Validate on the **durable statistical oracles**, not the byte-identity goldens (those must be regenerated once — any unification reorders rng draws).

- **Step 0 — pin the untested defensive paths first.** Add (a) a QuaSSE test that forces `total > rate_bound` and asserts the `ValueError` (`sse.py:620`, today only `rate_bound<=0` is checked at construction), and (b) an exact-sampler test that drives `total_rate<=0` in n_tips mode and asserts the `RuntimeError` (the DiversityDependent-at-K, μ=0 stall, today unreachable via the public API). Pin both *before* the code moves.
- **Step 1 — extract the core + `Action`/`Sampler` vocabulary, species-only.** Reimplement `_grow` as `_ForwardRates.resolve_event` over `ThinningSampler`, and `_grow_gillespie` as `ClaDSView`/`DDView`/`ShiftView` over `ExactSampler`. Keep `_finalize_present`/`_build_schedule` as-is. **Merges loops 1+2.** Run the full `test_species_forward.py`.
- **Step 2 — bring discrete SSE onto the core.** Add the payload hooks + the `Anagenesis` action; reimplement `_simulate_sse` as `MuSSEView` (covers BiSSE/MuSSE/HiSSE/CID via product states) over `ThinningSampler`, with the crown `born()` kick and the two-independent-draw `born()`. Run the `test_sse.py` discrete group incl. the stochastic-map tiling test.
- **Step 3 — bring QuaSSE on behind `on_interval` (gated).** Add `QuaSSEView` + the retained bound guard inside `ThinningSampler`. **Gate:** if `on_interval` forces continuous-only conditionals into the core, abort the loop-4 merge and keep `_simulate_quasse` as a documented sibling reusing `Lineage`/`finalize`/retry/`Action`. Run the QuaSSE + ClaSSE groups, especially `classe_continuous_jumps_are_normal_zero_jump_sigma2`.
- **Step 4 — unify dispatch + delete duplicates.** Route `simulate_forward` (`forward.py:500-517`) and `simulate_sse` (`sse.py:696-735`) through the single `forward_grow` entry, delete the four old loops, and re-confirm `feeds_gene_sim` / ghost-transfer partners still work off the complete tree.

Two design decisions to settle before Step 1 (flagged by the design pass): whether the core carries a **stateless schedule index** vs. the View popping a mutable `schedule` — the latter mutates across `max_attempts` retries, so the stateless index is safer; and keeping SSE **root-state resolution in the dispatch** (not `seed_crown`) to preserve the per-attempt resolution semantics.

---

## Part D — Rust event parity (Tax #2)

## 13. The three-place hardcoding

The gene-family event kinds are represented as `u8` opcodes in Rust and mapped to the Python `EventType` in **three independent places**, bound only by index convention — nothing checks them at build or import time:

- **Rust opcodes** `EV_O=0, EV_D=1, EV_T=2, EV_L=3, EV_S=4` (`rust/src/lib.rs:734-738`), plus code-indexed tables `EV_CHAR` and `ROLES` (`lib.rs:1126-1133`) used by the Rust file-writer.
- **Python bridge** `_EV = (ORIGINATION, DUPLICATION, TRANSFER, LOSS, SPECIATION)` and `_ROLES` (`zombi2/_rust.py:295-303`) — same index order, re-typed by hand.
- **`EventType` enum** string values (`zombi2/genomes/events.py:30-61`) — the ultimate meaning, but not referenced by Rust.

**Only the integer code crosses the FFI boundary.** `EventType.value` never enters Rust; Rust re-hardcodes the identical chars in `EV_CHAR`. The parity contract `_EV[i].value == EV_CHAR[i] == roles` is convention only. Adding an event touches ~9 sites across both languages and **requires a wheel rebuild** (`zombi2_core` is a non-editable compiled dependency; per the project's build rule, `maturin build --release` + force-reinstall) — the most common parity break in practice is a stale wheel.

There are also **two disjoint opcode spaces**: the nucleotide engine `NEngine::fire` (`lib.rs:1774`) uses an unrelated internal order (0=inv, 1=loss, 2=dup, 3=transfer, 4=transposition, _=orig) that never surfaces to Python. A parity scheme for the gene-family opcodes does not cover it.

## 14. Recommendation: three guards, not codegen

The build is `maturin` + `pyo3` `abi3-py310` with **no `build.rs`**. Codegen (a `build.rs` emitting `$OUT_DIR/opcodes.rs` from a shared spec) is feasible but heavy — it adds the crate's first `build.rs` and moves the source of truth into a data file. Prefer three cheap tests, because the opcode table alone does **not** cover the actual failure mode.

**(1) The real anti-silent-drop: make Rust eligibility *default-deny*.** The dangerous failure isn't the opcode table — it's that `eligible()` (`_rust.py:61-75`) is a hand-maintained *deny*-list. Add a new rate field (say `wgd`) and `SharedRates(wgd=…)` with the default `UnorderedGenome` returns `eligible()==True`, routes to the Rust engine, which never fires WGD → **a silent no-op on the fast path most users run.** The opcode parity test passes (the tables are internally consistent); nothing forces the new rate into the exclusion set. Fix: derive eligibility from Rust's *declared* firing capability:

```python
RUST_FIRES = frozenset({EventType.ORIGINATION, EventType.DUPLICATION,
                        EventType.TRANSFER, EventType.LOSS})     # exposed from Rust, single source
def eligible(rates, genome_factory, sampler):
    if genome_factory is not UnorderedGenome or sampler is not None: return False
    if type(rates) is not SharedRates: return False              # KEEP exact-type identity
    return all(w.event in RUST_FIRES                             # default-DENY: any other event -> Python
               for w in rates.event_weights(...) if w.rate > 0)
```

Now any nonzero rate mapping to an event Rust doesn't fire flips selection to the Python engine automatically. (Keep the `type(rates) is SharedRates` / `genome_factory is UnorderedGenome` *identity* checks — the registry must not soften these to `isinstance`, or the fast path is lost.)

**(2) Opcode parity test — against the enum, not the mirror.** Expose the Rust table via a tiny `#[pyfunction] opcode_chars()` / `opcode_roles()` and assert against the **authoritative `EventType`**, not the hand-typed `_rust._EV` copy:

```python
from zombi2.genomes.events import EventType
assert _rust.opcode_chars() == [EventType(c).value for c in "ODTLS"]   # source of truth
```

Make the "move the char/role table into `events.py` so `_rust.py` imports it" step **required, not optional** — otherwise a change made consistently to Rust *and* `_EV` but wrong vs `EventType` still passes. This test also fails loudly on a stale wheel (the exposed table won't match).

**(3) Cross-engine equivalence test — its own line item, mandatory per event.** Run the same event on both engines and assert identical profiles, so a new event can't ship on one engine only. This is what catches guard (1)'s failure mode end-to-end.

**Scoped out (tracked follow-up):** the nucleotide engine `NEngine::fire` uses a *disjoint* internal opcode order (`lib.rs:1774`) that never reaches Python; none of the three guards cover it. Either add a separate `nengine_opcode_chars()` assertion or explicitly leave it out of scope with a tracked issue — don't let §13's own enumerated third drift mode go silently unmitigated.

**This does not remove Tax #2** — a new genome event is still implemented twice. It removes every *silent* failure mode: default-deny eligibility (no fast-path no-op), enum-sourced parity (no table drift / stale wheel), per-event equivalence (no one-engine-only ship).

---

## 15. Sequencing & priority

Split into the lean slice (do these) and optional extensions (decide via §17).

**Recommended lean slice — high value, low risk, no shared framework:**

1. **Rust guards (§14):** default-deny `eligible()` + enum-sourced opcode parity test + a per-event cross-engine equivalence test. Tiny, and the only thing that makes adding a genome event *safe*. Do first.
2. **Species per-class capability descriptor (§8, lean form):** `SpeciesCaps` on each model class + `resolve_species_name` + `active_forward_features`. Collapses ~30 dispatch sites, kills the silent engine misroute, de-risks the "add more diversification models" roadmap. Bulk of the correctness win, no registry required.
3. **Sequences `make_model` warn-on-inapplicable-param (§6):** a few lines; kills the silent-param-swallowing wart. Uses the existing dict, no new framework.
4. **Forward loops 1+2 merge via shared helpers (§12, Steps 0–1):** first pin the two untested defensive paths (Step 0), then extract the crown-init, guard block, `_finalize_present`, and retry wrapper shared by `_grow`/`_grow_gillespie`.

**Optional extensions — adopt only if §17 resolves in favor:**

5. **The shared `_registry.py` + `Param`/`build_from_args` (§5)** and **retrofitting traits onto it.** Buys uniform CLI-flag declaration and cross-level consistency; costs a mandatory shared import for a "shared pattern, not shared code" goal. Justified only if a second level independently wants the identical machinery. *(Note: the trait CLI is not trivial — `Mk` alone has three sub-constructors and `threshold` parses a list, so this needs per-model build closures, not bare classes.)*
6. **Full forward-engine unification (§10):** the `Sampler`/`View`/`Action` framework absorbing loop 3 (SSE) and, gated, loop 4 (QuaSSE). High reuse but the most delicate change (7 correctness traps, full golden regeneration). Only after 1+2 lands cleanly; keep the QuaSSE gate.
7. **Auto-generated capability matrix (§5.2)** — a nice-to-have; skip unless Table 1 actually drifts.
8. **Genomes registry retrofit — recommend *not* doing** (§17 Q4): `_fire` stays hardcoded, so a genome caps row is decorative. If Table 1 must include genomes, tie the row to an *enforced* check (registered event set == union of `supported_events()` across shipped genome types), not a hand-maintained caps row.

## 16. Risks & mitigations

| Risk | Mitigation |
|---|---|
| `build_from_args` detects "foreign flag" by comparing to the declared default; explicitly passing the default reads as "not set" | Same imperfection the current `!= 1.0` guards already have — not a regression; start as a *warning* for sequences, escalate to error only at a major version |
| Capabilities are per-*type* but some rejections are per-*instance* value checks | Keep them separate by construction: `caps` for type-level questions, `active_forward_features` for instance-level. Call this out in the module docstring |
| `_rust.eligible` uses **exact-type identity** (`type(rates) is SharedRates`, `genome_factory is UnorderedGenome`) | The registry may only *record* the `rust_capable` predicate; it must never mediate the identity check, or the fast path is silently lost |
| Forward-loop unification reorders rng draws | Validate on statistical oracles, regenerate byte-identity golden files deliberately; add the two missing guard tests first |
| A new selection flag isn't threaded into `resolve_species_name` → model unreachable from CLI | CI test asserting every registry entry is reachable via some resolver input |
| Scope creep into "unify all event dispatch" / refactor `_fire` | Explicit non-goal (§3). `_fire`'s cross-genome orchestration stays put |
| "Add a model, done" overstated: a model needing a new event menu (protracted speciation, WGD) needs engine/`_fire` work, not just a caps row | §3 and §8 caveat this; the capability descriptor gates *selection*, not the *simulation path* |
| Rust `eligible()` deny-list silently routes a Python-only genome event to the fast engine (a no-op) | Default-deny eligibility (§14 guard 1) + per-event cross-engine equivalence test |

## 17. Open questions (for Adrián)

1. **Public or private registry?** Expose `z.species.SPECIES` so third parties can register models, or keep `_registry` strictly internal on a v0.2.0 package? *(Lean: private until an external use-case appears.)*
2. **One `Capabilities` struct or per-level?** Per-level (`SpeciesCaps`/`SeqCaps`/…) is leaner and avoids a god-struct. A *single* combined Table 1 would need either that god-struct or a projection layer (§5.2). Since the matrix is now a nice-to-have, the lean answer is **per-level caps + a manually maintained Table 1** — auto-generating one combined table isn't worth a god-struct. Comfortable with that?
3. **Where do entries live?** One `_registry.py` per level package (matches the 4-subpackage structure) or a central `zombi2/models_registry.py`?
4. **Register genomes at all?** Recommendation: **no** (§15 item 8). `_fire` stays hardcoded, so a genome caps row is decorative; the one concrete want — the `rust_capable` predicate — is a standalone function next to `_rust.eligible` (§14). Flagged only because it interacts with whether Table 1 auto-generates.
5. **Forward-engine unification: full or partial?** Commit to one `forward_grow`, or accept two cores (species-side + SSE-side) if unification makes the core leaky?

---

## Appendix — verified reference points

- **Sequences registry:** `models.py:192-224` (registry + `make_model`), `:88-108` (`_reversible_model`), `:202` (`is_protein_model`), `:48-85` (`SubstitutionModel`).
- **Species dispatch:** `sim.py:94, 102-105, 110, 123-126`; `forward.py:57-81, 439, 474-484, 500-510`; `ghosts.py:38-66`; `cli.py:410-483` (incl. arity heuristic `:428`).
- **Forward loops:** `forward.py:88` (`_grow`), `:300` (`_grow_gillespie`), `:167` (`_finalize_present`), `:191-198` (`_build_schedule`), `:289-297` (`_apply_scheduled`); `coevolve/sse.py:481` (`_simulate_sse`), `:570` (`_simulate_quasse`), `:696-711` (dispatch).
- **Rust parity:** `rust/src/lib.rs:734-738` (opcodes), `:1126-1133` (`EV_CHAR`/`ROLES`), `:1233-1252` (consumer match), `:1774` (nucleotide match); `zombi2/_rust.py:295-303` (`_EV`/`_ROLES`), `:61-75` (`eligible`); `zombi2/genomes/events.py:30-61` (`EventType`); `genomes/simulation.py:578-595` (engine selection).
- **Genome seam (the template):** `genome_sim.py:304-309` (`supported_events` ∩ `event_weights`), `:361-421` (`_fire`, left as-is); `genome.py:184, 200, 360, 392`; `rates.py:154-198, 223, 288, 402`.
