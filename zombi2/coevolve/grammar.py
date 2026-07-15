"""The coevolution grammar — one sentence for every coupled model.

Every coupling is the single sentence::

    driver  →  target-variable  :  response

with three slots, each drawn from a closed menu:

* **driver** — whose state pushes. Two kinds: a *state* driver (a value along the tree — a trait
  value, a gene count) or an *event* driver (a set of instants — speciation).
* **target-variable** — the variable of the target level that gets bent, from the fixed
  per-level menu (:data:`TARGET_VARIABLES`). Two kinds: a *rate* (a continuous flow — λ, μ, gene
  loss/gain, OU pull/diffusion, substitution speed, dN/dS) whose effect is a **modulation**; or a
  *state* (the trait value, gene presence, an OU optimum) whose effect is a **jump**.
* **response** — how the driver's value maps to the size of the effect. :class:`Scalar` (the
  default one-knob exp-link), :class:`Table` (per discrete state — recovers MuSSE / a
  ``theta_present``-vs-``theta_absent`` switch), or :class:`Curve` (a nonlinear continuous map —
  recovers QuaSSE). A :class:`Scalar` of strength ``0`` (see :func:`null_response`) is the matched
  **null**: every driver value maps to a multiplier of ``1``.

Two rules read the *graph* of couplings (:class:`CouplingGraph`):

* **solve rule** — a coupling is *directional* (run its driver first, then the target — **layer**)
  unless it sits in a cycle, in which case driver and target co-evolve and must be integrated
  together (**fuse**). Crucially, any arrow *into a substrate* (e.g. ``traits → species``) closes a
  cycle with that substrate's implicit downstream edge (``species → traits``), which is exactly why
  all into-species coupling grows the tree.
* **topology rule** — a coupling may only connect levels **within one tier** of each other. The
  ``species ↔ sequences`` diagonal is **forbidden** — a sequence rides *gene* trees, not the
  species tree.

This module is the **declarative core**: pure Python with no rate-engine dependency. It names the
levels / target-variables / responses, validates a coupling and a coupling graph, and classifies
each coupling as *layer* or *fuse*. Compiling a rate-target coupling down onto a
:class:`zombi2.genomes.rates.Modifier` (the ``CouplingModifier`` bridge) is a separate step, added
once the rate-model refactor (the ``per=`` opportunity knob) settles — see
``docs/design/coevolve-grammar.md`` §4.3.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field

# ═══════════════════════════════════════════════════════════════════════════════
# The diamond: levels, tiers, and the implicit downstream (substrate) edges
# ═══════════════════════════════════════════════════════════════════════════════
#: The four levels, top (substrate) to bottom (sink).
LEVELS: tuple[str, ...] = ("species", "traits", "genomes", "sequences")

#: Dependency tier of each level: species (substrate) → {traits, genomes} → sequences.
_TIER: dict[str, int] = {"species": 0, "traits": 1, "genomes": 1, "sequences": 2}

#: The implicit *downstream* edges — each level rides the tree its parent hands down. Sequences
#: ride GENE trees (genomes), **not** the species tree, so there is no ``species → sequences`` edge;
#: this is what makes the ``species ↔ sequences`` coupling both forbidden (topology) and, were it
#: allowed, a two-tier skip. These edges close the cycles that make into-substrate couplings *fuse*.
_SUBSTRATE_EDGES: frozenset[tuple[str, str]] = frozenset(
    {("species", "traits"), ("species", "genomes"), ("genomes", "sequences")}
)

_RATE, _STATE = "rate", "state"

#: The closed per-level target-variable menu: ``level → {variable → kind}``. A coupling may bend
#: only a variable on this list; the kind (``"rate"`` → modulation, ``"state"`` → jump) is read off
#: here, never supplied by the user.
TARGET_VARIABLES: dict[str, dict[str, str]] = {
    "species":   {"speciation": _RATE, "extinction": _RATE},
    "traits":    {"optimum": _STATE, "pull": _RATE, "diffusion": _RATE, "value": _STATE},
    "genomes":   {"loss": _RATE, "gain": _RATE, "duplication": _RATE, "transfer": _RATE,
                  "presence": _STATE},
    "sequences": {"substitution_speed": _RATE, "selection": _RATE, "residues": _STATE},
}

#: Clamp on the exp-link exponent so an extreme ``strength · driver`` never overflows ``exp`` (or
#: drives a Gillespie loop into an instant hot cycle). Matches ``trait_coupling._MAX_EXPONENT``.
_MAX_EXPONENT = 40.0


def _clamp(x: float) -> float:
    return max(-_MAX_EXPONENT, min(_MAX_EXPONENT, x))


# ═══════════════════════════════════════════════════════════════════════════════
# Response: how a driver value maps to the size of the effect
# ═══════════════════════════════════════════════════════════════════════════════
class Response(ABC):
    """How a driver's value maps to the size of the effect.

    Applied to a **rate** target-variable it returns a multiplier via :meth:`rate_multiplier`
    (``rate = base · multiplier``); applied to a **state**/optimum target it returns a value or
    offset via :meth:`state_offset`. A response is a **null** (:attr:`is_null`) when it does not
    depend on the driver at all — the matched null of every edge.
    """

    @abstractmethod
    def rate_multiplier(self, driver_value) -> float:
        """The multiplier on a *rate* target-variable for this driver value."""

    def state_offset(self, driver_value) -> float:
        """The value/offset for a *state* target-variable (e.g. an OU optimum). Optional; a
        response used only on rate targets need not define it."""
        raise NotImplementedError(
            f"{type(self).__name__} does not define a state (jump/optimum) effect")

    @property
    def is_null(self) -> bool:
        """Whether the response is independent of the driver (the matched null)."""
        return False


@dataclass
class Scalar(Response):
    """The default one-knob response: ``multiplier = exp(strength · driver)`` (a log link), and a
    linear optimum shift ``strength · driver`` for state targets.

    ``strength`` is the single interpretable coupling coefficient — positive raises the target
    rate with the driver, negative lowers it, and ``strength = 0`` is the null (multiplier ``1``
    for every driver value).
    """

    strength: float

    def rate_multiplier(self, driver_value) -> float:
        return math.exp(_clamp(self.strength * float(driver_value)))

    def state_offset(self, driver_value) -> float:
        return self.strength * float(driver_value)

    @property
    def is_null(self) -> bool:
        return self.strength == 0.0


@dataclass
class Table(Response):
    """A per-discrete-state response — recovers MuSSE free per-state rates and a
    ``theta_present``/``theta_absent`` switch. Maps a driver **state** (any hashable key) to its
    multiplier (rate targets) or its value (state targets); unlisted states fall back to
    ``default``."""

    per_state: dict
    default: float = 1.0

    def __post_init__(self) -> None:
        self.per_state = {k: float(v) for k, v in dict(self.per_state).items()}
        self.default = float(self.default)

    def rate_multiplier(self, driver_value) -> float:
        return self.per_state.get(driver_value, self.default)

    def state_offset(self, driver_value) -> float:
        return self.per_state.get(driver_value, self.default)

    @property
    def is_null(self) -> bool:
        # No dependence on the driver state → no signal.
        return len(set(self.per_state.values()) | {self.default}) == 1


@dataclass
class Curve(Response):
    """A nonlinear continuous response — recovers QuaSSE. ``fn(driver_value)`` gives the multiplier
    (or state value), optionally capped at ``bound`` (the rate ceiling a Gillespie thinner needs)."""

    fn: Callable[[float], float]
    bound: float | None = None

    def rate_multiplier(self, driver_value) -> float:
        y = float(self.fn(float(driver_value)))
        if self.bound is not None:
            y = min(y, self.bound)
        return y

    def state_offset(self, driver_value) -> float:
        return float(self.fn(float(driver_value)))


def null_response() -> Scalar:
    """The matched null: ``Scalar(0.0)`` — every driver value maps to a multiplier of ``1``."""
    return Scalar(0.0)


# ═══════════════════════════════════════════════════════════════════════════════
# The sentence: driver → target-variable : response
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class Driver:
    """Whose state pushes. ``kind`` is ``"state"`` (a value along the tree — a trait value, a gene
    count) or ``"event"`` (a set of instants — speciation)."""

    level: str
    kind: str = "state"

    def __post_init__(self) -> None:
        if self.level not in LEVELS:
            raise ValueError(f"unknown driver level {self.level!r}; expected one of {LEVELS}")
        if self.kind not in ("state", "event"):
            raise ValueError(f"driver kind must be 'state' or 'event', got {self.kind!r}")


@dataclass(frozen=True)
class TargetVariable:
    """The variable of a level that a coupling bends, from the closed per-level menu
    (:data:`TARGET_VARIABLES`). Its ``kind`` (``"rate"`` or ``"state"``) is looked up, not supplied."""

    level: str
    variable: str

    def __post_init__(self) -> None:
        menu = TARGET_VARIABLES.get(self.level)
        if menu is None:
            raise ValueError(f"unknown target level {self.level!r}; expected one of {LEVELS}")
        if self.variable not in menu:
            raise ValueError(
                f"{self.level!r} has no target-variable {self.variable!r}; the closed menu is "
                f"{sorted(menu)}")

    @property
    def kind(self) -> str:
        """``"rate"`` (a modulation) or ``"state"`` (a jump)."""
        return TARGET_VARIABLES[self.level][self.variable]


@dataclass(frozen=True)
class Coupling:
    """One grammar sentence: ``driver → target-variable : response``.

    Rejects a tier-skipping pair on construction (the topology rule); whether it *layers* or
    *fuses* is a property of the whole :class:`CouplingGraph`, not of a coupling in isolation.
    """

    driver: Driver
    target: TargetVariable
    response: Response

    def __post_init__(self) -> None:
        if not _adjacent(self.driver.level, self.target.level):
            raise ValueError(
                f"forbidden coupling {self.driver.level}→{self.target.level}: it skips a tier. "
                f"The species↔sequences diagonal is forbidden — a sequence rides gene trees, not "
                f"the species tree (to recapitulate the species tree, simulate one gene family "
                f"with no events).")

    @property
    def is_null(self) -> bool:
        """Whether this edge's response is decoupled (the matched null)."""
        return self.response.is_null


def _adjacent(a: str, b: str) -> bool:
    """The topology rule: two levels may couple only if within one tier of each other."""
    return abs(_TIER[a] - _TIER[b]) <= 1


def couple(driver: str, target: str, variable: str, response, *, driver_kind: str = "state") -> Coupling:
    """Sugar for building a :class:`Coupling` from level names — e.g.
    ``couple("traits", "genomes", "loss", -0.8)``. A bare number ``response`` becomes a
    :class:`Scalar` of that strength."""
    resp = response if isinstance(response, Response) else Scalar(float(response))
    return Coupling(Driver(driver, driver_kind), TargetVariable(target, variable), resp)


# ═══════════════════════════════════════════════════════════════════════════════
# The coupling graph: layer vs fuse
# ═══════════════════════════════════════════════════════════════════════════════
def _reachable(start: str, adj: dict[str, set[str]]) -> set[str]:
    seen: set[str] = set()
    stack = [start]
    while stack:
        u = stack.pop()
        for v in adj.get(u, ()):
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return seen


class _SCC:
    """Strongly-connected components of a tiny directed graph, by mutual reachability (the graph
    has ≤ 4 nodes, so the O(V·(V+E)) reachability method is both simplest and plenty fast)."""

    def __init__(self, nodes: Iterable[str], edges: Iterable[tuple[str, str]]):
        adj: dict[str, set[str]] = {}
        for u, v in edges:
            adj.setdefault(u, set()).add(v)
        nodes = list(nodes)
        reach = {n: _reachable(n, adj) | {n} for n in nodes}  # m in reach[n] ⇔ n reaches m
        self.component: dict[str, int] = {}
        self.members: dict[int, frozenset[str]] = {}
        cid = 0
        for n in nodes:
            if n in self.component:
                continue
            grp = frozenset(m for m in nodes if m in reach[n] and n in reach[m])  # mutually reachable
            for m in grp:
                self.component[m] = cid
            self.members[cid] = grp
            cid += 1

    def same(self, a: str, b: str) -> bool:
        """Whether ``a`` and ``b`` sit in one SCC of size > 1 (a genuine cycle, not a lone node)."""
        ca = self.component[a]
        return ca == self.component[b] and len(self.members[ca]) > 1


class CouplingGraph:
    """A set of couplings plus the two rules that read them.

    Topology is validated per-coupling on construction (:class:`Coupling`). This class classifies
    each coupling — and the whole run — as **directional** (layer, run in order) or **bidirectional**
    (fuse, run together). A coupling *fuses* when its two levels sit in one strongly-connected
    component of the combined graph of the *coupling* edges plus the implicit *substrate* edges
    (:data:`_SUBSTRATE_EDGES`). That one rule captures both cases:

    * an explicit ``A → B`` and ``B → A`` cycle (feedback, e.g. the trait↔gene loop);
    * a single arrow *into a substrate* (e.g. ``traits → species``), which closes a cycle with the
      substrate edge ``species → traits`` — the reason all into-species coupling grows the tree.
    """

    def __init__(self, couplings: Iterable[Coupling]):
        self.couplings: list[Coupling] = list(couplings)
        for c in self.couplings:
            if not isinstance(c, Coupling):
                raise TypeError(f"expected Coupling instances, got {type(c).__name__}")
        self._scc = _SCC(LEVELS, self._edges())

    def _edges(self) -> set[tuple[str, str]]:
        edges = set(_SUBSTRATE_EDGES)
        for c in self.couplings:
            edges.add((c.driver.level, c.target.level))
        return edges

    # --- classification ----------------------------------------------------
    def is_fused(self, coupling: Coupling) -> bool:
        """Whether ``coupling`` must be co-integrated (bidirectional) rather than layered."""
        return self._scc.same(coupling.driver.level, coupling.target.level)

    def fused_groups(self) -> list[list[Coupling]]:
        """The fused couplings grouped by the cycle (SCC) they belong to — each group is
        co-integrated together."""
        by_component: dict[int, list[Coupling]] = {}
        for c in self.couplings:
            if self.is_fused(c):
                by_component.setdefault(self._scc.component[c.driver.level], []).append(c)
        return list(by_component.values())

    def layered(self) -> list[Coupling]:
        """The directional couplings — simulated in order, driver before target."""
        return [c for c in self.couplings if not self.is_fused(c)]

    @property
    def grows_tree(self) -> bool:
        """Whether any coupling targets a *species* rate — i.e. the species tree is an OUTPUT
        (grown) rather than an INPUT (an overlay on a given tree)."""
        return any(c.target.level == "species" for c in self.couplings)

    @property
    def mode(self) -> str:
        """``"bidirectional"`` if any coupling fuses, else ``"directional"``."""
        return "bidirectional" if any(self.is_fused(c) for c in self.couplings) else "directional"
