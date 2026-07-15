"""Bridge: grammar couplings onto the sequence engine (the diamond's bottom tier).

Sequences evolve by molecular clocks (per-lineage substitution *rate*) and codon models (selection,
dN/dS ω) — a different engine from the genome rate machinery, so the sequence tier needs its own
bridge, not the :class:`~zombi2.coevolve.rate_bridge.CouplingModifier`.

Two target-variables live here:

* **``substitution_speed``** — how fast a lineage's sequences evolve. :class:`DriverClock` is a
  :class:`~zombi2.sequences.clocks.Clock` whose per-branch rate is a grammar
  :class:`~zombi2.coevolve.grammar.Response` of a :class:`~zombi2.coevolve.grammar.DriverSignal`, so
  it drops into :class:`~zombi2.SequenceEvolution`'s ``lineage_segments`` contract unchanged.
* **``selection``** (ω = dN/dS) — how tightly selection squeezes a coding sequence.
  :class:`OmegaSelector` (the **T→Σ** edge: a trait sets per-lineage ω) and :class:`GeneEventOmega`
  (the **G→Σ** edge: a gene event relaxes selection, e.g. post-duplication) both produce a
  per-gene-node codon model — an ``ω``-class cache — for
  :func:`~zombi2.sequences.models.evolve_on_tree`'s ``model_for`` hook.

Neither the clock nor the codon machinery lives in ``zombi2.genomes.rates``, so this whole tier is
independent of the rate rename. See ``docs/design/coevolve-grammar.md`` §5.
"""

from __future__ import annotations

from zombi2.coevolve.grammar import DriverSignal, Response
from zombi2.genomes.events import EventType
from zombi2.sequences.clocks import Clock
from zombi2.sequences.codon_models import gy94
from zombi2.tree import Tree


class DriverClock(Clock):
    """A molecular clock whose per-branch substitution rate is set by a grammar coupling on
    ``sequences.substitution_speed``.

    The rate on lineage ``b`` at time ``t`` is ``base_rate · response.rate_multiplier(driver_value)``,
    where the driver value is ``driver.value(b, t)``. Each branch is sub-segmented at the driver's
    interior change points (:meth:`DriverSignal.refresh_times`), so the rate tracks a within-branch
    driver change exactly. A null response (``Scalar(0)``) reduces this to a strict clock at
    ``base_rate``.

    Deterministic given the (already-simulated) driver — ``lineage_segments`` ignores its ``rng``.
    Satisfies the :class:`~zombi2.sequences.clocks.Clock` contract, so it is used exactly like any
    other clock (``.scale(tree)``, or as the shared clock in :class:`~zombi2.SequenceEvolution`).
    """

    def __init__(self, driver: DriverSignal, response: Response, *, base_rate: float = 1.0):
        if base_rate <= 0:
            raise ValueError(f"base_rate must be > 0, got {base_rate}")
        self.driver = driver
        self.response = response
        self.base_rate = float(base_rate)
        self.root_rate = self.base_rate

    def lineage_segments(self, tree: Tree, rng):
        segments: dict = {}
        avg: dict = {}
        for node in tree.nodes_preorder():
            if node.parent is None:
                segments[node.name] = []
                avg[node.name] = self.root_rate
                continue
            b0, b1 = node.parent.time, node.time
            cuts = sorted(t for (t, br) in self.driver.refresh_times(b0, b1)
                          if br == node.name and b0 < t < b1)
            bounds = [b0, *cuts, b1]
            segs = []
            rate_time = 0.0
            for s0, s1 in zip(bounds[:-1], bounds[1:]):
                if s1 <= s0:
                    continue
                r = self.base_rate * self.response.rate_multiplier(self.driver.value(node.name, s0))
                segs.append((r, s0, s1))
                rate_time += r * (s1 - s0)
            segments[node.name] = segs
            span = b1 - b0
            # A zero-length branch contributes no substitution length; report its instantaneous
            # driver-scaled rate (not a bare base_rate) so branch_rate matches sibling branches.
            avg[node.name] = (rate_time / span) if span > 0 else (
                self.base_rate * self.response.rate_multiplier(self.driver.value(node.name, b0)))
        return segments, avg


# ═══════════════════════════════════════════════════════════════════════════════
# selection (ω = dN/dS): a per-gene-node codon model, driven by the grammar
# ═══════════════════════════════════════════════════════════════════════════════
class _OmegaCoupling:
    """Shared machinery for a ``selection`` (ω) coupling: a lazily-built, cached codon model per
    **rounded** ω — the ω-class cache the docs call for, so a continuous driver does not mint a
    fresh 61×61 eigendecomposition per branch. Subclasses supply :meth:`omega_for` (the ω on a gene
    node); :meth:`model_for` and :attr:`base_model` plug into
    :func:`~zombi2.sequences.models.evolve_on_tree`.

    All ω-class models share the same stationary distribution / alphabet / ``k`` (they differ only in
    the non-synonymous rates), which is exactly what ``model_for`` requires.
    """

    def __init__(self, *, base_omega, kappa, freqs, builder, resolution):
        if base_omega < 0:
            raise ValueError(f"base_omega must be >= 0, got {base_omega}")
        if resolution <= 0:
            raise ValueError(f"resolution must be > 0, got {resolution}")
        self.base_omega = float(base_omega)
        self.kappa = float(kappa)
        self.freqs = freqs
        self.builder = builder
        self.resolution = float(resolution)
        self._cache: dict = {}
        #: a representative model for the root draw — shares π / alphabet / k with every ω class.
        self.base_model = self._model(self.base_omega)

    def _model(self, omega):
        key = round(max(0.0, float(omega)) / self.resolution)   # bucket ω → an ω class
        model = self._cache.get(key)
        if model is None:
            model = self.builder(kappa=self.kappa, omega=key * self.resolution, freqs=self.freqs)
            self._cache[key] = model
        return model

    def omega_for(self, node) -> float:
        raise NotImplementedError

    def model_for(self, node):
        """The codon model for ``node``'s ω class — pass as ``evolve_on_tree(..., model_for=…)``."""
        return self._model(self.omega_for(node))


class OmegaSelector(_OmegaCoupling):
    """**T→Σ selection**: a trait sets each gene lineage's dN/dS (ω).

    For every gene node, ``ω = base_omega · response.rate_multiplier(s)`` where ``s`` is the driver
    (trait) value on that node's **species branch** (``node.branch`` or, for a survivor,
    ``node.species``) at the branch midpoint — the same gene→species mapping
    :func:`~zombi2.sequences.evolution._annotate` uses. A positive-``strength`` ``Scalar`` makes a
    high trait *relax* selection (higher ω); the null (``Scalar(0)``) gives a uniform ω.
    """

    def __init__(self, driver: DriverSignal, response: Response, *, base_omega: float = 0.2,
                 kappa: float = 2.0, freqs=None, builder=gy94, resolution: float = 0.02):
        self.driver = driver
        self.response = response
        super().__init__(base_omega=base_omega, kappa=kappa, freqs=freqs, builder=builder,
                         resolution=resolution)

    def omega_for(self, node) -> float:
        name = node.branch if node.branch is not None else node.species
        s = self.driver.value(name, 0.5 * (node.birth + node.end)) if name is not None else 0.0
        return self.base_omega * self.response.rate_multiplier(s)


class GeneEventOmega(_OmegaCoupling):
    """**G→Σ selection**: a gene-content event relaxes (or tightens) selection on the branch it
    terminates — e.g. post-duplication relaxed selection.

    A gene node whose terminating event (``node.kind``) is one of ``events`` evolves at
    ``ω = base_omega · response.rate_multiplier(1)``; every other branch at
    ``base_omega · response.rate_multiplier(0) = base_omega``. So ``response=Scalar(β)`` relaxes an
    event branch's selection by ``e^β``.
    """

    def __init__(self, response: Response, *, base_omega: float = 0.2,
                 events=(EventType.DUPLICATION,), kappa: float = 2.0, freqs=None, builder=gy94,
                 resolution: float = 0.02):
        self.response = response
        self._events = frozenset(events)
        super().__init__(base_omega=base_omega, kappa=kappa, freqs=freqs, builder=builder,
                         resolution=resolution)

    def omega_for(self, node) -> float:
        driver_value = 1.0 if node.kind in self._events else 0.0
        return self.base_omega * self.response.rate_multiplier(driver_value)
