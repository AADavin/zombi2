"""Tests for the grammar → rate-engine bridge (:mod:`zombi2.coevolve.rate_bridge`).

Two layers:

1. **Unit — the factor** — a :class:`CouplingModifier` applies its :class:`Response` to its target
   event only, keys per-family weights the trait→gene way, and forwards the driver's change points.
2. **Integration — the engine** — composed through :class:`~zombi2.genomes.rates.ModifiedRates` over
   a real :class:`~zombi2.genomes.rates.Rates` (per='copy') base, it reproduces the trait→gene loss formula
   ``base_loss · cn · exp(-strength · w_f · s)`` per family, and a null coupling leaves the base
   weights byte-identical.
"""

import math

import pytest

from zombi2.coevolve.grammar import Scalar
from zombi2.coevolve.rate_bridge import EVENT_FOR_VARIABLE, CouplingModifier
from zombi2.genomes.events import EventType
from zombi2.genomes.genome import Gene, IdManager, UnorderedGenome
from zombi2.genomes.rates import ModifiedRates, Rates


class _ConstSignal:
    """A trivial DriverSignal: a constant value plus optional (time, branch) change points."""

    def __init__(self, value, changes=()):
        self._v = float(value)
        self._changes = list(changes)

    def value(self, lineage, time):
        return self._v

    def refresh_times(self, t0, t1):
        return [c for c in self._changes if t0 < c[0] < t1]


# ── 1. Unit: the factor ───────────────────────────────────────────────────────
def test_factor_applies_response_to_the_target_event_only():
    cm = CouplingModifier(_ConstSignal(1.0), Scalar(-0.8), EventType.LOSS)
    assert cm.factor(EventType.LOSS, None, "b", 0.0) == pytest.approx(math.exp(-0.8))
    assert cm.factor(EventType.DUPLICATION, None, "b", 0.0) == 1.0    # untargeted → no effect
    assert cm.factor(EventType.TRANSFER, None, "b", 0.0) == 1.0


def test_null_response_gives_factor_one():
    cm = CouplingModifier(_ConstSignal(5.0), Scalar(0.0), EventType.LOSS)
    assert cm.factor(EventType.LOSS, None, "b", 0.0) == 1.0


def test_per_family_weights_scale_the_driver_value():
    cm = CouplingModifier(_ConstSignal(2.0), Scalar(-1.0), EventType.LOSS, weights={"F1": 0.5})
    assert cm.keys_on_family is True
    # F1 sees weight·s = 0.5·2.0 = 1.0 → exp(-1.0·1.0)
    assert cm.factor(EventType.LOSS, "F1", "b", 0.0) == pytest.approx(math.exp(-1.0))
    # F2 is inert (weight 0) → no effect
    assert cm.factor(EventType.LOSS, "F2", "b", 0.0) == 1.0


def test_uniform_coupling_does_not_key_on_family():
    cm = CouplingModifier(_ConstSignal(1.0), Scalar(0.3), EventType.TRANSFER)
    assert cm.keys_on_family is False


def test_refresh_times_are_forwarded_from_the_signal():
    cm = CouplingModifier(_ConstSignal(1.0, changes=[(0.5, "b")]), Scalar(1.0), EventType.LOSS)
    assert cm.refresh_times(0.0, 1.0) == [(0.5, "b")]
    assert cm.refresh_times(0.6, 1.0) == []                            # change point is outside


def test_event_for_variable_maps_the_unambiguous_targets():
    assert EVENT_FOR_VARIABLE["loss"] is EventType.LOSS
    assert EVENT_FOR_VARIABLE["duplication"] is EventType.DUPLICATION
    assert EVENT_FOR_VARIABLE["transfer"] is EventType.TRANSFER


# ── 2. Integration: through ModifiedRates over a Rates (per='copy') base ──────
def _genome(copies):
    """An UnorderedGenome with ``copies = {family: n}``."""
    ids = IdManager()
    g = UnorderedGenome(ids)
    for fam, n in copies.items():
        for _ in range(n):
            g._add(Gene(ids.new_gene(), fam))
    return g


def _by_target(weights):
    return {(w.event, w.family): w.rate for w in weights}


def test_reproduces_trait_gene_loss_formula_per_family():
    # base_loss=1, F1 responsive (weight 1), F2 inert; constant trait s=1, strength -0.5.
    g = _genome({"F1": 2, "F2": 1})
    base = Rates(loss=1.0)
    cm = CouplingModifier(_ConstSignal(1.0), Scalar(-0.5), EventType.LOSS, weights={"F1": 1.0})
    weights = _by_target(ModifiedRates(base, [cm]).event_weights(g, "b", 0.0))

    # LOSS is expanded per family: rate = base_loss·cn, then ×exp(-strength·w_f·s).
    assert weights[(EventType.LOSS, "F1")] == pytest.approx(2.0 * math.exp(-0.5))   # cn=2, w=1
    assert weights[(EventType.LOSS, "F2")] == pytest.approx(1.0)                    # cn=1, inert


def test_null_coupling_is_byte_identical_to_the_base():
    g = _genome({"F1": 2, "F2": 3})
    base = Rates(duplication=0.4, transfer=0.6, loss=1.0, origination=0.1)
    cm = CouplingModifier(_ConstSignal(9.0), Scalar(0.0), EventType.LOSS)   # null: strength 0
    assert (ModifiedRates(base, [cm]).event_weights(g, "b", 0.0)
            == base.event_weights(g, "b", 0.0))


def test_only_the_targeted_channel_is_scaled():
    g = _genome({"F1": 2})
    base = Rates(transfer=0.5, loss=1.0)
    cm = CouplingModifier(_ConstSignal(1.0), Scalar(1.0), EventType.TRANSFER)   # scales TRANSFER only
    weights = _by_target(ModifiedRates(base, [cm]).event_weights(g, "b", 0.0))
    assert weights[(EventType.TRANSFER, None)] == pytest.approx(0.5 * 2 * math.exp(1.0))  # scaled
    assert weights[(EventType.LOSS, None)] == pytest.approx(1.0 * 2)                        # untouched
