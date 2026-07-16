"""Tests for the declarative coevolution grammar (:mod:`zombi2.coevolve.grammar`).

Three layers:

1. **Response** — the one-knob :class:`Scalar` exp-link (and its ``strength=0`` null), the
   per-state :class:`Table`, and the nonlinear :class:`Curve`.
2. **The sentence** — building and validating ``driver → target-variable : response``: the closed
   target-variable menu and the topology rule (the forbidden ``species↔sequences`` diagonal).
3. **The graph** — the solve rule: directional (layer) vs bidirectional (fuse), including the
   "into a substrate is a cycle" case that makes into-species coupling grow the tree.
"""

import math

import pytest

from zombi2.coevolve.grammar import (
    LEVELS, TARGET_VARIABLES, Coupling, CouplingGraph, Curve, Driver, DriverSignal, Jump, Scalar, Table,
    TargetVariable, couple, legal_null_kinds, make_null, null_response,
)


# ── 1. Response ───────────────────────────────────────────────────────────────
def test_scalar_is_exp_link():
    assert Scalar(0.8).rate_multiplier(1.0) == pytest.approx(math.exp(0.8))
    assert Scalar(-0.5).rate_multiplier(2.0) == pytest.approx(math.exp(-1.0))


def test_scalar_zero_is_the_null():
    r = Scalar(0.0)
    assert r.is_null
    for x in (-3.0, 0.0, 5.0, 1e6):
        assert r.rate_multiplier(x) == 1.0          # no dependence on the driver
    assert null_response().is_null
    assert null_response().rate_multiplier(42.0) == 1.0


def test_scalar_clamps_the_exponent():
    # An extreme strength·driver must not overflow exp().
    huge = Scalar(1e3).rate_multiplier(1e3)
    assert math.isfinite(huge)
    assert huge == pytest.approx(math.exp(40.0))     # clamped at _MAX_EXPONENT
    assert Scalar(-1e3).rate_multiplier(1e3) == pytest.approx(math.exp(-40.0))


def test_scalar_optimum_offset_is_linear():
    assert Scalar(2.0).state_offset(3.0) == pytest.approx(6.0)


def test_table_lookup_default_and_null():
    t = Table({"marine": 1.0, "fresh": 2.5, "soil": 0.6})
    assert t.rate_multiplier("fresh") == 2.5
    assert t.rate_multiplier("desert") == 1.0        # falls back to default
    assert not t.is_null
    assert Table({0: 3.0, 1: 3.0}, default=3.0).is_null   # uniform → no signal
    assert not Table({0: 1.0, 1: 3.0}).is_null


def test_curve_evaluates_and_caps_at_bound():
    hump = Curve(lambda x: 4.0 - (x - 2.0) ** 2, bound=5.0)
    assert hump.rate_multiplier(2.0) == pytest.approx(4.0)     # peak
    assert hump.rate_multiplier(0.0) == pytest.approx(0.0)     # 4 - 4
    capped = Curve(lambda x: 100.0, bound=5.0)
    assert capped.rate_multiplier(1.0) == 5.0


# ── 2. The sentence: validation ───────────────────────────────────────────────
def test_target_variable_kind_is_looked_up():
    assert TargetVariable("genomes", "loss").kind == "rate"
    assert TargetVariable("traits", "optimum").kind == "state"
    assert TargetVariable("species", "speciation").kind == "rate"


def test_unknown_target_variable_is_rejected():
    with pytest.raises(ValueError, match="no target-variable"):
        TargetVariable("genomes", "speciation")      # speciation is a species variable, not genomes
    with pytest.raises(ValueError, match="unknown target level"):
        TargetVariable("populations", "loss")


def test_unknown_driver_level_or_kind_is_rejected():
    with pytest.raises(ValueError, match="unknown driver level"):
        Driver("populations")
    with pytest.raises(ValueError, match="driver kind"):
        Driver("traits", kind="continuous")


def test_topology_rule_forbids_species_sequence_diagonal():
    with pytest.raises(ValueError, match="forbidden"):
        couple("species", "sequences", "residues", 1.0, driver_kind="event")
    with pytest.raises(ValueError, match="forbidden"):
        couple("sequences", "species", "speciation", 1.0)


@pytest.mark.parametrize("driver,target,variable", [
    ("traits", "species", "speciation"),      # SSE
    ("genomes", "species", "extinction"),     # key innovation
    ("species", "traits", "value"),           # cladogenetic
    ("traits", "genomes", "loss"),            # trait-linked
    ("genomes", "traits", "optimum"),         # gene-conditioned
    ("traits", "sequences", "selection"),     # T→Σ
    ("genomes", "sequences", "substitution_speed"),  # G→Σ
    ("sequences", "genomes", "loss"),         # Σ→G (concerted)
])
def test_all_diamond_edges_are_constructible(driver, target, variable):
    c = couple(driver, target, variable, 0.7)
    assert isinstance(c, Coupling)
    assert c.target.level == target


def test_couple_sugar_wraps_a_bare_number_as_scalar():
    c = couple("traits", "genomes", "loss", -0.8)
    assert isinstance(c.response, Scalar)
    assert c.response.strength == -0.8


def test_coupling_is_null_when_response_is_null():
    assert couple("traits", "genomes", "loss", 0.0).is_null
    assert not couple("traits", "genomes", "loss", 0.8).is_null


# ── 3. The graph: layer vs fuse ───────────────────────────────────────────────
def test_directional_edge_layers():
    g = CouplingGraph([couple("traits", "genomes", "loss", -0.8)])
    assert g.mode == "directional"
    assert not g.grows_tree
    assert len(g.layered()) == 1
    assert g.fused_groups() == []


def test_trait_gene_feedback_fuses():
    a = couple("traits", "genomes", "loss", -0.8)
    b = couple("genomes", "traits", "optimum", 0.8)
    g = CouplingGraph([a, b])
    assert g.mode == "bidirectional"
    assert g.is_fused(a) and g.is_fused(b)
    groups = g.fused_groups()
    assert len(groups) == 1 and len(groups[0]) == 2    # the T↔G pair co-integrates together


def test_into_species_edge_fuses_via_the_substrate_cycle():
    # A single arrow into species is bidirectional in disguise: it closes a cycle with the
    # implicit substrate edge species→traits, so the tree is GROWN.
    sse = couple("traits", "species", "speciation", 1.2)
    g = CouplingGraph([sse])
    assert g.is_fused(sse)
    assert g.grows_tree
    assert g.mode == "bidirectional"


def test_gene_to_sequence_layers_but_sequence_to_gene_fuses():
    # G→Σ rides the substrate downstream → directional; Σ→G closes the genomes↔sequences cycle
    # (concerted evolution) → fuse.
    down = CouplingGraph([couple("genomes", "sequences", "selection", 0.5)])
    assert not down.is_fused(down.couplings[0])

    up = couple("sequences", "genomes", "loss", 0.5)
    assert CouplingGraph([up]).is_fused(up)


def test_trait_to_sequence_is_directional():
    t2s = couple("traits", "sequences", "selection", 0.6)
    g = CouplingGraph([t2s])
    assert not g.is_fused(t2s)
    assert g.mode == "directional"


def test_menu_and_levels_are_the_diamond():
    assert set(TARGET_VARIABLES) == set(LEVELS)
    # species exposes only its two diversification rates; no state jump target.
    assert set(TARGET_VARIABLES["species"]) == {"speciation", "extinction"}


# ── 4. The null layer ─────────────────────────────────────────────────────────
def _state_edge():
    return couple("traits", "genomes", "loss", -0.8)                 # state driver


def _event_edge():
    return couple("species", "traits", "value", 1.0, driver_kind="event")   # event driver (clado)


def test_null_legality_matrix_by_driver_archetype():
    assert legal_null_kinds(_state_edge()) == frozenset({"neutral", "cid"})
    assert legal_null_kinds(_event_edge()) == frozenset({"neutral", "timing"})


def test_neutral_null_zeroes_the_response():
    n = make_null(_state_edge(), "neutral")
    assert n.kind == "neutral"
    assert n.coupling.is_null and isinstance(n.coupling.response, Scalar)
    assert n.coupling.response.strength == 0.0
    # neutral is legal for an event driver too
    assert make_null(_event_edge(), "neutral").coupling.is_null


def test_cid_is_a_hidden_driver_transform_not_a_decoupling():
    c = _state_edge()
    n = make_null(c, "cid")
    assert n.kind == "cid"
    assert n.coupling.driver.hidden is True           # observed driver swapped for a hidden twin
    assert n.coupling.driver.level == c.driver.level and n.coupling.driver.kind == "state"
    # the response (the signal) is UNCHANGED — cid keeps the heterogeneity, just hides its cause
    assert n.coupling.response.strength == c.response.strength
    assert not n.coupling.is_null


def test_cid_is_illegal_for_an_event_driver():
    with pytest.raises(ValueError, match="cid.*not legal"):
        make_null(_event_edge(), "cid")


def test_timing_is_legal_only_for_event_drivers():
    n = make_null(_event_edge(), "timing")
    assert n.kind == "timing"
    assert n.coupling is n.original                    # a marker; engine spreads it at sim time
    with pytest.raises(ValueError, match="timing.*not legal"):
        make_null(_state_edge(), "timing")


def test_unknown_null_kind_is_rejected():
    with pytest.raises(ValueError, match="unknown null kind"):
        make_null(_state_edge(), "shuffle")


# ── 5. The execution plan ─────────────────────────────────────────────────────
def _flatten(plan):
    out = []
    for mode, payload in plan:
        out.extend(payload if mode == "fuse" else [payload])
    return out


def test_plan_orders_directional_edges_by_dependency():
    a = couple("traits", "genomes", "loss", -0.8)
    b = couple("genomes", "sequences", "selection", 0.5)
    plan = CouplingGraph([b, a]).solve_plan()          # deliberately reversed input
    assert [m for m, _ in plan] == ["layer", "layer"]
    # traits→genomes must precede genomes→sequences (its driver, genomes, is produced first)
    assert [p for _, p in plan] == [a, b]


def test_plan_emits_a_cycle_as_one_fused_step():
    a = couple("traits", "genomes", "loss", -0.8)
    b = couple("genomes", "traits", "optimum", 0.8)
    plan = CouplingGraph([a, b]).solve_plan()
    assert len(plan) == 1
    mode, group = plan[0]
    assert mode == "fuse"
    assert len(group) == 2 and a in group and b in group    # membership by ==, not hashing


def test_plan_includes_a_layered_edge_feeding_into_a_cycle():
    # Regression: a layered (directional) edge whose TARGET sits inside a fused cycle must still
    # appear in the plan. species→traits is layered; traits↔genomes is a fused cycle containing
    # traits — the layered edge must not be dropped.
    lay = couple("species", "traits", "value", 0.3, driver_kind="event")
    a = couple("traits", "genomes", "loss", -0.8)
    b = couple("genomes", "traits", "optimum", 0.5)
    plan = CouplingGraph([a, b, lay]).solve_plan()
    flat = _flatten(plan)
    assert len(flat) == 3 and lay in flat and a in flat and b in flat   # every coupling, once
    assert [m for m, _ in plan].count("fuse") == 1                      # the cycle is one step
    assert ("layer", lay) in plan


def test_plan_covers_every_coupling_exactly_once():
    sse = couple("traits", "species", "speciation", 1.2)       # fused (into species)
    down = couple("genomes", "sequences", "selection", 0.4)    # layered
    g = CouplingGraph([sse, down])
    flat = _flatten(g.solve_plan())
    assert sorted(map(id, flat)) == sorted(map(id, [sse, down]))
    # the into-species fuse comes before the downstream layer step
    assert g.solve_plan()[0][0] == "fuse"


# ── 6. The driver-signal seam ─────────────────────────────────────────────────
def test_trait_trajectory_satisfies_the_driver_signal_protocol():
    # The grammar's DriverSignal contract must match the existing TraitTrajectory, so the eventual
    # Driver↔trajectory bridge needs only an adapter, not a new signal type.
    from zombi2.coevolve.trait_coupling import TraitTrajectory
    traj = TraitTrajectory({}, {}, [], default=0.0)            # a trivial (empty) trajectory
    assert isinstance(traj, DriverSignal)
    assert traj.value("any_lineage", 1.0) == 0.0              # falls back to the default value
    assert traj.refresh_times(0.0, 1.0) == []


# ── 7. The Jump response (event-driver state change) ──────────────────────────
def test_jump_response_null_and_validation():
    assert Jump().is_null                                     # no jump
    assert not Jump(scale=1.0).is_null
    assert not Jump(probability=0.5).is_null
    assert not Jump(gain=2.0).is_null
    with pytest.raises(ValueError, match="scale"):
        Jump(scale=-1.0)
    with pytest.raises(ValueError, match="probability"):
        Jump(probability=1.5)
    with pytest.raises(ValueError, match="gain"):
        Jump(gain=-1.0)


def test_jump_is_not_a_rate_multiplier():
    with pytest.raises(NotImplementedError):
        Jump(scale=1.0).rate_multiplier(0.5)                  # a jump is a state change, not a rate


# --- the two null encodings must be kept in step BY HAND (2026-07-16 audit) ---------------------

def test_grammar_null_matrix_agrees_with_the_models_that_actually_run():
    """``legal_null_kinds`` and the production models' ``.null()`` ladders must not drift apart.

    The grammar's null layer is **not wired up** (see the module note): every model still carries the
    per-edge if/error ladder ``make_null``/``legal_null_kinds`` was written to replace, and the CLI
    calls *those*. So the same rule is encoded twice and kept in step by hand — which is precisely
    what this test guards until P3 adopts the matrix (docs/design/coevolve-grammar.md §4.4).
    """
    import zombi2 as z
    from zombi2.coevolve.grammar import couple, legal_null_kinds

    def accepted(model, kinds):
        out = set()
        for k in kinds:
            try:
                model.null(k)
            except (ValueError, TypeError):
                continue
            out.add(k)
        return out

    edge = couple("traits", "species", "speciation", 1.0)          # a *state* driver
    assert legal_null_kinds(edge) == frozenset({"neutral", "cid"})

    musse = z.MuSSE(birth=[1.0, 2.0], death=[0.2, 0.2], Q=[[0.0, 0.1], [0.1, 0.0]])
    assert accepted(musse, ("neutral", "cid", "timing")) == legal_null_kinds(edge)


def test_grammar_null_matrix_knowingly_over_permits_a_continuous_driver():
    """The one place the matrix and the models *do* disagree, pinned so it stays deliberate.

    ``legal_null_kinds`` keys on the driver **archetype** (``"state"``), which does not yet carry the
    discrete-vs-continuous distinction — so it offers ``cid`` for QuaSSE, whose ``.null()`` correctly
    refuses it (CID is a discrete-character null). The grammar's own comment admits this. If P3 wires
    the matrix in without refining the archetype, QuaSSE + ``--null cid`` would be offered and then
    fail; this test is the reminder.
    """
    import zombi2 as z
    from zombi2.coevolve.grammar import couple, legal_null_kinds

    edge = couple("traits", "species", "speciation", 1.0)
    assert "cid" in legal_null_kinds(edge)                          # the matrix offers it ...

    quasse = z.QuaSSE(z.QuaSSE.sigmoid(0.5, 2.0), lambda x: 0.2, sigma2=0.4, rate_bound=2.2)
    quasse.null("neutral")                                          # ... neutral is fine ...
    with pytest.raises(TypeError, match="discrete-character null"):
        quasse.null("cid")                                          # ... but the model refuses it
