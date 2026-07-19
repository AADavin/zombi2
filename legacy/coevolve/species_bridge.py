"""Bridge: grammar couplings onto the into-species (tree-growing) engines.

The overlay edges compile to a rate :class:`~zombi2.genomes.rates.Modifier`, an OU walk, or a clock
on a **given** tree. The **into-species** edges are different — they *grow* the tree: a driver's
state sets the diversification rates, so the tree and the driver co-evolve in one forward Gillespie.
That is the grammar's **fuse** path (an arrow into ``species`` closes a cycle with the substrate), and
it cannot be an overlay.

Unlike the overlay rate-models, the into-species forward Gillespie is irreducible — there is no
bespoke duplication to delete — so this bridge is a thin, grammar-native **front-end** over the
existing SSE engine (:mod:`zombi2.coevolve.sse`), reused unchanged: it reads a grammar
:class:`~zombi2.coevolve.grammar.Response` on ``species.speciation`` / ``species.extinction`` as the
per-state birth/death rates (a free per-state :class:`~zombi2.coevolve.grammar.Table` = MuSSE) and
runs the tree-growing simulation. This realizes the **traits:species** edge end-to-end through the
grammar. See ``docs/design/coevolve-grammar.md`` §3.1 (traits:species) and §4.2.
"""

from __future__ import annotations

import numpy as np

from zombi2.coevolve.gene_diversification import (
    GeneDiversification, simulate_co_diversification, simulate_gene_diversification,
)
from zombi2.coevolve.grammar import Jump, Response, Scalar
from zombi2.coevolve.sse import MuSSE, simulate_sse


def musse_from_responses(states, transition, speciation: Response, extinction: Response) -> MuSSE:
    """Build a :class:`~zombi2.coevolve.sse.MuSSE` whose per-state birth/death rates are read from
    grammar responses.

    ``states`` are the discrete driver states; ``transition`` is their ``k×k`` rate matrix ``Q`` (the
    driver trait's own Mk dynamics). ``speciation`` / ``extinction`` are grammar
    :class:`~zombi2.coevolve.grammar.Response` s giving each state's rate — a
    :class:`~zombi2.coevolve.grammar.Table` is a free per-state MuSSE (``Table({0: λ0, 1: λ1, …})``);
    any response works, evaluated per state via :meth:`Response.rate_multiplier`.
    """
    states = list(states)
    birth = [float(speciation.rate_multiplier(s)) for s in states]
    death = [float(extinction.rate_multiplier(s)) for s in states]
    return MuSSE(birth=birth, death=death, Q=np.asarray(transition, dtype=float), states=states)


def simulate_trait_driven_diversification(states, transition, speciation: Response,
                                          extinction: Response, *, cladogenesis: Jump | None = None,
                                          age=None, n_tips=None, root_state=None, seed=None, rng=None):
    """Grow a tree whose diversification is driven by a discrete trait — the **traits:species** edge,
    run through the grammar. With ``cladogenesis`` it becomes the **ClaSSE** joint (``traits:species``
    *and* ``species:traits``: the trait also jumps at each speciation).

    A driver trait over ``states`` (with transition matrix ``transition``) sets each lineage's
    speciation/extinction from the grammar ``speciation`` / ``extinction`` responses; the tree and
    the trait grow together in the SSE forward Gillespie. ``cladogenesis`` is an optional grammar
    :class:`~zombi2.coevolve.grammar.Jump` (``scale`` = Gaussian jump variance for a continuous
    trait, ``probability`` = discrete move probability). Give exactly one stopping condition
    (``age`` or ``n_tips``); returns the :class:`~zombi2.traits.TraitResult`.
    """
    model = musse_from_responses(states, transition, speciation, extinction)
    return simulate_sse(model, age=age, n_tips=n_tips, root_state=root_state,
                        cladogenesis=_cladogenesis_from_jump(cladogenesis), seed=seed, rng=rng)


def _scalar_coefficient(response: Response, which: str) -> float:
    """The per-driver log-rate coefficient β from a :class:`~zombi2.coevolve.grammar.Scalar` response.
    ``genes:species`` is an exp-link (``λ = λ0·exp(Σ β_d)``), so the response must be a ``Scalar``."""
    if not isinstance(response, Scalar):
        raise TypeError(
            f"genes:species uses an exp-link on species.{which}: give a Scalar response (its "
            f"strength is the per-driver log-rate coefficient β), got {type(response).__name__}")
    return response.strength


def _jump_gain_to_probability(gain: float, k: int) -> float:
    """Realise a :class:`~zombi2.coevolve.grammar.Jump`'s ``gain`` (an *expected count* of families
    gained at a split) as the per-driver probability ``GeneDiversification`` wants.

    ``GeneDiversification`` bursts over a **finite** panel of ``k`` binary drivers: each *absent*
    driver is gained with probability ``p`` (see ``_burst_drivers``), so ``E[gained] = n_absent·p``.
    Setting ``p = gain/k`` makes ``E[gained] = gain`` for a lineage carrying none of the panel, and
    scales it down in proportion to what the lineage already carries — you cannot gain a driver you
    already have. This is the finite-panel analogue of the unbounded Poisson mean the same ``gain``
    field denotes for :func:`simulate_cladogenetic_genomes`; it is clamped to 1.0, so a ``gain``
    above ``k`` saturates at "gain every absent driver" rather than raising.
    """
    if k <= 0 or gain <= 0.0:
        return 0.0
    return min(1.0, gain / k)


def simulate_gene_driven_diversification(n_drivers, *, speciation: Response,
                                         extinction: Response | None = None,
                                         cladogenesis: Jump | None = None,
                                         lambda0: float = 1.0, mu0: float = 0.2,
                                         loss: float = 0.1, origination: float = 0.05,
                                         transfer: float = 0.5, root_drivers=0,
                                         age=None, n_tips=None, seed=None, rng=None):
    """Grow a tree whose diversification is driven by gene content — the **genes:species**
    (key-innovation) edge, run through the grammar. With ``cladogenesis`` it becomes the
    **co-diversification** joint (``genes:species`` *and* ``species:genomes``: a founder burst
    reshuffles the drivers at each split).

    ``K = n_drivers`` binary driver families set each lineage's speciation/extinction through an
    **exp-link** (``λ(S) = λ0·exp(Σ_{d∈S} βλ_d)``), so ``speciation`` / ``extinction`` are grammar
    :class:`~zombi2.coevolve.grammar.Scalar` responses whose ``strength`` is the per-driver
    coefficient βλ / βμ. ``cladogenesis`` is an optional grammar
    :class:`~zombi2.coevolve.grammar.Jump` giving the at-split burst (``probability`` = per-driver
    drop probability, ``gain`` = **mean drivers gained**, realised over the finite ``K``-driver panel
    by :func:`_jump_gain_to_probability`). Give exactly one stopping condition; returns a
    :class:`~zombi2.coevolve.gene_diversification.GeneDiversificationResult`. The engines are reused
    unchanged.
    """
    burst = cladogenesis if (cladogenesis is not None and not cladogenesis.is_null) else None
    model = GeneDiversification(
        n_drivers, lambda0=lambda0, mu0=mu0,
        driver_speciation=_scalar_coefficient(speciation, "speciation"),
        driver_extinction=(_scalar_coefficient(extinction, "extinction") if extinction is not None
                           else 0.0),
        loss=loss, origination=origination, transfer=transfer, root_drivers=root_drivers,
        cladogenetic_loss=(burst.probability if burst else 0.0),
        # `gain` is an expected COUNT everywhere in the grammar; the engine wants a per-driver
        # probability over its finite panel, so convert rather than pass it through raw.
        cladogenetic_gain=(_jump_gain_to_probability(burst.gain, n_drivers) if burst else 0.0))
    engine = simulate_co_diversification if burst is not None else simulate_gene_diversification
    return engine(model, age=age, n_tips=n_tips, seed=seed, rng=rng)


# ═══════════════════════════════════════════════════════════════════════════════
# species:X (cladogenetic) — speciation reshapes a character on a GIVEN tree
# ═══════════════════════════════════════════════════════════════════════════════
# The reverse direction: species drives a character. The driver is the speciation EVENT (not a
# state) and the effect is a grammar JUMP (:class:`~zombi2.coevolve.grammar.Jump`) at each split — an
# overlay on a given tree (a *layer* edge, unlike the into-species edges above). These reuse the
# existing cladogenetic engines unchanged.
def _cladogenesis_from_jump(jump: Jump | None):
    """A traits ``Cladogenesis`` kernel from a grammar :class:`~zombi2.coevolve.grammar.Jump`
    (``scale`` → Gaussian jump variance, ``probability`` → discrete shift), or ``None`` for no jump."""
    if jump is None or jump.is_null:
        return None
    from zombi2.traits.models import Cladogenesis
    return Cladogenesis(jump_sigma2=jump.scale, shift=jump.probability)


def simulate_cladogenetic_trait(tree, model, jump: Jump, *, root_state=None, seed=None, rng=None):
    """The **species:traits** edge: speciation reshapes the trait. At each branching a daughter's
    trait jumps per the grammar ``jump`` (``scale`` = Gaussian variance for a continuous trait,
    ``probability`` = move probability for a discrete one), layered on the anagenetic ``model``. An
    overlay on a given ``tree``. Reuses :func:`~zombi2.simulate_traits` with a ``Cladogenesis``.
    """
    from zombi2.traits.models import simulate_traits
    return simulate_traits(tree, model, cladogenesis=_cladogenesis_from_jump(jump),
                           root_state=root_state, seed=seed, rng=rng)


def simulate_cladogenetic_genomes(tree, jump: Jump, *, initial_families: int, loss: float = 0.0,
                                  origination: float = 0.0, seed=None, rng=None):
    """The **species:genomes** edge: speciation reshuffles the genome. At each split every daughter
    independently drops each family it carries with probability ``jump.probability`` and gains
    ``Poisson(jump.gain)`` new families (a founder burst), plus anagenetic ``loss`` / ``origination``
    along branches. An overlay on a given ``tree``. Reuses the
    :class:`~zombi2.coevolve.cladogenetic_genome.CladogeneticGenome` engine.
    """
    from zombi2.coevolve.cladogenetic_genome import (
        CladogeneticGenome, simulate_cladogenetic_genome,
    )
    model = CladogeneticGenome(initial_families=initial_families, loss=loss, origination=origination,
                               cladogenetic_loss=jump.probability, cladogenetic_gain=jump.gain)
    return simulate_cladogenetic_genome(tree, model, seed=seed, rng=rng)
