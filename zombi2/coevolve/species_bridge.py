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

from zombi2.coevolve.gene_diversification import GeneDiversification, simulate_gene_diversification
from zombi2.coevolve.grammar import Response, Scalar
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
                                          extinction: Response, *, age=None, n_tips=None,
                                          root_state=None, seed=None, rng=None):
    """Grow a tree whose diversification is driven by a discrete trait — the **traits:species** edge,
    run through the grammar.

    A driver trait over ``states`` (with transition matrix ``transition``) sets each lineage's
    speciation/extinction from the grammar ``speciation`` / ``extinction`` responses; the tree and
    the trait grow together in the SSE forward Gillespie. Give exactly one stopping condition
    (``age`` or ``n_tips``). Returns the :class:`~zombi2.traits.TraitResult` (the complete tree plus
    the realized state history); ``z.prune(result.tree)`` gives the survivors-only tree.
    """
    model = musse_from_responses(states, transition, speciation, extinction)
    return simulate_sse(model, age=age, n_tips=n_tips, root_state=root_state, seed=seed, rng=rng)


def _scalar_coefficient(response: Response, which: str) -> float:
    """The per-driver log-rate coefficient β from a :class:`~zombi2.coevolve.grammar.Scalar` response.
    ``genes:species`` is an exp-link (``λ = λ0·exp(Σ β_d)``), so the response must be a ``Scalar``."""
    if not isinstance(response, Scalar):
        raise TypeError(
            f"genes:species uses an exp-link on species.{which}: give a Scalar response (its "
            f"strength is the per-driver log-rate coefficient β), got {type(response).__name__}")
    return response.strength


def simulate_gene_driven_diversification(n_drivers, *, speciation: Response,
                                         extinction: Response | None = None,
                                         lambda0: float = 1.0, mu0: float = 0.2,
                                         loss: float = 0.1, origination: float = 0.05,
                                         transfer: float = 0.5, root_drivers=0,
                                         age=None, n_tips=None, seed=None, rng=None):
    """Grow a tree whose diversification is driven by gene content — the **genes:species**
    (key-innovation) edge, run through the grammar.

    ``K = n_drivers`` binary driver families set each lineage's speciation/extinction through an
    **exp-link**: ``λ(S) = λ0·exp(Σ_{d∈S} βλ_d)`` (and likewise μ), so ``speciation`` / ``extinction``
    are grammar :class:`~zombi2.coevolve.grammar.Scalar` responses whose ``strength`` is the
    per-driver coefficient βλ / βμ. The drivers spread by loss / origination / frequency-dependent
    ``transfer`` and ride the forward Gillespie *with* the growing tree (the fuse path). Give exactly
    one stopping condition (``age`` or ``n_tips``); returns a
    :class:`~zombi2.coevolve.gene_diversification.GeneDiversificationResult`
    (``.tree``, ``.tip_prevalence()``). The key-innovation engine is reused unchanged.
    """
    model = GeneDiversification(
        n_drivers, lambda0=lambda0, mu0=mu0,
        driver_speciation=_scalar_coefficient(speciation, "speciation"),
        driver_extinction=(_scalar_coefficient(extinction, "extinction") if extinction is not None
                           else 0.0),
        loss=loss, origination=origination, transfer=transfer, root_drivers=root_drivers)
    return simulate_gene_diversification(model, age=age, n_tips=n_tips, seed=seed, rng=rng)
