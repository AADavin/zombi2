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

from zombi2.coevolve.grammar import Response
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
