"""ZOMBI2 — simulating the evolution of species, genomes, sequences and traits.

There are **no top-level re-exports**: one canonical path per name, reached through each
level's package (see ``docs/design/MAP.md``). ::

    from zombi2 import species
    result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

    from zombi2.genomes import simulate_genomes_unordered
    from zombi2.rates import scope, modifiers

ZOMBI2 is a **clean core grown from** ``docs/design/SPEC.md`` — pure Python, one concept
per name across the API, the CLI and a ``--params`` file: ``zombi2.rates`` · ``zombi2.species`` ·
``zombi2.genomes`` · ``zombi2.sequences`` · ``zombi2.traits`` · ``zombi2.joint`` (conditioned
coupling folds into the target level, so there is no ``coupling`` package — SPEC §2–4).
"""

from __future__ import annotations

__version__ = "0.4.0"
