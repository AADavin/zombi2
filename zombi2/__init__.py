"""ZOMBI2 — simulating the evolution of species, genomes, sequences and traits.

There are **no top-level re-exports**: one canonical path per name, reached through each
level's package (see ``docs/design/MAP.md``). ::

    from zombi2 import species
    result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

    from zombi2.genomes import simulate_genomes_unordered
    from zombi2.rates import scope, modifiers

ZOMBI2 is being rebuilt as a **clean core grown from** ``docs/design/SPEC.md``; the old
code is quarantined under ``legacy/`` at the repo root, read-only, not importable. The
active tree is the clean core: ``zombi2.rates`` · ``zombi2.species`` · ``zombi2.genomes``.
"""

from __future__ import annotations

__version__ = "0.2.0"
