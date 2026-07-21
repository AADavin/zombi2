"""The ``zombi2`` command line — one subcommand per level, mirroring the Python API.

The CLI is grown from the clean core, one level at a time, exactly like the packages it wraps
(see ``docs/design/MAP.md``). Today it exposes the three built levels:

- ``zombi2 species``  → :func:`zombi2.species.simulate_species_tree`
- ``zombi2 genomes``  → :func:`zombi2.genomes.simulate_genomes_unordered` /
  :func:`~zombi2.genomes.simulate_genomes_ordered` (chosen by ``--resolution``)
- ``zombi2 traits``   → :func:`zombi2.traits.simulate_continuous` /
  :func:`~zombi2.traits.simulate_discrete` (chosen by ``--kind``)

Each subcommand's long options **are** the API keyword names (one word per concept across the API,
the CLI, and a ``--params`` file), and rates are bare numbers using their natural scope — the
``scope(base) × modifiers`` richness (SPEC §5) stays in the Python API for now. The remaining
levels (sequences, coupling) land here as they are rebuilt; until then their old commands live
only in ``legacy/``.
"""
from __future__ import annotations
