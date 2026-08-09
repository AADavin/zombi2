"""The ``zombi2`` command line — one subcommand per level, mirroring the Python API.

Six commands, grouped as ``zombi2 -h`` groups them — four levels, the joint run, and the tools that
read a finished run:

- ``zombi2 species``    → `zombi2.species.simulate_species_tree()`
- ``zombi2 genomes``    → `zombi2.genomes.simulate_genomes_family()` /
  `simulate_genomes_ordered()` / `simulate_genomes_nucleotide()` (chosen by ``--resolution``)
- ``zombi2 sequences``  → `zombi2.sequences.simulate_sequences()`
- ``zombi2 traits``     → `zombi2.traits.simulate_continuous()` /
  `simulate_discrete()` (chosen by ``--kind``)
- ``zombi2 joint``      → `zombi2.joint.simulate_joint()`, a tree and the level driving it at once
- ``zombi2 tools``      → analyses that read a run someone has already made

Each subcommand's long options **are** the API keyword names (one word per concept across the API,
the CLI, and a ``--params`` file), and rates are bare numbers using their natural scope — a rate
written from its scope with verbs chained onto it (SPEC §5) goes in full form on the rate flags.
"""
from __future__ import annotations
