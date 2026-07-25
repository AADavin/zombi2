"""Auxiliary runtime helpers shared across the levels — cross-cutting plumbing, not domain code.

``progress`` is the progress-bar façade; ``parallel`` is the shared scaffolding for the opt-in
parallel engines (worker resolution, gene-tree pickling). The per-level parallel engines themselves
live with their level (``genomes/_perfamily``, ``sequences/_pergenetree``)."""
