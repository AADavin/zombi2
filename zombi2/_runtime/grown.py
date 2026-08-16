"""What an engine hands back, in one shape.

Six engines grow a pair of levels together, and each one used to return whatever it happened to have:
a four-tuple, a five-tuple, a dict, a small record of its own. Reading two of them side by side meant
counting positions, and adding a field to one meant finding every caller.

They return this instead. Every field is optional, because every engine fills only what it grew, and
the field names are the ones the run's own result uses — so an engine's output reads like a smaller
copy of `~zombi2.joint.JointResult` rather than like a private convention.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Grown:
    """One engine's output. ``tree`` is the species tree it grew, or the one it was handed."""

    tree: Any = None
    species_events: list = field(default_factory=list)

    #: the trait level: its state or value at every node, and its log
    trait_values: dict = field(default_factory=dict)
    trait_events: list = field(default_factory=list)

    #: the genome level: each lineage's genome, the log, the declared families, and the genome the
    #: run started from — the last one a tuple of gene copies, as the genome level keeps it
    genomes: dict = field(default_factory=dict)
    genome_events: list = field(default_factory=list)
    genome_names: dict = field(default_factory=dict)
    genome_initial: tuple = ()

    #: the sequence level, per gene **name**: the states at every gene-tree node, the founding
    #: sequence, and each branch's length in substitutions per site
    sequences: dict = field(default_factory=dict)
