"""Species-tree namespace (scikit-learn-style).

Re-exports the species-tree public API so users can write::

    from zombi2.species import BirthDeath, simulate_species_tree

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.tree`, :mod:`zombi2.species.model`,
:mod:`zombi2.species.sim`, :mod:`zombi2.species.forward`,
:mod:`zombi2.species.ghosts`); it does not redefine anything.
"""

from __future__ import annotations

from zombi2.tree import Tree, TreeNode, read_newick, prune
from zombi2.species.model import (
    BirthDeath, Yule, EpisodicBirthDeath, ClaDS, DiversityDependent,
    SharedBirthDeath,  # noqa: F401 (deprecated preset for BirthDeath(per="shared"); not in __all__)
    CladeShiftBirthDeath,
)
from zombi2.species.sim import simulate_species_tree
from zombi2.species.ghosts import add_ghost_lineages

__all__ = [
    "Tree", "TreeNode", "read_newick", "prune",
    "BirthDeath", "Yule", "EpisodicBirthDeath", "ClaDS", "DiversityDependent",
    "CladeShiftBirthDeath", "simulate_species_tree", "add_ghost_lineages",
]
