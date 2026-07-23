"""``zombi2.tools`` — analyses that read a finished run and derive a new view of it.

The levels *simulate*; the tools *read back* what a run wrote and turn it into something else. The
first tool is homology classification: the true ortholog / paralog / xenolog relation of every gene
pair, straight from each gene tree's recorded events (:mod:`.homology`)."""
from __future__ import annotations

from zombi2.tools.homology import homology_table, homology_tsv, write_homology

__all__ = ["homology_table", "homology_tsv", "write_homology"]
