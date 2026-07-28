"""``zombi2.tools`` — analyses that read a finished run and derive a new view of it.

The levels *simulate*; the tools *read back* what a run wrote and turn it into something else. Both
tools here re-express one gene tree's recorded events: **homology** classification, how every gene
pair diverged and whether transfer is in its history (`homology`), and **recPhyloXML**, the gene tree
written inside the species tree in the community format for that (`recphylo`)."""
from __future__ import annotations

from zombi2.tools.homology import homology_table, homology_tsv, write_homology
from zombi2.tools.recphylo import recphylo_xml, write_recphylo

__all__ = ["homology_table", "homology_tsv", "write_homology", "recphylo_xml", "write_recphylo"]
