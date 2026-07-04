"""Parser for eggNOG-mapper output -> a ZOMBI2 phylogenetic profile + COG categories.

eggNOG-mapper annotations are one row per gene, with an ``eggNOG_OGs`` column (comma-separated
``OGid@taxlevel|taxname`` assignments, broadest first) and a ``COG_category`` letter. We collapse
genes to **orthologous groups** (the broadest OG, i.e. the COG/ENOG at root level), count copies per
genome, and return a :class:`zombi2.ProfileMatrix` (OGs x genomes) plus each OG's COG functional
category. Genome ids (``project_id``) are optionally mapped to short names (tree leaf labels) via a
genomes table.

    from eggnog import load_eggnog
    pm, category = load_eggnog("eggnog_annotations.tsv", genomes_tsv="genomes.tsv")
    # pm : ProfileMatrix(families=OGs, species=genome short codes)
    # category : {OG -> COG letter}   (the primary category, most frequent across the OG's genes)
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict

import numpy as np

from zombi2 import ProfileMatrix


def _short_map(genomes_tsv, id_col, short_col):
    with open(genomes_tsv) as f:
        return {r[id_col]: r[short_col] for r in csv.DictReader(f, delimiter="\t")}


def load_eggnog(annotations_tsv, *, genomes_tsv=None, id_col="project_id",
                short_col="ShortCode", og_field="eggNOG_OGs", cat_field="COG_category"):
    """Load an eggNOG-mapper annotation table into a profile + COG-category map.

    Parameters
    ----------
    annotations_tsv:
        Path to the eggNOG-mapper ``*.annotations`` TSV (one row per gene).
    genomes_tsv, id_col, short_col:
        Optional genomes table mapping the annotation's ``id_col`` (default ``project_id``) to a
        short genome name in ``short_col`` (default ``ShortCode`` -- the tree leaf labels). Without
        it, the raw ``project_id`` is used as the genome name.
    og_field, cat_field:
        Column names for the OG assignments and the COG category.

    Returns
    -------
    (ProfileMatrix, dict)
        The copy-number profile (rows = orthologous groups, columns = genomes) and a
        ``{OG -> COG letter}`` map (the primary category = the most frequent letter across that
        OG's genes; ``"?"`` if unannotated).
    """
    short = _short_map(genomes_tsv, id_col, short_col) if genomes_tsv else None
    counts = defaultdict(int)                 # (genome, OG) -> copy number
    og_cats = defaultdict(Counter)            # OG -> Counter of COG letters
    genomes = set()
    with open(annotations_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = r[id_col]
            genome = short.get(gid) if short is not None else gid
            if genome is None:
                continue
            ogs = r.get(og_field)
            if not ogs or ogs == "-":
                continue
            og = ogs.split(",")[0].split("@")[0]        # broadest OG (root-level COG/ENOG)
            counts[(genome, og)] += 1
            genomes.add(genome)
            cat = r.get(cat_field)
            if cat and cat != "-":
                og_cats[og][cat[0]] += 1                 # primary (first) category letter
    genomes = sorted(genomes)
    ogs = sorted({og for _, og in counts})
    matrix = np.array([[counts[(g, og)] for g in genomes] for og in ogs], dtype=int)
    category = {og: (og_cats[og].most_common(1)[0][0] if og_cats[og] else "?") for og in ogs}
    return ProfileMatrix(families=ogs, species=genomes, matrix=matrix), category
