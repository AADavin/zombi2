"""The model-free heterogeneity target: the coefficient of variation of root-to-tip substitution
distances in the GTDB archaeal tree. Extant tips are contemporaneous, so this spread is pure
among-lineage rate variation — no dating or ultrametricising required (Parks et al. 2018; the same
tree Rinke et al. 2021 used). Port of the retired ``red_clock/measure_gtdb.py``."""
from __future__ import annotations

import pathlib
import numpy as np


def root_to_tip_depths(tree_path: str | pathlib.Path) -> np.ndarray:
    """Root-to-tip substitution distance for every tip of a Newick phylogram (branch lengths in
    substitutions/site). Uses ete3 to read GTDB's quoted node labels."""
    from ete3 import Tree
    t = Tree(str(tree_path), format=1, quoted_node_names=True)
    depths: list[float] = []
    for node in t.traverse("preorder"):
        node.add_feature("depth", 0.0 if node.up is None else node.up.depth + node.dist)
        if node.is_leaf():
            depths.append(node.depth)
    return np.asarray(depths)


def cv(depths: np.ndarray) -> float:
    """Coefficient of variation (s.d. / mean)."""
    return float(depths.std() / depths.mean())


if __name__ == "__main__":
    here = pathlib.Path(__file__).parent
    d = root_to_tip_depths(here / "data" / "ar53.tree")
    print(f"GTDB archaea: n_tips={len(d)}  root-to-tip mean={d.mean():.4f}  CV={cv(d):.4f}")
