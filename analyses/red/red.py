"""RED analysis helpers.

The RED estimator now **ships in the clean core** — `zombi2.tree.relative_evolutionary_divergence`
(exact relative age on an ultrametric tree; also on the CLI as `zombi2 tools tree --red`). This
module used to carry a local port of it; that is gone. What stays is the one analysis-side helper:
the nodes RED is graded on.
"""
from __future__ import annotations

from zombi2.tree import Tree


def internal_nodes(tree: Tree) -> list[int]:
    """Ids of the nodes RED is graded on: internal, non-root (leaves are trivially 1, the root 0)."""
    return [i for i, nd in tree.nodes.items()
            if nd.children is not None and nd.parent is not None]
