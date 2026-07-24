"""Shared scaffolding for the opt-in parallel execution paths.

ZOMBI2 runs single-threaded unless you explicitly ask for parallelism (SPEC: serial by default).
When you do, each independent unit — a gene tree at the sequences level, a gene family at the
unordered-genomes level — is evolved under its **own** spawned RNG stream, so the result is
bit-identical for any worker count (``parallel=2`` and ``parallel=8`` agree to the byte), though it
differs from the serial reference engine: parallel is a **separate** engine, chosen deliberately
(the "A" decision — the serial default stays untouched, no fixture re-blessing).

Parallelism here is **process**-backed, not thread-backed. Measured on this codebase, a thread pool
barely helps and often hurts: numpy releases the GIL too little for the array sizes these inner loops
touch. A process pool scales (≈3× on genome-scale runs) at the cost of pickling across the boundary —
which is why gene trees are shipped in a **flat, recursion-free** form here (a deep gene tree overflows
the default pickle recursion limit; the rest of the codebase is iterative for the same reason).
"""

from __future__ import annotations

import os

from .genomes.gene_trees import GeneNode, GeneTree


def resolve_workers(parallel) -> int:
    """Turn the public ``parallel`` knob into a worker count.

    ``False`` / ``None`` → ``1`` (serial, the default). ``True`` → every core
    (``os.cpu_count()``). A positive ``int`` → that many workers. ``1`` means the parallel engine
    run inline (no pool) — still the spawned-stream engine, so it matches any higher count byte for
    byte; it is the serial *reference* engine (``parallel=False``) that differs."""
    if parallel is None or parallel is False:
        return 1
    if parallel is True:
        return os.cpu_count() or 1
    if isinstance(parallel, bool) or not isinstance(parallel, int) or parallel < 1:
        raise ValueError(
            f"parallel must be False (serial), True (all cores), or a positive int (worker count), "
            f"got {parallel!r}")
    return parallel


def flatten_gene_tree(gt: GeneTree) -> tuple[int, float, list[tuple[int, str, int, float, int]]]:
    """A gene tree as a flat, picklable triple ``(family, origination, nodes)`` — no object graph,
    so it crosses a process boundary without tripping the pickle recursion limit a deep tree would.

    ``nodes`` is a pre-order list of ``(parent_index, kind, species, time, copy)``; the root is index
    0 with ``parent_index = -1``. Children keep their original left-to-right order (siblings appear in
    order), which the evolution walk and the Newick both depend on."""
    root = gt.complete
    nodes: list[tuple[int, str, int, float, int]] = []
    stack: list[tuple[GeneNode, int]] = [(root, -1)]
    while stack:
        node, parent = stack.pop()
        idx = len(nodes)
        nodes.append((parent, node.kind, node.species, node.time, node.copy))
        # push reversed so the leftmost child pops next and takes the next index — a plain pre-order
        for child in reversed(node.children):
            stack.append((child, idx))
    return gt.family, gt.origination, nodes


def rebuild_gene_tree(flat: tuple[int, float, list[tuple[int, str, int, float, int]]]) -> GeneTree:
    """Invert :func:`flatten_gene_tree` — rebuild the :class:`GeneTree` in the worker, iteratively.
    Children are appended in flat-list order, which is the order :func:`flatten_gene_tree` emitted
    them, so the reconstructed tree is identical node-for-node (and sibling order is preserved)."""
    family, origination, nodes = flat
    built = [GeneNode(kind, species, time, copy) for (_p, kind, species, time, copy) in nodes]
    for i, (parent, *_rest) in enumerate(nodes):
        if parent >= 0:
            built[parent].children.append(built[i])
    return GeneTree(family, built[0], origination)
