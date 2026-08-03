"""Shared scaffolding for the opt-in parallel execution paths.

ZOMBI2 runs single-threaded unless you explicitly ask for parallelism (SPEC: serial by default).
When you do, each independent unit — a gene tree at the sequences level, a gene family at the
family-genomes level — is evolved under its **own** spawned RNG stream, so the result is
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

import contextlib
import multiprocessing
import os
import sys
from concurrent.futures.process import BrokenProcessPool
from typing import TYPE_CHECKING, Iterator

# NB: gene-tree types are imported *lazily* inside rebuild_gene_tree, not here. This module is the
# low-level shared scaffolding; importing zombi2.genomes at top level would make genomes → _perfamily →
# _parallel → genomes a cycle. The annotations below are strings (``from __future__ import annotations``);
# this block gives the names to a type-checker without a runtime import.
if TYPE_CHECKING:
    from ..genomes.gene_trees import GeneNode, GeneTree  # noqa: F401


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


def _pool_would_fail_to_start() -> bool:
    """True when a process pool cannot start here because the workers would have to re-import a
    ``__main__`` that has no importable file — an interactive shell, ``python -c``, a stdin heredoc, or
    a Jupyter / IPython kernel. Only the re-importing start methods (``spawn`` / ``forkserver``) are
    affected; ``fork`` inherits the parent, so it is always safe. In the failing case a
    ``ProcessPoolExecutor`` dies at startup with a raw ``BrokenProcessPool``, so a caller checks this
    first and runs single-process instead. A real ``.py`` script, the ``zombi2`` console script, and a
    pytest run all have a ``__main__`` with a real file, so none of them are affected."""
    try:
        method = multiprocessing.get_start_method()
    except Exception:                       # a locked-down environment with no start method: play safe
        method = None
    if method == "fork":
        return False
    path = getattr(sys.modules.get("__main__"), "__file__", None)
    return not (path and os.path.exists(path))


def _in_a_worker() -> bool:
    """True when this process is a multiprocessing child rather than the program the user started.

    Read off ``current_process().name``, which is ``"MainProcess"`` in the parent and
    ``"SpawnProcess-N"`` / ``"ForkProcess-N"`` in a child — and, unlike `multiprocessing.parent_process`,
    is already set while the child re-imports the main module, which is exactly the moment that
    matters here."""
    return multiprocessing.current_process().name != "MainProcess"


def guard_pool_workers(workers: int, *, what: str = "--parallel") -> int:
    """Return ``workers`` unchanged, unless a process pool would crash at startup here (see
    `_pool_would_fail_to_start()`) — then fall back to ``1`` (single-process), warning once. This is
    what lets ``parallel=`` be called from a notebook, ``python -c``, or a stdin heredoc without a raw
    ``BrokenProcessPool``: it runs single-process instead of dying. A ``.py`` script (or the CLI) is
    unaffected and keeps every worker. ``what`` names the knob in the message."""
    if workers <= 1 or not _pool_would_fail_to_start():
        return workers
    # Two very different situations reach here, and they need opposite answers.
    if _in_a_worker():
        # A worker is re-importing the caller's script and running its top level again — the script
        # has no `__main__` guard. Left alone this is worse than a crash: every child silently
        # repeated the whole simulation and the run "succeeded", having done N× the work and printed
        # N copies of its output. Refusing kills the child, and the parent turns the resulting
        # BrokenProcessPool into `pool_errors`'s message, which names the guard.
        raise RuntimeError(
            f"{what} started worker processes and one of them has re-run your program from the top, "
            f"so the call is not under `if __name__ == \"__main__\":`. Add the guard.")
    # Genuinely no importable `__main__`: a notebook, `python -c`, a heredoc. Nothing is wrong with
    # the caller's code — there is simply no file for a worker to import — so run single-process.
    print(f"note: {what} needs worker processes, but this looks like an interactive session, a "
          f"notebook, or 'python -c', where they cannot re-import your program — a process pool "
          f"would crash there. Running single-process instead; run from a .py script to use all cores.")
    return 1


@contextlib.contextmanager
def pool_errors(what: str = "parallel=") -> Iterator[None]:
    """Turn a dead worker pool into an error that names the usual cause.

    `guard_pool_workers()` catches the interactive cases up front, but the commonest failure is one
    it cannot see: a plain ``.py`` script whose simulate call sits at the top level, with no
    ``if __name__ == "__main__":``. Under the ``spawn`` start method each worker re-imports the
    script, runs it from the top, tries to open a pool of its own inside a child, and dies — and all
    the parent sees is ``BrokenProcessPool``, sometimes trailed by a `FileNotFoundError` from a
    half-created stream directory. Nothing about that names the guard, so a user's first parallel
    script fails and tells them nothing. Detecting the missing guard beforehand would mean reading
    the caller's source; saying what happened afterwards costs one ``except``."""
    try:
        yield
    except BrokenProcessPool as exc:
        raise RuntimeError(
            f"{what} needs worker processes, and they died as they started. The usual cause is a "
            f"script whose simulate call is at the top level: a worker re-imports your program, runs "
            f"it again from the top, and crashes. Put the call under `if __name__ == \"__main__\":` "
            f"(the standard multiprocessing requirement — the zombi2 CLI already does). If it is "
            f"already guarded, the workers hit something else: run once with the default "
            f"{what}False to see the real error."
        ) from exc


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
    """Invert `flatten_gene_tree()` — rebuild the `GeneTree` in the worker, iteratively.
    Children are appended in flat-list order, which is the order `flatten_gene_tree()` emitted
    them, so the reconstructed tree is identical node-for-node (and sibling order is preserved)."""
    from ..genomes.gene_trees import GeneNode, GeneTree     # lazy: avoids a genomes ↔ _parallel cycle

    family, origination, nodes = flat
    built = [GeneNode(kind, species, time, copy) for (_p, kind, species, time, copy) in nodes]
    for i, (parent, *_rest) in enumerate(nodes):
        if parent >= 0:
            built[parent].children.append(built[i])
    return GeneTree(family, built[0], origination)
