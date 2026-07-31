"""The parallel engine for `simulate_sequences()` — one gene tree per process.

Reached only when ``parallel`` is truthy; the serial default lives in the package ``__init__``. Each
family is evolved under its own RNG stream (spawned from the run seed in the parent), so the result is
identical for any worker count — the process pool is purely an execution detail, and a run too small to
be worth spawning workers is evaluated inline with the *same* streams, giving the *same* bytes.

Why processes and not threads: measured on this codebase, a thread pool barely helps and often hurts —
numpy releases the GIL too little for the per-site arrays here. So each worker is a real process, the
shared read-only inputs (the run's partitions or its two nucleotide models, and the lineage clock) are
shipped once via an initializer, and the gene tree crosses the boundary in the flat, recursion-free form
of `zombi2._runtime.parallel` (a deep tree overflows the pickle recursion limit otherwise).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .._runtime.parallel import flatten_gene_tree, rebuild_gene_tree
from .._runtime.progress import progress_bar

#: Below this many gene trees the process-pool spawn + IPC costs more than it saves (the measured
#: floor is ~0.2 s), so the parallel engine evaluates inline instead. Same streams, same output.
_MIN_FAMILIES_FOR_POOL = 2

# Per-worker state, set once by the initializer so the models and clock are not re-pickled per task.
# In-process (inline) runs set it directly. The CDF caches give each worker the same run-wide cache
# the serial engine keeps (branch lengths recur massively across gene trees), keyed the same way —
# by model identity, so two matrices can never share one cache of branch-length-keyed matrices.
_MODELS: tuple = ()
_CLOCK = None
_CACHES: dict = {}
_NAMES: dict = {}   # the run's node names — e<id> for a lineage that died
#: the run-wide ``((model, sites), …)`` a family or ordered run evolves, or ``None`` on a nucleotide
#: run, where each block's model and length come from its own task instead
_PARTS: "tuple | None" = None


def _init_worker(models, clock, names, partitions=None) -> None:
    global _MODELS, _CLOCK, _CACHES, _NAMES, _PARTS
    _MODELS, _CLOCK, _NAMES, _PARTS = models, clock, names, partitions
    _CACHES = {}


def _evolve_one(task):
    """Evolve one family's gene tree and return its labelled outputs. Runs in a worker process (or
    inline); reads the shared partitions/models/clock from module state, keeps a per-worker CDF cache
    per model."""
    from . import _evolve_partitions, _gene_newick, _scaled_gene_tree   # lazy: package, no cycle

    family, flat, midx, length, rate, seed_states, seedseq = task
    # A family or ordered run's partitions are run-wide and came through the initializer; a
    # nucleotide block names its model by index and brings its own length in the task.
    parts = _PARTS if _PARTS is not None else ((_MODELS[midx], length),)
    rng = np.random.default_rng(seedseq)
    gt = rebuild_gene_tree(flat)
    aln, anc, fnd = _evolve_partitions(gt, parts, rate, _CLOCK, rng, _CACHES, _NAMES,
                                       founding=seed_states)
    scaled = _scaled_gene_tree(gt, rate, _CLOCK)             # branch lengths in subs/site
    ext = scaled.extant
    phylo = {"complete": _gene_newick(scaled.complete, _NAMES),
             "extant": _gene_newick(ext, _NAMES) if ext is not None else None}
    return family, aln, anc, fnd, phylo


def evolve_families(gene_trees, per_block, model, intergene_model, length, rate_base, clock,
                    founding_seed, family_seeds, workers, progress, names, sink=None,
                    partitions=None):
    """Evolve every family concurrently and assemble the four output maps.

    ``family_seeds[i]`` is the spawned RNG stream for the *i*-th family in sorted order, so the family
    a result belongs to is fixed before any worker runs — the assignment is independent of which worker
    finishes when, which is what makes the run worker-count invariant.

    Which model evolves what is a **run-wide** property on a family or ordered run — every family
    takes the same ``partitions`` — so the partitions ride the initializer alongside the clock and
    are pickled once per worker rather than once per family. A nucleotide run is the other shape:
    genes and spacer are its only two models (``per_block`` maps each block to one of them), so there
    the task carries the model's index and the block's own length, rate multiplier and founding seed.

    ``sink``, when given, is handed each family the moment its result arrives and the four maps are
    left empty: that is the streamed run, where nothing family-sized is kept. Results come back in
    family order either way, so streaming writes the same bytes the in-memory path would."""
    families = sorted(gene_trees)
    models = () if partitions is not None else (model, intergene_model)
    tasks = []
    for i, family in enumerate(families):
        flat = flatten_gene_tree(gene_trees[family])
        if per_block is None:
            # the partitions say which models and how many sites; there is no index to give
            midx, f_len, f_rate, seed_states = None, None, rate_base, None
        else:                            # a nucleotide block: its own length, model and speed
            f_len, f_model, speed = per_block[family]
            midx = 0 if f_model is model else 1
            f_rate, seed_states = rate_base * speed, founding_seed[family]
        tasks.append((family, flat, midx, f_len, f_rate, seed_states, family_seeds[i]))

    alignments: dict[int, dict[str, str]] = {}
    ancestral: dict[int, dict[str, str]] = {}
    founding: dict[int, str] = {}
    phylograms: dict[int, dict[str, str | None]] = {}
    bar = progress_bar(len(families), "sequences", unit="family", enabled=progress)

    def _collect(results):
        for family, aln, anc, fnd, phylo in results:
            if sink is None:
                alignments[family], ancestral[family] = aln, anc
                founding[family], phylograms[family] = fnd, phylo
            else:
                sink(family, aln, anc, fnd, phylo)   # straight to disk; nothing is kept
            bar.update()

    n = len(families)
    if workers > 1 and n >= _MIN_FAMILIES_FOR_POOL:
        w = min(workers, n)
        with ProcessPoolExecutor(max_workers=w, initializer=_init_worker,
                                 initargs=(models, clock, names, partitions)) as ex:
            _collect(ex.map(_evolve_one, tasks, chunksize=max(1, n // (w * 8))))
    else:
        # inline: the same worker + per-process caches, no pool
        _init_worker(models, clock, names, partitions)
        _collect(_evolve_one(t) for t in tasks)
    bar.close()
    return alignments, ancestral, founding, phylograms
