"""The parallel engine for :func:`~zombi2.sequences.simulate_sequences` — one gene tree per process.

Reached only when ``parallel`` is truthy; the serial default lives in the package ``__init__``. Each
family is evolved under its own RNG stream (spawned from the run seed in the parent), so the result is
identical for any worker count — the process pool is purely an execution detail, and a run too small to
be worth spawning workers is evaluated inline with the *same* streams, giving the *same* bytes.

Why processes and not threads: measured on this codebase, a thread pool barely helps and often hurts —
numpy releases the GIL too little for the per-site arrays here. So each worker is a real process, the
shared read-only inputs (the models and the lineage clock) are shipped once via an initializer, and the
gene tree crosses the boundary in the flat, recursion-free form of :mod:`zombi2._parallel` (a deep tree
overflows the pickle recursion limit otherwise).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .._parallel import flatten_gene_tree, rebuild_gene_tree
from ..progress import progress_bar
from .evolution import evolve_gene_tree
from .substitution_models import decode

#: Below this many gene trees the process-pool spawn + IPC costs more than it saves (the measured
#: floor is ~0.2 s), so the parallel engine evaluates inline instead. Same streams, same output.
_MIN_FAMILIES_FOR_POOL = 2

# Per-worker state, set once by the initializer so the models and clock are not re-pickled per task.
# In-process (inline) runs set it directly. A list of CDF caches, one per model, gives each worker the
# same run-wide cache the serial engine keeps (branch lengths recur massively across gene trees).
_MODELS: tuple = ()
_CLOCK = None
_CACHES: list = []


def _init_worker(models, clock) -> None:
    global _MODELS, _CLOCK, _CACHES
    _MODELS, _CLOCK = models, clock
    _CACHES = [dict() for _ in models]


def _evolve_one(task):
    """Evolve one family's gene tree and return its labelled outputs. Runs in a worker process (or
    inline); reads the shared models/clock from module state, keeps a per-worker CDF cache per model."""
    from . import _gene_newick, _scaled_gene_tree, _split   # package helpers; imported lazily, no cycle

    family, flat, midx, length, rate, seed_states, seedseq = task
    model = _MODELS[midx]
    rng = np.random.default_rng(seedseq)
    gt = rebuild_gene_tree(flat)
    states, founding_states = evolve_gene_tree(gt.complete, model, length, rate, _CLOCK, rng,
                                               gt.origination, founding=seed_states,
                                               cdf_cache=_CACHES[midx])
    aln, anc = _split(gt, states, model)
    scaled = _scaled_gene_tree(gt, rate, _CLOCK)             # branch lengths in subs/site
    ext = scaled.extant
    phylo = {"complete": _gene_newick(scaled.complete),
             "extant": _gene_newick(ext) if ext is not None else None}
    return family, aln, anc, decode(founding_states, model.alphabet), phylo


def evolve_families(gene_trees, per_block, model, intergene_model, length, rate_base, clock,
                    founding_seed, family_seeds, workers, progress):
    """Evolve every family concurrently and assemble the four output maps.

    ``family_seeds[i]`` is the spawned RNG stream for the *i*-th family in sorted order, so the family
    a result belongs to is fixed before any worker runs — the assignment is independent of which worker
    finishes when, which is what makes the run worker-count invariant. Genes and spacer are the only two
    models a run uses (``per_block`` maps each nucleotide block to one of them); everything else — the
    per-block length, rate multiplier and founding seed — travels inside the task."""
    families = sorted(gene_trees)
    models = (model,) if per_block is None else (model, intergene_model)
    tasks = []
    for i, family in enumerate(families):
        flat = flatten_gene_tree(gene_trees[family])
        if per_block is None:
            midx, f_len, f_rate, seed_states = 0, length, rate_base, None
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
            alignments[family], ancestral[family] = aln, anc
            founding[family], phylograms[family] = fnd, phylo
            bar.update()

    n = len(families)
    if workers > 1 and n >= _MIN_FAMILIES_FOR_POOL:
        w = min(workers, n)
        with ProcessPoolExecutor(max_workers=w, initializer=_init_worker,
                                 initargs=(models, clock)) as ex:
            _collect(ex.map(_evolve_one, tasks, chunksize=max(1, n // (w * 8))))
    else:
        _init_worker(models, clock)      # inline: the same worker + per-process caches, no pool
        _collect(_evolve_one(t) for t in tasks)
    bar.close()
    return alignments, ancestral, founding, phylograms
