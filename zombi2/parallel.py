"""Replicate-level parallelism: run many independent simulations across CPU cores.

A single simulation is a sequential Markov process and isn't parallelised here; instead we
exploit the fact that *replicate* simulations (independent seeds) are embarrassingly
parallel. Each worker process runs one full simulation and **writes its own output to
disk**, returning only a small summary — so the millions of log records never cross the
process boundary.

Because the default process start method is "spawn" (macOS/Windows), every argument must be
picklable: use the built-in models/distributions or ``functools.partial`` for a
``genome_factory``, not lambdas. And, as always with multiprocessing, call
:func:`run_replicates` from within an ``if __name__ == "__main__":`` block.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .simulation import simulate_genomes
from .species_sim import simulate_species_tree


def _run_replicate(cfg: dict) -> dict:
    """Worker: run one replicate (species tree + gene families) and write it to disk."""
    rng = np.random.default_rng(cfg["seed"])  # one stream threaded through both steps
    tree = simulate_species_tree(
        cfg["species_model"], n_tips=cfg["n_tips"], age=cfg["age"],
        age_type=cfg["age_type"], rng=rng,
    )
    kwargs = {"initial_families": cfg["initial_size"], "rng": rng}
    if cfg["transfers"] is not None:
        kwargs["transfers"] = cfg["transfers"]
    if cfg["max_family_size"] is not None:
        kwargs["max_family_size"] = cfg["max_family_size"]
    if cfg["genome_factory"] is not None:
        kwargs["genome_factory"] = cfg["genome_factory"]
    if cfg["rates"] is not None:
        genomes = simulate_genomes(tree, cfg["rates"], **kwargs)
    else:
        genomes = simulate_genomes(
            tree, duplication=cfg["duplication"], transfer=cfg["transfer"],
            loss=cfg["loss"], origination=cfg["origination"], **kwargs,
        )

    path = os.path.join(cfg["outdir"], cfg["label"])
    genomes.write(path)
    return {
        "replicate": cfg["replicate"],
        "seed": cfg["seed"],
        "path": path,
        "n_species": len(genomes.profiles.species),
        "n_families": len(genomes.profiles.families),
        "n_events": len(genomes.event_log),
    }


def run_replicates(
    n_replicates: int,
    outdir: str,
    species_model,
    *,
    n_tips: int,
    age: float,
    age_type: str = "crown",
    rates=None,
    duplication: float = 0.0,
    transfer: float = 0.0,
    loss: float = 0.0,
    origination: float = 0.0,
    initial_families: int = 20,
    transfers=None,
    max_family_size=None,
    genome_factory=None,
    seed: int = 0,
    processes: int | None = None,
) -> list[dict]:
    """Run ``n_replicates`` independent simulations in parallel, writing each to disk.

    ``initial_families`` is the number of gene families seeded at the root of each
    replicate (default 20).

    Each replicate ``i`` is written to ``outdir/replicate_<i>/`` (the full ZOMBI-1-style
    output) and gets an independent seed derived from ``seed``, so the whole batch is
    reproducible and independent of the number of processes. The species-tree and
    gene-family parameters mirror :func:`~zombi2.simulate_species_tree` and
    :func:`~zombi2.simulate_genomes`.

    Parameters
    ----------
    processes:
        Number of worker processes (default: all cores). ``processes=1`` runs serially in
        the current process (handy for debugging).

    Returns
    -------
    list[dict]
        One summary per replicate (``replicate``, ``seed``, ``path``, ``n_species``,
        ``n_families``, ``n_events``), ordered by replicate index.
    """
    os.makedirs(outdir, exist_ok=True)
    child_seeds = [int(s) for s in np.random.SeedSequence(seed).generate_state(n_replicates)]
    width = max(4, len(str(n_replicates - 1)))
    configs = [
        {
            "replicate": i, "seed": cs, "label": f"replicate_{i:0{width}d}",
            "outdir": outdir, "species_model": species_model,
            "n_tips": n_tips, "age": age, "age_type": age_type,
            "rates": rates, "duplication": duplication, "transfer": transfer,
            "loss": loss, "origination": origination, "initial_size": initial_families,
            "transfers": transfers, "max_family_size": max_family_size,
            "genome_factory": genome_factory,
        }
        for i, cs in enumerate(child_seeds)
    ]

    if processes == 1:  # serial path (no multiprocessing) — deterministic and debuggable
        return [_run_replicate(c) for c in configs]
    with ProcessPoolExecutor(max_workers=processes) as executor:
        return list(executor.map(_run_replicate, configs))
