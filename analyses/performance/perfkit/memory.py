"""Peak-memory measurement by subprocess isolation.

``resource.getrusage`` reports a *high-water mark* for the whole process, so
measuring several sizes in one interpreter would only ever report the largest.
Instead we run each (task, size) in a fresh subprocess and read its peak RSS back
— the only way to get a clean per-size memory curve.

Run as a module for the child role::

    python -m perfkit.memory <task> <n>    # prints "<rss_mb>\t<seconds>\t<json>"
"""

from __future__ import annotations

import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent  # .../analyses/performance


def _peak_rss_mb() -> float:
    """This process's peak resident set size in MB (bytes on macOS, KB on Linux)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1e6 if platform.system() == "Darwin" else rss / 1e3


def measure(task: str, n: int, *, timeout: float = 1200.0) -> dict | None:
    """Run one (task, size) in a fresh subprocess; return its peak RSS + timing.

    Returns ``{"n", "task", "rss_mb", "seconds", "info"}`` or ``None`` if the
    child failed (e.g. ran out of address space allocating a dense N×N matrix).
    """
    cmd = [sys.executable, "-m", "perfkit.memory", task, str(n)]
    try:
        out = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True,
                             timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    if out.returncode != 0:
        return None
    try:
        rss, secs, info = out.stdout.strip().split("\t", 2)
        return dict(n=n, task=task, rss_mb=float(rss), seconds=float(secs),
                    info=json.loads(info))
    except Exception:
        return None


def _child(task: str, n: int) -> None:
    """Child role: run the task, print peak RSS + wall time + info as a TSV line."""
    import zombi2 as z
    import config

    model, rates = config.model(), config.rate_model()
    t0 = time.perf_counter()
    info: dict = {}
    if task == "tree":
        tree = z.simulate_species_tree(model, n_tips=n, age=config.TREE_AGE, seed=1)
        info["nodes"] = len(tree.leaves())
    elif task == "profiles":
        tree = z.simulate_species_tree(model, n_tips=n, age=config.TREE_AGE, seed=1)
        pm = z.simulate_genomes(tree, rates=rates, initial_families=config.INITIAL_SIZE,
                                output="profiles", seed=3)
        info["families"] = len(pm.families)
    elif task == "genomes":
        tree = z.simulate_species_tree(model, n_tips=n, age=config.TREE_AGE, seed=1)
        g = z.simulate_genomes(tree, rates=rates, initial_families=config.INITIAL_SIZE, seed=3)
        info["families"] = len(g.profiles.families)
        info["events"] = sum(1 for _ in g.event_log)
    elif task == "trace":
        tree = z.simulate_species_tree(model, n_tips=n, age=config.TREE_AGE, seed=1)
        tr = z.simulate_genomes(tree, rates=rates, initial_families=config.INITIAL_SIZE,
                                output="trace", seed=3)
        info["families"] = len(tr.profiles.families)
    else:  # pragma: no cover
        raise SystemExit(f"unknown task {task!r}")
    secs = time.perf_counter() - t0
    print(f"{_peak_rss_mb():.1f}\t{secs:.4f}\t{json.dumps(info)}")


if __name__ == "__main__":
    _child(sys.argv[1], int(sys.argv[2]))
