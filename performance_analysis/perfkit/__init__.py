"""perfkit — the reusable core of the ZOMBI2 performance suite.

Four small, stable pieces, deliberately decoupled so the workspace survives
heavy iteration on the benchmarks and figures themselves:

    timing       measure(fn) -> raw per-repeat times; the Point record
    environment  describe()  -> git/interpreter/Rust snapshot for provenance
    io           Result      -> save/load self-describing JSON under results/
    style        apply()/save() -> the matplotlib publication house style

Benchmarks live in ``benchmarks.py``; measurement (``run.py``) and plotting
(``plot.py``) never import each other — they meet only through ``results/*.json``.
"""

from __future__ import annotations

from .timing import Point, measure
from .environment import describe, one_line
from .io import Result, load_all, ROOT, RESULTS_DIR, FIGURES_DIR
from . import style

__all__ = [
    "Point", "measure",
    "describe", "one_line",
    "Result", "load_all", "ROOT", "RESULTS_DIR", "FIGURES_DIR",
    "style",
]
