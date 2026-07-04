"""Read/write benchmark results as self-describing JSON.

One file per benchmark under ``results/``. A file is a header (the environment
snapshot + the benchmark's own metadata) plus the flat list of measured points.
Plotting reads these back and never touches the simulator, so restyling a figure
is free and re-running a benchmark only rewrites its own file.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from .environment import describe
from .timing import Point

# Repo-relative anchors so scripts work from any cwd.
ROOT = Path(__file__).resolve().parent.parent          # .../analyses/performance
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"


@dataclass
class Result:
    """A benchmark's full output: identity, environment, metadata, and points."""

    name: str
    title: str
    x_label: str
    points: list[Point]
    meta: dict = field(default_factory=dict)
    env: dict = field(default_factory=dict)

    # --- persistence ------------------------------------------------------
    def save(self, results_dir: Path | None = None) -> Path:
        results_dir = results_dir or RESULTS_DIR
        results_dir.mkdir(parents=True, exist_ok=True)
        path = results_dir / f"{self.name}.json"
        payload = {
            "name": self.name,
            "title": self.title,
            "x_label": self.x_label,
            "meta": self.meta,
            "env": self.env or describe(),
            "points": [p.to_dict() for p in self.points],
        }
        path.write_text(json.dumps(payload, indent=2) + "\n")
        return path

    @classmethod
    def load(cls, path: Path) -> "Result":
        d = json.loads(Path(path).read_text())
        return cls(
            name=d["name"],
            title=d.get("title", d["name"]),
            x_label=d.get("x_label", "x"),
            points=[Point.from_dict(p) for p in d["points"]],
            meta=d.get("meta", {}),
            env=d.get("env", {}),
        )

    # --- convenience ------------------------------------------------------
    def series(self) -> list[str]:
        """Distinct series labels, in first-seen order."""
        seen: dict[str, None] = {}
        for p in self.points:
            seen.setdefault(p.series, None)
        return list(seen)

    def by_series(self, series: str) -> list[Point]:
        pts = [p for p in self.points if p.series == series]
        return sorted(pts, key=lambda p: p.x)


def load_all(results_dir: Path | None = None) -> dict[str, Result]:
    """Load every ``*.json`` result file, keyed by benchmark name."""
    results_dir = results_dir or RESULTS_DIR
    out: dict[str, Result] = {}
    if not results_dir.exists():
        return out
    for path in sorted(results_dir.glob("*.json")):
        r = Result.load(path)
        out[r.name] = r
    return out
