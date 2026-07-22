"""Capture the machine + code state a benchmark ran on.

Every result file embeds this so a figure is never ambiguous about *what* was
measured or *when*: a timing curve is only meaningful next to the git commit,
the interpreter, and whether the Rust extension was actually loaded.
"""

from __future__ import annotations

import platform
import subprocess
import sys
from datetime import datetime, timezone


def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(["git", *args], capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None if out.returncode == 0 else None
    except Exception:
        return None


def describe() -> dict:
    """A JSON-serialisable snapshot of the environment for the result header."""
    import numpy

    try:
        import zombi2
        zombi2_version = getattr(zombi2, "__version__", "?")
        rust = bool(zombi2.rust_available())
    except Exception as exc:  # pragma: no cover - defensive
        zombi2_version = f"import-failed: {exc}"
        rust = False

    commit = _git("rev-parse", "--short", "HEAD")
    dirty = _git("status", "--porcelain")
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": commit,
        "git_dirty": bool(dirty),
        "zombi2_version": zombi2_version,
        "rust_available": rust,
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "platform": platform.platform(),
        "processor": platform.processor() or platform.machine(),
        "cpu_count": _cpu_count(),
    }


def _cpu_count() -> int | None:
    import os

    return os.cpu_count()


def one_line(env: dict) -> str:
    """A compact caption string for stamping onto figures."""
    bits = [f"zombi2 {env.get('zombi2_version', '?')}"]
    if env.get("git_commit"):
        bits.append(f"@{env['git_commit']}{'*' if env.get('git_dirty') else ''}")
    bits.append("Rust" if env.get("rust_available") else "pure-Python")
    proc = env.get("processor") or ""
    if proc:
        bits.append(proc)
    bits.append(f"Python {env.get('python', '?')}")
    date = (env.get("timestamp_utc") or "")[:10]
    if date:
        bits.append(date)
    return " · ".join(bits)
