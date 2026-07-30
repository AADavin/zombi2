"""``zombi2 report`` — (re)write ``run.zombi2``, the one-page human-readable report for a run.

Every level command writes it on the way out, so a normal pipeline never needs this. It is here for
the times that flow is broken: a run made before the report existed, a directory assembled by hand, or
one whose report you deleted. It only reads the level records already in the directory and rewrites the
report from them — it simulates nothing.
"""
from __future__ import annotations

import argparse

from zombi2._runtime.report import write_run_report
from zombi2.cli.framework import _add_run_arg


def _add_report_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "the finished run to report on")


def run(args, parser):
    path = write_run_report(args.run)
    if path is None:
        raise FileNotFoundError(
            f"{args.run}/ holds no level records to report on — run a level there first (a --flat run "
            f"keeps no per-level records, so it has no report).")
    print(f"wrote {path}")
    return 0
