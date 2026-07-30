"""``zombi2 report`` — (re)write and show ``run.zombi2``, the one-page report for a run.

Every level command writes it on the way out, so a normal pipeline never needs this. It is here for
the times that flow is broken (a run made before the report existed, a directory assembled by hand, or
one whose report you deleted), and to *show* the report on demand. It only reads the level records
already in the directory and rewrites the report from them — it simulates nothing.

``--check`` gates a pipeline on the run agreeing with itself: it exits non-zero when a downstream level
was computed on an upstream one that has since changed, so a script or CI job can trust ``$?``.
"""
from __future__ import annotations

import argparse
import sys

from zombi2._runtime.report import stale_warnings, write_run_report
from zombi2.cli.framework import _add_run_arg


def _add_report_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "the finished run to report on")
    p.add_argument("--check", action="store_true",
                   help="exit non-zero if the run is stale — a downstream level computed on an upstream "
                        "one that has since changed — so a pipeline or CI job can gate on $?. The report "
                        "is still written and printed; the warnings are already in it.")


def run(args, parser):
    path = write_run_report(args.run)
    if path is None:
        raise FileNotFoundError(
            f"{args.run}/ has no level records to report on — run a level there first. (A --flat run "
            f"writes its records into the one directory rather than per-level sub-directories, so the "
            f"report cannot tell its levels apart and none is written.)")
    with open(path, encoding="utf-8") as f:
        sys.stdout.write(f.read())              # show the report; the path goes to stderr so stdout is it
    print(f"wrote {path}", file=sys.stderr)
    if args.check and stale_warnings(args.run):
        return 1                                # the warnings are in the report already; signal via $?
    return 0
