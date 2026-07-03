#!/usr/bin/env python
"""Run ZOMBI2 benchmarks and write raw timings to ``results/``.

    python run.py                     # core benchmarks, standard profile
    python run.py --quick             # fast smoke test (small grids, few reps)
    python run.py --full              # push into the large-tree regime
    python run.py gene_families       # just one benchmark
    python run.py --list              # show available benchmarks

Measurement only — no plotting. Re-run a single benchmark and only its result
file changes; ``plot.py`` reads whatever is in ``results/``. Kept under a
``__main__`` guard because parallel benchmarks spawn worker processes.
"""

from __future__ import annotations

import argparse
import sys
import time

from perfkit import RESULTS_DIR, describe, one_line
from benchmarks import REGISTRY, PROFILES

# The default set = the four overview panels: two scaling curves, the memory
# footprint, and parallel scaling. `write_output` is opt-in (name it explicitly).
DEFAULT = ["species_tree", "gene_families", "memory_scaling", "parallel_scaling"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*", help="benchmarks to run (default: the core set)")
    ap.add_argument("--profile", choices=list(PROFILES), default="standard")
    ap.add_argument("--quick", action="store_const", const="quick", dest="profile")
    ap.add_argument("--full", action="store_const", const="full", dest="profile")
    ap.add_argument("--all", action="store_true", help="run every registered benchmark")
    ap.add_argument("--list", action="store_true", help="list benchmarks and exit")
    args = ap.parse_args(argv)

    if args.list:
        for name, fn in REGISTRY.items():
            print(f"  {name:22s} {(fn.__doc__ or '').strip().splitlines()[0]}")
        return 0

    if args.all:
        names = list(REGISTRY)
    elif args.names:
        names = args.names
    else:
        names = DEFAULT

    unknown = [n for n in names if n not in REGISTRY]
    if unknown:
        ap.error(f"unknown benchmark(s): {', '.join(unknown)}. "
                 f"Available: {', '.join(REGISTRY)}")

    env = describe()
    print(f"Environment: {one_line(env)}")
    print(f"Profile: {args.profile}   Benchmarks: {', '.join(names)}\n")

    for name in names:
        print(f"[{name}]")
        t0 = time.perf_counter()
        result = REGISTRY[name](args.profile)
        result.env = env
        path = result.save(RESULTS_DIR)
        print(f"  -> {path.relative_to(RESULTS_DIR.parent)}  "
              f"({time.perf_counter() - t0:.1f}s)\n")

    print("Done. Render figures with:  python plot.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
