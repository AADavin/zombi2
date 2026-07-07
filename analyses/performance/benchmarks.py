"""The ZOMBI2 benchmark registry.

Each benchmark measures one *task* across a grid of problem sizes and returns a
:class:`perfkit.Result`. Benchmarks are registered with :func:`benchmark`,
discovered by ``run.py``, and never plot anything.

The library runs a single engine by default (the built-in model is Rust; no
Python-vs-Rust comparison anymore), so these measure ZOMBI2 as a user actually
runs it and push each task to the scale where it stops being cheap.

All benchmarks share one regime (``config.py``): trees are built once per size,
only the timed call sits inside :func:`perfkit.measure`, setup is untimed.
"""

from __future__ import annotations

import gc
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable

import zombi2 as z

import config
from perfkit import Point, Result, measure
from perfkit import memory as memkit

# --- registry -------------------------------------------------------------

REGISTRY: dict[str, Callable] = {}


def benchmark(name: str):
    def deco(fn: Callable) -> Callable:
        REGISTRY[name] = fn
        fn.bench_name = name
        return fn
    return deco


# --- run profiles (how big / how many repeats) ----------------------------

PROFILES = {
    "quick":    dict(reps=3, budget=6.0),
    "standard": dict(reps=5, budget=25.0),
    "full":     dict(reps=7, budget=90.0),
}


def _grid(profile: str, quick, standard, full):
    return {"quick": quick, "standard": standard, "full": full}[profile]


def _warmup(n: int) -> int:
    """One warm-up for cheap points; none for the multi-second giants."""
    return 1 if n < 100_000 else 0


def _timed_grid(profile, sizes, fn, series, *, cap=None, work_of=None):
    """Measure ``fn(n)`` over ``sizes`` (<= cap), tolerating failures at the top
    end (e.g. a dense N×N matrix that overflows address space)."""
    cfg = PROFILES[profile]
    points: list[Point] = []
    for n in sizes:
        if cap is not None and n > cap:
            continue
        counter = {"i": 0}

        def call():
            counter["i"] += 1
            return fn(n, counter["i"])

        try:
            times, result = measure(call, repeat=cfg["reps"], warmup=_warmup(n),
                                    max_seconds=cfg["budget"])
        except (MemoryError, RuntimeError, ValueError):
            print(f"  {series:26s} n={n:<9d} FAILED (skipped)")
            break
        work = work_of(result) if work_of else {}
        points.append(Point(series=series, x=n, times=times, work=work))
        print(f"  {series:26s} n={n:<9d} best={min(times)*1e3:10.2f} ms  {work}")
    return points


# =========================================================================
# 1. Species-tree simulation — scales to millions of tips
# =========================================================================

@benchmark("species_tree")
def species_tree(profile: str) -> Result:
    """Wall-clock to simulate a species tree vs number of extant tips.

    The backward reconstructed sampler and the forward complete-tree sampler,
    both O(N) after the assembly fix — pushed into the millions-of-tips regime.
    """
    sizes = _grid(profile,
                  quick=[100, 10_000, 100_000],
                  standard=[100, 1000, 10_000, 100_000, 1_000_000, 3_000_000],
                  full=[100, 1000, 10_000, 100_000, 1_000_000, 3_000_000, 10_000_000])
    fwd_cap = {"quick": 10_000, "standard": 100_000, "full": 1_000_000}[profile]

    points = []
    for direction, cap in (("backward", max(sizes)), ("forward", fwd_cap)):
        def make(direction):
            def fn(n, i):
                return z.simulate_species_tree(
                    config.model(), n_tips=n,
                    age=(config.TREE_AGE if direction == "backward" else None),
                    direction=direction, seed=1000 + i)
            return fn
        points += _timed_grid(
            profile, sizes, make(direction), f"Species tree · {direction}", cap=cap,
            work_of=lambda t: dict(n_leaves=len(t.leaves()),
                                   n_extant=sum(1 for x in t.leaves() if x.is_extant)))

    return Result(
        name="species_tree",
        title="Species-tree simulation scales to millions of tips",
        x_label="Number of extant tips",
        points=points,
        meta=dict(regime=config.label(), profile=profile),
    )


# =========================================================================
# 2. Gene-family simulation — full genomes vs counts-only, to the memory wall
# =========================================================================

@benchmark("gene_families")
def gene_families(profile: str) -> Result:
    """Gene-family (D/T/L/O) simulation vs tip count, three output modes.

    ``Rust · full genomes`` materialises the whole event log (one Python object per
    event → memory grows with the event count → caps around 10^5 tips). ``Rust · event
    trace`` keeps the genealogy as the engine's compact columns (``output="trace"``),
    deferring the per-event objects, so gene trees stay reconstructable yet it scales to
    millions like the counts path (~2x its cost). ``Rust · profiles only`` returns just
    the sparse counts matrix (the ABC / large-dataset fast path) — the cheapest, no
    genealogy at all.
    """
    sizes = _grid(profile,
                  quick=[1000, 10_000, 100_000],
                  standard=[1000, 3000, 10_000, 30_000, 100_000, 300_000, 1_000_000],
                  full=[1000, 10_000, 100_000, 1_000_000, 3_000_000])
    # Full event log materialises one object per event → caps around 10^5 tips.
    full_cap = {"quick": 10_000, "standard": 100_000, "full": 100_000}[profile]
    # Event trace keeps the genealogy as compact columns (no per-event objects) → reaches
    # the millions like the counts path, just carrying the full event columns (~2x cost).
    trace_cap = {"quick": 100_000, "standard": 1_000_000, "full": 3_000_000}[profile]
    # Counts-only output is SPARSE (COO) → O(N), so it scales to millions of tips.
    prof_cap = {"quick": 100_000, "standard": 1_000_000, "full": 3_000_000}[profile]

    trees = {n: z.simulate_species_tree(config.model(), n_tips=n, age=config.TREE_AGE,
                                        seed=100 + n) for n in sizes}
    rates = config.rate_model()

    def genomes_fn(n, i):
        return z.simulate_genomes(trees[n], rates=rates,
                                  initial_families=config.INITIAL_SIZE, seed=7000 + i)

    def trace_fn(n, i):
        return z.simulate_genomes(trees[n], rates=rates, initial_families=config.INITIAL_SIZE,
                                  output="trace", seed=7000 + i)

    def profiles_fn(n, i):
        return z.simulate_genomes(trees[n], rates=rates, initial_families=config.INITIAL_SIZE,
                                  output="profiles", seed=7000 + i)

    points = _timed_grid(profile, sizes, genomes_fn, "Rust · full genomes",
                         cap=full_cap, work_of=_genomes_work)
    points += _timed_grid(profile, sizes, trace_fn, "Rust · event trace",
                          cap=trace_cap, work_of=lambda tr: dict(n_families=len(tr.profiles.families)))
    points += _timed_grid(profile, sizes, profiles_fn, "Rust · profiles only",
                          cap=prof_cap, work_of=lambda pm: dict(n_families=len(pm.families)))

    return Result(
        name="gene_families",
        title="Gene-family simulation: full log vs event trace vs counts-only",
        x_label="Number of extant tips",
        points=points,
        meta=dict(regime=config.label(), profile=profile),
    )


def _genomes_work(g) -> dict:
    try:
        n_ev = sum(1 for _ in g.event_log)
    except Exception:
        n_ev = None
    return dict(n_families=len(g.profiles.families), n_events=n_ev)


# =========================================================================
# 3. Memory footprint — how big a run fits in RAM
# =========================================================================

@benchmark("memory_scaling")
def memory_scaling(profile: str) -> Result:
    """Peak resident memory vs tip count, measured in isolated subprocesses.

    The species tree is featherweight (linear, sub-GB per million tips) so it
    scales to tens of millions. The gene-family outputs cost more per tip: the
    full log grows with events, and the dense counts matrix is the ceiling.
    """
    tree_sizes = _grid(profile,
                       quick=[10_000, 100_000],
                       standard=[10_000, 100_000, 1_000_000, 3_000_000],
                       full=[10_000, 100_000, 1_000_000, 3_000_000, 10_000_000])
    # Full event log materialises one object per event → caps ~10^5 tips.
    full_sizes = _grid(profile,
                       quick=[10_000],
                       standard=[10_000, 30_000, 100_000],
                       full=[10_000, 30_000, 100_000])
    # Event trace keeps compact columns (no per-event objects) → far lighter than the
    # full log, reaching the millions with the tree.
    trace_sizes = _grid(profile,
                        quick=[10_000, 100_000],
                        standard=[10_000, 100_000, 1_000_000],
                        full=[10_000, 100_000, 1_000_000, 3_000_000])
    # Sparse counts output → O(N), scales with the tree into the millions.
    prof_sizes = _grid(profile,
                       quick=[10_000, 100_000],
                       standard=[10_000, 100_000, 1_000_000],
                       full=[10_000, 100_000, 1_000_000, 3_000_000])

    plan = ([("tree", n, "Species tree") for n in tree_sizes]
            + [("genomes", n, "Gene families · full") for n in full_sizes]
            + [("trace", n, "Gene families · trace") for n in trace_sizes]
            + [("profiles", n, "Gene families · profiles") for n in prof_sizes])

    points: list[Point] = []
    for task, n, series in plan:
        m = memkit.measure(task, n)
        if m is None:
            print(f"  {series:26s} n={n:<9d} FAILED (out of memory / timeout)")
            continue
        # store RSS in the "times" slot so the same Point machinery carries it;
        # the plotter reads work['rss_mb'] for the y-axis.
        points.append(Point(series=series, x=n, times=[m["seconds"]],
                            work=dict(rss_mb=m["rss_mb"], **m["info"])))
        print(f"  {series:26s} n={n:<9d} peakRSS={m['rss_mb']:8.0f} MB  "
              f"({m['seconds']:.1f}s)  {m['info']}")

    return Result(
        name="memory_scaling",
        title="Peak memory footprint vs run size",
        x_label="Number of extant tips",
        points=points,
        meta=dict(regime=config.label(), profile=profile,
                  note="peak RSS from resource.getrusage in an isolated subprocess"),
    )


# =========================================================================
# 4. Parallel scaling — embarrassingly-parallel replicates (compute-bound)
# =========================================================================

def _par_work(args: tuple[int, int]) -> int:
    """Module-level worker (picklable for spawn): one full simulation.

    The tip count travels *in the argument*, not a module global — a spawned
    worker re-imports this module fresh and would not see a runtime-set global.
    """
    seed, n_tips = args
    tree = z.simulate_species_tree(config.model(), n_tips=n_tips,
                                   age=config.TREE_AGE, seed=seed)
    g = z.simulate_genomes(tree, rates=config.rate_model(),
                           initial_families=config.INITIAL_SIZE, seed=seed + 1)
    return len(g.profiles.families)


@benchmark("parallel_scaling")
def parallel_scaling(profile: str) -> Result:
    """Wall-clock of N independent simulations vs worker-process count.

    A fixed batch of independent full simulations run through a process pool
    (compute-bound — no disk writes, unlike ``run_replicates``). The plotter
    turns this into a speed-up curve against the ideal-linear reference.

    Note (macOS): the default start method is *spawn*, so each worker re-imports
    zombi2 — a fixed per-worker cost that caps the achievable speed-up. Linux /
    Euler default to *fork* (no re-import) and scale closer to ideal.
    """
    import os
    from concurrent.futures import ProcessPoolExecutor

    n_reps = {"quick": 8, "standard": 20, "full": 40}[profile]
    n_tips = {"quick": 3000, "standard": 10_000, "full": 10_000}[profile]
    ncores = os.cpu_count() or 4
    proc_grid = _grid(profile,
                      quick=[1, 2, 4],
                      standard=sorted({1, 2, 4, 6, 8, ncores}),
                      full=sorted({1, 2, 4, 6, 8, ncores}))
    reps = 2 if profile != "full" else 3

    points: list[Point] = []
    for p in proc_grid:
        def call(p=p):
            work = [(s, n_tips) for s in range(1, n_reps + 1)]
            if p == 1:
                return [_par_work(w) for w in work]
            with ProcessPoolExecutor(max_workers=p) as ex:
                return list(ex.map(_par_work, work))

        times, _ = measure(call, repeat=reps, warmup=0,
                           max_seconds=PROFILES[profile]["budget"] * 4)
        points.append(Point(series="run batch", x=p, times=times,
                            work=dict(n_replicates=n_reps, n_tips=n_tips)))
        print(f"  processes={p:<3d} best={min(times):7.2f} s  "
              f"({n_reps} sims of {n_tips} tips)")

    return Result(
        name="parallel_scaling",
        title="Replicate-level parallelism (compute-bound)",
        x_label="Worker processes",
        points=points,
        meta=dict(regime=config.label(), profile=profile, cpu_count=ncores,
                  n_replicates=n_reps, n_tips=n_tips),
    )


# =========================================================================
# 5. Output writing — serialise a full simulation to a ZOMBI-1 folder
# =========================================================================

@benchmark("write_output")
def write_output(profile: str) -> Result:
    """Time ``Genomes.write()`` — reconstruct every gene tree + write the folder."""
    sizes = _grid(profile,
                  quick=[100, 1000],
                  standard=[100, 1000, 10_000, 30_000],
                  full=[100, 1000, 10_000, 30_000, 100_000])
    trees = {n: z.simulate_species_tree(config.model(), n_tips=n, age=config.TREE_AGE,
                                        seed=100 + n) for n in sizes}
    rates = config.rate_model()
    tmp = Path(tempfile.mkdtemp(prefix="zombi2_write_bench_"))

    points: list[Point] = []
    try:
        for n in sizes:
            genomes = z.simulate_genomes(trees[n], rates=rates,
                                         initial_families=config.INITIAL_SIZE, seed=9000 + n)
            outdir = tmp / f"n{n}"

            def call():
                if outdir.exists():
                    shutil.rmtree(outdir)
                genomes.write(outdir)
                return outdir

            times, _ = measure(call, repeat=PROFILES[profile]["reps"], warmup=1,
                               max_seconds=PROFILES[profile]["budget"])
            n_files = sum(1 for _ in outdir.rglob("*") if _.is_file())
            points.append(Point(series="write()", x=n, times=times,
                                work=dict(n_families=len(genomes.profiles.families),
                                          n_files=n_files)))
            print(f"  write() n={n:<9d} best={min(times)*1e3:10.2f} ms  files={n_files}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    return Result(
        name="write_output",
        title="Serialising a full simulation to disk (reconstruct + write)",
        x_label="Number of extant tips",
        points=points,
        meta=dict(regime=config.label(), profile=profile),
    )


# =========================================================================
# 6. Cross-tool comparison — ZOMBI2 (Rust) vs ZOMBI 1 (legacy pure-Python)
# =========================================================================

_VS_ZOMBI1_DIR = Path(__file__).resolve().parent / "vs_zombi1"


def _zombi1_dir() -> Path:
    return Path(os.environ.get("ZOMBI1_DIR", "/Users/aadria/Desktop/Github/ZOMBI"))


def _zombi1_gm(n: int, seed: int, timeout: float) -> dict | None:
    """Time ZOMBI 1's ``Gm`` genome step at ~n extant tips via the vendored harness
    (a subprocess to the ZOMBI 1 checkout). Returns its JSON dict, or None on failure."""
    proc = subprocess.run(
        [sys.executable, str(_VS_ZOMBI1_DIR / "harness.py"),
         "--n", str(n), "--seed", str(seed), "--timeout", str(timeout)],
        capture_output=True, text=True,
    )
    lines = [ln for ln in proc.stdout.strip().splitlines() if ln.startswith("{")]
    return json.loads(lines[-1]) if lines else None


@benchmark("vs_zombi1")
def vs_zombi1(profile: str) -> Result:
    """ZOMBI2 (Rust) vs ZOMBI 1 (legacy pure-Python) on the *same* gene-family task.

    Both evolve D/T/L/O gene families over a pure-birth (Yule λ=1) tree grown to a matched
    extant-tip count, in the same regime (D=0.2 T=0.1 L=0.25 O=0.5, initial genome size 20,
    size-neutral replacement transfers). Only the genome step is timed (the species tree is
    built untimed in both). ZOMBI 1's ``Gm`` is measured through ``vs_zombi1/harness.py`` and
    is strongly super-linear — its practical ceiling is ~1200 tips in a couple of minutes —
    whereas ZOMBI2 runs the same task in milliseconds and keeps scaling to millions.

    Opt-in (not a core benchmark): needs the ZOMBI 1 checkout ($ZOMBI1_DIR, or the default
    clone path). If it is absent, only the ZOMBI2 curve is produced.
    """
    z2_sizes = _grid(profile,
                     quick=[20, 100, 1000],
                     standard=[20, 100, 300, 1000, 3000, 10_000, 30_000, 100_000],
                     full=[20, 100, 300, 1000, 3000, 10_000, 30_000, 100_000, 300_000])
    z1_sizes = _grid(profile,
                     quick=[20, 100],
                     standard=[20, 100, 300, 1000, 1200],
                     full=[20, 100, 300, 1000, 1200, 1500])
    z1_timeout = {"quick": 60.0, "standard": 180.0, "full": 240.0}[profile]

    # Same regime as ZOMBI 1's Gm harness: size-neutral (replacement) transfers so the two
    # tools do comparable work; a pure-birth Yule tree grown forward to a matched tip count.
    rates = config.rate_model()
    transfers = z.TransferModel(replacement=1.0)
    trees = {n: z.simulate_species_tree(z.Yule(config.BIRTH), n_tips=n,
                                        direction="forward", seed=200 + n) for n in z2_sizes}

    def z2_fn(n, i):
        return z.simulate_genomes(trees[n], rates=rates, initial_families=config.INITIAL_SIZE,
                                  transfers=transfers, seed=8000 + i)

    points = _timed_grid(profile, z2_sizes, z2_fn, "ZOMBI2 · Rust",
                         work_of=lambda g: dict(n_families=len(g.profiles.families)))

    # --- ZOMBI 1 (opt-in; skipped if the checkout is absent) --------------------------
    have_z1 = (_zombi1_dir() / "Zombi.py").is_file()
    if not have_z1:
        print(f"  ZOMBI 1 · Python           SKIPPED (no checkout at {_zombi1_dir()}; "
              f"set $ZOMBI1_DIR)")
    else:
        for n in z1_sizes:
            res = _zombi1_gm(n, seed=1, timeout=z1_timeout)
            if res is None:
                print(f"  ZOMBI 1 · Python           n={n:<9d} FAILED (harness error)")
                break
            if res["status"] != "ok":
                print(f"  ZOMBI 1 · Python           n={n:<9d} {res['status'].upper()} "
                      f"(> {z1_timeout:.0f}s) — practical ceiling reached")
                break
            points.append(Point(series="ZOMBI 1 · Python", x=res["n_tips"],
                                times=[res["seconds"]],
                                work=dict(n_families=res["n_gene_families"])))
            print(f"  ZOMBI 1 · Python           n={res['n_tips']:<9d} "
                  f"best={res['seconds']*1e3:10.2f} ms  "
                  f"{{'n_families': {res['n_gene_families']}}}")

    return Result(
        name="vs_zombi1",
        title="ZOMBI2 (Rust) vs ZOMBI 1 (Python): same task, gene-family simulation",
        x_label="Number of extant tips",
        points=points,
        meta=dict(regime=config.label(), profile=profile,
                  zombi1_available=have_z1,
                  note="Yule(λ=1) tree grown to matched tips; genome step only is timed; "
                       "ZOMBI 1 Gm regime D=0.2 T=0.1 L=0.25 O=0.5, initial size 20, "
                       "replacement transfers (see vs_zombi1/NOTES.md)"),
    )


# =========================================================================
# 6b. Head-to-head at ONE fixed size on the SAME species tree (boxplot data)
# =========================================================================

def _z1_build_tree(work: Path, n_tips: int, seed: int) -> tuple[str | None, str]:
    """Run ZOMBI 1's ``T`` mode once to build a species tree of ~``n_tips`` extant tips.

    Uses the vendored template (pure-birth: SPECIATION=1, EXTINCTION=0, STOPPING_RULE=1),
    so the complete tree *is* the extant tree — no extinct lineages — and both engines can
    later run their genome step over the identical tree. Returns ``(ExtantTree.nwk path, "")``
    or ``(None, error)``. Not timed."""
    tmpl = (_VS_ZOMBI1_DIR / "SpeciesTreeParameters_template.tsv").read_text()
    tmpl = tmpl.replace("__NTIPS__", str(int(n_tips))).replace("__SEED__", str(int(seed)))
    work.mkdir(parents=True, exist_ok=True)
    sp_params = work / "sp_params.tsv"
    sp_params.write_text(tmpl)
    proc = subprocess.run([sys.executable, str(_zombi1_dir() / "Zombi.py"), "T",
                           str(sp_params), str(work)],
                          cwd=str(_zombi1_dir()), capture_output=True, text=True)
    extant = work / "T" / "ExtantTree.nwk"
    if proc.returncode != 0 or not extant.is_file():
        return None, (proc.stderr or proc.stdout)[-500:]
    return str(extant), ""


def _z1_gm_once(work: Path, seed: int, timeout: float) -> tuple[float | None, int]:
    """Time ONE ZOMBI 1 ``Gm`` genome step on the tree already in ``work`` (varying the
    genome seed only). Returns ``(seconds, n_families)`` or ``(None, 0)`` on timeout/error."""
    params = (_VS_ZOMBI1_DIR / "GenomeParameters_bench.tsv").read_text()
    lines = [f"SEED {seed}" if ln.strip().startswith("SEED") else ln
             for ln in params.splitlines()]
    params_r = work / f"genome_params_s{seed}.tsv"
    params_r.write_text("\n".join(lines) + "\n")
    t0 = time.perf_counter()
    try:
        proc = subprocess.run([sys.executable, str(_zombi1_dir() / "Zombi.py"), "Gm",
                               str(params_r), str(work)],
                              cwd=str(_zombi1_dir()), capture_output=True, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, 0
    if proc.returncode != 0:
        print(f"  ZOMBI 1 · Python   seed={seed} Gm FAILED: "
              f"{(proc.stderr or proc.stdout)[-300:]}")
        return None, 0
    dt = time.perf_counter() - t0
    gf_dir = work / "G" / "Gene_families"
    nfam = sum(1 for x in gf_dir.iterdir() if x.name.endswith("_events.tsv")) \
        if gf_dir.is_dir() else 0
    return dt, nfam


@benchmark("vs_zombi1_fixedtree")
def vs_zombi1_fixedtree(profile: str) -> Result:
    """Head-to-head at a FIXED size on ONE shared species tree — boxplot of per-run times.

    Builds a single ZOMBI 1 species tree of ~1000 extant tips (pure-birth, so the extant
    tree *is* the complete tree), then runs the *genome* step of BOTH engines on that
    identical tree, once per seed. Sharing the tree removes every confound except the
    engine — same topology, same branch-length/time integral, same D/T/L/O regime — so the
    per-seed spread across replicates is an honest wall-clock distribution. ZOMBI2 reads the
    tree via ``read_newick``; ZOMBI 1 runs its ``Gm`` step on the same ``T/`` folder. Each
    Point carries every replicate time in ``times`` (the plotter box-plots them).

    Opt-in (needs the ZOMBI 1 checkout at $ZOMBI1_DIR). If it is absent, an empty result is
    returned and the overview falls back to two panels.
    """
    n_tips = {"quick": 200, "standard": 1000, "full": 1000}[profile]
    reps = {"quick": 3, "standard": 10, "full": 20}[profile]
    z1_timeout = {"quick": 60.0, "standard": 180.0, "full": 240.0}[profile]

    have_z1 = (_zombi1_dir() / "Zombi.py").is_file()
    if not have_z1:
        print(f"  SKIPPED (no ZOMBI 1 checkout at {_zombi1_dir()}; set $ZOMBI1_DIR)")
        return Result(name="vs_zombi1_fixedtree",
                      title=f"ZOMBI2 vs ZOMBI 1 on one shared {n_tips}-tip tree",
                      x_label="Wall-clock time (one full genome simulation)",
                      points=[], meta=dict(regime=config.label(), profile=profile,
                                           zombi1_available=False))

    work = Path(tempfile.mkdtemp(prefix="zombi2_fixedtree_bench_"))
    points: list[Point] = []
    try:
        # --- one shared species tree (built by ZOMBI 1's T mode, untimed) -------------
        extant, err = _z1_build_tree(work, n_tips, seed=1)
        if extant is None:
            print(f"  species-tree build FAILED: {err}")
            raise RuntimeError("tree build failed")
        nwk = Path(extant).read_text().strip()
        tree = z.read_newick(nwk)
        actual = len(tree.extant_leaves())
        print(f"  shared tree: {actual} extant tips (age {tree.total_age:.2f})")

        # --- ZOMBI2 (Rust): one timed run per seed on the shared tree -----------------
        rates = config.rate_model()
        transfers = z.TransferModel(replacement=1.0)
        z.simulate_genomes(tree, rates=rates, initial_families=config.INITIAL_SIZE,
                           transfers=transfers, seed=999)  # warm-up, untimed
        z2_times: list[float] = []
        g = None
        for r in range(1, reps + 1):
            gc.collect()
            t0 = time.perf_counter()
            g = z.simulate_genomes(tree, rates=rates, initial_families=config.INITIAL_SIZE,
                                   transfers=transfers, seed=r)
            z2_times.append(time.perf_counter() - t0)
        points.append(Point(series="ZOMBI2 · Rust", x=actual, times=z2_times,
                            work=dict(n_families=len(g.profiles.families))))
        print(f"  ZOMBI2 · Rust      {reps} runs  median={_med(z2_times)*1e3:8.2f} ms")

        # --- ZOMBI 1 (Python): one Gm per seed on the SAME T/ folder ------------------
        z1_times: list[float] = []
        nfam = 0
        for r in range(1, reps + 1):
            dt, nfam = _z1_gm_once(work, seed=r, timeout=z1_timeout)
            if dt is None:
                print(f"  ZOMBI 1 · Python   seed={r} did not finish (> {z1_timeout:.0f}s)")
                break
            z1_times.append(dt)
            print(f"  ZOMBI 1 · Python   seed={r:<2d} {dt*1e3:9.1f} ms  fams={nfam}")
        if z1_times:
            points.append(Point(series="ZOMBI 1 · Python", x=actual, times=z1_times,
                                work=dict(n_families=nfam)))
            print(f"  ZOMBI 1 · Python   {len(z1_times)} runs  "
                  f"median={_med(z1_times):8.2f} s")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    return Result(
        name="vs_zombi1_fixedtree",
        title=f"ZOMBI2 vs ZOMBI 1 on one shared {n_tips}-tip tree",
        x_label="Wall-clock time (one full genome simulation)",
        points=points,
        meta=dict(regime=config.label(), profile=profile, zombi1_available=True,
                  n_tips=n_tips, reps=reps,
                  note="ONE ZOMBI 1 pure-birth tree (~%d extant tips); both engines run their "
                       "genome step on it, one run per seed. D=0.2 T=0.1 L=0.25 O=0.5, initial "
                       "size 20, replacement transfers (see vs_zombi1/NOTES.md)" % n_tips),
    )


def _med(xs: list[float]) -> float:
    from statistics import median
    return median(xs) if xs else float("nan")
