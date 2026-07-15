"""Byte-identity guardrail for the chromosome-tier refactor (Stage 0).

Freezes leaf-genome fingerprints for a fixed seed x config matrix, captured from the
*pre-refactor* engine, and asserts every later refactor stage reproduces them exactly. This is
the safety net the whole refactor leans on: promoting the ordered/nucleotide genomes onto a
first-class ``Chromosome`` tier must not change a single simulated genome until we deliberately
turn on a new event.

The fingerprint **includes gene ids**, on purpose: it therefore also catches id-counter
contamination — e.g. a chromosome id accidentally minted from the *gene* counter would shift
every gid without changing structure, and a family/orientation-only fingerprint would miss it.

Regenerate the baseline ONLY when a change to the fingerprints is intended and understood:

    ZOMBI2_UPDATE_FINGERPRINTS=1 python -m pytest tests/test_chromosome_fingerprints.py

Never regenerate to "make the test pass" during a refactor stage — drift is the signal.
"""
import hashlib
import json
import os
from pathlib import Path

import pytest

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    Rates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes

BASELINE_PATH = Path(__file__).parent / "data" / "chromosome_fingerprints.json"


# --- fingerprints (canonical, gid-inclusive) ------------------------------------------

def _fingerprint_ordered(genomes) -> str:
    """Every leaf's chromosomes, order preserved, as ``leaf=gid:family:+//gid:family:-...``.

    NOTE: reads ``genome.chromosomes`` as an iterable of per-chromosome gene sequences. When the
    representation changes to ``dict[chrom_id, Chromosome]`` (Stage 1) this extraction adapts but
    MUST still emit the identical string — that identity is exactly what byte-identity means here.
    """
    parts = []
    for leaf, genome in sorted(genomes.leaf_genomes.items(), key=lambda kv: kv[0].name):
        chroms = ["|".join(f"{g.gid}:{g.family}:{g.orientation:+d}" for g in chrom.genes)
                  for chrom in genome.chromosomes.values()]
        parts.append(f"{leaf.name}=" + "//".join(chroms))
    return ";".join(parts)


def _fingerprint_nucleotide(res) -> str:
    """Every leaf's nucleotide layout as ``leaf=source:pos:+,source:pos:-,...`` (``to_cells``)."""
    parts = []
    for leaf, genome in sorted(res.leaf_genomes.items(), key=lambda kv: kv[0].name):
        cells = ",".join(f"{src}:{pos}:{st:+d}" for (src, pos, st) in genome.to_cells())
        parts.append(f"{leaf.name}=" + cells)
    return ";".join(parts)


# --- config matrix --------------------------------------------------------------------

def _ordered_run(seed, n_chromosomes, circular):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=seed)
    rates = Rates(duplication=0.4, loss=0.3, transfer=0.2, origination=0.1,
                        inversion=0.3, transposition=0.3)
    return simulate_genomes(
        tree, rates, initial_families=20, seed=seed,
        genome_factory=lambda ids: OrderedGenome(ids, extension=0.6,
                                                 n_chromosomes=n_chromosomes, circular=circular),
    )


def _nucleotide_run(seed, initial_chromosomes=1, *, inversion=0.02, duplication=0.005, loss=0.005,
                    transposition=0.01, transfer=0.0, origination=0.0, insertion=0.0, deletion=0.0,
                    root_length=200, extension=0.9, gene_intervals=None, pseudogenization=0.0,
                    replacement=0.0):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=seed)
    return simulate_nucleotide_genomes(
        tree, inversion=inversion, duplication=duplication, loss=loss, transposition=transposition,
        transfer=transfer, origination=origination, insertion=insertion, deletion=deletion,
        root_length=root_length, extension=extension, initial_chromosomes=initial_chromosomes,
        gene_intervals=gene_intervals, pseudogenization=pseudogenization, replacement=replacement,
        seed=seed,
    )


#: name -> (fingerprint kind, zero-arg thunk running the sim). Ordered covers the byte-identity
#: surface for Stages 1-5; nucleotide is the Stage-6 surface. Single-chromosome ordered is the
#: one that must NEVER change, so it is covered at two seeds.
CONFIGS = {
    "ordered/single-circular/s1": ("ordered", lambda: _ordered_run(1, 1, True)),
    "ordered/single-circular/s7": ("ordered", lambda: _ordered_run(7, 1, True)),
    "ordered/multi3-circular/s4": ("ordered", lambda: _ordered_run(4, 3, True)),
    "ordered/multi4-circular/s4": ("ordered", lambda: _ordered_run(4, 4, True)),
    "ordered/multi8-linear/s5":   ("ordered", lambda: _ordered_run(5, 8, False)),
    "ordered/multi3-linear/s6":   ("ordered", lambda: _ordered_run(6, 3, False)),
    "nucleotide/single/s7":       ("nucleotide", lambda: _nucleotide_run(7, 1)),
    "nucleotide/multi3/s9":       ("nucleotide", lambda: _nucleotide_run(9, 3)),
    # Stage-6 surface: cover every event the nucleotide migration touches.
    "nucleotide/transfer/s3":     ("nucleotide", lambda: _nucleotide_run(3, 1, transfer=0.02, inversion=0.01)),
    "nucleotide/origination/s4":  ("nucleotide", lambda: _nucleotide_run(4, 1, origination=0.05)),
    "nucleotide/indels/s5":       ("nucleotide", lambda: _nucleotide_run(5, 1, insertion=0.01, deletion=0.01)),
    "nucleotide/genic/s6":        ("nucleotide", lambda: _nucleotide_run(
        6, 1, gene_intervals=[(20, 60), (100, 150)], loss=0.01, transfer=0.02,
        pseudogenization=0.5, replacement=0.5)),
    "nucleotide/rich-multi/s2":   ("nucleotide", lambda: _nucleotide_run(
        2, 2, inversion=0.02, duplication=0.01, loss=0.01, transposition=0.01, transfer=0.02,
        origination=0.03, insertion=0.005, deletion=0.005)),
}

_FINGERPRINT = {"ordered": _fingerprint_ordered, "nucleotide": _fingerprint_nucleotide}


def compute_fingerprints() -> dict:
    return {name: hashlib.sha256(_FINGERPRINT[kind](thunk()).encode()).hexdigest()
            for name, (kind, thunk) in CONFIGS.items()}


# --- tests ----------------------------------------------------------------------------

def test_fingerprints_match_baseline():
    """Every config reproduces the frozen pre-refactor fingerprint (or --update regenerates)."""
    current = compute_fingerprints()

    if os.environ.get("ZOMBI2_UPDATE_FINGERPRINTS"):
        BASELINE_PATH.parent.mkdir(parents=True, exist_ok=True)
        BASELINE_PATH.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")
        pytest.skip(f"regenerated {len(current)} fingerprints -> {BASELINE_PATH}")

    assert BASELINE_PATH.exists(), (
        f"no baseline at {BASELINE_PATH}; generate it once from the current (pre-refactor) tree "
        f"with ZOMBI2_UPDATE_FINGERPRINTS=1 python -m pytest {Path(__file__).name}")
    baseline = json.loads(BASELINE_PATH.read_text())

    assert set(current) == set(baseline), (
        f"config set changed (added {sorted(set(current) - set(baseline))}, "
        f"removed {sorted(set(baseline) - set(current))}); regenerate the baseline intentionally")
    drifted = [name for name in sorted(baseline) if current[name] != baseline[name]]
    assert not drifted, (
        "BYTE-IDENTITY DRIFT — a refactor changed simulation output for: " + ", ".join(drifted) +
        ". If this change is intended, regenerate with ZOMBI2_UPDATE_FINGERPRINTS=1.")


def test_fingerprints_are_reproducible():
    """Determinism, independent of the frozen baseline: same configs, same fingerprints twice."""
    assert compute_fingerprints() == compute_fingerprints()
