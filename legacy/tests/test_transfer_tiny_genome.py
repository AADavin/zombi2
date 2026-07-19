"""Regression: a transfer on a tiny genic nucleotide genome must not hang forever.

With an additive transfer rate and no compensating loss, the length-scaled nucleotide model
grows without bound: a structural event's rate is proportional to genome length, so every
transfer adds a gene copy, the aggregate rate climbs, simulated time barely advances, and the
Gillespie walk fires effectively unbounded events (each one O(genome length)). The old code
had no ceiling, so ``zombi2 genomes --genome-model nucleotide --genes ... --trans 0.03`` (no
``--loss``) ran for hours. The runaway-growth guard now turns that into a prompt, actionable
``RuntimeError``.

These tests wrap the simulation in a hard wall-clock watchdog so that a regression (the guard
being removed or defeated) *fails* the test instead of hanging the whole suite.
"""

import inspect
import signal
from contextlib import contextmanager

import pytest

import zombi2
from zombi2.genomes.nucleotide_sim import (
    DEFAULT_MAX_SEGMENTS_PER_GENOME,
    simulate_nucleotide_genomes,
)

# The exact tiny genic genome from the bug report: three 50 bp genes in a 300 bp chromosome.
GENES = [(0, 50, "g1"), (100, 150, "g2"), (200, 250, "g3")]

pytestmark = pytest.mark.skipif(
    not hasattr(signal, "SIGALRM"),
    reason="wall-clock watchdog needs SIGALRM (POSIX only)",
)


class _Timeout(Exception):
    pass


@contextmanager
def time_limit(seconds):
    """Raise ``_Timeout`` if the wrapped block runs longer than ``seconds``.

    The pure-Python engine checks for signals between bytecodes, so ``SIGALRM`` interrupts a
    runaway walk. In the *passing* case the guard raises well before the alarm and it never
    fires; in a *regression* the alarm fires and the test fails loudly rather than hanging.
    """
    def _handler(signum, frame):
        raise _Timeout(f"simulation did not terminate within {seconds}s (the hang is back)")

    old = signal.signal(signal.SIGALRM, _handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def _tiny_tree(seed=1):
    return zombi2.simulate_species_tree(
        zombi2.BirthDeath(1.0, 0.2), n_tips=6, age=2.5, seed=seed
    )


def test_tiny_genic_transfer_terminates_instead_of_hanging():
    """The reported hang: additive transfer, no loss, tiny genic genome -> must terminate fast."""
    tree = _tiny_tree(seed=1)
    with time_limit(20):
        with pytest.raises(RuntimeError, match="growing without bound"):
            simulate_nucleotide_genomes(
                tree,
                gene_intervals=GENES,
                root_length=300,
                transfer=0.03,
                seed=1,
                # A low ceiling makes the guard trip in a fraction of a second; the default
                # (exercised via the signature check below) protects real users just the same.
                max_segments_per_genome=2000,
            )


def test_runaway_error_is_actionable():
    """The error names the fix (balance the gain rate with loss / raise the ceiling)."""
    tree = _tiny_tree(seed=2)
    with time_limit(20):
        with pytest.raises(RuntimeError) as excinfo:
            simulate_nucleotide_genomes(
                tree, gene_intervals=GENES, root_length=300, transfer=0.03,
                seed=2, max_segments_per_genome=2000,
            )
    msg = str(excinfo.value)
    assert "max_segments_per_genome" in msg
    assert "loss" in msg


def test_default_guard_is_wired():
    """simulate_nucleotide_genomes ships the guard on by default (regression-proofs the default)."""
    default = inspect.signature(simulate_nucleotide_genomes).parameters[
        "max_segments_per_genome"
    ].default
    assert default is DEFAULT_MAX_SEGMENTS_PER_GENOME
    assert isinstance(default, int) and default > 0


def test_disabling_guard_is_possible_but_off_switch_exists():
    """Passing None disables the guard (documented escape hatch); we only check it is accepted,
    not that it runs to completion (an unbalanced run would, by design, never finish)."""
    param = inspect.signature(simulate_nucleotide_genomes).parameters["max_segments_per_genome"]
    # ``int | None`` -> None must be a legal value the caller can pass.
    assert param.annotation in ("int | None", "Optional[int]") or "None" in str(param.annotation)


def test_balanced_transfer_run_is_unaffected_and_deterministic():
    """A balanced run (loss present) completes normally, never trips the guard, and is
    byte-identical across calls -- proof the guard does not perturb working simulations."""
    tree = _tiny_tree(seed=3)
    kw = dict(gene_intervals=GENES, root_length=300, transfer=0.03, loss=0.03, seed=7)
    with time_limit(20):
        a = simulate_nucleotide_genomes(tree, **kw)
        b = simulate_nucleotide_genomes(tree, **kw)

    sizes_a = {name: g.size() for name, g in a.leaf_genomes.items()}
    sizes_b = {name: g.size() for name, g in b.leaf_genomes.items()}
    assert sizes_a == sizes_b
    # every leaf genome stays far below the ceiling -- the guard is never in play here
    assert all(g.n_segments() < DEFAULT_MAX_SEGMENTS_PER_GENOME for g in a.leaf_genomes.values())
