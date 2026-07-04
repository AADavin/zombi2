"""The one canonical simulation regime, shared by every benchmark.

Kept in one place (not scattered across benchmarks) so a scaling curve, a memory
probe, and a parallel run all describe the *same* workload — and so changing the
regime is a one-line edit that re-labels every figure consistently.
"""

from __future__ import annotations

import zombi2 as z

# Species-tree model: a mild net-diversification birth–death.
BIRTH, DEATH = 1.0, 0.3
TREE_AGE = 2.0

# Gene-family rate regime: loss > duplication so families never run away;
# moderate transfer + origination. Per-branch origination seeds new families.
RATES = dict(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5)
INITIAL_SIZE = 20


def model() -> "z.BirthDeath":
    return z.BirthDeath(BIRTH, DEATH)


def rate_model() -> "z.UniformRates":
    return z.UniformRates(**RATES)


def label() -> str:
    """One-line description of the regime for figure captions / metadata."""
    r = RATES
    return (f"BirthDeath(λ={BIRTH}, μ={DEATH}), age {TREE_AGE}; "
            f"D={r['duplication']} T={r['transfer']} L={r['loss']} O={r['origination']}, "
            f"initial_size={INITIAL_SIZE}")
