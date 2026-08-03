"""One random stream per level, from one seed.

A seed names a *run*, not a sequence of numbers. Every level used to open
``np.random.default_rng(seed)`` directly, which meant two levels handed the same integer replayed the
**same** PCG64 stream from the same state — so ``zombi2 species --seed 42`` and ``zombi2 genomes
--seed 42`` drew the same underlying variates, and a tree's height and its genome's copy count came
out correlated (Spearman −0.79 over 6000 seeds) with nothing in the output to say so. SPEC §2 calls
two levels that do not read each other *independent*; that is a claim about two random streams, and
they were one stream.

`stream()` derives each level's generator from a `numpy.random.SeedSequence` under a per-level spawn
key, so the levels are independent by construction while ``--seed 1`` still means the same
reproducible run it always did. The keys below are part of the output contract: **never renumber
them** — a level's key is what makes last year's seed reproduce last year's run.
"""

from __future__ import annotations

import numpy as np

#: The spawn key each level draws its stream under. Fixed forever: changing a number here changes
#: every run that level has ever produced. New levels take the next unused integer.
#:
#: ``species`` is ``None`` — the **root** sequence, unkeyed, which is exactly the stream
#: ``default_rng(seed)`` used to open. Independence needs the levels' streams to differ from one
#: another, not to all be new, so one of them can stay put; species is the one worth keeping, being
#: the level every other rides on and the one whose seeds are most often published. Every species
#: tree grown under an earlier version therefore still reproduces byte for byte, and the levels that
#: were colliding with it have moved off.
_LEVEL_KEYS: dict[str, int | None] = {
    "species": None,
    "genomes": 1,
    "sequences": 2,
    "traits": 3,
    "joint": 4,          # species and trait grown together — one engine, so one stream
}


def draw_seed() -> int:
    """A fresh seed from the OS, small enough that a reader can retype it as ``--seed N``.

    Called wherever ``seed=None`` reaches a level, so that a run nobody seeded still writes down the
    seed it used. A run that wants fresh randomness still gets it — this draws from the OS every
    time. What changes is only that the draw is written down."""
    entropy = np.random.SeedSequence().entropy
    assert isinstance(entropy, int)          # the default source is one integer, not a sequence
    return int(entropy % (2 ** 31))


def resolve_seed(seed: int | None) -> int:
    """``seed`` if it was given, else a freshly drawn one. For a level that picks between two engines
    after the seed is settled: it resolves once, up front, so both engines record the same number."""
    return draw_seed() if seed is None else seed


def seed_sequence(level: str, seed: int | None) -> tuple[np.random.SeedSequence, int]:
    """The `numpy.random.SeedSequence` for ``level``, and the seed it was built from.

    What the parallel engines want: they ``.spawn()`` a stream per unit off this, so a worker count
    changes nothing. Serial callers want `stream()` instead."""
    if level not in _LEVEL_KEYS:
        raise KeyError(f"unknown level {level!r}; levels are {sorted(_LEVEL_KEYS)}")
    if seed is None:
        seed = draw_seed()
    key = _LEVEL_KEYS[level]
    seq = np.random.SeedSequence(seed) if key is None \
        else np.random.SeedSequence(seed, spawn_key=(key,))
    return seq, seed


def stream(level: str, seed: int | None) -> tuple[np.random.Generator, int]:
    """The generator for ``level``, and the seed it was built from.

    Returns the **resolved** seed, so a caller writes ``rng, seed = stream("species", seed)`` and
    stores it on the result: ``seed=None`` then leaves a run that can still be reproduced, instead of
    a `None` nobody can replay.
    """
    seq, seed = seed_sequence(level, seed)
    return np.random.default_rng(seq), seed


__all__ = ["stream", "seed_sequence", "resolve_seed", "draw_seed"]
