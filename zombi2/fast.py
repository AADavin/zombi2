"""Optional Rust fast-path engine (increment 1: profiles-only, independent families).

The pure-Python engine (:func:`~zombi2.simulate_genomes`) stays the default and remains the
only path that produces the full event log and gene trees. This module offers a *much*
faster route to the **presence/copy-number profile matrix** — the σ dataset — for the
built-in ``UnorderedGenome`` + ``UniformRates`` model at large scale, by running the forward
Gillespie in Rust over per-family *counts* only (no gene ids / trees / log).

It requires the compiled ``zombi2_core`` extension. Build it once from the repo::

    pip install maturin
    cd rust && maturin build --release -i python3
    pip install --force-reinstall rust/target/wheels/*.whl

If the extension isn't built, :func:`rust_available` returns ``False`` and
:func:`simulate_profiles_fast` raises a clear error.
"""

from __future__ import annotations

import numpy as np

from .genome_sim import resolve_max_family_size
from .profiles import ProfileMatrix
from .rates import UniformRates

try:  # optional native extension
    import zombi2_core as _core
except ImportError:  # pragma: no cover - depends on whether the wheel is built
    _core = None


def rust_available() -> bool:
    """True if the compiled ``zombi2_core`` extension is importable."""
    return _core is not None


def _resolve_rates(rates, duplication, transfer, loss, origination):
    """Return (d, t, l, o) from a UniformRates object or the keyword shorthand, rejecting
    features the fast path does not implement yet."""
    if rates is None:
        return float(duplication), float(transfer), float(loss), float(origination)
    # The engine implements UniformRates' per-copy semantics specifically. Other models
    # (GenomeWiseRates, FamilySampledRates, BranchRates) differ and must use simulate_genomes.
    if type(rates) is not UniformRates:
        raise TypeError(
            f"the Rust fast path only supports UniformRates or the keyword shorthand, "
            f"not {type(rates).__name__}; use simulate_genomes for that model"
        )
    # Even a UniformRates with soft carrying capacity / rearrangements isn't supported yet.
    unsupported = []
    if rates.carrying_capacity is not None:
        unsupported.append("carrying_capacity")
    if rates.inversion:
        unsupported.append("inversion")
    if rates.transposition:
        unsupported.append("transposition")
    if unsupported:
        raise ValueError(
            f"the Rust fast path does not support {', '.join(unsupported)}; "
            f"use simulate_genomes instead"
        )
    return (float(rates.duplication), float(rates.transfer),
            float(rates.loss), float(rates.origination))


def simulate_profiles_fast(
    species_tree,
    rates=None,
    *,
    duplication: float = 0.0,
    transfer: float = 0.0,
    loss: float = 0.0,
    origination: float = 0.0,
    initial_size: int = 20,
    max_family_size=None,
    seed=None,
) -> ProfileMatrix:
    """Simulate gene families along ``species_tree`` in Rust and return the profile matrix.

    Same independent-families model as :func:`~zombi2.simulate_genomes` with
    ``UnorderedGenome`` + ``UniformRates``: per-copy duplication/transfer/loss and per-branch
    origination, with additive uniform-recipient transfers and an optional hard
    ``max_family_size`` (absolute int, or a fraction of the number of species). Because the
    RNG differs from the Python engine, results are **statistically** equivalent, not
    bit-identical; a given ``seed`` is reproducible within this engine.

    Returns only a :class:`~zombi2.ProfileMatrix` — no event log or gene trees (use
    :func:`~zombi2.simulate_genomes` for those).
    """
    if _core is None:
        raise RuntimeError(
            "zombi2_core (the Rust extension) is not built. Build it with:\n"
            "  pip install maturin\n"
            "  cd rust && maturin build --release -i python3\n"
            "  pip install --force-reinstall rust/target/wheels/*.whl"
        )

    d, t, l, o = _resolve_rates(rates, duplication, transfer, loss, origination)

    nodes = list(species_tree.nodes_preorder())
    index = {n: i for i, n in enumerate(nodes)}
    parent = [index[n.parent] if n.parent is not None else -1 for n in nodes]
    times = [float(n.time) for n in nodes]
    extant_leaf = [(not n.children) and n.is_extant for n in nodes]
    root = next(i for i, n in enumerate(nodes) if n.parent is None)

    n_extant = sum(extant_leaf)
    cap = -1 if max_family_size is None else int(resolve_max_family_size(max_family_size, n_extant))
    seed_val = (int(np.random.SeedSequence(seed).generate_state(1)[0])
                if seed is None else int(seed))

    result = _core.simulate_profiles(
        len(nodes), parent, times, extant_leaf, root,
        d, t, l, o, int(initial_size), cap, seed_val,
    )
    return _assemble(result, nodes)


def _assemble(result, nodes) -> ProfileMatrix:
    """Build a ProfileMatrix from the engine's per-leaf (family_id, count) lists."""
    # column order: extant leaves by natural name order (match ProfileMatrix.from_leaf_genomes)
    from .profiles import _natkey

    leaf_cols = {leaf_idx: cols for leaf_idx, cols in result}
    species_nodes = sorted((nodes[i] for i in leaf_cols), key=lambda n: _natkey(n.name))
    species = [n.name for n in species_nodes]

    famset = {fam for cols in leaf_cols.values() for fam, _ in cols}
    families_int = sorted(famset)
    families = [str(f + 1) for f in families_int]  # 1-based labels, ZOMBI-style
    frow = {f: i for i, f in enumerate(families_int)}
    col_of = {n: j for j, n in enumerate(species_nodes)}

    matrix = np.zeros((len(families), len(species)), dtype=int)
    for leaf_idx, cols in leaf_cols.items():
        j = col_of[nodes[leaf_idx]]
        for fam, count in cols:
            matrix[frow[fam], j] = count
    return ProfileMatrix(families=families, species=species, matrix=matrix)
