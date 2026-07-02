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

from collections import Counter

from .events import EventLog, EventRecord, EventType, GeneOp
from .genome_sim import resolve_max_family_size
from .profiles import ProfileMatrix
from .rates import UniformRates
from .simulation import Genomes

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

    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    cap, seed_val = _cap_and_seed(max_family_size, sum(extant_leaf), seed)

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


# --- shared tree/cap/seed plumbing -----------------------------------------------

def _tree_arrays(species_tree):
    """Flatten a Tree into the index-based arrays the Rust engines consume."""
    nodes = list(species_tree.nodes_preorder())
    index = {n: i for i, n in enumerate(nodes)}
    parent = [index[n.parent] if n.parent is not None else -1 for n in nodes]
    times = [float(n.time) for n in nodes]
    extant_leaf = [(not n.children) and n.is_extant for n in nodes]
    root = next(i for i, n in enumerate(nodes) if n.parent is None)
    return nodes, parent, times, extant_leaf, root


def _cap_and_seed(max_family_size, n_extant, seed):
    cap = -1 if max_family_size is None else int(resolve_max_family_size(max_family_size, n_extant))
    seed_val = (int(np.random.SeedSequence(seed).generate_state(1)[0])
                if seed is None else int(seed))
    return cap, seed_val


# --- full-log fast path ----------------------------------------------------------

class _FastGene:
    """Minimal Gene stand-in for gid->species lookup (int gid, str family)."""
    __slots__ = ("gid", "family")

    def __init__(self, gid, family):
        self.gid = gid
        self.family = family


class _FastGenome:
    """Minimal leaf genome exposing just what Genomes/ProfileMatrix need."""
    __slots__ = ("_pairs", "_counts")

    def __init__(self, pairs):  # pairs: list[(gid:int, family:int)]
        self._pairs = pairs
        self._counts = Counter(str(f + 1) for _, f in pairs)  # 1-based labels, ZOMBI-style

    def genes(self):
        return [_FastGene(gid, str(fam + 1)) for gid, fam in self._pairs]

    def families(self):
        return list(self._counts.keys())

    def copy_number(self, family):
        return self._counts.get(family, 0)


# event code (Rust) -> EventType, and per-event GeneOp roles for genes[0], genes[1], genes[2]
_EV = (EventType.ORIGINATION, EventType.DUPLICATION, EventType.TRANSFER,
       EventType.LOSS, EventType.SPECIATION)
_ROLES = (
    ("origin",),                            # O
    ("parent", "left", "right"),            # D
    ("parent", "donor_copy", "transfer_copy"),  # T
    ("lost",),                              # L
    ("parent", "child", "child"),           # S
)


def _build_event_log(cols, nodes) -> EventLog:
    ev, br, tm, dn, rc, fm, g0, g1, g2 = cols
    names = [n.name for n in nodes]
    log = EventLog()
    add = log.add
    for k in range(len(ev)):
        code = ev[k]
        family = str(fm[k] + 1)  # 1-based labels, ZOMBI-style
        roles = _ROLES[code]
        genes = [GeneOp(g0[k], family, roles[0])]
        if len(roles) == 3:
            genes.append(GeneOp(g1[k], family, roles[1]))
            genes.append(GeneOp(g2[k], family, roles[2]))
        donor = names[dn[k]] if dn[k] >= 0 else None
        recipient = names[rc[k]] if rc[k] >= 0 else None
        add(EventRecord(_EV[code], names[br[k]], tm[k], genes,
                        donor=donor, recipient=recipient))
    return log


def simulate_genomes_fast(
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
) -> Genomes:
    """Simulate gene families in Rust and return a full :class:`~zombi2.Genomes` result.

    Unlike :func:`simulate_profiles_fast`, this tracks individual gene lineages and emits the
    complete event genealogy, so the returned ``Genomes`` supports ``.event_log``,
    ``.gene_trees()`` and ``.write()`` exactly like :func:`~zombi2.simulate_genomes`.

    Same built-in model as :func:`~zombi2.simulate_genomes` with ``UnorderedGenome`` +
    ``UniformRates`` and the **default** ``TransferModel`` (additive, uniform recipient, no
    self-transfer); a hard ``max_family_size`` forces over-cap transfers to replacements.
    Because the RNG differs from the Python engine, results are statistically equivalent, not
    bit-identical (gene ids are integers rather than ``g``-prefixed strings). For custom
    transfer mechanics, other rate models, or ordered genomes, use
    :func:`~zombi2.simulate_genomes`.
    """
    if _core is None:
        raise RuntimeError(
            "zombi2_core (the Rust extension) is not built. Build it with:\n"
            "  pip install maturin\n"
            "  cd rust && maturin build --release -i python3\n"
            "  pip install --force-reinstall rust/target/wheels/*.whl"
        )

    d, t, l, o = _resolve_rates(rates, duplication, transfer, loss, origination)
    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    cap, seed_val = _cap_and_seed(max_family_size, sum(extant_leaf), seed)

    cols, leaves = _core.simulate_log(
        len(nodes), parent, times, extant_leaf, root,
        d, t, l, o, int(initial_size), cap, seed_val,
    )

    leaf_genomes = {nodes[li]: _FastGenome(pairs) for li, pairs in leaves}
    event_log = _build_event_log(cols, nodes)
    profiles = ProfileMatrix.from_leaf_genomes(leaf_genomes)
    return Genomes(species_tree=species_tree, leaf_genomes=leaf_genomes,
                   event_log=event_log, profiles=profiles)
