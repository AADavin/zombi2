"""Internal bridge to the compiled ``zombi2_core`` Rust engine.

This module is **not** part of the public API. It is called by
:func:`~zombi2.simulate_genomes` and :func:`~zombi2.simulate_nucleotide_genomes`, which
route the *built-in* model (order-free ``UnorderedGenome`` + ``UniformRates``) to Rust
automatically — there is no separate "fast" function and no engine switch. Flexible models
(family/genome-wise/branch rates, ordered genomes, carrying capacity, custom samplers) run
on the pure-Python engine instead.

The extension must be compiled. Build it once from the repo::

    pip install maturin
    cd rust && maturin build --release -i python3
    pip install --force-reinstall rust/target/wheels/*.whl

If it isn't built, :func:`available` returns ``False`` and the built-in model raises a clear
error telling the user to build it (there is deliberately no Python fallback for the
built-in model, so a given ``seed`` is reproducible against a single engine).
"""

from __future__ import annotations

import numpy as np

from collections import Counter

from .events import EventLog, EventRecord, EventType, GeneOp
from .genome_sim import resolve_max_family_size
from .nucleotide_sim import NucleotideResult
from .profiles import ProfileMatrix
from .rates import UniformRates
from .simulation import Genomes

try:  # optional native extension
    import zombi2_core as _core
except ImportError:  # pragma: no cover - depends on whether the wheel is built
    _core = None


_BUILD_HINT = (
    "the built-in gene-family model runs on the compiled zombi2_core extension, which is not "
    "built. Build it with:\n"
    "  pip install maturin\n"
    "  cd rust && maturin build --release -i python3\n"
    "  pip install --force-reinstall rust/target/wheels/*.whl\n"
    "or use a flexible rate model / ordered genome (which run on the Python engine)."
)


def available() -> bool:
    """True if the compiled ``zombi2_core`` extension is importable."""
    return _core is not None


def require() -> None:
    """Raise a clear build error if the extension is missing."""
    if _core is None:
        raise RuntimeError(_BUILD_HINT)


def eligible(rates, genome_factory, sampler) -> bool:
    """True if this model is the built-in one the Rust engine implements: the default
    ``UnorderedGenome``, a plain ``UniformRates`` (no soft carrying capacity, no
    rearrangements), and no custom event sampler. Everything else runs on Python."""
    from .genome import UnorderedGenome

    return (
        genome_factory is UnorderedGenome
        and sampler is None
        and type(rates) is UniformRates
        and rates.carrying_capacity is None
        and not rates.inversion
        and not rates.transposition
    )


class _FastNucleotideResult(NucleotideResult):
    """A NucleotideResult from the Rust profile path (no event log → no gene trees)."""

    def block_gene_trees(self):
        raise NotImplementedError(
            "the Rust nucleotide profile path emits only leaf segments (profile / mosaic / "
            "trace-back); per-block gene trees need the event log — use "
            "simulate_nucleotide_genomes(...) (the default) for those."
        )

    def block_histories(self):
        raise NotImplementedError(
            "block_histories needs the event log; use simulate_nucleotide_genomes(...)."
        )


# --- rate / tree / cap / transfer plumbing ---------------------------------------

def _resolve_rates(rates):
    """Return (d, t, l, o) from a UniformRates, rejecting features Rust does not implement."""
    if type(rates) is not UniformRates:
        raise TypeError(
            f"the Rust engine only supports UniformRates, not {type(rates).__name__}; "
            f"use simulate_genomes for that model"
        )
    unsupported = []
    if rates.carrying_capacity is not None:
        unsupported.append("carrying_capacity")
    if rates.inversion:
        unsupported.append("inversion")
    if rates.transposition:
        unsupported.append("transposition")
    if unsupported:
        raise ValueError(
            f"the Rust engine does not support {', '.join(unsupported)}; "
            f"use a flexible model on the Python engine"
        )
    return (float(rates.duplication), float(rates.transfer),
            float(rates.loss), float(rates.origination))


def _tree_arrays(species_tree):
    """Flatten a Tree into the index-based arrays the Rust engines consume.

    A single pre-order pass: because a node is always visited before its children, each parent's
    index is already assigned when its children are reached, so the parent pointers, times and
    extant-leaf mask are filled in one sweep. The previous version made ~six passes over the ~2N
    nodes (materialise, build an object-keyed index dict, then a comprehension each for parents /
    times / mask / root) — this is ~1.6x faster at a million tips. It stays Python-bound because
    it walks the ``TreeNode`` object graph; a Rust-speed tree would have to be parsed/built as
    arrays in Rust rather than materialised as Python objects first."""
    nodes = []
    parent = []
    times = []
    extant_leaf = []
    pos: dict[int, int] = {}       # id(node) -> pre-order index (parent seen before child)
    root = -1
    for i, n in enumerate(species_tree.nodes_preorder()):
        nodes.append(n)
        pos[id(n)] = i
        p = n.parent
        if p is None:
            parent.append(-1)
            root = i
        else:
            parent.append(pos[id(p)])
        times.append(float(n.time))
        extant_leaf.append(n.is_extant and not n.children)
    return nodes, parent, times, extant_leaf, root


def _cap_and_seed(max_family_size, n_extant, seed):
    cap = -1 if max_family_size is None else int(resolve_max_family_size(max_family_size, n_extant))
    seed_val = (int(np.random.SeedSequence(seed).generate_state(1)[0])
                if seed is None else int(seed))
    return cap, seed_val


def _transfer_params(transfers):
    """(replacement, distance_decay, allow_self) for Rust; decay = -1.0 means uniform."""
    if transfers is None:
        return 0.0, -1.0, False
    decay = transfers.distance_decay
    return (float(transfers.replacement),
            -1.0 if decay is None else float(decay),
            bool(transfers.allow_self))


# --- profiles (counts only; the fast ABC / σ-dataset path) -----------------------

def profiles(species_tree, rates, *, initial_size, transfers, max_family_size, seed):
    """Simulate gene families in Rust over per-family *counts* only and return the profile
    matrix — no gene ids, event log or trees. This is the path behind
    ``simulate_genomes(..., output="profiles")`` (the fast route for ABC / large σ datasets)."""
    require()
    d, t, l, o = _resolve_rates(rates)
    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    cap, seed_val = _cap_and_seed(max_family_size, sum(extant_leaf), seed)
    rep, dec, aself = _transfer_params(transfers)

    result = _core.simulate_profiles(
        len(nodes), parent, times, extant_leaf, root,
        d, t, l, o, int(initial_size), cap, seed_val, rep, dec, aself,
    )
    return _assemble_profiles(result, nodes)


def profiles_parallel(species_tree, rates, *, initial_size, transfers, max_family_size, seed,
                      threads):
    """Parallel counts-only path (behind ``simulate_genomes(..., output="profiles", threads=N)``).

    Runs ``threads`` **independent** copies of the engine, each with origination rate ``o/threads``
    and a ``1/threads`` share of the founding families, and sums the profiles. Gene families are
    independent and a Poisson process splits, so this is distributionally identical to one serial
    run — a different but equivalent realization. The result depends on ``(seed, threads)`` but not
    on scheduling (thread-count-independent given the copy count), and every transfer mode is
    supported (the recipient choice depends only on the tree + donor, so it decomposes unchanged)."""
    require()
    d, t, l, o = _resolve_rates(rates)
    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    cap, seed_val = _cap_and_seed(max_family_size, sum(extant_leaf), seed)
    rep, dec, aself = _transfer_params(transfers)

    # copies = threads (one balanced Poisson-thinned copy per worker); pool size = threads
    result = _core.simulate_profiles_parallel(
        len(nodes), parent, times, extant_leaf, root,
        d, t, l, o, int(initial_size), cap, seed_val, rep, dec, aself,
        int(threads), int(threads),
    )
    return _assemble_profiles(result, nodes)


def _assemble_profiles(result, nodes) -> ProfileMatrix:
    """Build a ProfileMatrix from the engine's flat COO byte buffers (``simulate_profiles`` /
    ``simulate_trace``).

    ``result`` is ``(leaf_nodes, col, fam, cnt)`` of packed native-endian bytes: ``leaf_nodes[c]``
    is the node index of column ``c`` (columns in the engine's emission order), then one cell per
    present family as parallel ``u32`` column / ``u64`` family id / ``u32`` count. Everything is
    ``np.frombuffer``-d (one memcpy, no per-cell Python objects) and mapped to sorted family rows
    and natkey-ordered species columns with vectorised numpy — the per-cell Python loop this
    replaces was ~40% of the wall-clock at millions of tips. Stays sparse COO throughout, so the
    dense O(N²) array is never built."""
    from .profiles import _natkey

    leaf_nodes_b, col_b, fam_b, cnt_b = result
    leaf_idx = np.frombuffer(leaf_nodes_b, dtype=np.uint32)
    col = np.frombuffer(col_b, dtype=np.uint32)
    fam = np.frombuffer(fam_b, dtype=np.uint64)
    cnt = np.frombuffer(cnt_b, dtype=np.uint32)

    # species (columns): extant-leaf names in natkey order; remap the engine's emission-order
    # column index to that order via the inverse permutation.
    leaf_names = [nodes[i].name for i in leaf_idx.tolist()]
    order = sorted(range(len(leaf_names)), key=lambda j: _natkey(leaf_names[j]))
    species = [leaf_names[j] for j in order]
    inv = np.empty(len(order), dtype=np.int64)
    inv[np.asarray(order, dtype=np.int64)] = np.arange(len(order), dtype=np.int64)
    cols = inv[col.astype(np.int64)] if col.size else np.zeros(0, dtype=np.int64)

    # families (rows): sorted unique family ids, 1-based ZOMBI-style labels.
    uniq = np.unique(fam) if fam.size else np.zeros(0, dtype=np.uint64)
    families = [str(int(f) + 1) for f in uniq.tolist()]
    rows = np.searchsorted(uniq, fam).astype(np.int64) if fam.size else np.zeros(0, dtype=np.int64)

    return ProfileMatrix(families=families, species=species,
                         coo=(rows, cols, cnt.astype(np.int64)))


# --- full genealogy (event log + gene trees), materialized as a Genomes ----------

class _RustGene:
    """Minimal Gene stand-in for gid->species lookup (str gid, str family)."""
    __slots__ = ("gid", "family")

    def __init__(self, gid, family):
        self.gid = gid
        self.family = family


class _RustGenome:
    """Minimal leaf genome exposing just what Genomes/ProfileMatrix need."""
    __slots__ = ("_pairs", "_counts")

    def __init__(self, pairs):  # pairs: list[(gid:int, family:int)]
        self._pairs = pairs
        self._counts = Counter(str(f + 1) for _, f in pairs)  # 1-based labels, ZOMBI-style

    def genes(self):
        # g-prefixed gids so the schema matches the Python engine (only the RNG differs)
        return [_RustGene(f"g{gid}", str(fam + 1)) for gid, fam in self._pairs]

    def families(self):
        return list(self._counts.keys())

    def copy_number(self, family):
        return self._counts.get(family, 0)

    def size(self):
        return len(self._pairs)


# event code (Rust) -> EventType, and per-event GeneOp roles for genes[0], genes[1], genes[2]
_EV = (EventType.ORIGINATION, EventType.DUPLICATION, EventType.TRANSFER,
       EventType.LOSS, EventType.SPECIATION)
_ROLES = (
    ("origin",),                                # O
    ("parent", "left", "right"),                # D
    ("parent", "donor_copy", "transfer_copy"),  # T
    ("lost",),                                  # L
    ("parent", "child", "child"),               # S
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
        genes = [GeneOp(f"g{g0[k]}", family, roles[0])]  # g-prefixed gids (schema parity)
        if len(roles) == 3:
            genes.append(GeneOp(f"g{g1[k]}", family, roles[1]))
            genes.append(GeneOp(f"g{g2[k]}", family, roles[2]))
        donor = names[dn[k]] if dn[k] >= 0 else None
        recipient = names[rc[k]] if rc[k] >= 0 else None
        add(EventRecord(_EV[code], names[br[k]], tm[k], genes,
                        donor=donor, recipient=recipient))
    return log


def events_trace_tsv(columns, nodes) -> str:
    """Serialise the raw Rust event columns as the compact ``Events_trace.tsv`` text — the fast
    path behind ``output="trace"``. This never builds an :class:`~zombi2.events.EventRecord`;
    it formats the columns directly, which is the whole point (object construction is the wall
    on large trees). The format matches the record-based
    :func:`zombi2.simulation.events_trace_from_log` exactly (g-prefixed gids, 1-based families)."""
    from .simulation import EVENTS_TRACE_HEADER

    ev, br, tm, dn, rc, fm, g0, g1, g2 = columns
    n = len(br)
    names = np.array([nd.name for nd in nodes], dtype=object)
    ev_a = np.frombuffer(ev, dtype=np.uint8)
    br_a = np.asarray(br, dtype=np.int64)
    dn_a = np.asarray(dn, dtype=np.int64)
    rc_a = np.asarray(rc, dtype=np.int64)
    # resolve names / event chars vectorially; -1 (no donor/recipient) becomes ""
    evchar = np.array([e.value for e in _EV], dtype=object)[ev_a].tolist()
    brn = names[br_a].tolist()
    dnn = np.where(dn_a >= 0, names[np.clip(dn_a, 0, None)], "").tolist()
    rcn = np.where(rc_a >= 0, names[np.clip(rc_a, 0, None)], "").tolist()
    tm_l = list(tm)
    fm_l = (np.asarray(fm, dtype=np.int64) + 1).tolist()  # 1-based family labels
    g0_l, g1_l, g2_l = list(g0), list(g1), list(g2)

    def gid(x):
        return f"g{x}" if x >= 0 else ""

    rows = [EVENTS_TRACE_HEADER]
    rows.extend(
        f"{tm_l[k]:.10g}\t{evchar[k]}\t{brn[k]}\t{dnn[k]}\t{rcn[k]}\t{fm_l[k]}\t"
        f"g{g0_l[k]}\t{gid(g1_l[k])}\t{gid(g2_l[k])}"
        for k in range(n)
    )
    return "\n".join(rows) + "\n"


def trace(species_tree, rates, *, initial_size, transfers, max_family_size, seed):
    """Simulate the compact **event trace** in Rust and return a :class:`~zombi2.GenomeTrace`.

    Uses the ``simulate_trace`` engine: identical D/T/L/O dynamics to :func:`genomes`, but
    speciations neither re-mint gene ids nor emit records (a lineage keeps its id across
    speciations). The returned columns are therefore O/D/T/L only — ~6x smaller than the full
    log in the low-rate regime — and the per-event Python objects are never built. The
    genealogy is reconstructed on demand by replaying the trace against the species tree
    (:func:`zombi2.reconciliation.expand_trace`). This is the path behind
    ``simulate_genomes(..., output="trace")``."""
    require()
    d, t, l, o = _resolve_rates(rates)
    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    cap, seed_val = _cap_and_seed(max_family_size, sum(extant_leaf), seed)
    rep, dec, aself = _transfer_params(transfers)

    cols, leaf_coo = _core.simulate_trace(
        len(nodes), parent, times, extant_leaf, root,
        d, t, l, o, int(initial_size), cap, seed_val, rep, dec, aself,
    )
    # leaf_coo is the flat COO profile buffers (per-family counts, no per-gene objects); assemble
    # straight from it. The compact trace's leaf identities are recovered by replaying against the
    # species tree, not read off leaf genomes, so leaf_genomes stays empty for this path.
    profs = _assemble_profiles(leaf_coo, nodes)

    from .simulation import GenomeTrace
    return GenomeTrace(species_tree=species_tree, leaf_genomes={},
                       profiles=profs, _columns=cols, _nodes=nodes)


def genomes(species_tree, rates, *, initial_size, transfers, max_family_size, seed):
    """Simulate the full genealogy in Rust and return a materialized :class:`~zombi2.Genomes`
    (event log, gene trees, profiles, write). This is the built-in path behind the default
    ``simulate_genomes(...)`` (``output="genomes"``)."""
    require()
    d, t, l, o = _resolve_rates(rates)
    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    cap, seed_val = _cap_and_seed(max_family_size, sum(extant_leaf), seed)
    rep, dec, aself = _transfer_params(transfers)

    cols, leaves = _core.simulate_log(
        len(nodes), parent, times, extant_leaf, root,
        d, t, l, o, int(initial_size), cap, seed_val, rep, dec, aself,
    )

    leaf_genomes = {nodes[li]: _RustGenome(pairs) for li, pairs in leaves}
    event_log = _build_event_log(cols, nodes)
    profs = ProfileMatrix.from_leaf_genomes(leaf_genomes)
    return Genomes(species_tree=species_tree, leaf_genomes=leaf_genomes,
                   event_log=event_log, profiles=profs)


# --- nucleotide-genome profile path (leaf segments -> blocks/profile/mosaics) ------

class _FastNucGenome:
    """Minimal leaf genome exposing just what NucleotideResult reads (segments + to_cells)."""
    __slots__ = ("_segments", "_length")

    def __init__(self, segments):
        self._segments = segments
        self._length = sum(s.length for s in segments)

    def n_segments(self):
        return len(self._segments)

    def size(self):
        return self._length

    def to_cells(self):
        out = []
        for seg in self._segments:
            if seg.strand == 1:
                out.extend((seg.source, p, 1) for p in range(seg.src_start, seg.src_end))
            else:
                out.extend((seg.source, p, -1) for p in range(seg.src_end - 1, seg.src_start - 1, -1))
        return out


def nucleotide(species_tree, *, inversion, loss, duplication, transfer, transposition,
               origination, root_length, extension, initial_size, transfers, seed):
    """Simulate the nucleotide structural model in Rust (leaf segments only) and return a
    :class:`NucleotideResult` supporting ``profile_matrix()`` / ``leaf_mosaic()`` /
    ``trace_back()`` but not the event-log products. This is the path behind
    ``simulate_nucleotide_genomes(..., output="profiles")``."""
    require()
    from .events import EventLog
    from .nucleotide_genome import Segment, SegmentRegistry
    from .nucleotide_sim import _build_blocks

    nodes, parent, times, extant_leaf, root = _tree_arrays(species_tree)
    _, seed_val = _cap_and_seed(None, sum(extant_leaf), seed)
    rep, dec, aself = _transfer_params(transfers)

    leaves = _core.simulate_nucleotide(
        len(nodes), parent, times, extant_leaf, root,
        float(inversion), float(loss), float(duplication), float(transfer),
        float(transposition), float(origination), int(root_length), float(extension),
        int(initial_size), seed_val, rep, dec, aself,
    )

    leaf_genomes = {}
    for leaf_idx, segs in leaves:
        objs = [Segment("", str(src + 1), start, end, strand)  # 1-based source labels
                for (src, start, end, strand) in segs]
        leaf_genomes[nodes[leaf_idx]] = _FastNucGenome(objs)

    blocks = _build_blocks(leaf_genomes, root_length)
    return _FastNucleotideResult(
        species_tree=species_tree, leaf_genomes=leaf_genomes,
        event_log=EventLog(), registry=SegmentRegistry(), blocks=blocks, root_length=root_length,
    )
