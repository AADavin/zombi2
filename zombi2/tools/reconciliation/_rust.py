"""Optional Rust fast-path for the dated ALElite likelihood.

Mirrors the pure-Python :mod:`zombi2.tools.reconciliation.dated` engine bit-for-bit (same slicing, same
backward-Euler sweep). The Python engine stays the reference; this is only invoked when the
compiled ``zombi2_core`` extension is present. Kept separate so ALElite stays liftable.

Build the extension with::

    cd rust && maturin build --release -i python3
    pip install --force-reinstall rust/target/wheels/*.whl
"""

from __future__ import annotations

try:  # optional native extension (shared with the simulator)
    import zombi2_core as _core

    _HAS = hasattr(_core, "dated_joint_loglik")
    _HAS_UNDATED = hasattr(_core, "undated_joint_loglik")
    _HAS_FAMILY = hasattr(_core, "dated_family_loglik") and hasattr(_core, "undated_family_loglik")
except ImportError:  # pragma: no cover - depends on whether the wheel is built
    _core = None
    _HAS = False
    _HAS_UNDATED = False
    _HAS_FAMILY = False

_ORIG = {"root": 0, "uniform": 1}
_TRANSFERS = {"global": 0, "dated": 1}


def available() -> bool:
    """True if the compiled dated-likelihood kernel is importable."""
    return _HAS


def available_undated() -> bool:
    """True if the compiled undated/reldated-likelihood kernel is importable."""
    return _HAS_UNDATED


def available_family() -> bool:
    """True if the compiled per-family (per-tree) kernels are importable."""
    return _HAS_FAMILY


def _flat_species(sp):
    parent = [b.parent if b.parent is not None else -1 for b in sp.branches]
    left = [b.left if b.left is not None else -1 for b in sp.branches]
    right = [b.right if b.right is not None else -1 for b in sp.branches]
    return (parent, left, right,
            [b.is_leaf for b in sp.branches],
            [b.time for b in sp.branches],
            [b.parent_time for b in sp.branches],
            sp.root)


def _flat_genes(gene_trees, sp):
    offsets = [0]
    is_leaf: list[bool] = []
    left: list[int] = []
    right: list[int] = []
    species: list[int] = []
    for gt in gene_trees:
        for g in gt.nodes:
            is_leaf.append(g.is_leaf)
            left.append(g.left if g.left is not None else -1)
            right.append(g.right if g.right is not None else -1)
            if g.is_leaf:
                si = sp.leaf_index.get(g.species)
                if si is None:
                    raise KeyError(f"gene tip species {g.species!r} is not a species-tree leaf")
                species.append(si)
            else:
                species.append(-1)
        offsets.append(len(is_leaf))
    return offsets, is_leaf, left, right, species


def dated_joint_loglik(gene_trees, sp, dup, transfer, loss, origination, n_extinct, n_steps):
    """Call the Rust dated kernel. Assumes :func:`available` is True."""
    sp_parent, sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root = _flat_species(sp)
    gt_off, gt_leaf, gt_left, gt_right, gt_species = _flat_genes(gene_trees, sp)
    return _core.dated_joint_loglik(
        sp_parent, sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root,
        gt_off, gt_leaf, gt_left, gt_right, gt_species,
        float(dup), float(transfer), float(loss), int(n_steps),
        _ORIG[origination], int(n_extinct),
    )


def undated_joint_loglik(gene_trees, sp, dup, transfer, loss, origination, transfers, n_extinct):
    """Call the Rust undated/reldated kernel. Assumes :func:`available_undated` is True."""
    _, sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root = _flat_species(sp)
    gt_off, gt_leaf, gt_left, gt_right, gt_species = _flat_genes(gene_trees, sp)
    return _core.undated_joint_loglik(
        sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root,
        gt_off, gt_leaf, gt_left, gt_right, gt_species,
        float(dup), float(transfer), float(loss),
        _TRANSFERS[transfers], _ORIG[origination], int(n_extinct),
    )


def dated_family_loglik(gene_trees, sp, dup, transfer, loss, origination, n_steps):
    """Per-family dated log-liks (one per tree). Assumes :func:`available_family` is True."""
    sp_parent, sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root = _flat_species(sp)
    gt_off, gt_leaf, gt_left, gt_right, gt_species = _flat_genes(gene_trees, sp)
    return _core.dated_family_loglik(
        sp_parent, sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root,
        gt_off, gt_leaf, gt_left, gt_right, gt_species,
        float(dup), float(transfer), float(loss), int(n_steps), _ORIG[origination],
    )


def undated_family_loglik(gene_trees, sp, dup, transfer, loss, origination, transfers):
    """Per-family undated/reldated log-liks (one per tree). Assumes :func:`available_family`."""
    _, sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root = _flat_species(sp)
    gt_off, gt_leaf, gt_left, gt_right, gt_species = _flat_genes(gene_trees, sp)
    return _core.undated_family_loglik(
        sp_left, sp_right, sp_is_leaf, sp_time, sp_ptime, sp_root,
        gt_off, gt_leaf, gt_left, gt_right, gt_species,
        float(dup), float(transfer), float(loss), _TRANSFERS[transfers], _ORIG[origination],
    )
