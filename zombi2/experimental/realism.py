"""Realism metric (**P4**): Fréchet-ESM distance between two protein sets in a PLM's embedding space.

An FID-style distance -- ``||mu_a - mu_b||^2 + tr(Sigma_a + Sigma_b - 2 (Sigma_a Sigma_b)^{1/2})`` over
mean-pooled ESM2 embeddings. It is ~0 for a set versus itself and grows the more the two distributions
differ. The intended use is to check that *simulated* proteins occupy the same region of embedding space
as *real* ones -- neutral drift should score far worse than selection, quantifying "do the simulated
sequences look real". Needs a critic implementing ``embed`` (e.g. :class:`~zombi2.experimental.selection.ESM2Critic`)
and scipy; both are imported lazily.
"""
from __future__ import annotations

import numpy as np

from zombi2.experimental import warn_experimental
from zombi2.experimental.selection import Critic

__all__ = ["frechet_esm_distance"]


def frechet_esm_distance(seqs_a, seqs_b, critic: Critic) -> float:
    """Fréchet distance between the ``critic`` embeddings of two protein sets (lower = more similar;
    ~0 for a set versus itself)."""
    warn_experimental("frechet_esm_distance")
    try:
        from scipy.linalg import sqrtm
    except ImportError as exc:  # pragma: no cover - only without the extra installed
        raise ImportError("frechet_esm_distance needs scipy; pip install 'zombi2[selection]'") from exc
    A = np.asarray(critic.embed(list(seqs_a)), dtype=float)
    B = np.asarray(critic.embed(list(seqs_b)), dtype=float)
    if A.ndim != 2 or B.ndim != 2 or A.shape[1] == 0 or A.shape[1] != B.shape[1]:
        raise ValueError("critic.embed must return matching (n, dim) arrays with dim > 0")
    if len(A) < 2 or len(B) < 2:
        raise ValueError("need at least 2 sequences per set to estimate a covariance")
    dmu = A.mean(0) - B.mean(0)
    Sa = np.atleast_2d(np.cov(A, rowvar=False))
    Sb = np.atleast_2d(np.cov(B, rowvar=False))
    d = Sa.shape[0]
    ridge = 1e-10 * (np.trace(Sa) + np.trace(Sb)) / (2.0 * d)     # keep Sa @ Sb strictly PSD -> real sqrt
    Sa = Sa + ridge * np.eye(d)
    Sb = Sb + ridge * np.eye(d)
    cov = sqrtm(Sa @ Sb)
    if np.iscomplexobj(cov):
        imag = float(np.abs(cov.imag).max())
        if imag > 1e-6 * (1.0 + float(np.abs(cov.real).max())):
            raise ValueError(f"Fréchet sqrtm has a large imaginary part ({imag:.2e}); the embedding "
                             "covariances look degenerate (too few or collinear sequences?)")
        cov = cov.real
    return float(dmu @ dmu + np.trace(Sa + Sb - 2.0 * cov))
