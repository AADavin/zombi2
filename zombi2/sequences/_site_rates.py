"""Across-site rate variation — the discretised Gamma (Yang 1994).

Real genes do not evolve at one speed along their length: some positions are held nearly fixed by
function and others race ahead. The standard way to say so is a **Gamma distribution of rates across
sites** — every site draws a multiplier of the substitution rate from a Gamma with mean 1, so one
parameter (its **shape**) sets how unequal the sites are, and a small shape means a few fast sites
among many slow ones.

**Why it is discretised.** Drawing a continuous rate per site would give every site its own branch
length, and a branch length is what the engine's transition matrix is computed from: one matrix
exponential per site per branch, where the whole level's speed rests on `_cdf_for()`'s
cache of a few hundred matrices reused across millions of sites. Yang's (1994) discretisation cuts
the Gamma into ``categories`` **equal-probability** classes and represents each by its mean, so a
site takes one of a handful of rates and the sites sharing a class share a branch length — and
therefore a cached matrix. Four classes is the field's default and is what this defaults to.

**The mean-1 contract.** The classes are normalised so the mean rate over all sites is exactly 1.
That is not cosmetic: it is what keeps a branch length in the phylograms equal to the **mean**
substitutions per site, so a run with rate variation and a run without have the same tree at the same
rate and differ only in how the change is spread across columns. `SubstitutionModel.across_sites()`
preserves it when an invariant class is added, by scaling the Gamma classes up to make room for the
sites that never change.

Everything here is `math`/`numpy` only. The discretisation needs the regularised lower incomplete
gamma and its inverse, which numpy does not carry and which scipy would — but scipy is a
**test-only** dependency here, never a runtime one, and `substitution_models` states
that as a promise. So the two functions are written out, following the standard series /
continued-fraction split; they run once per model, never per site, so the implementation is chosen
for being obviously right rather than for being fast. ``tests/test_sequences_gamma.py`` checks them
against scipy and against the category values published for PAML.
"""

from __future__ import annotations

import math

#: convergence controls for `_incomplete_gamma` — it is called a handful of times per model, so the
#: iteration caps are generous and the tolerance is at the edge of double precision
_ITMAX = 1000
_EPS = 3.0e-16
_FPMIN = 1.0e-300


def _incomplete_gamma(a: float, x: float) -> float:
    """The regularised lower incomplete gamma ``P(a, x) = γ(a, x) / Γ(a)`` — the CDF of a
    ``Gamma(shape=a, rate=1)`` at ``x``.

    Two expansions, split where each converges: the series for ``x < a + 1`` and the continued
    fraction (which computes ``Q = 1 - P``) beyond it. This is the standard treatment; the split
    point is where the series starts needing more terms than the fraction.
    """
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        # series: P(a,x) = e^{-x} x^a / Γ(a) · Σ_{n≥0} x^n / (a(a+1)…(a+n))
        ap, term, total = a, 1.0 / a, 1.0 / a
        for _ in range(_ITMAX):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * _EPS:
                break
        return total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    # continued fraction (modified Lentz) for Q(a,x); P = 1 - Q
    b, c = x + 1.0 - a, 1.0 / _FPMIN
    d = 1.0 / b
    h = d
    for i in range(1, _ITMAX + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _FPMIN:
            d = _FPMIN
        c = b + an / c
        if abs(c) < _FPMIN:
            c = _FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    return 1.0 - math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def _gamma_quantile(p: float, a: float) -> float:
    """The inverse of `_incomplete_gamma` in ``x``: the ``p``-quantile of ``Gamma(shape=a, rate=1)``.

    Bracket-and-bisect. It runs ``categories - 1`` times per model and never per site, so there is no
    reason to reach for a Newton step that would need guarding against the shape parameters where
    the density is unbounded at zero (``a < 1``, which is exactly the interesting case here).
    """
    if not 0.0 < p < 1.0:
        raise ValueError(f"a quantile needs a probability strictly inside (0, 1), got {p!r}")
    hi = max(a, 1.0)
    while _incomplete_gamma(a, hi) < p:      # push the upper bracket out until it straddles p
        hi *= 2.0
    lo = 0.0
    for _ in range(200):                     # 200 halvings takes the bracket below any float gap
        mid = 0.5 * (lo + hi)
        if _incomplete_gamma(a, mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def discrete_gamma(shape: float, categories: int) -> tuple[float, ...]:
    """The ``categories`` rate classes of a mean-1 Gamma of the given ``shape`` (Yang 1994).

    The Gamma is cut into equal-probability slices and each is represented by **its own mean**, so
    the classes are unequally spaced (crowded near zero for a small shape, where most of the
    probability sits) while every class holds the same share of sites. Returned in ascending order,
    renormalised so the mean over classes is exactly 1.0 — the contract the phylograms rest on.

    For ``X ~ Gamma(a, rate=a)`` the partial mean is ``∫₀^q x f(x) dx = P(a+1, a·q)``, so a class
    bounded by consecutive quantiles has mean ``K · [P(a+1, xᵢ) − P(a+1, xᵢ₋₁)]`` where ``xᵢ`` is the
    ``i/K`` quantile of ``Gamma(a, rate=1)`` — which is ``a·qᵢ``, so the rate parameter never has to
    be carried around. The last class's upper bound is infinity, where ``P`` is 1.
    """
    if isinstance(shape, bool) or not isinstance(shape, (int, float)):
        raise TypeError(f"gamma shape must be a real number, got {shape!r}")
    if not math.isfinite(shape) or shape <= 0:
        raise ValueError(f"gamma shape must be finite and positive, got {shape!r}")
    if isinstance(categories, bool) or not isinstance(categories, int) or categories < 2:
        raise ValueError(f"a Gamma needs at least 2 rate categories to vary anything, got {categories!r}")
    a = float(shape)
    # the cumulative partial means at each class boundary; 0 at the bottom, 1 at the top
    edges = [0.0]
    for i in range(1, categories):
        edges.append(_incomplete_gamma(a + 1.0, _gamma_quantile(i / categories, a)))
    edges.append(1.0)
    rates = [categories * (edges[i + 1] - edges[i]) for i in range(categories)]
    mean = sum(rates) / categories
    return tuple(r / mean for r in rates)     # pin the mean to 1 against the quantile's round-off


__all__ = ["discrete_gamma"]
