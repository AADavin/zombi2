"""Fit gene-family rates to an empirical profile table by rejection ABC.

Given a fixed species tree and an **empirical copy-number profile** (families x extant
species), :func:`match_profiles` searches for rates under which ZOMBI2's gene-family model
reproduces that profile. It is a plain **Approximate Bayesian Computation** (ABC) rejection
sampler:

1. draw a parameter set from the priors,
2. simulate a profile matrix under it (the built-in model runs on the fast Rust engine),
3. reduce both the simulated and the empirical matrix to a vector of **summary
   statistics**, and
4. keep the draws whose summaries land closest to the empirical one.

You cannot match a profile *table* directly — a simulation produces a different set of
families, in a different order, with different labels. So matching is done on
**permutation-invariant** summaries of the whole matrix (see :func:`default_summary`):

* the **gene frequency spectrum** — how many families are present in exactly 1, 2, ... S
  species (the pangenome core/shell/cloud curve; the most informative single statistic);
* the per-species **genome sizes** (total copies per species); and
* the **copy-number spectrum** — how many present cells hold 1, 2, 3, >=4 copies (this is
  what separates duplication from transfer).

Each summary component is scaled by its standard deviation across the simulated batch, so
the distance is scale-free without hand-tuned weights (Prangle 2015).

**Models.** ``model="uniform"`` (default) fits the four scalar D/T/L/O rates shared by
every family (the built-in model — runs on the fast Rust engine). ``model="family"`` fits
the same four values as the *means* of per-family rate distributions
(:class:`~zombi2.FamilySampledRates`), so families are heterogeneous; this runs on the
Python engine. You can also pass any callable ``params_dict -> RateModel`` as ``model`` for
full generality.

**On identifiability.** From presence/copy number alone the rates are only partly
separable — gain routes (origination vs transfer) and the gain/loss balance trade off. So
the result is a *posterior sample* (a ridge), not a point estimate; read ``.summary()``
credible intervals and ``.plot_spectra()``, not ``.best`` alone.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .distributions import Distribution, Fixed, Gamma, Uniform, as_distribution
from .events import EventType
from .profiles import ProfileMatrix
from .rates import FamilySampledRates
from .tree import Tree

#: The rate parameters the built-in models fit (the independent-family model).
RATE_PARAMS = ("duplication", "transfer", "loss", "origination")


def _wquantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Weighted quantile of ``values`` (weights need not be normalised)."""
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w) - 0.5 * w
    cum /= w.sum()
    return float(np.interp(q, cum, v))


# --- summary statistics ----------------------------------------------------------

def frequency_spectrum(pm: ProfileMatrix, n_species: int) -> np.ndarray:
    """Counts of families present in exactly ``k`` species, for ``k = 1 .. n_species``.

    Computed off the sparse profile (per-family presence counts), never the dense array."""
    if n_species == 0 or pm.shape[0] == 0:
        return np.zeros(n_species)
    present = pm.presence_per_family()
    return np.bincount(present, minlength=n_species + 1)[1:n_species + 1].astype(float)


def genome_sizes(pm: ProfileMatrix, species_order: list[str]) -> np.ndarray:
    """Total copy number per species, aligned to ``species_order`` (missing species -> 0)."""
    col_sum = pm.copies_per_species()
    by_name = {s: float(col_sum[j]) for j, s in enumerate(pm.species)}
    return np.array([by_name.get(s, 0.0) for s in species_order])


def copy_number_spectrum(pm: ProfileMatrix, max_copies: int = 4) -> np.ndarray:
    """Counts of present cells with copy number 1, 2, ..., ``>=max_copies`` (off the
    sparse non-zero values, no dense materialisation)."""
    out = np.zeros(max_copies)
    vals = pm.copy_values()
    if vals.size == 0:
        return out
    for c in range(1, max_copies):
        out[c - 1] = np.count_nonzero(vals == c)
    out[max_copies - 1] = np.count_nonzero(vals >= max_copies)
    return out


class _DefaultSummary:
    """The default summary function as a picklable object (so it survives to workers).

    Concatenates the frequency spectrum, per-species genome sizes, and copy-number
    spectrum, all computed on the empirical species axes it is built with. The frequency
    spectrum occupies the leading ``len(species_order)`` entries — which ``ABCFit``'s
    spectrum diagnostics rely on.
    """

    def __init__(self, species_order: list[str], max_copies: int = 4):
        self.species_order = list(species_order)
        self.n_species = len(species_order)
        self.max_copies = max_copies

    def __call__(self, pm: ProfileMatrix) -> np.ndarray:
        return np.concatenate([
            frequency_spectrum(pm, self.n_species),
            genome_sizes(pm, self.species_order),
            copy_number_spectrum(pm, self.max_copies),
        ])


def default_summary(species_order: list[str], max_copies: int = 4) -> _DefaultSummary:
    """The default summary function: frequency spectrum + genome sizes + copy spectrum.

    Returns a picklable callable ``pm -> 1-D float vector`` closed over the empirical
    species order, so simulated matrices are summarised on exactly the same axes as the
    empirical one.
    """
    return _DefaultSummary(species_order, max_copies)


def event_count_summary(genomes) -> np.ndarray:
    """Total number of duplication, transfer, and loss events in a :class:`~zombi2.Genomes`.

    Gene-tree-derived counts: they carry information the copy-number profile alone cannot,
    and pin the *gain*-side rates (duplication/transfer/origination) sharply. (Loss stays
    the hardest rate to identify — a fully lost lineage leaves little observable trace.)
    """
    dup = tr = ls = 0
    for r in genomes.event_log:
        if r.event is EventType.DUPLICATION:
            dup += 1
        elif r.event is EventType.TRANSFER:
            tr += 1
        elif r.event is EventType.LOSS:
            ls += 1
    return np.array([dup, tr, ls], dtype=float)


class _GeneTreeSummary:
    """Default summary when gene trees are available: the profile summary followed by the
    three event counts, with per-feature weights that give the (tiny) event block the same
    total weight as the (much larger) profile block — otherwise the counts are drowned out.

    Consumes a :class:`~zombi2.Genomes` (needs the event log). Picklable. The frequency
    spectrum still leads the vector, so the spectrum diagnostics keep working.
    """

    def __init__(self, species_order: list[str], max_copies: int = 4, event_weight=None):
        self.profile = _DefaultSummary(species_order, max_copies)
        self._pw = 2 * self.profile.n_species + max_copies      # profile-summary length
        # balance the two blocks: 3 * w^2 == profile length  ->  each block contributes equally
        w = (self._pw / 3.0) ** 0.5 if event_weight is None else float(event_weight)
        self.event_weight = w
        self.feature_weights = np.concatenate([np.ones(self._pw), np.full(3, w)])

    def __call__(self, genomes) -> np.ndarray:
        return np.concatenate([self.profile(genomes.profiles), event_count_summary(genomes)])


def default_gene_tree_summary(species_order: list[str], max_copies: int = 4, event_weight=None):
    """Picklable summary combining the profile summary with weighted event counts.

    ``event_weight`` defaults to a value that balances the event block against the profile
    block; pass a float to override (larger = more weight on the gene-tree counts).
    """
    return _GeneTreeSummary(species_order, max_copies, event_weight)


def cooccurrence_features(pm: ProfileMatrix, *, threshold: float = 0.35) -> np.ndarray:
    """Module structure of the gene-family co-occurrence graph — the Pellegrini signal.

    Builds the family x family presence-correlation matrix (over species), thresholds it into
    an undirected graph (an edge where two families co-occur with correlation above
    ``threshold``), and returns three graph statistics that detect *clusters* of co-occurring
    families (functional modules): ``[edge_count, triangle_count, transitivity]``.

    Unlike the mean or the tail of the pairwise correlations, the triangle / transitivity
    terms isolate *cliques* — coupling makes families that co-occur also co-occur with each
    other (closing triangles), whereas the spurious pairwise correlations left by noise or
    shared ancestry stay scattered and rarely close a triangle. Blind to family identity /
    order (permutation-invariant) and fixed-length, so it plugs into ABC as a summary block.
    Families present in every species or none (no variance) are dropped.
    """
    P = (pm.matrix > 0).astype(float)
    keep = P.std(axis=1) > 0
    if int(keep.sum()) < 3:
        return np.zeros(3)
    C = np.corrcoef(P[keep])
    A = (C > threshold).astype(float)
    np.fill_diagonal(A, 0.0)
    edges = A.sum() / 2.0
    triangles = float(np.trace(A @ A @ A)) / 6.0
    deg = A.sum(axis=1)
    triples = float(np.sum(deg * (deg - 1))) / 2.0
    transitivity = (3.0 * triangles / triples) if triples > 0 else 0.0
    return np.array([edges, triangles, transitivity])


class _CooccurrenceSummary:
    """The default marginal summary followed by the co-occurrence module features.

    Concatenates :class:`_DefaultSummary` (frequency spectrum + genome sizes + copy spectrum)
    with :func:`cooccurrence_features`, and carries ``feature_weights`` that give the small
    (3-feature) co-occurrence block the same total weight as the large profile block —
    otherwise it is drowned out in the distance. The frequency spectrum still leads the vector,
    so ``ABCFit``'s spectrum diagnostics keep working. Picklable.
    """

    def __init__(self, species_order: list[str], max_copies: int = 4, threshold: float = 0.35):
        self.profile = _DefaultSummary(species_order, max_copies)
        self.threshold = threshold
        self._pw = 2 * self.profile.n_species + max_copies      # profile-block length
        w = (self._pw / 3.0) ** 0.5                             # balance vs the 3 co-occ features
        self.feature_weights = np.concatenate([np.ones(self._pw), np.full(3, w)])

    def __call__(self, pm: ProfileMatrix) -> np.ndarray:
        return np.concatenate([self.profile(pm),
                               cooccurrence_features(pm, threshold=self.threshold)])


def cooccurrence_summary(species_order: list[str], max_copies: int = 4, threshold: float = 0.35):
    """Picklable summary = the default marginal summary + co-occurrence module features.

    Use with ``match_profiles(..., statistics=cooccurrence_summary(species_order))`` or
    :func:`match_coupled` to make ABC see gene-family *co-occurrence* (module structure), not
    just marginal prevalence — needed to fit coupling / non-independence, where the marginal
    alone can be blind. See :func:`cooccurrence_features`.
    """
    return _CooccurrenceSummary(species_order, max_copies, threshold)


# --- rate models to fit ----------------------------------------------------------

class _FamilyModel:
    """Build a :class:`~zombi2.FamilySampledRates` from per-family-rate **means**.

    Each of duplication/transfer/loss becomes a ``Gamma(shape, mean/shape)`` per-family
    distribution (its mean is the fitted value; ``shape`` sets the dispersion — larger is
    tighter, ``shape -> inf`` approaches a single shared rate). Origination stays a
    per-branch scalar, as in ``FamilySampledRates``. Picklable, so it survives to workers.
    """

    def __init__(self, shape: float = 2.0):
        if shape <= 0:
            raise ValueError("family_shape must be > 0")
        self.shape = float(shape)

    def _dist(self, mean: float):
        return Gamma(self.shape, mean / self.shape) if mean > 0 else Fixed(0.0)

    def __call__(self, params: dict) -> FamilySampledRates:
        return FamilySampledRates(
            duplication=self._dist(params.get("duplication", 0.0)),
            transfer=self._dist(params.get("transfer", 0.0)),
            loss=self._dist(params.get("loss", 0.0)),
            origination=params.get("origination", 0.0),
        )


def _resolve_model(model, family_shape):
    """Return the ``builder``: ``None`` for the uniform scalar model (built-in, runs on Rust
    via ``simulate_genomes``), otherwise a callable ``params -> RateModel`` (Python engine)."""
    if model is None or model == "uniform":
        return None
    if model == "family":
        return _FamilyModel(family_shape)
    if callable(model):
        return model
    raise ValueError("model must be 'uniform', 'family', or a callable params -> RateModel")


# --- the fit result --------------------------------------------------------------

@dataclass
class ABCFit:
    """Result of :func:`match_profiles` — the accepted parameter posterior and diagnostics."""

    param_names: list[str]
    samples: np.ndarray            # (n_sims, n_params) every drawn parameter set
    distances: np.ndarray          # (n_sims,) summary-space distance to the empirical target
    accepted: np.ndarray           # (n_accept,) indices into ``samples`` that were kept
    tolerance: float               # accepted distance threshold (largest accepted distance)
    empirical_summary: np.ndarray
    priors: dict[str, Distribution]
    accepted_summaries: np.ndarray  # (n_accept, summary_len) summaries of the kept draws
    summary_sd: np.ndarray          # (summary_len,) across-batch std used to scale the distance
    n_species: int
    empirical: ProfileMatrix
    uses_default_summary: bool
    sample_weights: np.ndarray | None = None  # per-sample importance weights (SMC); None = uniform

    def __post_init__(self):
        self._adjusted = None          # cache for regression_adjust()
        self._adjust_weights = None

    @property
    def posterior(self) -> dict[str, np.ndarray]:
        """Accepted values per parameter (the ABC posterior sample)."""
        return {n: self.samples[self.accepted, i] for i, n in enumerate(self.param_names)}

    @property
    def best(self) -> dict[str, float]:
        """The single closest-matching parameter set (one point on the ridge — see ``.summary``)."""
        i = int(np.argmin(self.distances))
        return {n: float(self.samples[i, j]) for j, n in enumerate(self.param_names)}

    def summary(self, adjusted: bool = False) -> dict[str, dict[str, float]]:
        """Per-parameter posterior mean, median, and 95% credible interval.

        With ``adjusted=True`` the summary is computed from the regression-adjusted
        posterior (:meth:`regression_adjust`), using its kernel weights.
        """
        if adjusted:
            post = self.regression_adjust()
            w = self._adjust_weights
            return {name: {
                "mean": float(np.average(v, weights=w)),
                "median": _wquantile(v, w, 0.5),
                "lo95": _wquantile(v, w, 0.025),
                "hi95": _wquantile(v, w, 0.975),
            } for name, v in post.items()}
        w = None if self.sample_weights is None else self.sample_weights[self.accepted]
        out: dict[str, dict[str, float]] = {}
        for name, vals in self.posterior.items():
            if w is None:
                out[name] = {
                    "mean": float(vals.mean()),
                    "median": float(np.median(vals)),
                    "lo95": float(np.quantile(vals, 0.025)),
                    "hi95": float(np.quantile(vals, 0.975)),
                }
            else:
                out[name] = {
                    "mean": float(np.average(vals, weights=w)),
                    "median": _wquantile(vals, w, 0.5),
                    "lo95": _wquantile(vals, w, 0.025),
                    "hi95": _wquantile(vals, w, 0.975),
                }
        return out

    def regression_adjust(self, ridge: float = 1.0) -> dict[str, np.ndarray]:
        """Local-linear **regression adjustment** of the accepted posterior (Beaumont 2002).

        Rejection ABC keeps draws whose summary is merely *close* to the empirical target;
        any residual dependence of the parameters on that summary discrepancy biases the
        posterior (this is what pulls loss up the ridge). This corrects it *post hoc, with
        no new simulations*: it regresses the accepted parameters on their (scaled) summary
        residuals ``s_i - s*`` — weighted by an Epanechnikov kernel of the distance — and
        subtracts the fitted slope, i.e. projects every accepted draw to what it would have
        been had its summary hit the target exactly. Returns ``{param: adjusted values}``
        (clipped at 0, since rates are non-negative) and caches the result.

        ``ridge`` regularises the regression: larger values shrink the correction toward
        the raw posterior (safe but mild), smaller values correct more aggressively but can
        overfit and narrow the intervals when the summary is high-dimensional relative to
        the number of accepted draws. The correction can only exploit information the
        summaries actually carry — it sharpens well-identified rates but cannot rescue one
        that sits on a ridge.
        """
        if self._adjusted is not None:
            return self._adjusted
        res = (self.accepted_summaries - self.empirical_summary) / self.summary_sd  # (k, L)
        # Drop summary columns that don't vary across the accepted set (e.g. the all-zero
        # tail bins of the frequency spectrum): they carry no information and make the
        # normal equations singular / invite overfitting.
        keep = res.std(axis=0) > 1e-9
        res = res[:, keep]
        y = self.samples[self.accepted]                                             # (k, P)
        d = self.distances[self.accepted]
        dmax = d.max() if d.max() > 0 else 1.0
        w = np.clip(1.0 - (d / dmax) ** 2, 1e-12, None)                             # Epanechnikov

        x = np.hstack([np.ones((len(res), 1)), res])                               # intercept + residuals
        sw = np.sqrt(w)[:, None]
        xw, yw = x * sw, y * sw
        if ridge > 0:
            reg = ridge * np.eye(x.shape[1]); reg[0, 0] = 0.0                      # don't penalise intercept
            beta = np.linalg.solve(xw.T @ xw + reg, xw.T @ yw)                     # (kept+1, P)
        else:
            beta = np.linalg.lstsq(xw, yw, rcond=None)[0]                          # min-norm (rank-safe)
        adjusted = np.clip(y - res @ beta[1:], 0.0, None)                          # remove slope, keep intercept

        self._adjust_weights = w
        self._adjusted = {n: adjusted[:, j] for j, n in enumerate(self.param_names)}
        return self._adjusted

    # --- diagnostics ---------------------------------------------------------
    def spectra_data(self) -> dict[str, np.ndarray]:
        """Data for the frequency-spectrum posterior-predictive check.

        Returns ``{"k", "empirical", "accepted"}``: the species counts ``k = 1..S``, the
        empirical gene frequency spectrum, and the ``(n_accept, S)`` spectra of the accepted
        simulations. Requires the default summary (which places the spectrum first).
        """
        if not self.uses_default_summary:
            raise ValueError(
                "spectrum diagnostics assume the default summary (statistics=None); "
                "with a custom summary, compute frequency_spectrum() yourself."
            )
        s = self.n_species
        return {
            "k": np.arange(1, s + 1),
            "empirical": self.empirical_summary[:s],
            "accepted": self.accepted_summaries[:, :s],
        }

    def plot_spectra(self, ax=None, *, band=(2.5, 97.5), draws: bool = False):
        """Plot the empirical gene frequency spectrum against the accepted simulations.

        Overlays the empirical spectrum (points) on the accepted simulations' spectra (a
        median line and a percentile band). If the empirical curve sits inside the band the
        fitted model can reproduce the data; if it falls outside, no setting of these rates
        reproduces the profile — a sign the model is missing something. Needs matplotlib;
        pass an existing ``ax`` or one is created. Returns the axes.
        """
        import matplotlib.pyplot as plt  # lazy — matplotlib is not a hard dependency

        d = self.spectra_data()
        k, emp, acc = d["k"], d["empirical"], d["accepted"]
        if ax is None:
            _, ax = plt.subplots()
        lo, med, hi = np.percentile(acc, [band[0], 50, band[1]], axis=0)
        ax.fill_between(k, lo, hi, alpha=0.25, color="C0",
                        label=f"accepted {band[0]:g}–{band[1]:g}%")
        ax.plot(k, med, lw=1.5, color="C0", label="accepted median")
        if draws:
            for row in acc:
                ax.plot(k, row, color="gray", alpha=0.12, lw=0.5)
        ax.plot(k, emp, "o-", color="k", label="empirical")
        ax.set_xlabel("number of species a family is present in")
        ax.set_ylabel("number of gene families")
        ax.set_title("Gene frequency spectrum: empirical vs accepted")
        ax.legend()
        return ax

    def __repr__(self) -> str:
        parts = [f"{n}={s['median']:.3g} [{s['lo95']:.3g}, {s['hi95']:.3g}]"
                 for n, s in self.summary().items()]
        return (f"ABCFit(n_accept={len(self.accepted)}, tol={self.tolerance:.3g}, "
                f"{', '.join(parts)})")


# --- driver ----------------------------------------------------------------------

def _normalize_priors(priors: dict, restrict: bool) -> dict[str, Distribution]:
    out: dict[str, Distribution] = {}
    for name, spec in priors.items():
        if restrict and name not in RATE_PARAMS:
            raise ValueError(f"unknown parameter {name!r}; the built-in models fit {RATE_PARAMS}")
        if isinstance(spec, (tuple, list)) and len(spec) == 2:
            out[name] = Uniform(*spec)          # (low, high) shorthand
        else:
            out[name] = as_distribution(spec)   # Distribution / float / scipy / callable
    if not out:
        raise ValueError("provide a prior for at least one parameter")
    return out


def _simulate(tree, vals, *, builder, initial_size, max_family_size, transfers,
              seed, return_genomes=False):
    """One ABC simulation. ``simulate_genomes`` picks the engine automatically: the uniform
    scalar model (``builder is None``) runs on Rust; family/custom builders run on Python.
    ``output="profiles"`` takes the fast counts-only path; ``output="genomes"`` returns the
    full result (needed for the gene-tree summary)."""
    from .simulation import simulate_genomes

    kw = dict(initial_families=initial_size, max_family_size=max_family_size,
              transfers=transfers, seed=seed,
              output="genomes" if return_genomes else "profiles")
    if builder is None:
        return simulate_genomes(tree, **vals, **kw)
    return simulate_genomes(tree, builder(vals), **kw)


def _simulate_and_summarize(draw, tree, builder, initial_size,
                            max_family_size, transfers, summarize, return_genomes):
    vals, seed = draw
    out = _simulate(tree, vals, builder=builder, initial_size=initial_size,
                    max_family_size=max_family_size, transfers=transfers, seed=seed,
                    return_genomes=return_genomes)
    return summarize(out)


# Per-worker cache of the fixed simulation config (set once per process via the pool
# initializer, so only the tiny per-draw payload crosses the process boundary each task).
_WORKER_CFG: tuple = ()


def _worker_init(*cfg) -> None:
    global _WORKER_CFG
    _WORKER_CFG = cfg


def _worker_run(draw) -> np.ndarray:
    return _simulate_and_summarize(draw, *_WORKER_CFG)


def match_profiles(
    tree: Tree,
    empirical,
    priors: dict,
    *,
    model=None,
    family_shape: float = 2.0,
    statistics: Callable[[ProfileMatrix], np.ndarray] | None = None,
    gene_trees: bool = False,
    feature_weights=None,
    n_sims: int = 1000,
    accept=0.05,
    initial_families: int = 20,
    max_family_size=None,
    transfers=None,
    processes: int | None = None,
    seed: int | None = None,
) -> ABCFit:
    """Fit gene-family rates to an empirical profile by rejection ABC.

    ``initial_families`` is the number of gene families seeded at the root of each
    simulation (default 20).

    Parameters
    ----------
    tree:
        The species tree the profiles were observed on (fixed; branch lengths set the rate
        scale). The same tree is used for every simulation.
    empirical:
        The target :class:`~zombi2.ProfileMatrix`, or a path / TSV text loadable by
        :meth:`ProfileMatrix.from_tsv`. Its species become the axes every simulation is
        summarised on.
    priors:
        ``{param: prior}`` over any of ``duplication/transfer/loss/origination``. A prior is
        a :class:`~zombi2.Distribution`, a bare float (fixed), a ``(low, high)`` tuple
        (uniform), or any scipy frozen dist / ``rng -> float`` callable. Parameters omitted
        here are held at 0.
    model:
        ``"uniform"`` / ``None`` (default) fits the four scalar rates shared by all families
        (the built-in model — runs on the Rust engine, so it is fast). ``"family"`` fits them
        as the *means* of per-family rate distributions (:class:`~zombi2.FamilySampledRates`,
        dispersion set by ``family_shape``; Python engine). Or pass a callable
        ``params -> RateModel``.
    statistics:
        A summary function ``pm -> 1-D array``; defaults to :func:`default_summary` over the
        empirical species. Must be picklable if ``processes`` is used. May expose a
        ``feature_weights`` attribute to weight its components in the distance.
    gene_trees:
        If True, use gene-tree information: ``empirical`` must be a :class:`~zombi2.Genomes`
        (not just a profile), the default summary becomes :func:`default_gene_tree_summary`
        (profile + weighted duplication/transfer/loss counts), and simulations return the full
        genealogy (``output="genomes"``) rather than the counts-only path. Sharpens the
        gain-side rates; loss stays the hardest to identify.
    feature_weights:
        Optional per-summary-component weights applied to the (scaled) distance, so a small
        informative block of statistics is not drowned by a large one. Overrides any
        ``feature_weights`` attribute on the summary.
    n_sims:
        Number of prior draws to simulate.
    accept:
        A float in ``(0, 1]`` keeps that fraction of the closest draws; an int keeps that
        many.
    processes:
        ``None`` (default) or ``1`` runs serially in-process. An int ``> 1`` distributes the
        simulations across that many worker processes (results are identical regardless of
        the count). As with any multiprocessing, call from a ``__main__`` guard, and any
        custom ``model`` / ``statistics`` must be picklable.
    seed:
        Seeds the whole procedure (parameter draws and per-simulation seeds); results are
        reproducible and independent of ``processes``.

    Returns
    -------
    ABCFit
        Holding the accepted posterior sample (``.posterior``/``.summary()``), the single
        closest draw (``.best``), per-draw distances, and the spectrum diagnostic
        (``.plot_spectra()``).
    """
    if gene_trees:
        if not (hasattr(empirical, "profiles") and hasattr(empirical, "event_log")):
            raise TypeError("gene_trees=True requires an empirical Genomes (with gene trees), "
                            "not a bare profile matrix")
        profile_mat = empirical.profiles
        species_order = list(profile_mat.species)
        summarize = statistics or default_gene_tree_summary(species_order)
    else:
        profile_mat = (empirical if isinstance(empirical, ProfileMatrix)
                       else ProfileMatrix.from_tsv(empirical))
        empirical = profile_mat
        species_order = list(profile_mat.species)
        summarize = statistics or default_summary(species_order)
    uses_default = statistics is None
    target = summarize(empirical)

    builder = _resolve_model(model, family_shape)
    priors = _normalize_priors(priors, restrict=builder is None or isinstance(builder, _FamilyModel))
    names = list(priors.keys())
    weights = feature_weights if feature_weights is not None else getattr(summarize, "feature_weights", None)
    rng = np.random.default_rng(seed)

    # Pre-draw every (parameters, sim-seed) up front, so the result is identical whether we
    # then simulate serially or across processes (the draw order is fixed by the master rng).
    draws = []
    samples = np.empty((n_sims, len(names)))
    for i in range(n_sims):
        vals = {n: max(0.0, priors[n].sample(rng)) for n in names}
        sim_seed = int(rng.integers(1, 2**63 - 1))
        draws.append((vals, sim_seed))
        samples[i] = [vals[n] for n in names]

    cfg = (tree, builder, initial_families, max_family_size, transfers, summarize, gene_trees)
    if processes is None or processes == 1:
        summaries = np.array([_simulate_and_summarize(d, *cfg) for d in draws])
    else:
        chunk = max(1, n_sims // (processes * 4))
        with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init,
                                 initargs=cfg) as ex:
            summaries = np.array(list(ex.map(_worker_run, draws, chunksize=chunk)))

    # scale-free distance: normalise each summary component by its across-batch spread, then
    # apply optional per-feature weights (so a small informative block is not drowned out).
    sd = summaries.std(axis=0)
    sd[sd == 0] = 1.0
    resid = (summaries - target) / sd
    if weights is not None:
        resid = resid * np.asarray(weights)
    distances = np.sqrt((resid ** 2).sum(axis=1))

    order = np.argsort(distances, kind="stable")
    if isinstance(accept, float):
        if not 0.0 < accept <= 1.0:
            raise ValueError("float accept must be in (0, 1]")
        k = max(1, round(accept * n_sims))
    else:
        k = min(int(accept), n_sims)
    accepted = np.sort(order[:k])

    return ABCFit(
        param_names=names,
        samples=samples,
        distances=distances,
        accepted=accepted,
        tolerance=float(distances[order[k - 1]]),
        empirical_summary=target,
        priors=priors,
        accepted_summaries=summaries[accepted],
        summary_sd=sd,
        n_species=len(species_order),
        empirical=profile_mat,
        uses_default_summary=uses_default,
    )


def match_coupled(
    tree: Tree,
    empirical,
    spec_builder,
    priors: dict,
    *,
    statistics: Callable[[ProfileMatrix], np.ndarray] | None = None,
    feature_weights=None,
    n_sims: int = 1000,
    accept=0.05,
    transfers=None,
    initial_presence=None,
    seed: int | None = None,
) -> ABCFit:
    """Fit parameters of a **coupled (Potts) gene-family model** to an empirical profile by ABC.

    The likelihood-free counterpart of :func:`match_profiles` for the non-independence model
    (see :func:`~zombi2.simulate_coupled`): the simulator is ``simulate_coupled`` over a
    **fixed family panel**, so it can fit coupling parameters — e.g. the coupling strength
    ``J`` — for which no likelihood exists.

    Parameters
    ----------
    tree:
        The species tree the profiles were observed on (fixed across simulations).
    empirical:
        Target :class:`~zombi2.ProfileMatrix` (or a path / TSV text). Its species set the axes
        every simulation is summarised on and must match the tree's extant leaves.
    spec_builder:
        A callable ``params -> CouplingSpec`` that turns a drawn parameter dict into a coupling
        specification. This is where you encode what is *fixed* (the coupling graph / pathway
        blocks, the field ``h``, ``base_loss``, ``transfer`` ...) and what is *fitted* (the keys
        of ``priors``). E.g. ``lambda p: pathway_blocks(sizes, within=p["within"], between=0,
        h=-0.5, base_loss=1.0, transfer=0.4)`` fits the within-pathway coupling ``within``.
    priors:
        ``{param: prior}`` over the ``spec_builder`` parameters (any names). A prior is a
        :class:`~zombi2.Distribution`, a bare float (fixed), a ``(low, high)`` uniform tuple,
        or a scipy dist / ``rng -> float`` callable. Values are **not** clamped to be positive
        (coupling parameters may be negative).
    statistics:
        Summary ``pm -> 1-D array``; defaults to :func:`cooccurrence_summary` (marginal + module
        structure), which is what makes coupling visible — the marginal alone can be blind to it.
    feature_weights, n_sims, accept, seed:
        As in :func:`match_profiles`.
    transfers, initial_presence:
        Passed through to :func:`~zombi2.simulate_coupled`.

    Returns
    -------
    ABCFit
        Accepted posterior over the fitted coupling parameters, with the usual diagnostics.
    """
    from .coupling import simulate_coupled           # local import avoids any import cycle

    profile_mat = (empirical if isinstance(empirical, ProfileMatrix)
                   else ProfileMatrix.from_tsv(empirical))
    species_order = list(profile_mat.species)
    summarize = statistics or cooccurrence_summary(species_order)
    uses_default = statistics is None
    target = summarize(profile_mat)

    priors = _normalize_priors(priors, restrict=False)     # arbitrary spec-builder param names
    names = list(priors.keys())
    weights = feature_weights if feature_weights is not None else getattr(summarize, "feature_weights", None)
    weights = None if weights is None else np.asarray(weights, float)
    rng = np.random.default_rng(seed)

    samples = np.empty((n_sims, len(names)))
    summaries = np.empty((n_sims, len(target)))
    for i in range(n_sims):
        vals = {n: priors[n].sample(rng) for n in names}   # no non-negativity clamp here
        sim_seed = int(rng.integers(1, 2**63 - 1))
        res = simulate_coupled(tree, spec_builder(vals), seed=sim_seed,
                               transfers=transfers, initial_presence=initial_presence)
        summaries[i] = summarize(res.profiles)
        samples[i] = [vals[n] for n in names]

    sd = summaries.std(axis=0)
    sd[sd == 0] = 1.0
    resid = (summaries - target) / sd
    if weights is not None:
        resid = resid * weights
    distances = np.sqrt((resid ** 2).sum(axis=1))

    order = np.argsort(distances, kind="stable")
    if isinstance(accept, float):
        if not 0.0 < accept <= 1.0:
            raise ValueError("float accept must be in (0, 1]")
        k = max(1, round(accept * n_sims))
    else:
        k = min(int(accept), n_sims)
    accepted = np.sort(order[:k])

    return ABCFit(
        param_names=names,
        samples=samples,
        distances=distances,
        accepted=accepted,
        tolerance=float(distances[order[k - 1]]),
        empirical_summary=target,
        priors=priors,
        accepted_summaries=summaries[accepted],
        summary_sd=sd,
        n_species=len(species_order),
        empirical=profile_mat,
        uses_default_summary=uses_default,
    )


# --- ABC-SMC (sequential Monte Carlo) --------------------------------------------

def _prepare(tree, empirical, statistics, gene_trees, model, family_shape,
             priors_spec, feature_weights):
    """Shared setup: resolve the summary/target, model builder, and priors."""
    if gene_trees:
        if not (hasattr(empirical, "profiles") and hasattr(empirical, "event_log")):
            raise TypeError("gene_trees=True requires an empirical Genomes (with gene trees)")
        profile_mat = empirical.profiles
        species_order = list(profile_mat.species)
        summarize = statistics or default_gene_tree_summary(species_order)
        target = summarize(empirical)
    else:
        profile_mat = (empirical if isinstance(empirical, ProfileMatrix)
                       else ProfileMatrix.from_tsv(empirical))
        species_order = list(profile_mat.species)
        summarize = statistics or default_summary(species_order)
        target = summarize(profile_mat)
    builder = _resolve_model(model, family_shape)
    priors = _normalize_priors(priors_spec, restrict=builder is None or isinstance(builder, _FamilyModel))
    names = list(priors.keys())
    weights = feature_weights if feature_weights is not None else getattr(summarize, "feature_weights", None)
    weights = None if weights is None else np.asarray(weights, float)
    return dict(summarize=summarize, target=target, weights=weights, profile_mat=profile_mat,
                species_order=species_order, builder=builder,
                return_genomes=gene_trees, uses_default=statistics is None,
                priors=priors, names=names)


def _uniform_bounds(priors, names):
    lo, hi = [], []
    for n in names:
        d = priors[n]
        if isinstance(d, Uniform):
            lo.append(d.low); hi.append(d.high)
        elif isinstance(d, Fixed):
            lo.append(d.value); hi.append(d.value)
        else:
            raise ValueError(f"match_profiles_smc needs Uniform (low, high) or fixed priors; "
                             f"{n!r} is {type(d).__name__}")
    return np.array(lo, float), np.array(hi, float)


def _smc_kernel(theta, parts, tau, active):
    """Gaussian perturbation-kernel density K(theta | each particle), up to a shared constant."""
    if not active.any():
        return np.ones(len(parts))
    d = (parts[:, active] - theta[active]) / tau[active]
    return np.exp(-0.5 * (d ** 2).sum(1) - np.log(tau[active]).sum())


def match_profiles_smc(
    tree: Tree,
    empirical,
    priors: dict,
    *,
    rounds: int = 5,
    n_particles: int = 200,
    quantile: float = 0.5,
    model=None,
    family_shape: float = 2.0,
    statistics=None,
    gene_trees: bool = False,
    feature_weights=None,
    initial_families: int = 20,
    max_family_size=None,
    transfers=None,
    seed: int | None = None,
    max_attempts_factor: int = 100,
) -> ABCFit:
    """Fit gene-family rates by **sequential Monte Carlo** ABC (Toni et al. 2009).

    ``initial_families`` is the number of gene families seeded at the root of each
    simulation (default 20).

    Plain rejection wastes most simulations far from the data. SMC instead evolves a
    population of ``n_particles`` over ``rounds`` with a shrinking tolerance: round 0 samples
    the prior; each later round resamples good particles, perturbs them with a Gaussian
    kernel (variance = twice the population variance), and keeps those under the new
    tolerance (the ``quantile`` of the previous round's distances). Importance weights keep
    the population an unbiased posterior sample. The result is a **sharper** posterior than
    rejection for a similar simulation budget.

    Takes the same model/summary/gene_trees options as :func:`match_profiles`, but
    **priors must be uniform** (``(low, high)`` / :class:`~zombi2.Uniform`) or fixed — the
    perturbation and weighting need bounded support. Returns an :class:`ABCFit` whose final
    population is weighted (``.summary()`` reports weighted credible intervals).
    """
    p = _prepare(tree, empirical, statistics, gene_trees, model, family_shape,
                 priors, feature_weights)
    summarize, target, w_feat = p["summarize"], p["target"], p["weights"]
    builder, return_genomes = p["builder"], p["return_genomes"]
    names, priors = p["names"], p["priors"]
    lo, hi = _uniform_bounds(priors, names)
    active = hi > lo
    P, N = len(names), n_particles
    rng = np.random.default_rng(seed)

    def sim_summ(theta):
        vals = {n: max(0.0, float(v)) for n, v in zip(names, theta)}
        out = _simulate(tree, vals, builder=builder, initial_size=initial_families,
                        max_family_size=max_family_size, transfers=transfers,
                        seed=int(rng.integers(1, 2**63 - 1)), return_genomes=return_genomes)
        return summarize(out)

    total_sims = N
    # round 0 — sample the prior; fix the distance scale from this batch
    part = np.array([rng.uniform(lo, hi) for _ in range(N)])
    summaries = np.array([sim_summ(th) for th in part])
    sd = summaries.std(0); sd[sd == 0] = 1.0

    def dist_of(s):
        r = (s - target) / sd
        return float(np.sqrt(((r * w_feat if w_feat is not None else r) ** 2).sum()))

    distances = np.array([dist_of(s) for s in summaries])
    W = np.full(N, 1.0 / N)

    for _ in range(1, rounds):
        eps = float(np.quantile(distances, quantile))
        mean = np.average(part, axis=0, weights=W)
        var = np.average((part - mean) ** 2, axis=0, weights=W)
        tau = np.sqrt(2 * var)
        tau[active] = np.maximum(tau[active], 1e-9)
        new_part = np.empty((N, P)); new_s = []; new_d = np.empty(N); new_W = np.empty(N)
        i = attempts = 0
        cap = max_attempts_factor * N
        while i < N and attempts < cap:
            attempts += 1
            theta = part[rng.choice(N, p=W)].copy()
            theta[active] += tau[active] * rng.standard_normal(int(active.sum()))
            if np.any(theta < lo) or np.any(theta > hi):
                continue
            s = sim_summ(theta); d = dist_of(s)
            if d < eps:
                new_part[i] = theta; new_s.append(s); new_d[i] = d
                new_W[i] = 1.0 / np.sum(W * _smc_kernel(theta, part, tau, active))  # uniform prior
                i += 1
        total_sims += attempts
        # Stop (keep the previous population) if the tolerance has hit the noise floor:
        # too few acceptances would just degenerate the population into a single point.
        if i < N or attempts >= cap:
            break
        part, summaries, distances = new_part, np.array(new_s), new_d
        W = new_W / new_W.sum()

    fit = ABCFit(
        param_names=names,
        samples=part,
        distances=distances,
        accepted=np.arange(len(part)),
        tolerance=float(distances.max()),
        empirical_summary=target,
        priors=priors,
        accepted_summaries=summaries,
        summary_sd=sd,
        n_species=len(p["species_order"]),
        empirical=p["profile_mat"],
        uses_default_summary=p["uses_default"],
        sample_weights=W,
    )
    fit.n_simulations = total_sims
    return fit
