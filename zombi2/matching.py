"""Fit gene-family rates to an empirical profile table by rejection ABC.

Given a fixed species tree and an **empirical copy-number profile** (families x extant
species), :func:`match_profiles` searches for rates under which ZOMBI2's gene-family model
reproduces that profile. It is a plain **Approximate Bayesian Computation** (ABC) rejection
sampler:

1. draw a parameter set from the priors,
2. simulate a profile matrix under it (the fast Rust path when available),
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
every family (Rust fast path). ``model="family"`` fits the same four values as the *means*
of per-family rate distributions (:class:`~zombi2.FamilySampledRates`), so families are
heterogeneous; this runs on the Python engine. You can also pass any callable
``params_dict -> RateModel`` as ``model`` for full generality.

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
    """Counts of families present in exactly ``k`` species, for ``k = 1 .. n_species``."""
    if pm.matrix.size == 0:
        return np.zeros(n_species)
    present = (pm.matrix > 0).sum(axis=1)
    return np.bincount(present, minlength=n_species + 1)[1:n_species + 1].astype(float)


def genome_sizes(pm: ProfileMatrix, species_order: list[str]) -> np.ndarray:
    """Total copy number per species, aligned to ``species_order`` (missing species -> 0)."""
    col_sum = pm.matrix.sum(axis=0) if pm.matrix.size else np.zeros(len(pm.species))
    by_name = {s: float(col_sum[j]) for j, s in enumerate(pm.species)}
    return np.array([by_name.get(s, 0.0) for s in species_order])


def copy_number_spectrum(pm: ProfileMatrix, max_copies: int = 4) -> np.ndarray:
    """Counts of present cells with copy number 1, 2, ..., ``>=max_copies``."""
    out = np.zeros(max_copies)
    if pm.matrix.size == 0:
        return out
    vals = pm.matrix[pm.matrix > 0]
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
    """Return ``(builder, fast_ok)``. ``builder`` is ``None`` for the uniform scalar model
    (which uses the Rust fast path); otherwise a callable ``params -> RateModel``."""
    if model is None or model == "uniform":
        return None, True
    if model == "family":
        return _FamilyModel(family_shape), False
    if callable(model):
        return model, False
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
        out: dict[str, dict[str, float]] = {}
        for name, vals in self.posterior.items():
            out[name] = {
                "mean": float(vals.mean()),
                "median": float(np.median(vals)),
                "lo95": float(np.quantile(vals, 0.025)),
                "hi95": float(np.quantile(vals, 0.975)),
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


def _resolve_engine(engine: str, fast_ok: bool) -> str:
    from .fast import rust_available
    if engine == "auto":
        return "fast" if (fast_ok and rust_available()) else "python"
    if engine not in ("fast", "python"):
        raise ValueError("engine must be 'auto', 'fast', or 'python'")
    if engine == "fast":
        if not fast_ok:
            raise ValueError("engine='fast' only supports model='uniform'; "
                             "use engine='python' for the family/custom models")
        if not rust_available():
            raise RuntimeError("engine='fast' requested but the zombi2_core Rust extension "
                               "is not built; use engine='python' or 'auto'")
    return engine


def _simulate(tree, vals, *, engine, builder, initial_size, max_family_size, transfers,
              seed, return_genomes=False):
    kw = dict(initial_size=initial_size, max_family_size=max_family_size,
              transfers=transfers, seed=seed)
    if builder is None and engine == "fast":   # uniform scalar model, Rust fast path
        from .fast import simulate_profiles_fast
        return simulate_profiles_fast(tree, **vals, **kw)
    from .simulation import simulate_genomes
    genomes = (simulate_genomes(tree, **vals, **kw) if builder is None
               else simulate_genomes(tree, builder(vals), **kw))
    return genomes if return_genomes else genomes.profiles


def _simulate_and_summarize(draw, tree, engine, builder, initial_size,
                            max_family_size, transfers, summarize, return_genomes):
    vals, seed = draw
    out = _simulate(tree, vals, engine=engine, builder=builder, initial_size=initial_size,
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
    initial_size: int = 20,
    max_family_size=None,
    transfers=None,
    engine: str = "auto",
    processes: int | None = None,
    seed: int | None = None,
) -> ABCFit:
    """Fit gene-family rates to an empirical profile by rejection ABC.

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
        (Rust fast path). ``"family"`` fits them as the *means* of per-family rate
        distributions (:class:`~zombi2.FamilySampledRates`, dispersion set by
        ``family_shape``; Python engine). Or pass a callable ``params -> RateModel``.
    statistics:
        A summary function ``pm -> 1-D array``; defaults to :func:`default_summary` over the
        empirical species. Must be picklable if ``processes`` is used. May expose a
        ``feature_weights`` attribute to weight its components in the distance.
    gene_trees:
        If True, use gene-tree information: ``empirical`` must be a :class:`~zombi2.Genomes`
        (not just a profile), the default summary becomes :func:`default_gene_tree_summary`
        (profile + weighted duplication/transfer/loss counts), and the Python engine is used
        (the Rust fast path yields no gene trees). Sharpens the gain-side rates; loss stays
        the hardest to identify.
    feature_weights:
        Optional per-summary-component weights applied to the (scaled) distance, so a small
        informative block of statistics is not drowned by a large one. Overrides any
        ``feature_weights`` attribute on the summary.
    n_sims:
        Number of prior draws to simulate.
    accept:
        A float in ``(0, 1]`` keeps that fraction of the closest draws; an int keeps that
        many.
    engine:
        ``"auto"`` (Rust fast path if available and the model allows it, else Python),
        ``"fast"``, or ``"python"``.
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
        if engine == "fast":
            raise ValueError("gene_trees=True requires the Python engine (the Rust fast path "
                             "produces no gene trees); use engine='auto' or 'python'")
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

    builder, fast_ok = _resolve_model(model, family_shape)
    if gene_trees:
        fast_ok = False                     # gene trees need the full (Python) engine
    priors = _normalize_priors(priors, restrict=builder is None or isinstance(builder, _FamilyModel))
    names = list(priors.keys())
    engine = _resolve_engine(engine, fast_ok)
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

    cfg = (tree, engine, builder, initial_size, max_family_size, transfers, summarize, gene_trees)
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
