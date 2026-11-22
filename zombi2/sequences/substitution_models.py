"""Substitution models — the **menu**.

A substitution model is the *chemistry* of a sequence: a ``K×K`` rate matrix ``Q`` (normalised to
one expected substitution per site per unit branch length) and its stationary frequencies ``π``.
Different models are genuinely different matrices — Jukes–Cantor, K80, HKY85, GTR differ in their
transition/transversion structure and base composition — so, unlike the clock, they do **not**
collapse to one grammar: they stay a menu of constructors, each taking its own physical parameters.

The four nucleotide models use four states, ``ACGT``: `jc69()`, `k80()`, `hky85()` and `gtr()`.
The five protein models use 20 states, `AMINO_ACIDS`: `poisson()`, `jtt()`, `dayhoff()`, `wag()`
and `lg()`. Researchers estimated the protein models' exchangeabilities and frequencies once from
large alignments. The protein models read those estimates from the published matrices in
``_aa_matrices``. They take no free parameters. `gy94()` implements Goldman and Yang's (1994) codon
model with codon states excluding stop codons. Runs write DNA with three bases per codon. Required
``omega`` sets selection on the protein.

**A matrix of your own** goes through `reversible()`: a symmetric exchangeability matrix ``S`` and
frequencies ``π``, normalised the same way, over any alphabet. It is the door the menu itself uses,
so `gtr()` is its four-state case and `lg()` its twenty-state one. It takes ``S`` and ``π`` rather
than ``Q`` because the engine requires a **time-reversible** model and ``S·π`` cannot be anything
else; a general ``Q`` is refused rather than evaluated wrongly (see `SubstitutionModel`).

**Across-site rate variation decorates any model on the menu** rather than adding entries to it:
`SubstitutionModel.across_sites()` returns the same chemistry with its sites sorted into
rate classes — a discretised Gamma (``+Γ``), a class that never changes (``+I``), or both. That is
where the field puts it (``HKY85+I+G4``) and, more to the point, it is not a **rate modifier**:
``SPEC §5``'s modifiers multiply *one* rate by a context factor, while this splits the sites of one
rate into classes. Putting it on the model keeps the rate grammar to one job and lets the spacer of
a nucleotide run carry its own (or no) variation.

Every model here is time-reversible, so the transition matrix over a branch of length ``t`` (in
substitutions/site), ``P(t) = exp(Q·t)``, is computed by eigendecomposition of the *symmetric*
matrix ``B = diag(√π)·Q·diag(1/√π)`` (numpy only, no scipy):
``P(t) = diag(1/√π)·V·exp(Λt)·Vᵀ·diag(√π)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from ._aa_matrices import (
    _DAYHOFF_EXCH, _DAYHOFF_PI, _JTT_EXCH, _JTT_PI, _LG_EXCH, _LG_PI, _WAG_EXCH, _WAG_PI,
)
from ._site_rates import discrete_gamma

#: the nucleotide alphabet, in the order ``Q`` and ``stationary`` follow.
BASES = "ACGT"

#: the 20-letter amino-acid alphabet, in the PAML column order every empirical protein matrix is
#: published in (``A R N D C Q E G H I L K M F P S T W Y V``) — the order ``Q`` and ``stationary``
#: follow for the protein models, and the order `decode()` reads them back in.
AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"

#: The genetic codes a codon model can be built on, by NCBI translation-table number: the amino
#: acid of each of the 64 codons, in the tables' own ``TCAG`` order (``TTT, TTC, TTA, TTG, TCT, …``),
#: with ``*`` for a stop. **11** is the bacterial, archaeal and plant-plastid code, whose amino
#: acids are the standard code's; **4** is the mycoplasma and spiroplasma code, where ``TGA`` is
#: tryptophan rather than a stop.
GENETIC_CODES: dict[int, str] = {
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
_TABLE_ORDER = tuple(a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG")


@dataclass(frozen=True)
class SubstitutionModel:
    """A ``K``-state reversible model: a normalised ``K×K`` rate matrix ``Q``, its stationary
    frequencies, and the ordered ``alphabet`` whose order ``Q`` / ``stationary`` follow.

    Built through the menu constructors (`jc69()`, `k80()`, `hky85()`, `gtr()`) or, for a matrix of
    your own, through `reversible()` — which takes a symmetric exchangeability matrix and the
    frequencies rather than ``Q`` itself, so a non-reversible model cannot be spelled. The reversible
    eigendecomposition behind `p_matrix()` is precomputed once in ``__post_init__``, and so is the
    check that the matrix it is being applied to is reversible at all: the transform is *only* a
    similarity transform of ``Q`` under detailed balance, so a matrix that violates it would give
    transition probabilities that are wrong without looking wrong. That check lives here rather than
    only in the constructors because this class is public, and every route into it has to end up in a
    state the engine can honour.

    ``site_rates`` and ``site_shares`` carry **across-site rate variation**: the sites of a sequence
    are sorted into classes, ``site_rates[c]`` multiplying the branch length for the sites in class
    ``c`` and ``site_shares[c]`` giving the proportion of sites that fall there. A plain model has
    the one class ``(1.0,)`` / ``(1.0,)`` and every site evolves at the same speed; `across_sites()`
    builds the rest. The two are dimensionless — a class scales a branch length, it is not a rate in
    ``SPEC §5``'s ``time⁻¹`` sense — and they obey one invariant:

        ``sum(rate × share) == 1``

    which is what makes a branch length the **mean** substitutions per site. Every phylogram in the
    result is drawn with that meaning, so the invariant is checked here rather than trusted.
    """

    name: str
    Q: np.ndarray
    stationary: np.ndarray
    alphabet: str = BASES
    #: the branch-length multiplier of each rate class, ascending; ``(1.0,)`` = no variation
    site_rates: tuple[float, ...] = (1.0,)
    #: the proportion of sites in each class, in the same order; sums to 1
    site_shares: tuple[float, ...] = (1.0,)
    #: a **profile mixture** (`profile_mixture`): one model per published profile, this model's
    #: exchangeabilities and rate classes over that profile's frequencies. Empty for a plain model.
    profile_models: tuple["SubstitutionModel", ...] = ()
    #: the weight of each profile, in the same order; sums to 1
    profile_weights: tuple[float, ...] = ()
    #: for a **codon** model (`gy94`), the triplet each state is written as, in the order ``Q`` and
    #: ``stationary`` follow. A state is then a codon and ``alphabet`` is the letters those triplets
    #: are spelled in, so the sequences a run writes are DNA, three letters to a state. Empty for a
    #: nucleotide or a protein model, where a state is one letter of ``alphabet``.
    codons: tuple[str, ...] = ()

    #: the reversible eigendecomposition `p_matrix()` runs on, filled in by ``__post_init__``
    #: (through ``object.__setattr__``, since the dataclass is frozen). Not part of the public
    #: surface: ``P(t) = _left · exp(_eigvals·t) · _right``.
    _eigvals: np.ndarray = field(init=False, repr=False, compare=False)
    _left: np.ndarray = field(init=False, repr=False, compare=False)
    _right: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self._check_reversible()
        if self.codons:
            ok = (len(self.codons) == self.Q.shape[0] and len(set(self.codons)) == len(self.codons)
                  and all(isinstance(c, str) and len(c) == 3 and set(c) <= set(self.alphabet)
                          for c in self.codons))
            if not ok:
                raise ValueError(
                    f"codons must name each of the model's {self.Q.shape[0]} states once, as "
                    f"three letters of its alphabet {self.alphabet!r}; build a codon model with "
                    f"gy94(), which does this.")
        rates, shares = self.site_rates, self.site_shares
        if len(rates) != len(shares) or not rates:
            raise ValueError(f"site_rates and site_shares must be the same non-empty length, got "
                             f"{len(rates)} and {len(shares)}")
        if any(r < 0 for r in rates) or not all(np.isfinite(rates)):
            raise ValueError(f"site rates must be finite and non-negative, got {rates}")
        if any(s <= 0 for s in shares) or not np.isclose(sum(shares), 1.0):
            raise ValueError(f"site shares must be positive and sum to 1, got {shares}")
        if not np.isclose(sum(r * s for r, s in zip(rates, shares)), 1.0):
            # the classes are what a branch length means; an unnormalised set would silently rescale
            # every phylogram in the run rather than raise anywhere near where it was built
            raise ValueError(
                f"site rates must average 1 over the sites (Σ rate×share), got "
                f"{sum(r * s for r, s in zip(rates, shares)):.6g} — a branch length is the mean "
                f"substitutions per site, so the classes have to be normalised to it. Build them "
                f"with across_sites(), which does this.")
        self._check_profiles()
        # Precompute the reversible eigendecomposition once (numpy only) for fast, scipy-free exp(Qt).
        pi = self.stationary
        sq = np.sqrt(pi)
        B = (sq[:, None] * self.Q) / sq[None, :]   # symmetric similarity transform of Q
        B = (B + B.T) / 2.0                         # kill round-off asymmetry before eigh
        w, V = np.linalg.eigh(B)
        # the pieces of P(t) = diag(1/√π) · V · exp(Λt) · Vᵀ · diag(√π)
        object.__setattr__(self, "_eigvals", w)
        object.__setattr__(self, "_left", V / sq[:, None])       # diag(1/√π) · V
        object.__setattr__(self, "_right", (V * sq[:, None]).T)   # Vᵀ · diag(√π)

    def _check_profiles(self) -> None:
        """A mixture's profile models must be this model over other frequencies: same alphabet, same
        rate classes, no mixture of their own. Built by `profile_mixture`, which guarantees it; the
        check is here because the fields are public and a hand-built mixture would otherwise evolve
        sites under a model the name does not claim."""
        models, weights = self.profile_models, self.profile_weights
        if len(models) != len(weights):
            raise ValueError(f"profile_models and profile_weights must be the same length, got "
                             f"{len(models)} and {len(weights)}")
        if not models:
            return
        if any(w <= 0 for w in weights) or not np.isclose(sum(weights), 1.0):
            raise ValueError(f"profile weights must be positive and sum to 1, got {weights}")
        for m in models:
            if (not isinstance(m, SubstitutionModel) or m.alphabet != self.alphabet
                    or m.site_rates != self.site_rates or m.site_shares != self.site_shares
                    or m.profile_models):
                raise ValueError(
                    f"every profile model of {self.name} must be a model over the same alphabet with "
                    f"the same rate classes and no mixture of its own; build the mixture with "
                    f"profile_mixture(), which does this.")

    def _check_reversible(self) -> None:
        """Refuse a matrix `p_matrix()` would silently mis-evaluate.

        `p_matrix()` does not exponentiate ``Q``. It eigendecomposes ``B = diag(√π)·Q·diag(1/√π)``,
        which is symmetric — and therefore cheap and stable to decompose with ``eigh`` — **only**
        when ``π_i·Q_ij = π_j·Q_ji`` for every pair, the detailed-balance condition that makes a
        model time-reversible. Given a ``Q`` that violates it, ``B`` is not symmetric, the
        symmetrisation ``(B + B.T)/2`` above quietly replaces it with a *different* matrix, and
        every ``P(t)`` that comes out is a well-formed stochastic matrix for the wrong model. The
        run would finish, the alignments would look like alignments, and nothing would say so. That
        is the reason this is checked rather than assumed, and checked here rather than only in the
        constructors: the class is public and in ``__all__``, so this is the one door every route
        goes through, `dataclasses.replace()` included.

        Unpickling restores the instance dict directly and does **not** re-run ``__post_init__``, so
        the parallel engine's workers pay nothing for this.
        """
        Q = np.asarray(self.Q, dtype=float)
        pi = np.asarray(self.stationary, dtype=float)
        k = Q.shape[0] if Q.ndim == 2 else -1
        if Q.ndim != 2 or Q.shape[1] != k or pi.shape != (k,):
            raise ValueError(
                f"{self.name} is not a K-state model: a rate matrix is K×K with K stationary "
                f"frequencies beside it, got Q of shape {Q.shape} and stationary of shape "
                f"{pi.shape}.")
        # A codon model writes each state as three letters, so its alphabet is the letters and its
        # `codons` are what has to match Q — checked in `__post_init__`, right after this.
        if self.codons:
            pass
        elif len(self.alphabet) != k:
            raise ValueError(
                f"{self.name}'s alphabet {self.alphabet!r} has {len(self.alphabet)} characters but "
                f"its Q is {k}×{k}. The alphabet is the *order* Q and stationary are written in and "
                f"the letters the sequences are decoded back into, so one character per state is "
                f"what makes those two agree. Use BASES (the 4 nucleotides), AMINO_ACIDS (the 20 "
                f"residues), or {k} distinct characters of your own.")
        if not self.codons and (len(set(self.alphabet)) != k or not self.alphabet.isascii()):
            raise ValueError(
                f"{self.name}'s alphabet {self.alphabet!r} must be {k} distinct ASCII "
                f"characters: a repeated letter would make two states decode the same and could "
                f"never be read back apart, and a non-ASCII one cannot go into a FASTA record this "
                f"level writes.")
        if pi.min() <= 0 or not np.isclose(pi.sum(), 1.0):
            raise ValueError(
                f"{self.name}'s stationary frequencies must be strictly positive and sum to 1, got "
                f"{pi} (summing to {pi.sum():.6g}). A zero-frequency state is not a rare state: the "
                f"reversible transform divides by √π, so it is a state the model cannot express at "
                f"all — drop it from the alphabet instead.")
        rows = Q.sum(axis=1)
        if not np.allclose(rows, 0.0):
            worst = int(np.argmax(np.abs(rows)))
            raise ValueError(
                f"{self.name}'s rate matrix has rows that do not sum to 0 — row {worst} sums to "
                f"{rows[worst]:.6g}. A rate matrix's diagonal is minus the rest of its row (a state "
                f"leaves exactly as fast as it goes anywhere else), so a non-zero row sum is not a "
                f"model with extra flux, it is a matrix that is not a generator. Build the model "
                f"with reversible(exchangeabilities, frequencies), which sets the diagonal itself.")
        flux = pi[:, None] * Q                 # π_i·Q_ij, the rate of i→j substitutions in a lineage
        gap = np.abs(flux - flux.T)
        if not np.allclose(flux, flux.T):
            i, j = np.unravel_index(int(np.argmax(gap)), gap.shape)
            raise ValueError(
                f"{self.name}'s rate matrix does not obey detailed balance: π_i·Q_ij must equal "
                f"π_j·Q_ji for every pair, and the worst mismatch here is {gap[i, j]:.6g} at states "
                f"{i} and {j} ({flux[i, j]:.6g} against {flux[j, i]:.6g}). That identity is what "
                f"makes a model time-reversible, and this engine needs it rather than merely liking "
                f"it: p_matrix computes exp(Qt) by eigendecomposing the symmetric matrix "
                f"diag(√π)·Q·diag(1/√π), which is similar to Q only under detailed balance. Given a "
                f"non-reversible Q it would return something that still looks like a transition "
                f"matrix and is not, so the run would produce plausible, wrong sequences instead of "
                f"failing. Build the model with reversible(exchangeabilities, frequencies): a symmetric "
                f"exchangeability matrix times π cannot express a non-reversible model. "
                f"Non-reversible models (UNREST and its relatives) are not implemented here — "
                f"that is a statement about this code, not about the models (SPEC §5).")

    @property
    def k(self) -> int:
        """Number of states in the alphabet."""
        return self.Q.shape[0]

    def p_matrix(self, t: float) -> np.ndarray:
        """Transition probabilities over branch length ``t`` (substitutions/site).

        ``P(t) = exp(Qt)`` via the reversible eigendecomposition; clipped to ``[0, ∞)`` to scrub tiny
        negative round-off so every row is a valid probability distribution.
        """
        # The BLAS matmul kernel can raise spurious FP flags on larger matrices even when every input
        # is finite and the result is a valid stochastic matrix; silence them — the clip is the guard.
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            P = (self._left * np.exp(self._eigvals * t)) @ self._right
        return np.clip(P, 0.0, None)

    def across_sites(self, *, gamma_shape: float | None = None, invariant: float = 0.0,
                     rate_categories: int = 4) -> "SubstitutionModel":
        """The same model with its sites sorted into **rate classes** — ``+Γ``, ``+I``, or both.

        Every site of a gene evolving at exactly one speed is a model no real gene obeys: some
        positions are held nearly fixed by function, others race ahead. Two standard knobs say so,
        and they compose:

        - ``gamma_shape`` (α) draws each site's rate from a **Gamma** with mean 1, discretised into
          ``rate_categories`` equal-probability classes (Yang 1994; see `_site_rates`). A small shape
          is strong variation — 0.5 gives a few fast sites among many slow ones — and a large one is
          nearly flat. There is no default: asking for ``+Γ`` means choosing how unequal the sites are.
        - ``invariant`` (0 by default) is the proportion of sites that **never** change: one more
          class, at rate 0. Real alignments have columns that are constant because the site cannot
          change, not because it happened not to, and a Gamma alone fits those badly.

        A site's class is drawn once and holds for the whole tree — that is what across-site
        variation means, as against a rate that varies along a branch. Sites stay independent, so
        this is a knob on the existing engine and not a new one (``SPEC §9``).

        **The classes are normalised to mean 1**, invariant sites included: the Gamma classes are
        scaled by ``1 / (1 - invariant)`` to make room for the sites that never change. So a branch
        length stays the **mean** substitutions per site, a run with variation and one without have
        the same phylogram at the same rate, and only the spread of change across columns differs.

        This is **not a rate modifier** (``SPEC §5``), and `simulate_sequences` will still refuse one
        on ``substitution``: a modifier multiplies one rate by a context factor, while this splits
        one rate's sites into classes. It belongs to the model, which is where the field puts it —
        the decorated ``name`` reads ``HKY85+I+G4``, and it is what the run reports.

        Returns a new model; the original is unchanged (they are frozen).
        """
        if len(self.site_rates) > 1:
            # Composing two calls would have to decide what a Gamma over an already-classed model
            # means, and there is exactly one sensible spelling of every combination this offers —
            # both knobs in one call. Refusing is the cheap half of that; the expensive half would be
            # a second set of classes replacing the first while the name claimed both.
            raise ValueError(
                f"{self.name} already varies across sites, and across_sites() replaces the classes "
                f"rather than layering on them. Give every knob in one call — "
                f"across_sites(gamma_shape=…, invariant=…) — starting from the plain model.")
        if isinstance(invariant, bool) or not isinstance(invariant, (int, float)):
            raise TypeError(f"invariant must be a real proportion of sites, got {invariant!r}")
        if not 0.0 <= invariant < 1.0:
            raise ValueError(
                f"invariant is the proportion of sites that never change, so it must be in [0, 1), "
                f"got {invariant!r}" + (" — at 1.0 no site could ever change and there would be no "
                                        "sequence evolution left to simulate" if invariant == 1.0 else ""))
        # named before the nothing-to-vary case below: someone who typed a category count meant to
        # ask for a Gamma, and the useful reply names the argument they left out
        if gamma_shape is None and rate_categories != 4:
            raise ValueError(
                f"rate_categories={rate_categories} counts the classes of the Gamma, but no "
                f"gamma_shape was given — add gamma_shape=…, or drop rate_categories (invariant= "
                f"on its own is a single never-changing class and takes no count).")
        if gamma_shape is None and not invariant:
            raise ValueError(
                "across_sites() has nothing to vary — give gamma_shape=… for a Gamma of rates "
                "across sites, invariant=… for a class of sites that never change, or both.")

        rates: list[float] = []
        shares: list[float] = []
        suffix = ""
        if invariant:
            rates.append(0.0)
            shares.append(float(invariant))
            suffix += "+I"
        # Everything not in the invariant class carries the whole mean of 1 between it, because the
        # invariant sites contribute nothing — hence the 1/(1 - invariant) scaling. Without it, +I
        # would quietly slow the whole sequence down instead of concentrating its change.
        varying = 1.0 - invariant
        scale = 1.0 / varying
        if gamma_shape is not None:
            rates.extend(r * scale for r in discrete_gamma(gamma_shape, rate_categories))
            shares.extend([varying / rate_categories] * rate_categories)
            suffix += f"+G{rate_categories}"
        else:
            # +I on its own: one class for every site that can change at all
            rates.append(scale)
            shares.append(varying)
        # replace() re-runs __post_init__, so the eigendecomposition is rebuilt (once) and the
        # mean-1 invariant above is checked on the way out rather than assumed. A mixture's profile
        # models take the same classes: a profile says which residues, the Gamma how fast.
        rates_t, shares_t = tuple(rates), tuple(shares)
        return replace(self, name=self.name + suffix, site_rates=rates_t, site_shares=shares_t,
                       profile_models=tuple(replace(m, site_rates=rates_t, site_shares=shares_t)
                                            for m in self.profile_models))

    def profile_mixture(self, mixture: str) -> "SubstitutionModel":
        """The same exchangeabilities over a published **mixture of site profiles**: ``C10`` to
        ``C60`` (Quang, Gascuel & Lartillot 2008), the site-heterogeneous models the field writes
        ``LG+C60``.

        A profile is one set of amino-acid frequencies. A mixture is ``K`` of them, each with a
        weight. When a gene evolves under the mixture, every site draws one profile by weight, once,
        and evolves under this model's exchangeabilities over that profile's frequencies. So a
        buried site can be hydrophobic and a surface site charged without a table of your own: that
        table (`simulate_sequences(profiles=…)`) is the other way to say the same thing, site by
        site, and the two are refused together.

        Composes with `across_sites` in either order: ``lg().profile_mixture("C60")
        .across_sites(gamma_shape=0.5)`` is ``LG+C60+G4``, a profile saying *which* residues and
        the Gamma *how fast*. Protein models only: the profiles are amino-acid frequencies. The
        profiles are ``profile_models`` and their weights ``profile_weights``; the run records
        which profile each site drew (``profile_assignments``).

        Returns a new model; the original is unchanged.
        """
        from ._profile_mixtures import MIXTURES
        if self.profile_models:
            raise ValueError(
                f"{self.name} already carries a profile mixture, and a second one would have to "
                f"say which of the two a site draws from. Start again from the plain model.")
        if self.alphabet != AMINO_ACIDS:
            raise ValueError(
                f"a profile mixture is a set of amino-acid profiles, and {self.name} is a "
                f"{self.k}-state model written over {self.alphabet!r}. Start from a protein "
                f"model: lg(), wag(), jtt(), dayhoff() or poisson().")
        label = mixture.upper() if isinstance(mixture, str) else mixture
        if label not in MIXTURES:
            raise ValueError(
                f"unknown profile mixture {mixture!r}; the published ones are "
                f"{', '.join(MIXTURES)}.")
        weights, rows = MIXTURES[label]
        w = np.asarray(weights, dtype=float)
        w = w / w.sum()
        models = tuple(_with_frequencies(self, row, name=f"{self.name}+{label}[{i}]")
                       for i, row in enumerate(rows, start=1))
        return replace(self, name=f"{self.name}+{label}", profile_models=models,
                       profile_weights=tuple(float(x) for x in w))


#: the published profile mixtures `SubstitutionModel.profile_mixture` knows, by label
PROFILE_MIXTURES: tuple[str, ...] = ("C10", "C20", "C30", "C40", "C50", "C60")


def _with_frequencies(model: SubstitutionModel, frequencies, *, name: str) -> SubstitutionModel:
    """The same model — same exchangeabilities, same rate classes — over different equilibrium
    frequencies.

    What a **site profile** needs (`simulate_sequences(profiles=…)`): the chemistry of which amino
    acids swap easily is a property of the substitution process and is kept, while *which* residues
    belong at this position is the profile's to say. ``S`` is recovered from the model as
    ``Q_ij / pi_j``, symmetric by construction because that is how ``Q`` was built; the halving
    kills round-off, not asymmetry."""
    S = model.Q / model.stationary[None, :]
    S = (S + S.T) / 2.0
    np.fill_diagonal(S, 0.0)
    out = _reversible_model(name, S, frequencies, model.alphabet, model.codons)
    # the rate classes ride along: a profile says which residues, `+G` says how fast, and the two
    # are independent (see `across_sites`)
    return replace(out, site_rates=model.site_rates, site_shares=model.site_shares)


def _reversible_model(name: str, S: np.ndarray, pi, alphabet: str = BASES,
                      codons: tuple[str, ...] = ()) -> SubstitutionModel:
    """Build a normalised reversible model from a symmetric exchangeability matrix ``S`` and ``pi``.

    ``Q_ij = S_ij · pi_j`` (i≠j), scaled so the expected rate ``-Σ pi_i Q_ii = 1`` (branch lengths
    are then in substitutions/site). ``S`` must be symmetric with a zero diagonal.
    """
    S = np.asarray(S, dtype=float)
    pi = np.asarray(pi, dtype=float)
    k = pi.shape[0]
    if pi.shape != (k,) or pi.min() <= 0 or not np.isclose(pi.sum(), 1.0):
        raise ValueError(f"stationary frequencies must be strictly positive and sum to 1, got {pi} "
                         "(a zero-frequency state makes the rate matrix degenerate)")
    if S.shape != (k, k) or (S < 0).any() or not np.allclose(S, S.T):
        raise ValueError(
            f"exchangeabilities must be a symmetric non-negative K×K matrix, and with K = {k} "
            f"(the number of frequencies) this one is {S.shape} with minimum entry "
            f"{S.min() if S.size else float('nan'):.6g}"
            + ("" if S.shape != (k, k) or np.allclose(S, S.T)
               else f" and asymmetric by up to {np.abs(S - S.T).max():.6g}")
            + ". Symmetry is not tidiness: S_ij = S_ji is exactly what makes Q_ij = S_ij·π_j "
              "reversible, which is what the engine's exp(Qt) requires.")
    if S.shape == (k, k) and np.any(np.diag(S) != 0.0):
        # Q's diagonal is set from its row sums below, so a diagonal given here would be dropped on
        # the floor. Say so rather than accept a number and ignore it: someone who wrote one meant
        # something by it, and what they meant is not expressible.
        raise ValueError(
            f"the exchangeability matrix must have a zero diagonal, got S_ii = {np.diag(S)}. S_ii "
            "is not a rate a state leaves itself at — the diagonal of Q is set from the row sums, "
            "so anything written here would be silently discarded.")
    pi = pi / pi.sum()          # renormalise (published frequencies round to 1 only to ~1e-6)
    Q = S * pi[None, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    scale = -(pi * np.diag(Q)).sum()
    if scale <= 0:
        raise ValueError("degenerate substitution model (zero rate)")
    # One substitution per STATE is the scale, and a codon model's state is three bases. Every
    # change it makes is one base, so three per codon is one per base — the unit every other model
    # on the menu writes a branch length in.
    per_state = 3.0 if codons else 1.0
    return SubstitutionModel(name, Q * (per_state / scale), pi, alphabet, codons=codons)


def _gtr_model(name: str, exch, pi) -> SubstitutionModel:
    """Build a GTR-family nucleotide model from 6 exchangeabilities ``[AC,AG,AT,CG,CT,GT]`` and frequencies.

    ``Q_ij = exch_ij · pi_j`` (i≠j), scaled so ``-Σ pi_i Q_ii = 1`` (branch lengths in subs/site).
    """
    pi = np.asarray(pi, dtype=float)
    if pi.shape != (4,):
        raise ValueError(f"stationary frequencies must be 4 values, got {pi}")
    ac, ag, at, cg, ct, gt = (float(x) for x in exch)
    if min(ac, ag, at, cg, ct, gt) < 0:
        raise ValueError("exchangeability rates must be non-negative")
    S = np.array([[0, ac, ag, at], [ac, 0, cg, ct], [ag, cg, 0, gt], [at, ct, gt, 0]], dtype=float)
    return _reversible_model(name, S, pi, BASES)


def jc69() -> SubstitutionModel:
    """Jukes–Cantor (1969): equal rates, equal base frequencies. No free parameters."""
    return _gtr_model("JC69", [1, 1, 1, 1, 1, 1], [0.25] * 4)


def k80(kappa: float = 2.0) -> SubstitutionModel:
    """Kimura 2-parameter (1980): transition/transversion ratio ``kappa``, equal frequencies."""
    return _gtr_model("K80", [1, kappa, 1, 1, kappa, 1], [0.25] * 4)


def hky85(kappa: float = 2.0, frequencies=(0.25, 0.25, 0.25, 0.25)) -> SubstitutionModel:
    """HKY85 (Hasegawa–Kishino–Yano 1985): transition bias ``kappa`` with unequal base ``frequencies`` (A,C,G,T)."""
    return _gtr_model("HKY85", [1, kappa, 1, 1, kappa, 1], frequencies)


def gtr(exchangeabilities=(1, 1, 1, 1, 1, 1),
        frequencies=(0.25, 0.25, 0.25, 0.25)) -> SubstitutionModel:
    """General time-reversible: 6 ``exchangeabilities`` ``[AC,AG,AT,CG,CT,GT]`` and ``frequencies``
    (A,C,G,T). This is `reversible()` restricted to four states, and the two name the numbers the
    same thing on purpose — they are the same numbers."""
    return _gtr_model("GTR", exchangeabilities, frequencies)


def reversible(exchangeabilities, frequencies, *, name: str = "Custom",
               alphabet: str = BASES) -> SubstitutionModel:
    """A model of **your own matrix**: a symmetric ``K×K`` exchangeability matrix and ``K``
    stationary frequencies, over any alphabet.

    ``Q_ij = S_ij · π_j`` for ``i ≠ j``, the diagonal set from the row sums, and the whole matrix
    scaled so that ``-Σ π_i Q_ii = 1`` — one expected substitution per site per unit branch length,
    exactly as every model on the menu is scaled. That normalisation is what keeps a branch length
    meaning **substitutions per site** whatever matrix produced it, so a phylogram from a model built
    here is comparable with one from `hky85()` without a conversion.

    This is the same door the menu already goes through, so the menu is its special cases: `gtr()` is
    this with ``K = 4`` and the six nucleotide exchangeabilities spelled out, and `lg()` is this with
    ``K = 20`` and a published triangle. Give ``alphabet=AMINO_ACIDS`` for a 20-state protein matrix
    of your own, or any ``K`` distinct ASCII characters for an alphabet that is neither — though the
    rest of the level only makes sense for those two: a **nucleotide** genome run reads its blocks on
    both strands and so refuses anything but ``ACGT``, and the FASTA it writes is whatever letters
    you chose.

    **You give S and π rather than Q, and that is the point.** The engine's `p_matrix()` computes
    ``exp(Qt)`` by eigendecomposing ``diag(√π)·Q·diag(1/√π)``, which is a symmetric matrix similar to
    ``Q`` only when ``π_i Q_ij = π_j Q_ji`` — detailed balance. A symmetric ``S`` times ``π`` obeys
    that by construction, so there is no way to spell a model here that the engine would evaluate
    wrongly. A general ``Q`` handed straight to `SubstitutionModel` could be, which is why the
    dataclass checks the same identity and refuses. Non-reversible models (UNREST and its relatives)
    are **not implemented** — a statement about this code, not about the models (SPEC §5).

    ``S`` must be symmetric, non-negative, and zero on the diagonal; ``π`` strictly positive and
    summing to 1. ``name`` is what the run prints and what error messages quote, so it is worth
    setting to whatever the matrix actually is.

    There is deliberately **no command-line flag** for this. SPEC §5 asks for one written form
    everywhere, and a ``K×K`` matrix has no written form a command line can carry: at ``K = 20`` it
    is 190 numbers, which is a file format, not an argument. It stays a Python constructor; the CLI's
    ``--model`` offers the menu, which is what a flag can honestly hold.

    Experimental (SPEC §9): Python API only, and the spelling may still move.
    """
    return _reversible_model(name, exchangeabilities, frequencies, alphabet)


# --- the codon model: 61 or 62 states, mechanistic ------------------------------------------------

def _codon_table(code: int) -> tuple[tuple[str, ...], dict[str, str]]:
    """``(sense codons in ACGT order, {codon: amino acid})`` for one genetic code."""
    if isinstance(code, bool) or code not in GENETIC_CODES:
        raise ValueError(
            f"code is an NCBI translation table, and the ones built in are "
            f"{', '.join(map(str, sorted(GENETIC_CODES)))} — 11 for bacteria and archaea (the "
            f"standard amino acids), 4 where TGA is tryptophan (mycoplasmas). Got {code!r}.")
    amino = dict(zip(_TABLE_ORDER, GENETIC_CODES[code]))
    return tuple(sorted(c for c in _TABLE_ORDER if amino[c] != "*")), amino


def gy94(*, omega: float, kappa: float = 2.0, frequencies=None, code: int = 11) -> SubstitutionModel:
    """``gy94`` implements the codon model of Goldman and Yang (1994). Each codon state contributes
    three bases to the DNA sequences that runs write. The codon model excludes stop codons, so runs
    never write them. Changes to two or three bases at once have substitution rate 0.

    Required ``omega`` multiplies substitution rates for amino acid changes relative to changes
    preserving the amino acid: 1 means no selection on the protein, 0.1 permits one such change in
    ten. ``kappa`` multiplies transition substitution rates (A with G, C with T) relative to
    transversions, with default 2.0, as for ``hky85``. ``frequencies`` accepts four base
    frequencies in A, C, G, T order or one frequency per codon excluding stop codons in
    alphabetical order. Four base frequencies give each codon's frequency as the product of its
    three base frequencies, renormalised over codons excluding stop codons. The default gives equal
    codon frequencies.

    ``code`` selects the NCBI genetic code, defaulting to 11 for bacteria and archaea, with 61
    codons excluding stop codons. Code 4 supports mycoplasmas, assigns TGA to tryptophan, and has
    62 codons excluding stop codons. A branch length measures substitutions per base, preserving
    the ``hky85`` meanings of ``divergence``, phylograms and the saturation warning. PAML and
    IQ-TREE report substitutions per codon, giving branch lengths three times larger."""
    for label, value in (("omega", omega), ("kappa", kappa)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) \
                or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be a positive number, got {value!r}")
    codons, amino = _codon_table(code)
    k = len(codons)
    if frequencies is None:
        pi = np.full(k, 1.0 / k)
    else:
        given = np.asarray(frequencies, dtype=float)
        if given.shape == (4,):
            # base frequencies, in ACGT order: a codon's is the product of its three bases',
            # renormalised over the codons that are not stops
            if given.min() <= 0 or not np.isclose(given.sum(), 1.0):
                raise ValueError(f"the four base frequencies must be positive and sum to 1, got "
                                 f"{tuple(given)}")
            base = dict(zip(BASES, given))
            pi = np.array([base[c[0]] * base[c[1]] * base[c[2]] for c in codons])
            pi = pi / pi.sum()
        elif given.shape == (k,):
            pi = given
        else:
            raise ValueError(
                f"frequencies takes 4 base frequencies (A, C, G, T) or {k} codon frequencies, one "
                f"per codon of code {code} that is not a stop, in alphabetical order (AAA, AAC, "
                f"AAG, …). Got {given.size} numbers.")
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    S = np.zeros((k, k))
    for i, x in enumerate(codons):
        for j in range(i + 1, k):
            y = codons[j]
            changed = [(a, b) for a, b in zip(x, y) if a != b]
            if len(changed) != 1:               # two or three bases apart: not one event
                continue
            rate = (kappa if changed[0] in transitions else 1.0) \
                * (1.0 if amino[x] == amino[y] else omega)
            S[i, j] = S[j, i] = rate
    return _reversible_model("GY94", S, pi, BASES, codons)


# --- the protein models: 20 states, empirical exchangeabilities + frequencies ----------------------

def _lower_triangle(tri, k: int) -> np.ndarray:
    """Expand a flat lower triangle (entry ``(i, j)`` for ``i = 1..k-1``, ``j < i``, row by row — the
    PAML layout of `_aa_matrices`) into the symmetric ``k×k`` exchangeability matrix."""
    S = np.zeros((k, k))
    it = iter(tri)
    for i in range(1, k):
        for j in range(i):
            S[i, j] = S[j, i] = next(it)
    return S


def _empirical_protein(name: str, tri, pi) -> SubstitutionModel:
    """Build a 20-state protein model from published lower-triangular exchangeabilities and frequencies,
    both in `AMINO_ACIDS` order — normalised, like every model here, to one expected
    substitution per site per unit branch length."""
    return _reversible_model(name, _lower_triangle(tri, 20), pi, AMINO_ACIDS)


def poisson() -> SubstitutionModel:
    """Poisson: equal exchangeabilities, equal frequencies — the JC69 of proteins. No free parameters."""
    S = np.ones((20, 20)) - np.eye(20)
    return _reversible_model("Poisson", S, np.full(20, 1.0 / 20.0), AMINO_ACIDS)


def jtt() -> SubstitutionModel:
    """JTT (Jones, Taylor & Thornton 1992): the empirical matrix from close protein homologues."""
    return _empirical_protein("JTT", _JTT_EXCH, _JTT_PI)


def dayhoff() -> SubstitutionModel:
    """Dayhoff (Dayhoff, Schwartz & Orcutt 1978): the original PAM matrix, in PAML's values."""
    return _empirical_protein("Dayhoff", _DAYHOFF_EXCH, _DAYHOFF_PI)


def wag() -> SubstitutionModel:
    """WAG (Whelan & Goldman 2001): estimated by maximum likelihood over a wide protein database."""
    return _empirical_protein("WAG", _WAG_EXCH, _WAG_PI)


def lg() -> SubstitutionModel:
    """LG (Le & Gascuel 2008): WAG's successor, fitted with across-site rate variation — the
    default protein model of modern phylogenetics."""
    return _empirical_protein("LG", _LG_EXCH, _LG_PI)


#: ASCII lookup tables for `decode()`, one per alphabet, built once on first use. The table is
#: tiny and read-only, so it is safe to reuse across calls — and `decode()` runs once per node of
#: every gene tree, so rebuilding it each time was pure waste.
_DECODE_LUT: dict[str, np.ndarray] = {}


def decode(states: np.ndarray, alphabet: str = BASES) -> str:
    """Map an array of integer states back to a string over ``alphabet`` — ``ACGT`` by default, or
    `AMINO_ACIDS` for a protein model (callers pass ``model.alphabet``).

    ``states`` are indices into ``alphabet``, so the whole array is one numpy gather into an ASCII
    lookup table — ``lut[states]`` — read out in a single ``.tobytes().decode()`` rather than one
    Python step per site. This is called once per node of every gene tree, so the per-site loop it
    replaces was the dominant cost of a sequence run; the result is byte-for-byte the same string. The
    lookup table itself is cached per alphabet (`_DECODE_LUT`), built once instead of per call.

    ``states`` may be multi-dimensional: a 2-D ``(rows, length)`` array decodes to the rows' strings
    concatenated back to back (row-major), which lets a caller decode a whole gene tree's nodes in one
    gather + one ASCII decode and slice the fixed-length rows out — see `_split()`."""
    lut = _DECODE_LUT.get(alphabet)
    if lut is None:
        lut = np.frombuffer(alphabet.encode("ascii"), dtype=np.uint8)
        _DECODE_LUT[alphabet] = lut
    return lut[np.asarray(states)].tobytes().decode("ascii")


def encode(seq: str, alphabet: str = BASES) -> np.ndarray:
    """The inverse of `decode()`: a string over ``alphabet`` to its integer states. Used to found a
    run's blocks from a real ``fasta=`` — the supplied DNA becomes a block's founding states. A character not in
    ``alphabet`` raises (the FASTA reader already rejects non-``ACGT``, so this is a second guard)."""
    index = {c: i for i, c in enumerate(alphabet)}
    try:
        return np.fromiter((index[c] for c in seq), dtype=np.int8, count=len(seq))
    except KeyError as e:
        raise ValueError(f"sequence has {e.args[0]!r}, not in the model's alphabet {alphabet!r}") from None


#: one lookup table per set of codons: ``lut[state]`` is the three ASCII bytes of that state
_CODON_LUT: dict[tuple[str, ...], np.ndarray] = {}


def letters_of(model: SubstitutionModel, states: np.ndarray) -> str:
    """`decode()` for a **model**: the letters its states are written as — one letter a state for a
    nucleotide or a protein model, three for a codon model. ``states`` may be multi-dimensional, as
    for `decode()`; the rows come back end to end, each ``width × length`` letters."""
    if not model.codons:
        return decode(states, model.alphabet)
    lut = _CODON_LUT.get(model.codons)
    if lut is None:
        lut = np.frombuffer("".join(model.codons).encode("ascii"), dtype=np.uint8).reshape(-1, 3)
        _CODON_LUT[model.codons] = lut
    return lut[np.asarray(states)].tobytes().decode("ascii")


def states_of(model: SubstitutionModel, seq: str) -> np.ndarray:
    """`encode()` for a **model**: a sequence in the model's letters to its integer states. For a
    codon model the sequence is read three letters at a time, so it has to be a whole number of
    codons, and none of them a stop — a stop is not a state the model has."""
    if not model.codons:
        return encode(seq, model.alphabet)
    if len(seq) % 3:
        raise ValueError(f"the sequence has {len(seq)} letters, which is not a whole number of "
                         f"codons")
    index = {c: i for i, c in enumerate(model.codons)}
    try:
        return np.fromiter((index[seq[i:i + 3]] for i in range(0, len(seq), 3)),
                           dtype=np.int8, count=len(seq) // 3)
    except KeyError as e:
        codon = e.args[0]
        if set(codon) <= set(model.alphabet):
            raise ValueError(
                f"codon {codon!r} is a stop codon in this model's genetic code, so it is not a "
                f"state the model has. Leave a gene's stop codon out of its sequence") from None
        raise ValueError(f"sequence has {codon!r}, which is not three letters of the model's "
                         f"alphabet {model.alphabet!r}") from None


def width_of(model: SubstitutionModel) -> int:
    """How many letters one state is written as: 3 for a codon model, 1 otherwise."""
    return 3 if model.codons else 1


__all__ = ["SubstitutionModel", "jc69", "k80", "hky85", "gtr", "gy94", "reversible",
           "poisson", "jtt", "dayhoff", "wag", "lg", "decode", "encode", "BASES", "AMINO_ACIDS",
           "GENETIC_CODES", "PROFILE_MIXTURES"]
