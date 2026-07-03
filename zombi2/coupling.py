"""Non-independent gene families — a Potts/Ising-style coupling rate model.

Today every other rate model evolves each family independently, so the phylogenetic
profile (:class:`~zombi2.ProfileMatrix`) correlates rows *only* through the shared species
tree. Real genomes correlate through **function**: families in the same pathway or complex
are present or absent together (Pellegrini 1999), a signal that (inverse) Potts / DCA
methods exploit (Croce 2019; Fukunaga & Iwasaki 2022). This module is the *forward*,
generative counterpart: it injects prescribed pairwise couplings ``J`` (and fields ``h``)
so that simulated profiles carry a **known ground-truth** coupling structure to benchmark
those inference methods against.

The design note ``docs/non_independence.tex`` derives three variants; this is the
implemented one, ``docs/coupling_model.md`` documents the choices. In one paragraph:

**The model.** Presence/absence of a fixed panel of ``N`` families inside one genome is an
Ising vector ``σ ∈ {0,1}^N``. Fields ``h_i`` and couplings ``J_ij`` define a local field

    f_i(σ) = h_i + Σ_j J_ij · σ_j          (partners only; J has a zero diagonal)

**Coupling enters through loss.** A *present* family is lost at rate

    loss_i = base_loss · exp(-β · f_i),

so present partners with ``J_ij > 0`` protect ``i`` (lower loss → co-occurrence) and
``J_ij < 0`` promote its loss (avoidance). **Gain is horizontal transfer**: the stock,
field-blind ``TRANSFER`` event carries a family into a recipient, and the coupled loss then
*selectively retains* it — kept where its partners are present, purged where they are not.
That differential retention of horizontally-acquired genes is what writes ``J`` into the
profiles.

**Faithfulness.** Because the gain channel is field-blind (HGT, not a detailed-balance
Glauber gain), this is an *approximate* Potts generator: recovered couplings track the
injected ``J`` in **sign and rank**, not as an exact Boltzmann constant. That is the
deliberate price of keeping regain mechanistic — a lost family returns only from a donor
that still has it, exactly as in real genomes.

**Architecture.** :class:`PottsRates` is an ordinary :class:`~zombi2.RateModel`: the
simulator already hands it the whole genome via ``event_weights(genome, branch, time)``, so
coupling needs **no change to the simulator, sampler, genome or output**. A custom rate
model is automatically ineligible for the Rust fast path, so a coupled run uses the
pure-Python engine (the coupling breaks per-family independence, so the fast path could not
apply anyway). Cost is ``O(N + nnz(J))`` per event — fine at benchmark scale.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

from .events import EventType
from .genome import Gene, UnorderedGenome
from .genome_sim import GenomeSimulator
from .profiles import ProfileMatrix
from .rates import EventWeight, RateModel
from .transfers import TransferModel

#: Clamp on the loss exponent ``-β·f_i`` so extreme fields never overflow ``exp`` (or drive
#: the Gillespie loop into a gain/instant-loss hot cycle). ``exp(±40)`` already spans ~35
#: orders of magnitude, far beyond any useful coupling.
_MAX_EXPONENT = 40.0


def _natkey(name: str) -> tuple[int, str]:
    """Natural-ish sort key (matches :mod:`zombi2.profiles`)."""
    digits = re.sub(r"\D", "", name)
    return (int(digits) if digits else 0, name)


# ═══════════════════════════════════════════════════════════════════════════════
# Coupling specification
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class CouplingSpec:
    """A fixed panel of ``n_families`` families with pairwise couplings ``J`` and fields ``h``.

    Build it from a dense matrix (:meth:`from_dense`), a sparse edge list
    (:meth:`from_edges`), or pathway blocks (:func:`pathway_blocks`). The panel's families
    are named ``F0 .. F{n-1}`` by default (override with ``prefix``); ``J[i, j]`` couples
    family ``i`` and ``j``.

    Parameters
    ----------
    n_families : the panel size ``N``.
    adjacency  : per-family neighbour list ``adjacency[i] = [(j, J_ij), ...]`` (symmetric,
                 zero diagonal). Prefer the constructors below to building this by hand.
    h          : length-``N`` field vector (intrinsic retention bias); ``None`` → zeros.
    beta       : global coupling strength (inverse temperature); scales the whole field.
    base_loss  : baseline per-family loss rate (the loss at ``f_i = 0``).
    transfer   : per-copy HGT rate — the (field-blind) gain channel.
    origination: background rate of brand-new, *uncoupled* families (0 → closed panel).
    prefix     : family-id prefix; family ``i`` is ``f"{prefix}{i}"``.
    """

    n_families: int
    adjacency: list[list[tuple[int, float]]]
    h: np.ndarray
    beta: float = 1.0
    base_loss: float = 1.0
    transfer: float = 0.5
    origination: float = 0.0
    prefix: str = "F"

    def __post_init__(self) -> None:
        if self.n_families <= 0:
            raise ValueError("n_families must be positive")
        self.h = np.asarray(self.h, dtype=float)
        if self.h.shape != (self.n_families,):
            raise ValueError(f"h must have shape ({self.n_families},), got {self.h.shape}")
        for name in ("beta", "base_loss", "transfer", "origination"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be >= 0")
        self.panel_ids: list[str] = [f"{self.prefix}{i}" for i in range(self.n_families)]
        self.index: dict[str, int] = {fam: i for i, fam in enumerate(self.panel_ids)}

    # --- constructors ------------------------------------------------------
    @classmethod
    def from_dense(cls, J, *, h=None, **kw) -> "CouplingSpec":
        """From a symmetric ``N×N`` coupling matrix (diagonal ignored)."""
        J = np.asarray(J, dtype=float)
        if J.ndim != 2 or J.shape[0] != J.shape[1]:
            raise ValueError("J must be a square 2-D matrix")
        n = J.shape[0]
        adjacency: list[list[tuple[int, float]]] = [[] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j and J[i, j] != 0.0:
                    adjacency[i].append((j, float(J[i, j])))
        return cls(n_families=n, adjacency=adjacency,
                   h=np.zeros(n) if h is None else h, **kw)

    @classmethod
    def from_edges(cls, n_families: int, edges, *, h=None, **kw) -> "CouplingSpec":
        """From a sparse edge list: ``edges`` maps ``(i, j) -> J_ij`` (or is an iterable of
        ``(i, j, J_ij)``). Edges are symmetrised; a repeated ``(i, j)`` is summed."""
        adjacency: list[list[tuple[int, float]]] = [[] for _ in range(n_families)]
        acc: dict[tuple[int, int], float] = {}
        norm = (((i, j), v) for (i, j), v in edges.items()) if isinstance(edges, dict) \
            else (((i, j), v) for i, j, v in edges)
        for (i, j), v in norm:
            if i == j:
                raise ValueError(f"self-coupling ({i},{j}) is not allowed (zero diagonal)")
            if not (0 <= i < n_families and 0 <= j < n_families):
                raise ValueError(f"edge ({i},{j}) out of range for n_families={n_families}")
            key = (i, j) if i < j else (j, i)
            acc[key] = acc.get(key, 0.0) + float(v)
        for (i, j), v in acc.items():
            adjacency[i].append((j, v))
            adjacency[j].append((i, v))
        return cls(n_families=n_families, adjacency=adjacency,
                   h=np.zeros(n_families) if h is None else h, **kw)

    # --- introspection -----------------------------------------------------
    def dense_J(self) -> np.ndarray:
        """Materialise the dense ``N×N`` coupling matrix (for inspection / plots)."""
        J = np.zeros((self.n_families, self.n_families), dtype=float)
        for i, nbrs in enumerate(self.adjacency):
            for j, v in nbrs:
                J[i, j] = v
        return J


def pathway_blocks(pathway_sizes, *, within: float = 3.0, between: float = 0.0,
                   h=0.0, **kw) -> CouplingSpec:
    """Build a :class:`CouplingSpec` whose ``J`` has pathway block structure — the hybrid
    (c) of the design note.

    Families are laid out in contiguous blocks (``pathway_sizes = [s0, s1, ...]``). Within a
    block every pair is coupled by ``within`` (positive → co-occurring pathway members).
    ``between`` couples every *cross-block* pair (use a negative value for mutually-exclusive
    "rival" pathways; the default ``0`` leaves blocks independent). ``h`` is a scalar field
    applied to every family (or a length-``N`` vector).
    """
    sizes = list(pathway_sizes)
    n = sum(sizes)
    block_of = np.empty(n, dtype=int)
    start = 0
    for b, s in enumerate(sizes):
        block_of[start:start + s] = b
        start += s
    J = np.full((n, n), between, dtype=float)
    for i in range(n):
        for j in range(n):
            if block_of[i] == block_of[j]:
                J[i, j] = within
    np.fill_diagonal(J, 0.0)
    hv = np.full(n, float(h)) if np.isscalar(h) else np.asarray(h, dtype=float)
    return CouplingSpec.from_dense(J, h=hv, **kw)


# ═══════════════════════════════════════════════════════════════════════════════
# The rate model
# ═══════════════════════════════════════════════════════════════════════════════
class PottsRates(RateModel):
    """Coupled loss + horizontal-transfer gain over a fixed family panel.

    For each *present* family ``i`` emits a loss weight ``base_loss·exp(-β·f_i)`` (so
    partners modulate retention), plus one field-blind ``TRANSFER`` gain channel and an
    optional background ``ORIGINATION``. Families outside the panel (e.g. originated ones)
    are treated as uncoupled — field ``0``, loss ``base_loss``.

    Weights depend only on the discrete presence vector, so they change solely at events;
    ``time_dependent`` stays ``False`` and the simulator refreshes only the branch that just
    changed — genuine Glauber-style dynamics along the tree.
    """

    def __init__(self, spec: CouplingSpec):
        self.spec = spec

    # --- local field -------------------------------------------------------
    def _loss_rate(self, idx: int, present: set[int]) -> float:
        """``base_loss · exp(-β·f_i)`` for panel family ``idx`` given the present set."""
        s = self.spec
        f_i = s.h[idx]
        for j, j_ij in s.adjacency[idx]:
            if j in present:
                f_i += j_ij
        expo = max(-_MAX_EXPONENT, min(_MAX_EXPONENT, -s.beta * f_i))
        return s.base_loss * math.exp(expo)

    def event_weights(self, genome, branch, time):
        s = self.spec
        index = s.index
        fams = genome.families()

        # presence vector σ restricted to the panel (non-panel families do not enter fields)
        present = {index[f] for f in fams if f in index}

        out: list[EventWeight] = []
        for fam in fams:
            idx = index.get(fam)
            rate = self._loss_rate(idx, present) if idx is not None else s.base_loss
            if rate > 0.0:
                out.append(EventWeight(EventType.LOSS, fam, rate))

        n = genome.size()
        if s.transfer > 0.0 and n > 0:
            out.append(EventWeight(EventType.TRANSFER, None, s.transfer * n))
        if s.origination > 0.0:
            out.append(EventWeight(EventType.ORIGINATION, None, s.origination))
        return out


# ═══════════════════════════════════════════════════════════════════════════════
# Driving a coupled simulation
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class CoupledResult:
    """Output of :func:`simulate_coupled`."""

    profiles: ProfileMatrix              # panel families × extant species (all panel rows kept)
    leaf_genomes: dict                   # extant leaf TreeNode -> final genome
    event_log: object                    # the full DTL event log
    spec: CouplingSpec = field(repr=False, default=None)


def _panel_factory(seed_families):
    """A ``genome_factory`` that seeds a genome with exactly ``seed_families`` present."""
    def factory(ids):
        g = UnorderedGenome(ids)
        for fam in seed_families:
            g._add(Gene(ids.new_gene(), fam))
        return g
    return factory


def _panel_profile(leaf_genomes, spec: CouplingSpec) -> ProfileMatrix:
    """Build the profile over the *full* fixed panel, keeping all-absent panel rows.

    (``ProfileMatrix.from_leaf_genomes`` would drop families absent from every species, but
    a benchmark wants the known ``N`` rows to survive so ground-truth couplings line up.)
    """
    species_nodes = sorted(leaf_genomes, key=lambda n: _natkey(n.name))
    species = [n.name for n in species_nodes]
    index = spec.index
    rows, cols, data = [], [], []
    for j, node in enumerate(species_nodes):
        genome = leaf_genomes[node]
        for fam in genome.families():
            k = index.get(fam)
            if k is None:
                continue  # background (originated) family — not part of the panel
            cn = genome.copy_number(fam)
            if cn:
                rows.append(k); cols.append(j); data.append(cn)
    return ProfileMatrix(families=list(spec.panel_ids), species=species,
                         coo=(rows, cols, data))


def simulate_coupled(tree, spec: CouplingSpec, *, seed=None, rng=None,
                     transfers: TransferModel | None = None,
                     initial_presence=None) -> CoupledResult:
    """Simulate coupled gene families along ``tree`` under ``spec``.

    The root genome is seeded with the panel (all families present by default; pass
    ``initial_presence`` as a length-``N`` 0/1 mask to start from a chosen configuration).
    Transfers default to full **replacement** (``TransferModel(replacement=1.0)``) so a
    re-acquired family does not stack copies — keeping the state cleanly presence/absence.

    Returns a :class:`CoupledResult`; ``.profiles`` is the ``N × species`` panel matrix.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    if initial_presence is None:
        seed_families = list(spec.panel_ids)
    else:
        mask = np.asarray(initial_presence)
        if mask.shape != (spec.n_families,):
            raise ValueError(f"initial_presence must have shape ({spec.n_families},)")
        seed_families = [fam for fam, on in zip(spec.panel_ids, mask) if on]

    tm = transfers if transfers is not None else TransferModel(replacement=1.0)
    result = GenomeSimulator().simulate(
        tree, PottsRates(spec), rng,
        initial_size=0, transfers=tm, genome_factory=_panel_factory(seed_families),
    )
    return CoupledResult(
        profiles=_panel_profile(result.leaf_genomes, spec),
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        spec=spec,
    )
