"""High-level gene-family simulation on a fixed species tree."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._sampling import EventSampler
from .events import EventType
from .genome import UnorderedGenome
from .genome_sim import GenomeSimulator
from .profiles import ProfileMatrix
from .reconciliation import build_gene_trees
from .rates import RateModel, UniformRates
from .tree import Tree


#: Header of the compact single-file event trace (``Events_trace.tsv``): one row per event
#: (O/D/T/L/S), the scalable alternative to the per-family ``gene_family_events/`` directory.
EVENTS_TRACE_HEADER = ("time\tevent\tbranch\tdonor\trecipient\tfamily\tparent\tchild1\tchild2")


def _write_tree_and_nodes(out: Path, tree: Tree) -> None:
    """Write the always-present ``species_tree.nwk`` + ``species_nodes.tsv`` pair."""
    (out / "species_tree.nwk").write_text(tree.to_newick() + "\n")
    node_lines = ["name\ttime\tis_leaf\tis_extant"]
    for n in tree.nodes_preorder():
        node_lines.append(f"{n.name}\t{n.time:.10g}\t{int(n.is_leaf())}\t{int(n.is_extant)}")
    (out / "species_nodes.tsv").write_text("\n".join(node_lines) + "\n")


def _write_profiles(out: Path, profiles: ProfileMatrix, sparse: bool) -> None:
    """Write the copy-number profile: one sparse long table, or the dense matrix pair."""
    if sparse:
        (out / "Profiles_sparse.tsv").write_text(profiles.to_coo_tsv())
    else:
        (out / "Profiles.tsv").write_text(profiles.to_tsv())
        (out / "Presence.tsv").write_text(profiles.to_tsv(presence=True))


def events_trace_from_log(event_log) -> str:
    """Serialise an :class:`~zombi2.events.EventLog` as the compact ``Events_trace.tsv`` text.

    One row per event; the gene-lineage ids are the record's ``parent -> child1[,child2]``
    gids (empty when the event has no children). This is the record-based writer; the Rust
    fast path emits the same format straight from its columns via
    :func:`zombi2._rust.events_trace_tsv` (no per-event Python objects)."""
    rows = [EVENTS_TRACE_HEADER]
    for r in event_log:
        g = r.genes
        parent = g[0].gid
        child1 = g[1].gid if len(g) > 1 else ""
        child2 = g[2].gid if len(g) > 2 else ""
        rows.append(f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t{r.donor or ''}\t"
                    f"{r.recipient or ''}\t{r.family}\t{parent}\t{child1}\t{child2}")
    return "\n".join(rows) + "\n"


def read_events_trace(text: str, species_tree=None) -> dict:
    """Parse ``Events_trace.tsv`` text back into ``{family: [EventRecord]}`` (time-ordered).

    The inverse of :func:`events_trace_from_log`. Each row rebuilds one
    :class:`~zombi2.events.EventRecord`; only the gene-lineage genealogy is recovered (``GeneOp``
    roles are not stored in the trace, and gene-tree / sequence reconstruction does not need
    them). This is the from-disk entry point behind ``zombi2 sequence`` — gene trees stay
    reconstructable on demand from the compact trace.

    The compact ``output="trace"`` file carries **no speciation rows** (a lineage keeps its id
    across speciations). When ``species_tree`` is supplied and the file is compact, the trace is
    replayed against it (:func:`~zombi2.reconciliation.expand_trace`) so the returned records are
    the full O/S/D/T/L genealogy the ordinary reconstruction expects. A file that already
    contains speciation rows (a full log written as a trace) is returned as-is.
    """
    from .events import EventRecord, EventType, GeneOp

    families: dict[str, list] = {}
    has_speciation = False
    lines = text.splitlines()
    if not lines:
        return families
    if lines[0] != EVENTS_TRACE_HEADER:
        raise ValueError(f"not an Events_trace.tsv (header was {lines[0]!r})")
    for line in lines[1:]:
        if not line.strip():
            continue
        time, event, branch, donor, recipient, family, parent, child1, child2 = line.split("\t")
        genes = [GeneOp(parent, family, "parent")]
        if child1:
            genes.append(GeneOp(child1, family, "child"))
        if child2:
            genes.append(GeneOp(child2, family, "child"))
        ev = EventType(event)
        has_speciation = has_speciation or ev is EventType.SPECIATION
        rec = EventRecord(ev, branch, float(time), genes,
                          donor=donor or None, recipient=recipient or None)
        families.setdefault(family, []).append(rec)

    if species_tree is not None and not has_speciation:
        from .reconciliation import expand_trace
        families = expand_trace(families, species_tree)
    return families


@dataclass
class Genomes:
    """Result of :func:`simulate_genomes`."""

    species_tree: Tree
    leaf_genomes: dict  # extant leaf TreeNode -> its final genome
    event_log: object   # EventLog
    profiles: ProfileMatrix
    # Optional precomputed gid -> extant species. Used when the genealogy was reconstructed by
    # replaying a compact trace (no leaf genomes to read it off); otherwise it is derived from
    # ``leaf_genomes``. See :meth:`GenomeTrace.genomes`.
    gid2species: dict | None = None

    @property
    def gene_families(self):
        """Per-family event lists (family id -> list[EventRecord])."""
        return self.event_log.by_family()

    def _gid_to_species(self) -> dict[str, str]:
        if self.gid2species is not None:
            return self.gid2species
        out: dict[str, str] = {}
        for leaf, genome in self.leaf_genomes.items():
            for g in genome.genes():
                out[g.gid] = leaf.name
        return out

    def gene_trees(self, annotate_species: bool = False) -> dict[str, tuple[str, str | None]]:
        """Reconstruct ``{family: (complete_newick, extant_newick)}`` from the event log.

        ``extant_newick`` is ``None`` for families with no surviving copies. With
        ``annotate_species=True`` internal gene nodes are labelled ``<gid>|<species-branch>``.
        """
        gid2species = self._gid_to_species()
        total_age = self.species_tree.total_age
        return {
            fam: build_gene_trees(records, gid2species, total_age, annotate_species)
            for fam, records in self.gene_families.items()
        }

    def reconciliations(self) -> dict:
        """Reconcile each family against the species tree.

        Returns ``{family: Reconciliation(complete, extant, events)}`` — the complete gene
        tree annotated with its species mapping (real events + losses), the extant (pruned)
        tree annotated (cherries, no losses), and the S/D/T/L event list. See
        :func:`~zombi2.reconciliation.reconcile`.
        """
        from .reconciliation import reconcile

        gid2species = self._gid_to_species()
        total_age = self.species_tree.total_age
        return {
            fam: reconcile(records, gid2species, total_age)
            for fam, records in self.gene_families.items()
        }

    # Selectable components for write(include=...). species_tree.nwk + species_nodes.tsv
    # are always written; the CLI's --output maps onto these names.
    WRITE_PARTS = ("profiles", "trace", "trees", "events", "transfers", "summary")

    # --- output ------------------------------------------------------------
    def write(self, outdir: str | Path, *, include=None, sparse: bool = False,
              annotate_species: bool = False) -> None:
        """Write the ZOMBI-1-style output folder.

        ``include`` selects which components to write — any subset of :attr:`WRITE_PARTS`
        (``"profiles"``, ``"trace"``, ``"trees"``, ``"events"``, ``"transfers"``,
        ``"summary"``); ``None`` (the default) writes them all. ``species_tree.nwk`` and
        ``species_nodes.tsv`` are always written. Omitted components do no work — notably
        ``"trees"`` drives the (expensive) gene-tree reconstruction.

        ``"trace"`` writes the compact single-file event log ``Events_trace.tsv`` (one row per
        event); it is the scalable alternative to ``"events"`` (per-family files) and is the
        record from which gene trees can be reconstructed later on demand. See
        :class:`GenomeTrace` for the fast path that produces it without materialising the log.

        With ``sparse=True`` the copy-number profile is written as a single sparse long
        table (``Profiles_sparse.tsv``) instead of the dense ``Profiles.tsv`` /
        ``Presence.tsv`` pair — the latter is ``families × species`` (O(N²) in tip count)
        and is the one output that does not scale; the sparse form is O(present cells).
        """
        want = set(self.WRITE_PARTS) if include is None else set(include)
        unknown = want - set(self.WRITE_PARTS)
        if unknown:
            raise ValueError(f"unknown write component(s) {sorted(unknown)}; "
                             f"choose from {self.WRITE_PARTS}")

        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        _write_tree_and_nodes(out, self.species_tree)

        if "trace" in want:
            (out / "Events_trace.tsv").write_text(events_trace_from_log(self.event_log))

        # per-family event lists — needed only by the events and summary tables
        families = self.gene_families if (want & {"events", "summary"}) else {}

        if "events" in want:
            # per-family event tables (O/D/T/S/L with from -> to lineage ids)
            gdir = out / "gene_family_events"
            gdir.mkdir(exist_ok=True)
            for family, records in families.items():
                lines = ["time\tevent\tbranch\tdonor\trecipient\tnodes"]
                for r in records:
                    nodes = ";".join(f"{op.role}={op.gid}" for op in r.genes)
                    lines.append(
                        f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t"
                        f"{r.donor or ''}\t{r.recipient or ''}\t{nodes}"
                    )
                (gdir / f"{family}_events.tsv").write_text("\n".join(lines) + "\n")

        if "trees" in want:
            # gene trees (complete + extant)
            tdir = out / "gene_trees"
            tdir.mkdir(exist_ok=True)
            for family, (complete, extant) in self.gene_trees(annotate_species).items():
                if complete:
                    (tdir / f"{family}_complete.nwk").write_text(complete + "\n")
                if extant:
                    (tdir / f"{family}_extant.nwk").write_text(extant + "\n")

        if "transfers" in want:
            # all transfers, one row each
            tr_lines = ["time\tfamily\tdonor_branch\trecipient_branch\tparent_id\tdonor_copy_id\ttransfer_id"]
            for r in self.event_log:
                if r.event is EventType.TRANSFER:
                    p, dc, tc = (op.gid for op in r.genes)
                    tr_lines.append(f"{r.time:.10g}\t{r.family}\t{r.donor}\t{r.recipient}\t{p}\t{dc}\t{tc}")
            (out / "Transfers.tsv").write_text("\n".join(tr_lines) + "\n")

        if "summary" in want:
            # per-family summary
            counts_of = lambda recs, ev: sum(1 for r in recs if r.event is ev)
            sum_lines = ["family\torigin_time\torigin_branch\tn_dup\tn_transfer\tn_loss"
                         "\tn_speciation\textant_copies\tspecies_present"]
            pmat = self.profiles
            fam_row = {f: i for i, f in enumerate(pmat.families)}
            # Per-family totals off the sparse profile, computed once (no dense N² array).
            copies_by_row = pmat.copies_per_family()
            present_by_row = pmat.presence_per_family()
            for family, records in families.items():
                origin = next((r for r in records if r.event is EventType.ORIGINATION), None)
                ot = f"{origin.time:.10g}" if origin else ""
                ob = origin.branch if origin else ""
                if family in fam_row:
                    i = fam_row[family]
                    extant_copies = int(copies_by_row[i])
                    species_present = int(present_by_row[i])
                else:
                    extant_copies = species_present = 0
                sum_lines.append(
                    f"{family}\t{ot}\t{ob}\t{counts_of(records, EventType.DUPLICATION)}\t"
                    f"{counts_of(records, EventType.TRANSFER)}\t{counts_of(records, EventType.LOSS)}\t"
                    f"{counts_of(records, EventType.SPECIATION)}\t{extant_copies}\t{species_present}"
                )
            (out / "Gene_family_summary.tsv").write_text("\n".join(sum_lines) + "\n")

        if "profiles" in want:
            _write_profiles(out, self.profiles, sparse)


@dataclass
class GenomeTrace:
    """A gene-family simulation kept as a compact *event trace*, not a materialised log.

    This is what ``simulate_genomes(..., output="trace")`` returns. The genealogy is held in
    its cheapest form — the engine's raw event columns (Rust path) or an already-built
    :class:`~zombi2.events.EventLog` (Python path) — so writing the scalable single-file
    ``Events_trace.tsv`` and the copy-number profile costs almost nothing: the per-event
    Python objects and the gene trees are built **only if you ask** (via :meth:`event_log` /
    :meth:`genomes` / :meth:`gene_trees`). It is the intermediate between the counts-only
    ``ProfileMatrix`` and the full :class:`Genomes`.
    """

    species_tree: Tree
    leaf_genomes: dict
    profiles: ProfileMatrix
    _columns: tuple | None = None   # raw Rust event columns (fast path); None on the Python path
    _nodes: list | None = None      # preorder node list aligned with the column indices
    _event_log: object | None = None  # prebuilt EventLog (Python path) or lazily materialised
    _genomes_cache: object = None   # promoted Genomes (see genomes())

    @property
    def event_log(self):
        """The :class:`~zombi2.events.EventLog` — materialised from the columns on first access
        (the deferred cost the trace exists to avoid). Cached thereafter.

        For the Rust ``output="trace"`` path this is the **compact** log: O/D/T/L only, no
        speciation records (see :meth:`genomes` for the reconstructed full genealogy)."""
        if self._event_log is None:
            if self._columns is None:
                raise RuntimeError("GenomeTrace has neither raw columns nor an event log")
            from ._rust import _build_event_log
            self._event_log = _build_event_log(self._columns, self._nodes)
        return self._event_log

    def genomes(self) -> "Genomes":
        """Promote to a full :class:`Genomes` (the deferred reconstruction cost). Cached.

        A compact trace (Rust ``output="trace"``: O/D/T/L only) is first **expanded** — replayed
        against the species tree to re-insert the implied speciations and remint per-instance
        gene ids (:func:`zombi2.reconciliation.expand_trace`) — so the resulting ``Genomes`` has a
        full O/S/D/T/L log and behaves exactly like a directly-simulated one. A log that already
        carries speciations (the Python-engine path) is wrapped as-is."""
        if self._genomes_cache is not None:
            return self._genomes_cache
        from .events import EventType
        log = self.event_log
        if any(r.event is EventType.SPECIATION for r in log):   # full log (Python path)
            self._genomes_cache = Genomes(self.species_tree, self.leaf_genomes, log, self.profiles)
        else:                                                    # compact trace → expand
            from .events import EventLog
            from .reconciliation import expand_trace, extant_species_from_records
            full = expand_trace(log.by_family(), self.species_tree)
            g2s = extant_species_from_records(full, self.species_tree)
            exp = EventLog()
            for recs in full.values():
                for r in recs:
                    exp.add(r)
            self._genomes_cache = Genomes(self.species_tree, {}, exp, self.profiles, gid2species=g2s)
        return self._genomes_cache

    def gene_trees(self, annotate_species: bool = False):
        """Reconstruct ``{family: (complete, extant)}`` on demand (see :meth:`Genomes.gene_trees`)."""
        return self.genomes().gene_trees(annotate_species)

    def reconciliations(self):
        """Reconcile each family on demand (see :meth:`Genomes.reconciliations`)."""
        return self.genomes().reconciliations()

    def write(self, outdir: str | Path, *, include=None, sparse: bool = False,
              annotate_species: bool = False) -> None:
        """Write the output folder for a trace run.

        ``include`` defaults to ``("trace", "profiles")`` — the two cheap components. Any of
        the heavier :attr:`Genomes.WRITE_PARTS` (``"trees"``, ``"events"``, ``"transfers"``,
        ``"summary"``) is honoured too, but reconstructing them materialises the event log
        first (:meth:`genomes`), so the fast path is ``include`` ⊆ ``{"trace", "profiles"}``.
        """
        want = {"trace", "profiles"} if include is None else set(include)
        unknown = want - set(Genomes.WRITE_PARTS)
        if unknown:
            raise ValueError(f"unknown write component(s) {sorted(unknown)}; "
                             f"choose from {Genomes.WRITE_PARTS}")

        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        _write_tree_and_nodes(out, self.species_tree)

        if "trace" in want:
            if self._columns is not None:  # fast: straight from the Rust columns, no objects
                from ._rust import events_trace_tsv
                (out / "Events_trace.tsv").write_text(events_trace_tsv(self._columns, self._nodes))
            else:
                (out / "Events_trace.tsv").write_text(events_trace_from_log(self.event_log))

        if "profiles" in want:
            _write_profiles(out, self.profiles, sparse)

        heavy = want - {"trace", "profiles"}
        if heavy:  # trees / events / transfers / summary — needs the materialised log
            self.genomes().write(out, include=heavy, sparse=sparse,
                                 annotate_species=annotate_species)


def simulate_genomes(
    species_tree: Tree,
    rates: RateModel | None = None,
    *,
    duplication: float = 0.0,
    transfer: float = 0.0,
    loss: float = 0.0,
    origination: float = 0.0,
    initial_size: int = 20,
    transfers=None,
    max_family_size=None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
    genome_factory=UnorderedGenome,
    output: str = "genomes",
):
    """Simulate gene families forward along ``species_tree``.

    Provide either a rate model (``rates=z.UniformRates(...)`` /
    ``z.FamilySampledRates(...)``) or the convenience shorthand
    (``duplication=..., transfer=..., loss=..., origination=...``), which builds a
    :class:`~zombi2.UniformRates`. ``transfers`` (a :class:`~zombi2.TransferModel`) sets
    transfer mechanics; ``max_family_size`` (int absolute, or float as a multiple of the
    number of species) bounds family growth across duplication and transfer.

    Engine (automatic, not a user choice): the **built-in** model — the default
    ``UnorderedGenome`` with a plain :class:`~zombi2.UniformRates` — runs on the compiled
    ``zombi2_core`` Rust extension and **requires** it (a clear error asks you to build it if
    it is missing). Flexible models (``FamilySampledRates`` / ``GenomeWiseRates`` /
    ``BranchRates``, ordered genomes, soft carrying capacity, custom samplers) run on the
    pure-Python engine.

    ``output``:
        ``"genomes"`` (default) returns a full :class:`Genomes` (event log, gene trees,
        ``write()``). ``"profiles"`` returns only a :class:`~zombi2.ProfileMatrix`; for the
        built-in model this takes the Rust counts-only path (no gene ids / log / trees) — the
        fast route for ABC and large presence/absence datasets. ``"trace"`` returns a
        :class:`GenomeTrace` — the compact event trace (``Events_trace.tsv``) plus the profile,
        without materialising the per-event Python objects or gene trees (built lazily, only if
        asked). It is the intermediate for very large datasets: near counts-only speed, yet
        gene trees remain reconstructable on demand.
    """
    if output not in ("genomes", "profiles", "trace"):
        raise ValueError(f"output must be 'genomes', 'profiles' or 'trace', got {output!r}")

    shorthand = any((duplication, transfer, loss, origination))
    if rates is None:
        rates = UniformRates(duplication, transfer, loss, origination)
    elif shorthand:
        raise ValueError(
            "pass a rate model OR the duplication/transfer/loss/origination shorthand, not both"
        )

    # Resolve an integer seed for the Rust engine when only an rng is supplied.
    if seed is None and rng is not None:
        seed = int(rng.integers(0, 2**63 - 1))

    from . import _rust

    # The Rust engine assumes a strictly binary species tree. Degree-two nodes (FBD sampled
    # ancestors, from forward simulation with removal < 1) are handled by the Python engine's
    # pass-through instead — Rust cannot process them, so this is a capability boundary, not a
    # silent engine preference (each tree has exactly one engine that can run it).
    binary_tree = all(len(n.children) != 1 for n in species_tree.nodes_preorder())
    if _rust.eligible(rates, genome_factory, sampler) and binary_tree:
        _rust.require()  # one engine for the built-in model on a binary tree; no Python fallback
        kw = dict(initial_size=initial_size, transfers=transfers,
                  max_family_size=max_family_size, seed=seed)
        if output == "profiles":
            return _rust.profiles(species_tree, rates, **kw)
        if output == "trace":
            return _rust.trace(species_tree, rates, **kw)
        return _rust.genomes(species_tree, rates, **kw)

    # --- flexible models: pure-Python engine -----------------------------------
    if rng is None:
        rng = np.random.default_rng(seed)
    result = GenomeSimulator(sampler).simulate(
        species_tree, rates, rng, initial_size=initial_size, transfers=transfers,
        max_family_size=max_family_size, genome_factory=genome_factory,
    )
    profiles = ProfileMatrix.from_leaf_genomes(result.leaf_genomes)
    if output == "profiles":
        return profiles
    if output == "trace":
        # the Python engine already built the log during simulation, so there is no object cost
        # to defer; wrap it so the return type and Events_trace.tsv output match the Rust path.
        return GenomeTrace(species_tree=species_tree, leaf_genomes=result.leaf_genomes,
                           profiles=profiles, _event_log=result.event_log)
    return Genomes(
        species_tree=species_tree,
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        profiles=profiles,
    )
