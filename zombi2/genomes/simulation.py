"""High-level gene-family simulation on a fixed species tree."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from zombi2._sampling import EventSampler
from zombi2.genomes.events import EventType
from zombi2.genomes.genome import UnorderedGenome
from zombi2.genomes.genome_sim import GenomeSimulator
from zombi2.genomes.profiles import ProfileMatrix
from zombi2.genomes.reconciliation import build_gene_trees
from zombi2.genomes.rates import RateModel, SharedRates
from zombi2.tree import Tree


#: Header of the compact single-file event trace (``Events_trace.tsv``): one row per event
#: (O/D/T/L/S), the scalable alternative to the per-family ``gene_family_events/`` directory.
EVENTS_TRACE_HEADER = ("time\tevent\tbranch\tdonor\trecipient\tfamily\tparent\tchild1\tchild2")

#: Header of the per-species-branch event summary (``Branch_events.tsv``): one row per branch of
#: the species tree, with the count of each gene-family event that fired on it. Transfers are split
#: into ``transfer_out`` (this branch is the donor) and ``transfer_in`` (the recipient). ``is_extant``
#: flags the branches that survive to the present, so a table restricted to the *extant* tree is
#: just this file filtered to ``is_extant == 1``.
BRANCH_EVENTS_HEADER = (
    "branch\ttime\tis_leaf\tis_extant\torigination\tduplication\ttransfer_in\ttransfer_out"
    "\tloss\tinversion\ttransposition\ttotal"
)

#: Event kinds that are counted, per firing branch, in ``Branch_events.tsv`` (transfers are handled
#: separately because they touch two branches). The value is the column name.
_BRANCH_EVENT_COLUMNS = {
    EventType.ORIGINATION: "origination",
    EventType.DUPLICATION: "duplication",
    EventType.LOSS: "loss",
    EventType.INVERSION: "inversion",
    EventType.TRANSPOSITION: "transposition",
}

#: Every count column (order-free); ``_TOTAL_COLUMNS`` is the subset summed into ``total`` — the
#: events that fired *on* the branch (``transfer_in`` is excluded: it fired on the donor).
_COUNTED_COLUMNS = ("origination", "duplication", "transfer_in", "transfer_out", "loss",
                    "inversion", "transposition")
_TOTAL_COLUMNS = ("origination", "duplication", "transfer_out", "loss", "inversion", "transposition")


def _extant_branches(species_tree) -> set:
    """Names of the branches on the *extant* species tree — those ancestral to a present-day leaf.

    A leaf is present-day when it sits at the tree's maximum depth (``total_age``); an internal
    branch is extant when any descendant leaf is. This is derived from node *times*, so it is
    correct even for a tree loaded from a plain Newick (which carries no extant flag) — a leaf that
    died before the present sits above ``total_age``'s depth and is not counted. In an ultrametric
    (fully-sampled) tree every leaf is present-day, so every branch is extant.
    """
    total_age = species_tree.total_age
    tol = 1e-9 * (abs(total_age) or 1.0)
    extant: set = set()
    # reversed preorder visits every node after all its descendants (a valid post-order)
    for node in reversed(list(species_tree.nodes_preorder())):
        if node.is_leaf():
            if abs(node.time - total_age) <= tol:
                extant.add(node.name)
        elif any(c.name in extant for c in node.children):
            extant.add(node.name)
    return extant


def branch_events_table(event_log, species_tree) -> str:
    """Aggregate an :class:`~zombi2.events.EventLog` into per-species-branch event counts.

    One row per branch of ``species_tree`` — *every* branch, including those on which nothing
    happened — carrying how many of each event *fired on* that branch: originations,
    duplications, losses and (for ordered genomes) inversions and transpositions, plus
    ``transfer_out`` (the branch was the transfer donor). ``transfer_in`` counts transfers the
    branch *received* (it was the recipient). ``total`` is the number of events that fired on the
    branch — it therefore includes ``transfer_out`` but not ``transfer_in`` (that event fired on
    the donor). ``is_extant`` marks branches reaching the present, so the per-branch table for the
    *extant* species tree is this file filtered to ``is_extant == 1``.

    Counts follow the event log's own granularity, so a multi-gene transfer on an ordered genome
    contributes one row per gene — matching ``Transfers.tsv``. The compact trace path routes here
    through :meth:`Genomes.write`, so speciation markers are already reinstated.
    """
    counts: dict[str, dict[str, int]] = {}

    def bump(branch: str, column: str) -> None:
        row = counts.get(branch)
        if row is None:
            row = counts[branch] = dict.fromkeys(_COUNTED_COLUMNS, 0)
        row[column] += 1

    for r in event_log:
        if r.event is EventType.TRANSFER:
            bump(r.branch, "transfer_out")
            if r.recipient is not None:
                bump(r.recipient, "transfer_in")
        else:
            column = _BRANCH_EVENT_COLUMNS.get(r.event)
            if column is not None:
                bump(r.branch, column)

    zero = dict.fromkeys(_COUNTED_COLUMNS, 0)
    extant = _extant_branches(species_tree)
    rows = [BRANCH_EVENTS_HEADER]
    for n in species_tree.nodes_preorder():
        c = counts.get(n.name, zero)
        total = sum(c[k] for k in _TOTAL_COLUMNS)
        rows.append(
            f"{n.name}\t{n.time:.10g}\t{int(n.is_leaf())}\t{int(n.name in extant)}\t"
            f"{c['origination']}\t{c['duplication']}\t{c['transfer_in']}\t{c['transfer_out']}\t"
            f"{c['loss']}\t{c['inversion']}\t{c['transposition']}\t{total}"
        )
    return "\n".join(rows) + "\n"


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


def _write_reconciliations(out: Path, recons: dict) -> None:
    """Write annotated reconciled gene trees: ``Reconciled_complete.nwk`` / ``Reconciled_extant.nwk``
    (one family per line — tips ``<species>|<gid>``, internal labels ``<species-branch>|<EVENT>``,
    the format ``tools recon-accuracy`` and ``tools reconcile`` read) plus a flat S/D/T/L event
    table. Mirrors ``tools simulate``'s ground-truth output so a simulated truth can be scored."""
    complete, extant = [], []
    ev = ["family\tevent\tspecies\trecipient\ttime\tgene"]
    for fam, recon in recons.items():
        if recon.complete is not None:
            complete.append(recon.complete)
        if recon.extant is not None:
            extant.append(recon.extant)
        for e in recon.events:
            ev.append(f"{fam}\t{e.event}\t{e.species}\t{e.recipient or ''}\t"
                      f"{e.time:.10g}\t{e.gene or ''}")
    (out / "Reconciled_complete.nwk").write_text("\n".join(complete) + ("\n" if complete else ""))
    (out / "Reconciled_extant.nwk").write_text("\n".join(extant) + ("\n" if extant else ""))
    (out / "Reconciliation_events.tsv").write_text("\n".join(ev) + "\n")


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


#: Header of ``Geneorder_events.tsv`` — the structural-event log with *positional* columns, one
#: row per event across every branch (filter on ``branch`` for the per-branch view). ``chrom /
#: start / length / strand`` are the physical arc the event acted on (half-open ``[start,
#: start+length)`` on the acting genome, circular), empty for events with no region (e.g. a plain
#: origination). Emitted only by the ordered / nucleotide models, which carry an event ``region``.
GENEORDER_EVENTS_HEADER = (
    "time\tevent\tbranch\tfamily\tchrom\tstart\tlength\tstrand\tdonor\trecipient")


def geneorder_events_from_log(event_log) -> str:
    """Serialise an :class:`~zombi2.events.EventLog` as ``Geneorder_events.tsv`` text.

    Native zombi2 coordinates: the ``region`` on each record is the physical arc of the operation
    (``[start, start+length)``, half-open, circular). Inversions and losses are fully described by
    their arc; transposition/duplication record the *source* arc (the paste destination is not yet
    logged — see docs/design/geneorder-export.md). Records without a region (origination, transfer
    on the family path) serialise with empty positional cells.
    """
    rows = [GENEORDER_EVENTS_HEADER]
    for r in event_log:
        if r.event.value in ("S", "F"):
            continue  # speciation / leaf markers are not gene-order events (tree is in the .nwk)
        reg = r.region
        chrom, start, length, strand = (
            ("", "", "", "") if reg is None
            else (reg.chromosome, reg.start, reg.length, reg.strand))
        family = r.genes[0].family if r.genes else ""
        rows.append(f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t{family}\t"
                    f"{chrom}\t{start}\t{length}\t{strand}\t{r.donor or ''}\t{r.recipient or ''}")
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
    from zombi2.genomes.events import EventRecord, EventType, GeneOp

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
        from zombi2.genomes.reconciliation import expand_trace
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
        from zombi2.genomes.reconciliation import reconcile

        gid2species = self._gid_to_species()
        total_age = self.species_tree.total_age
        return {
            fam: reconcile(records, gid2species, total_age)
            for fam, records in self.gene_families.items()
        }

    def reconciliation_likelihoods(self, duplication: float, transfer: float, loss: float, *,
                                   models=("dated", "undated"), origination: str = "root",
                                   n_steps: int = 100, backend: str = "auto"):
        """ALE reconciliation log-likelihood of every extant family's gene tree (ALElite).

        Scores each family's **extant** reconciled gene tree under each ``model`` (any of
        ``"dated"``, ``"undated"``, ``"reldated"``) at the given DTL rates — for ``dated`` the
        rates are per-unit-time δ/τ/λ (as simulated); for the undated models they are per-branch
        odds. Returns a list of :class:`~zombi2.tools.reconciliation.FamilyScore` (family, extant
        copy count, ``{model: loglik}``). Pass the simulation's own rates to score each family
        under the generating model. This is a convenience wrapper over the
        :mod:`zombi2.tools.reconciliation` tool (ALElite); the tool's own
        ``score_reconciliations`` / ``reconciliation_likelihood`` are the primary surface.
        """
        from zombi2.tools.reconciliation import score_reconciliations

        return score_reconciliations(
            self.species_tree, self.reconciliations(), duplication, transfer, loss,
            models=models, origination=origination, n_steps=n_steps, backend=backend,
        )

    # Selectable components for write(include=...). species_tree.nwk + species_nodes.tsv
    # are always written; the CLI's --write maps onto these names.
    WRITE_PARTS = ("profiles", "trace", "trees", "events", "transfers", "summary", "branch_events",
                   "reconciliations", "layout", "karyotype")
    #: Written only when explicitly requested (never by ``include=None``): the ordered-genome
    #: ``layout``/``karyotype`` are meaningful only for ordered genomes, and ``reconciliations``
    #: reconstructs the (expensive) gene-tree genealogy, so they stay out of a run's output unless
    #: asked for — default (single-chromosome) output is therefore unchanged.
    _OPT_IN_PARTS = frozenset({"layout", "karyotype", "reconciliations"})

    # --- output ------------------------------------------------------------
    def write(self, outdir: str | Path, *, include=None, sparse: bool = False,
              annotate_species: bool = False) -> None:
        """Write the ZOMBI-1-style output folder.

        ``include`` selects which components to write — any subset of :attr:`WRITE_PARTS`.
        ``None`` (the default) writes them all **except** the opt-in ordered-genome parts
        ``"layout"`` (``Gene_order.tsv`` — the per-leaf chromosome layout) and ``"karyotype"``
        (``Karyotype_trace.tsv`` — the fission/fusion/origination/loss genealogy), which are written
        only when asked for, so single-chromosome output is unchanged. ``species_tree.nwk`` and
        ``species_nodes.tsv`` are always written. Omitted components do no work — notably ``"trees"``
        drives the (expensive) gene-tree reconstruction.

        ``"trace"`` writes the compact single-file event log ``Events_trace.tsv`` (one row per
        event); it is the scalable alternative to ``"events"`` (per-family files) and is the
        record from which gene trees can be reconstructed later on demand. See
        :class:`GenomeTrace` for the fast path that produces it without materialising the log.

        ``"branch_events"`` writes ``Branch_events.tsv`` — the per-species-branch event counts
        (D/T/L/O and, for ordered genomes, inversion/transposition), one row per branch with an
        ``is_extant`` flag so the extant-tree view is a filter (see :func:`branch_events_table`).

        With ``sparse=True`` the copy-number profile is written as a single sparse long
        table (``Profiles_sparse.tsv``) instead of the dense ``Profiles.tsv`` /
        ``Presence.tsv`` pair — the latter is ``families × species`` (O(N²) in tip count)
        and is the one output that does not scale; the sparse form is O(present cells).
        """
        want = (set(self.WRITE_PARTS) - self._OPT_IN_PARTS) if include is None else set(include)
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

        if "reconciliations" in want:
            # annotated reconciled gene trees + event table — the input `tools recon-accuracy` and
            # `tools reconcile` read (tips <species>|<gid>). Same format as `tools simulate`, so a
            # simulated ground truth can be scored against an inferred reconciliation.
            _write_reconciliations(out, self.reconciliations())

        if "transfers" in want:
            # all transfers, one row each
            tr_lines = ["time\tfamily\tdonor_branch\trecipient_branch\tparent_id\tdonor_copy_id\ttransfer_id"]
            for r in self.event_log:
                if r.event is EventType.TRANSFER:
                    p, dc, tc = (op.gid for op in r.genes)
                    tr_lines.append(f"{r.time:.10g}\t{r.family}\t{r.donor}\t{r.recipient}\t{p}\t{dc}\t{tc}")
            (out / "Transfers.tsv").write_text("\n".join(tr_lines) + "\n")

        if "branch_events" in want:
            # per-species-branch event counts (see branch_events_table); is_extant column
            # makes the extant-tree view a filter.
            (out / "Branch_events.tsv").write_text(
                branch_events_table(self.event_log, self.species_tree))

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

        if "layout" in want:
            # per-leaf ordered layout: which chromosome each gene sits on, and in what order. The
            # karyotype is state that is otherwise unserialised; only ordered genomes have a layout,
            # so other genome models contribute no rows.
            lay = ["species\tchromosome\tposition\tfamily\tgid\torientation"]
            for leaf, genome in sorted(self.leaf_genomes.items(), key=lambda kv: kv[0].name):
                chroms = getattr(genome, "chromosomes", None)
                if not isinstance(chroms, dict):
                    continue
                for chrom in chroms.values():
                    for pos, g in enumerate(chrom.genes):
                        strand = "+" if getattr(g, "orientation", 1) >= 0 else "-"
                        lay.append(f"{leaf.name}\t{chrom.chrom_id}\t{pos}\t"
                                   f"{g.family}\t{g.gid}\t{strand}")
            (out / "Gene_order.tsv").write_text("\n".join(lay) + "\n")

        if "karyotype" in want:
            # the chromosome-tier genealogy: fission / fusion / chromosome origination / loss, each
            # with its source (parents) and resulting (children) chromosome ids. Empty (header only)
            # when the karyotype never changed.
            kar = ["time\tevent\tbranch\tparents\tchildren"]
            for r in self.event_log.chromosome_records:
                parents = ";".join(str(p) for p in r.parents)
                children = ";".join(str(c) for c in r.children)
                kar.append(f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t{parents}\t{children}")
            (out / "Karyotype_trace.tsv").write_text("\n".join(kar) + "\n")


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
            from zombi2._rust import _build_event_log
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
        from zombi2.genomes.events import EventType
        log = self.event_log
        if any(r.event is EventType.SPECIATION for r in log):   # full log (Python path)
            self._genomes_cache = Genomes(self.species_tree, self.leaf_genomes, log, self.profiles)
        else:                                                    # compact trace → expand
            from zombi2.genomes.events import EventLog
            from zombi2.genomes.reconciliation import expand_trace, extant_species_from_records
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
                from zombi2._rust import events_trace_tsv
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
    conversion: float = 0.0,
    initial_families: int = 20,
    transfers=None,
    conversions=None,
    max_family_size=None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
    genome_factory=UnorderedGenome,
    output: str = "genomes",
    threads: int = 1,
):
    """Simulate gene families forward along ``species_tree``.

    ``initial_families`` is the number of gene families seeded at the root (default 20).

    Provide either a rate model (``rates=z.SharedRates(...)`` /
    ``z.FamilySampledRates(...)``) or the convenience shorthand
    (``duplication=..., transfer=..., loss=..., origination=..., conversion=...``), which builds a
    :class:`~zombi2.SharedRates`. ``transfers`` (a :class:`~zombi2.TransferModel`) sets
    transfer mechanics; ``max_family_size`` (int absolute, or float as a multiple of the
    number of species) bounds family growth across duplication and transfer.

    ``conversion`` is the per-copy intra-genome gene-conversion rate (a copy of a family overwrites
    another copy of the same family; concerted evolution), and ``conversions`` sets its donor
    directionality — any object with a ``bias`` attribute (a :class:`~zombi2.ConversionModel`).
    Conversion has an effect only when ``conversion > 0``; ``conversions`` without it is inert.

    Engine (automatic, not a user choice): the **built-in** model — the default
    ``UnorderedGenome`` with a plain :class:`~zombi2.SharedRates` — runs on the compiled
    ``zombi2_core`` Rust extension and **requires** it (a clear error asks you to build it if
    it is missing). Flexible models (``FamilySampledRates`` / ``PerGenomeRates`` /
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

    ``threads`` (only with ``output="profiles"`` on the built-in Rust model, binary tree): run the
    counts-only engine on ``threads`` cores by Poisson-thinning the gene families across ``threads``
    independent copies and summing the profiles. Distributionally identical to the serial run (a
    different but equivalent realization); the output depends on ``(seed, threads)`` but not on
    scheduling. ``threads=1`` (default) is the exact serial engine.
    """
    if output not in ("genomes", "profiles", "trace"):
        raise ValueError(f"output must be 'genomes', 'profiles' or 'trace', got {output!r}")
    if threads > 1 and output != "profiles":
        raise ValueError("threads>1 (parallel simulation) is only supported for output='profiles'")

    # Tree-shape contract: every engine (the Rust kernel and the pure-Python engine) assumes a
    # bifurcating species tree. Degree-two nodes (FBD sampled ancestors) are handled by the Python
    # engine's pass-through, but a polytomy (a node with >2 children) is supported by no engine: the
    # Rust log/trace kernel underflows its alive-list index (a process-killing panic) and the Python
    # engine's _speciate unpacks exactly two children. Reject it here with a clear, actionable error
    # instead of a deep crash.
    _polytomies = [n.name for n in species_tree.nodes_preorder() if len(n.children) > 2]
    if _polytomies:
        shown = ", ".join(str(x) for x in _polytomies[:5]) + ("..." if len(_polytomies) > 5 else "")
        raise ValueError(
            f"simulate_genomes requires a bifurcating species tree, but found {len(_polytomies)} "
            f"polytomous node(s) with >2 children: {shown}. Resolve polytomies into bifurcations "
            "(e.g. as zero-length internal branches) before simulating."
        )

    shorthand = any((duplication, transfer, loss, origination, conversion))
    if rates is None:
        rates = SharedRates(duplication, transfer, loss, origination, conversion=conversion)
    elif shorthand:
        raise ValueError(
            "pass a rate model OR the duplication/transfer/loss/origination/conversion shorthand, "
            "not both"
        )

    from zombi2 import _rust

    # The Rust engine assumes a strictly binary species tree. Degree-two nodes (FBD sampled
    # ancestors, from forward simulation with removal < 1) are handled by the Python engine's
    # pass-through instead — Rust cannot process them, so this is a capability boundary, not a
    # silent engine preference (each tree has exactly one engine that can run it).
    binary_tree = all(len(n.children) != 1 for n in species_tree.nodes_preorder())
    if _rust.eligible(rates, genome_factory, sampler) and binary_tree:
        _rust.require()  # one engine for the built-in model on a binary tree; no Python fallback
        # Resolve an integer seed for the Rust engine when only an rng is supplied. Do this
        # inside the Rust branch (not before engine selection): the Python path below draws
        # directly from ``rng``, so deriving a seed there would consume a draw it never uses,
        # making ``seed=x`` and ``rng=default_rng(x)`` diverge for Python-engine models.
        if seed is None and rng is not None:
            seed = int(rng.integers(0, 2**63 - 1))
        kw = dict(initial_size=initial_families, transfers=transfers,
                  max_family_size=max_family_size, seed=seed)
        if output == "profiles":
            if threads > 1:
                return _rust.profiles_parallel(species_tree, rates, threads=threads, **kw)
            return _rust.profiles(species_tree, rates, **kw)
        if output == "trace":
            return _rust.trace(species_tree, rates, **kw)
        return _rust.genomes(species_tree, rates, **kw)

    if threads > 1:
        raise ValueError("threads>1 (parallel profiles) requires the built-in model "
                         "(UnorderedGenome + SharedRates) on a binary species tree")

    # --- flexible models: pure-Python engine -----------------------------------
    if rng is None:
        rng = np.random.default_rng(seed)
    result = GenomeSimulator(sampler).simulate(
        species_tree, rates, rng, initial_size=initial_families, transfers=transfers,
        conversions=conversions, max_family_size=max_family_size, genome_factory=genome_factory,
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
