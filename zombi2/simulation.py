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


@dataclass
class Genomes:
    """Result of :func:`simulate_genomes`."""

    species_tree: Tree
    leaf_genomes: dict  # extant leaf TreeNode -> its final genome
    event_log: object   # EventLog
    profiles: ProfileMatrix

    @property
    def gene_families(self):
        """Per-family event lists (family id -> list[EventRecord])."""
        return self.event_log.by_family()

    def _gid_to_species(self) -> dict[str, str]:
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

    # --- output ------------------------------------------------------------
    def write(self, outdir: str | Path, annotate_species: bool = False) -> None:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)

        (out / "species_tree.nwk").write_text(self.species_tree.to_newick() + "\n")

        node_lines = ["name\ttime\tis_leaf\tis_extant"]
        for n in self.species_tree.nodes_preorder():
            node_lines.append(f"{n.name}\t{n.time:.10g}\t{int(n.is_leaf())}\t{int(n.is_extant)}")
        (out / "species_nodes.tsv").write_text("\n".join(node_lines) + "\n")

        families = self.gene_families

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

        # gene trees (complete + extant)
        tdir = out / "gene_trees"
        tdir.mkdir(exist_ok=True)
        for family, (complete, extant) in self.gene_trees(annotate_species).items():
            if complete:
                (tdir / f"{family}_complete.nwk").write_text(complete + "\n")
            if extant:
                (tdir / f"{family}_extant.nwk").write_text(extant + "\n")

        # all transfers, one row each
        tr_lines = ["time\tfamily\tdonor_branch\trecipient_branch\tparent_id\tdonor_copy_id\ttransfer_id"]
        for r in self.event_log:
            if r.event is EventType.TRANSFER:
                p, dc, tc = (op.gid for op in r.genes)
                tr_lines.append(f"{r.time:.10g}\t{r.family}\t{r.donor}\t{r.recipient}\t{p}\t{dc}\t{tc}")
        (out / "Transfers.tsv").write_text("\n".join(tr_lines) + "\n")

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

        (out / "Profiles.tsv").write_text(self.profiles.to_tsv())
        (out / "Presence.tsv").write_text(self.profiles.to_tsv(presence=True))


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
        fast route for ABC and large presence/absence datasets.
    """
    if output not in ("genomes", "profiles"):
        raise ValueError(f"output must be 'genomes' or 'profiles', got {output!r}")

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
    return Genomes(
        species_tree=species_tree,
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        profiles=profiles,
    )
