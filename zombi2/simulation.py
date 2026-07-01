"""High-level gene-family simulation on a fixed species tree."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._sampling import EventSampler
from .genome import UnorderedGenome
from .genome_sim import GenomeSimulator
from .profiles import ProfileMatrix
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

    def write(self, outdir: str | Path) -> None:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)

        (out / "species_tree.nwk").write_text(self.species_tree.to_newick() + "\n")

        node_lines = ["name\ttime\tis_leaf\tis_extant"]
        for n in self.species_tree.nodes_preorder():
            node_lines.append(f"{n.name}\t{n.time:.10g}\t{int(n.is_leaf())}\t{int(n.is_extant)}")
        (out / "species_nodes.tsv").write_text("\n".join(node_lines) + "\n")

        gdir = out / "gene_family_events"
        gdir.mkdir(exist_ok=True)
        for family, records in self.gene_families.items():
            lines = ["time\tevent\tbranch\tdonor\trecipient\tgenes"]
            for r in records:
                genes = ";".join(f"{op.role}:{op.gid}" for op in r.genes)
                lines.append(
                    f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t"
                    f"{r.donor or ''}\t{r.recipient or ''}\t{genes}"
                )
            (gdir / f"{family}_events.tsv").write_text("\n".join(lines) + "\n")

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
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
    genome_factory=UnorderedGenome,
) -> Genomes:
    """Simulate gene families forward along ``species_tree``.

    Provide either a rate model (``rates=z.UniformRates(...)`` /
    ``z.FamilySampledRates(...)``) or the convenience shorthand
    (``duplication=..., transfer=..., loss=..., origination=...``), which builds a
    :class:`~zombi2.UniformRates`.
    """
    shorthand = any((duplication, transfer, loss, origination))
    if rates is None:
        rates = UniformRates(duplication, transfer, loss, origination)
    elif shorthand:
        raise ValueError(
            "pass a rate model OR the duplication/transfer/loss/origination shorthand, not both"
        )

    if rng is None:
        rng = np.random.default_rng(seed)

    result = GenomeSimulator(sampler).simulate(
        species_tree, rates, rng, initial_size=initial_size, genome_factory=genome_factory
    )
    profiles = ProfileMatrix.from_leaf_genomes(result.leaf_genomes)
    return Genomes(
        species_tree=species_tree,
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        profiles=profiles,
    )
