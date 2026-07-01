"""Top-level driver that ties the species tree, gene families and profiles together."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._sampling import EventSampler
from .genome_sim import GenomeResult, GenomeSimulator
from .profiles import ProfileMatrix
from .rates import RateModel
from .species_model import SpeciesTreeModel
from .species_sim import SpeciesTreeSimulator
from .tree import Tree


@dataclass
class SimulationResult:
    """Everything a run produces."""

    species_tree: Tree
    genomes: GenomeResult
    profiles: ProfileMatrix

    @property
    def event_log(self):
        return self.genomes.event_log

    @property
    def gene_families(self):
        """Per-family event lists (family id -> list[EventRecord])."""
        return self.genomes.event_log.by_family()

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


class Simulation:
    """Run a full ZOMBI2 v1 simulation: backward species tree, then forward genomes."""

    def __init__(
        self,
        species_model: SpeciesTreeModel,
        rate_model: RateModel,
        *,
        seed: int | None = None,
        initial_size: int = 20,
        sampler: EventSampler | None = None,
        genome_factory=None,
    ):
        self.species_model = species_model
        self.rate_model = rate_model
        self.seed = seed
        self.initial_size = initial_size
        self.sampler = sampler
        self.genome_factory = genome_factory

    def run(self) -> SimulationResult:
        rng = np.random.default_rng(self.seed)
        tree = SpeciesTreeSimulator().simulate(self.species_model, rng)
        kwargs = {"initial_size": self.initial_size}
        if self.genome_factory is not None:
            kwargs["genome_factory"] = self.genome_factory
        genomes = GenomeSimulator(self.sampler).simulate(tree, self.rate_model, rng, **kwargs)
        profiles = ProfileMatrix.from_leaf_genomes(genomes.leaf_genomes)
        return SimulationResult(species_tree=tree, genomes=genomes, profiles=profiles)
