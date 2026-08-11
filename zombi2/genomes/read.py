"""Read a written genome run back in — the Python side of ``zombi2 sequences --from DIR``.

A run's files are the handoff between levels: the CLI writes them at one level and reopens them at
the next, and ``stream_to=`` exists precisely so a million families never have to fit in memory at
once. Reopening was CLI-only, so from Python that handoff was a dead end — you could write a streamed
run and then not feed it to the sequence level without shelling out. This is that reader, made
public.
"""

from __future__ import annotations

import os

from ..tree import Tree, read_newick
from .events import edges_from_tsv, gene_from_label
from .family import FamilyGenomesResult, GeneCopy

#: Where a genomes run's files sit inside a run directory: the grouped layout first, then ``--flat``.
_LAYOUTS = ("genomes", "")


def _resolve(directory: str) -> str:
    """The directory that actually holds ``genome_events.tsv``, in either layout."""
    for candidate in _LAYOUTS:
        full = os.path.join(directory, candidate) if candidate else directory
        if os.path.exists(os.path.join(full, "genome_events.tsv")):
            return full
    raise FileNotFoundError(
        f"{directory} holds no genomes run — looked for genome_events.tsv in "
        f"{os.path.join(directory, 'genomes')}/ and in {directory} itself. Run the genomes level "
        f"there first, and keep 'events' among its outputs (it is written by default).")


def _tree(handoff: str, directory: str) -> Tree:
    """The species tree the events index against. A genomes run writes its own canonical copy beside
    them, so the run is self-contained; the species level's own file is the fallback for a
    hand-assembled directory."""
    for candidate in (os.path.join(handoff, "species_complete.nwk"),
                      os.path.join(directory, "species", "species_complete.nwk"),
                      os.path.join(directory, "species_complete.nwk")):
        if os.path.exists(candidate):
            with open(candidate, encoding="utf-8") as f:
                return read_newick(f.read())[0]
    raise FileNotFoundError(
        f"{handoff} has an event log but no species_complete.nwk beside it, and the events name "
        f"branches of a tree that is not there. Re-run the genomes level (it writes the tree it ran "
        f"on), or put the species tree back.")


def _genomes(handoff: str, tree: Tree) -> dict[int, tuple[GeneCopy, ...]]:
    """Every node's gene content from ``genomes.tsv``, or ``{}`` when that output was not written.

    Read rather than re-derived: the run already computed it, and ``profiles`` is a table people
    publish. An empty map is honest — `FamilyGenomesResult.profiles` says what is missing — and the
    gene trees, which the sequence level actually runs on, come from the events either way."""
    path = os.path.join(handoff, "genomes.tsv")
    if not os.path.exists(path):
        return {}
    ids = {label: i for i, label in tree.labels().items()}   # n<id> / e<id> → node id
    out: dict[int, list[GeneCopy]] = {i: [] for i in tree.nodes}
    with open(path, encoding="utf-8") as f:
        for lineno, raw in enumerate(f, 1):
            row = raw.rstrip("\n").split("\t")
            if lineno == 1 and row[:1] == ["lineage"]:
                continue
            if len(row) != 3:
                raise ValueError(f"{path}:{lineno}: expected 'lineage<TAB>family<TAB>copy', "
                                 f"got {raw.rstrip()!r}")
            lineage, family, copy = row
            if lineage not in ids:
                raise ValueError(f"{path}:{lineno}: lineage {lineage!r} is not a node of the species "
                                 f"tree beside it — are the two files from the same run?")
            out[ids[lineage]].append(GeneCopy(gene_from_label(copy), int(family)))
    return {i: tuple(copies) for i, copies in out.items()}


def read_run(directory) -> FamilyGenomesResult:
    """Reopen a genome run written to ``directory`` — the run object, from its files.

    Takes a run directory in either layout (``out/`` with a ``genomes/`` inside it, or a ``flat=True``
    directory), or a `StreamedRun` handle. The event log is the source of truth, so ``gene_trees``,
    ``events`` and the whole genealogy come back exactly; ``genomes`` and ``initial_genome`` come from
    their own tables when the run wrote them, and are empty otherwise.

    This is what ``simulate_sequences`` calls when handed a path, so the two-command CLI pipeline and
    a two-step Python one are the same pipeline::

        st = genomes.simulate_genomes_family(sp, ..., stream_to="run/")
        sequences.simulate_sequences(st, model=hky85(), length=1000, seed=1)   # reads it back

    The ``seed`` of a reopened run is the seed recorded on the handle when there is one, and ``None``
    from a bare directory: the files are the run, and a number that did not produce them would be a
    lie. It is in the run's own ``run.zombi2`` / ``genomes.log`` either way.
    """
    seed = getattr(directory, "seed", None)                 # a StreamedRun carries its own
    directory = getattr(directory, "directory", directory)
    directory = os.fspath(directory)
    handoff = _resolve(directory)
    if os.path.exists(os.path.join(handoff, "blocks.tsv")):
        raise NotImplementedError(
            f"{handoff} is a nucleotide genome run (it has blocks.tsv), which reads back through "
            f"zombi2.genomes.nucleotide.read_nucleotide_genomes(directory, tree) — it needs the "
            f"blocks as well as the events, so it takes the tree explicitly.")
    tree = _tree(handoff, directory)
    with open(os.path.join(handoff, "genome_events.tsv"), encoding="utf-8") as f:
        edges = edges_from_tsv(f.read())
    initial: tuple[GeneCopy, ...] = ()
    initial_path = os.path.join(handoff, "initial_genome.tsv")
    if os.path.exists(initial_path):
        with open(initial_path, encoding="utf-8") as f:
            rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
        initial = tuple(GeneCopy(gene_from_label(copy), int(family))
                        for family, copy in rows if family != "family")
    return FamilyGenomesResult(complete_tree=tree, node_genomes=_genomes(handoff, tree), edges=edges,
                               seed=seed, initial_genome=initial)


__all__ = ["read_run"]
