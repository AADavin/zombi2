"""Cladogenetic gene-family dynamics — the ``species:genes`` edge (speciation drives the genome).

Ordinary gene-family evolution changes a genome *along* its branches (anagenesis: gains and losses
accumulate gradually in time). Here **speciation itself** reorganises the genome: at every branching
event each daughter undergoes a *burst* of gene loss and gain — a founder effect / genome upheaval
concentrated at the split. This is the gene-content counterpart of the cladogenetic trait model
([`species:traits`](coevolution_models)), i.e. **punctuational genome evolution**: the ``species:genes``
arrow, where *speciation drives gene content*.

Because gene content does **not** feed back into diversification here, the tree is an *input* (a given
Newick) and this is an **overlay**: a genome is evolved down the fixed tree with

* **anagenetic** change along each branch — each present family is lost at rate ``loss``; brand-new
  families originate at rate ``origination`` (a Gillespie process on the branch);
* a **cladogenetic** burst at each speciation — every daughter independently drops each family it
  carries with probability ``cladogenetic_loss`` and gains ``Poisson(cladogenetic_gain)`` new
  families.

Set the anagenetic rates to zero for *pure* punctuational genomes (constant along branches, changing
only at speciation); leave them on for both. The result carries the genome at **every** node, from
which the extant presence/absence profile is built.

    import zombi2 as z
    tree = z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=40, age=5, seed=1)
    m = z.CladogeneticGenome(initial_families=30, cladogenetic_loss=0.15, cladogenetic_gain=3)
    res = z.simulate_cladogenetic_genome(tree, m, seed=2)
    res.profile_matrix().to_tsv(presence=True)   # families x extant species (0/1)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from zombi2.genomes.profiles import ProfileMatrix
from zombi2.tree import Tree


class CladogeneticGenome:
    """Punctuational gene-family dynamics on a given tree (the ``species:genes`` edge).

    Parameters
    ----------
    initial_families:
        Number of distinct families in the root genome (``>= 0``).
    loss:
        Anagenetic per-family loss rate along a branch (``>= 0``).
    origination:
        Anagenetic origination rate of brand-new families, per lineage (``>= 0``).
    cladogenetic_loss:
        Probability, in ``[0, 1]``, that a daughter drops each family it carries **at speciation**
        (the founder-effect burst).
    cladogenetic_gain:
        Mean number of brand-new families a daughter gains **at speciation** (Poisson; ``>= 0``).
    """

    def __init__(self, initial_families: int = 30, *, loss: float = 0.0, origination: float = 0.0,
                 cladogenetic_loss: float = 0.1, cladogenetic_gain: float = 2.0):
        if initial_families < 0:
            raise ValueError(f"initial_families must be >= 0, got {initial_families}")
        if loss < 0 or origination < 0:
            raise ValueError("loss and origination must be >= 0")
        if not (0.0 <= cladogenetic_loss <= 1.0):
            raise ValueError(f"cladogenetic_loss must be a probability in [0, 1], got {cladogenetic_loss}")
        if cladogenetic_gain < 0:
            raise ValueError(f"cladogenetic_gain must be >= 0, got {cladogenetic_gain}")
        self.initial_families = int(initial_families)
        self.loss = float(loss)
        self.origination = float(origination)
        self.cladogenetic_loss = float(cladogenetic_loss)
        self.cladogenetic_gain = float(cladogenetic_gain)

    def __repr__(self):
        return (f"CladogeneticGenome(initial_families={self.initial_families}, "
                f"cladogenetic_loss={self.cladogenetic_loss:g}, "
                f"cladogenetic_gain={self.cladogenetic_gain:g})")


@dataclass
class CladogeneticGenomeResult:
    """The outcome of :func:`simulate_cladogenetic_genome`.

    ``node_genomes`` maps **every** node to the (frozen) set of family ids present there — the tips
    are the observable profile; internal nodes are the exact ancestral gene content.
    """

    tree: Tree
    model: CladogeneticGenome
    node_genomes: dict

    def genome_sizes(self) -> dict:
        """Number of families at each node."""
        return {node: len(g) for node, g in self.node_genomes.items()}

    def profile_matrix(self) -> ProfileMatrix:
        """A families × extant-species presence :class:`~zombi2.ProfileMatrix` (copies are 0/1).

        Only families present in at least one extant tip are kept (a family that is lost before the
        present, or lives only in an extinct lineage, is not observable)."""
        tips = self.tree.extant_leaves()
        present = sorted(set().union(*[self.node_genomes[t] for t in tips])) if tips else []
        index = {fam: i for i, fam in enumerate(present)}
        matrix = np.zeros((len(present), len(tips)), dtype=np.int64)
        for j, tip in enumerate(tips):
            for fam in self.node_genomes[tip]:
                matrix[index[fam], j] = 1
        return ProfileMatrix([f"F{fam}" for fam in present], [t.name for t in tips], matrix)


# --------------------------------------------------------------------------- engine
def _evolve_branch(genome: set, dt: float, loss: float, origination: float,
                   next_id: list, rng) -> set:
    """Anagenetic loss + origination along a branch of duration ``dt`` (exact Gillespie).

    ``next_id`` is a one-element list used as a shared, monotonically increasing family-id counter.
    """
    if loss <= 0.0 and origination <= 0.0:
        return genome
    elapsed = 0.0
    while True:
        n = len(genome)
        rate = loss * n + origination
        if rate <= 0.0:
            break
        elapsed += rng.exponential(1.0 / rate)
        if elapsed >= dt:
            break
        if n > 0 and rng.random() < (loss * n) / rate:          # a loss
            fams = list(genome)
            genome.discard(fams[int(rng.integers(n))])
        else:                                                   # a new family originates
            genome.add(next_id[0])
            next_id[0] += 1
    return genome


def _cladogenetic_burst(genome: set, clado_loss: float, clado_gain: float,
                        next_id: list, rng) -> set:
    """A daughter's genome at its birth: drop each family w.p. ``clado_loss`` and gain
    ``Poisson(clado_gain)`` brand-new families (the speciation burst)."""
    if clado_loss > 0.0:
        genome = {f for f in genome if rng.random() >= clado_loss}
    else:
        genome = set(genome)
    if clado_gain > 0.0:
        for _ in range(int(rng.poisson(clado_gain))):
            genome.add(next_id[0])
            next_id[0] += 1
    return genome


def simulate_cladogenetic_genome(
    tree: Tree,
    model: CladogeneticGenome,
    *,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
) -> CladogeneticGenomeResult:
    """Evolve a genome down ``tree`` with anagenetic change along branches and a cladogenetic burst
    at each speciation (the ``species:genes`` edge — speciation drives gene content).

    ``tree`` is any :class:`~zombi2.tree.Tree` (a simulated species tree or one read with
    :func:`~zombi2.read_newick`). Returns a :class:`CladogeneticGenomeResult` holding the genome at
    every node; ``.profile_matrix()`` is the observable extant presence/absence profile.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    root = tree.root
    next_id = [model.initial_families]                          # next fresh family id
    node_genomes = {root: frozenset(range(model.initial_families))}
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        start = set(node_genomes[node.parent])                 # parent's end genome, inherited
        start = _cladogenetic_burst(start, model.cladogenetic_loss, model.cladogenetic_gain,
                                    next_id, rng)               # burst as this branch is born
        end = _evolve_branch(start, node.branch_length(), model.loss, model.origination,
                             next_id, rng)                      # then gradual change along it
        node_genomes[node] = frozenset(end)

    return CladogeneticGenomeResult(tree=tree, model=model, node_genomes=node_genomes)
