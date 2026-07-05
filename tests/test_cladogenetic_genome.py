"""Tests for cladogenetic gene-family dynamics (``zombi2.cladogenetic_genome``: species:genes).

Change is concentrated **at speciation**, so the model is checked by its expected signal:

* no change (all rates 0) -> every tip carries exactly the root genome;
* pure cladogenetic *loss* -> every genome is a shrinking subset of the root;
* a cladogenetic burst -> sister (cherry) tips differ, because change happens at their split;
* gains grow the family universe beyond the root; reproducible under a fixed seed.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.cli import main


def _tree(seed=1, tips=40):
    return z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=tips, age=5, seed=seed)


# --------------------------------------------------------------------------- validation
def test_cladogenetic_genome_validation():
    with pytest.raises(ValueError):
        z.CladogeneticGenome(initial_families=-1)
    with pytest.raises(ValueError):
        z.CladogeneticGenome(30, loss=-0.1)
    with pytest.raises(ValueError):
        z.CladogeneticGenome(30, cladogenetic_loss=1.5)           # not a probability
    with pytest.raises(ValueError):
        z.CladogeneticGenome(30, cladogenetic_gain=-1)


def test_cladogenetic_genome_reproducible():
    tree = _tree()
    m = z.CladogeneticGenome(30, cladogenetic_loss=0.15, cladogenetic_gain=3)
    a = z.simulate_cladogenetic_genome(tree, m, seed=9)
    b = z.simulate_cladogenetic_genome(tree, m, seed=9)
    assert a.profile_matrix().to_tsv(presence=True) == b.profile_matrix().to_tsv(presence=True)


# --------------------------------------------------------------------------- deterministic corners
def test_no_change_gives_root_genome_at_every_tip():
    tree = _tree()
    m = z.CladogeneticGenome(20, loss=0, origination=0, cladogenetic_loss=0, cladogenetic_gain=0)
    res = z.simulate_cladogenetic_genome(tree, m, seed=1)
    root = res.node_genomes[tree.root]
    assert root == frozenset(range(20))
    assert all(res.node_genomes[t] == root for t in tree.extant_leaves())


def test_pure_cladogenetic_loss_only_shrinks():
    """No gains, no anagenesis: every node genome is a subset of the root and never larger."""
    tree = _tree()
    m = z.CladogeneticGenome(30, loss=0, origination=0, cladogenetic_loss=0.2, cladogenetic_gain=0)
    res = z.simulate_cladogenetic_genome(tree, m, seed=3)
    root = res.node_genomes[tree.root]
    assert all(res.node_genomes[n] <= root for n in tree.nodes())
    sizes = res.genome_sizes()
    assert all(sizes[n] <= len(root) for n in tree.nodes())


def test_gains_grow_the_family_universe():
    """With cladogenetic gain, far more families exist than the root held."""
    tree = _tree()
    m = z.CladogeneticGenome(20, cladogenetic_loss=0.1, cladogenetic_gain=3)
    res = z.simulate_cladogenetic_genome(tree, m, seed=2)
    assert len(res.profile_matrix().families) > 20


# --------------------------------------------------------------------------- the cladogenetic signature
def test_sister_tips_differ_because_change_is_at_the_split():
    tree = _tree()
    m = z.CladogeneticGenome(30, cladogenetic_loss=0.15, cladogenetic_gain=3)
    res = z.simulate_cladogenetic_genome(tree, m, seed=2)
    cherries = [n for n in tree.internal_nodes()
                if len(n.children) == 2 and all(not c.children for c in n.children)]
    differ = sum(res.node_genomes[n.children[0]] != res.node_genomes[n.children[1]] for n in cherries)
    assert cherries and differ >= 0.7 * len(cherries)          # nearly all sisters differ


def test_anagenetic_only_still_evolves():
    """With no cladogenetic change but non-zero anagenetic rates, genomes still turn over."""
    tree = _tree()
    m = z.CladogeneticGenome(20, loss=0.4, origination=0.4, cladogenetic_loss=0, cladogenetic_gain=0)
    res = z.simulate_cladogenetic_genome(tree, m, seed=4)
    tips = tree.extant_leaves()
    assert not all(res.node_genomes[t] == frozenset(range(20)) for t in tips)   # not frozen


# --------------------------------------------------------------------------- result views
def test_profile_matrix_is_presence_only():
    tree = _tree(tips=20)
    res = z.simulate_cladogenetic_genome(tree, z.CladogeneticGenome(15, cladogenetic_gain=2), seed=1)
    pm = res.profile_matrix()
    assert pm.shape[1] == len(tree.extant_leaves())            # one column per extant tip
    assert set(np.unique(pm.copy_values())) <= {1}            # presence/absence -> copies are 1


# --------------------------------------------------------------------------- CLI (species:genes edge)
def test_cli_species_genes(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--tips", "30", "--age", "4", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "sg"
    rc = main(["coevolve", "--couple", "species:genes", "-t", str(sp / "species_tree.nwk"),
               "--genome-size", "25", "--clado-gene-loss", "0.15", "--clado-gene-gain", "3",
               "--seed", "2", "-o", str(out)])
    assert rc == 0
    for f in ("species_tree.nwk", "Profiles.tsv", "Presence.tsv", "genome_sizes.tsv"):
        assert (out / f).exists()
    # genome_sizes covers every node; the root has exactly --genome-size families
    rows = dict(line.split("\t") for line in (out / "genome_sizes.tsv").read_text().splitlines()[1:])
    assert rows["root"] == "25"


def test_cli_species_genes_needs_tree(tmp_path):
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "species:genes", "-o", str(tmp_path / "a")])


def test_cli_species_genes_rejects_age(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "species:genes", "-t", str(sp / "species_tree.nwk"),
              "--age", "3", "-o", str(tmp_path / "a")])
