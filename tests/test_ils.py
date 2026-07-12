"""Incomplete lineage sorting via the multispecies coalescent (zombi2.experimental.ils).

Pure numpy -- the whole suite, including the end-to-end CLI run, works on a bare install.
The sharp correctness check is the rooted-triple concordance against the analytic
``1 - (2/3) e^{-T/N}`` (Hudson 1983), plus the two no-/full-ILS limits.
"""
from __future__ import annotations

import ast
import math

import numpy as np
import pytest

import zombi2.experimental as ex
from zombi2 import BirthDeath, simulate_species_tree
from zombi2.cli import main
from zombi2.experimental.ils import (
    MultispeciesCoalescent, _coalesce_bounded, _g_lineages_to_one,
    expected_triple_concordance, is_concordant, rooted_clades,
)
from zombi2.genomes.simulation import simulate_genomes
from zombi2.tree import TreeNode, read_newick

# rooted triple ((A,B),C) with internal branch length T and tips at age T + 0.5
_T = 0.5
_TRIPLE = read_newick(f"((A:0.5,B:0.5):{_T},C:{_T + 0.5});")


def _cherry(gene_tree):
    """The pair of leaves forming a cherry in a 3-taxon gene tree -- its topology label."""
    for n in gene_tree.nodes():
        if not n.is_leaf() and all(c.is_leaf() for c in n.children):
            return frozenset(c.name for c in n.children)
    return frozenset()


# --------------------------------------------------------------------------- #
# experimental lifecycle + conventions
# --------------------------------------------------------------------------- #
def test_multispeciescoalescent_is_experimental():
    ex._warned.discard("MultispeciesCoalescent")
    with pytest.warns(ex.ExperimentalWarning, match="MultispeciesCoalescent"):
        MultispeciesCoalescent(population_size=1.0)


def test_exports_stay_in_the_experimental_namespace():
    import zombi2
    for name in ("MultispeciesCoalescent", "expected_triple_concordance",
                 "is_concordant", "rooted_clades"):
        assert name in ex.__all__ and hasattr(ex, name)
        assert not hasattr(zombi2, name), f"{name} leaked into the top-level zombi2 namespace"


def test_module_pulls_no_heavy_dependencies():
    """The 'pure numpy, no optional deps' promise: the module imports only numpy/stdlib/zombi2."""
    import zombi2.experimental.ils as m
    with open(m.__file__) as f:
        tree = ast.parse(f.read())
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".")[0])
    assert roots <= {"__future__", "math", "numpy", "zombi2"}, f"unexpected imports: {roots}"


# --------------------------------------------------------------------------- #
# correctness: the multispecies coalescent
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("N", [0.5, 1.0, 2.0])
def test_triple_concordance_matches_theory(N):
    msc = MultispeciesCoalescent(population_size=N)
    trees = msc.sample_gene_trees(_TRIPLE, 6000, rng=np.random.default_rng(0))
    empirical = np.mean([is_concordant(g, _TRIPLE) for g in trees])
    assert abs(empirical - expected_triple_concordance(_T, N)) < 0.03


def test_discordant_topologies_are_equifrequent():
    msc = MultispeciesCoalescent(population_size=1.0)
    trees = msc.sample_gene_trees(_TRIPLE, 8000, rng=np.random.default_rng(1))
    ac = np.mean([_cherry(g) == frozenset("AC") for g in trees])
    bc = np.mean([_cherry(g) == frozenset("BC") for g in trees])
    assert abs(ac - bc) < 0.03                      # the two minor resolutions are equally likely


def test_no_ils_limit_recovers_the_species_tree():
    # N tiny relative to the branches: coalescence is forced at the nodes -> gene tree == species tree
    msc = MultispeciesCoalescent(population_size=1e-3)
    trees = msc.sample_gene_trees(_TRIPLE, 200, rng=np.random.default_rng(2))
    assert all(is_concordant(g, _TRIPLE) for g in trees)


def test_full_ils_limit_approaches_one_third():
    # N huge relative to the branches: the three rooted topologies become equiprobable
    msc = MultispeciesCoalescent(population_size=100.0)
    trees = msc.sample_gene_trees(_TRIPLE, 8000, rng=np.random.default_rng(3))
    concordant = np.mean([is_concordant(g, _TRIPLE) for g in trees])
    assert 0.30 < concordant < 0.37


# --------------------------------------------------------------------------- #
# the gene tree is a well-formed genealogy
# --------------------------------------------------------------------------- #
def test_gene_tree_is_a_wellformed_genealogy():
    species = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=6, age=1.0, seed=7)
    msc = MultispeciesCoalescent(population_size=0.5)
    g = msc.sample_gene_tree(species, rng=np.random.default_rng(0))

    sp_leaves = {lf.name for lf in species.leaves()}
    assert {lf.name for lf in g.leaves()} == sp_leaves          # one copy per species
    internal = g.internal_nodes()
    assert all(len(n.children) == 2 for n in internal)          # strictly bifurcating
    assert len(internal) == len(sp_leaves) - 1                  # n-1 coalescences
    assert all(n.branch_length() >= 0 for n in g.nodes() if n.parent)
    leaf_times = [lf.time for lf in g.leaves()]
    assert max(leaf_times) - min(leaf_times) < 1e-9            # contemporaneous tips (ultrametric)
    assert g.root.time < species.root.time                     # TMRCA predates species root (deep coalescence)


def test_samples_multiple_copies_per_tip():
    species = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=3, age=1.0, seed=5)
    msc = MultispeciesCoalescent(population_size=1.0)
    g = msc.sample_gene_tree(species, samples=2, rng=np.random.default_rng(0))
    names = sorted(lf.name for lf in g.leaves())
    assert len(names) == 6                                      # 2 copies x 3 species
    for lf in species.leaves():
        assert f"{lf.name}_1" in names and f"{lf.name}_2" in names


def test_reproducible_under_a_seed():
    msc = MultispeciesCoalescent(population_size=1.0)
    a = msc.sample_gene_tree(_TRIPLE, rng=np.random.default_rng(42)).to_newick()
    b = msc.sample_gene_tree(_TRIPLE, rng=np.random.default_rng(42)).to_newick()
    assert a == b


# --------------------------------------------------------------------------- #
# input validation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bad", [0.0, -1.0, float("inf"), float("nan")])
def test_population_size_must_be_positive_and_finite(bad):
    with pytest.raises(ValueError, match="population_size"):
        MultispeciesCoalescent(population_size=bad)


def test_sampling_zero_copies_everywhere_errors():
    msc = MultispeciesCoalescent(population_size=1.0)
    with pytest.raises(ValueError, match="at least one"):
        msc.sample_gene_tree(_TRIPLE, samples=0, rng=np.random.default_rng(0))


def test_rooted_clades_capture_topology():
    assert rooted_clades(_TRIPLE) == {frozenset("AB"), frozenset("ABC")}


# --------------------------------------------------------------------------- #
# CLI (torch-free -- no optional deps to guard)
# --------------------------------------------------------------------------- #
def test_ils_help_exits_zero():
    with pytest.raises(SystemExit) as e:
        main(["experimental", "ils", "-h"])
    assert e.value.code == 0


def test_ils_missing_population_size_errors(tmp_path):
    t = tmp_path / "sp.nwk"
    t.write_text(_TRIPLE.to_newick())
    with pytest.raises(SystemExit):                             # -N is required -> argparse exit 2
        main(["experimental", "ils", "-t", str(t), "-o", str(tmp_path / "o")])


def test_ils_replicates_must_be_positive(tmp_path, capsys):
    t = tmp_path / "sp.nwk"
    t.write_text(_TRIPLE.to_newick())
    rc = main(["experimental", "ils", "-t", str(t), "-N", "1.0", "-n", "0", "-o", str(tmp_path / "o")])
    assert rc == 1 and "replicates" in capsys.readouterr().err   # clean one-line error, not a traceback


def test_ils_end_to_end_writes_gene_trees(tmp_path, capsys):
    species = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=5, age=1.0, seed=1)
    t = tmp_path / "species_tree.nwk"
    t.write_text(species.to_newick())
    out = tmp_path / "out"
    rc = main(["experimental", "ils", "-t", str(t), "-N", "0.5", "-n", "50", "--seed", "1", "-o", str(out)])
    assert rc == 0
    gene_trees = (out / "gene_trees.nwk").read_text().splitlines()
    assert len(gene_trees) == 50
    assert (out / "species_tree.nwk").exists() and (out / "ils.log").exists()
    # each line is a parseable gene tree on the 5 species
    assert {lf.name for lf in read_newick(gene_trees[0]).leaves()} == {lf.name for lf in species.leaves()}
    assert "match the species-tree topology" in capsys.readouterr().out   # concordance diagnostic


# --------------------------------------------------------------------------- #
# v2: DTL + ILS -- the coalescent within each family's locus tree
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def dtl_run():
    """A gene-family simulation with real duplication / transfer / loss to exercise the locus tree."""
    sp = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=8, age=1.0, seed=3)
    g = simulate_genomes(sp, duplication=0.4, transfer=0.25, loss=0.4, origination=0.1,
                         initial_families=40, seed=7, output="genomes")
    return sp, g


# -- the bounded ("founder") coalescent primitive -----------------------------------------
@pytest.mark.parametrize("t", [0.2, 0.7, 1.5])
def test_bounded_coalescent_g_matches_tavare(t):
    n = 1.3
    assert abs(_g_lineages_to_one(2, t, n) - (1 - math.exp(-t / n))) < 1e-12   # g_{2,1}(t)=1-e^{-t/n}
    assert _g_lineages_to_one(4, 0.0, n) < 1e-12                               # can't coalesce in 0 time
    assert _g_lineages_to_one(4, 1e6, n) == pytest.approx(1.0)                 # always coalesced by t=inf


def test_bounded_coalescent_returns_a_single_founder():
    rng = np.random.default_rng(0)
    for k in range(1, 7):
        out = _coalesce_bounded([TreeNode(f"x{i}", 1.0) for i in range(k)], 1.0, 0.3, 1.0, rng)
        assert len(out) == 1                                   # always bottlenecks to one founder
        stack = [out[0]]                                       # every coalescence within the branch
        while stack:
            nd = stack.pop()
            assert nd.time >= 1.0 - 0.3 - 1e-9
            stack.extend(nd.children)


def test_bounded_coalescent_m2_time_is_truncated_exponential():
    s, n = 0.4, 1.0
    rng = np.random.default_rng(1)
    times = [1.0 - _coalesce_bounded([TreeNode("a", 1.0), TreeNode("b", 1.0)], 1.0, s, n, rng)[0].time
             for _ in range(20000)]
    lam = 1.0 / n
    theory = 1 / lam - s * math.exp(-lam * s) / (1 - math.exp(-lam * s))       # E[X | X <= s], X~Exp(lam)
    assert abs(float(np.mean(times)) - theory) < 0.01
    assert max(times) <= s + 1e-9                              # never exceeds the deadline


# -- the locus-tree coalescent ------------------------------------------------------------
def test_dtl_ils_reproduces_reconstruction_at_zero_ils(dtl_run):
    """The key check: N->0 collapses the locus-tree coalescent onto ZOMBI2's OWN deterministic
    gene-tree reconstruction -- across duplications, transfers, losses and originations."""
    _sp, g = dtl_run
    extant = {fam: ext for fam, (_c, ext) in g.gene_trees().items() if ext}
    mine = MultispeciesCoalescent(1e-4).sample_family_gene_trees(g, rng=np.random.default_rng(0))
    compared = 0
    for fam, ext in extant.items():
        z = read_newick(ext)
        if len(z.leaves()) < 2 or fam not in mine:            # single-copy families are trivially equal
            continue
        compared += 1
        assert rooted_clades(mine[fam]) == rooted_clades(z), f"family {fam} differs from reconstruction"
    assert compared >= 10                                     # the run really did exercise multi-copy families


def test_dtl_ils_creates_discordance_when_ils_is_on(dtl_run):
    _sp, g = dtl_run
    extant = {fam: ext for fam, (_c, ext) in g.gene_trees().items() if ext}
    mine = MultispeciesCoalescent(0.5).sample_family_gene_trees(g, rng=np.random.default_rng(0))
    discord = sum(1 for fam, ext in extant.items()
                  if fam in mine and len(read_newick(ext).leaves()) >= 3
                  and rooted_clades(mine[fam]) != rooted_clades(read_newick(ext)))
    assert discord > 0                                        # ILS makes gene trees disagree with the locus tree


def test_family_gene_trees_are_wellformed(dtl_run):
    _sp, g = dtl_run
    mine = MultispeciesCoalescent(0.5).sample_family_gene_trees(g, rng=np.random.default_rng(0))
    assert mine
    for tree in mine.values():
        assert all(len(nd.children) == 2 for nd in tree.internal_nodes())     # binary genealogy
        assert all("_" in lf.name for lf in tree.leaves())                    # named "<species>_<gid>"
        assert all(nd.branch_length() >= 0 for nd in tree.nodes() if nd.parent)


# -- the standalone CLI DTL + ILS mode (--events-trace) -----------------------------------
def test_ils_dtl_cli_end_to_end(tmp_path):
    sp = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=6, age=1.0, seed=2)
    t = tmp_path / "sp.nwk"
    t.write_text(sp.to_newick())
    grun = tmp_path / "grun"
    rc = main(["genomes", "-t", str(t), "--dup", "0.3", "--trans", "0.2", "--loss", "0.3",
               "--orig", "0.1", "--initial-families", "20", "--seed", "5", "--write", "trace", "-o", str(grun)])
    assert rc == 0 and (grun / "Events_trace.tsv").exists()

    out = tmp_path / "ilsout"
    rc = main(["experimental", "ils", "-t", str(t), "--events-trace", str(grun / "Events_trace.tsv"),
               "-N", "0.5", "-n", "2", "--seed", "1", "-o", str(out)])
    assert rc == 0
    files = list((out / "gene_trees").glob("*.nwk"))
    assert files                                              # a gene tree per surviving family
    lines = files[0].read_text().splitlines()
    assert len(lines) == 2                                    # -n 2 replicates per family
    assert read_newick(lines[0]).leaves()                    # parseable gene tree
