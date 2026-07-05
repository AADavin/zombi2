"""The built-in gene-family model runs on the compiled Rust engine (``zombi2_core``).

There is no separate "fast" function: ``simulate_genomes`` routes the built-in model
(``UnorderedGenome`` + ``UniformRates``) to Rust automatically — ``output="profiles"`` for
the counts-only path, the default ``output="genomes"`` for the full genealogy. Flexible
models (genome-wise rates, soft carrying capacity, …) fall through to the Python engine.
Skipped entirely if the extension isn't built.
"""

import numpy as np
import pytest

import zombi2 as z

pytestmark = pytest.mark.skipif(not z.rust_available(),
                                reason="zombi2_core (Rust extension) not built")

RATES = dict(duplication=0.15, transfer=0.1, loss=0.2, origination=0.5,
             initial_families=20, max_family_size=0.5)
FULL = dict(RATES)


def _tree(n=30, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n, age=5.0, seed=seed)


def _extant_leaves(newick):
    # gene trees are strictly binary, so #leaves = #commas + 1
    return 0 if newick is None else newick.count(",") + 1


def _read_profiles(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    return rows[0][1:], {r[0]: list(map(int, r[1:])) for r in rows[1:]}


def _python_profiles(tree, seed, *, initial_families, max_family_size, **rates):
    """Reference profiles from the pure-Python engine (the built-in model otherwise routes
    to Rust), for the cross-engine statistical-equivalence check."""
    from zombi2.genome_sim import GenomeSimulator

    rng = np.random.default_rng(seed)
    res = GenomeSimulator().simulate(tree, z.UniformRates(**rates), rng,
                                     initial_size=initial_families, max_family_size=max_family_size)
    return z.ProfileMatrix.from_leaf_genomes(res.leaf_genomes)


# --- output="profiles" (Rust counts-only path) -----------------------------------

def test_profiles_reproducible_same_seed():
    tree = _tree()
    a = z.simulate_genomes(tree, seed=7, output="profiles", **RATES)
    b = z.simulate_genomes(tree, seed=7, output="profiles", **RATES)
    assert np.array_equal(a.matrix, b.matrix)


def test_profiles_shape_and_species_are_extant_leaves():
    tree = _tree(n=30)
    pm = z.simulate_genomes(tree, seed=3, output="profiles", **RATES)
    assert pm.matrix.shape[1] == len(tree.extant_leaves())
    assert set(pm.species) == {n.name for n in tree.extant_leaves()}
    assert (pm.matrix >= 0).all()


def test_profiles_hard_cap_respected():
    tree = _tree(n=40)
    cap = 5
    pm = z.simulate_genomes(tree, duplication=0.6, transfer=0.2, loss=0.05,
                            origination=0.4, initial_families=20, max_family_size=cap,
                            seed=11, output="profiles")
    assert pm.matrix.max() <= cap


def test_profiles_accepts_uniform_rates_object():
    tree = _tree()
    obj = z.simulate_genomes(tree, z.UniformRates(0.15, 0.1, 0.2, 0.5),
                             initial_families=20, max_family_size=0.5, seed=7, output="profiles")
    kw = z.simulate_genomes(tree, seed=7, output="profiles", **RATES)
    assert np.array_equal(obj.matrix, kw.matrix)


def test_flexible_models_route_to_python_for_profiles():
    # These are no longer "rejected" — simulate_genomes runs them on the Python engine and
    # still returns a ProfileMatrix for output="profiles".
    tree = _tree()
    a = z.simulate_genomes(tree, z.UniformRates(0.2, 0, 0.1, 0.3, carrying_capacity=10),
                           seed=1, output="profiles")
    b = z.simulate_genomes(tree, z.GenomeWiseRates(0.2, 0.1, 0.2, 0.5), seed=1, output="profiles")
    for pm in (a, b):
        assert isinstance(pm, z.ProfileMatrix)
        assert pm.matrix.shape[1] == len(tree.extant_leaves())


# --- output="genomes" (Rust full-genealogy path) ---------------------------------

def test_full_log_reproducible():
    tree = _tree()
    a = z.simulate_genomes(tree, seed=7, **FULL)
    b = z.simulate_genomes(tree, seed=7, **FULL)
    assert len(a.event_log) == len(b.event_log)
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)


def test_full_log_gene_tree_invariant():
    # reconciliation invariant: extant gene-tree leaves == family's extant copy count
    tree = _tree(n=40, seed=3)
    g = z.simulate_genomes(tree, seed=5, **FULL)
    fam_row = {f: i for i, f in enumerate(g.profiles.families)}
    checked = 0
    for fam, (_complete, extant) in g.gene_trees().items():
        copies = int(g.profiles.matrix[fam_row[fam]].sum()) if fam in fam_row else 0
        assert _extant_leaves(extant) == copies, (fam, _extant_leaves(extant), copies)
        checked += 1
    assert checked > 0


def test_full_log_event_types_present():
    tree = _tree(n=40, seed=2)
    g = z.simulate_genomes(tree, seed=9, **FULL)
    kinds = {r.event.value for r in g.event_log}
    assert {"O", "S", "D", "L"}.issubset(kinds)  # transfers may be absent by chance


def test_full_log_gene_ids_are_g_prefixed():
    # schema parity with the Python engine: gene ids look like "g0", "g1", …
    tree = _tree()
    g = z.simulate_genomes(tree, seed=1, **FULL)
    gids = [op.gid for r in g.event_log for op in r.genes]
    assert gids and all(isinstance(x, str) and x.startswith("g") for x in gids)


def test_flexible_model_full_log_runs_on_python():
    # GenomeWiseRates is not the built-in model → Python engine, full Genomes (no raise).
    tree = _tree()
    g = z.simulate_genomes(tree, z.GenomeWiseRates(0.2, 0.1, 0.2, 0.5), seed=1)
    assert isinstance(g, z.Genomes)
    assert len(g.event_log) > 0


def test_profiles_vs_genomes_family_counts_agree():
    # the two Rust paths run the same model; mean family count should agree
    tree = _tree(n=60, seed=2)
    a = np.mean([len(z.simulate_genomes(tree, duplication=0.2, loss=0.1, origination=0.4,
                                        initial_families=30, max_family_size=0.5,
                                        seed=100 + s).profiles.families) for s in range(8)])
    b = np.mean([z.simulate_genomes(tree, duplication=0.2, loss=0.1, origination=0.4,
                                    initial_families=30, max_family_size=0.5, seed=200 + s,
                                    output="profiles").matrix.shape[0] for s in range(8)])
    assert abs(a - b) / b < 0.2


# --- write() on a Rust-backed Genomes --------------------------------------------

def test_write_produces_all_files(tmp_path):
    tree = _tree()
    z.simulate_genomes(tree, seed=7, **FULL).write(tmp_path)
    species, _ = _read_profiles(tmp_path / "Profiles.tsv")
    assert len(species) == len(tree.extant_leaves())
    for f in ["species_tree.nwk", "species_nodes.tsv", "Transfers.tsv",
              "Gene_family_summary.tsv", "Profiles.tsv", "Presence.tsv"]:
        assert (tmp_path / f).exists()
    assert (tmp_path / "gene_trees").is_dir()
    assert (tmp_path / "gene_family_events").is_dir()


def test_write_gene_tree_invariant(tmp_path):
    tree = _tree(n=40, seed=3)
    z.simulate_genomes(tree, seed=5, **FULL).write(tmp_path)
    _, prof = _read_profiles(tmp_path / "Profiles.tsv")
    checked = 0
    for fam, row in prof.items():
        rowsum = sum(row)
        ep = tmp_path / "gene_trees" / f"{fam}_extant.nwk"
        n_leaves = (ep.read_text().count(",") + 1) if ep.exists() else 0
        assert n_leaves == rowsum
        checked += 1
    assert checked > 0


def test_write_presence_matches_profiles(tmp_path):
    z.simulate_genomes(_tree(), seed=1, **FULL).write(tmp_path)
    _, prof = _read_profiles(tmp_path / "Profiles.tsv")
    _, pres = _read_profiles(tmp_path / "Presence.tsv")
    for fam in prof:
        assert pres[fam] == [1 if x > 0 else 0 for x in prof[fam]]


def test_write_reproducible(tmp_path):
    tree = _tree()
    z.simulate_genomes(tree, seed=9, **FULL).write(tmp_path / "a")
    z.simulate_genomes(tree, seed=9, **FULL).write(tmp_path / "b")
    assert (tmp_path / "a" / "Profiles.tsv").read_text() == (tmp_path / "b" / "Profiles.tsv").read_text()
    assert (tmp_path / "a" / "Transfers.tsv").read_text() == (tmp_path / "b" / "Transfers.tsv").read_text()


# --- transfer mechanics ----------------------------------------------------------

def _self_transfer_rows(path):
    with open(path) as f:
        next(f)
        return sum(1 for line in f if (lambda c: c[2] == c[3])(line.split("\t")))


def test_transfer_self(tmp_path):
    tree = _tree(n=30, seed=2)
    base = dict(duplication=0.0, transfer=0.4, loss=0.1, origination=0.4,
                initial_families=25, max_family_size=0.9, seed=5)
    z.simulate_genomes(tree, transfers=z.TransferModel(allow_self=False), **base).write(tmp_path / "no")
    z.simulate_genomes(tree, transfers=z.TransferModel(allow_self=True), **base).write(tmp_path / "yes")
    assert _self_transfer_rows(tmp_path / "no" / "Transfers.tsv") == 0
    assert _self_transfer_rows(tmp_path / "yes" / "Transfers.tsv") > 0


def test_transfer_replacement_reduces_growth():
    tree = _tree(n=30, seed=2)
    base = dict(transfer=0.4, loss=0.02, origination=0.4, initial_families=25, max_family_size=0.9)
    add = np.mean([z.simulate_genomes(tree, seed=200 + s, output="profiles",
                   transfers=z.TransferModel(replacement=0.0), **base).matrix.sum() for s in range(4)])
    rep = np.mean([z.simulate_genomes(tree, seed=200 + s, output="profiles",
                   transfers=z.TransferModel(replacement=1.0), **base).matrix.sum() for s in range(4)])
    assert rep < add


def test_transfer_distance_decay_prefers_close(tmp_path):
    tree = _tree(n=40, seed=2)
    nodes = list(tree.nodes_preorder())
    name2i = {n.name: i for i, n in enumerate(nodes)}
    par = [name2i[n.parent.name] if n.parent else -1 for n in nodes]
    tm = [n.time for n in nodes]

    def mean_dist(path):
        ds = []
        with open(path) as f:
            next(f)
            for line in f:
                c = line.split("\t")
                t, db, rb = float(c[0]), c[2], c[3]
                if db == rb:
                    ds.append(0.0)
                    continue
                anc, n = set(), name2i[db]
                while n >= 0:
                    anc.add(n)
                    n = par[n]
                n = name2i[rb]
                while n not in anc:
                    n = par[n]
                ds.append(2 * (t - tm[n]))
        return sum(ds) / len(ds) if ds else 0.0

    base = dict(transfer=0.4, loss=0.1, origination=0.4, initial_families=40, max_family_size=0.9, seed=11)
    z.simulate_genomes(tree, **base).write(tmp_path / "u")
    z.simulate_genomes(tree, transfers=z.TransferModel(distance_decay=5.0), **base).write(tmp_path / "d")
    assert mean_dist(tmp_path / "d" / "Transfers.tsv") < mean_dist(tmp_path / "u" / "Transfers.tsv")


def test_transfer_mechanics_invariant():
    # full-log invariant still holds with replacement + distance decay + self combined
    tree = _tree(n=40, seed=3)
    g = z.simulate_genomes(tree, duplication=0.05, transfer=0.25, loss=0.15, origination=0.5,
                           initial_families=30, max_family_size=0.5, seed=7,
                           transfers=z.TransferModel(replacement=0.3, distance_decay=2.0, allow_self=True))
    fr = {f: i for i, f in enumerate(g.profiles.families)}
    for fam, (_c, e) in g.gene_trees().items():
        leaves = 0 if e is None else e.count(",") + 1
        copies = int(g.profiles.matrix[fr[fam]].sum()) if fam in fr else 0
        assert leaves == copies


def test_rust_matches_python_engine():
    # the Rust engine should match the pure-Python reference within Monte-Carlo error
    tree = _tree(n=60, seed=2)
    kw = dict(duplication=0.2, loss=0.1, origination=0.4, initial_families=30, max_family_size=0.5)
    r = np.mean([z.simulate_genomes(tree, seed=1000 + s, output="profiles", **kw).matrix.mean()
                 for s in range(15)])
    p = np.mean([_python_profiles(tree, seed=2000 + s, **kw).matrix.mean() for s in range(15)])
    assert abs(r - p) / p < 0.15
