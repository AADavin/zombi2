"""SequenceEvolution — the gene x lineage substitution clock over reconciled gene trees.

rate(family g, species branch b) = R_b (shared lognormal lineage clock) x s_g (family speed).
Built on GenomeWiseRates so the fixture runs on the pure-Python engine (no Rust needed).
"""

import pytest

import zombi2 as z


def _sim():
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=1)
    g = z.simulate_genomes(
        tree, z.GenomeWiseRates(duplication=0.3, transfer=0.15, loss=0.2, origination=0.6),
        initial_size=12, seed=2)
    return tree, g


def _total_length(newick: str) -> float:
    t = z.read_newick(newick)
    return sum(n.branch_length() for n in t.nodes_preorder() if n.parent is not None)


def test_rate1_identity_is_the_chronogram():
    """sigma=0 and unit family speed -> the phylogram IS the input timetree, byte for byte."""
    _, g = _sim()
    ph = z.SequenceEvolution(branch_sigma=0.0, family_speed=1.0).scale(g, seed=5)
    chrono = g.gene_trees()
    checked = 0
    for fam, (_, extant) in chrono.items():
        if extant is None:
            continue
        checked += 1
        assert ph.extant[fam] == extant
    assert checked > 0


def test_global_scaling_multiplies_chronogram():
    """A flat lineage rate r and a fixed family speed c scale every branch by r*c."""
    _, g = _sim()
    ph = z.SequenceEvolution(branch_sigma=0.0, root_rate=2.0, family_speed=3.0).scale(g, seed=5)
    chrono = g.gene_trees()
    checked = 0
    for fam, (_, extant) in chrono.items():
        if extant is None:
            continue
        checked += 1
        assert _total_length(ph.extant[fam]) == pytest.approx(6.0 * _total_length(extant), rel=1e-4)
    assert checked > 0


def test_lineage_rates_shared_and_vary():
    """One lineage clock shared by all families (a rate per species branch); speeds differ."""
    tree, g = _sim()
    ph = z.SequenceEvolution(branch_sigma=0.7, family_speed=z.LogNormal(0.0, 0.4)).scale(g, seed=6)
    assert set(ph.branch_rate) == {n.name for n in tree.nodes_preorder()}
    assert len({round(v, 6) for v in ph.branch_rate.values()}) > 1   # lineage rates vary
    assert len({round(v, 6) for v in ph.family_speed.values()}) > 1  # family speeds vary


def test_autocorrelation_relates_parent_and_child_rates():
    """Adjacent species branches have correlated rates (the relaxed clock is autocorrelated)."""
    import numpy as np
    tree, g = _sim()
    ph = z.SequenceEvolution(branch_sigma=0.8).scale(g, seed=8)
    parent, child = [], []
    for node in tree.nodes_preorder():
        if node.parent is not None and node.parent.parent is not None:
            parent.append(ph.branch_rate[node.parent.name])
            child.append(ph.branch_rate[node.name])
    assert np.corrcoef(parent, child)[0, 1] > 0.3


def test_strict_clock_gives_flat_rates():
    _, g = _sim()
    ph = z.SequenceEvolution(branch_sigma=0.0, root_rate=1.5).scale(g, seed=1)
    assert all(v == 1.5 for v in ph.branch_rate.values())


def test_reproducible_given_seed():
    _, g = _sim()
    make = lambda: z.SequenceEvolution(branch_sigma=0.5, family_speed=z.LogNormal(0.0, 0.3))  # noqa: E731
    a = make().scale(g, seed=11)
    b = make().scale(g, seed=11)
    assert a.extant == b.extant
    assert a.family_speed == b.family_speed
    assert a.branch_rate == b.branch_rate


def test_discrete_bin_lineage_clock_runs():
    """A RateVariation (GTDB discrete-bin) lineage clock is selectable and integrates."""
    _, g = _sim()
    rv = z.RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=1.0)
    ph = z.SequenceEvolution(lineage=rv, family_speed=z.LogNormal(0.0, 0.4)).scale(g, seed=5)
    extant = [e for e in ph.extant.values() if e]
    assert extant and all(e.endswith(";") for e in extant)
    assert len({round(v, 6) for v in ph.branch_rate.values()}) > 1   # per-branch avg rates vary


def test_single_bin_rate1_is_the_chronogram():
    """A one-bin [1.0] strict clock with unit speed reproduces the timetree — segments integrate."""
    _, g = _sim()
    rv = z.RateVariation(bins=[1.0], switch_rate=0.0)
    ph = z.SequenceEvolution(lineage=rv, family_speed=1.0).scale(g, seed=3)
    chrono = g.gene_trees()
    checked = 0
    for fam, (_, extant) in chrono.items():
        if extant is None:
            continue
        checked += 1
        assert ph.extant[fam] == extant
    assert checked > 0


def test_two_lineage_clocks_rejected():
    with pytest.raises(ValueError):
        z.SequenceEvolution(branch_sigma=0.5, lineage=z.RateVariation(bins=[1.0], switch_rate=0.0))


def test_complete_phylograms_are_valid_newick():
    _, g = _sim()
    ph = z.SequenceEvolution(branch_sigma=0.3, family_speed=z.LogNormal(0.0, 0.3)).scale(g, seed=2)
    complete = [c for c in ph.complete.values() if c]
    assert complete and all(c.endswith(";") for c in complete)


def test_bad_parameters_rejected():
    with pytest.raises(ValueError):
        z.SequenceEvolution(branch_sigma=-0.1)
    with pytest.raises(ValueError):
        z.SequenceEvolution(root_rate=0.0)


# --- replaying from a written Events_trace.tsv (what `zombi2 sequence` does) ----------------

def test_scale_families_equals_scale():
    """The low-level scale_families entry point matches scale() byte-for-byte at full precision."""
    _, g = _sim()
    se = z.SequenceEvolution(branch_sigma=0.5, family_speed=z.LogNormal(0.0, 0.4))
    a = se.scale(g, seed=9)
    b = se.scale_families(g.species_tree, g.gene_families, g._gid_to_species(), seed=9)
    assert a.extant == b.extant and a.complete == b.complete
    assert a.family_speed == b.family_speed and a.branch_rate == b.branch_rate


def test_gid_to_species_recovered_from_trace():
    """Reconstructing gid->species from the written trace matches the live leaf-genome mapping."""
    from zombi2.reconciliation import extant_species_from_records
    from zombi2.simulation import events_trace_from_log, read_events_trace

    for seed in range(1, 5):
        tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.25), n_tips=14, age=4.0, seed=seed)
        g = z.simulate_genomes(tree, z.GenomeWiseRates(duplication=0.4, transfer=0.2, loss=0.3,
                                                       origination=0.7), initial_size=12, seed=seed)
        families = read_events_trace(events_trace_from_log(g.event_log))
        assert extant_species_from_records(families, g.species_tree) == g._gid_to_species()


def test_replay_from_trace_matches_live_within_precision():
    """A phylogram replayed from the trace matches the live one to the trace's time precision."""
    from zombi2.reconciliation import extant_species_from_records
    from zombi2.simulation import events_trace_from_log, read_events_trace

    tree, g = _sim()
    families = read_events_trace(events_trace_from_log(g.event_log))
    g2s = extant_species_from_records(families, g.species_tree)
    se = z.SequenceEvolution(branch_sigma=0.5, family_speed=z.LogNormal(0.0, 0.4))
    live = se.scale(g, seed=13)
    replay = se.scale_families(g.species_tree, families, g2s, seed=13)
    assert live.family_speed == replay.family_speed          # speeds are exact
    checked = 0
    for fam, tree_live in live.extant.items():
        if tree_live is None:
            continue
        checked += 1
        # same topology/labels; branch lengths agree to ~the trace's 10-sig-fig time precision
        assert _total_length(replay.extant[fam]) == pytest.approx(_total_length(tree_live), rel=1e-5)
    assert checked > 0
