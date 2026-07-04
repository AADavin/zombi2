"""The event-trace output — ``simulate_genomes(..., output="trace")`` / ``--output trace``.

A trace keeps the genealogy in its cheapest form (raw engine columns on the Rust path, a
prebuilt log on the Python path) and defers the per-event objects and gene trees until asked.
It must produce **exactly** the same profile, the same ``Events_trace.tsv``, and the same gene
trees as the full ``output="genomes"`` path — it is only a lazier route to them.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2 import GenomeTrace
from zombi2.simulation import EVENTS_TRACE_HEADER, events_trace_from_log

RATES = dict(duplication=0.2, transfer=0.12, loss=0.3, origination=1.0)


def _tree(n=200, seed=1):
    return z.simulate_species_tree(z.Yule(birth=1.0), n_tips=n, age=8.0, seed=seed)


def _profile_key(p):
    r, c, d = p.coo
    return (p.families, p.species, sorted(zip(r.tolist(), c.tolist(), d.tolist())))


def test_output_trace_returns_a_genome_trace():
    tr = z.simulate_genomes(_tree(), output="trace", seed=5, **RATES)
    assert isinstance(tr, GenomeTrace)


def test_invalid_output_value_rejected():
    with pytest.raises(ValueError, match="output must be"):
        z.simulate_genomes(_tree(), output="banana", seed=1, **RATES)


def test_trace_profile_matches_full_genomes():
    tree = _tree()
    g = z.simulate_genomes(tree, seed=9, **RATES)
    tr = z.simulate_genomes(tree, output="trace", seed=9, **RATES)
    assert _profile_key(g.profiles) == _profile_key(tr.profiles)


def test_trace_gene_trees_and_reconciliations_match_full():
    tree = _tree()
    g = z.simulate_genomes(tree, seed=11, **RATES)
    tr = z.simulate_genomes(tree, output="trace", seed=11, **RATES)

    gt_full, gt_trace = g.gene_trees(), tr.gene_trees()
    assert gt_full.keys() == gt_trace.keys()
    assert all(gt_full[f] == gt_trace[f] for f in gt_full)

    rc_full, rc_trace = g.reconciliations(), tr.reconciliations()
    assert all(rc_full[f] == rc_trace[f] for f in rc_full)


def test_events_trace_file_matches_the_log(tmp_path):
    tree = _tree()
    g = z.simulate_genomes(tree, seed=3, **RATES)
    tr = z.simulate_genomes(tree, output="trace", seed=3, **RATES)

    tr.write(tmp_path / "trace", include={"trace"})
    text = (tmp_path / "trace" / "Events_trace.tsv").read_text()

    # header + one row per event; the fast (column) writer equals the record-based writer
    assert text.splitlines()[0] == EVENTS_TRACE_HEADER
    assert text.count("\n") == len(g.event_log) + 1  # +1 for the header line
    assert text == events_trace_from_log(g.event_log)


def test_trace_write_default_and_heavy_parts(tmp_path):
    tree = _tree()
    tr = z.simulate_genomes(tree, output="trace", seed=4, **RATES)

    # default include is trace + profiles
    tr.write(tmp_path / "a")
    assert (tmp_path / "a" / "Events_trace.tsv").exists()
    assert (tmp_path / "a" / "Profiles.tsv").exists()
    assert (tmp_path / "a" / "species_tree.nwk").exists()

    # a heavy part (trees) is honoured by promoting to the full log on demand
    tr2 = z.simulate_genomes(tree, output="trace", seed=4, **RATES)
    tr2.write(tmp_path / "b", include={"trace", "trees"})
    assert (tmp_path / "b" / "Events_trace.tsv").exists()
    assert (tmp_path / "b" / "gene_trees").is_dir()


def test_trace_write_rejects_unknown_part(tmp_path):
    tr = z.simulate_genomes(_tree(), output="trace", seed=1, **RATES)
    with pytest.raises(ValueError, match="unknown write component"):
        tr.write(tmp_path / "x", include={"bogus"})


def test_python_engine_supports_trace():
    """A flexible model routes to the Python engine; output="trace" still yields a GenomeTrace
    (wrapping the log it already built) and writes the same trace file."""
    tree = _tree(n=60)
    tr = z.simulate_genomes(tree, rates=z.GenomeWiseRates(0.2, 0.0, 0.3, 1.0),
                            output="trace", seed=2)
    assert isinstance(tr, GenomeTrace)
    assert tr._columns is None and tr._event_log is not None  # Python path: prebuilt log
    # trace file writes fine off the prebuilt log
    txt = events_trace_from_log(tr.event_log)
    assert txt.splitlines()[0] == EVENTS_TRACE_HEADER


def test_genomes_write_supports_trace_part(tmp_path):
    """The 'trace' component is also available on a full Genomes (record-based writer)."""
    g = z.simulate_genomes(_tree(), seed=6, **RATES)
    g.write(tmp_path / "g", include={"trace"})
    assert (tmp_path / "g" / "Events_trace.tsv").exists()
