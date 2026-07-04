"""The event-trace output — ``simulate_genomes(..., output="trace")`` / ``--output trace``.

A trace keeps the genealogy in its cheapest form (raw engine columns on the Rust path, a
prebuilt log on the Python path). The Rust ``output="trace"`` engine drops speciation records
and keeps gene ids across speciations, so the on-disk trace is compact (O/D/T/L only); the full
genealogy is recovered by replaying it against the species tree. It must reconstruct the *same*
species-level structure as the full ``output="genomes"`` path (gene ids differ between the two
minting schemes, so comparisons are gid-independent).
"""

import re

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


def _recon_signature(genlike):
    """Species-level reconciliation signature (gid-independent): the multiset of
    (family, event, species, recipient, time) over every reconciled event."""
    from collections import Counter
    sig = Counter()
    for fam, r in genlike.reconciliations().items():
        for e in r.events:
            sig[(fam, e.event, e.species, e.recipient, round(e.time, 6))] += 1
    return sig


def _extant_leafsets(genlike):
    """Per-family multiset of extant-tree leaf *species* (gid-independent)."""
    from collections import Counter
    out = {}
    for fam, (_complete, extant) in genlike.gene_trees().items():
        if extant:
            leaves = re.findall(r"[(,]([^(),:;]+):", extant) or \
                     [re.match(r"([^(),:;]+):", extant).group(1)]
            out[fam] = Counter(lbl.rsplit("_", 1)[0] for lbl in leaves)
    return out


def test_trace_reconstruction_matches_full_log():
    """The no-remint trace engine sees the same D/T/L/O draws as the full-log engine (same seed),
    so replaying it must reconstruct the *same* genealogy — identical reconciliation events and
    extant-tree species leaves. Gene ids differ (the two schemes mint them differently), so we
    compare gid-independent, species-level structure."""
    tree = _tree()
    g = z.simulate_genomes(tree, seed=11, **RATES)
    tr = z.simulate_genomes(tree, output="trace", seed=11, **RATES)

    assert _recon_signature(g) == _recon_signature(tr)
    assert _extant_leafsets(g) == _extant_leafsets(tr)


def test_events_trace_file_is_compact(tmp_path):
    """The trace file carries only O/D/T/L rows — no speciation ``S`` rows (a lineage keeps its
    id across speciations), so it is far smaller than the full log."""
    tree = _tree()
    g = z.simulate_genomes(tree, seed=3, **RATES)
    tr = z.simulate_genomes(tree, output="trace", seed=3, **RATES)

    tr.write(tmp_path / "trace", include={"trace"})
    lines = (tmp_path / "trace" / "Events_trace.tsv").read_text().splitlines()
    assert lines[0] == EVENTS_TRACE_HEADER
    events = [ln.split("\t")[1] for ln in lines[1:] if ln.strip()]
    assert "S" not in events                       # no speciation rows
    assert set(events) <= {"O", "D", "T", "L"}
    assert len(events) < len(g.event_log)          # strictly smaller than the full log


def test_events_trace_file_roundtrips(tmp_path):
    """A written compact trace, read back with the species tree, expands to a full genealogy that
    reconstructs the same species-level structure as the in-memory trace."""
    from zombi2 import read_events_trace
    from zombi2.reconciliation import extant_species_from_records, reconcile
    from collections import Counter

    tree = _tree()
    tr = z.simulate_genomes(tree, output="trace", seed=7, **RATES)
    tr.write(tmp_path / "trace", include={"trace"})
    text = (tmp_path / "trace" / "Events_trace.tsv").read_text()

    families = read_events_trace(text, tree)               # expands (compact → full)
    g2s = extant_species_from_records(families, tree)
    total_age = tree.total_age

    # compare gid- and time-independently: on-disk times are 10-sig-fig, so a 6th-decimal
    # rounding boundary would spuriously differ. Event type + species branch + recipient is
    # what a downstream reconciliation consumer cares about.
    def sig_no_time(events):
        return Counter((fam, e.event, e.species, e.recipient) for fam, es in events for e in es)

    from_disk = sig_no_time((fam, reconcile(recs, g2s, total_age).events)
                            for fam, recs in families.items())
    in_memory = sig_no_time((fam, r.events) for fam, r in tr.reconciliations().items())
    assert from_disk == in_memory


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
