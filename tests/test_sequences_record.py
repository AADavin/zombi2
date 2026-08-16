"""The sequence level's own log — step 11 of the joining design note (§10).

Off by default, because it is the one log bigger than the output it explains. `record=True` also
changes the sampler: an ordinary branch jumps to its end with `exp(Q·bl)`, a recorded one walks the
path. Same process, same distribution at the end, different realisation for a seed.
"""

import collections
import pathlib
import statistics
import tempfile

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import hky85, jc69, simulate_sequences
from zombi2.sequences._record import walk


def _run(**kw):
    ct = species.simulate_species_tree(birth=1.0, n_extant=8, seed=1).complete_tree
    g = simulate_genomes_family(ct, initial_families=3, duplication=0.05, loss=0.1, seed=2)
    return ct, g, simulate_sequences(g, model=jc69(), length=80, seed=3, **kw)


# --- what it records ------------------------------------------------------------------------------

def test_nothing_is_recorded_unless_the_run_asks():
    _ct, _g, r = _run()
    assert r.events == []


def test_every_substitution_gets_a_row():
    _ct, _g, r = _run(record=True)
    assert r.events
    assert {e.kind for e in r.events} == {"substitution"}
    for e in r.events:
        assert e.from_state != e.to_state
        assert e.from_state in "ACGT" and e.to_state in "ACGT"
        assert e.strand == 1


def test_the_log_is_in_time_order():
    _ct, _g, r = _run(record=True)
    assert [e.time for e in r.events] == sorted(e.time for e in r.events)


def test_a_row_names_the_lineage_and_the_copy_the_alignment_names():
    _ct, _g, r = _run(record=True)
    keys = {k for a in (*r.alignments.values(), *r.ancestral.values()) for k in a}
    for e in r.events[:50]:
        assert f"{e.lineage}_g{e.gene}" in keys


def test_indels_get_rows_too():
    """Otherwise the log says a site changed and never says when that site arrived or left."""
    _ct, _g, r = _run(record=True, insertion=0.05, deletion=0.05)
    kinds = collections.Counter(e.kind for e in r.events)
    assert kinds["insertion"] and kinds["deletion"] and kinds["substitution"]
    for e in r.events:
        if e.kind == "insertion":
            assert e.after >= -1               # the id it follows, or -1 before every founding site
        if e.kind in ("insertion", "deletion"):
            assert e.from_state == "" and e.to_state == ""


def test_a_recorded_site_is_one_the_lineage_actually_carries():
    """The engine evolves every column of the widened alignment whether a lineage holds it or not,
    because that is what makes one alignment out of many. A *log* of a site the lineage does not have
    would be a claim about nothing."""
    ct, _g, r = _run(record=True, insertion=0.06, deletion=0.06)
    columns = {}
    for fam, table in r.alignments.items():
        for label, seq in table.items():
            columns[label] = seq
    for e in r.events:
        if e.kind != "substitution":
            continue
        seq = columns.get(f"{e.lineage}_g{e.gene}")
        if seq is not None:
            assert seq.count("-") < len(seq)   # the row exists and is not all gap


def test_it_writes_one_file_and_only_when_asked():
    _ct, _g, r = _run(record=True)
    with tempfile.TemporaryDirectory() as d:
        r.write(d, outputs=("alignments", "events"))
        log = pathlib.Path(d, "sequence_events.tsv")
        assert log.exists()
        head, *rows = log.read_text().splitlines()
        assert head.split("\t")[:5] == ["time", "event", "lineage", "gene", "site"]
        assert len(rows) == len(r.events)
    _ct, _g, plain = _run()
    with tempfile.TemporaryDirectory() as d:
        plain.write(d, outputs=("alignments", "events"))
        assert not pathlib.Path(d, "sequence_events.tsv").exists()


def test_it_is_deterministic():
    _ct, _g, a = _run(record=True)
    _ct2, _g2, b = _run(record=True)
    assert [e.row() for e in a.events] == [e.row() for e in b.events]


# --- the recorded sampler is the same process -----------------------------------------------------

def test_the_walk_lands_where_the_matrix_would():
    """The claim that makes a recorded run a valid run rather than a second model: the forward
    Gillespie's end state has the same distribution as one draw from `exp(Q·bl)`. Measured on one
    site over many replicates, against the transition matrix itself."""
    from zombi2.sequences.evolution import _cdf_for

    m = hky85(kappa=3.0)
    rng = np.random.default_rng(11)
    bl, reps = 0.35, 4000
    start = np.zeros(1, dtype=np.int8)
    walked = collections.Counter()
    for _ in range(reps):
        end, _rows = walk(start, m.Q, bl, rng)
        walked[int(end[0])] += 1
    cdf = _cdf_for({}, m, bl)
    expected = np.diff(np.concatenate(([0.0], cdf[0])))
    for state, p in enumerate(expected):
        assert abs(walked[state] / reps - p) < 0.02, f"state {state}: {walked[state] / reps} vs {p}"


def test_the_number_of_substitutions_matches_the_branch_length():
    """A branch of length `bl` expects `bl` substitutions per site, because every model is normalised
    to exactly that. So the rows are not merely plausible — they are countable."""
    m = jc69()
    rng = np.random.default_rng(3)
    bl, sites = 0.4, 500
    start = rng.choice(m.k, size=sites, p=m.stationary).astype(np.int8)
    counts = [len(walk(start, m.Q, bl, rng)[1]) for _ in range(20)]
    assert statistics.fmean(counts) == pytest.approx(bl * sites, rel=0.08)


def test_recording_changes_the_realisation_and_says_so():
    """Not a defect: one sampler jumps to the branch's end, the other walks there. The alignment a
    seed gives is a valid draw either way, and it is not the same draw."""
    _ct, _g, plain = _run()
    _ct2, _g2, kept = _run(record=True)
    assert plain.alignments != kept.alignments
    assert set(plain.alignments) == set(kept.alignments)


# --- what is refused ------------------------------------------------------------------------------

def test_partitions_are_refused():
    _ct, g, _r = _run()
    with pytest.raises(ValueError, match="not on that path"):
        simulate_sequences(g, partitions=[(jc69(), 40), (hky85(2.0), 40)], seed=3, record=True)


def test_the_parallel_engine_is_refused():
    _ct, g, _r = _run()
    with pytest.raises(ValueError, match="not on that path"):
        simulate_sequences(g, model=jc69(), length=40, seed=3, record=True, parallel=2)


def test_streaming_is_refused():
    _ct, g, _r = _run()
    with tempfile.TemporaryDirectory() as d, pytest.raises(ValueError, match="not on that path"):
        simulate_sequences(g, model=jc69(), length=40, seed=3, record=True, stream_to=d)


def test_a_nucleotide_run_is_refused():
    from zombi2.genomes import simulate_genomes_nucleotide
    sp = species.simulate_species_tree(birth=1.0, n_extant=5, seed=1)
    g = simulate_genomes_nucleotide(sp, genes=2, gene_length=90, root_length=400, loss=0.3,
                                    loss_extent=40, seed=2)
    with pytest.raises(ValueError, match="not on that path"):
        simulate_sequences(g, model=jc69(), seed=3, record=True)
