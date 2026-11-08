"""The published profile mixtures C10 to C60 (issue #479): `SubstitutionModel.profile_mixture`.

Every site of every family draws one profile by weight and evolves under the model's
exchangeabilities over that profile's frequencies. The draw follows the seed on both engines; the
profile each site drew is kept (`profile_assignments`) and written (`profile_assignments.tsv`). That
the sites of a profile really carry its frequencies is a validation check, in ``test_validation.py``;
here the plumbing: composition with `across_sites`, both engines, the table, recording, refusals."""

import collections
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences
from zombi2.sequences._profile_mixtures import MIXTURES
from zombi2.sequences.lineage_models import Models
from zombi2.sequences.substitution_models import (PROFILE_MIXTURES, AMINO_ACIDS, hky85, lg, poisson,
                                                  wag)
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def genomes():
    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1).complete_tree
    return simulate_genomes_family(tree, duplication=0.1, loss=0.2, origination=0.4,
                                   initial_families=12, seed=1)


def test_the_six_published_mixtures_are_whole():
    assert PROFILE_MIXTURES == tuple(MIXTURES) == ("C10", "C20", "C30", "C40", "C50", "C60")
    for label, (weights, rows) in MIXTURES.items():
        k = int(label[1:])
        assert len(weights) == len(rows) == k
        assert sum(weights) == pytest.approx(1.0, abs=1e-8)
        assert all(len(r) == 20 and sum(r) == pytest.approx(1.0, abs=1e-8) and min(r) > 0 for r in rows)
    # one number against the source: the weight of C10's first profile, as IQ-TREE's table has it
    assert MIXTURES["C10"][0][0] == pytest.approx(0.1191344178)


def test_the_mixture_decorates_a_protein_model_and_composes_with_rate_classes():
    m = lg().profile_mixture("C60")
    assert m.name == "LG+C60" and len(m.profile_models) == 60 and len(m.profile_weights) == 60
    assert sum(m.profile_weights) == pytest.approx(1.0)
    # each profile model: LG's exchangeabilities over that profile's frequencies, normalised
    pi12 = np.asarray(MIXTURES["C60"][1][11])
    assert np.allclose(m.profile_models[11].stationary, pi12 / pi12.sum())
    assert m.profile_models[11].name == "LG+C60[12]"
    S_lg = lg().Q / lg().stationary[None, :]
    S_12 = m.profile_models[11].Q / m.profile_models[11].stationary[None, :]
    off = ~np.eye(20, dtype=bool)
    ratio = S_12[off] / S_lg[off]
    assert np.allclose(ratio, ratio[0])            # the same chemistry, up to the normalisation
    # either order of the two decorations gives every profile model the rate classes
    g1 = m.across_sites(gamma_shape=0.5)
    g2 = lg().across_sites(gamma_shape=0.5).profile_mixture("c60")
    assert g1.name == "LG+C60+G4" and g2.name == "LG+G4+C60"
    for g in (g1, g2):
        assert all(p.site_rates == g.site_rates and p.site_shares == g.site_shares
                   for p in g.profile_models)
    assert wag().profile_mixture("C10").name == "WAG+C10"


def test_a_mixture_is_checked():
    with pytest.raises(ValueError, match="amino-acid profiles"):
        hky85().profile_mixture("C60")
    with pytest.raises(ValueError, match="unknown profile mixture"):
        lg().profile_mixture("C70")
    with pytest.raises(ValueError, match="already carries"):
        lg().profile_mixture("C10").profile_mixture("C20")
    from dataclasses import replace
    with pytest.raises(ValueError, match="same length"):
        replace(lg(), profile_models=(poisson(),), profile_weights=())
    with pytest.raises(ValueError, match="same alphabet"):
        replace(lg(), profile_models=(hky85(),), profile_weights=(1.0,))


@pytest.mark.parametrize("parallel", [False, 2])
def test_every_site_of_every_family_draws_a_profile_by_the_seed(genomes, parallel):
    model = lg().profile_mixture("C20").across_sites(gamma_shape=0.8)
    kw = dict(model=model, length=150, seed=3, parallel=parallel)
    r = simulate_sequences(genomes, **kw)
    again = simulate_sequences(genomes, **kw)
    assert r.alignments == again.alignments and r.profile_assignments == again.profile_assignments
    assert set(r.profile_assignments) == set(r.alignments)
    for fam, drawn in r.profile_assignments.items():
        assert len(drawn) == 150 and 1 <= min(drawn) and max(drawn) <= 20
        for seq in list(r.alignments[fam].values()) + list(r.ancestral[fam].values()):
            assert len(seq) == 150 and set(seq) <= set(AMINO_ACIDS)
    pooled = collections.Counter(p for drawn in r.profile_assignments.values() for p in drawn)
    assert len(pooled) >= 15                 # with 150 sites × 12+ families, nearly every profile shows
    if parallel:
        more = simulate_sequences(genomes, **{**kw, "parallel": 4})
        assert more.alignments == r.alignments and more.profile_assignments == r.profile_assignments


def test_a_plain_model_draws_no_profiles(genomes):
    r = simulate_sequences(genomes, model=lg(), length=50, seed=3)
    assert r.profile_assignments == {}


def test_the_assignments_are_written_and_streamed_alike(genomes, tmp_path):
    model = lg().profile_mixture("C10")
    r = simulate_sequences(genomes, model=model, length=40, seed=2)
    r.write(tmp_path / "mem")
    lines = (tmp_path / "mem" / "profile_assignments.tsv").read_text().splitlines()
    assert lines[0] == "family\tsite\tprofile_number"
    assert len(lines) - 1 == 40 * len(r.profile_assignments)
    fam0 = sorted(r.profile_assignments)[0]
    assert lines[1] == f"{fam0}\t0\t{r.profile_assignments[fam0][0]}"
    # a streamed run writes the same table as the in-memory run of the same engine (the parallel
    # engine is its own realisation for a seed, so it is compared with its own in-memory run)
    for engine in ({}, {"parallel": 2}):
        ref = r if not engine else simulate_sequences(genomes, model=model, length=40, seed=2, **engine)
        ref.write(tmp_path / ("ref" + ("_p" if engine else "")))
        expected = (tmp_path / ("ref" + ("_p" if engine else "")) / "profile_assignments.tsv").read_text()
        out = tmp_path / ("stream" + ("_p" if engine else ""))
        simulate_sequences(genomes, model=model, length=40, seed=2, stream_to=out, **engine)
        assert (out / "profile_assignments.tsv").read_text() == expected
        for fam, aln in ref.alignments.items():
            if aln:
                text = (out / "alignments" / f"fam{fam}.fasta").read_text().replace("\n", "")
                assert all(seq in text for seq in aln.values())
    # the token: left out, no file; a plain model writes none either way
    r.write(tmp_path / "without", outputs=("alignments", "summary"))
    assert not (tmp_path / "without" / "profile_assignments.tsv").exists()
    simulate_sequences(genomes, model=lg(), length=40, seed=2).write(tmp_path / "plain")
    assert not (tmp_path / "plain" / "profile_assignments.tsv").exists()


def test_a_recorded_run_names_the_sites_the_user_sees():
    """The sites of one profile evolve as one block and the columns go back in site order
    afterwards; the substitution log has to name the column the user sees, not the block's."""
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1).complete_tree
    g = simulate_genomes_family(tree, duplication=0.0, loss=0.0, origination=0.0,
                                initial_families=1, seed=1)
    r = simulate_sequences(g, model=lg().profile_mixture("C20"), length=120, seed=5, record=True)
    seqs = {**r.alignments[0], **r.ancestral[0]}
    last: dict = {}
    for e in sorted(r.events, key=lambda e: e.time):
        last[(f"{e.lineage}_g{e.gene}", e.site)] = e.to_state
    assert last and all(seqs[label][site] == state for (label, site), state in last.items())
    touched = {site for _label, site in last}
    founding = r.founding[0]
    assert not [i for seq in r.alignments[0].values() for i, ch in enumerate(seq)
                if ch != founding[i] and i not in touched]


def test_what_a_mixture_is_refused_with(genomes):
    model = lg().profile_mixture("C10")
    with pytest.raises(ValueError, match="two answers"):
        simulate_sequences(genomes, model=model, length=20, seed=1,
                           profiles={0: np.full((20, 20), 0.05)})
    with pytest.raises(ValueError, match="indels are refused"):
        simulate_sequences(genomes, model=model, length=20, seed=1, deletion=0.05)
    with pytest.raises(ValueError, match="partitions cannot hold"):
        simulate_sequences(genomes, partitions=[(model, 10), (lg(), 10)], seed=1)
    clades = Models(models={"a": model, "rest": lg()})      # a set_by, before it is painted
    with pytest.raises(ValueError, match="per-clade model set"):
        simulate_sequences(genomes, model=clades, length=20, seed=1)


def test_the_command_line_takes_the_mixture(tmp_path):
    def zombi2(*args):
        return subprocess.run([sys.executable, "-m", "zombi2.cli.main", *args],
                              capture_output=True, text=True)
    out = tmp_path / "run"
    assert zombi2("species", str(out), "--birth", "1.0", "--n-extant", "6", "--seed", "1").returncode == 0
    assert zombi2("genomes", str(out), "--duplication", "0.1", "--loss", "0.2", "--origination", "0.3",
                  "--initial-families", "5", "--seed", "1").returncode == 0
    done = zombi2("sequences", str(out), "--model", "lg", "--mixture", "c20", "--gamma-shape", "0.7",
                  "--length", "30", "--seed", "2")
    assert done.returncode == 0, done.stderr
    assert "LG+C20+G4" in done.stdout
    rows = (out / "sequences" / "profile_assignments.tsv").read_text().splitlines()
    assert rows[0] == "family\tsite\tprofile_number" and len(rows) > 1
    log = (out / "sequences" / "sequences.log").read_text()
    assert "mixture\tC20" in log                    # the log reproduces the run
    bad = zombi2("sequences", str(out), "--model", "hky85", "--mixture", "C10", "--length", "30",
                 "--seed", "2", "--force")
    assert bad.returncode != 0 and "needs a protein model" in bad.stderr
    typo = zombi2("sequences", str(out), "--model", "lg", "--mixture", "C70", "--length", "30",
                  "--seed", "2", "--force")
    assert typo.returncode != 0 and "invalid choice" in typo.stderr


def test_a_run_without_a_mixture_is_unchanged(genomes):
    """The byte-identity promise: nothing here touches a run whose model has no mixture, so its
    draws are the ones of before. Pinned on an export of main (fe005543)."""
    import hashlib
    r = simulate_sequences(genomes, model=lg().across_sites(gamma_shape=0.7), length=60, seed=11)
    text = "".join(f"{f}:{k}:{v}\n" for f in sorted(r.alignments)
                   for k, v in sorted(r.alignments[f].items()))
    assert hashlib.sha256(text.encode()).hexdigest()[:12] == _BEFORE_PLAIN


#: taken on an export of main (fe005543); the branch gives the same bytes
_BEFORE_PLAIN = "fdb35d0dcca7"
