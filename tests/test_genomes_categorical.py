"""The `Categorical` law (issue #478): one category per family, drawn by weight and shared by every
rate that reads the law; each rate maps the category to its own multiplier, normalised to mean 1
under the weights. Runs that use today's draws are byte-identical to main: the cache that shares a
category changes nothing for them."""

import collections
import hashlib

import pytest

from zombi2.genomes import read_run, simulate_genomes_family, simulate_genomes_ordered
from zombi2.params import Categorical, Drift, LogNormal, PerCopy, Random
from zombi2.params.parse import parse_rate
from zombi2.species import simulate_species_tree


def _digest(r) -> str:
    # times rounded to nine decimals: a lognormal draw goes through exp(), whose last bit differs
    # between platforms' maths libraries, and that noise is not what this digest is about
    return hashlib.sha256(repr([(round(e.time, 9), e.kind, e.family, e.copy)
                                for e in r.edges]).encode()).hexdigest()[:12]


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.0, n_extant=20, seed=1).complete_tree


WEIGHTS = {"stable": 0.1, "mobile": 0.9}
TABLES = {"duplication": {"stable": 0.2, "mobile": 1.1}, "transfer": {"stable": 0.05, "mobile": 1.1},
          "loss": {"stable": 0.05, "mobile": 1.1}}


def _rates():
    kind = Random("families", Categorical(WEIGHTS))
    base = {"duplication": 0.05, "transfer": 0.15, "loss": 0.5}
    return {k: PerCopy(base[k]).varying_among(kind, TABLES[k]) for k in base}


def _expected(key, category):
    mean = sum(WEIGHTS[c] * TABLES[key][c] for c in WEIGHTS)
    return TABLES[key][category] / mean


@pytest.mark.parametrize("label", ["family", "parallel", "ordered"])
def test_one_category_per_family_read_by_every_rate(tree, label):
    fn = simulate_genomes_ordered if label == "ordered" else simulate_genomes_family
    extra = {"parallel": True} if label == "parallel" else {}
    r = fn(tree, origination=1.0, initial_families=300, seed=2, **_rates(), **extra)
    cats = r.family_categories
    assert len(cats) >= 300
    share = collections.Counter(cats.values())["stable"] / len(cats)
    se = (0.1 * 0.9 / len(cats)) ** 0.5
    assert abs(share - 0.1) < 4 * se, share
    # every rate's multiplier is the table's number for the family's category, over the weighted mean
    for fam, cat in cats.items():
        for key in ("duplication", "transfer", "loss"):
            assert r.family_multipliers[fam][key] == pytest.approx(_expected(key, cat))


def test_the_category_is_written_read_back_and_streamed(tree, tmp_path):
    r = simulate_genomes_family(tree, origination=1.0, initial_families=50, seed=2, **_rates())
    r.write(tmp_path / "mem")
    lines = (tmp_path / "mem" / "family_multipliers.tsv").read_text().splitlines()
    assert lines[0] == "family\tduplication\ttransfer\tloss\tcategory"
    assert all(ln.split("\t")[-1] in WEIGHTS for ln in lines[1:])
    back = read_run(tmp_path / "mem")
    assert back.family_categories == r.family_categories
    s = simulate_genomes_family(tree, origination=1.0, initial_families=50, seed=2, stream_to=tmp_path / "s",
                                **_rates())
    streamed = (tmp_path / "s" / "family_multipliers.tsv").read_text().splitlines()
    assert streamed[0] == lines[0] and all(ln.split("\t")[-1] in WEIGHTS for ln in streamed[1:])
    assert s.n_families >= 50
    # without the law, no column
    r0 = simulate_genomes_family(tree, origination=1.0, initial_families=10, seed=2,
                                 loss=PerCopy(0.5).varying_among("families", LogNormal(0.0, 0.5)))
    r0.write(tmp_path / "plain")
    assert (tmp_path / "plain" / "family_multipliers.tsv").read_text().splitlines()[0] == \
        "family\tduplication\ttransfer\tloss"
    assert r0.family_categories == {}


def test_the_text_form_parses_for_one_rate():
    r = parse_rate("PerCopy(0.5).varying_among(Random('families', Categorical({'stable': 0.1, "
                   "'mobile': 0.9})), {'stable': 0.05, 'mobile': 1.1})")
    m = r.modifiers[0]
    assert m.factor_for("stable") == pytest.approx(0.05 / (0.1 * 0.05 + 0.9 * 1.1))


def test_a_categorical_law_is_checked():
    with pytest.raises(ValueError, match="two or more"):
        Categorical({"only": 1.0})
    with pytest.raises(ValueError, match="positive"):
        Categorical({"a": 1.0, "b": -1.0})
    kind = Random("families", Categorical(WEIGHTS))
    with pytest.raises(TypeError, match="multiplier per category"):
        PerCopy(0.5).varying_among(kind)
    with pytest.raises(TypeError, match="name the draw"):
        PerCopy(0.5).varying_among("families", Categorical(WEIGHTS))
    with pytest.raises(ValueError, match="exactly the law's categories"):
        PerCopy(0.5).varying_among(kind, {"stable": 1.0})
    with pytest.raises(ValueError, match="no meaning"):
        Random("families", Drift(Categorical(WEIGHTS)))


def test_two_laws_in_one_run_are_refused(tree):
    a = Random("families", Categorical(WEIGHTS))
    b = Random("families", Categorical({"x": 0.5, "y": 0.5}))
    with pytest.raises(ValueError, match="two Categorical laws"):
        simulate_genomes_family(tree, initial_families=5, seed=1,
                                duplication=PerCopy(0.1).varying_among(a, TABLES["duplication"]),
                                loss=PerCopy(0.5).varying_among(b, {"x": 0.5, "y": 1.5}))


#: digests of runs that use today's draws, taken on an export of main (9d37d4f0): the shared cache
#: the category law needs must not move any of them
_BEFORE = {"family shared": "a4cf3e70028f", "parallel shared": "dbd099d54791",
           "family lineages": "3710a80059b3", "ordered shared": "2464fec2b300"}


def test_runs_with_the_draws_of_before_are_byte_identical():
    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1).complete_tree
    speed = Random("families", LogNormal(0.0, 0.6))
    shared = dict(duplication=PerCopy(0.2).varying_among(speed), loss=PerCopy(0.3).varying_among(speed),
                  transfer=PerCopy(0.2).varying_among("families", LogNormal(0.0, 0.4)), origination=0.5,
                  initial_families=40, seed=3)
    assert _digest(simulate_genomes_family(tree, **shared)) == _BEFORE["family shared"]
    assert _digest(simulate_genomes_family(tree, parallel=True, **shared)) == _BEFORE["parallel shared"]
    lineages = dict(duplication=PerCopy(0.2).varying_among("lineages", Drift(LogNormal(0.0, 0.3))),
                    loss=PerCopy(0.3).varying_among("lineages", LogNormal(0.0, 0.3)), origination=0.5,
                    initial_families=40, seed=3)
    assert _digest(simulate_genomes_family(tree, **lineages)) == _BEFORE["family lineages"]
    assert _digest(simulate_genomes_ordered(tree, **shared)) == _BEFORE["ordered shared"]
