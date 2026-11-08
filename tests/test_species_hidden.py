"""Hidden lineages for a fixed sampled tree (issue #480): `add_hidden_lineages` and
`fit_birth_death_rates`. The maths is in the module's docstring; here the guarantees: pruning the
result gives the input back exactly, the degenerate case adds nothing (the attachments and the
unsampled tips against their closed forms are a validation check, in ``test_validation.py``), and the fit recovers the rates it was given."""

import re
import statistics

import pytest

from zombi2.species import add_hidden_lineages, fit_birth_death_rates, simulate_species_tree
from zombi2.tree import read_newick


def _shape(tree, labels=None) -> str:
    """Topology and branch lengths, labels dropped: what pruning must give back."""
    def emit(i):
        n = tree.nodes[i]
        bl = round(n.end_time - n.birth_time, 9)
        return f"{bl}" if not n.children else "(" + ",".join(emit(c) for c in n.children) + f"):{bl}"
    return emit(tree.root)


@pytest.fixture(scope="module")
def sampled():
    full = simulate_species_tree(birth=1.0, death=0.4, n_extant=60, sampling=0.5, seed=1)
    text = full.extant_tree.to_newick()
    labels = re.findall(r"n\d+", text)
    relabelled = text
    for k, lab in enumerate(labels):                      # a user's own tip labels, not ZOMBI2's
        relabelled = relabelled.replace(lab + ":", f"sp{k}:", 1)
    return read_newick(relabelled)


def test_pruning_the_result_recovers_the_input_exactly(sampled):
    tree, names = sampled
    r = add_hidden_lineages(tree, birth=1.0, death=0.4, sampling=0.5, seed=1, names=names)
    back = r.extant_tree
    assert _shape(back) == _shape(tree)                   # topology, lengths and children order
    # the labels travel: every sampled tip of the result maps to the input's label of the same tip
    tips_in = {i: names[i] for i in tree.leaves()}
    assert sorted(r.names[i] for i in back.leaves()) == sorted(tips_in.values())
    assert set(n.fate for n in r.complete_tree.nodes.values() if not n.children) >= {"extant"}
    assert r.complete_tree.nodes[r.complete_tree.root].birth_time == 0.0


def test_a_stem_adds_time_above_the_crown_and_attachments_on_it(sampled):
    tree, _ = sampled
    r = add_hidden_lineages(tree, birth=1.0, death=0.4, sampling=0.5, stem=2.0, seed=2)
    root = r.complete_tree.nodes[r.complete_tree.root]
    crown = r.extant_tree.nodes[r.extant_tree.root]
    assert crown.end_time - crown.birth_time == pytest.approx(2.0)      # the stem, as pruned
    # a hidden lineage born on the stem starts before the crown
    assert any(n.birth_time < crown.end_time and n.fate in ("extinct", "unsampled")
               for n in r.complete_tree.nodes.values() if not n.children)
    assert root.birth_time == 0.0


def test_no_death_and_full_sampling_add_nothing():
    cherry, _ = read_newick("(a:1.0,b:1.0);")
    r = add_hidden_lineages(cherry, birth=1.0, death=0.0, sampling=1.0, seed=3)
    assert len(r.complete_tree.nodes) == 3
    assert r.complete_tree.unsampled_leaves() == [] and r.complete_tree.extinct_leaves() == []


def test_full_sampling_leaves_only_extinct_hidden_lineages():
    cherry, _ = read_newick("(a:1.0,b:1.0);")
    seen_hidden = False
    for s in range(20):
        r = add_hidden_lineages(cherry, birth=1.5, death=0.8, sampling=1.0, seed=s)
        assert r.complete_tree.unsampled_leaves() == []
        seen_hidden |= bool(r.complete_tree.extinct_leaves())
    assert seen_hidden


def test_the_input_is_checked():
    cherry, _ = read_newick("(a:1.0,b:1.0);")
    with pytest.raises(ValueError, match="positive"):
        add_hidden_lineages(cherry, birth=0.0, death=0.0, sampling=1.0)
    with pytest.raises(ValueError, match="sampling"):
        add_hidden_lineages(cherry, birth=1.0, sampling=0.0)
    skewed, _ = read_newick("(a:1.0,b:0.5);", tip_fates={"a": "extant", "b": "extinct"})
    with pytest.raises(ValueError, match="extant"):
        add_hidden_lineages(skewed, birth=1.0, death=0.2, sampling=1.0)
    with pytest.raises(ValueError, match="stem"):
        add_hidden_lineages(cherry, birth=1.0, stem=-1.0)


def test_the_fit_recovers_the_rates_it_was_given():
    births, both = [], []
    for s in range(20):
        t = simulate_species_tree(birth=1.0, death=0.4, n_extant=200, sampling=0.5, seed=100 + s).extant_tree
        births.append(fit_birth_death_rates(t, sampling=0.5, death=0.4)[0])
        both.append(fit_birth_death_rates(t, sampling=0.5))
    assert abs(statistics.median(births) - 1.0) < 0.15, statistics.median(births)
    assert abs(statistics.median(b for b, _ in both) - 1.0) < 0.2
    assert abs(statistics.median(d for _, d in both) - 0.4) < 0.2


def test_the_result_writes_the_species_files_and_the_names(sampled, tmp_path):
    tree, names = sampled
    r = add_hidden_lineages(tree, birth=1.0, death=0.4, sampling=0.5, seed=1, names=names)
    r.write(tmp_path)
    written = sorted(p.name for p in tmp_path.iterdir())
    assert {"species_complete.nwk", "species_extant.nwk", "species_fates.tsv",
            "species_events.tsv", "names.tsv"} <= set(written)
    fates = [ln.split("\t")[1] for ln in (tmp_path / "species_fates.tsv").read_text().splitlines()[1:]]
    assert set(fates) >= {"extant"} and set(fates) <= {"extant", "extinct", "unsampled"}
    names_rows = (tmp_path / "names.tsv").read_text().splitlines()
    assert names_rows[0] == "node\tname" and len(names_rows) - 1 == len(names)   # tips and labelled splits
