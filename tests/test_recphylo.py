"""Tests for recPhyloXML output (zombi2.tools.recphylo)."""

import collections
import xml.etree.ElementTree as ET

import pytest

from zombi2.genomes import simulate_genomes_family
from zombi2.species import simulate_species_tree
from zombi2.tools.recphylo import recphylo_xml, write_recphylo


@pytest.fixture(scope="module")
def run():
    tree = simulate_species_tree(birth=1.0, death=0.4, n_extant=10, seed=3)
    g = simulate_genomes_family(tree, duplication=0.3, transfer=0.3, loss=0.3, origination=0.3,
                                initial_families=8, seed=5)
    return g


def _clades(elem):
    """Every ``<clade>`` under ``elem``, in document order."""
    return elem.findall(".//clade")


def _terminal(clade):
    """The tag that ends this clade's ``<eventsRec>`` — the event the gene-tree node *is*."""
    return list(clade.find("eventsRec"))[-1].tag


def test_the_document_is_well_formed_and_has_the_expected_skeleton(run):
    root = ET.fromstring(recphylo_xml(run.gene_trees, run.complete_tree))
    assert root.tag == "recPhylo"
    assert [c.tag for c in root][:1] == ["spTree"]                  # the species tree comes first
    assert len(root.findall("recGeneTree")) == len(run.gene_trees)  # then one block per family
    for rgt in root.findall("recGeneTree"):
        phy = rgt.find("phylogeny")
        assert phy.get("rooted") == "true"
        assert phy.find("name").text.startswith("fam")


def test_the_species_tree_holds_every_node_a_gene_points_at(run):
    root = ET.fromstring(recphylo_xml(run.gene_trees, run.complete_tree))
    named = {c.find("name").text for c in _clades(root.find("spTree"))}
    # complete: extinct species too — and they are named e<id>, which is how the file says so
    assert named == set(run.complete_tree.labels().values())
    assert any(lab[0] == "e" for lab in named), "the run should have lost a lineage"
    pointed_at = set()
    for rgt in root.findall("recGeneTree"):
        for clade in _clades(rgt):
            for ev in clade.find("eventsRec"):
                pointed_at.add(ev.get("speciesLocation") or ev.get("destinationSpecies"))
    assert pointed_at <= named                                     # no dangling speciesLocation


def test_every_gene_tree_node_becomes_one_clade_with_the_right_event(run):
    # the mapping is the whole point: one clade per node, and its terminal tag is that node's event
    tag_of = {"duplication": "duplication", "speciation": "speciation", "transfer": "branchingOut",
              "loss": "loss", "extant": "leaf", "extinct": "leaf", "unsampled": "leaf"}
    for fam, gt in run.gene_trees.items():
        nodes, stack = [], [gt.complete]
        while stack:
            n = stack.pop()
            nodes.append(n)
            stack.extend(n.children)
        rgt = ET.fromstring(recphylo_xml({fam: gt}, run.complete_tree)).find("recGeneTree")
        clades = _clades(rgt)
        assert len(clades) == len(nodes)
        assert (collections.Counter(_terminal(c) for c in clades)
                == collections.Counter(tag_of[n.kind] for n in nodes))


def test_a_transfer_is_the_format_s_two_step(run):
    # branchingOut on the node the copy left from; the child that arrived opens with transferBack,
    # naming where it landed — and the other child stays put, with no transferBack at all
    seen = 0
    for fam, gt in run.gene_trees.items():
        rgt = ET.fromstring(recphylo_xml({fam: gt}, run.complete_tree)).find("recGeneTree")
        for clade in _clades(rgt):
            events = list(clade.find("eventsRec"))
            if events[-1].tag != "branchingOut":
                assert [e.tag for e in events[:-1]] in ([], ["transferBack"])
                continue
            donor = events[-1].get("speciesLocation")
            kids = clade.findall("clade")
            assert len(kids) == 2
            backs = [k for k in kids if list(k.find("eventsRec"))[0].tag == "transferBack"]
            assert len(backs) == 1                              # exactly one child arrived
            stayed = next(k for k in kids if k not in backs)
            assert list(stayed.find("eventsRec"))[0].tag != "transferBack"
            assert list(stayed.find("eventsRec"))[-1].get("speciesLocation") == donor
            seen += 1
    assert seen > 5                                             # the run really did transfer


def test_transfers_and_losses_are_counted_as_the_event_log_counts_them(run):
    # a transfer writes two event rows for one gene-tree node, so the format must show ONE
    # branchingOut per transfer, not two (the error class the `event` column exists to kill)
    xml = recphylo_xml(run.gene_trees, run.complete_tree)
    tags = collections.Counter(ev.tag for rgt in ET.fromstring(xml).findall("recGeneTree")
                               for c in _clades(rgt) for ev in c.find("eventsRec"))
    log = collections.Counter(e.kind for e in run.events)
    assert tags["branchingOut"] == tags["transferBack"] == log["transfer_additive"]
    assert tags["duplication"] == log["duplication"]
    assert tags["loss"] == log["loss"]
    assert tags["speciation"] == log["speciation"]


def test_losses_are_in_the_file_at_all(run):
    # the reason the COMPLETE gene tree goes in: a loss has nothing to hang on in the extant tree
    xml = recphylo_xml(run.gene_trees, run.complete_tree)
    assert "<loss " in xml and xml.count("<loss ") > 5


def test_write_recphylo_writes_one_file_per_family(run, tmp_path):
    what = write_recphylo(run.gene_trees, run.complete_tree, tmp_path)
    files = sorted(p.name for p in tmp_path.iterdir())
    assert files == sorted(f"recphylo_fam{f}.xml" for f in run.gene_trees)
    assert what == f"{len(run.gene_trees)} file(s)"
    for p in tmp_path.iterdir():                                # each stands alone: tree + one family
        root = ET.fromstring(p.read_text(encoding="utf-8"))
        assert root.find("spTree") is not None
        assert len(root.findall("recGeneTree")) == 1


def test_the_cli_writes_recphylo_for_a_run(tmp_path):
    from zombi2.cli.main import main

    run_dir = tmp_path / "r"
    assert main(["species", str(run_dir), "--birth", "1", "--death", "0.3", "--n-extant", "8",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run_dir), "--duplication", "0.2", "--transfer", "0.2", "--loss",
                 "0.2", "--initial-families", "5", "--seed", "1", "--quiet"]) == 0
    assert main(["tools", "format", str(run_dir), "--format", "recphylo", "--quiet"]) == 0
    out = run_dir / "genomes" / "recphylo"
    assert out.is_dir() and list(out.glob("recphylo_fam*.xml"))
    for p in out.iterdir():
        ET.fromstring(p.read_text(encoding="utf-8"))                            # every file parses
