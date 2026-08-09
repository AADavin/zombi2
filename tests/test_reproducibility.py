"""What a seed guarantees.

ZOMBI2's one promise is that a seed names a run: give the same version the same seed and you get the
same dataset back, on any machine, any supported Python, any operating system. Everything downstream
rests on it — a paper's "seed 42" is not a citation unless the run behind it can be reproduced, and
a benchmark whose ground truth shifts between machines is not a benchmark.

Nothing else in the suite checks that. The other tests assert *properties* of a run — that a tree has
the right number of tips, that a rate does what it says — and every one of them would still pass if
the random stream shifted by a single draw and every number came out different.

So this file records **digests**. A short run of each level is dumped to a canonical text form and
hashed, and the hash is written down below. CI runs it on Linux, macOS and Windows across Python 3.10
to 3.13, which is what turns the promise into a checked fact rather than a hope: six configurations,
six independent Python hash seeds, one number.

**When one of these fails**, do not re-record it to make the suite green. The digest changing means
the run changed, which is either a bug or a deliberate change to the model. If it is deliberate, it
is a **breaking** change — every seed anyone published now means something different — so it belongs
in a minor release with a `CHANGELOG.md` entry saying so, and *then* the digest is updated.
"""

from __future__ import annotations

import hashlib

import pytest

from zombi2.genomes import simulate_genomes_family, simulate_genomes_ordered
from zombi2.rates import ScaledBy
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

# Floats are written at 9 significant digits, not repr. A double carries ~16, and the last one or two
# are where a platform's own libm can differ — `exp` and `log` are permitted a rounding error, and
# LAPACK is a different implementation on macOS than on Linux. Nine digits is far more precision than
# any result of ours is claimed to, and far enough from the noise to be stable.
_G = "%.9g"


def _fmt(v: object) -> str:
    return _G % v if isinstance(v, float) else str(v)


def _digest(*blocks: str) -> str:
    h = hashlib.sha256()
    for b in blocks:
        h.update(b.encode("utf-8"))
        h.update(b"\x00")
    return h.hexdigest()[:16]


def _tree_text(tree) -> str:
    """A tree as its Newick, at the precision every ZOMBI2 file is written to."""
    return tree.to_newick()


def _events_text(events) -> str:
    """A log as the columns that decide what happened, in the order it happened.

    The genome runs hash their `edges` — one entry per gene-tree edge — rather than `events`, which
    is the grouped view. Not a detail: the edge list is the finer record, so a change the grouping
    would hide still moves the digest, and it is what these digests hashed before `events` came to
    mean one row per event. They are unchanged across that rename, which is what says it moved no
    randomness."""
    return "\n".join("\t".join(_fmt(getattr(e, f, "")) for f in
                               ("time", "kind", "lineage", "family", "copy", "parent",
                                "donor", "recipient"))
                     for e in events)


# --- the runs ---------------------------------------------------------------------------------
#
# Small and quick, but each exercises the whole path down to its own level: the species engine feeds
# the genome engine feeds the sequence engine. A drift anywhere upstream changes every digest below
# it, which is the point — one of these going red localises the change by which ones stayed green.

def _species():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=42)
    return _digest(_tree_text(r.complete_tree), _tree_text(r.extant_tree),
                   "\n".join(f"{n.id}\t{n.fate}\t{_G % n.birth_time}\t{_G % n.end_time}"
                             for n in sorted(r.complete_tree.nodes.values(), key=lambda n: n.id)))


def _genomes_family():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=42)
    g = simulate_genomes_family(sp, initial_families=8, duplication=0.2, transfer=0.15, loss=0.25,
                                origination=0.3, seed=7)
    trees = "\n".join(f"{fam}\t{gt.to_newick('complete')}" for fam, gt in sorted(g.gene_trees.items()))
    return _digest(_events_text(g.edges), trees)


def _genomes_ordered():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=42)
    g = simulate_genomes_ordered(sp, initial_families=12, duplication=0.15, loss=0.15,
                                 inversion=0.3, transposition=0.1, seed=7)
    order = "\n".join(f"{s}\t{c.id}\t{c.topology}\t"
                      + ",".join(f"{gene.family}:{gene.strand}" for gene in c.genes)
                      for s in sorted(g.genomes) for c in g.genomes[s])
    return _digest(_events_text(g.edges), order)


def _sequences():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=42)
    g = simulate_genomes_family(sp, initial_families=4, duplication=0.15, loss=0.15, seed=7)
    s = simulate_sequences(g, model=hky85(kappa=2.5), length=120, divergence=0.4, seed=11)
    aln = "\n".join(f"{fam}\t{k}\t{v}" for fam, a in sorted(s.alignments.items())
                    for k, v in sorted(a.items()))
    return _digest(aln)


def _traits_continuous():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=42)
    t = simulate_continuous(sp.complete_tree, start=0.0, rate=1.5, seed=3)
    return _digest("\n".join(f"{k}\t{_fmt(v)}" for k, v in sorted(t.values.items())))


def _traits_discrete():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=42)
    t = simulate_discrete(sp.complete_tree, states=["a", "b", "c"], start="a", switch=0.4, seed=3)
    return _digest(_events_text(t.events),
                   "\n".join(f"{k}\t{v}" for k, v in sorted(t.values.items())))


def _driven():
    """A conditioned run: the trait is grown first and a genome rate reads it. Worth its own digest
    because the driver threading is a second stream through the same engine, and a change there would
    leave every undriven digest above untouched."""
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=15, seed=42)
    hab = simulate_discrete(sp.complete_tree, states=["wet", "dry"], start="wet", switch=0.3, seed=3)
    g = simulate_genomes_family(sp, initial_families=6, duplication=0.1,
                                loss=0.2 * ScaledBy(hab, {"dry": 5.0, "wet": 1.0}), seed=7)
    return _digest(_events_text(g.edges))


#: run name → the digest that run has produced since it was recorded. Read the module docstring
#: before changing one.
RUNS = {
    "species": (_species, "4f33df380040b900"),
    "genomes_family": (_genomes_family, "ea19d2eb7751ed4e"),
    "genomes_ordered": (_genomes_ordered, "60f87e542d759973"),
    "sequences": (_sequences, "5d6ae687bfb8c219"),
    "traits_continuous": (_traits_continuous, "719b66f31de35cd7"),
    "traits_discrete": (_traits_discrete, "4d46d145dd81174a"),
    "driven": (_driven, "421c647ee645e8c0"),
}


@pytest.mark.parametrize("name", sorted(RUNS))
def test_a_seed_still_names_the_same_run(name):
    run, expected = RUNS[name]
    assert run() == expected, (
        f"the {name!r} run no longer produces what it produced when this digest was recorded. A seed "
        f"names a run, so this means every published seed for this level now means something else. "
        f"If that was deliberate it is a breaking change: say so in CHANGELOG.md under a minor "
        f"version, then update the digest. If it was not, something upstream shifted the random "
        f"stream — find it rather than re-recording."
    )


@pytest.mark.parametrize("name", sorted(RUNS))
def test_a_seed_reproduces_a_run_within_one_process(name):
    """The same call twice, in one process. The digests above cannot catch a run that depends on
    something *other* than its seed — a leaked global, a cached result, an iteration order set by
    object identity — because each of them runs only once. This one does."""
    run, _ = RUNS[name]
    assert run() == run()


def test_two_seeds_produce_two_runs():
    """The other half of the promise, and the reason a digest test is not enough on its own: a
    simulator that ignored its seed entirely would satisfy every assertion above."""
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=42)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=43)
    assert a.complete_tree.to_newick() != b.complete_tree.to_newick()
