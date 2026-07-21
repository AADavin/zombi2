"""A byte-level pin on the nucleotide engine's seeded output.

The engine draws every choice from one rng in one order, so *any* change to how many values it
draws, or when, silently reshuffles every seeded run — the results stay valid but stop matching
what was published, replayed or reported. The ordinary determinism tests compare a run to itself,
which cannot see that: both halves move together.

This compares against a stored fixture instead, so it holds the engine to the draw order it had
when the fixture was recorded. Regenerate it deliberately, never to make a red test go green:

    python tests/test_genomes_nucleotide_golden.py

It was written to protect the extraction of the mutators' apply halves (`duplicate`, `delete`,
`transpose`, `originate` alongside `invert`), whose whole risk was disturbing that order.
"""

import pathlib

from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.species import read_newick

GOLDEN = pathlib.Path(__file__).parent / "data" / "nucleotide_golden.txt"

#: a fixed tree, so the fixture does not depend on the species engine's draws either
TREE = "(((A:0.4,B:0.4):0.3,C:0.7):0.3,(D:0.5,E:0.5):0.5);"

#: runs chosen to fire every mutator: content events, rearrangements, and the chromosome tier
RUNS = {
    "content": dict(duplication=0.6, loss=0.5, origination=0.7, transfer=0.4, root_length=400,
                    genes=3),
    "rearrangements": dict(inversion=0.8, transposition=0.6, translocation=0.5, chromosomes=3,
                           inversion_probability=0.5, root_length=300, genes=2),
    "tier": dict(fission=2.0, fusion=0.4, chromosome_origination=0.3, chromosome_loss=0.2,
                 chromosomes=2, root_length=300, genes=2),
    "everything": dict(duplication=0.4, loss=0.3, origination=0.5, transfer=0.3, inversion=0.5,
                       transposition=0.4, translocation=0.3, fission=1.0, fusion=0.4,
                       chromosomes=3, root_length=400, genes=3),
}


def _render(name, params):
    """One run, flattened to text: every node's blocks, then every event and rearrangement."""
    tree, _ = read_newick(TREE)
    r = simulate_genomes_nucleotide(tree, seed=11, **params)
    lines = [f"## {name}"]
    for node in sorted(r.genomes):
        for chrom in r.genomes[node].chromosomes:
            blocks = " ".join(f"{b.source}:{b.start}-{b.end}:{b.strand}:{b.copy}:{b.gene}"
                              for b in chrom.blocks)
            lines.append(f"node {node} chrom {chrom.id} {chrom.topology} {blocks}")
    lines += [f"event {e}" for e in r.events]
    lines += [f"rearrangement {e}" for e in r.rearrangements]
    lines += [f"chromosome_event {e}" for e in r.chromosome_events]
    return lines


def _render_all():
    out = []
    for name, params in RUNS.items():
        out += _render(name, params)
    return "\n".join(out) + "\n"


def test_seeded_runs_match_the_recorded_output():
    assert GOLDEN.exists(), f"missing fixture — regenerate with: python {__file__}"
    assert _render_all() == GOLDEN.read_text(), (
        "a seeded nucleotide run no longer matches the recorded output. If the engine's sampling "
        "genuinely changed, regenerate the fixture; if not, an rng draw was added, removed or "
        "reordered.")


def test_the_fixture_covers_every_mutator():
    # a guard on the fixture: if the runs above stopped firing an event kind, the pin above would
    # still pass while protecting nothing
    text = GOLDEN.read_text()
    for kind in ("Duplication", "Loss", "Origination", "Transfer", "Inversion", "Transposition",
                 "Translocation"):
        assert kind in text, f"the fixture fires no {kind} — it is not pinning that mutator"
    for tier in ("fission", "fusion"):
        assert tier in text, f"the fixture fires no {tier}"


if __name__ == "__main__":
    GOLDEN.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN.write_text(_render_all())
    print(f"wrote {GOLDEN} ({len(GOLDEN.read_text().splitlines())} lines)")
