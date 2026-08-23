"""Who trades genes with whom: the simulation and its truth tables.

One conditioned run. A two-state habitat is grown first on the complete tree; genomes then
evolve 1,000 families under constant duplication, transfer and loss, with the transfer
recipients weighted toward the donor's habitat:

    transfer_to = Recipients().weighted_by(habitat,
                      Between({("A", "A"): 10.0, ("B", "B"): 10.0}, default=1.0))

The run knows every transfer (donor lineage, recipient lineage, time), the habitat of every
lineage at every instant, and which gene copies survive to the present. This script writes:

    data/species_extant.nwk     the extant species tree (internal nodes labelled n<id>)
    data/gene_trees/fam_*.nwk   the extant gene tree of every family in >= MIN_SPECIES species
    data/transfers.tsv          every true transfer with its visibility and projections
    data/branches.tsv           every extant-tree branch with its true habitat and clade
    data/manifest.json          parameters, seeds, versions, counts

`oracle.py` reads these tables and reports what a perfect reconciliation could recover.
Everything is determined by the seeds; rerunning regenerates the data byte-identically.

    python experiment.py
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

sys.setrecursionlimit(1_000_000)

import zombi2
from zombi2.genomes import simulate_genomes_family
from zombi2.params import Between, Recipients
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

HERE = pathlib.Path(__file__).parent
DATA = HERE / "data"

BIRTH, DEATH, N_TIPS, SEED_TREE = 1.0, 0.5, 100, 11
STATES, SWITCH, SEED_TRAIT = ("A", "B"), 0.05, 1
N_FAMILIES, DUPLICATION, TRANSFER, LOSS, SEED_GENOMES = 1000, 0.02, 0.05, 0.15, 7
SAME_HABITAT_WEIGHT = 10.0
MIN_SPECIES = 4          # a gene tree needs a few leaves before reconciling it means anything


def build_world():
    sp = simulate_species_tree(birth=BIRTH, death=DEATH, n_extant=N_TIPS, seed=SEED_TREE)
    ct = sp.complete_tree
    tr = simulate_discrete(ct, states=list(STATES), start="A", seed=SEED_TRAIT,
                           switch={"A->B": SWITCH, "B->A": SWITCH})
    weights = Recipients().weighted_by(
        tr, Between({("A", "A"): SAME_HABITAT_WEIGHT, ("B", "B"): SAME_HABITAT_WEIGHT},
                    default=1.0))
    g = simulate_genomes_family(ct, initial_families=N_FAMILIES, duplication=DUPLICATION,
                                transfer=TRANSFER, loss=LOSS, transfer_to=weights,
                                seed=SEED_GENOMES)
    return sp, ct, tr, g


class Projection:
    """Everything about mapping complete-tree lineages onto the extant tree.

    An extant-tree branch is named by its tipward node, whose id is also its label in
    ``species_extant.nwk``. A living lineage projects to the branch it lies on; a dead one
    to the branch carrying its divergence from the sampled tree.
    """

    def __init__(self, ct, tr):
        self.ct, self.tr = ct, tr
        self.extant = set(ct.extant_leaves())
        self._flag, self._rep, self._chain = {}, {}, {}

    def _children(self, n):
        return [c if isinstance(c, int) else c.id for c in self.ct.nodes[n].children]

    def _parent(self, n):
        p = self.ct.nodes[n].parent
        return p if (p is None or isinstance(p, int)) else p.id

    def has_extant(self, n):
        if n not in self._flag:
            self._flag[n] = n in self.extant or any(
                self.has_extant(c) for c in self._children(n))
        return self._flag[n]

    def is_branch_node(self, n):
        return self.has_extant(n) and (
            n in self.extant
            or sum(1 for c in self._children(n) if self.has_extant(c)) >= 2)

    def tipward_rep(self, n):
        if n not in self._rep:
            self._rep[n] = n if self.is_branch_node(n) else self.tipward_rep(
                next(c for c in self._children(n) if self.has_extant(c)))
        return self._rep[n]

    def project(self, n):
        while not self.has_extant(n):
            n = self._parent(n)
            if n is None:
                return None
        return self.tipward_rep(n)

    def chain(self, rep):
        """The complete-tree lineages an extant-tree branch absorbs, tipward first."""
        if rep not in self._chain:
            members, n = [rep], self._parent(rep)
            while n is not None and self.has_extant(n) and not self.is_branch_node(n):
                members.append(n)
                n = self._parent(n)
            self._chain[rep] = members
        return self._chain[rep]

    def habitat(self, lineage, t):
        """Habitat of a lineage at ABSOLUTE time t (segment ends are branch-local)."""
        t_local = t - self.ct.nodes[lineage].birth_time
        for state, end in self.tr.history[lineage]:
            if t_local <= end + 1e-12:
                return state
        return self.tr.history[lineage][-1][0]

    def habitat_of_branch_at(self, rep, t):
        for n in self.chain(rep):
            node = self.ct.nodes[n]
            if node.birth_time - 1e-12 <= t <= node.end_time + 1e-12:
                return self.habitat(n, t)
        oldest = self.chain(rep)[-1]
        return self.habitat(oldest, max(t, self.ct.nodes[oldest].birth_time))

    def branch_truth(self, rep):
        """Duration-majority habitat over the branch, and whether it switched inside."""
        from collections import Counter
        dur = Counter()
        for n in self.chain(rep):
            node = self.ct.nodes[n]
            t0 = node.birth_time
            for state, end in self.tr.history[n]:
                end_abs = min(node.birth_time + end, node.end_time)
                if end_abs > t0:
                    dur[state] += end_abs - t0
                    t0 = end_abs
        return dur.most_common(1)[0][0], sum(1 for d in dur.values() if d > 0) > 1


def surviving_copies(gene_trees):
    """{family: set of copy ids alive at the present}, read off the extant gene trees."""
    out = {}
    for fam, tree in gene_trees.items():
        nw = tree.to_newick("extant", annotate=False) or ""
        out[fam] = {int(m) for m in re.findall(r"_g(\d+)[,:)]", nw + ")")}
    return out


def main() -> int:
    sp, ct, tr, g = build_world()
    proj = Projection(ct, tr)
    node_ids = list(range(len(ct.nodes))) if not isinstance(ct.nodes, dict) else list(ct.nodes)
    for n in node_ids:
        proj.has_extant(n)

    DATA.mkdir(exist_ok=True)
    (DATA / "gene_trees").mkdir(exist_ok=True)

    # --- the extant species tree, whose internal labels are the branch names below ---
    nwk = sp.extant_tree.to_newick()
    (DATA / "species_extant.nwk").write_text(nwk + "\n")
    # ALE's newick parser refuses internal node labels, so it gets an unlabeled copy;
    # its nodes are mapped back to ours by clade content (the clades are in branches.tsv).
    (DATA / "species_extant_ale.nwk").write_text(
        re.sub(r"\)n\d+", ")", nwk) + "\n")

    # --- the extant gene trees, one file per family with enough species ---
    lab = ct.labels()
    profiles = g.profiles
    import numpy as np
    M = np.asarray(profiles.matrix)
    fam_index = {fam: i for i, fam in enumerate(profiles.families)}
    written = []
    for fam, tree in g.gene_trees.items():
        row = fam_index.get(fam)
        if row is None or int((M[row] > 0).sum()) < MIN_SPECIES:
            continue
        nw = tree.to_newick("extant", annotate=False)
        if nw:
            (DATA / "gene_trees" / f"fam_{fam}.nwk").write_text(nw + "\n")
            written.append(fam)

    # --- gene-side visibility of every transfer ---
    by_family = {}
    for e in g.edges:
        by_family.setdefault(e.family, []).append(e)
    surviving = surviving_copies(g.gene_trees)

    def descends_from(c, anc, parent_of):
        while c is not None:
            if c == anc:
                return True
            c = parent_of.get(c)
        return False

    rows = []
    for fam, edges in by_family.items():
        parent_of = {e.copy: e.parent for e in edges if e.kind != "loss"}
        alive = surviving.get(fam, set())
        has_desc = set()
        for c in alive:
            while c is not None and c not in has_desc:
                has_desc.add(c)
                c = parent_of.get(c)
        for e in edges:
            if e.kind != "transfer" or e.recipient is None:
                continue
            trace = e.copy in has_desc
            detectable = trace and any(
                not descends_from(s, e.copy, parent_of) for s in alive)
            pd, pr = proj.project(e.donor), proj.project(e.recipient)
            rows.append({
                "family": fam, "time": f"{e.time:.6f}",
                "donor": e.donor, "recipient": e.recipient,
                "donor_alive": int(proj.has_extant(e.donor)),
                "recipient_alive": int(proj.has_extant(e.recipient)),
                "trace": int(trace), "detectable": int(detectable),
                "donor_proj": "" if pd is None else pd,
                "recipient_proj": "" if pr is None else pr,
                "donor_habitat": proj.habitat(e.donor, e.time),
                "recipient_habitat": proj.habitat(e.recipient, e.time),
                "donor_proj_habitat":
                    "" if pd is None else proj.habitat_of_branch_at(pd, e.time),
                "recipient_proj_habitat":
                    "" if pr is None else proj.habitat_of_branch_at(pr, e.time),
            })
    cols = list(rows[0])
    with open(DATA / "transfers.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # --- every extant-tree branch: truth habitat, mid-branch switching, clade ---
    reps = [n for n in node_ids if proj.is_branch_node(n)]
    leaves_below = {}

    def clade(n):
        if n not in leaves_below:
            leaves_below[n] = ({lab[n]} if n in proj.extant else
                               set().union(*(clade(c) for c in proj._children(n)
                                             if proj.has_extant(c))))
        return leaves_below[n]

    with open(DATA / "branches.tsv", "w") as fh:
        fh.write("branch\tkind\ttruth_habitat\tswitched_mid_branch\ttip_state\tclade\n")
        for rep in sorted(reps):
            hab, switched = proj.branch_truth(rep)
            terminal = rep in proj.extant
            tip_state = tr.history[rep][-1][0] if terminal else ""
            fh.write(f"{rep}\t{'terminal' if terminal else 'internal'}\t{hab}"
                     f"\t{int(switched)}\t{tip_state}\t{';'.join(sorted(clade(rep)))}\n")

    n_transfers = len(rows)
    manifest = {
        "zombi2": zombi2.__version__,
        "tree": {"birth": BIRTH, "death": DEATH, "n_extant": N_TIPS, "seed": SEED_TREE},
        "habitat": {"states": list(STATES), "switch": SWITCH, "seed": SEED_TRAIT},
        "genomes": {"initial_families": N_FAMILIES, "duplication": DUPLICATION,
                    "transfer": TRANSFER, "loss": LOSS, "seed": SEED_GENOMES,
                    "same_habitat_weight": SAME_HABITAT_WEIGHT},
        "min_species_per_gene_tree": MIN_SPECIES,
        "counts": {"extinct_lineages": len(ct.extinct_leaves()),
                   "surviving_families": int(M.shape[0]),
                   "gene_trees_written": len(written),
                   "true_transfers": n_transfers},
    }
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {len(written)} gene trees, {n_transfers} transfers, "
          f"{len(reps)} branches -> {DATA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
