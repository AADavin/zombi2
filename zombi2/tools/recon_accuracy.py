"""Accuracy of an inferred reconciliation against a known (simulated) one.

The mirror of ALElite (:mod:`zombi2.tools.reconciliation`): ALElite asks *how probable* a
reconciliation is; this asks *how close an inferred reconciliation is to the truth*. Both read
ZOMBI2's reconciled gene trees; this one takes two of them — the true reconciliation a
simulation emits and an inferred reconciliation for the **same gene tree** — and reports, node
by node, how much of the event history was recovered.

Scope: **fixed topology, node-by-node.** Both reconciliations must annotate the *same* extant
gene tree (identical tip labels, identical branching), as when an inference method (ALE,
GeneRax, ecceTERA, …) is run on the fixed simulated gene tree. The tool aligns the two trees
structurally — matching children by the set of tips beneath them, so it is insensitive to
Newick child order and handles the unary pseudogenization nodes ZOMBI2 keeps in the extant
tree — and compares the annotation at each internal node:

* **event type** — is each node's inferred event (S/D/T, plus G pseudogenization / C
  conversion) the true one? Reported as overall accuracy and per-class precision / recall,
  the standard way to score D/T/L recovery against a simulated truth.
* **species mapping** — does the node map to the right species branch? This is the LCA / MRCA
  mapping accuracy.
* **transfer donor & recipient** — for each true transfer, was it detected, and were its donor
  and recipient branches recovered?

Losses leave no node in the extant tree, so they fall outside node-by-node scoring (a
per-species-branch loss-count comparison would be a separate, topology-agnostic tool). Inputs
are ZOMBI2 annotated reconciled Newicks — a :class:`~zombi2.genomes.reconciliation.Reconciliation`
(its ``.extant`` string is used) or the string itself, with internal labels
``"<species branch>|<EVENT>"`` (``"<donor>|T>recipient"`` for transfers) and tips
``"<species>|<gid>"``. Converting a third-party method's output into that format is the glue
the caller supplies; this tool scores, it does not run inference.
"""

from __future__ import annotations

from collections import Counter, namedtuple

from zombi2.tree import Tree, read_newick

__all__ = [
    "reconciliation_accuracy",
    "ReconAccuracy",
    "EventPR",
    "TransferRecovery",
]

#: Per-event-class precision / recall. ``support_true`` / ``support_pred`` are the node counts
#: labelled this event in the truth / the inference; ``tp`` is how many agree.
EventPR = namedtuple("EventPR", ["precision", "recall", "f1", "tp", "support_true", "support_pred"])

#: Transfer donor/recipient recovery. All fields are node counts: ``n_true`` true transfers,
#: ``detected`` of them also called a transfer, and — among detected — ``donor_correct`` /
#: ``recipient_correct`` / ``both_correct`` with the right donor / recipient / both branches.
TransferRecovery = namedtuple(
    "TransferRecovery",
    ["n_true", "detected", "donor_correct", "recipient_correct", "both_correct"],
)

#: Full result. Accuracies are fractions in ``[0, 1]`` over the ``n_nodes`` internal gene-tree
#: nodes; ``per_event`` maps each event char to an :class:`EventPR`; ``transfer`` is a
#: :class:`TransferRecovery` (``n_true == 0`` if the true tree has no transfers).
ReconAccuracy = namedtuple(
    "ReconAccuracy",
    ["n_nodes", "event_accuracy", "mapping_accuracy", "joint_accuracy", "per_event", "transfer"],
)

#: Internal-node annotation parsed from a reconciled-tree label.
_Ann = namedtuple("_Ann", ["event", "species", "recipient"])


def _extant_tree(x) -> Tree:
    """Accept a :class:`~zombi2.genomes.reconciliation.Reconciliation`, or an annotated Newick."""
    if isinstance(x, str):
        return read_newick(x)
    if hasattr(x, "extant"):
        if x.extant is None:                          # a Reconciliation of a fully extinct family
            raise ValueError("reconciliation has no extant tree (the family is fully extinct)")
        return read_newick(x.extant)
    raise TypeError(f"expected a Reconciliation or an annotated Newick string, got {type(x).__name__}")


def _parse_label(name: str) -> _Ann:
    """``"n5|D"`` -> event D; ``"n5|T>n8"`` -> transfer, recipient n8; ``"n5|G"`` -> pseudogenization."""
    if "|" not in name:
        return _Ann("?", name, None)                  # unannotated internal node
    species, spec = name.split("|", 1)
    if spec.startswith("T>"):
        return _Ann("T", species, spec[2:])
    return _Ann(spec, species, None)


def _clades(tree: Tree) -> dict:
    """``{node: frozenset(leaf names beneath it)}`` (post-order); leaf atoms are ``species|gid``."""
    clades: dict = {}
    for node in reversed(tree.nodes()):
        if node.is_leaf():
            clades[node] = frozenset((node.name,))
        else:
            acc: set = set()
            for c in node.children:
                acc |= clades[c]
            clades[node] = frozenset(acc)
    return clades


def _aligned_internal_pairs(truth: Tree, inferred: Tree) -> list:
    """Pairs of ``(_Ann, _Ann)`` for every aligned internal node, matched structurally.

    Children are matched by their tip-set, so Newick child order is irrelevant and unary
    (pseudogenization) nodes align to unary nodes. Raises if the two trees are not the same
    topology over the same tips — this tool's fixed-topology contract.
    """
    tc, ic = _clades(truth), _clades(inferred)
    if tc[truth.root] != ic[inferred.root]:
        only_t = sorted(tc[truth.root] - ic[inferred.root])[:5]
        only_i = sorted(ic[inferred.root] - tc[truth.root])[:5]
        raise ValueError(
            "the two reconciliations must annotate the same gene tree (same tip labels); "
            f"tips only in truth e.g. {only_t}, only in inferred e.g. {only_i}"
        )
    pairs: list = []

    def rec(tn, in_):
        if tn.is_leaf():
            return
        pairs.append((_parse_label(tn.name), _parse_label(in_.name)))
        t_by = {tc[c]: c for c in tn.children}
        i_by = {ic[c]: c for c in in_.children}
        if set(t_by) != set(i_by):
            raise ValueError(
                "the two reconciliations differ in gene-tree topology (a node's children do not "
                "match); this tool assumes a shared fixed topology"
            )
        for clade, tchild in t_by.items():
            rec(tchild, i_by[clade])

    rec(truth.root, inferred.root)
    return pairs


def _pr(tp: int, support_true: int, support_pred: int) -> EventPR:
    precision = tp / support_pred if support_pred else 0.0
    recall = tp / support_true if support_true else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return EventPR(precision, recall, f1, tp, support_true, support_pred)


def reconciliation_accuracy(truth, inferred) -> ReconAccuracy:
    """Node-by-node accuracy of an inferred reconciliation against a true one.

    ``truth`` and ``inferred`` are each a :class:`~zombi2.genomes.reconciliation.Reconciliation`
    (its ``.extant`` reconciled tree is scored) or an annotated reconciled-Newick string. They
    must annotate the **same** extant gene tree (same tip labels, same topology); a mismatch is
    an error, not a partial score. Returns a :class:`ReconAccuracy`.

    A gene tree with no internal nodes (a lone surviving copy) yields zero comparisons and
    all-zero accuracies.
    """
    pairs = _aligned_internal_pairs(_extant_tree(truth), _extant_tree(inferred))
    n = len(pairs)
    if n == 0:
        return ReconAccuracy(0, 0.0, 0.0, 0.0, {}, TransferRecovery(0, 0, 0, 0, 0))

    event_correct = mapping_correct = joint_correct = 0
    support_true: Counter = Counter()
    support_pred: Counter = Counter()
    tp: Counter = Counter()
    n_true_T = detected = donor_ok = recip_ok = both_ok = 0

    for t, i in pairs:
        ev_ok = t.event == i.event
        map_ok = t.species == i.species
        event_correct += ev_ok
        mapping_correct += map_ok
        joint_correct += ev_ok and map_ok
        support_true[t.event] += 1
        support_pred[i.event] += 1
        if ev_ok:
            tp[t.event] += 1
        if t.event == "T":
            n_true_T += 1
            if i.event == "T":
                detected += 1
                d_ok = i.species == t.species
                r_ok = i.recipient == t.recipient
                donor_ok += d_ok
                recip_ok += r_ok
                both_ok += d_ok and r_ok

    events = set(support_true) | set(support_pred)
    per_event = {ev: _pr(tp[ev], support_true[ev], support_pred[ev]) for ev in sorted(events)}
    return ReconAccuracy(
        n_nodes=n,
        event_accuracy=event_correct / n,
        mapping_accuracy=mapping_correct / n,
        joint_accuracy=joint_correct / n,
        per_event=per_event,
        transfer=TransferRecovery(n_true_T, detected, donor_ok, recip_ok, both_ok),
    )
