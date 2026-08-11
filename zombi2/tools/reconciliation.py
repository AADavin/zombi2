"""The extant-only reconciliation — the gene tree a dataset holds, drawn inside the species tree a
dataset holds.

`zombi2.tools.recphylo` writes the **complete** history: every lineage that ever lived, every gene
copy that ever existed. That is the truth, and it is not the thing a reconciliation method is scored
against, because a method sees only the survivors. Scoring needs the same history *projected* onto
what a dataset can contain — and the projection is not mechanical, which is why it lives in its own
module rather than as a flag deep inside the writer.

**What survives the projection.** A branch of the complete species tree is visible when it has an
extant descendant; a gene event is visible when it happened on one. Everything else is not merely
unobserved, it is unobservable, and putting it in a benchmark's answer key would score methods on
evidence that never existed. So:

- a **speciation** where the gene followed one daughter is kept, with a ``loss`` on the other, when
  that other daughter is itself visible — that is the one loss anyone can infer, and it stands for
  however many losses really happened down inside it;
- the same speciation is **suppressed** when the abandoned daughter has no extant descendant: the
  split is not in the extant species tree at all, so there is nothing to lose the gene on;
- a **duplication** whose second copy died is suppressed. Both copies sat on one branch, so the
  extant data shows one copy and no method could ever say otherwise;
- a **transfer** whose donor side died is suppressed, which is what leaves the arriving copy looking
  as though it appeared from nowhere — see ``entered_by`` below;
- a surviving transfer's **donor** is the nearest visible ancestor of the branch it really left. That
  branch had already split by then, so the donation is attributed to a branch that no longer existed
  at that instant. It is where an inference method would put it, which is the honest place for it.

**Where the family enters, and the one real choice.** Two reconciliations are defensible and they
differ in exactly one thing — the branch the family is rooted on:

- ``"true"`` roots it where the family really originated. A family present at the root and surviving
  in a scattered handful of genomes is *recorded as* present at the root, with the losses that
  narrowed it. This is what ancestral gene-content reconstruction is scored against: it is the answer
  key, and it is the only one of the two that contains the answer.
- ``"recoverable"`` roots it at the surviving copies' common ancestor — the most any method could
  recover from extant data alone. A family that left no trace above that point cannot be placed
  higher by any amount of cleverness, so a method graded against this one is never marked wrong for
  missing something invisible.

``"true"`` contains ``"recoverable"``: trimming the ancestral chain back to the surviving MRCA and
dropping the losses that go with it is mechanical, and going the other way is impossible. Both are
written, because the difference between them is itself worth seeing — it is exactly the part of the
history that a perfect method still cannot reach.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..genomes.gene_trees import GeneNode, GeneTree
from ..tree import Tree, prune

#: how a family enters the visible tree: it really began there, or it arrived from a lineage nobody
#: can see. The second is invisible in the reconciliation itself — the arrival looks like an
#: origination — so it is recorded beside it rather than faked into the XML, which has no tag for
#: either (a recPhyloXML gene tree simply starts mid-branch).
ORIGINATION, TRANSFER = "origination", "transfer"


@dataclass(frozen=True)
class Reconciliation:
    """One family's extant-only reconciliation. ``root`` is a gene tree in the ordinary `GeneNode`
    shape whose every ``species`` is a node of the **extant** species tree, so the recPhyloXML writer
    serialises it with no special case. ``entered_by`` is `ORIGINATION` or `TRANSFER` and ``branch``
    is where it entered; ``losses`` counts the synthesised loss leaves, which is what the two scopes
    differ by."""

    root: GeneNode
    entered_by: str
    branch: int
    losses: int


def visible_branches(tree: Tree) -> tuple[set[int], dict[int, int]]:
    """``(kept, image)`` for a complete species tree: ``kept`` is the extant tree's node set, and
    ``image`` sends every complete-tree branch that has an extant descendant to the extant-tree branch
    containing it. A branch with no extant descendant is absent from ``image`` — it has no image,
    because a dataset has no way to refer to it.

    A chain of complete-tree branches collapses to one extant-tree branch when the splits along it
    led nowhere, so several keys share an image. That is the whole reason a projected event can need
    moving: its own branch may not be a thing the extant tree names."""
    et = prune(tree, keep="extant")
    if et is None:
        return set(), {}
    kept = set(et.nodes)
    surviving: dict[int, bool] = {}
    for i in sorted(tree.nodes, reverse=True):        # children have higher ids than their parent
        nd = tree.nodes[i]
        surviving[i] = (nd.fate == "extant" if not nd.children
                        else any(surviving[c] for c in nd.children))
    image: dict[int, int] = {}
    for i in sorted(tree.nodes, reverse=True):        # children first, so a chain resolves in one pass
        if not surviving[i]:
            continue
        if i in kept:
            image[i] = i
        else:                                          # suppressed: exactly one child survives
            live = [c for c in (tree.nodes[i].children or ()) if surviving[c]]
            image[i] = image[live[0]]
    return kept, image


def _nearest_visible(branch: int, tree: Tree, image: dict[int, int]) -> int | None:
    """The branch's own image, or the nearest ancestor that has one — where a transfer out of a
    lineage nobody can see gets attributed. ``None`` only when nothing above it is visible either,
    which cannot happen below a surviving root."""
    i: int | None = branch
    while i is not None:
        if i in image:
            return image[i]
        i = tree.nodes[i].parent
    return None


def _loss(branch: int, time: float, copy: int) -> GeneNode:
    """A synthesised ``loss`` leaf on ``branch``. It carries the copy id of the lineage that failed to
    reach it, so the leaf still names a real gene rather than an invented one."""
    return GeneNode("loss", branch, time, copy)


def _moved_child(n: GeneNode) -> GeneNode | None:
    """Of a transfer's two children, the copy that **left** — the one now on another branch. A
    self-transfer puts both on one branch, and there it is the second, since the engine logs the
    donor's continuation first."""
    if n.kind != "transfer" or not n.children:
        return None
    return next((c for c in n.children if c.species != n.species), n.children[-1])


def _project(node: GeneNode, tree: Tree, image: dict[int, int],
             counter: list[int]) -> tuple[GeneNode | None, int | None]:
    """The projection proper, bottom-up over the complete gene tree.

    Returns ``(subtree, arrived_from)``. ``arrived_from`` is the visible donor branch when the link
    from this subtree to its parent is a **horizontal** one — which happens when a transfer out of a
    lineage nobody can see got suppressed on the way up. Carrying it is what stops the arriving copy
    being grafted under an ancestor as though it had descended from it: the node where it rejoins its
    visible relatives becomes the transfer, attributed to the nearest visible ancestor of the branch
    it really left. That is where an inference method would put it, and it is the only placement the
    extant tree can express."""
    order: list[GeneNode] = []
    stack = [node]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(n.children)

    done: dict[int, tuple[GeneNode | None, int | None]] = {}
    for n in reversed(order):                          # children before parents
        if n.is_leaf:
            done[id(n)] = ((GeneNode(n.kind, image[n.species], n.time, n.copy), None)
                           if n.kind == "extant" else (None, None))
            continue
        here = _nearest_visible(n.species, tree, image)
        assert here is not None          # a surviving root is visible, so something above is
        moved = _moved_child(n)
        live = [(c, *done[id(c)]) for c in n.children]
        alive = [(c, built, came) for c, built, came in live if built is not None]
        if not alive:
            done[id(n)] = (None, None)
            continue
        if len(alive) == 1:
            (c, built, came), = alive
            # this edge is horizontal if the survivor is the copy that left, or if it already was
            horizontal = here if c is moved else came
            dead = [d for d, b, _ in live if b is None]
            # a loss is visible only where the abandoned branch is one the extant tree names and is
            # not the branch its sibling still sits on: a dead duplicate leaves nothing to see, and a
            # dead donor branch is not in the extant tree at all
            visible_dead = [d for d in dead if d.species in image and image[d.species] != here]
            if n.kind == "speciation" and visible_dead and horizontal is None:
                counter[0] += len(visible_dead)
                out = GeneNode(n.kind, here, n.time, n.copy)
                out.children = [built] + [_loss(image[d.species], d.time, d.copy)
                                          for d in visible_dead]
                done[id(n)] = (out, None)
            else:
                done[id(n)] = (built, horizontal)      # suppressed; carry the jump up
            continue
        # two or more survivors: a real divergence. If either arrived horizontally, this is where the
        # transfer becomes visible, so the node IS the transfer — sitting on the donor's branch.
        # The donor is `here`, not the branch the copy really left. The branch it left is invisible,
        # and its nearest visible ancestor need not sit below this node — attributing the donation
        # there would put a transfer above its own parent and break the tree. `here` is the visible
        # common ancestor of both surviving groups, so the gene demonstrably WAS on it, and the
        # donation came from a descendant of it. Same smear in time as any ghost donor, and the only
        # placement that is also well formed.
        arrived = any(came is not None for _, _, came in alive)
        out = GeneNode("transfer" if arrived else n.kind, here, n.time, n.copy)
        out.children = [built for _, built, _ in alive]
        done[id(n)] = (out, None)
    return done[id(node)]


def _trim_to_recoverable(root: GeneNode) -> tuple[GeneNode, int]:
    """Drop the family's ancestral presence back to the surviving copies' common ancestor, and the
    losses that go with it. Returns ``(root, losses_dropped)``.

    The projection walks the **complete** gene tree, so it already carries every ancestral speciation
    whose abandoned daughter is visible — that is the ``"true"`` scope, built for free. The
    ``"recoverable"`` one is this trim: while the root's only real child is one subtree and everything
    beside it is a loss, the family's presence up there is a thing no method could recover, so cut it.
    One direction only — the ancestry, once cut, cannot be put back."""
    dropped = 0
    while True:
        real = [c for c in root.children if c.kind != "loss"]
        if len(real) != 1:
            return root, dropped
        dropped += len(root.children) - 1
        root = real[0]


def extant_reconciliation(gene_tree: GeneTree, tree: Tree, *, scope: str = "true",
                          visible: tuple[set[int], dict[int, int]] | None = None
                          ) -> Reconciliation | None:
    """One family's extant-only reconciliation, or ``None`` when the family left no extant copy.

    ``scope`` is ``"true"`` (rooted where the family really originated) or ``"recoverable"`` (rooted
    at the surviving copies' common ancestor) — see the module docstring for why both exist. Pass
    ``visible`` (from `visible_branches`) to share one species-tree analysis across a whole run."""
    if scope not in ("true", "recoverable"):
        raise ValueError(f"scope is 'true' or 'recoverable', got {scope!r}")
    if gene_tree.extant is None:
        return None
    _, image = visible if visible is not None else visible_branches(tree)
    counter = [0]
    root, arrived = _project(gene_tree.complete, tree, image, counter)
    if root is None:                                   # cannot happen with an extant tree, but say so
        return None
    # `arrived` survives all the way up exactly when the family's whole visible history hangs off one
    # horizontal jump — nothing it could have rejoined, because there was nothing else left. Then the
    # family did not begin here at all; it landed here, which is what a reader of the extant data
    # would call an origination and would be wrong about.
    entered_by = TRANSFER if arrived is not None else ORIGINATION
    losses = counter[0]
    if scope == "recoverable":
        root, dropped = _trim_to_recoverable(root)
        losses -= dropped
    return Reconciliation(root, entered_by, root.species, losses)


def origins_tsv(recs: dict[int, Reconciliation]) -> str:
    """The companion table: how each family entered the visible tree, and where.

    It is a table and not an attribute in the XML because recPhyloXML has no notion of origination at
    all — a gene tree there simply starts part-way down a branch — so there is nowhere legal to put
    it, and bending the format would cost us the readers that validate against it."""
    rows = ["family\tentered_by\tbranch\tlosses"]
    for fam, r in sorted(recs.items()):
        rows.append(f"{fam}\t{r.entered_by}\tn{r.branch}\t{r.losses}")
    return "\n".join(rows) + "\n"


__all__ = ["Reconciliation", "extant_reconciliation", "visible_branches", "origins_tsv",
           "ORIGINATION", "TRANSFER"]
