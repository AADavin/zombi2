"""Transfer mechanics shared across genome resolutions — the ``transfer_to`` weighting.

A transfer's *rate* is an ordinary rate; what is special is **who receives** once it fires. That
mechanic is the same whether the genome is a multiset of families or an ordered set of chromosomes,
so it lives here, imported by every resolution. ``transfer_to`` **chooses who receives** (SPEC §5) —
the numbers in it are per-candidate weights, normalised across the contemporaneous lineages, so they
change neither how fast nor how many transfers happen, only **who** receives. Four rules:

- ``"uniform"`` — every contemporaneous lineage gets equal weight;
- `Distance` — weight by relatedness (closer relatives likelier), which needs the tree's mean
  root-to-tip time to stay scale-free;
- `Clades` — weight by the **pair** (donor's named clade, recipient's named clade), a fact read from
  the tree;
- `ScaledBy` — weight by **another level**: candidate ``k``'s weight is
  the mapping of the driver's value on lineage ``k`` at this instant (a trait that makes a lineage
  competent to take DNA up), and with a `Between` mapping the donor's value too.

All four rules work at **every** resolution — family, ordered and nucleotide. That is what
`resolve_transfer_to()` and `prepare_transfer_to()` are for: it is validated and prepared here,
once, so the three engines cannot drift apart in what they accept or in what they say when they
refuse.
"""

from __future__ import annotations

import math

from ..rates.choice import Clades, Distance
from ..rates.mapping import Between, check_kernel_fires
from ..rates.modifiers import Driven, SetBy
from ..rates.rate import Rate
from .._runtime.draw import weighted_index as _weighted_index
from ..tree import node_label


def _resolve_node(tree, spec) -> int:
    """A lineage reference — an ``int`` node id or an ``"n<id>"`` / ``"<id>"`` string — to a node id in
    ``tree`` (the form ``species_events.tsv`` and ``to_newick`` write). Raises if it is neither, or not
    a node of this tree."""
    if isinstance(spec, bool):
        raise ValueError(f"a clade reference must be a node id or 'n<id>' label, got {spec!r}")
    if isinstance(spec, int):
        nid = spec
    elif isinstance(spec, str):
        s = spec.strip()
        s = s[1:] if s[:1] == "n" else s
        try:
            nid = int(s)
        except ValueError:
            raise ValueError(
                f"clade reference {spec!r} is not a node id — name a lineage by its integer id or its "
                f"'n<id>' label (as species_events.tsv / to_newick write it)") from None
    else:
        raise ValueError(f"a clade reference must be a node id (int) or an 'n<id>' label, got {spec!r}")
    if nid not in tree.nodes:
        raise ValueError(
            f"clade reference {node_label(nid, tree.nodes[nid].fate if nid in tree.nodes else None)}"
            f" is not a lineage of this tree")
    return nid


def _mrca(tree, node_ids) -> int:
    """The most recent common ancestor of ``node_ids`` in ``tree`` — the deepest node on every one of
    their root-ward paths. A single id is its own MRCA."""
    def ancestors(i):
        chain = []
        while i is not None:
            chain.append(i)
            i = tree.nodes[i].parent
        return chain

    common = ancestors(node_ids[0])  # deepest-first: self, parent, …, root
    for other in node_ids[1:]:
        others = set(ancestors(other))
        common = [a for a in common if a in others]  # keep depth order, drop non-shared
    return common[0]


def _subtree(tree, root_id) -> set:
    """Every node in the clade below ``root_id`` (inclusive), extinct and internal lineages included."""
    out, stack = set(), [root_id]
    while stack:
        i = stack.pop()
        out.add(i)
        kids = tree.nodes[i].children
        if kids is not None:
            stack.extend(kids)
    return out


def resolve_groups(tree, groups) -> dict:
    """Paint every node of the complete ``tree`` with its clade label — ``{node_id: label}``, ``"rest"``
    for a lineage in no named clade. A clade named by a list of tips is the subtree below their MRCA; a
    clade named by a single node id is that node's subtree. Clades must be **disjoint** (an overlap —
    one clade nested in another — is refused). Computed once per run; membership is constant along a
    branch, so unlike a driver it adds no Gillespie breakpoints."""
    group_of = {i: "rest" for i in tree.nodes}
    claimed: dict[int, str] = {}
    for label, spec in groups.items():
        root = _mrca(tree, [_resolve_node(tree, t) for t in spec]) if isinstance(spec, (list, tuple)) \
            else _resolve_node(tree, spec)
        for i in _subtree(tree, root):
            if i in claimed and claimed[i] != label:
                raise ValueError(
                    f"clades {claimed[i]!r} and {label!r} overlap at "
                    f"{node_label(i, tree.nodes[i].fate)}; groups must be disjoint — "
                    f"is one clade nested inside the other?")
            claimed[i] = label
            group_of[i] = label
    return group_of


# --- who receives, in one place: validate it once, prepare it once ------------------------------------
# Three engines take ``transfer_to``, and the words they use to refuse a bad one — and the work they
# do before the first transfer fires — must be the same words and the same work, or the resolutions
# quietly become three slightly different models. So both live here and each engine calls them.

def resolve_transfer_to(transfer_to):
    """Validate ``transfer_to`` — who receives — and return the rule the engine will run on:
    ``"uniform"``, a `Distance`, a `Clades` or a `Driven` — with the ``"distance"`` shorthand
    coerced to ``Distance()``.

    This is not a rate (SPEC §5). The numbers in it are per-candidate **weights**, normalised
    across the contemporaneous candidates, so they change neither how fast nor how many transfers
    happen — only **who** receives. Two of the four messages below exist because that distinction is
    exactly what a user coming from the rate grammar gets wrong, and the generic "must be one of …"
    message would list the alternatives without saying why what they wrote is a different kind of
    thing:

    - ``transfer_to = 0.1 * Weights(...)`` is the *rate* spelling. There is no base here, because
      there is no rate, so the modifier is given on its own.
    - ``transfer_to = (Distance(), Weights(...))`` asks for two rules at once. Composing a
      topological weight with a driven one is a later slice, not a thing that is refused on principle.
    """
    if transfer_to == "distance":
        return Distance()
    if isinstance(transfer_to, Rate):
        raise ValueError(
            "transfer_to takes the Weights modifier on its own, not a rate — write "
            "transfer_to=Weights(driver, {...}) with no base number. Here the mapping's "
            "numbers are relative WEIGHTS over the candidate recipients (normalised), not a rate "
            "multiplier: they change which lineage receives, never how often transfer happens."
        )
    if isinstance(transfer_to, (list, tuple)):
        raise ValueError(
            "transfer_to takes one recipient rule, not several — combining Distance (relatedness) "
            "with a Weights rule is a later slice. Give 'uniform', 'distance' / "
            "Distance(decay=), or ScaledBy(driver, {...})."
        )
    if isinstance(transfer_to, Driven) and type(transfer_to) is Driven and transfer_to.verb == "ScaledBy":
        # Now that the verb carries the meaning, a mismatched one is catchable: `ScaledBy` says the
        # number multiplies a base, and a choice has none. Same object, wrong word — and saying so
        # is cheaper than a reader wondering why their "factor" behaved like a weight.
        raise ValueError(
            "transfer_to takes Weights(driver, {...}), not ScaledBy: its numbers are weights over "
            "the candidate recipients, compared against each other and normalised, rather than a "
            "factor on a base. The same driver and the same mapping — only the verb changes.")
    if isinstance(transfer_to, SetBy):
        # `SetBy` is a `Driven`, so the test below admits it — and admitting it runs the model the
        # user did not write. A choice has no base: only the ratios between the candidates are read,
        # so there is nothing here for a replaced base to mean, and it ran as an ordinary weighting.
        raise ValueError(
            "transfer_to cannot be SetBy: it weights the candidate recipients against each other, "
            "so there is no base for a driver to replace. The same numbers, spelled as what they "
            "are, is Weights(driver, {...}).")
    if transfer_to != "uniform" and not isinstance(transfer_to, (Distance, Driven, Clades)):
        raise ValueError(
            f"transfer_to must be 'uniform', 'distance' / Distance(decay=), "
            f"Clades({{...}}, Between({{...}})) (weight by named clade), or "
            f"ScaledBy(driver, {{...}}) (a recipient weight driven by an evolved value), "
            f"got {transfer_to!r}")
    return transfer_to


def prepare_transfer_to(tree, transfer_to, resolved=None, *, level=None):
    """Everything a run must work out **once**, before its first transfer, for the ``transfer_to``
    rule it was given — returned as the pair ``(group_of, to_traj)`` that `recipient_index()` takes
    as its ``groups`` and ``to_traj``. ``(None, None)`` for ``"uniform"`` and `Distance`, which need
    nothing prepared.

    - A `Clades` rule paints every lineage with its clade label. Membership is a fact about the
      **tree**, so it is constant along a branch and adds no Gillespie breakpoints — which is the
      whole reason it can be computed here and then never touched again.
    - A `ScaledBy` weight resolves its driver into a driver trajectory. ``resolved`` is the caller's
      ``{driver key: trajectory}`` cache, mutated in place, so a driver shared with a driven *rate*
      is loaded once and the two read the very same trajectory.

    **The trajectory returned here must not join the engine's ``trajs``.** A driven ``transfer_to``
    moves no rate: its weights are read at the instant a transfer fires, and a weight that is not a
    rate can neither change the total hazard nor make the loop per-lineage. Putting it in ``trajs``
    would add a Gillespie breakpoint at every one of the driver's switches, which changes the draws —
    and so the run — while every assertion about who received still passes. Call this **after**
    ``trajs`` is built, and hand the trajectory straight to ``_do_transfer``.

    Both branches also refuse a mapping that could never fire, for the same reason a driven rate
    does: a kernel or table naming only states the driver never takes leaves every candidate on the
    default weight, so the recipient is drawn uniformly while the run's log records it as steered.

    ``level`` names the calling engine (``"genomes.ordered"``, …), passed on so a driver this level may
    not read at all is refused here too — who receives is as much a read of the driver as a rate is.
    """
    if isinstance(transfer_to, Clades):
        group_of = resolve_groups(tree, transfer_to.groups)
        check_kernel_fires(transfer_to.between, set(group_of.values()), driver_label="clades")
        return group_of, None
    if isinstance(transfer_to, Driven):
        # imported here, not at module scope, so a run with no driver anywhere never pays for the
        # driver machinery (the same lazy import the family engine makes for its rate drivers)
        from ..rates.driver import check_mapping_fires, resolve_driver
        if resolved is None:
            resolved = {}
        if transfer_to.key not in resolved:
            resolved[transfer_to.key] = resolve_driver(transfer_to.driver, tree,
                                                       step=transfer_to.step, level=level)
        to_traj = resolved[transfer_to.key]
        label = transfer_to.driver if isinstance(transfer_to.driver, str) \
            else f"<{type(transfer_to.driver).__name__}>"
        if isinstance(transfer_to.mapping, Between):
            # a donor-conditioned kernel. The two checks are an either/or, not a pair:
            # check_mapping_fires only knows Table and would return without looking at a kernel.
            check_kernel_fires(transfer_to.mapping, to_traj.states(), driver_label=label)
        else:
            check_mapping_fires(transfer_to.mapping, to_traj.states(), driver_label=label)
        return None, to_traj
    return None, None


def recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth, to_traj=None, groups=None):
    """Pick a recipient lineage index (into ``alive``) from the candidate indices ``cand`` by the
    ``transfer_to`` rule: ``"uniform"`` gives every contemporaneous lineage equal weight; a
    `Distance` weights by relatedness (closer relatives likelier); a `Clades` weights by
    the kernel on (donor's clade, candidate's clade), read from the precomputed ``groups`` map; a
    `ScaledBy` weights by the driver's value on each candidate, read
    from ``to_traj`` (the trajectory the engine resolved for that driver) — and, with a
    `Between` mapping, by the donor's value too.

    Returns ``None`` — "nobody can receive" — when a driven weighting gives **every** candidate a
    weight of 0. The caller must then make the event a **no-op**: leaving it unrecorded is exactly the
    model in which the transfer rate itself drops to zero while no eligible recipient exists, because
    rejecting an event whose acceptance depends only on the current state is Poisson thinning, and a
    rejected event changes nothing (see `_do_transfer()`)."""
    if transfer_to == "uniform":
        return cand[int(rng.integers(len(cand)))]
    if isinstance(transfer_to, Clades):
        # topological, donor-conditioned: candidate k's weight is the kernel on (donor's clade, k's
        # clade), read from the precomputed membership map. A weight of 0 means "cannot receive".
        g_d = groups[donor]
        weights = [transfer_to.between.weight(g_d, groups[alive[k]]) for k in cand]
        total = sum(weights)
        if total <= 0.0:
            return None
        return cand[_weighted_index(rng, weights, total)]
    if isinstance(transfer_to, Driven):
        # who receives: candidate k's weight is the mapping of the driver on lineage k right now,
        # normalised over the candidates. A weight of 0 means "cannot receive". A Between mapping is
        # donor-conditioned — the weight reads the driver on the DONOR too — so a trait can steer
        # transfer between guilds exactly as Clades does between clades.
        if isinstance(transfer_to.mapping, Between):
            g_d = to_traj.value(donor, t)
            weights = [transfer_to.mapping.weight(g_d, to_traj.value(alive[k], t)) for k in cand]
        else:
            weights = [transfer_to.mapping.multiplier(to_traj.value(alive[k], t)) for k in cand]
        total = sum(weights)
        if total <= 0.0:
            return None
        return cand[_weighted_index(rng, weights, total)]
    # Distance: patristic distance d(donor, x) = 2·(t − t_mrca); scale-free in the tree depth. Mark
    # the donor's ancestor end-times once, then climb each candidate to its first marked ancestor.
    anc = {}
    p = tree.nodes[donor].parent
    while p is not None:
        anc[p] = tree.nodes[p].end_time
        p = tree.nodes[p].parent
    dists = []
    for k in cand:
        x = alive[k]
        if x == donor:
            dists.append(0.0)  # self (only reachable under self_transfer): closest
            continue
        q = x
        while q not in anc:
            q = tree.nodes[q].parent
        dists.append(2.0 * (t - anc[q]))
    dmin = min(dists)
    weights = [math.exp(-transfer_to.decay * (d - dmin) / depth) for d in dists]  # dmin: softmax-stable
    return cand[_weighted_index(rng, weights, sum(weights))]


def mean_root_to_tip(tree) -> float:
    """The tree's mean root-to-tip time — the timescale that makes `Distance` decay scale-free.
    Over the extant tips (all leaves if none survive); 1.0 for a degenerate zero-height tree."""
    root_t = tree.nodes[tree.root].birth_time
    tips = tree.extant_leaves() or tree.leaves()
    depth = sum(n.end_time - root_t for n in tips) / len(tips)
    return depth if depth > 0 else 1.0
