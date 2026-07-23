"""Species trees — the forward birth-death engine.

Per-lineage birth and death grow a tree forward in time and record every
speciation and extinction. ``birth`` and ``death`` are full **rate specs** — a number,
a scope wrapper, or a product — so ``birth = scope.Global(1.0)`` gives one shared
tree-wide budget (linear growth), ``birth = 1.0 * mod.OnTotalDiversity(cap=100)`` slows the
tree as it fills up, ``birth = 1.0 * mod.OnTime({...})`` runs a skyline (the interval-aware
sampler steps to each breakpoint), and ``birth = 1.0 * mod.FromParent(spread=0.2)`` lets
the rate drift down the tree (clade drift): each lineage threads its own inherited factor and
the lineage that speciates or dies is drawn **weighted** by its effective rate.

Still to come: the full result spine, the CLI, and the move to ``zombi2.species``. This
lives here for now so the old package is untouched.
"""

from __future__ import annotations

import functools
import math
import pathlib
import re
from dataclasses import dataclass, field

import numpy as np

from ..rates.modifiers import FromParent, OnTime, OnTotalDiversity
from ..progress import progress_bar
from ..rates.rate import as_rate
from ..rates.scope import Global, PerLineage

#: The rate grammar this level wires (SPEC §5). Both the engine's gate below and the CLI's help read
#: this, so a modifier can never be advertised without being implemented — or silently ignored.
WIRED_SCOPES = (PerLineage, Global)
WIRED_MODIFIERS = (OnTime, OnTotalDiversity, FromParent)


@dataclass
class Node:
    """One lineage segment: born at ``birth_time``, ended at ``end_time`` by a split, a
    death, or reaching the present. A split has two ``children``; a leaf has none."""

    id: int
    parent: int | None
    birth_time: float
    end_time: float | None = None
    children: tuple[int, int] | None = None
    fate: str = "alive"  # alive → "extant" | "extinct" | "unsampled"; internal splits are "speciation"


@dataclass(frozen=True)
class Event:
    """A recorded event in the true history: a speciation (with its two children) or an extinction."""

    time: float
    kind: str  # "speciation" | "extinction"
    node: int
    children: tuple[int, int] | None = None


@dataclass
class Tree:
    """The complete tree: every lineage that ever lived, keyed by id, rooted at ``root``."""

    nodes: dict[int, Node]
    root: int

    def leaves(self) -> list[Node]:
        """Every lineage with no descendants — extant **and** extinct."""
        return [n for n in self.nodes.values() if n.children is None]

    def extant(self) -> list[Node]:
        """The lineages alive at the present."""
        return [n for n in self.nodes.values() if n.fate == "extant"]

    def extinct(self) -> list[Node]:
        """The lineages that died before the present."""
        return [n for n in self.nodes.values() if n.fate == "extinct"]

    def unsampled(self) -> list[Node]:
        """Survivors not observed under incomplete ``sampling`` — kept in the complete tree (told
        apart by their fate) but pruned from the extant tree."""
        return [n for n in self.nodes.values() if n.fate == "unsampled"]

    def to_newick(self) -> str:
        """Serialise to Newick (matching ``tree.to_newick()`` elsewhere in the codebase). Each
        branch length is ``end_time - birth_time`` and every node — leaves and internals — is named
        ``n<id>``.

        The root carries a branch length like any other node: its **stem**, the time from the origin
        to the first split. A forward birth–death run starts from one lineage, so that stem is real
        simulated time in which events happen, and writing ``)n0;`` would silently discard it — for a
        tree whose crown comes late, a large fraction of its history. It is emitted as ``)n0:<stem>;``
        and :func:`read_newick` reads it back."""

        def emit(i: int) -> str:
            node = self.nodes[i]
            bl = node.end_time - node.birth_time
            if node.children is None:
                return f"n{i}:{bl:.6g}"
            inner = ",".join(emit(c) for c in node.children)
            return f"({inner})n{i}:{bl:.6g}"

        root = self.nodes[self.root]
        stem = root.end_time - root.birth_time
        if root.children is None:
            return f"n{self.root}:{stem:.6g};"
        return f"({','.join(emit(c) for c in root.children)})n{self.root}:{stem:.6g};"


_WRITE_OUTPUTS = ("complete", "extant", "events", "fossils", "fates")  # the write vocabulary the CLI reuses


@dataclass
class SpeciesResult:
    """What ``simulate_species_tree`` returns: the ``complete_tree`` (with the dead) and the derived
    ``extant_tree`` (the observed survivors), the ``events`` log (the recorded true history), the
    ``seed``, and any ``fossils``. (The ``record=`` memory dial lands with the data-heavy levels.)"""

    complete_tree: Tree
    events: list[Event]
    seed: int | None
    #: recovered fossils as ``(lineage_id, time)`` pairs, sorted by time — a side output, present
    #: only when ``fossils`` was set; the fossil's lineage is not removed and is not in the extant tree
    fossils: list[tuple[int, float]] = field(default_factory=list)

    @property
    def n_extant(self) -> int:
        """The number of **observed** survivors — the extant tips. Under ``sampling < 1`` this is
        the sampled subset (the rest are ``unsampled``), so it matches the extant tree's tip count."""
        return len(self.complete_tree.extant())

    @functools.cached_property
    def extant_tree(self) -> Tree | None:
        """The survivors' tree — the complete tree pruned to extant lineages with the
        unifurcations suppressed (dated, bifurcating). ``None`` if nothing survived."""
        return prune(self.complete_tree, keep="extant")

    def write(self, directory, outputs=None) -> None:
        """Write outputs to ``directory``, each file prefixed ``species_``; ``outputs`` selects which
        (default = all applicable): ``"complete"`` → ``species_complete.nwk``, ``"extant"`` →
        ``species_extant.nwk`` (if any survived), ``"events"`` → ``species_events.tsv`` (the
        always-recorded true history), ``"fossils"`` → ``species_fossils.tsv`` (if any recovered),
        ``"fates"`` → ``species_fates.tsv`` (each tip's resolved fate).

        ``species_fates.tsv`` is the tip-fate table: one ``lineage<TAB>fate`` row per tip, with fate
        one of ``extant`` / ``extinct`` / ``unsampled``. Fate is resolved once, at the end of the run,
        on the same stable ``n<id>`` that keys every other file, so it never renames anything — it is
        a materialised view of information the run already holds. It exists because the ``.nwk`` records
        only branch lengths, from which a reader cannot tell an extinct tip from a survivor that sits at
        the present; this table says so directly, so a downstream level can build the extant set from
        fate rather than guessing from tip depth."""
        if outputs is None:
            outputs = _WRITE_OUTPUTS
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "complete" in outputs:
            (d / "species_complete.nwk").write_text(self.complete_tree.to_newick() + "\n")
        if "extant" in outputs and self.extant_tree is not None:
            (d / "species_extant.nwk").write_text(self.extant_tree.to_newick() + "\n")
        if "events" in outputs:
            rows = ["time\tkind\tlineage\tchildren"]
            for e in self.events:
                kids = ";".join(f"n{c}" for c in e.children) if e.children else ""
                rows.append(f"{e.time:.6g}\t{e.kind}\tn{e.node}\t{kids}")
            (d / "species_events.tsv").write_text("\n".join(rows) + "\n")
        if "fossils" in outputs and self.fossils:
            rows = ["lineage\ttime"] + [f"n{i}\t{t:.6g}" for i, t in self.fossils]
            (d / "species_fossils.tsv").write_text("\n".join(rows) + "\n")
        if "fates" in outputs:
            # one row per tip (extant / extinct / unsampled); internal nodes are always speciations
            rows = ["lineage\tfate"]
            for n in sorted(self.complete_tree.leaves(), key=lambda nd: nd.id):
                rows.append(f"n{n.id}\t{n.fate}")
            (d / "species_fates.tsv").write_text("\n".join(rows) + "\n")


def prune(tree: Tree, keep: str = "extant") -> Tree | None:
    """Prune the complete tree to a kept set (matching ``prune(tree, keep=...)`` in the codebase):
    drop the pruned subtrees and suppress the unifurcations they leave behind, giving a dated,
    bifurcating tree. Branch lengths merge across suppressed nodes; ``None`` if nothing is kept.

    ``keep="extant"`` (default) keeps the survivors — the extant tree. (``"sampled"``, the
    fossil/serially-sampled tree, arrives with the sampling and fossils slices; it raises for now.)"""
    if keep != "extant":
        raise ValueError(f"keep must be 'extant' until the sampling/fossils slices land, got {keep!r}")
    nodes = tree.nodes
    surviving: dict[int, bool] = {}
    for i in sorted(nodes, reverse=True):  # children have higher ids → processed before parents
        nd = nodes[i]
        surviving[i] = nd.fate == "extant" if nd.children is None else any(surviving[c] for c in nd.children)
    if not any(surviving.values()):
        return None

    def surv_children(i: int) -> list[int]:
        nd = nodes[i]
        return [] if nd.children is None else [c for c in nd.children if surviving[c]]

    # keep the extant leaves and the genuine bifurcations (≥2 surviving children)
    kept = {i for i in nodes
            if (nodes[i].children is None and nodes[i].fate == "extant") or len(surv_children(i)) >= 2}

    new: dict[int, Node] = {}
    ext_root: int | None = None
    for i in kept:
        p = nodes[i].parent  # walk up to the nearest kept ancestor
        while p is not None and p not in kept:
            p = nodes[p].parent
        branch_start = nodes[p].end_time if p is not None else 0.0  # merge the suppressed edges
        new[i] = Node(i, p, branch_start, nodes[i].end_time, None, nodes[i].fate)
        if p is None:
            ext_root = i
    for i in sorted(kept):  # rebuild children from parents, in id order for a stable Newick
        p = new[i].parent
        if p is not None:
            existing = new[p].children
            new[p].children = (i,) if existing is None else existing + (i,)

    return Tree(new, ext_root)


_ZOMBI_LABEL = re.compile(r"^n(\d+)$")  # the id-bearing label to_newick writes on every node


def _assign_fates_from_map(leaves: list[Node], labels: dict[int, str],
                           tip_fates: dict[str, str], *, source: str) -> None:
    """Set each leaf's fate from ``tip_fates`` (``{label: fate}``), joined to the leaves through
    ``labels`` (``{leaf id: label}``). Every leaf must be uniquely labelled and covered, and every
    value must be ``extant`` / ``extinct`` / ``unsampled``; anything off raises, naming ``source`` (the
    map it came from) so the message points at the right file. Used for both a ZOMBI tree keyed by its
    ``n<id>`` labels and an external tree keyed by the user's labels."""
    labelled = [labels.get(n.id) for n in leaves]  # every tip must be uniquely named to map a fate
    if any(lbl is None for lbl in labelled):
        raise ValueError(f"every tip must be named to declare its fate with {source}, but the "
                         "tree has unlabelled tips")
    if len(set(labelled)) != len(labelled):
        dups = sorted({lbl for lbl in labelled if labelled.count(lbl) > 1})
        raise ValueError(f"tip labels must be unique to map fates; repeated: {', '.join(dups)}")
    fates = {k: str(v).lower() for k, v in tip_fates.items()}
    missing = sorted(set(labelled) - set(fates))
    unknown = sorted(set(fates) - set(labelled))
    if missing:
        raise ValueError(f"{source} is missing a fate for: {', '.join(missing)}")
    if unknown:
        raise ValueError(f"{source} names tips not in the tree: {', '.join(unknown)}")
    bad = sorted({v for v in fates.values() if v not in ("extant", "extinct", "unsampled")})
    if bad:
        raise ValueError(f"tip fates must be 'extant', 'extinct' or 'unsampled', got: {', '.join(bad)}")
    for n in leaves:
        n.fate = fates[labels[n.id]]


def _assign_external_fates(leaves: list[Node], names: dict[int, str],
                           tip_fates: dict[str, str] | None, gap: float) -> None:
    """Set the fate on the tips of a **non-ultrametric** external tree from a user-supplied
    ``tip_fates`` (``{tip label: "extant" | "extinct" | "unsampled"}`` — the same vocabulary a species
    run writes to ``species_fates.tsv``). ZOMBI will not infer a tip's fate from its depth here (a
    shallow tip could be extinct *or* an early sample), so a missing or mismatched map raises rather
    than guesses."""
    if tip_fates is None:
        raise ValueError(
            f"input tree is not ultrametric (tip depths differ by {gap:.3g}); ZOMBI can't tell "
            "extinct lineages from early samples — declare each tip's fate with --tip-fates FILE "
            "(a 'tip_name<TAB>extant|extinct|unsampled' row per tip)")
    _assign_fates_from_map(leaves, names, tip_fates, source="--tip-fates")


def read_newick(newick: str, *, tip_fates: dict[str, str] | None = None) -> tuple[Tree, dict[int, str]]:
    """Parse a Newick string into a complete :class:`Tree` and a name-map ``{id: user label}``.

    This is how the CLI loads a species tree back for the downstream levels. Branch lengths are read
    as **durations**: the root sits at time 0 and each node's ``birth_time`` is its parent's
    ``end_time``, so ``end_time - birth_time`` is the parsed length. Two kinds of tree are accepted,
    told apart by the labels:

    - a **ZOMBI complete tree** (every node — internal ones too — is ``n<id>``, as ``to_newick``
      writes it): the ids come from the labels, and the name-map is empty (the labels *are* the ids).
      Fate comes from ``tip_fates`` when given — the run's ``species_fates.tsv``, keyed by the same
      ``n<id>`` — which is authoritative; without it a leaf is ``"extinct"`` if it ends before the
      tree's greatest depth, else ``"extant"`` (a fallback that cannot recover an ``"unsampled"`` tip).
    - any **external tree** (leaves named freely, internal nodes usually unlabelled): fresh ids are
      minted in traversal order (root 0, parents before children), the original labels are returned as
      the **name-map** (``{minted id: user label}``), and fates depend on whether the tree is
      **ultrametric** (all root-to-tip depths equal, within ``1e-6 × height``):

      - **ultrametric** → the tips are contemporaneous, so **every tip is ``"extant"``** (observed);
      - **not ultrametric** → the differing tip depths could mean extinct lineages *or* early
        samples, which ZOMBI cannot tell apart, so it **refuses to guess**: pass ``tip_fates`` — a
        ``{tip label: "extant" | "extinct" | "unsampled"}`` map covering every tip — or a
        :class:`ValueError` is raised. (The CLI fills ``tip_fates`` from ``--tip-fates FILE``, which
        reads the same format a species run writes to ``species_fates.tsv``.)

    A root branch length is read when present — ``to_newick`` writes one, so a ZOMBI tree round-trips
    with its stem intact. External trees usually have none, and then the root gets zero duration and
    the tree starts at its crown, which is all the file says.

    The ``.nwk`` records only branch lengths, so ``"unsampled"`` fate cannot be read from it alone —
    pass ``tip_fates`` (the run's ``species_fates.tsv``) to recover it; without one a survivor reads
    back ``"extant"``, which is still fine for evolving genomes/traits along the tree.

    Only bifurcating trees are supported (an internal node with other than two children raises).
    """
    s = newick.strip().rstrip(";").strip()
    if not s:
        raise ValueError("empty Newick string — is the tree file empty?")
    i = 0

    def skip_ws() -> None:
        # whitespace (incl. the newlines of a line-wrapped file) is insignificant between tokens
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def read_name() -> str:
        # a quoted label ('...' / "...") is taken verbatim (a doubled quote unwrapped to one); an
        # unquoted label runs to the next whitespace or structural char, so stray whitespace never leaks
        nonlocal i
        if i < len(s) and s[i] in "'\"":
            quote = s[i]
            i += 1
            chars: list[str] = []
            while i < len(s):
                if s[i] == quote:
                    if i + 1 < len(s) and s[i + 1] == quote:
                        chars.append(quote)
                        i += 2
                        continue
                    i += 1
                    break
                chars.append(s[i])
                i += 1
            return "".join(chars)
        start = i
        while i < len(s) and s[i] not in ",():;" and not s[i].isspace():
            i += 1
        return s[start:i]

    # parse into a lightweight tree of (name, length, children) so we can decide ids/fates in a
    # second pass, once we know whether every node is ``n<id>``-labelled and the tree's depth
    class _P:
        __slots__ = ("name", "length", "children")

        def __init__(self, name, length, children):
            self.name, self.length, self.children = name, length, children

    def parse() -> _P:
        nonlocal i
        children: list[_P] = []
        skip_ws()
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                children.append(parse())
                skip_ws()
                if i >= len(s):
                    raise ValueError("malformed Newick: unbalanced parentheses — the string ended "
                                     "while a clade was still open (is the tree file truncated?)")
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
                else:
                    raise ValueError(f"malformed Newick: expected ',' or ')' after a clade at "
                                     f"position {i}, got {s[i]!r}")
        skip_ws()
        name = read_name()
        length = 0.0
        skip_ws()
        if i < len(s) and s[i] == ":":
            i += 1
            skip_ws()
            start = i
            while i < len(s) and s[i] not in ",():;" and not s[i].isspace():
                i += 1
            try:
                length = float(s[start:i])
            except ValueError:
                raise ValueError(f"malformed Newick: expected a branch length after ':' at position "
                                 f"{start}, got {s[start:i]!r}") from None
        if children and len(children) != 2:
            raise ValueError(f"only bifurcating trees are supported: node {name or '(unnamed)'!r} has "
                             f"{len(children)} children (collapse polytomies / unifurcations first)")
        return _P(name, length, children)

    root_p = parse()
    skip_ws()
    if i < len(s):
        raise ValueError(f"malformed Newick: unexpected trailing text at position {i}: {s[i:]!r}")

    # ZOMBI complete tree ⟺ every node carries an ``n<id>`` label; then ids come from the labels,
    # otherwise we mint them ourselves and call every tip extant.
    all_labelled = True

    def _scan(p: _P) -> None:
        nonlocal all_labelled
        if not _ZOMBI_LABEL.match(p.name):
            all_labelled = False
        for c in p.children:
            _scan(c)

    _scan(root_p)

    nodes: dict[int, Node] = {}
    names: dict[int, str] = {}  # {minted id: user label} — for external trees; empty for ZOMBI ones
    counter = 0

    def _mint(p: _P) -> int:
        nonlocal counter
        if all_labelled:
            return int(_ZOMBI_LABEL.match(p.name).group(1))
        i_ = counter
        counter += 1
        return i_

    # first pass: assign ids (parents before children) and absolute times from durations
    def _build(p: _P, parent: int | None, birth: float) -> int:
        nid = _mint(p)
        if nid in nodes:
            raise ValueError(f"duplicate node id n{nid} in the Newick (labels must be unique)")
        end = birth + p.length
        child_ids = tuple(_build(c, nid, end) for c in p.children) or None
        nodes[nid] = Node(nid, parent, birth, end, child_ids)  # fate filled in below
        if not all_labelled and p.name:  # external tree: keep the user's label for the name-map
            names[nid] = p.name
        return nid

    root_id = _build(root_p, None, 0.0)

    # second pass: fates. Internal nodes are always speciations; the tips depend on the tree kind.
    for n in nodes.values():
        if n.children is not None:
            n.fate = "speciation"
    leaves = [n for n in nodes.values() if n.children is None]

    if all_labelled:
        if tip_fates is not None:
            # the run's own species_fates.tsv (or a --tip-fates file) states each tip's fate directly,
            # keyed by the same n<id> the tree carries. It is authoritative: depth cannot tell an
            # unsampled survivor from an extant one (both sit at the present), nor an extinct tip that
            # died just before the present, so when it is given we use it rather than guess.
            _assign_fates_from_map(leaves, {n.id: f"n{n.id}" for n in leaves}, tip_fates,
                                   source="tip fates")
            return Tree(nodes, root_id), names
        # no fate table: a tip is extinct if it ends before the present (the greatest end_time). The
        # tolerance is depth-relative — ``to_newick`` prints 6 significant figures, whose rounding
        # accumulates to ~5e-6·height along a root-to-tip path, so a tip at the present can fall a few
        # ×1e-5·height short of the max, far below any real extinction gap. This cannot recover an
        # unsampled tip (it sits at the present, so it reads back extant) — pass the fate table for that.
        present = max(n.end_time for n in nodes.values())
        tol = max(1e-9, 1e-4 * present)
        for n in leaves:
            n.fate = "extinct" if n.end_time < present - tol else "extant"
        return Tree(nodes, root_id), names

    # an external tree: ultrametric ⟺ every tip is contemporaneous (all extant); otherwise the
    # differing depths could be extinctions or early samples, which we refuse to guess (SPEC decision).
    depths = [n.end_time for n in leaves]  # root sits at 0, so a tip's depth is its end_time
    height = max(depths)
    if max(depths) - min(depths) <= max(1e-12, 1e-6 * height):  # ultrametric
        for n in leaves:
            n.fate = "extant"
        return Tree(nodes, root_id), names

    _assign_external_fates(leaves, names, tip_fates, max(depths) - min(depths))
    return Tree(nodes, root_id), names


_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant


def _drift(rate) -> FromParent | None:
    """The :class:`~zombi2.modifiers.FromParent` modifier a rate carries, or ``None``. When present
    the rate is *per-lineage*: the engine threads each lineage's own inherited factor (clade drift)."""
    for m in rate.modifiers:
        if isinstance(m, FromParent):
            return m
    return None


def _weighted_index(rng, weights: list[float], total: float) -> int:
    """Pick an index in proportion to ``weights`` (which must sum to ``total``)."""
    r = rng.random() * total
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r < acc:
            return i
    return len(weights) - 1  # floating-point guard: r == total lands on the last lineage


def _grow(rng, birth_rate, death_rate, n_extant: int | None, total_time: float | None,
          pulses: list[tuple[float, float]], progress: bool = False,
          max_lineages: int | None = None) -> tuple[Tree, list[Event]]:
    """Grow one forward birth-death tree until it reaches ``n_extant`` living lineages,
    reaches ``total_time``, or dies out. Returns the complete tree and the event log.

    When ``birth`` or ``death`` carries an :class:`~zombi2.modifiers.FromParent` modifier the rate
    is *per-lineage*: every lineage threads its own inherited factor (its parent's, nudged at the
    split), so the lineage that speciates or dies is drawn **weighted** by its effective rate rather
    than uniformly. Birth and death drift independently. A rate with no ``FromParent`` keeps a factor
    of 1 and picks uniformly, exactly as before.

    ``pulses`` are scheduled mass extinctions as ``(time, survival)`` pairs sorted by time (time runs
    forward from the crown): at each instant every standing lineage is kept with probability
    ``survival`` and otherwise becomes an extinct leaf. They sit at a point on the timeline, so the
    caller passes them only when ``total_time`` is set."""
    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent: int | None, t: float) -> int:
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    birth_drift = _drift(birth_rate)  # the FromParent modifier on each rate, or None
    death_drift = _drift(death_rate)

    root = new_node(None, 0.0)
    alive = [root]  # a list so picks are reproducible given the seed
    # each lineage's inherited factor, kept in lock-step with `alive` under swap-remove; a rate with
    # no FromParent stays at 1.0, so its total is just scope(base) × modifiers over n, picked uniform
    inh_b = [birth_drift.initial() if birth_drift else 1.0]
    inh_d = [death_drift.initial() if death_drift else 1.0]
    t = 0.0
    events: list[Event] = []
    pulse_idx = 0  # the next unfired mass extinction in `pulses`

    # a tree grows toward whichever stop condition was given: a tip count, or a time
    bar = progress_bar(n_extant if n_extant is not None else total_time, "species",
                       unit="tip" if n_extant is not None else "time", enabled=progress)
    # A run conditioned on time has no natural ceiling: standing diversity grows like
    # exp((birth - death) * t), so a rate a little too high or a time a little too long is the
    # difference between a thousand lineages and ten million. The guard RAISES rather than stopping
    # early — a tree cut off at a size is no longer a sample from the process asked for, and handing
    # one back as if it were would be worse than not running at all.
    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    while alive:
        bar.to(len(alive) if n_extant is not None else t)
        n = len(alive)
        if ceiling is not None and n > ceiling:
            bar.close()
            raise RuntimeError(
                f"the tree passed {ceiling} standing lineages at time {t:.3g} and is still growing "
                f"— birth exceeds death by enough that this run has no realistic end. Lower the "
                f"rates, shorten total_time, cap the growth with a modifier "
                f"(birth * OnTotalDiversity(cap=...)), or raise max_lineages if the size is what "
                f"you want (max_lineages=None removes the guard).")
        # standing diversity = the living lineages; OnTotalDiversity/OnTime read `diversity`/`time`
        ctx = {"diversity": n, "time": t}
        # a drifting rate's total is the sum over lineages of each lineage's effective rate —
        # scope(base) × modifiers evaluated per lineage through its inherited factor (lineages=1
        # is one lineage); a non-drifting rate is scope(base) × modifiers once, over all n lineages
        if birth_drift:
            w_b = [birth_rate.effective(lineages=1, inherited=x, **ctx) for x in inh_b]
            total_birth = sum(w_b)
        else:
            total_birth = birth_rate.effective(lineages=n, **ctx)
        if death_drift:
            w_d = [death_rate.effective(lineages=1, inherited=x, **ctx) for x in inh_d]
            total_death = sum(w_d)
        else:
            total_death = death_rate.effective(lineages=n, **ctx)
        total = total_birth + total_death
        # the total rate is constant until the next skyline breakpoint, mass extinction, or the total_time
        # limit — advance no further than the earliest of them before re-evaluating
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t))
        next_pulse = pulses[pulse_idx][0] if pulse_idx < len(pulses) else math.inf
        horizon = min(next_change, next_pulse)
        if total_time is not None:
            horizon = min(horizon, total_time)

        if total > 0.0:
            t_event = t + float(rng.exponential(1.0 / total))
            if t_event < horizon:  # an event fires before the rate changes
                t = t_event
                if n == n_extant:
                    # already at the target: stop at the time this next event WOULD fire (do not
                    # apply it), so the present sits *after* the last split and the two newest tips
                    # get a real, non-zero branch length.
                    break
                # birth vs death by their totals; then WHICH lineage — weighted by its effective
                # rate if that rate drifts (so faster lineages are likelier), else uniform
                speciates = rng.random() < total_birth / total
                if speciates:
                    i = _weighted_index(rng, w_b, total_birth) if birth_drift else int(rng.integers(n))
                else:
                    i = _weighted_index(rng, w_d, total_death) if death_drift else int(rng.integers(n))
                node = alive[i]
                parent_b, parent_d = inh_b[i], inh_d[i]
                alive[i] = alive[-1]  # swap-remove keeps picks O(1); the inherited factors move in step
                alive.pop()
                inh_b[i] = inh_b[-1]; inh_b.pop()
                inh_d[i] = inh_d[-1]; inh_d.pop()
                if speciates:
                    nodes[node].end_time = t
                    nodes[node].fate = "speciation"
                    c1, c2 = new_node(node, t), new_node(node, t)
                    nodes[node].children = (c1, c2)
                    alive.extend((c1, c2))
                    # each daughter inherits the parent's inherited factor, nudged (1.0 if no drift)
                    inh_b.extend((birth_drift.descend(parent_b, rng), birth_drift.descend(parent_b, rng))
                                 if birth_drift else (1.0, 1.0))
                    inh_d.extend((death_drift.descend(parent_d, rng), death_drift.descend(parent_d, rng))
                                 if death_drift else (1.0, 1.0))
                    events.append(Event(t, "speciation", node, (c1, c2)))
                else:
                    nodes[node].end_time = t
                    nodes[node].fate = "extinct"
                    events.append(Event(t, "extinction", node))
                continue

        # no stochastic event fired before the horizon
        if math.isinf(horizon):
            break  # nothing scheduled and the rate never changes again → nothing more can happen
        if total_time is not None and horizon == total_time:
            t = total_time
            break
        if pulse_idx < len(pulses) and horizon == next_pulse:
            # a mass extinction: each standing lineage is kept with probability `survival`, the rest
            # become extinct leaves at this instant (their inherited factors leave with them)
            t = next_pulse
            survival = pulses[pulse_idx][1]
            pulse_idx += 1
            kept_a: list[int] = []
            kept_b: list[float] = []
            kept_d: list[float] = []
            for k, node_id in enumerate(alive):
                if survival >= 1.0 or rng.random() < survival:
                    kept_a.append(node_id)
                    kept_b.append(inh_b[k])
                    kept_d.append(inh_d[k])
                else:
                    nodes[node_id].end_time = t
                    nodes[node_id].fate = "extinct"
                    events.append(Event(t, "extinction", node_id))
            alive[:] = kept_a
            inh_b[:] = kept_b
            inh_d[:] = kept_d
            continue
        t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    bar.close()
    for i in alive:  # whoever is still alive reached the present
        nodes[i].end_time = t
        nodes[i].fate = "extant"

    return Tree(nodes, root), events


def _mass_extinction_pulses(mass_extinctions, total_time: float | None) -> list[tuple[float, float]]:
    """Turn user ``(time, fraction_lost)`` pulses into the engine's ``(time, survival)`` pairs,
    sorted by time. Time runs **forward from the crown**, so ``(3.0, 0.75)`` = at time 3.0, 75% of
    the standing lineages die (survival 0.25). Empty when none are given. A pulse sits at a point on
    the timeline, so it needs a fixed end — ``total_time`` set — and must fall inside ``(0, total_time)``."""
    if not mass_extinctions:
        return []
    if total_time is None:
        raise ValueError(
            "mass_extinctions need a run with a fixed end — give total_time=..., not n_extant= "
            "(under n_extant= the run can stop before a pulse's time is reached)"
        )
    pulses: list[tuple[float, float]] = []
    for pulse in mass_extinctions:
        time, fraction = pulse
        if (isinstance(time, bool) or not isinstance(time, (int, float))
                or not math.isfinite(time) or not 0.0 < time < total_time):
            raise ValueError(
                f"each mass extinction time must be a number strictly between 0 and total_time ({total_time}), "
                f"got {time!r}"
            )
        if (isinstance(fraction, bool) or not isinstance(fraction, (int, float))
                or not 0.0 <= fraction <= 1.0):
            raise ValueError(f"each mass extinction fraction lost must be in [0, 1], got {fraction!r}")
        pulses.append((time, 1.0 - fraction))
    pulses.sort()
    return pulses


def _apply_sampling(tree: Tree, rho: float, rng) -> None:
    """Incomplete extant sampling: relabel each surviving lineage ``"unsampled"`` with probability
    ``1 - rho`` (in place). The extant tree then prunes to the sampled survivors, while the unsampled
    ones stay in the complete tree, told apart by their fate. ``rho = 1`` observes everyone."""
    if rho >= 1.0:
        return
    for i in sorted(tree.nodes):  # id order + one draw per survivor → reproducible given the seed
        node = tree.nodes[i]
        if node.fate == "extant" and float(rng.random()) >= rho:
            node.fate = "unsampled"


def _recover_fossils(tree: Tree, rate: float, rng) -> list[tuple[int, float]]:
    """Recover fossils along every branch of the complete tree: a branch of length ``L`` yields
    ``Poisson(rate × L)`` fossils, each at a uniform time on the branch. Returns ``(lineage_id,
    time)`` pairs sorted by time. A pure side output — no lineage is removed."""
    if rate <= 0.0:
        return []
    fossils: list[tuple[int, float]] = []
    for i in sorted(tree.nodes):  # id order, then Poisson + uniforms → reproducible given the seed
        node = tree.nodes[i]
        length = node.end_time - node.birth_time
        if length <= 0.0:
            continue
        for _ in range(int(rng.poisson(rate * length))):
            fossils.append((i, node.birth_time + float(rng.random()) * length))
    fossils.sort(key=lambda ft: ft[1])
    return fossils


def simulate_species_tree(birth, death=0.0, *, n_extant=None, total_time=None,
                          mass_extinctions=None, sampling=1.0, fossils=0.0, seed=None,
                          progress=False, max_lineages=100_000) -> SpeciesResult:
    """Grow a forward birth-death tree.

    ``birth`` and ``death`` are rate specs (a number, a ``scope`` wrapper, or a product
    with modifiers); the default scope is **per lineage** (each lineage speciates/dies at
    the base rate, so the tree grows exponentially). Yule = ``death=0``.

    Stop at exactly ``n_extant`` living lineages, **or** at ``total_time`` — give exactly
    one. ``n_extant`` is **conditioned on survival**: a birth-death tree can die out, so we
    restart (advancing the same generator) until one reaches ``n_extant``. ``total_time`` is not
    conditioned. Deterministic given ``seed``.

    ``mass_extinctions`` is a list of ``(time, fraction_lost)`` pulses — e.g. ``[(3.0, 0.75)]`` culls
    75% of the lineages alive at time 3.0 (time runs forward from the crown). It is a point-in-time
    intervention on the process (not a rate) placed on the timeline, so it needs a fixed end:
    give ``total_time`` (not ``n_extant``), with each time strictly inside ``(0, total_time)``.

    ``sampling`` (ρ, default 1.0) is incomplete extant sampling: each survivor is observed with
    probability ρ, the rest relabelled ``unsampled``. It prunes the **extant tree** to the sampled
    survivors (the unsampled ones remain only in the complete tree). ``n_extant`` still stops at that
    many *survivors*; sampling then thins what you observe, so ``result.n_extant`` can be smaller.

    ``fossils`` is a recovery rate along the branches: each branch of length ``L`` yields
    ``Poisson(fossils × L)`` fossils, returned as ``result.fossils`` = ``(lineage, time)`` pairs. A
    **side output** — the fossil's lineage is not removed and does not enter the extant tree.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        # a modifier this engine does not thread would return its default factor of 1.0 — a run that
        # is quietly not the model asked for — so reject it (SPEC §5, the genome engine's discipline)
        if not isinstance(rate.scope, WIRED_SCOPES):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the species engine counts "
                f"lineages — use PerLineage(...) (the default, so a bare number is enough) or "
                f"Global(...) for one shared budget."
            )
        for m in rate.modifiers:
            if not isinstance(m, WIRED_MODIFIERS):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the species engine does not "
                    f"support — OnTime (skyline), OnTotalDiversity (diversity-dependent) and "
                    f"FromParent (clade drift, ClaDS) are wired."
                )
        if _drift(rate) is not None and not isinstance(rate.scope, PerLineage):
            raise ValueError(
                f"{label} carries FromParent (per-lineage drift) but its scope is "
                f"{type(rate.scope).__name__}; a drifting rate must be per lineage — drop the "
                f"scope wrapper (per lineage is the default) or use PerLineage(...)"
            )
    if (n_extant is None) == (total_time is None):
        raise ValueError("give exactly one of n_extant or total_time")
    if n_extant is not None and (isinstance(n_extant, bool) or not isinstance(n_extant, int) or n_extant < 1):
        raise ValueError(f"n_extant must be a positive integer, got {n_extant!r}")
    if total_time is not None and (not isinstance(total_time, (int, float)) or not math.isfinite(total_time) or total_time <= 0):
        raise ValueError(f"total_time must be a positive finite number, got {total_time!r}")
    if isinstance(fossils, bool) or not isinstance(fossils, (int, float)) or not math.isfinite(fossils) or fossils < 0:
        raise ValueError(f"fossils must be a non-negative finite rate, got {fossils!r}")
    if isinstance(sampling, bool) or not isinstance(sampling, (int, float)) or not 0.0 < sampling <= 1.0:
        raise ValueError(f"sampling must be a fraction in (0, 1], got {sampling!r}")
    pulses = _mass_extinction_pulses(mass_extinctions, total_time)  # [] unless mass_extinctions given (needs total_time)

    rng = np.random.default_rng(seed)

    def _finish(tree: Tree, events: list[Event]) -> SpeciesResult:
        # observe (sampling relabels survivors) then recover fossils along the grown branches
        _apply_sampling(tree, sampling, rng)
        return SpeciesResult(tree, events, seed, _recover_fossils(tree, fossils, rng))

    if total_time is not None:
        tree, events = _grow(rng, birth_rate, death_rate, None, total_time, pulses, progress,
                             max_lineages)
        # A time-conditioned run is not conditioned on survival, so with death ≥ birth it can reach
        # total_time with nothing alive. An empty tree is not a sample anyone can use — the extant
        # tree is None and every downstream level would otherwise mistake the last-dying tip for a
        # survivor — so refuse it here rather than hand back a tree with no present.
        if not any(nd.fate == "extant" for nd in tree.nodes.values()):
            raise RuntimeError(
                f"the run went extinct before total_time={total_time:g}: no lineage is alive at the "
                f"present, so there is nothing to grow a genome, sequence or trait along. With death "
                f"close to or above birth, total extinction is likely — lower death, shorten "
                f"total_time, or use n_extant=... (which is conditioned on survival).")
        return _finish(tree, events)

    for _ in range(_MAX_ATTEMPTS):
        tree, events = _grow(rng, birth_rate, death_rate, n_extant, None, [], progress,
                             max_lineages)
        if sum(1 for nd in tree.nodes.values() if nd.fate == "extant") == n_extant:  # survivors (pre-sampling)
            return _finish(tree, events)
    raise RuntimeError(
        f"could not grow a tree to {n_extant} extant lineages in {_MAX_ATTEMPTS} attempts; "
        "birth must comfortably exceed death for large n_extant"
    )


__all__ = ["simulate_species_tree", "SpeciesResult", "Tree", "Node", "Event", "prune"]
