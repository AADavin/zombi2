"""Drawing one of several things in proportion to its weight — the one sampler every engine needs.

Whenever a rate differs across the living lineages — because a driver varies it, or because a
per-family draw weights it — the total is summed *with* those weights, so whatever the event lands
on has to be drawn with them too. Summing one way and drawing another would make the rate say one
thing and the run do another.

It lived three times over, once in `zombi2.species`, once in `zombi2.genomes._live` and once in
`zombi2.joint`, and the genome copy was imported across packages from the species one. Here it is
plumbing rather than domain code, which is what this package is for, and it can be imported from
anywhere without a cycle: `zombi2.genomes` imports `zombi2.species`, so the reverse could not.
"""

from __future__ import annotations

__all__ = ["WeightedIndex", "weighted_index"]


def weighted_index(rng, weights, total: float) -> int:
    """The index of one of ``weights``, drawn in proportion to it. ``total`` is their sum, passed in
    because the caller has already computed it to race the events."""
    r = float(rng.random()) * total
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r < acc:
            return i
    return len(weights) - 1     # floating-point guard: r == total lands on the last one


class WeightedIndex:
    """One weight per candidate, with their **total** and a **draw** from them, neither of which
    walks the candidates.

    `weighted_index` reads every weight to draw one, and its caller has already added every weight up
    to race the events. Where the candidates are the living lineages and an event changes one weight,
    that pair is the whole of a run's growth with the tree: every step pays for every lineage.

    The weights are the leaves of a binary tree, each node the sum of its two children. Setting one
    weight rewrites the path from its leaf to the root — ``log2(candidates)`` additions — the total
    *is* the root, and a draw walks down from the root choosing the side the draw falls in.

    Two things follow, and the second is a change a user can see:

    - Each node is **recomputed** from its children rather than moved by a difference, so the total is
      a sum of the weights as they stand. It cannot drift from them, however many changes it has
      taken, and the same seed reads the same number.
    - A tree adds in pairs where a list adds left to right, so the total is not bit-for-bit the number
      ``sum()`` gives, and neither is the boundary a draw is compared against. Runs stay the same
      process and the same draw count — each `pick` takes one ``rng.random()``, as `weighted_index`
      does — but event times differ in their last digit, and from there a long run diverges.

    The candidates are the live lineage set, so they enter at the end (`append`) and leave by
    swap-remove (`remove`), the two moves `zombi2.genomes._live` makes.
    """

    __slots__ = ("_tree", "_cap", "_n")

    def __init__(self) -> None:
        self._n = 0
        self._cap = 1
        self._tree = [0.0, 0.0]

    def __len__(self) -> int:
        return self._n

    def __getitem__(self, k: int) -> float:
        if not 0 <= k < self._n:
            raise IndexError(k)
        return self._tree[self._cap + k]

    def __iter__(self):
        return iter(self._tree[self._cap:self._cap + self._n])

    @property
    def total(self) -> float:
        """The weights summed, as they stand — the root of the tree."""
        return self._tree[1]

    def set(self, k: int, value: float) -> None:
        """Give candidate ``k`` this weight, and mend the path from it to the root."""
        tree = self._tree
        i = self._cap + k
        tree[i] = value
        i >>= 1
        while i:
            tree[i] = tree[2 * i] + tree[2 * i + 1]
            i >>= 1

    def append(self, value: float = 0.0) -> None:
        """A candidate enters at the end, holding this weight."""
        if self._n == self._cap:
            self._double()
        self._n += 1
        self.set(self._n - 1, value)

    def remove(self, k: int) -> None:
        """Retire candidate ``k``: the last one takes its place, mirroring `zombi2.genomes._live`'s
        swap-remove, and the tail empties out."""
        last = self._n - 1
        if k != last:
            self.set(k, self._tree[self._cap + last])
        self.set(last, 0.0)
        self._n -= 1

    def pick(self, rng) -> int:
        """One candidate drawn in proportion to its weight, walking down from the root. Takes the
        same single ``rng.random()`` draw `weighted_index` takes, so the two are interchangeable in
        the random stream."""
        r = float(rng.random()) * self._tree[1]
        tree, cap = self._tree, self._cap
        i = 1
        while i < cap:
            i *= 2
            left = tree[i]
            if r >= left:
                r -= left
                i += 1
        k = i - cap
        return k if k < self._n else self._n - 1   # the r == total guard, as `weighted_index` has

    def _double(self) -> None:
        """Twice the room, the leaves carried over and the tree above them laid again."""
        cap = self._cap * 2
        tree = [0.0] * (2 * cap)
        tree[cap:cap + self._n] = self._tree[self._cap:self._cap + self._n]
        for i in range(cap - 1, 0, -1):
            tree[i] = tree[2 * i] + tree[2 * i + 1]
        self._cap, self._tree = cap, tree
