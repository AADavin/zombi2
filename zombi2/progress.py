"""A progress bar for a long run, shown only when a caller asks for one.

A simulation of any size is a wait with nothing on screen, so every level takes ``progress=True`` and
reports how far along it is. The bar is off by default: a library function that printed to the
terminal on its own would be a nuisance inside a loop, a notebook, or another program. The ``zombi2``
commands turn it on, since a person is watching.

What a level counts differs — a species tree grows toward a tip count or a time, a genome run works
through the species tree's schedule, a sequence run goes family by family — so each passes its own
``total`` and unit, and the shape is the same for all of them::

    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    try:
        while ...:
            bar.to(si)
    finally:
        bar.close()

``tqdm`` is an optional dependency. Without it the bar is silently the no-op below, so a run behaves
the same either way and only the display is missing.
"""

from __future__ import annotations


class _Silent:
    """The bar when nobody is watching: every call is a no-op, so a level's loop needs no branch."""

    def to(self, value: float) -> None:
        pass

    def update(self, delta: float = 1) -> None:
        pass

    def close(self) -> None:
        pass


class _Tqdm:
    """A ``tqdm`` bar behind the same calls. ``to`` is an absolute position, which is what a level
    that can go backwards needs — extinctions take a growing tree's tip count down again."""

    def __init__(self, bar) -> None:
        self._bar = bar

    def to(self, value: float) -> None:
        self._bar.n = min(value, self._bar.total) if self._bar.total else value
        # update(0) redraws on tqdm's own schedule; refresh() would redraw on every call, which for
        # a Gillespie loop is a screenful per event and megabytes down a pipe
        self._bar.update(0)

    def update(self, delta: float = 1) -> None:
        self._bar.update(delta)

    def close(self) -> None:
        self._bar.close()


def progress_bar(total: float | None, desc: str, *, unit: str = "it", enabled: bool = False):
    """A progress bar over ``total`` units, or a silent stand-in when ``enabled`` is false or ``tqdm``
    is not installed. Close it when the loop ends — ``try/finally``, so an exception mid-run does not
    leave a half-drawn line behind."""
    if not enabled:
        return _Silent()
    try:
        from tqdm import tqdm
    except ImportError:                                   # optional dependency, never required
        return _Silent()
    # disable=None silences the bar when the output is not a terminal, so a redirected run writes
    # its summary and nothing else
    return _Tqdm(tqdm(total=total, desc=desc, unit=unit, leave=False, disable=None,
                      bar_format="{desc:>9} {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} {unit} "
                                 "[{elapsed}<{remaining}]"))


def track(items, desc: str, *, unit: str = "it", enabled: bool = False):
    """Iterate ``items`` behind a bar — the shape for a loop whose length is known up front, where
    :func:`progress_bar` would be three lines of bookkeeping around a ``for``."""
    if not enabled:
        return items
    bar = progress_bar(len(items), desc, unit=unit, enabled=True)

    def walk():
        try:
            for item in items:
                yield item
                bar.update()
        finally:
            bar.close()

    return walk()


__all__ = ["progress_bar", "track"]
