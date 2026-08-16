"""The slice contract, in one place.

Three engines slice: a diffusing trait driving speciation, two genes reading each other, and a trait
with a gene's sequence. Each one holds a driver fixed across a stretch of time and releases it at the
boundary, because the driver moves at every instant and so no rate reading it is ever constant.

``step`` is the length of that stretch. It rides on the **connection** rather than on the run, because
it belongs to that reading: a steep response curve needs a finer step than a flat one. It has no
default. Any number this module could invent would be a claim about a timescale only the model knows.

The two rules below were written three times, once per engine, which is how a rule stops being one
rule. `step_of` reads it off the connections and refuses a disagreement; `check_step` refuses a step
the tree is not long enough for.
"""

from __future__ import annotations

import math

from ..params.connection import Driven


def step_of(mods, *, what: str, how: str) -> float:
    """The one ``step`` a set of connections agree on.

    ``mods`` are the `~zombi2.params.connection.Driven` modifiers that read the live driver.
    ``what`` names who is reading — "birth and death", "the trait's switch", "gene 'chaperone'" —
    and ``how`` shows the spelling to fix it with. Both go into the refusals, which is why they are
    arguments: one rule, and each engine still names its own parts."""
    steps, missing = set(), False
    for m in mods:
        if not isinstance(m, Driven):
            continue
        if m.step is None:
            missing = True
        else:
            steps.add(float(m.step))
    if missing:
        raise ValueError(
            f"{what} reads a driver that moves at every instant, so it needs a step= — the stretch "
            f"of time that driver is held fixed across. There is no interval where this rate holds "
            f"still on its own, and so nothing exact to draw against; the run slices instead. Write "
            f"it on the connection, {how}, and pick it so the driver moves little within one slice. "
            f"Then halve it, rerun the same seed, and see whether the answer moves.")
    if len(steps) > 1:
        raise ValueError(
            f"{what} reads one live driver at two resolutions, step={sorted(steps)}. One walk "
            f"carries them and it has one set of slice boundaries, so the readings have to agree.")
    if not steps:
        raise ValueError(f"{what} reads no live driver, so there is no step to take one from.")
    return next(iter(steps))


def check_step(step, tallest: float) -> float:
    """``step`` as a number, checked against the tree it will be walked over."""
    if isinstance(step, bool) or not isinstance(step, (int, float)) or not math.isfinite(step) \
            or step <= 0:
        raise ValueError(
            f"step is the stretch of time a driver is held fixed across, in the tree's own units, "
            f"so it must be finite and positive; got {step!r}.")
    if step >= tallest:
        raise ValueError(
            f"step={step:g} is not shorter than the tree itself ({tallest:.3g}), so the driver would "
            f"be read once, at the start, and the run would be one conditioned run in each "
            f"direction. Pick a step the driver moves little within.")
    return float(step)
