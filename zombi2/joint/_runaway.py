"""The names and the ceiling every joint engine shares.

Three engines grow a tree whose birth rate reads a driver their own growth feeds. Each needs the
same guard, and it raises rather than truncating: a tree cut off at a size is no longer a sample
from the process asked for."""

from __future__ import annotations

#: The live gene-content driver a birth rate reads for a lineage's whole gene count.
GENOME_COUNT = "genomes:count"
#: What a trait driver is called when the run holds one and it was given no name.
#: Naming is what lets a run hold two, so the bare form stays for the common case.
ONE_TRAIT = "trait"


#: What ``max_lineages`` says when it fires. The species engine has had this guard since it grew a
#: tree conditioned on time; a joint run needs it more, not less, because its birth rate can read a
#: driver that its own growth feeds — gene content accumulates, birth rises, more lineages accumulate
#: more gene content — so a rate that looks calm on paper can have no realistic end. It RAISES rather
#: than stopping early, for the reason the species engine gives: a tree cut off at a size is no longer
#: a sample from the process asked for.
def runaway(ceiling: int, t: float) -> RuntimeError:
    return RuntimeError(
        f"the tree passed {ceiling} standing lineages at time {t:.3g} and is still growing — the "
        f"driven birth rate has no realistic end. Lower the rates, shorten total_time, flatten the "
        f"mapping the driver is read through, or raise max_lineages if the size is what you want "
        f"(max_lineages=None removes the guard).")



