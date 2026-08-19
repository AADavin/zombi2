# zombi2.joint

Two levels that cannot be run in order, grown in **one** run. A pair is *joint* when neither level
can be finished first — a trait that drives speciation, or gene content that does. When the driver
*can* be grown first, the run is **conditioned** instead: two ordinary runs in order, the finished
driver handed to the second as an object (in Python) or as its written log (across two commands).
See [Dependent runs](../guide/conditioning.md) and [Joint runs](../guide/joining.md).

A level joined to **itself** stays on that level's own function, with `joint=True`:
[`simulate_genomes_family`][zombi2.genomes.simulate_genomes_family] for a gene family driving the
rest of its genome, and [`simulate_traits`][zombi2.traits.simulate_traits] for two traits reading
each other.

::: zombi2.joint.simulate

::: zombi2.joint.JointResult
