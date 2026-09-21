# zombi2.traits

Level 4: a trait evolving on the species tree, alongside the genome chain rather than inside it.
Two functions, because a continuous trait and a discrete one take genuinely different arguments. A
third, `simulate_traits`, is for several traits that depend on each other — the case with no
order to grow them in. Two discrete traits are one chain over the product of their states, walked
exactly. A pair holding a continuous trait has no exact solution, so that run holds the traits
still over short steps (`zombi2.traits.stepped`).

::: zombi2.traits.simulate_continuous

::: zombi2.traits.simulate_discrete

::: zombi2.traits.simulate_traits

::: zombi2.traits.TraitsResult

::: zombi2.traits.Change

::: zombi2.traits.DiscreteTrait

::: zombi2.traits.ContinuousTrait
