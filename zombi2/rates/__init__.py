"""The cross-level rate grammar (SPEC §5): ``effective rate = scope(base) × modifiers``.

Shared by every level, so it lives in one place. Reach the pieces as submodules::

    from zombi2.rates import scope, modifiers
    birth = scope.Global(1.0)
    birth = 1.0 * modifiers.Time({0: 1.0, 3: 0.3})

- ``scope`` — *per what?* ``PerCopy`` · ``PerLineage`` · ``PerSite`` · ``Global``
- ``modifiers`` — *depends on what?* ``Time`` · ``Diversity`` · ``Inherited``
- ``rate`` — the internal ``Rate`` plumbing (users never build a ``Rate`` directly)
- ``distributions`` — value/length distributions
"""
