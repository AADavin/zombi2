"""Parameters for the backward species-tree simulation.

Follows the msprime idiom: an immutable, validated model object handed to a stateless
simulator.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SpeciesTreeModel:
    """Constant-rate birth-death parameters, conditioned on the number of extant tips.

    Parameters
    ----------
    birth, death:
        Speciation rate ``λ`` (> 0) and extinction rate ``μ`` (>= 0).
    n_tips:
        Number of extant species ``N`` to condition on (>= 2).
    age:
        Tree age. Interpreted per ``age_type``. v1 requires an explicit age
        (conditioning on ``N`` alone, sampling the age, is a future addition).
    age_type:
        ``"crown"`` (default) — ``age`` is the age of the root/MRCA; the reconstructed
        tree's root sits at time 0 and extant leaves at ``age``.
        ``"stem"`` — ``age`` is the time of origin (a stem lineage precedes the crown).
    """

    birth: float
    death: float
    n_tips: int
    age: float | None = None
    age_type: str = "crown"

    def validate(self) -> None:
        if self.birth <= 0:
            raise ValueError(f"birth rate must be > 0, got {self.birth}")
        if self.death < 0:
            raise ValueError(f"death rate must be >= 0, got {self.death}")
        if self.n_tips < 2:
            raise ValueError(f"n_tips must be >= 2, got {self.n_tips}")
        if self.age_type not in ("crown", "stem"):
            raise ValueError(f"age_type must be 'crown' or 'stem', got {self.age_type!r}")
        if self.age is None:
            raise NotImplementedError(
                "v1 requires an explicit `age`; conditioning on n_tips alone "
                "(sampling the age) is not yet implemented."
            )
        if self.age <= 0:
            raise ValueError(f"age must be > 0, got {self.age}")
