"""Which substitution model each lineage evolves under (SPEC §5, §9 — experimental).

A rate says how *fast* a lineage evolves; the model says what the change *looks like* — which
residues turn into which, and what composition the sequence settles at. Until now a run had one
model for every branch, so a clade could evolve faster but not differently, and the whole class of
compositional questions was out of reach: an endosymbiont clade drifting toward AT, a lineage with a
different GC bias, the compositional attraction that misleads a tree-builder in a way long-branch
attraction does not.

    model = Models().set_by(Clade({"endo": ["n7", "n9"]}),
                            {"endo": hky85(frequencies=(0.4, 0.1, 0.1, 0.4)),
                             "rest": hky85()})

A model is the fourth thing a driver can target, beside a rate, an extent and a choice — and the odd
one out, because the first three take a **factor** and are multiplied by it while a model takes none
and is **selected**. That is why the verb is ``set_by``: the driver supplies the whole thing rather
than a multiple of it.

Python only, no CLI flag, exactly as `substitution_models.reversible` is: a model is a K×K matrix,
which SPEC §5 already says has no written form a flag can carry.

What may vary along the tree is the matrix. What may not is the **alphabet** — the sequences of one
gene copy are one string, so there is no half-DNA half-protein run — and the **across-site rate
classes**, because a site's class is drawn once for the family and holds all the way down the tree,
so two clades cannot sort the same site into different classes.
"""

from __future__ import annotations

from ..params.driver import Clade
from ..params.parameter import RateCompositionError
from .substitution_models import SubstitutionModel

__all__ = ["Models"]

_A_MODEL = ("take one from the menu (jc69(), hky85(kappa=2.0), lg(), …) or build your own with "
            "substitution_models.reversible(S, frequencies)")


class Models:
    """A per-lineage substitution model set — see the module docstring.

    Written from its own entry point, with one legal verb, exactly as ``Recipients()`` is for
    ``transfer_to``: there is no base model to scale, so there is no scope in front of it and
    ``set_by`` is the only thing it takes."""

    def __init__(self, driver=None, models=None, per_species=None) -> None:
        self.driver = driver
        self.models = models
        self.per_species = per_species

    # --- the one verb ---------------------------------------------------------------------------

    def set_by(self, driver, mapping=None, *, step=None) -> "Models":
        """Give each named clade its own model. ``mapping`` is ``{clade label: model}``, and every
        label the driver paints — including ``"rest"``, the lineages in no named clade — needs one."""
        if self.models is not None:
            raise RateCompositionError(
                "Models() carries one set_by: the driver supplies the whole model, so a second "
                "would be a second answer to the same question.")
        if step is not None:
            raise ValueError(
                "step is the resolution a CONTINUOUS driver is read at; a clade never switches "
                "along a branch, so there is nothing to cut finer.")
        if not isinstance(driver, Clade):
            raise ValueError(
                "a per-lineage model is read off the tree today: "
                "Models().set_by(Clade({...}), {...}). A trait switches partway along a branch, and "
                "this engine samples ONE transition matrix per branch endpoint pair — a model that "
                "changed mid-branch would need the branch cut at the switch and sampled stretch by "
                "stretch, which is not implemented. A trait may still drive the substitution RATE: "
                "substitution=PerSite(0.05).scaled_by(trait, {...}).")
        if not isinstance(mapping, dict) or not mapping:
            raise ValueError(f"Models().set_by needs a non-empty {{clade label: model}} dict, got "
                             f"{mapping!r}")
        for label, m in mapping.items():
            if not isinstance(label, str) or not label.strip():
                raise ValueError(f"a clade label must be a non-empty string, got {label!r}")
            if not isinstance(m, SubstitutionModel):
                raise ValueError(f"the model for {label!r} is not a substitution model — {_A_MODEL}, "
                                 f"got {m!r}")
        first, *rest = mapping.items()
        for label, m in rest:
            if m.alphabet != first[1].alphabet:
                raise ValueError(
                    f"{first[0]!r} is {first[1].name} over {first[1].alphabet!r} and {label!r} is "
                    f"{m.name} over {m.alphabet!r}. One gene copy's sequence is one string, so a "
                    f"run has one alphabet: a clade may evolve under a different matrix, not over "
                    f"different residues.")
            if m.site_rates != first[1].site_rates or m.site_shares != first[1].site_shares:
                raise ValueError(
                    f"{first[0]!r} is {first[1].name} and {label!r} is {m.name}, whose across-site "
                    f"rate classes differ. A site's class is drawn ONCE for the family and holds "
                    f"down the whole tree, so two clades cannot sort the same site into different "
                    f"classes — give every model the same across_sites(...) call, or none.")
        return Models(driver, dict(mapping))

    def _no(self, verb: str):
        raise RateCompositionError(
            f"a model set takes set_by and nothing else: {verb} contributes a factor, and a model is "
            f"not a quantity to multiply — it is the thing itself. Write "
            f"Models().set_by(Clade({{...}}), {{...}}).")

    def scaled_by(self, *_a, **_k):
        self._no("scaled_by")

    def weighted_by(self, *_a, **_k):
        self._no("weighted_by")

    def varying_among(self, *_a, **_k):
        self._no("varying_among")

    # --- what the level reads -------------------------------------------------------------------

    @property
    def alphabet(self) -> str:
        """The one alphabet every model in the set shares — checked in `set_by`."""
        return next(iter(self.models.values())).alphabet

    @property
    def name(self) -> str:
        """For a message, not for parsing: ``'endo:HKY85 | rest:HKY85'``."""
        return " | ".join(f"{k}:{m.name}" for k, m in self.models.items())

    def resolve(self, tree) -> "Models":
        """Paint the tree and return the set with ``per_species`` filled — ``{node id: model}``.

        Painted with `resolve_groups`, the one function the engine paints clade membership with, so
        this cannot disagree with a rate scoped to the same clade."""
        from ..genomes._transfer import resolve_groups

        painted = resolve_groups(tree, self.driver.groups)
        labels = set(painted.values())
        missing = sorted(labels - set(self.models))
        if missing:
            raise ValueError(
                f"no model for {missing}: every lineage evolves under one, and a lineage in no named "
                f"clade is in 'rest'. Name a model for each of {sorted(labels)}.")
        stray = sorted(set(self.models) - labels)
        if stray:
            raise ValueError(
                f"{stray} name(s) no lineage of this tree, so those models would evolve nothing. "
                f"The tree's groups are {sorted(labels)}.")
        return Models(self.driver, self.models,
                      {i: self.models[painted[i]] for i in tree.nodes})

    def at(self, species: int) -> SubstitutionModel:
        """The model on species branch ``species`` — after `resolve`."""
        return self.per_species[species]

    def __repr__(self) -> str:
        # a description, not a written form: a model is a K×K matrix, which SPEC §5 says has no
        # spelling a flag can carry, so this does not round-trip and does not claim to.
        inner = ", ".join(f"{k!r}: {m.name}" for k, m in (self.models or {}).items())
        return f"Models().set_by({self.driver!r}, {{{inner}}})"
