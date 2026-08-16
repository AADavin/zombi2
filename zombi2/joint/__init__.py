"""Joint models — one run simulates two levels at once (SPEC §2–4).

A run is **joint** when neither level can be finished before the other starts, so one run has to
produce both. That is the whole of it, and it says nothing on its own about the species tree: two
levels can drive each other on a tree handed to the run. In the two models built here the species
tree **is** one of the two being simulated, so it comes out of the run rather than going into it:

- a **discrete trait** drives speciation (BiSSE / MuSSE), ``P(Species, Traits)`` — birth/death read the
  trait state on each lineage while the trait evolves by its own Mk process on the growing tree;
- **gene content** drives speciation, ``P(Species, Genomes)`` — birth/death read a summary of each
  lineage's live genome (its total gene count, or the presence of a named family) while the genome
  evolves by duplication/loss/origination on the growing tree;
- a **continuous trait** drives speciation (QuaSSE) — birth/death read a diffusing value on each
  lineage while it diffuses on the growing tree.

One Gillespie races the event classes over the living lineages at once: **speciation** and
**extinction** (per lineage, driver-read), plus the driver's own events — a **trait switch** (the CTMC
out-rate) or a genome **duplication/loss/origination**. A driver event changes a lineage's state
without touching the topology; a speciation hands the parent's driver state (its trait, its genome) to
both daughters. Because a discrete driver only changes at events, the rate is piecewise-constant
between them and the race is **exact** — no thinning. A diffusion is the exception: it moves at every
instant, so that run **slices**, holding the value fixed across a ``step`` the driven rate declares.

The mechanism is the same ``scaled_by`` as conditioning; only the ``driver`` differs — here a
**live level name** (``"trait"``, ``"genomes:count"``, ``"genomes:<family>"``) rather than a filename.
Driving *both* birth and death recovers full state-dependent diversification (BiSSE's λ and μ)."""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass
from typing import Any


from .._runtime.summary import write_summary
from ..params.mapping import check_not_a_kernel
from ..rng import stream
from ..params.driver import OnTime, OnTotalDiversity
from ..params.evaluate import DRAWN, INHERITED, describe, is_implemented
from ..params.connection import Driven
from ..genomes import FamilyGenomesResult, FamilyGenome
from ..species import SpeciesResult
from ..tree import Tree
from ..traits import ContinuousTrait, DiscreteTrait, TraitsResult

from .._runtime.slicing import step_of
from ..params.parameter import as_rate
from ..params.scope import PerLineage
from . import _genomes_traits, _species_continuous, _species_genomes, _species_traits

#: The rate grammar a joint run supports on ``birth`` / ``death`` (SPEC §5). Declared, like every
#: other level, so the gate below cannot fall behind what the engine threads: the loop passes ``time``
#: and ``diversity`` into every rate and steps its Gillespie at each ``next_change``, so the two
#: covariates are as real here as at the species level, and ``scaled_by`` is what makes the run joint
#: at all. What is missing is missing on purpose — see the rejections in `_simulate_joint()`.
IMPLEMENTED_MODIFIERS = (OnTime, OnTotalDiversity, Driven)

#: `JointResult.write`'s vocabulary. The tokens are the two **levels** and the run's own summary, not
#: their files: a joint run's whole claim is that it writes each level exactly as that level's own
#: command does, so each is written with its own default and there is nothing here to restate.
_WRITE_OUTPUTS = ("summary", "species", "driver")

_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant
from ._runaway import GENOME_COUNT as _GENOME_COUNT      # noqa: E402  (the shared names)
from ._runaway import ONE_TRAIT as _ONE_TRAIT            # noqa: E402


@dataclass
class JointResult:
    """What `simulate()` returns — **both** simulated levels of a joint run. ``species`` is the
    grown tree (a `SpeciesResult`: ``complete_tree``, ``extant_tree``, the
    speciation/extinction ``events``); the **driver** level that grew with it is either ``trait`` (a
    `TraitsResult`, for a trait→speciation run) or ``genome`` (a
    `FamilyGenomesResult`, for a gene-content→speciation run) — exactly one is set.
    The tree is an output, grown by the driver it carries, so the levels share one ``complete_tree``."""

    species: SpeciesResult
    seed: int | None
    trait: TraitsResult | None = None
    genome: FamilyGenomesResult | None = None
    #: the sequence level, when it is one of the two — a `~zombi2.sequences.SequencesResult` holding
    #: the one gene the run evolved. Set only by the trait-and-sequence model, whose tree and gene
    #: trees both come in. Typed loosely because importing that result here would be a cycle: the
    #: sequence level imports the genome level, which this module imports too.
    sequences: "Any | None" = None

    def __repr__(self) -> str:
        grown = [n for n, v in (("trait", self.trait), ("genome", self.genome),
                                ("sequences", self.sequences)) if v is not None]
        return (f"JointResult({self.n_extant} extant tips, {' and '.join(grown)}, "
                f"seed={self.seed})")

    @property
    def complete_tree(self) -> Tree:
        return self.species.complete_tree

    @property
    def extant_tree(self):
        return self.species.extant_tree

    @property
    def n_extant(self) -> int:
        return self.species.n_extant

    @property
    def events(self) -> list:
        """The species events (speciation / extinction). The driver level's own events are
        ``trait.events`` / ``genome.events``."""
        return self.species.events

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``joint_summary.json``.

        A joint run grew two levels at once, so this holds both of their summaries under one roof
        rather than inventing a third vocabulary: ``species`` is the tree that came out, and exactly
        one of ``trait`` / ``genome`` is the driver that shaped it. The tree is an *output* here, which
        is the whole point of the command, so its realised birth and death rates are the numbers worth
        reading — they are what the driver did."""
        out = {"level": "joint", "seed": self.seed,
               "driver": "trait" if self.trait is not None else "genome",
               "species": self.species.summary()}
        if self.sequences is not None:
            out["sequences"] = self.sequences.summary()
        if self.trait is not None:
            out["trait"] = self.trait.summary()
        if self.genome is not None:
            out["genome"] = self.genome.summary()
        return out

    def write(self, directory, outputs=_WRITE_OUTPUTS, *, flat: bool = False) -> None:
        """Write both levels to ``directory`` (created if needed), each exactly as its own command
        writes it: ``"species"`` → the `SpeciesResult` files (``species_complete.nwk`` /
        ``species_extant.nwk`` / ``species_events.tsv`` / ``species_fates.tsv`` /
        ``species_summary.json``); ``"driver"`` → the level that grew with it, a trait's
        ``trait_values.tsv`` / ``trait_events.tsv`` / ``trait_tree.nwk`` / ``trait_summary.json`` or
        a genome's ``genome_events.tsv`` / ``profiles.tsv`` / ``genomes.tsv`` /
        ``initial_genome.tsv`` / ``gene_trees/`` / ``genome_summary.json``; ``"summary"`` →
        ``joint_summary.json``, the one file that is the joint run's own.

        The tokens are the two **levels**, not their files, because each is written with that level's
        own default — which is what makes a joint run's directory the two runs it stands in for. Pick
        files *within* a level through the level itself: ``result.species.write(d, outputs=…)``,
        ``result.trait.write(d, outputs=…)``. ``flat`` is passed to the driver level, the only one of
        the two with a many-files-per-run output.

        Both levels land in the one directory named here; ``zombi2 joint`` groups them under
        ``species/`` and ``traits/`` / ``genomes/`` instead, and writes the same files."""
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "summary" in outputs:
            write_summary(d / "joint_summary.json", self.summary())
        if "species" in outputs:
            self.species.write(d)
        if "driver" in outputs:
            # two independent tests, not an if/else: a result carrying neither writes neither, rather
            # than reaching for `.write` on None
            if self.trait is not None:
                self.trait.write(d)
            if self.genome is not None:
                self.genome.write(d, flat=flat)
            if self.sequences is not None:
                self.sequences.write(d, flat=flat)


def _simulate_joint(*, birth, death=0.0, trait=None, genome=None, n_extant=None,
                    total_time=None, seed=None, max_lineages=100_000) -> JointResult:
    """Grow a tree **and** the driver that drives its speciation, in one run (SPEC §2–4).

    ``birth`` and ``death`` are rate specs (per lineage). Make either read the driver with
    ``.scaled_by(driver, mapping)``. Driving speciation has its own function because the driver
    cannot be grown first: it would have to grow on the tree it drives, so the two grow together.
    The driver is a **live level name** (``"trait"``, ``"genomes:count"``); it has to agree with the
    spec given below, and a filename is refused, because a driver read from a file is conditioning.
    Give **exactly one** driver:

    - ``trait = traits.discrete(...)`` — a discrete trait drives speciation (BiSSE / MuSSE), read as
      ``.scaled_by("trait", {"small": 1.0, "large": 2.0})``. Driving both birth and death gives
      state-dependent λ *and* μ.
    - ``trait = traits.continuous(...)`` — a **diffusing** trait drives speciation (QuaSSE), read
      through a `~zombi2.params.Curve` or a `~zombi2.params.Scalar` rather than a table, and with a
      ``step=`` on the connection: a diffusion moves at every instant, so this one run slices.
    - ``genome = genomes.genome(...)`` — **gene content** drives speciation (``P(Species, Genomes)``),
      read as the total gene count ``.scaled_by("genomes:count", curve)`` or the presence of a named
      family ``.scaled_by("genomes:toxin", {"present": 2.0, "absent": 1.0})`` (declare it with
      ``families=[family("toxin")]``).

    The engine behind `simulate()`, which is the way to call it::

        joint.simulate(
            species.birth_death(
                birth = PerLineage(1.0).scaled_by("genomes:toxin", {"present": 3.0, "absent": 1.0}),
                n_extant = 100),
            genomes.genome(origination=0.2, loss=0.1, families=[family("toxin")]), seed = 1)

    ``max_lineages`` (default 100000) stops a run that has no realistic end. A joint birth rate reads
    a driver the run itself grows, so it can feed itself — gene content accumulates, birth rises, and
    more lineages accumulate more gene content — and a rate that looks calm on paper then grows
    without bound. It **raises** rather than truncating, for the reason the species engine gives: a
    tree cut off at a size is no longer a sample from the process asked for. ``max_lineages=None``
    removes the guard.

    The driver is an **unexecuted** process spec, grown with the tree. Stop at exactly ``n_extant``
    living lineages (conditioned on survival — a birth-death tree can die out, so it restarts,
    advancing the same generator) **or** at ``total_time`` — give exactly one. Returns a
    `JointResult` carrying the grown tree and the driver level (``.trait`` or ``.genome``).
    Deterministic given ``seed``. Clade drift (an inherited value) combined with driving, and gene
    transfer in a joint run, are not available.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    if (trait is None) == (genome is None):
        raise TypeError(
            "give exactly one driver: trait=traits.discrete(...) OR genome=genomes.genome(...)."
        )
    # collect the Driven driver names on birth/death (a joint model's diversification must be per lineage)
    driver_names: list[str] = []
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        # `Rate.scope` holds the scope **class**, so this is an identity test rather than an
        # isinstance: a scope instance never exists (SPEC §5).
        if rate.scope is not PerLineage:
            assert rate.scope is not None      # `as_rate` above fills in the default scope
            raise ValueError(
                f"{label} has a {rate.scope.__name__} scope, but a joint diversification rate is "
                f"per lineage — write PerLineage(...) (the default, so a bare number is enough)."
            )
        for m in rate.modifiers:
            if m.reads == (DRAWN, "families"):
                # not a missing feature: there is nothing here for it to mean (the species level
                # says the same thing about the same modifier, for the same reason)
                raise ValueError(
                    f"{label} carries a value drawn among families, but a diversification rate has "
                    f"no gene families — varying_among('families', ...) belongs on a genomes rate. "
                    f"To make speciation depend on gene content, drive it: "
                    f"birth = PerLineage(1.0).scaled_by(\"genomes:count\", ...)."
                )
            if m.reads == (INHERITED, "lineages"):
                raise ValueError(
                    f"{label} carries a value inherited among lineages (clade drift); drift and a "
                    f"driven rate are not available together — use one or the other."
                )
            if m.reads == (DRAWN, "lineages"):
                # The species level takes this; the joint engine does not thread it, so accepting it
                # here would run the model without the rate variation the user asked for — and the
                # same `--birth` expression working on `zombi2 species` makes that a trap rather than
                # merely a gap (SPEC §5: reject, never silently ignore).
                raise ValueError(
                    f"{label} carries a value drawn among lineages (independent per-lineage rates); "
                    f"per-lineage rate variation and a driven rate are not available together in a "
                    f"joint run — use one or the other. On its own, "
                    f"varying_among('lineages', ...) works at the species level."
                )
            if not is_implemented(m, IMPLEMENTED_MODIFIERS, "joint"):
                # the backstop: anything this engine does not thread would come back as its default
                # factor of 1.0, which is a run quietly not the model that was asked for (SPEC §5).
                # Declared rather than enumerated here, so a modifier added later cannot slip through.
                raise ValueError(
                    f"{label} carries {describe(m)}, which a joint run does not support. It "
                    f"takes changing_at (skyline), scaled_by(TotalDiversity(cap=...)) "
                    f"(diversity-dependent) and scaled_by (the driver that makes the run joint)."
                )
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
                if not isinstance(m.driver, str):
                    raise TypeError(
                        f"{label} is driven by a {type(m.driver).__name__} object, but a joint model "
                        f"drives from a live level *name* (a string, e.g. \"trait\" / \"genomes:count\"). "
                        + ("A clade is read off a finished tree, and a joint run grows the tree as it "
                           "goes, so there is no clade to read yet."
                           if type(m.driver).__name__ == "Clade" else
                           "A grown result object is conditioning — pass it to the driven level's run.")
                    )
                driver_names.append(m.driver)
    if not driver_names:
        raise ValueError(
            "a joint model needs the driver to drive something: give birth (or death) a "
            "scaled_by(...). With neither driven, grow the two levels as independent runs instead."
        )
    # the driver spec must match the driver names
    if trait is not None:
        if not isinstance(trait, (DiscreteTrait, ContinuousTrait)):
            raise TypeError(
                "trait= must be traits.discrete(states=[...], switch=...) or "
                "traits.continuous(rate=...) — a trait process spec.")
        # a run holds one trait, so the bare "trait" always names it; a name additionally lets a rate
        # say which, which is what a run holding two will need
        trait_keys = {_ONE_TRAIT} | ({f"traits:{trait.name}"} if trait.name else set())
        bad = sorted({s for s in driver_names if s not in trait_keys})
        if bad:
            named = " or ".join(f'"{k}"' for k in sorted(trait_keys))
            raise ValueError(
                f"with a trait participant, drive from that trait — scaled_by({named}, ...); got "
                f"driver(s) {bad}. (A filename driver is conditioning, not a joint run.)"
            )
        if isinstance(trait, ContinuousTrait):
            step = step_of([m for r in (birth_rate, death_rate) for m in r.modifiers
                            if isinstance(m, Driven) and m.driver in trait_keys],
                           what="birth and death",
                           how='scaled_by("trait", Curve(f), step=0.05)')
    else:
        if not isinstance(genome, FamilyGenome):
            raise TypeError("genome= must be genomes.genome(...) — a family-genome process spec.")
        if genome.transfer:
            raise ValueError(
                "transfer is not available while the tree is being simulated: a transfer needs the "
                "set of lineages alive at that instant, and on a growing tree that set is still "
                "forming. It works on a tree handed to the run — genomes with traits.")
        for s in driver_names:
            if s == _GENOME_COUNT:
                continue
            if s.startswith("genomes:"):
                name = s.split(":", 1)[1]
                if name not in genome.family_names:
                    raise ValueError(
                        f'scaled_by("{s}", ...) names family {name!r}, but genomes.family was not '
                        f"declared with it — add families=[…, family({name!r})]."
                    )
                continue
            raise ValueError(
                f'with genome=, drive from gene content — "genomes:count" or "genomes:<family>"; '
                f"got {s!r}."
            )
    if (n_extant is None) == (total_time is None):
        raise ValueError("give exactly one of n_extant or total_time")
    if n_extant is not None and (isinstance(n_extant, bool) or not isinstance(n_extant, int) or n_extant < 1):
        raise ValueError(f"n_extant must be a positive integer, got {n_extant!r}")
    if total_time is not None and (not isinstance(total_time, (int, float))
                                   or not math.isfinite(total_time) or total_time <= 0):
        raise ValueError(f"total_time must be a positive finite number, got {total_time!r}")

    rng, seed = stream("joint", seed)       # own stream, and a drawn seed if none was given
    unique_driver_names = sorted(set(driver_names))

    def grow_once(target_n, tt) -> tuple[Tree, JointResult]:
        if isinstance(trait, ContinuousTrait):
            g = _species_continuous.grow(rng, birth_rate, death_rate, trait, step, target_n, tt,
                                       max_lineages, tuple(sorted(trait_keys)))
            g.trait_events.sort(key=lambda c: c.time)
            result = JointResult(SpeciesResult(g.tree, g.species_events, seed, []), seed,
                                 trait=TraitsResult(g.tree, g.trait_values, g.trait_events, seed))
        elif trait is not None:
            g = _species_traits.grow(rng, birth_rate, death_rate, trait, target_n, tt, max_lineages,
                            tuple(sorted(trait_keys)))
            g.trait_events.sort(key=lambda c: c.time)
            result = JointResult(SpeciesResult(g.tree, g.species_events, seed, []), seed,
                                 trait=TraitsResult(g.tree, g.trait_values, g.trait_events, seed,
                                                    kind="discrete"))
        else:
            g = _species_genomes.grow(rng, birth_rate, death_rate, genome, unique_driver_names,
                                   target_n, tt, max_lineages)
            result = JointResult(SpeciesResult(g.tree, g.species_events, seed, []), seed,
                                 genome=FamilyGenomesResult(g.tree, g.genomes, g.genome_events,
                                                            seed, g.genome_names, {}))
        return g.tree, result

    if total_time is not None:
        return grow_once(None, total_time)[1]

    for _ in range(_MAX_ATTEMPTS):
        tree, result = grow_once(n_extant, None)
        if sum(1 for nd in tree.nodes.values() if nd.fate == "extant") == n_extant:
            return result
    raise RuntimeError(
        f"could not grow a tree to {n_extant} extant lineages in {_MAX_ATTEMPTS} attempts; "
        "birth must comfortably exceed death for large n_extant"
    )


def _classify(participants):
    """Sort the things handed to `simulate()` by which level each belongs to.

    A participant is a **process spec** — `~zombi2.species.BirthDeath`,
    `~zombi2.traits.DiscreteTrait`, `~zombi2.genomes.FamilyGenome` — never a finished result. A
    finished result is a driver you already have, which is conditioning, and it belongs in the driven
    level's own run."""
    from ..species import BirthDeath

    kinds = {"species": [], "traits": [], "genomes": [], "sequences": []}
    for p in participants:
        if isinstance(p, BirthDeath):
            kinds["species"].append(p)
        elif isinstance(p, (DiscreteTrait, ContinuousTrait)):
            kinds["traits"].append(p)
        elif isinstance(p, FamilyGenome):
            kinds["genomes"].append(p)
        elif type(p).__name__ == "GeneSpec":
            kinds["sequences"].append(p)
        else:
            raise TypeError(
                f"joint.simulate takes process specs — species.birth_death(...), "
                f"traits.discrete(...), traits.continuous(...), genomes.genome(...), "
                f"sequences.gene(...) — and got {p!r}. A finished result is a "
                f"driver you already have, which is conditioning: pass it to the driven level's own "
                f"run instead.")
    return kinds


def _on_a_given_tree(kinds, *, tree, seed) -> JointResult:
    """The joint models whose tree is an **input** — today, a genome and a discrete trait.

    ``tree`` is required here for the reason it is refused when the species level is a participant:
    a run simulates what it is given specs for, and takes everything else."""
    from ..tree import as_tree

    n_traits, n_genomes = len(kinds["traits"]), len(kinds["genomes"])
    if not (n_traits == 1 and n_genomes == 1):
        raise NotImplementedError(
            f"a joint run on a tree you supply is built for one genome and one trait; got "
            f"{n_genomes} genome(s) and {n_traits} trait(s). A level joined to **itself** does not "
            f"come here at all — that is one level and one result, so it stays on that level's own "
            f"function with joint=True.")
    if tree is None:
        raise ValueError(
            "neither participant simulates the species tree, so this run needs one: pass tree=. "
            "Give what you are not simulating.")
    tree = as_tree(tree, level="joint")
    genome, trait = kinds["genomes"][0], kinds["traits"][0]
    if not isinstance(trait, DiscreteTrait):
        raise NotImplementedError(
            "a genome and a CONTINUOUS trait driving each other on a given tree is not built. The "
            "genome's race would have to be sliced against the diffusion, as a continuous trait "
            "driving speciation already is. Use traits.discrete(...) here, or grow one level first "
            "and condition the other on it.")

    # each level must actually read the other, or these are two independent runs wearing one call
    trait_keys = {_ONE_TRAIT} | ({f"traits:{trait.name}"} if trait.name else set())
    reads_trait, reads_genome = False, False
    for rate in (genome.duplication, genome.loss, genome.origination, genome.transfer):
        r = as_rate(rate, default_scope=PerLineage)
        for m in r.modifiers:
            if isinstance(m, Driven) and m.driver in trait_keys:
                reads_trait = True
    from ..traits.discrete import _switch_specs
    for spec in _switch_specs(trait.switch):
        if not isinstance(spec, (int, float)):
            for m in as_rate(spec, default_scope=PerLineage).modifiers:
                if isinstance(m, Driven) and str(m.driver).startswith("genomes:"):
                    reads_genome = True
    if not (reads_trait or reads_genome):
        raise ValueError(
            'neither level reads the other, so this is two independent runs rather than one joint '
            'one. Give a genome rate a scaled_by("trait", ...), or the trait\'s switch a '
            'scaled_by("genomes:<family>", ...) — or run the two levels separately.')

    rng, seed = stream("joint", seed)
    g = _genomes_traits.grow(rng, tree, genome, trait, tuple(sorted(trait_keys)), seed)
    g.trait_events.sort(key=lambda c: c.time)
    return JointResult(
        # the tree came in rather than out, so its own event log is not this run's to write
        SpeciesResult(tree, [], seed, []), seed,
        trait=TraitsResult(tree, g.trait_values, g.trait_events, seed, kind="discrete"),
        genome=FamilyGenomesResult(tree, g.genomes, g.genome_events, seed, g.genome_names, {},
                                   g.genome_initial, genome.max_family_size))


def _genomes_and_sequences(kinds, *, tree, genomes, seed, record=False) -> JointResult:
    """A genome and a gene's sequence, each driving the other — the last cell of the map.

    Both levels are participants, so both come out of the run. The **tree** is the one thing handed
    over, which is the rule every joint model follows: give what you are not simulating. The gene
    trees are not handed over either — the genome participant produces them.

    See ``docs/design/genomes-sequences.md`` for why this is the family resolution only, and for what
    a genome event does to a sequence."""
    from .._runtime.slicing import check_step, step_of
    from ..genomes import FamilyGenomesResult
    from ..genomes.gene_trees import gene_trees_from_edges
    from ..params.parameter import as_rate
    from ..params.scope import PerCopy
    from ..sequences import GeneSpec, SequencesResult, _gene_newick, _split
    from ..sequences._loop import scaled_tree
    from ..sequences._record import recorder_for
    from ..sequences.substitution_models import SubstitutionModel, decode
    from ..tree import as_tree
    from . import _genomes_sequences

    if genomes is not None:
        raise ValueError(
            "genomes= hands over a FINISHED genome run, and here the genome is one of the things "
            "being simulated. Pass the tree instead: joint.simulate(genomes.genome(...), "
            "sequences.gene(...), tree=ct).")
    if len(kinds["genomes"]) != 1 or len(kinds["sequences"]) != 1:
        raise NotImplementedError(
            f"this pair is built for one genome and one gene; got {len(kinds['genomes'])} genome(s) "
            f"and {len(kinds['sequences'])} gene(s).")
    if tree is None:
        raise ValueError(
            "neither participant simulates the species tree, so this run needs one: pass tree=. "
            "Give what you are not simulating.")
    tree = as_tree(tree, level="joint")
    genome, spec = kinds["genomes"][0], kinds["sequences"][0]
    if not isinstance(spec, GeneSpec):
        raise TypeError(f"the sequence participant is sequences.gene(...), got {spec!r}.")
    if not isinstance(spec.model, SubstitutionModel):
        raise TypeError(f"gene {spec.name!r} needs one substitution model — model=jc69() — and got "
                        f"{spec.model!r}.")
    if spec.offers is None:
        raise ValueError(
            f"the genome reads this gene, so the gene has to say what it publishes: "
            f"offers=sequences.composition('GC', absent=0.5) on sequences.gene({spec.name!r}, ...).")
    stray = sorted(set(spec.offers.letters) - set(spec.model.alphabet))
    if stray:
        raise ValueError(
            f"gene {spec.name!r} offers composition({spec.offers.letters!r}), which names {stray} — "
            f"not in this model's alphabet ({spec.model.alphabet}).")
    if spec.name not in genome.family_names:
        raise ValueError(
            f"gene {spec.name!r} names no family this genome spec declared. Declare it there — "
            f"genomes.genome(..., families=[family({spec.name!r})]).")

    gene_name = f"sequences:{spec.name}"
    connections, lookup = [], set()
    for label, rate in (("duplication", genome.duplication), ("loss", genome.loss),
                        ("transfer", genome.transfer), ("origination", genome.origination)):
        for m in as_rate(rate, default_scope=PerCopy).modifiers:
            if not isinstance(m, Driven):
                continue
            if m.driver != gene_name:
                raise ValueError(
                    f"{label} reads {m.driver!r}. Here a genome rate reads the gene in this same "
                    f'run: scaled_by("{gene_name}", Curve(f), step=0.05).')
            connections.append(m)
            lookup.add(m.key)
    if not connections:
        raise ValueError(
            "neither level reads the other, so this is two independent runs wearing one call. Give "
            f'a genome rate a scaled_by("{gene_name}", ...), or run the two levels in order — '
            "a genome run, then simulate_sequences over it.")
    tallest = max(n.end_time for n in tree.nodes.values())
    step = check_step(step_of(connections, what="the genome's rates",
                              how=f'scaled_by("{gene_name}", Curve(f), step=0.05)'), tallest)

    rng, seed = stream("joint", seed)
    events, recorder = recorder_for(record, tree.labels())
    g = _genomes_sequences.grow(rng, tree, genome, spec, tuple(lookup), step, record=recorder)
    events.sort(key=lambda e: e.time)

    labels = tree.labels()
    fam = g.genome_names[spec.name]
    trees = gene_trees_from_edges(g.genome_events, tree)
    by_copy, founding_states, length_by_copy = g.sequences[spec.name]
    gt = trees[fam]
    states = {id(n): by_copy[n.copy] for n in _walk_nodes(gt.complete)}
    lengths = {id(n): length_by_copy[n.copy] for n in _walk_nodes(gt.complete)}
    aln, anc = _split(gt, states, labels, spec.model)
    scaled = scaled_tree(gt, lengths)
    ext = scaled.extant
    seqs = SequencesResult(
        {fam: aln}, {fam: anc}, {fam: decode(founding_states, spec.model.alphabet)},
        {fam: {"complete": _gene_newick(scaled.complete, labels),
               "extant": _gene_newick(ext, labels) if ext is not None else None}},
        {"complete": None, "extant": None}, seed, {}, {}, "family", spec.model.alphabet,
        tuple(labels[i] for i in sorted(tree.extant_leaves())), (spec.name,), events)
    return JointResult(
        SpeciesResult(tree, [], seed, []), seed,
        genome=FamilyGenomesResult(tree, g.genomes, g.genome_events, seed, g.genome_names, {},
                                   g.genome_initial, genome.max_family_size),
        sequences=seqs)


def _walk_nodes(root):
    """Every node of a gene tree, in any order."""
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return out


def _traits_and_sequences(kinds, *, tree, genomes, seed, record=False) -> JointResult:
    """A trait and a gene's sequence, each driving the other — the cross-level join whose two ends
    are the furthest apart (design note §7).

    A sequence lives on a gene tree, which lives on the species tree, so this run needs a **genome
    run** handed over rather than a bare tree: that is where the gene trees and the tree they sit on
    both come from. Everything else follows the same rule as the other joint models — you give what
    you are not simulating."""
    from ..genomes import FamilyGenomesResult, OrderedGenomesResult
    from ..params.parameter import as_rate
    from ..params.scope import PerSite
    from ..sequences import SequencesResult, _gene_newick, _split
    from .._runtime.slicing import check_step, step_of
    from ..sequences._loop import scaled_tree
    from ..sequences._record import recorder_for
    from ..sequences.substitution_models import SubstitutionModel, decode
    from ..traits.discrete import _driven_entries, _switch_specs
    from . import _traits_sequences

    n_traits, n_genes = len(kinds["traits"]), len(kinds["sequences"])
    if kinds["genomes"] and not kinds["species"] and not kinds["traits"]:
        return _genomes_and_sequences(kinds, tree=tree, genomes=genomes, seed=seed, record=record)
    if kinds["species"] or kinds["genomes"]:
        raise ValueError(
            "a trait and a sequence drive each other on a tree the run is handed, so the species "
            "tree and the genome are inputs rather than participants: hand the genome run over with "
            "genomes=g, which carries the tree its gene trees sit on.")
    if not (n_traits == 1 and n_genes == 1):
        raise NotImplementedError(
            f"this pair is built for one trait and one gene; got {n_traits} trait(s) and "
            f"{n_genes} gene(s).")
    if not isinstance(genomes, (FamilyGenomesResult, OrderedGenomesResult)):
        raise ValueError(
            "a sequence needs the gene trees to evolve along, so this run takes the genome run that "
            "produced them: joint.simulate(traits.discrete(...), sequences.gene(...), genomes=g). "
            f"Got {type(genomes).__name__}." + (
                " A bare tree is not enough — a sequence lives on a gene tree, and the genome run is "
                "what has both." if tree is not None else ""))
    if tree is not None:
        raise ValueError(
            "the tree comes with the genome run here — its gene trees sit on it — so passing both "
            "leaves two answers to one question. Give genomes=g alone.")
    trait, spec = kinds["traits"][0], kinds["sequences"][0]
    if not isinstance(trait, DiscreteTrait):
        raise NotImplementedError(
            "a CONTINUOUS trait and a sequence driving each other is not built: the trait's own walk "
            "would have to be sliced against the sequence's, which the discrete one is not. Use "
            "traits.discrete(...) here.")
    if not isinstance(spec.model, SubstitutionModel):
        raise TypeError(f"gene {spec.name!r} needs one substitution model — model=lg() — and got "
                        f"{spec.model!r}.")
    declared = getattr(genomes, "family_names", {}) or {}
    if spec.name not in declared:
        raise ValueError(
            f"gene {spec.name!r} names no family of this genome run, which declared "
            f"{sorted(declared) or 'none'}. Declare it there — "
            f"simulate_genomes_family(..., families=[family({spec.name!r})]).")
    if spec.offers is None:
        raise ValueError(
            f"the trait reads this gene, so the gene has to say what it publishes: "
            f"offers=sequences.composition('KR', absent=0.08) on sequences.gene({spec.name!r}, ...).")
    stray = sorted(set(spec.offers.letters) - set(spec.model.alphabet))
    if stray:
        raise ValueError(
            f"gene {spec.name!r} offers composition({spec.offers.letters!r}), which names {stray} — "
            f"not in this model's alphabet ({spec.model.alphabet}).")

    # each side's connection: the trait reads "sequences:<name>", the gene reads the trait by name
    trait_keys = {_ONE_TRAIT} | ({f"traits:{trait.name}"} if trait.name else set())
    gene_name = f"sequences:{spec.name}"
    connections: list = []            # every reading, for the one place the step rule lives
    reads = 0
    gene_lookup, trait_lookup = set(), set()
    for sw in _switch_specs(trait.switch):
        if isinstance(sw, (int, float)):
            continue
        for m in as_rate(sw, default_scope=PerLineage).modifiers:
            if not isinstance(m, Driven):
                continue
            if m.driver != gene_name:
                raise ValueError(
                    f"the trait's switch rate reads {m.driver!r}. Here it reads the gene in this "
                    f'same run: scaled_by("{gene_name}", Curve(...), step=0.05).')
            gene_lookup.add(m.key)
            connections.append(m)
            reads += 1
    rate = as_rate(1.0 if spec.substitution is None else spec.substitution, default_scope=PerSite)
    if rate.scope is not PerSite or rate.base is None:
        raise ValueError(
            f"gene {spec.name!r}'s substitution rate is read per site and needs a base — write "
            f"PerSite(1.0).scaled_by(...).")
    factors = []
    for m in rate.modifiers:
        if not isinstance(m, Driven) or m.driver not in trait_keys:
            raise ValueError(
                f"gene {spec.name!r}'s substitution rate carries {describe(m)}; here it reads the "
                f'trait of this same run: scaled_by("{sorted(trait_keys)[0]}", {{...}}).')
        trait_lookup.add(m.key)
        factors.append((m.key, m.mapping))
        reads += 1
    if not reads:
        raise ValueError(
            "neither level reads the other, so this is two independent runs wearing one call. Give "
            f'the trait\'s switch a scaled_by("{gene_name}", ...), or the gene\'s substitution rate '
            f'a scaled_by("trait", ...).')
    tree = genomes.complete_tree
    gene_tree = genomes.gene_trees[declared[spec.name]]
    tallest = max(n.end_time for n in tree.nodes.values())
    step = check_step(step_of(connections, what="the trait's switch",
                              how=f'scaled_by("{gene_name}", Curve(f), step=0.05)'), tallest)
    rng, seed = stream("joint", seed)
    # `DiscreteTrait._resolve` settles the generator into one constant matrix, which is exactly what
    # a switch rate reading a composition cannot be. So the alphabet, the root state and the split
    # shift are taken from the spec here — its own checks along with them — and the generator is left
    # as rate specs, rebuilt per lineage inside the walk.
    states = list(trait.states)
    idx = {s: i for i, s in enumerate(states)}
    if trait.start is None:
        start_i = int(rng.integers(len(states)))
    elif trait.start in idx:
        start_i = idx[trait.start]
    else:
        raise ValueError(f"start must be one of states={states} (or None for a uniform draw), "
                         f"got {trait.start!r}")
    at_split = trait.at_speciation
    if at_split is not None and (isinstance(at_split, bool)
                                 or not isinstance(at_split, (int, float))
                                 or not 0.0 <= at_split <= 1.0):
        raise ValueError(f"at_speciation must be a probability in [0, 1] (the shift chance), "
                         f"got {at_split!r}")
    shift = 0.0 if at_split is None else float(at_split)
    founder = (spec.model if spec.start is None else spec.start)
    founding = rng.choice(spec.model.k, size=spec.length,
                          p=founder.stationary).astype("int8")
    events, recorder = recorder_for(record, tree.labels())
    grown = _traits_sequences.grow(
        rng, tree, gene_name=spec.name, gene_tree=gene_tree, model=spec.model, length=spec.length,
        founder=founding,
        letters=spec.offers.letters, absent=spec.offers.absent, gene_keys=tuple(gene_lookup),
        base_rate=float(rate.base), gene_factors=tuple(factors),
        trait_states=states, trait_entries=_driven_entries(list(states), trait.switch),
        trait_start=start_i, trait_shift=shift, trait_keys=tuple(trait_lookup), step=step,
        record=recorder)
    events.sort(key=lambda e: e.time)

    labels = tree.labels()
    fam = declared[spec.name]
    states, founding_states, length_of = grown.sequences[spec.name]
    aln, anc = _split(gene_tree, states, labels, spec.model)
    scaled = scaled_tree(gene_tree, length_of)
    ext = scaled.extant
    seqs = SequencesResult(
        {fam: aln}, {fam: anc}, {fam: decode(founding_states, spec.model.alphabet)},
        {fam: {"complete": _gene_newick(scaled.complete, labels),
               "extant": _gene_newick(ext, labels) if ext is not None else None}},
        {"complete": None, "extant": None}, seed, {}, {}, "family", spec.model.alphabet,
        tuple(labels[i] for i in sorted(tree.extant_leaves())), (spec.name,), events)
    return JointResult(
        SpeciesResult(tree, [], seed, []), seed,
        trait=TraitsResult(tree, grown.trait_values, grown.trait_events, seed, kind="discrete"),
        sequences=seqs)


def simulate(*participants, tree=None, genomes=None, seed=None, record: bool = False,
             max_lineages=100_000) -> JointResult:
    """Simulate two levels **at once**, because neither can be finished before the other starts
    (SPEC §2–4).

    Each participant is a process spec, and you **give what you are not simulating**::

        # the tree is one of the two, so it comes out of the run
        joint.simulate(species.birth_death(birth=faster_if_large, death=0.2, n_extant=100),
                       traits.discrete(name="size", states=["small", "large"], switch=0.1), seed=1)

        joint.simulate(species.birth_death(birth=faster_with_toxin, n_extant=100),
                       genomes.genome(origination=0.2, loss=0.1, families=[family("toxin")]), seed=1)

    A rate reads the other participant by **name**, ``"<level>:<handle>"`` — ``"traits:size"`` for a
    named trait, ``"genomes:toxin"`` for a declared family, ``"genomes:count"`` for a lineage's whole
    gene count. A run holding one unnamed trait also answers to ``"trait"``.

    `~zombi2.traits.continuous` is a participant too, and it is the one driver that does not race
    exactly: a diffusion moves at every instant, so the run holds it fixed across a ``step=`` written
    on the connection and releases it at each boundary.

    The species tree is an output exactly when `~zombi2.species.birth_death` is one of the
    participants; otherwise ``tree`` supplies it. A level driving **itself** does not come here at
    all — that is one level and one result, so it stays on that level's own function with
    ``joint=True``.

    Returns a `JointResult`. Deterministic given ``seed``.
    """
    kinds = _classify(participants)
    n_species, n_traits, n_genomes = (len(kinds[k]) for k in ("species", "traits", "genomes"))
    if kinds["sequences"]:
        return _traits_and_sequences(kinds, tree=tree, genomes=genomes, seed=seed,
                                     record=record)
    if record:
        raise ValueError(
            "record= keeps the SEQUENCE level's own history, and only a run holding that level has "
            "one to keep. The trait and genome levels record theirs always, in trait.events and "
            "genome.events.")
    if genomes is not None:
        raise ValueError(
            "genomes= hands over a finished genome run, and only a sequence participant needs one — "
            "a sequence lives on the gene trees it produced. Drop it, or add "
            "sequences.gene(name=..., model=..., length=...).")
    if n_species == 0:
        return _on_a_given_tree(kinds, tree=tree, seed=seed)
    if tree is not None:
        raise ValueError(
            "the species tree is one of the things this run simulates, so it comes out rather than "
            "going in: drop tree=, or drop species.birth_death(...) and hand the tree over.")
    if n_species > 1:
        raise ValueError("give one species.birth_death(...) — a run grows one tree.")
    if n_traits + n_genomes != 1:
        raise ValueError(
            "give exactly one level for the tree to be simulated with: traits.discrete(...) or "
            f"genomes.genome(...). Got {n_traits} trait(s) and {n_genomes} genome(s).")
    spec = kinds["species"][0]
    driver = kinds["traits"][0] if n_traits else kinds["genomes"][0]
    return _simulate_joint(birth=spec.birth, death=spec.death,
                          n_extant=spec.n_extant, total_time=spec.total_time,
                          seed=seed, max_lineages=max_lineages,
                          **({"trait": driver} if n_traits else {"genome": driver}))


def simulate_joint(**_):
    """Retired. A joint run is written as its **participants** now (SPEC §2–4)::

        joint.simulate(species.birth_death(birth=…, death=…, n_extant=100),
                       traits.discrete(name="size", states=[…], switch=0.1), seed=1)

    The rates and the stop condition move onto `~zombi2.species.birth_death`, which makes the tree
    one of the things being simulated rather than a keyword of the run; and the driver is a
    participant beside it rather than a ``trait=`` / ``genome=`` slot. That is what lets one function
    take every joint model instead of one per pair.
    """
    raise TypeError(
        "simulate_joint is no longer written — a joint run is its participants: "
        "joint.simulate(species.birth_death(birth=…, death=…, n_extant=100), "
        "traits.discrete(name='size', states=[…], switch=0.1), seed=1). The rates and the stop "
        "condition go on species.birth_death, and the driver is a participant beside it.")


__all__ = ["simulate", "JointResult"]
