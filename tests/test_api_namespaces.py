"""scikit-learn-style submodule namespaces (``zombi2.species``, ...).

These are thin, *additive* re-export modules layered over the implementation
modules. This suite pins the contract:

* each namespace imports and exposes every name it advertises in ``__all__``;
* a namespaced name is the *same object* as the top-level one
  (``zombi2.species.BirthDeath is zombi2.BirthDeath``);
* ``import zombi2`` still exposes every original public name; and
* the namespaces partition the top-level public API (every public name lives
  in exactly one namespace, and no namespace advertises a name absent from the
  top level).
"""

from __future__ import annotations

import importlib

import pytest

import zombi2 as z

# The scikit-learn-style namespaces and the exact public names each advertises
# (kept in sync, by hand, with the module ``__all__`` lists so a drift on either
# side is caught).
NAMESPACES = {
    "species": [
        "Tree", "TreeNode", "read_newick", "prune",
        "BirthDeath", "Yule", "EpisodicBirthDeath", "ClaDS", "DiversityDependent",
        "CladeShiftBirthDeath", "simulate_species_tree", "add_ghost_lineages",
    ],
    "genomes": [
        "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
        "NucleotideGenome", "Segment", "simulate_nucleotide_genomes",
        "NucleotideResult", "Block", "read_gff", "GffGenome",
        "RateModel", "SharedRates", "PerGenomeRates", "FamilySampledRates",
        "BranchRates", "EventWeight", "TransferModel", "ConversionModel",
        "GenomeSimulator", "GenomeResult", "ProfileMatrix",
        "simulate_genomes", "Genomes", "GenomeTrace", "read_events_trace",
        "build_gene_trees", "run_replicates",
    ],
    "traits": [
        "BrownianMotion", "OrnsteinUhlenbeck", "MultivariateBrownian",
        "MultivariateOU", "MultiOptimumOU", "ThresholdModel", "EarlyBurst", "Mk",
        "CorrelatedBinary", "CorrelatedBinaryK", "HiddenStateMk", "simulate_traits",
        "replicate_traits",
        "TraitResult", "pagel_lambda", "pagel_delta", "pagel_kappa",
        "DEC", "simulate_biogeography", "Cladogenesis",
    ],
    "sequences": [
        # relaxed molecular clocks (folded into sequences — the rate half of sequence evolution)
        "Clock", "RateScaledTree", "StrictClock", "UncorrelatedLogNormalClock",
        "UncorrelatedGammaClock", "WhiteNoiseClock", "AutocorrelatedLogNormalClock",
        "CIRClock", "RateVariation",
        "SequenceEvolution", "GenePhylograms",
        "SubstitutionModel", "GammaRates", "jc69", "k80", "hky85", "gtr",
        "poisson", "lg", "wag", "jtt", "dayhoff", "make_model", "is_protein_model",
        "DNA_MODELS", "PROTEIN_MODELS", "AMINO_ACIDS",
        "evolve_on_tree", "read_fasta", "write_fasta",
    ],
    "coevolve": [
        "TraitGeneCoupling", "TraitTrajectory", "TraitLinkedRates",
        "TraitLinkedResult", "simulate_trait_linked_genomes",
        "GeneDiversification", "GeneDiversificationResult",
        "simulate_gene_diversification", "simulate_co_diversification",
        "CladogeneticGenome", "CladogeneticGenomeResult",
        "simulate_cladogenetic_genome",
        "GeneConditionedTrait", "GeneConditionedTraitResult",
        "simulate_gene_conditioned_trait",
        "TraitGeneFeedback", "TraitGeneFeedbackResult", "simulate_trait_gene_feedback",
        "BiSSE", "MuSSE", "HiSSE", "QuaSSE", "simulate_sse",
    ],
    "distributions": [
        "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform",
        "as_distribution",
    ],
}

# Low-level primitives deliberately left at the top level only (no namespace).
LOW_LEVEL = {
    "__version__",
    "EventType", "GeneOp", "EventRecord", "Selection", "Region", "TargetParams",
    "rust_available",
}


@pytest.mark.parametrize("ns", sorted(NAMESPACES))
def test_namespace_imports_and_all_matches(ns):
    """The namespace imports and its ``__all__`` is exactly the expected set."""
    mod = importlib.import_module(f"zombi2.{ns}")
    assert set(mod.__all__) == set(NAMESPACES[ns])


@pytest.mark.parametrize(
    "ns,name",
    [(ns, n) for ns, names in NAMESPACES.items() for n in names],
)
def test_namespace_exposes_name(ns, name):
    """(a) every advertised name is importable as an attribute of the namespace."""
    mod = importlib.import_module(f"zombi2.{ns}")
    assert hasattr(mod, name), f"zombi2.{ns} is missing {name}"


@pytest.mark.parametrize(
    "ns,name",
    [(ns, n) for ns, names in NAMESPACES.items() for n in names],
)
def test_namespace_object_identity(ns, name):
    """(b) the namespaced object *is* the top-level object (same identity)."""
    mod = importlib.import_module(f"zombi2.{ns}")
    assert getattr(mod, name) is getattr(z, name), (
        f"zombi2.{ns}.{name} is not zombi2.{name}"
    )


def test_from_import_style_works():
    """The scikit-learn ``from zombi2.<ns> import <Name>`` style resolves."""
    from zombi2.species import DiversityDependent, simulate_species_tree
    from zombi2.traits import OrnsteinUhlenbeck
    from zombi2.sequences import lg, StrictClock, UncorrelatedLogNormalClock
    from zombi2.genomes import simulate_genomes
    from zombi2.coevolve import simulate_trait_linked_genomes
    from zombi2.distributions import Gamma

    assert DiversityDependent is z.DiversityDependent
    assert simulate_species_tree is z.simulate_species_tree
    assert OrnsteinUhlenbeck is z.OrnsteinUhlenbeck
    assert lg is z.lg
    assert StrictClock is z.StrictClock
    assert UncorrelatedLogNormalClock is z.UncorrelatedLogNormalClock
    assert simulate_genomes is z.simulate_genomes
    assert simulate_trait_linked_genomes is z.simulate_trait_linked_genomes
    assert Gamma is z.Gamma


def test_top_level_still_exposes_all_original_names():
    """(c) ``import zombi2`` still exposes every name in its ``__all__``."""
    assert len(z.__all__) == 129   # + ConversionModel (gene conversion promoted to core)
    missing = [n for n in z.__all__ if not hasattr(z, n)]
    assert missing == [], f"top-level zombi2 lost names: {missing}"


def test_namespaces_partition_public_api():
    """The namespaces cover exactly the top-level public API (minus primitives).

    Every top-level public name (excluding the deliberately un-namespaced
    low-level primitives) appears in exactly one namespace, and no namespace
    advertises a name that is absent from the top level.
    """
    top_public = set(z.__all__) - LOW_LEVEL

    seen: dict[str, str] = {}
    for ns, names in NAMESPACES.items():
        for n in names:
            assert n in set(z.__all__), f"zombi2.{ns} advertises {n}, absent up top"
            assert n not in seen, f"{n} is in both zombi2.{seen[n]} and zombi2.{ns}"
            seen[n] = ns

    uncovered = top_public - set(seen)
    assert uncovered == set(), f"public names in no namespace: {sorted(uncovered)}"


def test_traits_impl_moved_but_backward_compatible():
    """The trait implementation moved to ``_traits_impl`` without breaking imports.

    ``zombi2/traits.py`` is now a re-export namespace, but the historical
    ``from zombi2.traits import ...`` (including the private ``_expm`` helper)
    must keep working, and the implementation must be importable from its new
    home.
    """
    from zombi2.traits import OrnsteinUhlenbeck as ns_ou, _expm  # noqa: F401
    from zombi2.traits.models import OrnsteinUhlenbeck as impl_ou

    assert ns_ou is impl_ou is z.OrnsteinUhlenbeck
