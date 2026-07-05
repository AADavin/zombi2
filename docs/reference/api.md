# API reference

Auto-generated from the source docstrings.

## Species trees

::: zombi2.species_model.BirthDeath
::: zombi2.species_model.Yule
::: zombi2.species_model.EpisodicBirthDeath
::: zombi2.species_sim.simulate_species_tree
::: zombi2.ghosts.add_ghost_lineages
::: zombi2.tree.prune

## Tree

::: zombi2.tree.Tree
::: zombi2.tree.TreeNode
::: zombi2.tree.read_newick

## Relaxed molecular clocks

::: zombi2.rate_variation.Clock
::: zombi2.rate_variation.StrictClock
::: zombi2.rate_variation.UncorrelatedLogNormalClock
::: zombi2.rate_variation.UncorrelatedGammaClock
::: zombi2.rate_variation.WhiteNoiseClock
::: zombi2.rate_variation.AutocorrelatedLogNormalClock
::: zombi2.rate_variation.CIRClock
::: zombi2.rate_variation.RateVariation
::: zombi2.rate_variation.RateScaledTree
::: zombi2.sequence_evolution.SequenceEvolution
::: zombi2.sequence_evolution.GenePhylograms

## Trait evolution

::: zombi2.traits.simulate_traits
::: zombi2.traits.TraitResult
::: zombi2.traits.BrownianMotion
::: zombi2.traits.OrnsteinUhlenbeck
::: zombi2.traits.EarlyBurst
::: zombi2.traits.MultivariateBrownian
::: zombi2.traits.MultivariateOU
::: zombi2.traits.MultiOptimumOU
::: zombi2.traits.Mk
::: zombi2.traits.CorrelatedBinary
::: zombi2.traits.HiddenStateMk
::: zombi2.traits.ThresholdModel
::: zombi2.traits.pagel_lambda
::: zombi2.traits.pagel_delta
::: zombi2.traits.pagel_kappa

## Historical biogeography

::: zombi2.biogeography.simulate_biogeography
::: zombi2.biogeography.DEC

## Gene-family coupling (Potts)

::: zombi2.coupling.simulate_coupled
::: zombi2.coupling.CouplingSpec
::: zombi2.coupling.pathway_blocks
::: zombi2.coupling.PottsRates
::: zombi2.coupling.CoupledResult

## Trait-linked gene families

::: zombi2.trait_coupling.simulate_trait_linked_genomes
::: zombi2.trait_coupling.TraitGeneCoupling
::: zombi2.trait_coupling.TraitLinkedResult
::: zombi2.trait_coupling.TraitLinkedRates
::: zombi2.trait_coupling.TraitTrajectory

## Rate models

::: zombi2.rates.RateModel
::: zombi2.rates.SharedRates
::: zombi2.rates.PerGenomeRates
::: zombi2.rates.FamilySampledRates
::: zombi2.rates.BranchRates
::: zombi2.rates.EventWeight

## Distributions

::: zombi2.distributions.Distribution
::: zombi2.distributions.Fixed
::: zombi2.distributions.Exponential
::: zombi2.distributions.Gamma
::: zombi2.distributions.LogNormal
::: zombi2.distributions.Uniform
::: zombi2.distributions.as_distribution

## Transfers

::: zombi2.transfers.TransferModel

## Genomes

::: zombi2.genome.Gene
::: zombi2.genome.Genome
::: zombi2.genome.UnorderedGenome
::: zombi2.genome.OrderedGene
::: zombi2.genome.OrderedGenome

## Simulation driver

::: zombi2.simulation.simulate_genomes
::: zombi2.simulation.Genomes
::: zombi2.genome_sim.GenomeSimulator

## Profiles

::: zombi2.profiles.ProfileMatrix

<!-- ABC profile-matching inference (zombi2.matching) is experimental and withheld from the v1
     public API; its reference entries return when the module is stabilised. -->


## Parallel replicates

::: zombi2.parallel.run_replicates

## Rust engine

The built-in model runs on the Rust engine automatically via `simulate_genomes` (see the
[Rust engine](../guide/rust-engine.md) guide); `rust_available()` reports whether the compiled
extension is installed.

::: zombi2._rust.available

## Gene-tree reconstruction

::: zombi2.reconciliation.build_gene_trees

## Events

::: zombi2.events.EventType
::: zombi2.events.EventRecord
::: zombi2.events.GeneOp
::: zombi2.events.Selection
