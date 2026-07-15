# API reference

Auto-generated from the source docstrings.

## Species trees

::: zombi2.species.model.BirthDeath
::: zombi2.species.model.Yule
::: zombi2.species.model.EpisodicBirthDeath
::: zombi2.species.sim.simulate_species_tree
::: zombi2.species.ghosts.add_ghost_lineages
::: zombi2.tree.prune

## Tree

::: zombi2.tree.Tree
::: zombi2.tree.TreeNode
::: zombi2.tree.read_newick

## Relaxed molecular clocks

::: zombi2.sequences.clocks.Clock
::: zombi2.sequences.clocks.StrictClock
::: zombi2.sequences.clocks.UncorrelatedLogNormalClock
::: zombi2.sequences.clocks.UncorrelatedGammaClock
::: zombi2.sequences.clocks.WhiteNoiseClock
::: zombi2.sequences.clocks.AutocorrelatedLogNormalClock
::: zombi2.sequences.clocks.CIRClock
::: zombi2.sequences.clocks.RateVariation
::: zombi2.sequences.clocks.RateScaledTree
::: zombi2.sequences.evolution.SequenceEvolution
::: zombi2.sequences.evolution.GenePhylograms

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

::: zombi2.traits.biogeography.simulate_biogeography
::: zombi2.traits.biogeography.DEC

## Trait-linked gene families

::: zombi2.coevolve.trait_coupling.simulate_trait_conditioned_genomes
::: zombi2.coevolve.trait_coupling.TraitGeneCoupling
::: zombi2.coevolve.trait_coupling.TraitGeneResult
::: zombi2.coevolve.trait_coupling.TraitGeneRates
::: zombi2.coevolve.trait_coupling.TraitTrajectory

## Rate models

::: zombi2.genomes.rates.RateModel
::: zombi2.genomes.rates.Rates
::: zombi2.genomes.rates.FamilySampledRates
::: zombi2.genomes.rates.LineageRates
::: zombi2.genomes.rates.EventWeight
::: zombi2.genomes.read_rates.read_family_rates
::: zombi2.genomes.read_rates.read_lineage_rates

## Distributions

::: zombi2.distributions.Distribution
::: zombi2.distributions.Fixed
::: zombi2.distributions.Exponential
::: zombi2.distributions.Gamma
::: zombi2.distributions.LogNormal
::: zombi2.distributions.Uniform
::: zombi2.distributions.as_distribution

## Transfers

::: zombi2.genomes.transfers.TransferModel

## Gene conversion

::: zombi2.genomes.conversion.ConversionModel

## Genomes

::: zombi2.genomes.genome.Gene
::: zombi2.genomes.genome.Genome
::: zombi2.genomes.genome.UnorderedGenome
::: zombi2.genomes.genome.OrderedGene
::: zombi2.genomes.genome.OrderedGenome

## Simulation driver

::: zombi2.genomes.simulation.simulate_genomes
::: zombi2.genomes.simulation.Genomes
::: zombi2.genomes.genome_sim.GenomeSimulator

## Profiles

::: zombi2.genomes.profiles.ProfileMatrix

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

::: zombi2.genomes.reconciliation.build_gene_trees

## Events

::: zombi2.genomes.events.EventType
::: zombi2.genomes.events.EventRecord
::: zombi2.genomes.events.GeneOp
::: zombi2.genomes.events.Selection
