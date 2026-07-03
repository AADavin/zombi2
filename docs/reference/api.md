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

## Rate variation

::: zombi2.rate_variation.RateVariation
::: zombi2.rate_variation.RateScaledTree

## Rate models

::: zombi2.rates.RateModel
::: zombi2.rates.UniformRates
::: zombi2.rates.GenomeWiseRates
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

## Profile matching (ABC)

::: zombi2.matching.match_profiles
::: zombi2.matching.match_profiles_smc
::: zombi2.matching.ABCFit
::: zombi2.matching.default_summary
::: zombi2.matching.default_gene_tree_summary
::: zombi2.matching.event_count_summary
::: zombi2.matching.frequency_spectrum
::: zombi2.matching.genome_sizes
::: zombi2.matching.copy_number_spectrum
::: zombi2.profiles.ProfileMatrix

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
