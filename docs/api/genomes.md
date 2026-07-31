# zombi2.genomes

Level 2: genomes evolving along the species tree. The level has three **resolutions** — family ⊂
ordered ⊂ nucleotide — and one entry point each; more detail costs more compute, so the resolution
is a dial.

::: zombi2.genomes.simulate_genomes_family

::: zombi2.genomes.simulate_genomes_ordered

::: zombi2.genomes.simulate_genomes_nucleotide

## Results

One result type per resolution, each carrying the true history behind the dataset: the gene trees,
the event log, and the genomes themselves.

::: zombi2.genomes.FamilyGenomesResult

::: zombi2.genomes.OrderedGenomesResult

::: zombi2.genomes.NucleotideGenomesResult

::: zombi2.genomes.StreamedRun

## Gene trees

::: zombi2.genomes.GeneTree

::: zombi2.genomes.GeneCopy

## Who receives a transfer

`transfer_to` is a **choice slot**: it redistributes who receives a horizontal transfer without
changing how many transfers happen. A weight of 0 means "cannot receive", and when every candidate
weighs 0 the transfer does not fire.

::: zombi2.genomes.Distance

::: zombi2.genomes.Clades

::: zombi2.genomes.Between
