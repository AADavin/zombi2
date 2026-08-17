# Examples

ZOMBI2 is a tool to generate datasets in which the user knows everything about them. Let's
see how that is used in practice. Each example asks a question about a real method, builds
the dataset that can answer it, and reads the answer off. Every file behind each example,
the scripts, the results and the exact reproduction recipe, is in the repository under
[`analyses/`](https://github.com/AADavin/zombi2/tree/main/analyses).

**[Can RED be trusted?](red.md)** Relative Evolutionary Divergence turns a phylogram into
a relative divergence scale, and GTDB uses it to align taxonomic ranks across the tree of
life. That works only if branch length stands in for time, and real data cannot say
whether it does, because the true node ages are exactly what a real tree withholds. We
measure how rate-variable real archaea are with one model-free number, simulate trees
exactly that variable, and grade RED where the ages are known.

**[Can BiSSE find the gene that drives speciation?](bisse.md)** A gene family multiplies a
lineage's speciation rate while it is present: a genomic key innovation, simulated as one
joint run in which the genome shapes the very tree it evolves along. The dataset carries
its own controls, a matched null and an undriven family in the same genomes, so the
standard test for state-dependent diversification can be scored on data where the truth
is known: how often it fires when nothing drives, and how often it finds the gene that
does.
