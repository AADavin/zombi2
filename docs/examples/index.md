# Examples

ZOMBI2 is a tool to generate datasets in which the user knows everything about them. Let's
see how that is used in practice. Each example asks a question about a real method and
answers it with a dataset where the truth is known. Every file behind each example, the
scripts, the results and the exact reproduction recipe, is in the repository under
[`analyses/`](https://github.com/AADavin/zombi2/tree/main/analyses).

**[Can RED be trusted?](red.md)** GTDB uses Relative Evolutionary Divergence to align
taxonomic ranks across the tree of life; the measure assumes branch length tracks time. We
measure how rate-variable real archaea are, simulate trees with the same variability, and
test RED on trees whose node ages are known.

**[Can BiSSE find the gene that drives speciation?](bisse.md)** A gene family raises the
speciation rate of the lineages that carry it, in a joint run where the genome shapes its
own tree. The dataset includes a matched null and an undriven family, so the standard test
can be scored on known truth: how often it fires when nothing drives, and how often it
finds the gene that does.

**[Which rearrangement parameters can be recovered from gene order?](rearrangements.md)**
Genomes are simulated at known inversion and translocation rates and the parameters
inferred back. The rates are recovered, including the mix of the two event types. The
size of the events is not, and fixing it at a wrong value biases the rates that are
recoverable.

**[Can Pagel's test detect a feedback?](pagel.md)** A habitat drives gene loss across
the genome, and the absence of one gene family drives the switch rate into the
parasitic habitat, a feedback closed in one joint run. Pagel's test detects the
feedback and the switch-rate connection in about nine replicates of ten, largely misses the
loss connection, and rejects at the nominal rate with no connections and on a control
family carried on the same trees.

**[Who trades genes with whom?](transfers.md)** Gene transfers are drawn preferentially
between lineages sharing a habitat, and the habitats of ancestral lineages are read
back from the transfer network that ALE infers. A vote over the inferred transfers
classifies all ancestral branches at 92% accuracy, equal to the ceiling a perfect
reconciliation could reach; the cost of extinction, of dead donors and of the
reconciliation itself is measured at every step.
