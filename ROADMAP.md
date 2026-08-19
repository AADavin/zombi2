# Roadmap

What ZOMBI2 does not do yet, and where each gap stands. The manual documents only what ships,
so this page holds the rest: models we intend to build, models we may build, and models we have
decided against. It is a statement of intent, not a schedule — nothing here has a date.

Each entry carries one of three statuses:

- **planned** — we intend to build it.
- **considering** — a plausible extension with no commitment attached. An issue that argues for
  one of these, with the study it would serve, is how it moves up.
- **not planned** — a deliberate no, with the reason, listed so it is not asked again.

Many of these gaps are also named inside the code: asking for one at the entry point raises an
error saying the model is not built. Those errors are deliberate — each is a statement about the
code, not about the model.

## Species trees

What ships is in [Species trees](https://aadavin.github.io/zombi2/docs/guide/species-trees/) —
including mass extinctions, incomplete sampling, fossil recovery, piecewise time-varying rates,
and diversity-dependent rates.

| Feature | What it adds | Status |
| --- | --- | --- |
| Sampled trees for fossils | Today `fossils=` reports (lineage, time) pairs beside the tree. This would put fossils *on* the tree: a removal probability, sampled-ancestor nodes, and a pruned sampled tree — the fossilized birth–death process in full. | considering |
| Time-varying fossil recovery and sampling | The fossil recovery rate and the sampling fraction are plain numbers; the piecewise (skyline) form both have in the literature is not written yet. | considering |
| Exact conditioning on tip count | With extinction, `n_extant` stops the run the first time the count is reached, a biased sample of the conditioned distribution (the manual quantifies the bias). An exact sampler is a known, harder build. | considering |
| Age-dependent rates (non-exponential waiting times, e.g. Weibull) | A speciation or extinction rate that depends on the lineage's age. Needs the engine to integrate a rate between events rather than read it at a point. | considering |
| Clade rate shifts | A rate shift inherited by one clade. On a growing tree there is no clade to name yet, so this needs its own written form — a shift at a time, in a lineage drawn then. | considering |
| Backward simulation | Growing the reconstructed tree directly, with ghost lineages available on request. | considering |
| Protracted speciation | Speciation as a window rather than an instant: a split opens an incipient lineage, which becomes a species at a completion event — or dies unrecorded. The incipient/good state along a branch is a driver the connection grammar can already read, so gene flow, clock changes or trait jumps tied to the window come with it. | considering |
| Hybridization and species networks | A node with two parents: hybrid speciation, and allopolyploidy as whole-genome duplication by merger. A different object from a tree. Distance-biased transfer already approximates soft introgression. | considering |
| Stochastic rate shifts | Rate shifts arriving as their own process on branches and inherited by the descendants — no clade has to be named in advance. | considering |
| Diversified taxon sampling | Keeping one tip per clade rather than a uniform fraction — the sampling scheme dating studies condition on. | considering |

## Genomes

What ships is in the three genome chapters — the family, ordered and nucleotide resolutions;
duplication, transfer, loss and origination; inversions, translocations, fissions and fusions;
the chromosome network; and founding a nucleotide run from a real annotation (`gff=` + `fasta=`).

| Feature | What it adds | Status |
| --- | --- | --- |
| Real-genome founding at every resolution | `gff=`/`fasta=` found a nucleotide run today. The family and ordered resolutions start only from synthetic genomes; founding them from an annotation is the most requested item from ZOMBI v1 users. | considering |
| Homology at founding | Declaring two annotated genes copies of one family. Today every founding gene starts its own family. | considering |
| Whole-genome duplication | One event that duplicates every family — and, at the structured resolutions, every chromosome. | considering |
| Gene conversion | Ectopic conversion between copies of a family, with a bias toward the older copy. | considering |
| GC-biased gene conversion | Conversion resolved toward G and C, moving base composition — the model that connects gene conversion here with composition change at the sequence level. | considering |
| Pseudogenization | A loss at the nucleotide resolution that demotes the gene to intergene and keeps its sequence. Today a loss deletes. | considering |
| Replacing transfer at the nucleotide resolution | A nucleotide transfer is always additive today. | considering |
| Gene-order files at the ordered resolution | The nucleotide resolution writes GFF and BED; the ordered resolution writes `gene_order.tsv` only. | considering |
| A graph file for the chromosome network | `chromosome_events.tsv` holds the network as an edge list; a standard graph format (GraphML, DOT) would open it to graph tools. | considering |
| Random rate variation on genome rates | A drawn per-lineage or per-chromosome factor, and inherited clade drift. Today a genome rate varies by time, by driver, or by family. | considering |
| Per-family rates and self-reading runs beyond the family resolution | Named families with their own rates, and `joint=True`, are family-resolution only today. | considering |

## Sequences

What ships is in the Sequences chapter — the reversible menu (JC69 through GTR, the empirical
protein matrices), `reversible()` over any alphabet, +I/+G site rates, per-site frequency
profiles, per-clade models, partitions, indels, and the drawn, inherited and driven clocks.

| Feature | What it adds | Status |
| --- | --- | --- |
| Non-reversible substitution models | A model built from a raw rate matrix, its stationary frequencies derived from it rather than supplied — UNREST and strand-asymmetric mutation. | planned |
| Composition change along the tree | Stationary frequencies that move over time or differ between clades while exchangeabilities stay put — GC content drifting toward a new equilibrium. | planned |
| Epoch clock | The substitution rate as a piecewise function of time (`changing_at` at this level). | planned |
| White-noise clock | An uncorrelated clock whose branch factor has variance shrinking with branch duration — the one i.i.d. clock that stays the same model when a branch is split. | planned |
| Free site-rate categories | User-chosen rate categories and weights beside `gamma_shape`. The machinery exists; the written form does not. | considering |
| Codon models | GY94/MG94 with ω, the site-model series on top, and mutation–selection models (site-wise fitnesses) as the natural extension. | considering |
| Site-profile mixtures | `profiles=` takes per-site frequency rows today; a built-in mixture to draw them from is missing, and profiles do not run on the parallel engine yet. | considering |
| Covarion | Hidden rate-switching states over a doubled alphabet. | considering |
| Heterotachy | A site's rate class re-drawn within clades. Refused by name today: it breaks the once-per-family class draw the engine's speed rests on. | considering |
| Autocorrelated rates across sites | Neighbouring sites sharing rate classes instead of independent draws. | considering |
| Indels in the rate grammar | Indel rates and extents are plain numbers today — no modifiers, no CLI flags. | considering |
| Paired-site models | RNA stems as 16-state doublets, the two columns evolving together. Needs an engine where a site reads other sites, which does not exist yet. | considering |

## Traits

What ships is in the Traits chapter — Brownian motion and OU with the full rate grammar,
multi-optimum regimes, correlated traits, punctuated change at speciation, Mk with arbitrary
matrices, the threshold model, and rates driven by other levels in both directions.

| Feature | What it adds | Status |
| --- | --- | --- |
| A paired draw at speciation | Every at-speciation change today draws each daughter independently. Drawing the pair jointly is the one missing piece shared by range division (DEC), cladogenetic transition matrices, and one-daughter-changes models. | considering |
| DEC biogeography | Ranges as sets of areas, with dispersal, extirpation and cladogenetic division. Needs the paired draw above. | considering |
| Range-dependent diversification | GeoSSE: speciation and extinction read the range. Falls out of DEC and the paired draw, joined to the species tree. | considering |
| Hidden states as a first-class model | Compound state labels can spell a hidden-state model today, but the observed/hidden split, its projection in the result, and the character-independent null are all hand work. | considering |
| Along-branch jumps | Pulsed evolution: jumps at their own rate along a branch, not only at splits. | considering |
| A trend on Brownian motion | A directional drift term on the value. | considering |
| A driven OU optimum | Only rates can be driven today; the optimum is a plain argument. An optimum moved by time, an environment, or another trait is a new cell in the grammar — driving a parameter that is not a rate. | considering |
| Full multivariate OU | One trait's deviation pulling another — an off-diagonal pull matrix. Refused by name today. | considering |
| Bounded traits | A reflecting or absorbing bound on a continuous trait. | considering |
| Mk root options | The root state drawn from the stationary distribution or a supplied vector; today it is uniform or fixed. | considering |
| More modifiers on the switch rate | Piecewise time variation, a drawn per-lineage factor, diversity dependence — the discrete switch rate takes only a driver today. | considering |
| Tree transforms | Pagel's λ, δ and κ as tree-to-tree functions. | considering |

## Dependent runs

What ships is in the Dependent runs chapter and the
[connection reference](https://aadavin.github.io/zombi2/docs/guide/connection-reference/) —
all six level pairs, an engine for every joinable pair, and a gallery example for every
connection.

| Feature | What it adds | Status |
| --- | --- | --- |
| Structured genomes joined to the species tree | Gene content driving speciation is family-resolution only today. | considering |
| Genomes and sequences joined at the nucleotide resolution | Refused today because an event can fall inside a gene, so a copy's sequence stops being one string. | considering |
| Continuous traits in more joint pairs | A genome and a continuous trait, or a continuous trait and a sequence, driving each other. Both need the diffusion sliced against the other level's events, as the quantitative-trait speciation model already is. | considering |
| Transfer on a growing tree | A joint run refuses transfer because the set of live lineages is still forming. | considering |
| Drawn variation beside a driven rate in joint runs | Clade drift with a driver, or a relaxed clock with a driver, are refused together today. | considering |
| Three participants | A joint run holds exactly two parts, and a trait cycle holds exactly two traits. | considering |
| An environmental curve as a driver | An external time series — temperature, sea level — read by any level's rate: environment-dependent diversification, a climate-tracking optimum. The grammar reserves the entry point (a measured driver); nothing supplies one beyond time itself. | considering |
| New drivers | `Distance` split into a driver with a mapping of its own (designed, not built). | considering |

Not planned:

- **Species ↔ sequences.** Too far apart to connect: a sequence evolves along a gene tree,
  which evolves along the species tree, so connecting the two would require simulating the
  genome as well.
- **A genome driven by a second genome of the same lineage.** One genome per lineage is the
  model.

## Tools

What ships is in the Tools appendix — `format` (homology, markers, recphylo), `tree`, and
`treedist`.

| Feature | What it adds | Status |
| --- | --- | --- |
| Incomplete lineage sorting | A multispecies coalescent run inside the gene trees a genome run produces — discordance from ILS on top of duplication, transfer and loss. | considering |
| Reconciliation scoring | The other half of the benchmarking loop the recphylo answer keys open: a reconciliation likelihood, an accuracy scorer against the true events, and readers for ALE/AleRax output. Prototyped, not shipped. | considering |
| Quartet distance | A fourth `treedist` metric. | considering |
| Pruning foreign trees | `tree --prune` needs fates, which only a ZOMBI2 run has; a way to hand fates in from elsewhere. | considering |
| Benchmark realism | Missing data, fragmentary sequences and contamination applied to a finished run — the corruption a method meets in practice, added after the truth is written. | considering |
