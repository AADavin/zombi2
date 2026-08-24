# The connection reference
The tables of Chapter 8, in one place: every connection, every driver, every target, and the verbs and mappings of the link.

## Every connection

![What can connect to what. Rows drive, columns are driven, and the numbers are the rows of the table below. A shaded cell cannot be connected: three would need two genome runs for the same lineage, and the rest pair levels that cannot depend on each other. A cell with the arrow with two heads runs joint only, because neither level can be simulated first. A boxed cell is one part of a level driving another part of the same level. A starred cell can also be written on the command line.](figures/conditioning_map_print.png){width=95%}

| # | Driver | Target | What it says | Conditioned | Joint |
|---|---|---|---|---|---|
| **1** | a trait | a gene family | habitat sets the loss rate; all four rates take a driver | [Co1–Co8](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:genome_reduction--><!--gallery:genome_expansion--><!--gallery:hgt_uptake--><!--gallery:continuous_conditioning--><!--gallery:curve_saturating--><!--gallery:curve_optimum--><!--gallery:set_by_habitat--><!--gallery:scalar_response--> | [Jo5](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:cave_genomes--> |
| **2** | a trait | an ordered or nucleotide genome | eleven rates at the ordered resolution, thirteen at the nucleotide one, and the extents besides | [Co9](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_inversions--> | — |
| **3** | a trait | a sequence | habitat sets the substitution rate | [Co10](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_substitution--> | [Jo2](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:trait_and_sequence--> |
| **4** | a trait | a trait | one character sets another's `rate` or `switch` | [Co11–Co12](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:driven--><!--gallery:trait_drives_trait--> | [Jo4](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:trait_loop--> |
| **5** | a gene family | a gene family | a mobile element makes transfer likelier for the rest of the genome | [Co13](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:mobile_element--> | [Jo6](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:mobile_element_joint--> |
| **6** | a gene family | a sequence | lose the repair gene and evolve faster | [Co14](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:repair_gene--> | — |
| **7** | a gene family | a trait | carry the toxin family and turn pathogenic | [Co15](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gene_drives_trait--> | [Jo5](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:cave_genomes--> |
| **8** | an ordered or nucleotide genome | a sequence | as **6**, with coordinates in the genome run | [Co16](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:operon_substitution--> | — |
| **9** | an ordered or nucleotide genome | a trait | as **7**, with coordinates in the genome run | [Co17–Co18](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:module_drives_motility--><!--gallery:operon_trait--> | — |
| **10** | a sequence | a sequence | one gene's composition indexes something about the lineage, and that sets another gene's rate | [Co19–Co20](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gc_drives_sequence--><!--gallery:named_family_drives_sequence--> | [Jo3](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:sequence_loop--> |
| **11** | a sequence | a trait | GC content sets how fast a trait changes | [Co21](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gc_drives_trait--> | [Jo2](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:trait_and_sequence--> |
| **12** | a sequence | a genome | composition sets the loss rate of the genome that carries the gene; it can never be simulated first | — | [Jo1](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:genome_and_sequence--> |
| **13** | a trait | the species tree | speciation and extinction follow the state; the tree comes out of the run | — | [Jo8–Jo12](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:bisse--><!--gallery:state_extinction--><!--gallery:musse--><!--gallery:classe--><!--gallery:quasse--> |
| **14** | gene content | the species tree | carrying a key family sets the speciation rate | — | [Jo7](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:key_innovation--> |

## What can be a driver

In a conditioned run, a driver is a finished result:

| Driver | What it offers | Gallery |
|---|---|---|
| a discrete trait | one of its states | [Co1](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:genome_reduction--> |
| a continuous trait | a number, taken every `step` of time | [Co4](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:continuous_conditioning--> |
| a gene family | `present` or `absent` | [Co14–Co15](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:repair_gene--><!--gallery:gene_drives_trait--> |
| a module | a fraction, 0 to 1: how complete a declared group of families is | [Co17](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:module_drives_motility--> |
| a sequence's composition | a number, 0 to 1: the share of its letters in a given set | [Co19](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gc_drives_sequence--> |
| **one family's** composition | the same, on a run restricted to that family, plus an `absent=` | [Co20](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:named_family_drives_sequence--> |

In a joint run, a driver is named, and the name says which level it comes from:

| Driver | What it offers | Mapping |
|---|---|---|
| `"traits:<name>"`, discrete | that trait's current state | a table over the states |
| `"traits:<name>"`, continuous | that trait's current value | a curve, or a `Scalar` |
| `"genomes:<family>"` | whether that family is there | a table over `present` / `absent` |
| `"genomes:count"` | how many genes the lineage has | a curve, or a `Scalar` |
| `"sequences:<name>"` | how much of that gene is a given set of letters | a curve, or a `Scalar` |

## What can be a target

| Target | Kind | Level | Gallery |
|---|---|---|---|
| `duplication`, `transfer`, `loss`, `origination` | how often | genomes, every resolution | [Co1](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:genome_reduction--> |
| `inversion`, `transposition`, `translocation`, `fission`, `fusion`, `chromosome_origination`, `chromosome_loss` | how often | genomes, ordered and nucleotide | [Co9](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_inversions--> |
| `substitution` | how often | sequences | [Co10](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_substitution--> |
| `rate` (continuous), `switch` (discrete) | how often | traits | [Co11–Co12](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:driven--><!--gallery:trait_drives_trait--> |
| every event's extent | how much | genomes, ordered and nucleotide; Python only | [Co9](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_inversions--> |
| `transfer_to` | a choice | genomes, every resolution | [Co3](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:hgt_uptake--> |
| the substitution model | a model | sequences | [Sq2](https://aadavin.github.io/zombi2/gallery.html#sequences)<!--gallery:clade_own_model--> |
| `birth`, `death` | how often | species; joint runs only | [Jo8](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:bisse--> |

## The verbs

| Verb | What the number does | Written on | Gallery |
|---|---|---|---|
| `scaled_by` | multiplies the base in front of it | a rate, an extent | [Co1](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:genome_reduction--> |
| `set_by` | replaces the base, in the rate's own units, so nothing is written in front; on a substitution model it gives the model itself | a rate, a substitution model | [Co7](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:set_by_habitat--> |
| `weighted_by` | weighs the candidates against each other | `transfer_to`, from `Recipients()` | [Co3](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:hgt_uptake--> |

`set_by` is for when the literature states the rate itself (the loss rate is 1.0 in the water) rather than a multiple of a base you had to invent. `weighted_by` needs no base because its weights are normalised over whoever is alive when a transfer fires, so they choose the recipient and never change how many transfers happen.

## The mappings

| The driver gives | Mapping | What it is | Gallery |
|---|---|---|---|
| named states | `Table` | one factor per state; a bare dict is read as one, and a state left out is unchanged | [Co1](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:genome_reduction--> |
| a number | `Curve` | any function of it; a bare function is read as one, and `bound=` puts a ceiling on the factor, which an unbounded function needs | [Co5–Co6](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:curve_saturating--><!--gallery:curve_optimum--> |
| a number | `Scalar` | a strength, giving the factor `exp(strength × value)` | [Co8](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:scalar_response--> |
| a pair of states | `Between` | one weight per ordered (donor, recipient) pair, `default=` for the rest; `transfer_to` only | [Ge3](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_transfer_highway--> |

Whatever the shape, the number that comes out of a mapping has no units and cannot be negative: it multiplies, or weighs, what it lands on. `set_by` is the one exception: its mapping gives the rate itself, in the rate's own units.
