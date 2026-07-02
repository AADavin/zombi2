# Non-independence of gene families

The headline scientific goal of ZOMBI2 is to model the **non-independence of gene
families**: functionally coupled families (same pathway or complex) are present or absent
*together* across genomes. Today's rate models treat families independently; making that
correlation real — in a way that phylogenetic-profiling / inverse-Potts methods can recover
— is the next major extension.

A didactic write-up covers the background and the plan:

**:material-file-pdf-box: [non_independence.pdf](non_independence.pdf)** — from phylogenetic
profiles to Potts models.

It explains, intuitively and then formally:

1. **The papers** — Pellegrini et al. (1999), Croce et al. (2019, PhyDCA), and
   Fukunaga & Iwasaki (2022, inverse Potts): why correlated presence/absence predicts
   shared function, and how direct-coupling analysis separates true couplings from
   transitive and phylogenetic artifacts.
2. **The maths** — the Ising and Potts models, the energy `E(σ) = −Σ hᵢσᵢ − Σ Jᵢⱼσᵢσⱼ`
   and its Boltzmann distribution, the inverse problem (inferring `J` from data), and the
   two confounders (shared ancestry, transitivity).
3. **How to simulate it in ZOMBI2** — three designs (pathway-completeness rates; direct
   Potts couplings via a detailed-balance gain/loss process = Glauber dynamics; a hybrid),
   how they map onto a `PottsRates(RateModel)` reading `genome.presence_vector(σ)`, and an
   inject-and-recover validation experiment.

This lands through the existing [rate-model seam](guide/extending.md) — no change to the
simulation engine.
