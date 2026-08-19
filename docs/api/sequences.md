# zombi2.sequences

Level 3: a sequence evolving inside a gene, along its gene tree. The level takes a whole **genome
run**, not bare gene trees — a sequence sees the species tree only through its gene tree, but its
rate is set by the species branch the gene sits on.

::: zombi2.sequences.simulate_sequences

::: zombi2.sequences.SequencesResult

::: zombi2.sequences.StreamedSequences

::: zombi2.sequences.mean_pairwise_identity

## The substitution-model menu

A substitution model is a `K×K` rate matrix `Q`, normalised to one expected substitution per site
per unit branch length, and its stationary frequencies. Different models are genuinely different
matrices, so they stay a **menu** of constructors rather than one grammar.

::: zombi2.sequences.substitution_models
    options:
      members:
        - SubstitutionModel
        - jc69
        - k80
        - hky85
        - gtr
        - reversible
        - poisson
        - jtt
        - dayhoff
        - wag
        - lg
        - decode
        - encode
