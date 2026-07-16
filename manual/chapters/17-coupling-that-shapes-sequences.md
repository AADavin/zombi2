# Coupling that shapes sequences

The last tier of the diamond is **sequences** (Σ). Sequences are a *target only*: a trait or a gene
event can bend how a sequence evolves, but a sequence drives nothing back. And because a sequence
rides its *gene* tree — downstream of the genome layer, not the species tree — these couplings live on
the `zombi2 sequences` command rather than `coevolve`, replaying a genome run's event trace along its
gene trees. There is deliberately no species–sequence edge (that diagonal skips a tier).

Two things about a sequence can be driven, so three edges point into Σ:

- its **selection** — the strength of purifying/positive selection, read out as `dN/dS` ($\omega$) — driven
  by a trait (`traits:selection`) or by a gene event such as duplication (`genomes:selection`);
- its **substitution speed** — the molecular clock — driven by a trait (`traits:speed`).

Each is one `sequences --couple driver:target-variable` invocation. They all run on a prior genome
simulation whose event trace was written:

```bash
# a species tree, then gene families with the event trace the sequence sim replays
zombi2 species --tips 30 --age 5 --seed 1 -o run/
zombi2 genomes -t run/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace -o run/
```

Three flags shape any of these couplings: `--couple-strength` is the exp-link coefficient (how hard
the driver pushes), `--couple-base-omega` the baseline `dN/dS` a selection edge modulates, and
`--couple-trait-sigma` the Brownian variance of the driving trait for the `traits:*` edges.

## `traits:selection` — trait-driven dN/dS

A continuous trait diffuses along every lineage, and its value sets that lineage's selection strength:
`dN/dS` rides an exponential link $\omega = \omega_0\,\exp(\beta s)$ of the trait value $s$ — with
$\omega_0$ set by `--couple-base-omega` and $\beta$ by `--couple-strength` — so lineages that drift to
a high trait value fall under relaxed (or positive) selection and those that drift low are held under
tighter purifying selection. The edge builds its own `GY94` codon model, so you do **not** pass a
`--subst-model`:

```bash
zombi2 sequences --genomes run/ --couple traits:selection \
    --couple-strength 1.5 --couple-base-omega 0.2 --couple-trait-sigma 0.5 \
    --seq-length 300 --seed 2 -o run/
```

**What it recovers:** `dN/dS` varies from branch to branch in step with a latent phenotype — the
lineage-heterogeneous selection that codon models try to detect, here with the trait that caused it on
record. It writes the codon alignments, the substitution-unit gene trees, and `branch_rates.tsv`.

## `genomes:selection` — post-duplication relaxed selection

Here a **gene event** relaxes selection: after a duplication, the redundant copy is freed from
purifying selection, so its `dN/dS` rises for a while — the classic signature of neofunctionalisation.
The coupling reads the genome event trace and lifts $\omega$ on the branches following a duplication,
again on a `GY94` model built for you:

```bash
zombi2 sequences --genomes run/ --couple genomes:selection \
    --couple-strength 2 --couple-base-omega 0.15 \
    --seq-length 300 --seed 2 -o run/
```

**What it recovers:** duplicated copies show elevated `dN/dS` relative to the single-copy background —
the relaxed-selection burst that follows gene duplication, written into the codon alignments as a
known ground truth. `--couple-strength` sets how much selection relaxes; with it at `0` the edge is cut
and every branch keeps the baseline $\omega$.

## `traits:speed` — a trait-driven clock

The third edge drives the substitution **rate** rather than selection: a trait scales the molecular
clock, so lineages with a high trait value accumulate substitutions faster. Unlike the two selection
edges, this evolves ordinary sequences under a substitution model you choose, and the coupling replaces
the usual lineage clock with a trait-driven one:

```bash
zombi2 sequences --genomes run/ --couple traits:speed \
    --couple-strength 1.0 --couple-trait-sigma 0.6 --subst-model gtr \
    --seq-length 300 --seed 2 -o run/
```

**What it recovers:** among-lineage rate variation that *tracks a phenotype* — a relaxed clock with a
known cause, rather than the phenotype-blind random clocks of [Relaxed molecular
clocks](#molecular-clocks). The drawn per-branch rates are written to `branch_rates.tsv`, so the
trait–rate link that generated the phylograms is on record.

::: note
The sequence couplings are the newest tier of the diamond and the youngest part of ZOMBI2. Today they
are driven from the `zombi2 sequences` command as shown; a Python entry point and matched `--null`
runs are on the near-term roadmap, alongside a `genomes:speed` edge (a gene event scaling the clock).
:::
