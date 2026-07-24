# Analyses

Self-contained studies built on ZOMBI2. Each subfolder is independent — its own scripts, data,
figures, and write-up — and regenerates deterministically from fixed seeds. Run the scripts from the
subfolder (they import the installed `zombi2` and read/write paths relative to themselves).

| Study | Question | What it recovers | Regenerate |
|-------|----------|------------------|------------|
| [`red/`](red/) | Does **RED** (the GTDB tree-rescaling measure) recover relative node ages once uneven molecular rates distort branch lengths? | RED holds at real archaeal raggedness (CV = 0.23): Pearson r ≈ 0.94–0.95, nRMSE ≈ 6% of tree depth. | `python red/observable.py && python red/experiment.py && python red/figures.py` |
| [`synteny_inversions/`](synteny_inversions/) | What is the genome **inversion rate** in yeast, inferred by matching gene-order conservation between real genomes and a nucleotide-model simulation down a dated tree (ABC)? | ≈ 3–5×10⁻⁴ inversions per gene·Myr (*Lachancea* 2.7×10⁻⁴, *Kluyveromyces* 4.6×10⁻⁴); the rate is identifiable, the event size is not. | `python synteny_inversions/fit.py && python synteny_inversions/figures.py` |

Each study keeps a `REPORT.md` write-up beside its code. Where a study needs a capability the clean
core does not ship as a public API (RED's estimator, for instance), it carries a **local, faithful
port** in its own folder rather than un-quarantining the package `tools/` — the core stays lean; the
analysis stays reproducible.
