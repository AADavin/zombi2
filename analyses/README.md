# Analyses

Self-contained studies built on ZOMBI2. Each subfolder is independent — its own scripts, data,
figures, and write-up — and regenerates deterministically from fixed seeds. Run the scripts from the
subfolder (they import the installed `zombi2` and read/write paths relative to themselves).

| Study | Question | What it recovers | Regenerate |
|-------|----------|------------------|------------|
| [`red/`](red/) | Does **RED** (the GTDB tree-rescaling measure) recover relative node ages once uneven molecular rates distort branch lengths? | RED holds at real archaeal raggedness (CV = 0.23), under every arrangement of rate variation the core can put it under: Pearson r = 0.94–0.95 uncorrelated (nRMSE ≈ 6% of tree depth), 0.99 autocorrelated (2.3%). Published as a worked example on the docs site. | `python red/observable.py && python red/experiment.py && python red/figures.py` |

`red/` is also the docs site's worked example ([`docs/example-red.md`](../docs/example-red.md)); its
`figures.py` writes into `docs/assets/red/` as well as its own folder, so the two cannot drift.

Each study keeps a `REPORT.md` write-up beside its code. Where a study needs a capability the clean
core does not ship as a public API (RED's estimator, for instance), it carries a **local, faithful
port** in its own folder rather than un-quarantining the package `tools/` — the core stays lean; the
analysis stays reproducible.
