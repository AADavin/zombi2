# Analyses

Self-contained validation studies built on ZOMBI2. Each subfolder is independent — its own scripts,
data, figures and `REPORT.md` — and regenerates deterministically from fixed seeds. Run the scripts
from the subfolder; they import the installed `zombi2` and read and write paths relative to
themselves.

**[`red/`](red/) — can you trust RED?** Does the measure GTDB uses to normalise taxonomic ranks
across the tree of life still recover relative node ages once molecular rates vary? It is also the
docs site's worked example ([`docs/example-red.md`](../docs/example-red.md)), and its `figures.py`
writes into `docs/assets/red/` as well as its own folder, so the page and the study cannot drift
apart.

```bash
cd analyses/red && python observable.py && python experiment.py && python figures.py
```

Where a study needs a capability the clean core does not ship as a public API (RED's estimator, for
instance), it carries a **local, faithful port** in its own folder rather than un-quarantining the
package `tools/` — the core stays lean and the analysis stays reproducible.
