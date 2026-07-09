# Parsing reconciliation output (reconparser)

**reconparser** is the *interop* complement to [ALElite](reconciliation-likelihood.md). Where
ALElite **computes** the reconciliation likelihood of a ZOMBI2 gene family, reconparser **reads**
what the established reconciliation programs already wrote — so you can pull a real reconciliation
into Python and line it up against a ZOMBI2 simulation of the same system.

It parses two formats:

| Software | Versions | What it reads |
| --- | --- | --- |
| **ALE** (`ALEml`, `ALEml_undated`) | v0.4, v1.0 | `.ucons_tree` (consensus gene tree), `.uTs` (transfers), `.uml_rec` (species tree, reconciled gene trees, ML D/T/L rates, log-likelihood, summary + per-branch event tables) |
| **AleRax** | v1.2+ | a whole run directory — species trees, per-gene model parameters, per-family likelihoods, global and per-family transfers, per-species event counts, origination probabilities, coverage |

Trees come back as native [`zombi2.tree.Tree`](../reference/api.md) objects — parsed with
ZOMBI2's own `read_newick`, the same tree type the simulator produces — and tables as
[`pandas`](https://pandas.pydata.org/) DataFrames. That means a parsed reconciliation drops
straight into the rest of ZOMBI2 (the [reconciliation scorer](reconciliation-likelihood.md), the
tree utilities) with no foreign tree type to convert. The reconciliation annotations ALE and
AleRax bake into node names — `.T@donor->recipient` for transfers, `.D@…` for duplications,
consensus support values — survive verbatim as each node's `name`. AleRax parsing is **lazy** —
nothing is read until you ask for it, which matters for runs with thousands of families.

## Optional dependency

reconparser needs only `pandas` (for the table outputs), which is **not** in the base install;
trees use ZOMBI2's built-in parser, so there is no `ete3` dependency. Like the `selection` extra,
it ships behind an opt-in extra:

```bash
pip install 'zombi2[reconparser]'
```

Nothing is re-exported into the top-level `zombi2` namespace, so plain `import zombi2` never pulls
in `pandas` — you only pay for it when you import the subpackage.

## From Python

```python
from zombi2.tools.reconparser import ALEParser, AleRaxRun

# --- classic ALE ---
ale = ALEParser("results.ale")            # base path, or any of its .ucons_tree/.uTs/.uml_rec files
print(ale.get_ml_rates())                 # {'duplications': .., 'transfers': .., 'losses': ..}
print(ale.get_log_likelihood())
transfers = ale.get_transfers()           # DataFrame: from, to, freq
branch_stats = ale.get_branch_statistics()

# --- AleRax run directory ---
run = AleRaxRun("alerax_output/")
print(run.get_run_info())                 # version, num_families, num_species, model, ...
print(run.get_total_log_likelihood())
fam = run.get_family("K00192")            # lazy per-family view
fam.get_event_counts(sample=0)            # {'S': .., 'D': .., 'T': .., 'L': .., ...}
```

`ALEParser`, `AleRaxRun` and `AleRaxFamily` expose the full run/family surface (consensus and
sampled gene trees, model parameters, per-species events, coverage, origins, …); each getter is
cached after its first call. See the class docstrings for the complete method list.

## From the command line

`zombi2 tools parse` reads a run and prints a summary — the ML DTL rates, the log-likelihood, and
the top transfers. The tool is auto-detected from the path (a directory is an AleRax run,
anything else is classic ALE); override with `--tool`.

```bash
# summarize a classic ALE result (base path, without the .uml_rec extension)
zombi2 tools parse results.ale

# summarize an AleRax run directory, show the 20 strongest transfers,
# and also write the transfer / per-family tables as TSV into out/
zombi2 tools parse alerax_output/ --top 20 -o out/
```

```text
ALE reconciliation: results.ale
  files present: consensus_tree, transfers, reconciliation
  log-likelihood: -42.500000
  ML rates:  D=0.1  T=0.05  L=0.15
  total events:  D=3  T=2  L=5  S=4
  transfers: 2 edge(s) (top 10 by frequency)
     A -> B   0.500
     B -> C   0.300
```

With `-o DIR` the command also writes the transfer table (and, for ALE, the per-branch table; for
AleRax, the per-family likelihoods) as tab-separated files you can load elsewhere.

## Scope

reconparser is a **reader**, not an inference engine: it turns files into trees and tables and
stops there. To *score* a ZOMBI2 gene family under the ALE models, use
[ALElite](reconciliation-likelihood.md); to *run* the reconciliation in the first place, use ALE
or AleRax. reconparser is the bridge between them and ZOMBI2.

!!! note "Provenance"
    reconparser is vendored from the standalone
    [`reconparser`](https://github.com/AADavin/reconparser) library (MIT, same author) so it ships
    with ZOMBI2 rather than requiring a separate install. The parser modules are kept close to
    upstream to make future syncs easy.
