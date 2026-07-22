# Analyses

Self-contained studies built on ZOMBI2. Each subfolder is independent — its own
scripts, figures, and write-up — and regenerates deterministically from fixed seeds.
Run the scripts from the repository root.

| Study | Question | Regenerate |
|-------|----------|------------|
| [`ecoli_nt/`](ecoli_nt/) | Validation of the nucleotide genome model on a real *E. coli* genome (GFF-imported gene layout). | `python analyses/ecoli_nt/run_analysis.py` |

Each study keeps its own `README.md` (where present) and a LaTeX `report/` write-up.
See the individual folders for details and options.
