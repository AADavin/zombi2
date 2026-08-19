# The manual

Eight chapters and four appendices, one file each:

| file | chapter |
|---|---|
| `ch1.md` | Introduction (with the tour of the tool) |
| `ch2.md` | Species trees |
| `ch3.md` | Genomes I: gene families |
| `ch4.md` | Genomes II: ordered |
| `ch5.md` | Genomes III: nucleotide |
| `ch6.md` | Sequences |
| `ch7.md` | Traits |
| `ch8.md` | Dependent runs |
| `appendix-a.md` | Rates in detail, and the Gillespie algorithm |
| `appendix-b.md` | Output files |
| `appendix-c.md` | The connection reference |
| `appendix-d.md` | Tools |

`figures/` holds the chapter figures.

[`docs/design/SPEC.md`](../../docs/design/SPEC.md) is the constitution — the model, the words, the
chapter list. This directory is its exposition. When a chapter and SPEC disagree, SPEC wins.

Ch1–Ch8 are published to the docs site under `docs/guide/`, Appendix A as `docs/rates.md`,
Appendix B as `docs/output-files.md`, Appendix C as `docs/guide/connection-reference.md`, and
Appendix D as `docs/tools.md` plus one page per tool, all by snippet include.
Renaming a file here breaks those includes; CI's `mkdocs --strict` will fail if you do.

Run every example before trusting a chapter. Chapters drift behind the code in both directions.

*(The round-1 review apparatus — `dashboard.yml`, `render.py`, `proposed-index.md`, `coverage-audit.md`
— was retired on 2026-07-21, once all nine chapters had landed. It is in the git history.)*
