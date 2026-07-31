# The manual

Nine chapters and three appendices, one file each:

| file | chapter |
|---|---|
| `ch1.md` | Introduction |
| `ch2.md` | A tour of ZOMBI2 |
| `ch3.md` | Species trees |
| `ch4.md` | Genomes I: gene families |
| `ch5.md` | Genomes II: ordered |
| `ch6.md` | Genomes III: nucleotide |
| `ch7.md` | Sequence evolution |
| `ch8.md` | Trait evolution |
| `ch9.md` | Conditioning |
| `ch10.md` | Joining |
| `appendix-a.md` | Rates in detail, and the Gillespie algorithm |
| `appendix-b.md` | Output files |
| `appendix-c.md` | Tools |

`figures/` holds the chapter figures.

[`docs/design/SPEC.md`](../../docs/design/SPEC.md) is the constitution — the model, the words, the
chapter list. This directory is its exposition. When a chapter and SPEC disagree, SPEC wins.

Ch3 and Ch4 are also published to the docs site, and Appendix B is single-sourced into
`docs/reference/output-files.md`, all by snippet include. Renaming a file here breaks those includes;
CI's `mkdocs --strict` will fail if you do.

Run every example before trusting a chapter. Chapters drift behind the code in both directions.

*(The round-1 review apparatus — `dashboard.yml`, `render.py`, `proposed-index.md`, `coverage-audit.md`
— was retired on 2026-07-21, once all nine chapters had landed. It is in the git history.)*
