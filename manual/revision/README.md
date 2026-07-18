# Manual revision — round 1

Adrián read the manual slowly on 2026-07-17 and produced 134 comments. This directory is the
apparatus for working through them without losing any, which was his explicit worry:

> From previous experience, I see that in these cases of very long feedback, many comments tend to
> be just dropped or ignored.

## The files

| file | what it is |
|---|---|
| `dashboard.yml` | **The source of truth.** One entry per comment, plus the decisions, themes, contradictions, corrections and the wave plan. Hand-edit this. |
| `render.py` | Draws `dashboard.yml` as a single self-contained HTML page. Never edit the HTML. |
| `dashboard.html` | Generated. Open it in a browser, or read it as a published artifact. |

```bash
python manual/revision/render.py     # dashboard.yml -> dashboard.html
```

## The one rule

**Every id in `dashboard.yml` must reach a terminal `outcome`.** That is the whole point: the
register is a contract that nothing gets dropped silently. When a chapter lands, set each of its
comments to `done`, `dropped` or `deferred` — and if `dropped`, say why in the same entry. The diff
on this file is the evidence of what a chapter rewrite actually addressed.

`status` describes the comment as it was received (`actionable`, `needs-decision`, `stale-claim`,
`already-done`, `conflicts`). `outcome` describes what we did about it. They are different fields
and both matter.

## The shape of the work

The 134 comments are not 134 edits: 90 are downstream of 13 decisions (`D1`–`D13`), and 5 of those
are foundational. That is why the plan starts with decision sessions rather than with Chapter 1.
Only 34 comments are both unblocked and actionable today — those are Wave 1, which needs no meeting.

Chapter 2 is the bottleneck. It owns the vocabulary every other chapter quotes, so no chapter
downstream can be rewritten until D3 (the rate vocabulary) is settled.

## Two things to know before deciding anything

1. **Some comments reverse decisions ratified days earlier** (the ch4/ch5 split, the by-target
   coevolution cut, the "diamond", "resolution"). See `contradictions` in the YAML. Reverse them if
   you like, but knowingly — and update `../REVISION_NOTES.md` when you do, because that record
   currently contradicts the book.
2. **Some comments rest on a premise that is out of date** (see `stale_claims`). Several are wrong
   in Adrián's favour: `z.simulate_traits` already exists, traits already run on gene trees, the
   within-branch stochastic map is already exact, `--driver-extinction` is already signed. Those
   close with a few sentences of prose rather than with code.

## Provenance

Built by 25 agents reading every chapter against the code, then one synthesis pass. Every factual
claim in the register carries a `file:line`. Where a claim could not be verified, the entry says so.
The analysis is a snapshot of `main` at 2026-07-17 — re-verify anything load-bearing before acting
on it, especially the code claims.
