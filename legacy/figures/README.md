# Quarantined figures

The figure set of the manual as it stood before the rewrite. Nothing here is referenced by a
chapter, a docs page, the README or SPEC; most of it illustrates models that now live in `legacy/`
and have not been ported to the clean core yet.

It is kept for the same reason the rest of `legacy/` is: these are the drawings to **port from**, not
to delete and redraw. Many are good, and several will come back as their level is rebuilt.

    img/       the committed renders and hand-authored SVGs, straight out of `docs/img/`
    scripts/   the generators, straight out of `figures/scripts/`
    data/      the trees and other inputs some generators read, one directory per figure

## What stayed behind in the active tree

Only what something a reader reaches actually names:

| Figure | Named by |
| --- | --- |
| `manual/book/figures/fig-2-1-four-levels.svg` | Ch2, the docs home page, the README |
| `manual/book/figures/composition.svg`, `composition2.svg` | SPEC §7, as reference figures |
| `docs/img/age_crown.svg` | SPEC §7, as a reference figure |
| `docs/img/gillespie_{poisson,step,loop,everywhere}.svg` | Appendix A |
| `docs/img/event_levels.svg` | the genome resolution ladder — kept deliberately, not yet placed |
| `assets/logo.svg` | the README |

`figures/scripts/zombi_style.py` is the house style every generator imports, and `figures/STYLE.md`
describes it; both stayed.

## Porting one back

Move the script to `figures/scripts/` and its inputs to `figures/<name>/`, run it, and put the render
where the chapter that wants it can see it. Check it against SPEC §7 first — the figure conventions
changed in the rewrite, and most of these predate them.
