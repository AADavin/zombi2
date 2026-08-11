# ZOMBI2 examples gallery

A scikit-learn-style **examples gallery** for ZOMBI2. Each example **simulates** with ZOMBI2 and
**plots** with [Phylustrator](https://github.com/AADavin/Phylustrator) (trees via `ph.trees`, genomes
/ synteny / alignments via `ph.genomes`). Six sections: Species · Genomes · Sequences · Traits ·
Conditioning · Joining.

The built page is published at **`/gallery.html`** (the site deploy copies `web/` to the root), and
the landing-page level cards link to its per-level sections (`gallery.html#species`, `#genomes`, …).

## Build it

```bash
pip install -e ".[gallery]"      # phylustrator + cairosvg + matplotlib + pillow
cd gallery && python build.py     # renders figures/ and writes ../web/gallery.html
```

`figures/` (and its cached `_data/` runs, including a one-off NCBI GFF download for the *Mycoplasma*
example) is regenerated and git-ignored; the committed artefact is `web/gallery.html`. `BW` at the
top of [`helpers.py`](helpers.py) is the branch width `helpers.style()` uses; most examples set their
own `branch_width` in their `ph.Style(...)`.

## Layout

| File | Role |
|---|---|
| [`helpers.py`](helpers.py) | shared: `Example`, the house `style`, the extinct-dashing handoff, the cached CLI runs, and the **matplotlib companion-panel** compositing (`composite_below` / `composite_beside` / `composite_markov` / `composite_model_realization`) |
| [`species.py`](species.py), [`genomes.py`](genomes.py), [`sequences.py`](sequences.py), [`traits.py`](traits.py), [`joining.py`](joining.py) | one module per level — each exposes an `EXAMPLES` list of `Example(id, title, caption, tag, render, code)`. `joining.py` splits its list into `CONDITIONING` and `JOINING`, which are two sections |
| [`build.py`](build.py) | renders every example and writes `web/gallery.html` from the `LEVELS` registry (each entry carries the URL slug used as its section anchor) |

**Add an example:** write a `render(out_png)` function in the level module and append an `Example` to
its list. **Add a level:** new module with an `EXAMPLES` list, add it to `LEVELS` in `build.py` (with a
slug) — and add a matching card + link on the landing page (`web/index.html`).

## Design decisions (settled)

- **Six sections:** Species · Genomes · Sequences · Traits · **Conditioning** · **Joining**. The last
  two are relations rather than levels — conditioning grows the driver first and holds it fixed, joining
  grows both in one run. Validation **case studies** (RED, synteny) stay a separate group under
  [`analyses/`](../analyses) — those are validations, not demos.
- **Title outside the figure** (the caption carries it — no in-figure title), sklearn-style.
- **Mixed rendering:** Phylustrator for trees and genomes; matplotlib for companion panels (LTT
  skyline, Markov-chain insets, scatter), composited on a precise shared time axis (`px(t) == t`; see
  `helpers._extent`).
- Species trees show **extinct lineages dashed** via `plot(tree, dashed=…)`.
