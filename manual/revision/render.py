#!/usr/bin/env python
"""Render manual/revision/dashboard.yml to a single self-contained HTML page.

The YAML is the source of truth; this file only draws it. Re-run after editing
statuses so the page and the register can never disagree.

    python manual/revision/render.py            # -> manual/revision/dashboard.html
"""

from __future__ import annotations

import html
import pathlib
import re
import sys

import yaml

HERE = pathlib.Path(__file__).parent
DATA = HERE / "dashboard.yml"
OUT = HERE / "dashboard.html"

E = lambda s: html.escape(str(s if s is not None else ""))


def md(s: str) -> str:
    """The comment text carries `code`, **bold** and *italic* from the agents."""
    s = E(s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", s)
    return s


CSS = """
/* ---- tokens ------------------------------------------------------------
   Palette is the ZOMBI2 house style (figures/scripts/zombi_style.py): INK
   #1a1a1a, the muted teal STATE_ON #2f7d84 (which is also the docs site's
   teal primary), warm taupe STATE_OFF #b9b0a4, and the coupled-figure
   semantics COOCCUR/AVOID. Neutrals are biased a few degrees toward the
   teal so they read as chosen rather than inherited.                     */
:root {
  --ground:      #eef1f0;
  --panel-bg:    #f8faf9;
  --panel-line:  #d8ded9;
  --ink:         #1a1a1a;
  --ink-2:       #4a524f;
  --ink-3:       #7d857f;
  --accent:      #2f7d84;
  --accent-soft: #e2eeee;
  --taupe:       #b9b0a4;
  --ok:          #2f8f4e;
  --warn:        #c9762c;
  --crit:        #c0402f;
  --ok-bg:       #e6f1e8;
  --warn-bg:     #fbeedd;
  --crit-bg:     #f8e6e2;
  --shadow:      0 1px 2px rgba(26,26,26,.05), 0 4px 14px rgba(26,26,26,.045);
  --serif: "Iowan Old Style", "Charter", "Palatino Linotype", Palatino, Georgia, serif;
  --sans: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: ui-monospace, "SF Mono", SFMono-Regular, Menlo, monospace;
}
@media (prefers-color-scheme: dark) {
  :root {
    --ground:     #101413;
    --panel-bg:   #171c1b;
    --panel-line: #2a3230;
    --ink:        #e7ebe9;
    --ink-2:      #a6b0ac;
    --ink-3:      #79837f;
    --accent:     #63b3ba;
    --accent-soft:#1c2a2b;
    --taupe:      #8a8378;
    --ok:         #63b87c;
    --warn:       #e0a05c;
    --crit:       #e0705e;
    --ok-bg:      #16241a;
    --warn-bg:    #2a2013;
    --crit-bg:    #2b1815;
    --shadow:     0 1px 2px rgba(0,0,0,.3), 0 4px 14px rgba(0,0,0,.25);
  }
}
:root[data-theme="dark"] {
  --ground:#101413; --panel-bg:#171c1b; --panel-line:#2a3230; --ink:#e7ebe9;
  --ink-2:#a6b0ac; --ink-3:#79837f; --accent:#63b3ba; --accent-soft:#1c2a2b;
  --taupe:#8a8378; --ok:#63b87c; --warn:#e0a05c; --crit:#e0705e;
  --ok-bg:#16241a; --warn-bg:#2a2013; --crit-bg:#2b1815;
  --shadow:0 1px 2px rgba(0,0,0,.3), 0 4px 14px rgba(0,0,0,.25);
}
:root[data-theme="light"] {
  --ground:#eef1f0; --panel-bg:#f8faf9; --panel-line:#d8ded9; --ink:#1a1a1a;
  --ink-2:#4a524f; --ink-3:#7d857f; --accent:#2f7d84; --accent-soft:#e2eeee;
  --taupe:#b9b0a4; --ok:#2f8f4e; --warn:#c9762c; --crit:#c0402f;
  --ok-bg:#e6f1e8; --warn-bg:#fbeedd; --crit-bg:#f8e6e2;
  --shadow:0 1px 2px rgba(26,26,26,.05), 0 4px 14px rgba(26,26,26,.045);
}

* { box-sizing: border-box; }
body {
  margin: 0; background: var(--ground); color: var(--ink);
  font-family: var(--sans); font-size: 15px; line-height: 1.55;
  -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 1120px; margin: 0 auto; padding: 0 24px 96px; }
code { font-family: var(--mono); font-size: .88em; background: var(--accent-soft);
       padding: .08em .32em; border-radius: 3px; word-break: break-word; }
strong { font-weight: 650; }

/* ---- masthead ---- */
header.top { padding: 56px 0 30px; border-bottom: 1px solid var(--panel-line); }
.eyebrow { font-family: var(--mono); font-size: 11px; letter-spacing: .14em;
           text-transform: uppercase; color: var(--accent); margin: 0 0 14px; }
h1 { font-family: var(--serif); font-weight: 600; font-size: clamp(30px, 4.6vw, 46px);
     line-height: 1.1; margin: 0 0 6px; text-wrap: balance; letter-spacing: -.012em; }
.sub { font-family: var(--serif); font-size: 19px; font-style: italic;
       color: var(--ink-2); margin: 0 0 26px; }
.headline { font-size: 15.5px; line-height: 1.65; color: var(--ink-2);
            max-width: 66ch; border-left: 2px solid var(--accent); padding-left: 18px; }
.headline strong { color: var(--ink); }

/* ---- kpi strip ---- */
.kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
        gap: 1px; background: var(--panel-line); border: 1px solid var(--panel-line);
        border-radius: 8px; overflow: hidden; margin: 30px 0 8px; }
.kpi { background: var(--panel-bg); padding: 16px 18px; }
.kpi .n { font-family: var(--serif); font-size: 30px; font-weight: 600;
          font-variant-numeric: tabular-nums; line-height: 1; }
.kpi .l { font-size: 11.5px; color: var(--ink-3); margin-top: 6px; letter-spacing: .02em; }
.kpi.hot .n { color: var(--warn); }
.kpi.good .n { color: var(--ok); }

/* ---- sections ---- */
section { margin-top: 60px; scroll-margin-top: 20px; }
h2 { font-family: var(--serif); font-size: 27px; font-weight: 600; margin: 0 0 6px;
     letter-spacing: -.01em; }
.lede { color: var(--ink-2); max-width: 72ch; margin: 0 0 24px; font-size: 14.5px; }

/* ---- decision cards ---- */
.dcard { background: var(--panel-bg); border: 1px solid var(--panel-line);
         border-radius: 9px; margin-bottom: 12px; overflow: hidden; box-shadow: var(--shadow);
         border-left: 3px solid var(--taupe); }
.dcard[data-sev="foundational"] { border-left-color: var(--crit); }
.dcard[data-sev="major"]        { border-left-color: var(--warn); }
.dcard[data-sev="minor"]        { border-left-color: var(--taupe); }
.dhead { display: flex; align-items: flex-start; gap: 14px; padding: 15px 18px;
         cursor: pointer; }
.dhead:hover { background: var(--accent-soft); }
.dhead::-webkit-details-marker { display: none; }
.did { font-family: var(--mono); font-size: 13px; font-weight: 600; color: var(--accent);
       background: var(--accent-soft); border-radius: 4px; padding: 2px 7px;
       flex: 0 0 auto; margin-top: 1px; }
.dtitle { font-family: var(--serif); font-size: 17.5px; font-weight: 600; flex: 1;
          line-height: 1.35; }
.dq { display: block; color: var(--ink-2); font-size: 13.5px; margin-top: 5px;
      font-family: var(--sans); font-weight: 400; line-height: 1.5; }
.dmeta { display: flex; gap: 6px; flex: 0 0 auto; align-items: center; flex-wrap: wrap;
         justify-content: flex-end; max-width: 210px; }
.dbody { padding: 4px 18px 20px; border-top: 1px solid var(--panel-line); }
.dsec { margin-top: 16px; }
.dsec h4 { font-size: 11px; text-transform: uppercase; letter-spacing: .1em;
           color: var(--ink-3); margin: 0 0 7px; font-weight: 600; font-family: var(--mono); }
.opt { border: 1px solid var(--panel-line); border-radius: 6px; padding: 10px 12px;
       margin-bottom: 7px; background: var(--ground); }
.opt .ol { font-weight: 600; font-size: 14px; margin-bottom: 3px; }
.opt .oc { font-size: 12.5px; color: var(--ink-3); font-family: var(--mono); line-height: 1.5; }
.opt .op { font-size: 13px; color: var(--ink-2); }
.rec { background: var(--accent-soft); border-radius: 6px; padding: 13px 15px;
       font-size: 14px; line-height: 1.62; }
.rec::before { content: "Recommendation"; display: block; font-family: var(--mono);
               font-size: 10.5px; letter-spacing: .1em; text-transform: uppercase;
               color: var(--accent); margin-bottom: 6px; font-weight: 600; }
.rec.decided { background: var(--ok-bg); }
.rec.decided::before { content: "Decided · 17 Jul"; color: var(--ok); }
.pill.st-decided { color: var(--ok); background: var(--ok-bg); border-color: transparent; font-weight: 600; }
.dcard[data-status="decided"] > .dhead .dtitle,
.dcard[data-status="partly-decided"] > .dhead .dtitle { color: var(--ink-2); }
.dcard[data-status="decided"] { border-left-color: var(--ok); }
.blocks { display: flex; flex-wrap: wrap; gap: 4px; }

/* ---- pills ---- */
.pill { font-family: var(--mono); font-size: 10.5px; letter-spacing: .04em; padding: 2px 7px;
        border-radius: 99px; border: 1px solid var(--panel-line); color: var(--ink-3);
        white-space: nowrap; background: var(--ground); }
.pill.sev-foundational { color: var(--crit); background: var(--crit-bg); border-color: transparent; font-weight: 600; }
.pill.sev-major { color: var(--warn); background: var(--warn-bg); border-color: transparent; }
.pill.st-actionable { color: var(--ok); background: var(--ok-bg); border-color: transparent; }
.pill.st-needs-decision { color: var(--warn); background: var(--warn-bg); border-color: transparent; }
.pill.st-stale-claim, .pill.st-conflicts { color: var(--crit); background: var(--crit-bg); border-color: transparent; }
.pill.st-already-done { color: var(--ink-3); text-decoration: line-through; }
.pill.ref { color: var(--accent); background: var(--accent-soft); border-color: transparent;
            cursor: pointer; font-weight: 600; }
.pill.ref:hover { background: var(--accent); color: var(--panel-bg); }

/* ---- callouts (contradictions / stale) ---- */
.call { background: var(--panel-bg); border: 1px solid var(--panel-line); border-radius: 9px;
        padding: 0; margin-bottom: 10px; overflow: hidden; box-shadow: var(--shadow); }
.call.contra { border-left: 3px solid var(--warn); }
.call.stale  { border-left: 3px solid var(--crit); }
.call > summary { padding: 13px 18px; cursor: pointer; list-style: none;
                  display: flex; gap: 12px; align-items: baseline; }
.call > summary::-webkit-details-marker { display: none; }
.call > summary:hover { background: var(--accent-soft); }
.call .ct { font-family: var(--serif); font-size: 15.5px; font-weight: 600; flex: 1; line-height: 1.4; }
.call .cb { padding: 2px 18px 18px; font-size: 14px; line-height: 1.62; color: var(--ink-2); }
.call .cb .lab { font-family: var(--mono); font-size: 10.5px; letter-spacing: .09em;
                 text-transform: uppercase; color: var(--ink-3); display: block;
                 margin: 13px 0 4px; font-weight: 600; }
.call .cb .lab:first-child { margin-top: 4px; }
.between { font-family: var(--mono); font-size: 12px; color: var(--ink-3); margin-bottom: 8px; }
.between li { margin-bottom: 3px; }

/* ---- register ---- */
.controls { display: flex; gap: 8px; flex-wrap: wrap; align-items: center;
            margin-bottom: 14px; position: sticky; top: 0; background: var(--ground);
            padding: 12px 0; z-index: 5; border-bottom: 1px solid var(--panel-line); }
.controls select, .controls input {
  font-family: var(--sans); font-size: 13px; padding: 6px 9px; border-radius: 6px;
  border: 1px solid var(--panel-line); background: var(--panel-bg); color: var(--ink);
}
.controls input { flex: 1; min-width: 180px; }
.controls select:focus-visible, .controls input:focus-visible, .dhead:focus-visible,
.call > summary:focus-visible, .rowhead:focus-visible, button:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px;
}
.count { font-family: var(--mono); font-size: 12px; color: var(--ink-3); margin-left: auto; }
button.clear { font-family: var(--sans); font-size: 12.5px; padding: 6px 11px; border-radius: 6px;
               border: 1px solid var(--panel-line); background: var(--panel-bg);
               color: var(--ink-2); cursor: pointer; }
button.clear:hover { border-color: var(--accent); color: var(--accent); }

.row { background: var(--panel-bg); border: 1px solid var(--panel-line); border-radius: 8px;
       margin-bottom: 6px; overflow: hidden; }
.rowhead { display: grid; grid-template-columns: 62px 1fr auto; gap: 12px; align-items: baseline;
           padding: 11px 14px; cursor: pointer; list-style: none; }
.rowhead::-webkit-details-marker { display: none; }
.rowhead:hover { background: var(--accent-soft); }
.cid { font-family: var(--mono); font-size: 12.5px; font-weight: 600; color: var(--accent); }
.cask { font-size: 14px; line-height: 1.45; }
.cmeta { display: flex; gap: 5px; flex-wrap: wrap; justify-content: flex-end; }
.rowbody { padding: 4px 14px 16px 14px; border-top: 1px solid var(--panel-line);
           font-size: 13.5px; line-height: 1.62; color: var(--ink-2); }
.rowbody .lab { font-family: var(--mono); font-size: 10.5px; letter-spacing: .09em;
                text-transform: uppercase; color: var(--ink-3); display: block;
                margin: 12px 0 3px; font-weight: 600; }
.rowbody .lab:first-child { margin-top: 4px; }
.verb { font-family: var(--serif); font-style: italic; font-size: 14.5px; color: var(--ink);
        border-left: 2px solid var(--taupe); padding-left: 12px; white-space: pre-wrap; }
.anchor { font-family: var(--mono); font-size: 12px; color: var(--ink-3); word-break: break-all; }

/* ---- waves ---- */
.wave { display: grid; grid-template-columns: 30px 1fr; gap: 16px; }
.wave + .wave { margin-top: 0; }
.wrail { position: relative; }
.wrail::before { content: ""; position: absolute; left: 50%; top: 0; bottom: 0; width: 1px;
                 background: var(--panel-line); transform: translateX(-.5px); }
.wdot { position: relative; width: 11px; height: 11px; border-radius: 99px;
        background: var(--panel-bg); border: 2px solid var(--accent); margin: 15px auto 0; z-index: 1; }
.wave.w0 .wdot { background: var(--crit); border-color: var(--crit); }
.wcard { background: var(--panel-bg); border: 1px solid var(--panel-line); border-radius: 9px;
         padding: 14px 17px; margin-bottom: 8px; box-shadow: var(--shadow); }
.wtitle { font-family: var(--serif); font-size: 16.5px; font-weight: 600; margin-bottom: 2px; }
.wgoal { font-size: 13.5px; color: var(--ink-2); margin-bottom: 9px; }
.wcont { font-size: 13px; color: var(--ink-3); line-height: 1.6; }
.wblock { font-family: var(--mono); font-size: 11px; color: var(--ink-3); margin-top: 9px; }

.themegrid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 10px; }
.theme { background: var(--panel-bg); border: 1px solid var(--panel-line); border-radius: 9px;
         padding: 15px 17px; box-shadow: var(--shadow); }
.theme h3 { font-family: var(--serif); font-size: 16px; font-weight: 600; margin: 0 0 4px;
            line-height: 1.35; }
.theme .tn { font-family: var(--mono); font-size: 11px; color: var(--accent); margin-bottom: 8px; }
.theme p { margin: 0; font-size: 13px; line-height: 1.6; color: var(--ink-2); }

footer { margin-top: 70px; padding-top: 22px; border-top: 1px solid var(--panel-line);
         font-size: 12.5px; color: var(--ink-3); line-height: 1.7; }
.hide { display: none !important; }
@media (prefers-reduced-motion: reduce) { * { transition: none !important; animation: none !important; } }
@media (max-width: 720px) {
  .rowhead { grid-template-columns: 1fr; gap: 5px; }
  .cmeta { justify-content: flex-start; }
  .dhead { flex-wrap: wrap; }
  .dmeta { max-width: none; justify-content: flex-start; }
}
"""


def pill(text, cls=""):
    return f'<span class="pill {cls}">{E(text)}</span>'


def build(d: dict) -> str:
    m, dec, comments = d["meta"], d["decisions"], d["comments"]
    n = len(comments)
    n_dec_needed = sum(c["status"] == "needs-decision" for c in comments)
    n_blocked = sum(1 for c in comments if c["blocked_by"])
    n_free = sum(1 for c in comments if not c["blocked_by"] and c["status"] == "actionable")
    n_found = sum(1 for x in dec if x["severity"] == "foundational")

    P = []
    A = P.append

    # ---------- masthead ----------
    A('<div class="wrap"><header class="top">')
    A('<p class="eyebrow">Manual revision · round 1 · 17 July 2026</p>')
    A(f"<h1>{n} comments, {len(dec)} decisions</h1>")
    A('<p class="sub">Where the ZOMBI2 manual stands, before a word is rewritten.</p>')
    A(f'<div class="headline">{md(d["headline"])}</div>')
    A("</header>")

    A('<div class="kpis">')
    for num, lab, cls in [
        (n, "comments catalogued", ""),
        (len(dec), f"decisions ({n_found} foundational)", ""),
        (n_dec_needed, "blocked on you", "hot"),
        (n_blocked, "downstream of a decision", "hot"),
        (n_free, "free to start today", "good"),
    ]:
        A(f'<div class="kpi {cls}"><div class="n">{num}</div><div class="l">{E(lab)}</div></div>')
    A("</div>")

    # ---------- decisions ----------
    A('<section id="decisions"><h2>The decision spine</h2>')
    A(
        '<p class="lede">Ordered by how much each one blocks. These are the questions only you can '
        "answer, and every one of them is quoted by comments in chapters you have not reached yet — "
        "which is why they come before the rewrite, not during it. Each card opens to the options, "
        "their real migration cost, and a recommendation.</p>"
    )
    for x in dec:
        st = x.get("status", "open")
        settled = st in ("decided", "partly-decided")
        A(f'<details class="dcard" data-sev="{E(x["severity"])}" data-status="{E(st)}" id="{E(x["id"])}">')
        A('<summary class="dhead">')
        A(f'<span class="did">{E(x["id"])}</span>')
        A(f'<span class="dtitle">{md(x["title"])}<span class="dq">{md(x["question"])}</span></span>')
        A('<span class="dmeta">')
        if settled:
            A(pill("decided ✓" if st == "decided" else "partly decided", "st-decided"))
        A(pill(x["severity"], f'sev-{x["severity"]}'))
        A(pill(x["scope"]))
        A(pill(f'blocks {len(x["blocks"])}'))
        A("</span></summary>")
        A('<div class="dbody">')
        A(f'<div class="dsec"><h4>Why it gates</h4><div>{md(x["why_it_gates"])}</div></div>')
        if x.get("options"):
            A('<div class="dsec"><h4>Options</h4>')
            for o in x["options"]:
                A('<div class="opt">')
                A(f'<div class="ol">{md(o["label"])}</div>')
                if o.get("pros"):
                    A(f'<div class="op"><strong>+</strong> {md(o["pros"])}</div>')
                if o.get("cons"):
                    A(f'<div class="op"><strong>−</strong> {md(o["cons"])}</div>')
                A(f'<div class="oc">cost · {md(o["cost"])}</div>')
                A("</div>")
            A("</div>")
        if settled and x.get("decision"):
            A(f'<div class="dsec"><div class="rec decided">{md(x["decision"])}</div></div>')
        elif x.get("recommendation"):
            A(f'<div class="dsec"><div class="rec">{md(x["recommendation"])}</div></div>')
        if x["blocks"]:
            A(f'<div class="dsec"><h4>Blocks {len(x["blocks"])} comments</h4><div class="blocks">')
            for b in x["blocks"]:
                A(f'<span class="pill ref" data-jump="{E(b)}">{E(b)}</span>')
            A("</div></div>")
        A("</div></details>")
    A("</section>")

    # ---------- read first ----------
    A('<section id="before"><h2>Read before you decide</h2>')
    A(
        '<p class="lede">Two categories that are worth more than the rest of this page. '
        "<strong>Reversals</strong> are places where these comments undo something you ratified days ago, "
        "or where they contradict each other. <strong>Corrections</strong> are places where the premise "
        "of a comment is factually wrong. Several are wrong in your favour: the thing already works.</p>"
    )
    A(f"<h3 style='font-family:var(--serif);font-size:19px;margin:22px 0 10px'>Reversals and internal conflicts <span class='pill' style='vertical-align:middle'>{len(d['contradictions'])}</span></h3>")
    for c in d["contradictions"]:
        first = c["between"][0]
        A('<details class="call contra">')
        A(f'<summary><span class="ct">{md(first[:150])}{"…" if len(first) > 150 else ""}</span>{pill(f"{len(c['between'])}-way")}</summary>')
        A('<div class="cb"><span class="lab">The conflict</span><ul class="between">')
        for b in c["between"]:
            A(f"<li>{md(b)}</li>")
        A("</ul>")
        A(f'{md(c["explanation"])}<span class="lab">Resolution</span>{md(c["resolution"])}</div></details>')

    A(f"<h3 style='font-family:var(--serif);font-size:19px;margin:26px 0 10px'>Corrections to the premise <span class='pill' style='vertical-align:middle'>{len(d['stale_claims'])}</span></h3>")
    for c in d["stale_claims"]:
        A('<details class="call stale">')
        A(f'<summary><span class="did">{E(c["comment_id"])}</span><span class="ct">{md(c["claim"][:160])}{"…" if len(c["claim"]) > 160 else ""}</span></summary>')
        A(f'<div class="cb"><span class="lab">Actually</span>{md(c["truth"])}'
          f'<span class="lab">Evidence</span>{md(c["evidence"])}</div></details>')
    A("</section>")

    # ---------- themes ----------
    A('<section id="themes"><h2>What the comments are really about</h2>')
    A('<p class="lede">Eight patterns account for most of the 133. They are why chapter-by-chapter '
      "editing alone would not converge: the same defect surfaces in six chapters wearing six different "
      "faces.</p>")
    A('<div class="themegrid">')
    for t in d["themes"]:
        A('<div class="theme">')
        A(f'<h3>{md(t["name"])}</h3>')
        A(f'<div class="tn">{len(t["comment_ids"])} comments · {", ".join(E(i) for i in t["comment_ids"][:9])}</div>')
        A(f'<p>{md(t["description"])}</p></div>')
    A("</div></section>")

    # ---------- register ----------
    A('<section id="register"><h2>The register</h2>')
    A(
        '<p class="lede">Every comment you wrote, anchored to the line it is about and checked against '
        "the code. Nothing here is summarised away: this list is the contract that nothing gets dropped. "
        "Filter it, or click a decision id above to see what it blocks.</p>"
    )
    A('<div class="controls">')
    A('<input id="q" type="search" placeholder="Search comments, evidence, recommendations…" aria-label="Search the register">')
    chapters = sorted({c["chapter"] for c in comments}, key=lambda x: (x == "GLOBAL", int(x) if x.isdigit() else 0))
    A('<select id="fch" aria-label="Filter by chapter"><option value="">All chapters</option>')
    for ch in chapters:
        lab = "Global" if ch == "GLOBAL" else f"Chapter {ch}"
        A(f'<option value="{E(ch)}">{E(lab)}</option>')
    A("</select>")
    A('<select id="fst" aria-label="Filter by status"><option value="">Any status</option>')
    for st in ["actionable", "needs-decision", "stale-claim", "already-done", "conflicts"]:
        A(f'<option value="{E(st)}">{E(st)}</option>')
    A("</select>")
    A('<select id="fdf" aria-label="Filter by difficulty"><option value="">Any effort</option>')
    for df in ["trivial", "small", "medium", "large", "epic"]:
        A(f'<option value="{E(df)}">{E(df)}</option>')
    A("</select>")
    A('<select id="fd" aria-label="Filter by blocking decision"><option value="">Any decision</option>')
    A('<option value="__none__">Blocked by nothing</option>')
    for x in dec:
        A(f'<option value="{E(x["id"])}">{E(x["id"])} — {E(x["title"][:38])}</option>')
    A("</select>")
    A('<button class="clear" id="clear">Reset</button>')
    A('<span class="count" id="count"></span>')
    A("</div>")
    A('<div id="rows">')
    for c in comments:
        blocked = " ".join(c["blocked_by"]) or "__none__"
        hay = " ".join([c["id"], c["verbatim"], c["ask"], c["evidence"], c["recommendation"], c["anchor"]]).lower()
        A(f'<details class="row" id="{E(c["id"])}" data-ch="{E(c["chapter"])}" '
          f'data-st="{E(c["status"])}" data-df="{E(c["difficulty"])}" data-d="{E(blocked)}" '
          f'data-hay="{E(hay)}">')
        A('<summary class="rowhead">')
        A(f'<span class="cid">{E(c["id"])}</span>')
        A(f'<span class="cask">{md(c["ask"])}</span>')
        A('<span class="cmeta">')
        A(pill(c["status"], f'st-{c["status"]}'))
        A(pill(c["difficulty"]))
        A(pill(c["scope"]))
        for b in c["blocked_by"]:
            A(f'<span class="pill ref" data-jump="{E(b)}">{E(b)}</span>')
        A("</span></summary>")
        A('<div class="rowbody">')
        A(f'<span class="lab">You wrote · {E(c["kind"])}</span><div class="verb">{md(c["verbatim"])}</div>')
        A(f'<span class="lab">Where</span><span class="anchor">{md(c["anchor"])}</span>')
        A(f'<span class="lab">What we checked</span>{md(c["evidence"])}')
        A(f'<span class="lab">Recommendation</span>{md(c["recommendation"])}')
        if c["cross_refs"]:
            A('<span class="lab">Related</span>' + " · ".join(md(x) for x in c["cross_refs"]))
        A("</div></details>")
    A("</div></section>")

    # ---------- waves ----------
    A('<section id="waves"><h2>The order of work</h2>')
    A('<p class="lede">Wave 0 is four sessions with you; Wave 1 runs in parallel and is blocked on '
      "nothing, so it can start today. Everything after that waits on the spine. Chapter 2 is the "
      "bottleneck of the whole plan.</p>")
    for w in d["waves"]:
        w0 = "w0" if w["wave"].strip().lower().startswith("wave 0") else ""
        A(f'<div class="wave {w0}"><div class="wrail"><div class="wdot"></div></div><div class="wcard">')
        A(f'<div class="wtitle">{md(w["wave"])}</div>')
        A(f'<div class="wgoal">{md(w["goal"])}</div>')
        A(f'<div class="wcont">{md(w["contents"])}</div>')
        bb = ", ".join(w["blocked_by"]) if w["blocked_by"] else "nothing — can start now"
        A(f'<div class="wblock">blocked by · {E(bb)}</div>')
        A("</div></div>")
    A("</section>")

    A(
        "<footer>Generated from <code>manual/revision/dashboard.yml</code>, which is the source of truth — "
        "edit the YAML and re-run <code>python manual/revision/render.py</code>; do not hand-edit this page. "
        f'Every one of the {n} comments carries an id and must reach a terminal outcome, so none can be '
        "silently dropped. Evidence was gathered by 25 agents reading the chapters against the code at "
        f'<code>{E(m.get("generated", ""))}</code>.</footer>'
    )
    A("</div>")

    # ---------- js ----------
    A(
        """
<script>
(function () {
  var q = document.getElementById('q'), fch = document.getElementById('fch'),
      fst = document.getElementById('fst'), fdf = document.getElementById('fdf'),
      fd = document.getElementById('fd'), count = document.getElementById('count'),
      rows = Array.prototype.slice.call(document.querySelectorAll('#rows .row'));

  function apply() {
    var t = q.value.trim().toLowerCase(), ch = fch.value, st = fst.value,
        df = fdf.value, dv = fd.value, shown = 0;
    rows.forEach(function (r) {
      var ok = (!t || r.dataset.hay.indexOf(t) !== -1)
        && (!ch || r.dataset.ch === ch)
        && (!st || r.dataset.st === st)
        && (!df || r.dataset.df === df)
        && (!dv || (dv === '__none__' ? r.dataset.d === '__none__'
                                      : r.dataset.d.split(' ').indexOf(dv) !== -1));
      r.classList.toggle('hide', !ok);
      if (ok) shown++;
    });
    count.textContent = shown + ' of ' + rows.length;
  }
  [q, fch, fst, fdf, fd].forEach(function (el) {
    el.addEventListener('input', apply); el.addEventListener('change', apply);
  });
  document.getElementById('clear').addEventListener('click', function () {
    q.value = ''; fch.value = ''; fst.value = ''; fdf.value = ''; fd.value = ''; apply();
  });

  // A decision pill filters the register to what that decision blocks.
  document.addEventListener('click', function (e) {
    var p = e.target.closest('.pill.ref');
    if (!p) return;
    var id = p.dataset.jump;
    if (/^D\\d+$/.test(id)) {
      fd.value = id; q.value = ''; fch.value = ''; fst.value = ''; fdf.value = '';
      apply();
      document.getElementById('register').scrollIntoView({ behavior: 'smooth', block: 'start' });
    } else {
      var el = document.getElementById(id);
      if (el) { el.open = true; el.scrollIntoView({ behavior: 'smooth', block: 'center' }); }
    }
  });
  apply();
})();
</script>"""
    )
    return "\n".join(P)


def main() -> int:
    if not DATA.exists():
        print(f"missing {DATA}", file=sys.stderr)
        return 1
    d = yaml.safe_load(DATA.read_text())
    page = f"<title>ZOMBI2 manual revision · decision dashboard</title>\n<style>{CSS}</style>\n" + build(d)
    # Emit pure ASCII with numeric references for everything else. The page carries Σ, ×, ≤ and
    # arrows, and a host that serves it without `charset=utf-8` would mangle them into Latin-1
    # (Σ -> Î£). Entities are correct under any charset, so this is not worth gambling on.
    page = page.encode("ascii", "xmlcharrefreplace").decode("ascii")
    OUT.write_text(page, encoding="ascii")
    print(f"wrote {OUT}  ({len(page) / 1024:.0f} KB, {len(d['comments'])} comments, {len(d['decisions'])} decisions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
