"""Render every example and (re)build the published gallery page.

    cd gallery && python build.py

Renders figures/ (local, regenerated) and writes the page to ../web/gallery.html, which the site
deploy copies to the root so it publishes at /gallery.html. Adding a level = a new module with an
EXAMPLES list, added to LEVELS (each entry carries a URL slug used as the section anchor). One module
may feed more than one section: joining.py supplies both the conditioning and the joining lists.
"""

from __future__ import annotations

import inspect
import json
import os
import textwrap

from zombi2 import __version__ as __zombi2_version__

from PIL import Image

import genomes
import crosslevel
import joining
import sequences
import species
import traits

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
OUT = os.path.abspath(os.path.join(HERE, "..", "web", "gallery.html"))   # published at /gallery.html
#: Where the page's figures are written, beside it rather than inside it. They used to be embedded as
#: base64 data URIs, which made one self-contained file — and made every rebuild a fresh six-megabyte
#: blob in git, because base64 does not delta against its previous version. Fifty of those were the
#: bulk of a seventy-five-megabyte clone. As files, a rebuild costs only the figures that changed.
WEBFIG = os.path.abspath(os.path.join(HERE, "..", "web", "figures"))

#: Each section's reference prefix. An example is cited as `Ge3` — the section, then its position in
#: it. The number is **derived at build time**, never stored: the id (`genome_synteny_tree`) is the
#: identity, so an example can be inserted anywhere and every number after it simply follows on the
#: next build. Two letters because species and sequences both start with one S.
PREFIX = {"species": "Sp", "genomes": "Ge", "sequences": "Sq",
          "traits": "Tr", "conditioning": "Co", "joining": "Jo"}

# (slug, title, blurb, examples) — the slug is the section id the landing-page cards link to.
# The first section is the one that starts open; every other section starts folded.
LEVELS = [
    ("species", "Species trees",
     "Forward birth-death trees. The run keeps the whole history, survivors and extinctions, "
     "and the diversification model shows on the tree.",
     species.EXAMPLES),
    ("genomes", "Genomes",
     "Genes on chromosomes. A genome draws as a ring, and two genomes show their synteny. "
     "Gene-family events and copy number read against the species tree.",
     genomes.EXAMPLES),
    ("sequences", "Sequences",
     "The dated tree the sequences evolve down, and an alignment lined up row-for-row with its tips.",
     sequences.EXAMPLES),
    ("traits", "Traits",
     "A trait evolving down the tree. Branches take the colour of its value. Some examples add "
     "a companion panel.",
     traits.EXAMPLES),
    ("conditioning", "Conditioning",
     "Two runs, in order. The first run grows the driver on the tree and holds it fixed. The "
     "second run reads it. A driver is a trait, a gene family or a whole module. It drives a "
     "rate, or which lineage receives a transfer.",
     joining.CONDITIONING + crosslevel.EXAMPLES),
    ("joining", "Joining",
     "One run makes both. The trait sets the speciation or extinction rate of the lineage "
     "carrying it. The trait and the tree therefore come out together.",
     joining.JOINING),
]


def render_all():
    os.makedirs(FIGDIR, exist_ok=True)
    for _, _, _, examples in LEVELS:
        for ex in examples:
            ex.render(os.path.join(FIGDIR, f"{ex.id}.png"))
            print("rendered", ex.id)


#: The widest an embedded figure is stored at. The detail view and the click-to-zoom lightbox read
#: the same data URI as the card, so this is also the resolution anyone zooming in actually gets —
#: at 920 the conditioning figures went to mush on a tree with 45 tips.
MAXW = 1200


def _figure_url(ex, maxw=MAXW):
    """Write the card's figure to ``web/figures/`` and return the page-relative URL for it."""
    im = Image.open(os.path.join(FIGDIR, f"{ex.id}.png")).convert("RGB")
    if im.width > maxw:
        im = im.resize((maxw, round(im.height * maxw / im.width)), Image.LANCZOS)
    os.makedirs(WEBFIG, exist_ok=True)
    im.save(os.path.join(WEBFIG, f"{ex.id}.png"), format="PNG", optimize=True)
    return f"figures/{ex.id}.png"


def _code_for(ex):
    """The snippet shown on the detail view: the example's curated ``code`` if given, else the body
    of its render function (def line dropped, dedented)."""
    if ex.code:
        return ex.code
    src = textwrap.dedent(inspect.getsource(ex.render)).splitlines()
    i = 0
    while i < len(src) and not src[i].lstrip().startswith("def "):
        i += 1
    return textwrap.dedent("\n".join(src[i + 1:])).strip()


def _detail_data(examples, store, slug):
    for i, ex in enumerate(examples, start=1):
        store[ex.id] = {"num": f"{PREFIX[slug]}{i}", "title": ex.title, "caption": ex.caption,
                        "tag": ex.tag, "code": _code_for(ex)}


def _cards(examples, slug):
    out = []
    for i, ex in enumerate(examples, start=1):
        uri = _figure_url(ex)
        num = f"{PREFIX[slug]}{i}"
        out.append(f"""      <figure class="card" tabindex="0" role="button" data-id="{ex.id}" aria-label="Open: {ex.title}">
        <div class="thumb"><img loading="lazy" src="{uri}" alt="{ex.title}"></div>
        <figcaption><h3><span class="num">{num}</span>{ex.title}</h3><p>{ex.caption}</p><span class="tag">{ex.tag}</span></figcaption>
      </figure>""")
    return "\n".join(out)


def build_html():
    sections, detail = [], {}
    for i, (slug, name, blurb, examples) in enumerate(LEVELS):
        _detail_data(examples, detail, slug)
        # <details>/<summary>: folding works with JavaScript off, and the count and blurb stay in the
        # summary, so a folded section still says what is inside it. Only the first section opens.
        sections.append(f"""  <details class="level" id="{slug}"{" open" if i == 0 else ""}>
    <summary class="level-head"><h2>{name}</h2><span class="count">{len(examples)}</span><span class="blurb">{blurb}</span></summary>
    <div class="grid">
{_cards(examples, slug)}
    </div>
  </details>""")
    data = "<script>window.EX = " + json.dumps(detail) + ";</script>"
    # A minimal standards-mode head so the page stands on its own at /gallery.html.
    head = ('<!doctype html>\n<meta charset="utf-8">\n'
            '<title>ZOMBI2 — Examples gallery</title>\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n')
    html = head + _CSS + _PAGE_OPEN + "\n".join(sections) + _PAGE_CLOSE + data + _JS
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write(html)
    print(f"wrote {os.path.relpath(OUT, HERE)} ({len(html) // 1024} KB)")


_CSS = """<style>
:root{--bg:#f4f8f7;--surface:#fff;--ink:#16211f;--muted:#57655f;--faint:#7c8a85;--line:#e2ebe8;
 --accent:#0d7d74;--accent-ink:#0a625b;--mat:#fff;--code-bg:#eef4f2;--code-ink:#1d2a27;--shadow:0 1px 2px rgba(18,40,36,.05),0 10px 28px -14px rgba(18,40,36,.16);--radius:13px}
@media (prefers-color-scheme:dark){:root{--bg:#0d1513;--surface:#141e1b;--ink:#e9f0ee;--muted:#8ea09b;--faint:#748681;
 --line:#243330;--accent:#45bcae;--accent-ink:#63cec2;--mat:#f8fbfa;--code-bg:#0b120f;--code-ink:#cfe0db;--shadow:0 1px 2px rgba(0,0,0,.45),0 12px 32px -16px rgba(0,0,0,.7)}}
:root[data-theme="light"]{--bg:#f4f8f7;--surface:#fff;--ink:#16211f;--muted:#57655f;--faint:#7c8a85;--line:#e2ebe8;
 --accent:#0d7d74;--accent-ink:#0a625b;--mat:#fff;--code-bg:#eef4f2;--code-ink:#1d2a27;--shadow:0 1px 2px rgba(18,40,36,.05),0 10px 28px -14px rgba(18,40,36,.16)}
:root[data-theme="dark"]{--bg:#0d1513;--surface:#141e1b;--ink:#e9f0ee;--muted:#8ea09b;--faint:#748681;--line:#243330;
 --accent:#45bcae;--accent-ink:#63cec2;--mat:#f8fbfa;--code-bg:#0b120f;--code-ink:#cfe0db;--shadow:0 1px 2px rgba(0,0,0,.45),0 12px 32px -16px rgba(0,0,0,.7)}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;-webkit-font-smoothing:antialiased;transition:background .25s,color .25s}
.wrap{max-width:1160px;margin:0 auto;padding:36px 24px 88px}
.top{display:flex;align-items:center;justify-content:space-between;gap:16px;margin-bottom:34px}
.eyebrow{font:600 12px/1 ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:.14em;text-transform:uppercase;color:var(--accent-ink);display:flex;align-items:center;gap:9px}
.eyebrow::before{content:"";width:22px;height:2px;background:var(--accent);border-radius:2px}
.toggle{appearance:none;border:1px solid var(--line);background:var(--surface);color:var(--muted);font:600 12px/1 system-ui;padding:9px 13px;border-radius:9px;cursor:pointer;display:flex;align-items:center;gap:7px;transition:.18s}
.toggle:hover{color:var(--ink);border-color:var(--accent)}
.toggle:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.masthead h1{font-size:clamp(2rem,4.4vw,2.7rem);line-height:1.06;letter-spacing:-.02em;margin:0 0 12px;text-wrap:balance;font-weight:700}
.masthead .lede{margin:0;max-width:60ch;color:var(--muted);font-size:1.05rem}
.masthead .lede b{color:var(--ink);font-weight:600}
.masthead .gallery-note{margin-top:14px}
.masthead .gallery-install{display:inline-block;margin:10px 0 0;padding:8px 14px;border:1px solid var(--line);border-radius:8px;background:var(--code-bg);font-family:var(--mono,ui-monospace,monospace);font-size:.92rem;color:var(--code-ink)}
.note{margin:18px 0 0;display:inline-flex;gap:9px;align-items:center;font-size:.82rem;color:var(--faint);border:1px dashed var(--line);border-radius:8px;padding:7px 12px}
.note b{color:var(--accent-ink);font-weight:600}
.level{margin-top:52px}
summary.level-head{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap;padding-bottom:14px;border-bottom:1px solid var(--line);cursor:pointer;list-style:none;-webkit-tap-highlight-color:transparent;transition:border-color .18s}
summary.level-head::-webkit-details-marker{display:none}
summary.level-head::marker{content:""}
summary.level-head:hover{border-color:var(--accent)}
summary.level-head:focus-visible{outline:2px solid var(--accent);outline-offset:4px;border-radius:4px}
.level-head h2{margin:0;font-size:1.4rem;letter-spacing:-.01em;font-weight:650}
.level-head .count{font:600 12px/1 ui-monospace,monospace;color:var(--faint);padding:4px 9px;border:1px solid var(--line);border-radius:20px}
.level-head .blurb{margin:0;color:var(--muted);font-size:.95rem;flex:1;min-width:200px}
summary.level-head::after{content:"";flex:none;align-self:center;width:9px;height:9px;margin-left:4px;border-right:2px solid var(--faint);border-bottom:2px solid var(--faint);transform:rotate(45deg) translate(-2px,-2px);transition:transform .2s}
.level[open]>summary.level-head::after{transform:rotate(-135deg) translate(-2px,-2px)}
.level>.grid{margin-top:24px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(315px,1fr));gap:22px}
.card{margin:0;background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);overflow:hidden;cursor:pointer;transition:transform .2s,box-shadow .2s,border-color .2s;display:flex;flex-direction:column}
.card:hover,.card:focus-visible{transform:translateY(-3px);box-shadow:var(--shadow);border-color:var(--accent);outline:none}
.card:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.thumb{background:var(--mat);aspect-ratio:16/11;display:flex;align-items:center;justify-content:center;overflow:hidden;border-bottom:1px solid var(--line)}
.thumb img{width:100%;height:100%;object-fit:contain;display:block}
figcaption{padding:16px 17px 17px;display:flex;flex-direction:column;gap:6px}
figcaption h3{margin:0;font-size:1.03rem;font-weight:640;letter-spacing:-.01em}
figcaption p{margin:0;color:var(--muted);font-size:.88rem;line-height:1.5}
.num{display:inline-block;min-width:2.6em;color:var(--accent);font:600 .78rem/1 ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:.02em;vertical-align:.08em}
.version{font:500 .5em/1 ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--faint);vertical-align:.5em;letter-spacing:.02em}
.tag{margin-top:4px;font:600 11px/1 ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:.04em;color:var(--accent-ink);text-transform:lowercase}
footer{margin-top:64px;padding-top:22px;border-top:1px solid var(--line);color:var(--faint);font-size:.85rem;display:flex;justify-content:space-between;flex-wrap:wrap;gap:10px}
footer code{font:600 .82rem ui-monospace,monospace;color:var(--muted)}
.detail{position:fixed;inset:0;background:rgba(8,16,14,.86);backdrop-filter:blur(3px);display:none;align-items:flex-start;justify-content:center;padding:4vmin;z-index:50;overflow:auto}
.detail.open{display:flex}
.sheet{position:relative;background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);max-width:1000px;width:100%;box-shadow:0 24px 70px -20px rgba(0,0,0,.7);overflow:hidden;margin:auto}
.det-close{position:absolute;top:12px;right:12px;z-index:2;appearance:none;border:1px solid var(--line);background:var(--surface);color:var(--muted);width:34px;height:34px;border-radius:9px;cursor:pointer;font-size:15px;line-height:1;transition:.15s}
.det-close:hover{color:var(--ink);border-color:var(--accent)}
.det-fig{background:var(--mat);display:flex;align-items:center;justify-content:center;padding:22px;border-bottom:1px solid var(--line)}
.det-fig img{max-width:100%;max-height:58vh;object-fit:contain;display:block;cursor:zoom-in}
.det-fig{position:relative}
.det-zoomhint{position:absolute;bottom:10px;right:12px;font:600 10px/1 ui-monospace,monospace;letter-spacing:.06em;text-transform:uppercase;color:var(--faint);background:var(--surface);border:1px solid var(--line);border-radius:6px;padding:5px 8px;pointer-events:none;opacity:.85}
.det-meta{padding:20px 26px 4px}
.det-meta .tag{margin:0}
.det-meta h2{margin:.4rem 0 .5rem;font-size:1.4rem;letter-spacing:-.01em;font-weight:650}
.det-meta p{margin:0;color:var(--muted);font-size:.96rem;line-height:1.55;max-width:72ch}
.det-codewrap{padding:16px 26px 26px}
.det-codehead{display:flex;justify-content:space-between;align-items:center;font:600 11px/1 ui-monospace,monospace;letter-spacing:.1em;text-transform:uppercase;color:var(--faint);margin-bottom:9px}
.det-codehead button{appearance:none;border:1px solid var(--line);background:var(--surface);color:var(--muted);font:600 11px/1 ui-monospace,monospace;padding:6px 11px;border-radius:7px;cursor:pointer;text-transform:none;letter-spacing:0;transition:.15s}
.det-codehead button:hover{color:var(--ink);border-color:var(--accent)}
.det-code{margin:0;background:#0f1714;border:1px solid #20302c;border-radius:10px;padding:16px 18px;overflow-x:auto;font:12.5px/1.7 ui-monospace,SFMono-Regular,Menlo,monospace;color:#cdd9d4;white-space:pre;tab-size:4}
.hl-sec{color:#59d3c3;font-weight:700}
.hl-c{color:#6c8079;font-style:italic}
.hl-s{color:#93d2a1}
.hl-n{color:#e3ab60}
.hl-k{color:#6cb7d7}
.lightbox{position:fixed;inset:0;background:rgba(6,12,10,.95);display:none;z-index:80;overflow:auto;cursor:zoom-out}
.lightbox.open{display:block}
.lightbox img{display:block;margin:auto;max-width:100%;max-height:100vh;object-fit:contain}
.lightbox.actual{cursor:grab}
.lightbox.actual img{max-width:none;max-height:none;margin:24px auto;cursor:zoom-out}
.lb-close{position:fixed;top:14px;right:16px;z-index:81;appearance:none;border:1px solid rgba(255,255,255,.25);background:rgba(0,0,0,.35);color:#fff;width:38px;height:38px;border-radius:9px;cursor:pointer;font-size:16px;line-height:1}
.lb-close:hover{background:rgba(0,0,0,.6)}
@media (prefers-reduced-motion:reduce){*{transition:none!important}}
</style>"""

_PAGE_OPEN = """
<div class="wrap">
  <div class="top">
    <div class="eyebrow"><a href="./">ZOMBI2</a> · examples</div>
    <button class="toggle" id="toggle" aria-label="Toggle colour theme">◑ <span>Dark</span></button>
  </div>
  <header class="masthead">
    <h1>Examples gallery <span class="version">ZOMBI2 @VERSION@</span></h1>
    <p class="lede gallery-note">Click any figure for the code that produces it. Every example
      simulates with ZOMBI2 and plots with <a href="https://pypi.org/project/phylustrator/">Phylustrator</a>,
      a separate package &mdash; so to run one you need both:</p>
    <pre class="gallery-install">pip install zombi2 phylustrator</pre>
  </header>
""".replace("@VERSION@", __zombi2_version__)

_PAGE_CLOSE = """
  <footer>
    <span><a href="./">ZOMBI2</a> · <a href="docs/">docs</a> · <a href="https://github.com/AADavin/zombi2">GitHub</a></span>
    <span><code>plot(tree) + color_branches(…) + …</code></span>
  </footer>
</div>
<div class="detail" id="detail" aria-hidden="true">
  <div class="sheet" role="dialog" aria-modal="true" aria-labelledby="det-title">
    <button class="det-close" id="det-close" aria-label="Close">✕</button>
    <div class="det-fig"><img id="det-img" alt=""><span class="det-zoomhint">click to zoom</span></div>
    <div class="det-meta">
      <span class="tag" id="det-tag"></span>
      <h2 id="det-title"></h2>
      <p id="det-cap"></p>
    </div>
    <div class="det-codewrap">
      <div class="det-codehead"><span>example code</span><button id="det-copy">copy</button></div>
      <pre class="det-code"><code id="det-code"></code></pre>
    </div>
  </div>
</div>
<div class="lightbox" id="lightbox" aria-hidden="true">
  <button class="lb-close" id="lb-close" aria-label="Close zoom">✕</button>
  <img id="lb-img" alt="">
</div>
"""

_JS = """<script>
(function(){
  var root=document.documentElement, KEY="zombi2-gallery-theme";
  var saved=localStorage.getItem(KEY); if(saved) root.setAttribute("data-theme",saved);
  function cur(){return root.getAttribute("data-theme")||(matchMedia("(prefers-color-scheme:dark)").matches?"dark":"light");}
  var btn=document.getElementById("toggle");
  function label(){btn.querySelector("span").textContent=cur()==="dark"?"Light":"Dark";}
  label();
  btn.addEventListener("click",function(){var n=cur()==="dark"?"light":"dark";root.setAttribute("data-theme",n);localStorage.setItem(KEY,n);label();});
  var EX=window.EX||{};
  function esc(s){return s.replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");}
  var HLRE=/(###[^\\n]*)|(#[^\\n]*)|("[^"]*"|'[^']*')|\\b(\\d+(?:\\.\\d+)?)\\b|\\b(import|from|as|for|in|if|else|elif|def|return|and|or|not|None|True|False|lambda|with|print)\\b/g;
  function hl(code){
    var out="",last=0,m;
    while((m=HLRE.exec(code))){
      out+=esc(code.slice(last,m.index));
      var cls=m[1]?"sec":m[2]?"c":m[3]?"s":m[4]?"n":"k";
      out+='<span class="hl-'+cls+'">'+esc(m[0])+"</span>";
      last=m.index+m[0].length;
    }
    return out+esc(code.slice(last));
  }
  var det=document.getElementById("detail"),
      dImg=document.getElementById("det-img"), dTag=document.getElementById("det-tag"),
      dTitle=document.getElementById("det-title"), dCap=document.getElementById("det-cap"),
      dCode=document.getElementById("det-code"), dCopy=document.getElementById("det-copy");
  function open(c){
    var id=c.getAttribute("data-id"), meta=EX[id]||{}, img=c.querySelector("img");
    dImg.src=img.src; dImg.alt=img.alt;
    dTag.innerHTML=meta.tag||""; dTitle.textContent=(meta.num?meta.num+"  ":"")+(meta.title||img.alt);
    dCap.innerHTML=meta.caption||""; dCode.innerHTML=hl(meta.code||"");
    dCopy.textContent="copy";
    det.classList.add("open"); det.setAttribute("aria-hidden","false"); det.scrollTop=0;
  }
  function close(){det.classList.remove("open");det.setAttribute("aria-hidden","true");dImg.src="";}
  document.querySelectorAll(".card").forEach(function(c){
    c.addEventListener("click",function(){open(c);});
    c.addEventListener("keydown",function(e){if(e.key==="Enter"||e.key===" "){e.preventDefault();open(c);}});
  });
  document.getElementById("det-close").addEventListener("click",close);
  det.addEventListener("click",function(e){if(e.target===det)close();});
  dCopy.addEventListener("click",function(){
    navigator.clipboard&&navigator.clipboard.writeText(dCode.textContent).then(function(){
      dCopy.textContent="copied";setTimeout(function(){dCopy.textContent="copy";},1400);});
  });

  // --- click-to-zoom lightbox: fit-to-screen, click again for actual size (scrollable) ---
  var lb=document.getElementById("lightbox"), lbImg=document.getElementById("lb-img");
  function openLB(){lbImg.src=dImg.src;lbImg.alt=dImg.alt;lb.classList.remove("actual");
    lb.classList.add("open");lb.setAttribute("aria-hidden","false");lb.scrollTop=0;}
  function closeLB(){lb.classList.remove("open","actual");lb.setAttribute("aria-hidden","true");lbImg.src="";}
  dImg.addEventListener("click",openLB);
  lbImg.addEventListener("click",function(e){e.stopPropagation();  // toggle fit <-> actual size
    lb.classList.toggle("actual");lb.scrollTop=0;lb.scrollLeft=0;});
  lb.addEventListener("click",closeLB);   // click the backdrop closes
  document.getElementById("lb-close").addEventListener("click",function(e){e.stopPropagation();closeLB();});
  document.addEventListener("keydown",function(e){
    if(e.key!=="Escape")return;
    if(lb.classList.contains("open"))closeLB(); else close();   // zoom first, then the detail sheet
  });

  // --- deep links: gallery.html#genomes must unfold a section that starts closed ---
  function openHash(){
    var id=(location.hash||"").slice(1); if(!id)return;
    var el=document.getElementById(id); if(!el)return;
    var d=el.closest("details"); if(d)d.open=true;              // the section itself, or one above it
    el.scrollIntoView();
  }
  openHash();
  window.addEventListener("hashchange",openHash);
})();
</script>"""


if __name__ == "__main__":
    render_all()
    build_html()
