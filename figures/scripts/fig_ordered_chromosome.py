"""Figure: one circular chromosome as a signed gene order.

Chapter 5 (ordered genomes) reads a leaf as ``[ 0+ 1+ 3+ 3+ 4− ]`` — a run of gene tokens, each a
family with a strand. This draws that same chromosome as the ring it is: five genes spaced evenly
around a circle, each an arrow that points **the way its strand reads**. Colour is by family, so the
two copies of family ``3`` — a tandem duplication — share a colour and sit adjacent, and family ``4``,
left backwards by an inversion, is the one arrow pointing against the others.

The genome is the leaf ``n2`` of the chapter's first example (``seed=2``); the figure is generated
from that run, not drawn by hand, so it cannot drift from what the code produces.

House style: near-black ink, ColorBrewer identities, no title inside the figure (the manual captions
it).

Run:  python figures/scripts/fig_ordered_chromosome.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered
from zombi_style import save, FONT, INK, MUTED, FS_LABEL

# families are told apart by shade, not hue — a dark-to-light grey ramp. Two genes of the same
# family share a shade, so a tandem duplication reads as two neighbours of one grey.
GREYS = ["#2b2b2b", "#575757", "#838383", "#afafaf"]


def _leaf_order(node_label: str = "n2"):
    """The (family, strand) order of one extant leaf of the chapter's first example."""
    tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=2)
    g = simulate_genomes_ordered(tree, duplication=0.3, loss=0.2, origination=0.15, inversion=0.5,
                                 chromosomes=1, initial_families=5, seed=2)
    target = int(node_label[1:])
    chrom = g.genomes[target][0]
    return [(gene.family, gene.strand) for gene in chrom.genes]


W, H = 760, 720
CX, CY = W / 2, H / 2 - 6
R = 205                       # ring radius (gene centre-line)
GENE_W = 42                   # gene arc thickness
GAP_DEG = 15                  # blank arc between neighbouring genes
HEAD = 15                     # arrowhead length in degrees of extra reach beyond the arc
HEAD_HALF = GENE_W * 0.78     # half-width of the arrowhead base (a touch wider than the body)


def _xy(angle_deg: float, radius: float) -> tuple[float, float]:
    """A point at ``angle_deg`` clockwise from 12 o'clock, on ``radius``."""
    a = math.radians(angle_deg)
    return CX + radius * math.sin(a), CY - radius * math.cos(a)


def _gene_arc(d, a0: float, a1: float, strand: int, colour: str) -> None:
    """One gene as a thick arc from ``a0`` to ``a1`` (clockwise degrees), arrowhead at the end the
    strand reads toward: the high-angle end for ``+1``, the low-angle end for ``−1``."""
    lead, tail = (a1, a0) if strand == 1 else (a0, a1)
    body_lead = lead - math.copysign(HEAD, lead - tail)      # leave room for the head
    # the arc body (a stroked path along the ring)
    x0, y0 = _xy(tail, R)
    x1, y1 = _xy(body_lead, R)
    large = 1 if abs(body_lead - tail) > 180 else 0
    sweep = 1 if body_lead > tail else 0
    d.append(draw.Path(f"M {x0:.2f} {y0:.2f} A {R} {R} 0 {large} {sweep} {x1:.2f} {y1:.2f}",
                       stroke=colour, stroke_width=GENE_W, fill="none", stroke_linecap="butt"))
    # the arrowhead: apex on the ring at `lead`, base a short chord just wider than the body
    apex = _xy(lead, R)
    br = _xy(body_lead, R + HEAD_HALF)
    bl = _xy(body_lead, R - HEAD_HALF)
    d.append(draw.Lines(apex[0], apex[1], br[0], br[1], bl[0], bl[1], close=True, fill=colour))


def render() -> None:
    order = _leaf_order("n2")
    n = len(order)
    families = sorted({fam for fam, _ in order})
    colour_of = {fam: GREYS[i % len(GREYS)] for i, fam in enumerate(families)}

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    # the backbone ring, faint, so the chromosome reads as a closed loop behind the genes
    d.append(draw.Circle(CX, CY, R, stroke=MUTED, stroke_width=1.4, fill="none",
                         stroke_dasharray="2,5"))

    slot = 360.0 / n
    for i, (fam, strand) in enumerate(order):
        a0 = i * slot + GAP_DEG / 2
        a1 = (i + 1) * slot - GAP_DEG / 2
        _gene_arc(d, a0, a1, strand, colour_of[fam])
        # the family label, outside the ring at the gene's mid-angle
        lx, ly = _xy((a0 + a1) / 2, R + GENE_W / 2 + 34)
        d.append(draw.Text(str(fam), FS_LABEL, lx, ly, font_family=FONT, fill=INK,
                           font_weight="bold", text_anchor="middle", dominant_baseline="central"))

    save(d, "ordered_chromosome")


if __name__ == "__main__":
    render()
