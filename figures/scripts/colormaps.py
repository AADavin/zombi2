"""A handful of hand-tabulated colormaps (anchor stops), so we don't need matplotlib.

Each palette is a list of (position 0..1, (r,g,b)) stops; ``cmap(name)`` returns a
function t in [0,1] -> (r,g,b). Sequential ones are perceptually ordered; the two
'diverging' ones (coolwarm, RdBu) pass through a light middle — nice for a trait
that ranges below/above a central value.
"""

from __future__ import annotations

CMAPS: dict[str, list] = {
    "viridis":  [(0.0, (68, 1, 84)), (0.25, (59, 82, 139)), (0.5, (33, 145, 140)), (0.75, (94, 201, 98)), (1.0, (253, 231, 37))],
    "magma":    [(0.0, (0, 0, 4)), (0.25, (81, 18, 124)), (0.5, (183, 55, 121)), (0.75, (252, 137, 97)), (1.0, (252, 253, 191))],
    "plasma":   [(0.0, (13, 8, 135)), (0.25, (126, 3, 168)), (0.5, (204, 71, 120)), (0.75, (248, 149, 64)), (1.0, (240, 249, 33))],
    "cividis":  [(0.0, (0, 34, 78)), (0.25, (58, 66, 101)), (0.5, (122, 124, 116)), (0.75, (183, 180, 112)), (1.0, (255, 234, 70))],
    "mako":     [(0.0, (11, 4, 10)), (0.25, (43, 44, 110)), (0.5, (51, 109, 150)), (0.75, (65, 176, 164)), (1.0, (222, 245, 229))],
    "rocket":   [(0.0, (3, 5, 26)), (0.25, (110, 28, 70)), (0.5, (203, 65, 64)), (0.75, (242, 153, 120)), (1.0, (250, 235, 221))],
    "coolwarm": [(0.0, (59, 76, 192)), (0.25, (123, 145, 235)), (0.5, (221, 221, 221)), (0.75, (229, 138, 116)), (1.0, (180, 4, 38))],
    "RdBu":     [(0.0, (33, 102, 172)), (0.25, (146, 197, 222)), (0.5, (247, 247, 247)), (0.75, (244, 165, 130)), (1.0, (178, 24, 43))],
}


def ramp(stops, t):
    t = max(0.0, min(1.0, t))
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return tuple(int(round(c0[i] + (c1[i] - c0[i]) * f)) for i in range(3))
    return stops[-1][1]


def cmap(name):
    stops = CMAPS[name]
    return lambda t: ramp(stops, t)


def hexc(rgb):
    return "#%02x%02x%02x" % tuple(rgb)
