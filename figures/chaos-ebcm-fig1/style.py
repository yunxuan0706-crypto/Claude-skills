"""Shared CSF/Elsevier-style drawing helpers for Fig. 1."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Polygon
from matplotlib.path import Path

# ---------------------------------------------------------------- typography
plt.rcParams.update({
    "font.family": "Liberation Sans",          # Arial-metric, Elsevier house look
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
    "pdf.fonttype": 42,                        # embed as TrueType, keep text selectable
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

FS_PANEL = 9.0     # (a) (b) (c) (d)
FS_SUB   = 7.5     # panel subtitle
FS_LBL   = 7.2     # in-panel labels
FS_MATH  = 7.8     # display-ish math
FS_SM    = 7.0     # smallest text used anywhere

# ---------------------------------------------------------------- palette
# PALETTE=full  three hues + warm/cool layer tints (richest)
# PALETTE=duo   two hues, one neutral layer tint      (recommended for CSF)
# PALETTE=mono  single red accent, everything else neutral
_P = __import__("os").environ.get("PALETTE", "duo")

if _P == "full":
    C_S, C_S_E = "#A7C8E8", "#2F5C8A"
    C_I, C_I_E = "#E2726A", "#A32E26"
    C_R, C_R_E = "#4A4E52", "#1E2124"
    _F1, _A1, _F2, _A2 = "#5B7FA6", 0.085, "#B08A3E", 0.10
    C_TRAN, C_COUP = "#C0392B", "#00757F"
elif _P == "mono":
    C_S, C_S_E = "#C6C6C6", "#6E6E6E"
    C_I, C_I_E = "#C0392B", "#8C2A20"
    C_R, C_R_E = "#8E9296", "#3D4145"
    _F1, _A1, _F2, _A2 = "#8A8A8A", 0.06, "#8A8A8A", 0.06
    C_TRAN, C_COUP = "#C0392B", "#3D3D3D"
else:                                    # duo
    C_S, C_S_E = "#9DC3E6", "#2E5F8F"
    C_I, C_I_E = "#DE6B63", "#973630"
    C_R, C_R_E = "#43474A", "#1B1E20"
    _F1, _A1, _F2, _A2 = "#7A8794", 0.07, "#7A8794", 0.07
    C_TRAN, C_COUP = "#C0392B", "#3D3D3D"

C_U, C_U_E = "#FFFFFF", "#111111"
C_INK, C_MUTE = "#2B2B2B", "#9AA0A6"
L1 = dict(ec=C_INK, fc=_F1, alpha_f=_A1, ls="-")
L2 = dict(ec=C_INK, fc=_F2, alpha_f=_A2, ls=(0, (3, 2)))

NR = 1.55          # node radius (mm)

# Every drawing helper appends to REC so the same layout code can be replayed
# into a PowerPoint deck (see export_spec) as well as onto the matplotlib canvas.
REC = []
def _rec(**kw): REC.append(kw); return kw

# ---------------------------------------------------------------- geometry
def _hull(P):
    """Andrew monotone chain -> convex hull, counter-clockwise."""
    pts = sorted(set(map(tuple, np.asarray(P, float).tolist())))
    if len(pts) <= 2:
        return np.array(pts, float)
    def half(seq):
        st = []
        for q in seq:
            while len(st) >= 2:
                a = np.subtract(st[-1], st[-2]); b = np.subtract(q, st[-1])
                if a[0] * b[1] - a[1] * b[0] <= 0: st.pop()
                else: break
            st.append(q)
        return st
    lo, up = half(pts), half(pts[::-1])
    return np.array(lo[:-1] + up[:-1], float)

def _offset(P, pad, n_arc=18):
    """Minkowski sum of a convex polygon with a disk: rounded convex blob."""
    P = np.asarray(P, float); n = len(P)
    if n == 1:
        t = np.linspace(0, 2 * np.pi, 72)
        return P[0] + pad * np.c_[np.cos(t), np.sin(t)]
    out = []
    for i in range(n):
        p, prv, nxt = P[i], P[(i - 1) % n], P[(i + 1) % n]
        din, dout = p - prv, nxt - p
        din = din / (np.linalg.norm(din) or 1.0)
        dout = dout / (np.linalg.norm(dout) or 1.0)
        nin = np.array([din[1], -din[0]]); nout = np.array([dout[1], -dout[0]])
        a0 = np.arctan2(nin[1], nin[0]); a1 = np.arctan2(nout[1], nout[0])
        while a1 < a0 - 1e-9: a1 += 2 * np.pi
        for t in np.linspace(a0, a1, n_arc):
            out.append(p + pad * np.array([np.cos(t), np.sin(t)]))
    return np.array(out)

def blob(ax, pts, layer, pad=4.2, lw=0.9, z=1, alpha_scale=1.0):
    """One hyperedge, drawn as a rounded convex blob enclosing its member nodes."""
    curve = _offset(_hull(pts), pad)
    ax.add_patch(PathPatch(Path(curve, closed=True), facecolor=layer["fc"],
                           alpha=layer["alpha_f"] * alpha_scale, edgecolor="none", zorder=z))
    ax.add_patch(PathPatch(Path(curve, closed=True), facecolor="none",
                           edgecolor=layer["ec"], lw=lw, ls=layer["ls"],
                           zorder=z + 0.5, alpha=alpha_scale))
    _rec(kind="blob", pts=curve.tolist(), fill=layer["fc"],
         alpha=layer["alpha_f"] * alpha_scale, edge=layer["ec"], lw=lw,
         dashed=layer["ls"] != "-", edge_alpha=alpha_scale)
    return curve

def cross(ax, xy, s=1.6, color=C_TRAN, lw=1.3, z=8):
    x, y = xy
    ax.plot([x - s, x + s], [y - s, y + s], color=color, lw=lw, zorder=z, solid_capstyle="round")
    ax.plot([x - s, x + s], [y + s, y - s], color=color, lw=lw, zorder=z, solid_capstyle="round")
    _rec(kind="cross", x=x, y=y, s=s, color=color, lw=lw)

def rule(ax, x0, y0, x1, y1, color="#DCDCDC", lw=0.6, dashed=False, z=0):
    ax.plot([x0, x1], [y0, y1], color=color, lw=lw, zorder=z,
            ls=(0, (3, 2)) if dashed else "-")
    _rec(kind="rule", p0=(x0, y0), p1=(x1, y1), color=color, lw=lw, dashed=dashed)

def roundbox(ax, x, y, w, h, r=1.2, fc="#F2F4F6", ec="#B9C0C7", lw=0.7, z=3):
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x, y), w, h, zorder=z, facecolor=fc, edgecolor=ec, lw=lw,
                                boxstyle=f"round,pad=0,rounding_size={r}"))
    _rec(kind="roundbox", x=x, y=y, w=w, h=h, r=r, fill=fc, edge=ec, lw=lw)

def curve(ax, a, b, rad, color="#4A4A4A", lw=0.8, z=4, ms=5):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms, lw=lw,
                                 color=color, zorder=z,
                                 connectionstyle=f"arc3,rad={rad}"))
    ctl = ((a[0] + b[0]) / 2 + rad * (b[1] - a[1]),
           (a[1] + b[1]) / 2 - rad * (b[0] - a[0]))
    _rec(kind="arrow", p0=tuple(a), p1=tuple(b), ctl=ctl, color=color, lw=lw,
         dashed=False, alpha=1.0)

def node(ax, xy, state="S", r=NR, z=6, lw=0.75, alpha=1.0):
    fc, ec = {"S": (C_S, C_S_E), "I": (C_I, C_I_E),
              "R": (C_R, C_R_E), "U": (C_U, C_U_E)}[state]
    ax.add_patch(Circle(xy, r, facecolor=fc, edgecolor=ec,
                        lw=1.5 if state == "U" else lw, zorder=z, alpha=alpha))
    _rec(kind="node", x=xy[0], y=xy[1], r=r, fill=fc, edge=ec,
         lw=1.5 if state == "U" else lw, alpha=alpha)

def arrow(ax, a, b, color=C_TRAN, lw=0.85, z=4, rad=0.0, ls="-",
          sa=NR + 0.5, sb=NR + 1.4, ms=5.5, alpha=1.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms,
                                 shrinkA=sa * 2.83, shrinkB=sb * 2.83,
                                 lw=lw, color=color, zorder=z, alpha=alpha,
                                 linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}"))
    (x1, y1), (x2, y2) = a, b
    dx, dy = x2 - x1, y2 - y1
    L = (dx * dx + dy * dy) ** 0.5 or 1.0
    ux, uy = dx / L, dy / L
    p0 = (x1 + ux * sa, y1 + uy * sa)
    p1 = (x2 - ux * sb, y2 - uy * sb)
    ctl = None
    if rad:                       # matplotlib Arc3 control point
        ctl = ((p0[0] + p1[0]) / 2 + rad * (p1[1] - p0[1]),
               (p0[1] + p1[1]) / 2 - rad * (p1[0] - p0[0]))
    _rec(kind="arrow", p0=p0, p1=p1, ctl=ctl, color=color, lw=lw,
         dashed=ls != "-", alpha=alpha)

def txt(ax, x, y, s, size=FS_LBL, ha="left", va="center", color="#1A1A1A",
        weight="normal", style="normal", z=8, **kw):
    t = ax.text(x, y, s, fontsize=size, ha=ha, va=va, color=color,
                fontweight=weight, fontstyle=style, zorder=z, **kw)
    _rec(kind="text", x=x, y=y, s=s, size=size, ha=ha, va=va, color=color,
         weight=weight, style=style, _artist=t)
    return t

def panel_head(ax, x, y, letter, title):
    txt(ax, x, y, f"({letter})", size=FS_PANEL, weight="bold", va="top")
    txt(ax, x + 5.2, y, title, size=FS_SUB, va="top", color="#333333")


# ------------------------------------------------- PowerPoint spec export
def export_spec(path, fig, ax, W, H):
    """Dump REC as JSON with y measured from the top, so the same layout can be
    replayed as native PowerPoint shapes. Text extents are measured off the real
    render, which is what lets the deck place text boxes exactly."""
    import json
    import mathrun
    fig.canvas.draw()
    inv = ax.transData.inverted()
    out = []
    for r in REC:
        d = {k: v for k, v in r.items() if k != "_artist"}
        if d["kind"] == "text":
            bb = r["_artist"].get_window_extent(fig.canvas.get_renderer())
            (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
            d["w"], d["h"] = abs(x1 - x0), abs(y1 - y0)
            d["runs"] = mathrun.runs(d.pop("s"))
            d["y"] = H - d["y"]
        elif d["kind"] == "ellipse":
            d["cy"] = H - d["cy"]
            d["ang"] = -d["ang"]                  # y flips, so rotation mirrors
        elif d["kind"] == "blob":
            d["pts"] = [[x, H - y] for x, y in d["pts"]]
        elif d["kind"] in ("node", "cross"):
            d["y"] = H - d["y"]
        elif d["kind"] == "roundbox":
            d["y"] = H - d["y"] - d["h"]
        else:                                     # arrow, rule
            for k in ("p0", "p1", "ctl"):
                if d.get(k) is not None:
                    d[k] = [d[k][0], H - d[k][1]]
        out.append(d)
    json.dump({"w_mm": W, "h_mm": H, "shapes": out}, open(path, "w"), indent=1)
    print("wrote %s  (%d shapes)" % (path, len(out)))


# ---------------------------------- hyperedges as plain ellipses
# Everything here maps to one standard PowerPoint shape: Insert > Shapes > Oval,
# then set Size and Rotation. Nodes sit on the ellipse boundary.
def on_ell(c, a, b, ang, t):
    """Point at parameter t (deg) on an ellipse centred at c, rotated ang (deg)."""
    import math
    tr, ar = math.radians(t), math.radians(ang)
    x, y = a * math.cos(tr), b * math.sin(tr)
    return (c[0] + x * math.cos(ar) - y * math.sin(ar),
            c[1] + x * math.sin(ar) + y * math.cos(ar))


def hyperedge(ax, c, a, b, ang, ts, layer, lw=0.9, z=1, alpha_scale=1.0):
    """One hyperedge drawn as an ellipse; returns member positions at params ts."""
    from matplotlib.patches import Ellipse
    for fc, ec, al, zz in ((layer["fc"], "none", layer["alpha_f"] * alpha_scale, z),
                           ("none", layer["ec"], alpha_scale, z + 0.5)):
        ax.add_patch(Ellipse(c, 2 * a, 2 * b, angle=ang, facecolor=fc, edgecolor=ec,
                             lw=lw, ls=layer["ls"], alpha=al, zorder=zz))
    _rec(kind="ellipse", cx=c[0], cy=c[1], a=a, b=b, ang=ang, fill=layer["fc"],
         alpha=layer["alpha_f"] * alpha_scale, edge=layer["ec"], lw=lw,
         dashed=layer["ls"] != "-", edge_alpha=alpha_scale)
    return [on_ell(c, a, b, ang, t) for t in ts]
