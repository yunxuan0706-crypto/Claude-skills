#!/usr/bin/env python3
"""Schematic model figure for the multiplex-hypergraph SIR write-up.

Faithful vector redraw of the three-panel hand sketch:
  (a) a multiplex hypergraph with a shared cavity node u;
  (b) the intra-group cascade, contrasting the naive count (m-1)T with the
      cascade count C, and the extension of the group's active period;
  (c) group-level branching, showing the reciprocity subtraction -delta_ab that
      applies only within the same layer (a=b).

English labels throughout. STIX mathtext so the maths matches figs 1-2.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..", "figures")

DARK = "#555555"      # infected member / hub node
LIGHT = "#c9c9c9"     # source node in panel (c)
RED = "#b0302f"       # transmission arrow
GREY = "#8a8a8a"      # secondary labels

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "DejaVu Serif"],
    "font.size": 9,
    "mathtext.fontset": "stix",
    "mathtext.rm": "serif",
})


def node(ax, xy, kind="S", r=0.17, z=5):
    fc = {"S": "white", "I": DARK, "src": LIGHT}[kind]
    ax.add_patch(Circle(xy, r, facecolor=fc, edgecolor="black", lw=1.0, zorder=z))


def group(ax, c, r, dashed=False, z=2):
    ax.add_patch(Circle(c, r, fill=False, edgecolor="black", lw=1.2,
                        ls="--" if dashed else "-", zorder=z))


def arrow(ax, p0, p1, rad=0.0, z=6, lw=1.4, shrink=7):
    ax.add_patch(FancyArrowPatch(p0, p1, connectionstyle=f"arc3,rad={rad}",
                                 arrowstyle="-|>", mutation_scale=11,
                                 lw=lw, color=RED, shrinkA=shrink, shrinkB=shrink,
                                 zorder=z))


def ring(c, r, angles):
    return [(c[0] + r * math.cos(math.radians(a)),
             c[1] + r * math.sin(math.radians(a))) for a in angles]


# ------------------------------------------------------------- panel (a) ----
def panel_a(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    e1, e2, e3 = (3.0, 6.3), (6.0, 6.3), (4.5, 3.6)
    R = 1.8
    group(ax, e1, R); group(ax, e2, R); group(ax, e3, R, dashed=True)
    u = (4.5, 5.31)
    # e1 members: u + two others
    for p in ring(e1, R, [150, 205]): node(ax, p)
    # e2 members
    for p in ring(e2, R, [30, 335]): node(ax, p)
    # e3 members
    for p in ring(e3, R, [212, 270, 328]): node(ax, p)
    node(ax, u, "I", r=0.20, z=8)
    ax.text(0.2, 9.2, r"layer 1,  $m_1=3$", fontsize=9)
    ax.text(0.2, 0.9, r"layer 2,  $m_2=4$", fontsize=9)
    ax.text(2.2, 6.3, r"$e_1$", color=GREY, ha="center")
    ax.text(5.4, 6.6, r"$e_2$", color=GREY, ha="center")
    ax.text(3.6, 3.3, r"$e_3$", color=GREY, ha="center")
    ax.text(4.5, 4.72, r"$u$", ha="center", fontsize=10)
    ax.text(8.1, 5.7, r"$k^{(1)}(u)=2$", fontsize=9, ha="left")
    ax.text(8.1, 4.9, r"$k^{(2)}(u)=1$", fontsize=9, ha="left")


# ------------------------------------------------------------- panel (b) ----
def small_group(ax, c, states, arrows, r=0.85):
    group(ax, c, r)
    pos = ring(c, r, [90, 0, -90, 180])   # top, right, bottom, left
    for p, st in zip(pos, states):
        node(ax, p, st, r=0.15)
    for i, j, rad in arrows:
        arrow(ax, pos[i], pos[j], rad=rad, lw=1.2, shrink=5)
    return pos


def panel_b(ax):
    ax.set_xlim(0, 12); ax.set_ylim(0, 10); ax.axis("off")
    # naive (top)
    small_group(ax, (2.6, 8.5), ["I", "S", "S", "S"],
                [(0, 1, 0), (0, 2, 0), (0, 3, 0)])
    ax.text(1.0, 8.5, "naive", color=GREY, ha="right", va="center")
    ax.text(4.4, 8.5, r"$\Rightarrow\ (m-1)\,T$", va="center", fontsize=9.5)
    # cascade (middle): three snapshots, adjacent spread
    small_group(ax, (2.6, 5.4), ["I", "S", "S", "S"], [(0, 3, 0.2)])
    small_group(ax, (5.6, 5.4), ["I", "S", "S", "I"], [(3, 2, 0.2)])
    small_group(ax, (8.6, 5.4), ["I", "S", "I", "I"], [])
    ax.text(1.0, 5.4, "cascade", color=GREY, ha="right", va="center")
    ax.text(10.4, 5.4, r"$\Rightarrow\ C$", va="center", fontsize=9.5)
    for x, lab in ((2.6, "(1,3)"), (5.6, "(2,2)"), (8.6, "(3,1)")):
        ax.text(x, 4.1, r"$(i,s)=%s$" % lab, ha="center", fontsize=8.5)
    # active-period bars
    y1, y2 = 3.1, 2.5
    ax.text(2.2, y1, "naive", color=GREY, ha="right", va="center", fontsize=8.5)
    ax.text(2.2, y2, "cascade", color=GREY, ha="right", va="center", fontsize=8.5)
    ax.plot([2.5, 5.6], [y1, y1], color="black", lw=1.3, solid_capstyle="butt")
    for x in (2.5, 5.6): ax.plot([x, x], [y1-0.12, y1+0.12], color="black", lw=1.1)
    ax.plot([2.5, 5.6], [y2, y2], color="black", lw=1.3, solid_capstyle="butt")
    ax.plot([5.6, 8.9], [y2, y2], color=RED, lw=1.3, solid_capstyle="butt")
    for x in (2.5, 8.9): ax.plot([x, x], [y2-0.12, y2+0.12], color="black", lw=1.1)
    ax.text(7.25, 2.85, "extension", color=RED, ha="center", fontsize=8.5)
    ax.text(5.6, 1.8, "active period of the group", color=GREY, ha="center",
            fontsize=8.5)
    ax.text(6.0, 0.8, r"$C(m,\lambda,\theta)\ >\ (m-1)\,T$", ha="center",
            fontsize=10)


# ------------------------------------------------------------- panel (c) ----
def panel_c(ax):
    ax.set_xlim(0, 11); ax.set_ylim(0, 10); ax.axis("off")
    src_grp = (2.1, 5.2)
    ax.add_patch(Circle(src_grp, 1.45, fill=False, edgecolor=GREY, lw=1.2, zorder=2))
    node(ax, (2.9, 5.9), "src", r=0.18, z=6)
    node(ax, (1.4, 4.4), "src", r=0.18, z=6)
    # the big cross: this source group is the arrival route, not re-counted
    ax.plot([1.7, 2.5], [4.8, 5.6], color="black", lw=2.2, zorder=7)
    ax.plot([1.7, 2.5], [5.6, 4.8], color="black", lw=2.2, zorder=7)
    hub = (5.0, 5.0)
    node(ax, hub, "I", r=0.22, z=8)
    arrow(ax, (2.9, 5.9), hub, rad=-0.06)
    # two offspring groups, left of their labels
    ga, gb = (6.9, 7.0), (6.9, 3.0)
    group(ax, ga, 1.15)
    group(ax, gb, 1.15, dashed=True)
    for p in ring(ga, 1.15, [30, -30]): node(ax, p, "S", r=0.14)
    for p in ring(gb, 1.15, [30, -30]): node(ax, p, "S", r=0.14)
    arrow(ax, hub, (6.1, 6.4), rad=-0.1)
    arrow(ax, hub, (6.1, 3.6), rad=0.1)
    ax.text(6.9, 8.5, r"layer $a$", color=GREY, ha="center", fontsize=9)
    ax.text(8.5, 7.0, r"$C\,(X_{aa}-1)$", ha="left", fontsize=9)
    ax.text(6.9, 1.5, r"layer $b$", color=GREY, ha="center", fontsize=9)
    ax.text(8.5, 3.0, r"$C\,X_{ab}$", ha="left", fontsize=9)
    ax.text(2.1, 7.0, r"layer $a$", color=GREY, ha="center", fontsize=9)
    ax.text(1.6, 3.1, r"$-\,\delta_{ab}\ \ (a=b\ \mathrm{only})$", ha="center",
            fontsize=8.5)
    ax.text(5.2, 0.6, r"$N_{ab}=C\,(X_{ab}-\delta_{ab})$", ha="center",
            fontsize=10)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.1))
    panel_a(axes[0]); panel_b(axes[1]); panel_c(axes[2])
    titles = ["multiplex hypergraph", "intra-group cascade", "group-level branching"]
    for ax, letter, tt in zip(axes, "abc", titles):
        ax.set_aspect("equal")
        ax.text(0.0, 1.02, f"$\\bf{{{letter}}}$  {tt}", transform=ax.transAxes,
                fontsize=10, va="bottom")
    # thin separators between panels
    for x in (0.35, 0.665):
        fig.add_artist(plt.Line2D([x, x], [0.05, 0.95], color="#cccccc",
                                  lw=0.8, transform=fig.transFigure))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.03, wspace=0.12)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig_model.{ext}"), dpi=600,
                    bbox_inches="tight")
    print("wrote figures/fig_model.{pdf,png}")
    plt.close(fig)


if __name__ == "__main__":
    main()
