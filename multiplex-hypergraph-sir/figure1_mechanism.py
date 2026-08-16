"""
Figure 1 (schematic): the mechanism behind the group closure.

Conventions follow the reference figures: thin black strokes, open circles for
susceptible nodes and filled for infectious, one red accent reserved for
transmission, serif italic symbols, and no prose inside the drawing -- the
annotation is mathematical and the explanation belongs to the caption.

Hyperedges are thin circular contours, solid for layer 1 and dashed for
layer 2, arranged so that groups genuinely share the focal node.

  (a) the multiplex hypergraph and the layer degree;
  (b) the intra-group cascade, labelled by the (i,s) counts the recursion uses;
  (c) the group-level branching step and the excess subtraction.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.lines import Line2D

BLK, GRY, LGY, RED = "#141414", "#7d7d7d", "#b9b9b9", "#a8291f"
F_I, F_S, F_R = "#8f8f8f", "#ffffff", "#ececec"

mpl.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 360,
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9,
    "text.color": BLK,
})

RN, LWC, LWN = 0.092, 0.8, 0.8


def group(ax, c, r, m, a0, dashed=False, color=BLK, lw=LWC, z=1):
    ax.add_patch(Circle(c, r, facecolor="none", edgecolor=color, lw=lw,
                        ls=(0, (3.2, 2.0)) if dashed else "-", zorder=z))
    return [np.array([c[0] + r * np.cos(a0 + 2 * np.pi * j / m),
                      c[1] + r * np.sin(a0 + 2 * np.pi * j / m)])
            for j in range(m)]


def node(ax, xy, fill=F_S, r=RN, lw=LWN, z=5, ec=BLK):
    ax.add_patch(Circle(xy, r, facecolor=fill, edgecolor=ec, lw=lw, zorder=z))


def arrow(ax, a, b, color=BLK, lw=0.95, ms=7.5, z=6, rad=0.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=color, zorder=z,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=5.5, shrinkB=5.5))


def panel_label(ax, letter, title):
    ax.text(0.0, 1.0, f"{letter}. ", transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")
    ax.text(0.035, 1.0, title, transform=ax.transAxes, fontsize=10,
            va="top", ha="left")


fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(12.4, 3.55))
fig.subplots_adjust(left=0.012, right=0.988, bottom=0.03, top=0.90, wspace=0.10)
for ax in (axA, axB, axC):
    ax.set_aspect("equal"); ax.axis("off")
for x in (0.3455, 0.6725):                      # thin rules between panels
    fig.add_artist(Line2D([x, x], [0.06, 0.92], color=LGY, lw=0.7))

# ================================================================ (a)
rg = 0.62
cA = rg * np.array([np.cos(np.deg2rad(150)), np.sin(np.deg2rad(150))])
cB = rg * np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])
cC = rg * np.array([0.0, -1.0])
gA = group(axA, cA, rg, 3, np.deg2rad(-30))
gB = group(axA, cB, rg, 3, np.deg2rad(210))
gC = group(axA, cC, rg, 4, np.deg2rad(90), dashed=True)
for g in (gA, gB, gC):
    for p in g[1:]:
        node(axA, p)
node(axA, np.zeros(2), fill=F_I, r=0.112, lw=1.15)
axA.text(0.0, -0.245, r"$u$", ha="center", va="top", fontsize=10,
         bbox=dict(fc="white", ec="none", pad=0.8), zorder=8)
axA.text(*(cA + [0, 0.30]), r"$e_1$", ha="center", fontsize=9, color=GRY)
axA.text(*(cB + [0, 0.30]), r"$e_2$", ha="center", fontsize=9, color=GRY)
axA.text(*(cC + [0.34, -0.20]), r"$e_3$", ha="center", fontsize=9, color=GRY)

axA.text(-1.30, 1.16, r"layer 1,  $m_1=3$", fontsize=9, ha="left")
axA.text(-1.30, -1.46, r"layer 2,  $m_2=4$", fontsize=9, ha="left")
axA.text(1.34, 0.24, r"$k^{(1)}(u)=2$", fontsize=9.5, ha="left")
axA.text(1.34, -0.06, r"$k^{(2)}(u)=1$", fontsize=9.5, ha="left")
axA.set_xlim(-1.55, 2.55); axA.set_ylim(-1.62, 1.42)
panel_label(axA, "a", "multiplex hypergraph")

# ================================================================ (b)
rb, yb = 0.455, 0.34


def snap(ax, cx, n_inf):
    pts = group(ax, (cx, yb), rb, 4, np.deg2rad(90))
    for j, p in enumerate(pts):
        node(ax, p, fill=F_I if j < n_inf else F_S)
    return pts


p1, p2, p3 = snap(axB, -1.28, 1), snap(axB, 0.10, 2), snap(axB, 1.48, 3)
arrow(axB, p1[0], p1[1], color=RED)
arrow(axB, p2[1], p2[2], color=RED)
for cx, tri in ((-1.28, "(1,3)"), (0.10, "(2,2)"), (1.48, "(3,1)")):
    axB.text(cx, yb - 0.72, rf"$(i,s)={tri}$", ha="center", va="top",
             fontsize=9, color=BLK)
axB.text(-0.59, yb + 0.05, r"$\lambda s$", fontsize=9, color=RED, ha="center")
axB.text(0.79, yb + 0.05, r"$\lambda s$", fontsize=9, color=RED, ha="center")

ya = -1.05
axB.plot([-1.28, 0.42], [ya, ya], "-", color=BLK, lw=1.3, solid_capstyle="butt")
axB.plot([0.42, 1.86], [ya, ya], "-", color=RED, lw=1.3, solid_capstyle="butt")
for x in (-1.28, 0.42, 1.86):
    axB.plot([x, x], [ya - 0.075, ya + 0.075], "-", color=BLK, lw=0.85)
axB.text(-0.43, ya + 0.13, r"$\mathrm{Exp}(\mu)$", fontsize=8.6, color=GRY,
         ha="center")
axB.text(1.14, ya + 0.13, "extension", fontsize=8.6, color=RED, ha="center")
axB.text(0.29, ya - 0.30, "active period of the group", fontsize=8.6,
         color=GRY, ha="center", va="top")

axB.text(0.10, 1.16, r"$C(m,\lambda,\theta)\;>\;(m-1)\,T$", ha="center",
         fontsize=10.5)
axB.set_xlim(-1.90, 2.10); axB.set_ylim(-1.62, 1.46)
panel_label(axB, "b", "intra-group cascade")

# ================================================================ (c)
rc = 0.58
cIn = rc * np.array([-1.0, 0.0])
c1 = rc * np.array([np.cos(np.deg2rad(58)), np.sin(np.deg2rad(58))])
c2 = rc * np.array([np.cos(np.deg2rad(-58)), np.sin(np.deg2rad(-58))])
gIn = group(axC, cIn, rc, 3, np.deg2rad(0), color=GRY)
g1 = group(axC, c1, rc, 3, np.deg2rad(58 - 180))
g2 = group(axC, c2, rc, 4, np.deg2rad(-58 - 180), dashed=True)
for p in gIn[1:]:
    node(axC, p, fill=F_R, ec=GRY, lw=0.75)
for p in g1[1:] + g2[1:]:
    node(axC, p)
arrow(axC, gIn[1], np.zeros(2), color=RED, lw=1.05)
node(axC, np.zeros(2), fill=F_I, r=0.112, lw=1.15)

axC.plot([cIn[0] - 0.095, cIn[0] + 0.095], [-0.095, 0.095], "-",
         color=GRY, lw=1.25, zorder=7)
axC.plot([cIn[0] - 0.095, cIn[0] + 0.095], [0.095, -0.095], "-",
         color=GRY, lw=1.25, zorder=7)
axC.text(cIn[0], -0.80, r"$-\,\delta_{ab}$", fontsize=10, color=GRY,
         ha="center", va="center")
axC.text(*(c1 + 0.80 * c1 / np.linalg.norm(c1)), r"$C$", fontsize=10.5,
         ha="center", va="center", color=RED)
axC.text(*(c2 + 0.80 * c2 / np.linalg.norm(c2)), r"$C$", fontsize=10.5,
         ha="center", va="center", color=RED)
axC.text(*(c1 + [0.05, 0.30]), r"layer $b$", fontsize=8.6, color=GRY, ha="center")
axC.text(*(c2 + [0.05, -0.34]), r"layer $b'$", fontsize=8.6, color=GRY, ha="center")
axC.text(cIn[0], 0.80, r"layer $a$", fontsize=8.6, color=GRY, ha="center")
axC.set_xlim(-1.85, 1.85); axC.set_ylim(-1.62, 1.46)
panel_label(axC, "c", "group-level branching")

fig.savefig("figure1_mechanism.pdf", bbox_inches="tight")
fig.savefig("figure1_mechanism.png", bbox_inches="tight", dpi=210)
print("saved figure1_mechanism.pdf / .png")
