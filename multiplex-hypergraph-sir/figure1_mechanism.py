"""
Figure 1 (schematic): the mechanism behind the group closure.

Drawn in the line-and-circle idiom of the reference figures: thin black
strokes, open circles for susceptible nodes, filled for infectious, a single
red accent reserved for transmission events, and no fills or shading. Groups
(hyperedges) are thin closed contours -- solid for layer 1, dashed for layer 2.

  (a) what a multiplex hypergraph is, and what the layer degree counts;
  (b) the intra-group cascade -- why C(m,lambda,theta) exceeds (m-1)T;
  (c) the group-level branching process, with the excess subtraction.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

BLK, GRY, RED = "#1a1a1a", "#8c8c8c", "#b5342a"
FILL_I, FILL_S, FILL_R = "#9a9a9a", "#ffffff", "#e6e6e6"

mpl.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 340,
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
    "text.color": BLK,
})

R_NODE = 0.105


def group(ax, c, r, m, u_angle, dashed=False, lw=0.9, color=BLK):
    """Hyperedge as a thin circular contour; returns its m member positions,
    the first at angle u_angle, so groups can be made to share a node."""
    ax.add_patch(Circle(c, r, facecolor="none", edgecolor=color, lw=lw,
                        ls=(0, (3.5, 2.2)) if dashed else "-", zorder=1))
    return [np.array([c[0] + r * np.cos(u_angle + 2 * np.pi * j / m),
                      c[1] + r * np.sin(u_angle + 2 * np.pi * j / m)])
            for j in range(m)]


def node(ax, xy, fill=FILL_S, r=R_NODE, lw=0.9, z=5):
    ax.add_patch(Circle(xy, r, facecolor=fill, edgecolor=BLK, lw=lw, zorder=z))


def arrow(ax, a, b, color=BLK, lw=1.0, ms=8, z=6, rad=0.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=color, zorder=z,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=6.5, shrinkB=6.5))


def tag(ax, s):
    ax.text(0.0, 1.0, s, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left")


fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(12.6, 3.9))
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.04, top=0.88, wspace=0.05)
for ax in (axA, axB, axC):
    ax.set_aspect("equal"); ax.axis("off")

# ============================================================ (a) structure
rg = 0.66
cA = rg * np.array([np.cos(np.deg2rad(148)), np.sin(np.deg2rad(148))])
cB = rg * np.array([np.cos(np.deg2rad(32)), np.sin(np.deg2rad(32))])
cC = rg * np.array([0.0, -1.0])
gA = group(axA, cA, rg, 3, np.deg2rad(148 - 180))
gB = group(axA, cB, rg, 3, np.deg2rad(32 - 180))
gC = group(axA, cC, rg, 4, np.deg2rad(90), dashed=True)
for g in (gA, gB, gC):
    for p in g[1:]:
        node(axA, p)
node(axA, np.zeros(2), fill=FILL_I, r=0.125, lw=1.3)
axA.text(0.0, -0.30, r"$u$", ha="center", va="top", fontsize=10.5,
         bbox=dict(fc="white", ec="none", pad=1.2), zorder=8)

axA.text(-1.40, 1.28, "layer 1", fontsize=9.5, ha="right")
axA.text(-1.40, 1.02, r"$m_1=3$", color=GRY, fontsize=9, ha="right")
axA.text(-1.40, -1.10, "layer 2", fontsize=9.5, ha="right")
axA.text(-1.40, -1.36, r"$m_2=4$", color=GRY, fontsize=9, ha="right")
axA.text(1.45, 0.46, r"$k^{(1)}(u)=2$", fontsize=9.5, ha="left")
axA.text(1.45, 0.16, r"$k^{(2)}(u)=1$", fontsize=9.5, ha="left")
axA.text(1.45, -0.24, "one state per node,", fontsize=8.6, color=GRY, ha="left")
axA.text(1.45, -0.48, "shared by all layers", fontsize=8.6, color=GRY, ha="left")
axA.set_xlim(-2.15, 2.75); axA.set_ylim(-1.75, 1.55)
axA.set_title("multiplex hypergraph", fontsize=10.5, pad=4)
tag(axA, "a")

# ============================================================ (b) cascade
rb, yb = 0.50, 0.32


def snap(ax, cx, n_inf):
    pts = group(ax, (cx, yb), rb, 4, np.deg2rad(90))
    for j, p in enumerate(pts):
        node(ax, p, fill=FILL_I if j < n_inf else FILL_S)
    return pts


p1, p2, p3 = snap(axB, -1.35, 1), snap(axB, 0.10, 2), snap(axB, 1.55, 3)
arrow(axB, p1[0], p1[1], color=RED, lw=1.1)
arrow(axB, p2[1], p2[2], color=RED, lw=1.1)
for cx, lab in ((-1.35, "one seed"), (0.10, "infects a member"),
                (1.55, "which infects another")):
    axB.text(cx, yb - 0.84, lab, ha="center", va="top", fontsize=8.6, color=GRY)

ya = -1.30
axB.plot([-1.35, 0.55], [ya, ya], "-", color=BLK, lw=1.6, solid_capstyle="butt")
axB.plot([0.55, 1.95], [ya, ya], "-", color=RED, lw=1.6, solid_capstyle="butt")
for x in (-1.35, 0.55, 1.95):
    axB.plot([x, x], [ya - 0.09, ya + 0.09], "-", color=BLK, lw=1.0)
axB.text(-0.40, ya + 0.15, "seed infectious", fontsize=8.2, color=GRY, ha="center")
axB.text(1.25, ya + 0.15, "clock extended", fontsize=8.2, color=RED, ha="center")
axB.text(-1.35, ya - 0.28, "the group transmits while any member is infectious",
         fontsize=8.4, color=GRY, ha="left", va="top")

axB.text(0.10, 1.58, r"$C(m,\lambda,\theta)\;>\;(m-1)\,T$", ha="center",
         fontsize=11.5)
axB.set_xlim(-2.00, 2.20); axB.set_ylim(-1.95, 1.80)
axB.set_title("intra-group cascade", fontsize=10.5, pad=4)
tag(axB, "b")

# ============================================================ (c) branching
v = np.zeros(2)
rc = 0.62
cIn = rc * np.array([np.cos(np.deg2rad(180)), np.sin(np.deg2rad(180))])
c1 = rc * np.array([np.cos(np.deg2rad(56)), np.sin(np.deg2rad(56))])
c2 = rc * np.array([np.cos(np.deg2rad(-56)), np.sin(np.deg2rad(-56))])
gIn = group(axC, cIn, rc, 3, np.deg2rad(0), color=GRY)
g1 = group(axC, c1, rc, 3, np.deg2rad(56 - 180))
g2 = group(axC, c2, rc, 4, np.deg2rad(-56 - 180), dashed=True)
for p in gIn[1:]:
    node(axC, p, fill=FILL_R, lw=0.8)
for p in g1[1:] + g2[1:]:
    node(axC, p)
arrow(axC, gIn[1], v, color=RED, lw=1.2)
node(axC, v, fill=FILL_I, r=0.125, lw=1.3)

xm = cIn.copy()                      # cross sits on the arriving group itself
axC.plot([xm[0] - 0.10, xm[0] + 0.10], [xm[1] - 0.10, xm[1] + 0.10], "-",
         color=GRY, lw=1.5, zorder=7)
axC.plot([xm[0] - 0.10, xm[0] + 0.10], [xm[1] + 0.10, xm[1] - 0.10], "-",
         color=GRY, lw=1.5, zorder=7)
axC.text(-1.40, -1.34, "arriving group,", fontsize=8.5, color=GRY,
         ha="center", va="center")
axC.text(-1.40, -1.58, r"not counted again $(-\delta_{ab})$", fontsize=8.5,
         color=GRY, ha="center", va="center")
for cc in (c1, c2):
    d = cc / np.linalg.norm(cc)
    axC.text(*(cc + 0.98 * d), r"$C$", fontsize=11, ha="center", va="center",
             color=RED)
axC.text(0.0, 1.80,
         r"$N_{ab}=C(m_b,\lambda_b,\theta_b)\!\left[\dfrac{\langle k^{(a)}k^{(b)}\rangle}"
         r"{\langle k^{(a)}\rangle}-\delta_{ab}\right]$", ha="center",
         va="center", fontsize=10.5)
axC.set_xlim(-2.15, 2.15); axC.set_ylim(-2.05, 2.15)
axC.set_title("group-level branching", fontsize=10.5, pad=4)
tag(axC, "c")

fig.savefig("figure1_mechanism.pdf", bbox_inches="tight")
fig.savefig("figure1_mechanism.png", bbox_inches="tight", dpi=210)
print("saved figure1_mechanism.pdf / .png")
