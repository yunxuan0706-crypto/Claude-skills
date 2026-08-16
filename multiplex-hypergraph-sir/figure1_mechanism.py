"""
Figure 1 (schematic): the mechanism behind the group closure.

Three things the quantitative figures presuppose but never draw:
  (a) what a multiplex hypergraph is, and what the layer degree counts;
  (b) the intra-group cascade -- the reason C(m,lambda,theta) exceeds (m-1)T,
      and the source of everything new in this model;
  (c) the group-level branching process, including the excess subtraction.
Panels (b) and (c) are exactly the mechanisms whose absence from the degree
mean field Figure 3 quantifies.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
L1, L2 = "#3E7CB1", "#5FA88D"          # layer 1 / layer 2 groups
INF, SUS, REC = "#E2703A", "#ffffff", "#d9d8d2"   # infectious / susceptible / recovered

mpl.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 340,
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
    "text.color": INK, "axes.linewidth": 0,
})


# ------------------------------------------------------------------ helpers
def blob(ax, pts, color, pad=30, alpha=0.16, z=1):
    """Draw a hyperedge as a rounded convex blob around its member points."""
    pts = np.asarray(pts, float)
    c = pts.mean(0)
    order = np.argsort(np.arctan2(pts[:, 1] - c[1], pts[:, 0] - c[0]))
    p = pts[order]
    p = np.vstack([p, p[:1]])
    ax.plot(p[:, 0], p[:, 1], "-", color=color, lw=pad, alpha=alpha,
            solid_capstyle="round", solid_joinstyle="round", zorder=z)
    ax.fill(p[:, 0], p[:, 1], color=color, alpha=alpha, lw=0, zorder=z)


def node(ax, xy, state=SUS, r=0.115, z=5, lw=1.1, ec=None):
    ax.add_patch(Circle(xy, r, facecolor=state, edgecolor=ec or SEC,
                        lw=lw, zorder=z))


def arrow(ax, a, b, color=INK, lw=1.2, style="-|>", ms=9, z=6, rad=0.0, ls="-"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=ms,
                                 lw=lw, color=color, zorder=z, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=7, shrinkB=7))


def tag(ax, s):
    ax.text(0.005, 1.0, s, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left", color=INK)


fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.2, 4.05))
fig.subplots_adjust(left=0.012, right=0.99, bottom=0.03, top=0.90, wspace=0.06)
for ax in (axA, axB, axC):
    ax.set_aspect("equal"); ax.axis("off")

# =============================================================== (a) structure
u = np.array([0.0, 0.0])
gA = [u, np.array([-1.30, 0.62]), np.array([-0.78, 1.28])]          # layer 1
gB = [u, np.array([1.30, 0.62]), np.array([0.78, 1.28])]            # layer 1
gC = [u, np.array([-0.95, -0.95]), np.array([0.0, -1.36]),
      np.array([0.95, -0.95])]                                      # layer 2
for g in (gA, gB):
    blob(axA, g, L1)
blob(axA, gC, L2)
for g in (gA, gB, gC):
    for p in g[1:]:
        node(axA, p)
node(axA, u, state="#f3f2ec", r=0.145, lw=1.7, ec=INK)
axA.text(0.0, 0.0, r"$u$", ha="center", va="center", fontsize=10, zorder=7)

axA.text(-1.62, 1.62, r"layer 1  ($m_1\!=\!3$)", color=L1, fontsize=9.5, ha="left")
axA.text(-1.62, -1.72, r"layer 2  ($m_2\!=\!4$)", color=L2, fontsize=9.5, ha="left")
axA.text(2.52, 0.34, r"$k^{(1)}(u)=2$", color=L1, fontsize=9.5, ha="right")
axA.text(2.52, -0.06, r"$k^{(2)}(u)=1$", color=L2, fontsize=9.5, ha="right")
axA.text(0.0, -2.06, "one state per node, shared by every layer", fontsize=8.8,
         color=SEC, ha="center")
axA.set_xlim(-1.95, 2.60); axA.set_ylim(-2.25, 1.95)
axA.set_title("multiplex hypergraph", fontsize=10.5, pad=6, color=INK)
tag(axA, "a")

# =============================================================== (b) cascade
def group_snapshot(ax, cx, cy, states, R=0.52):
    ang = np.linspace(90, 450, 5)[:4] * np.pi / 180
    pts = [np.array([cx + R * np.cos(t), cy + R * np.sin(t)]) for t in ang]
    blob(ax, pts, L1, pad=26, alpha=0.15)
    for p, s in zip(pts, states):
        node(ax, p, state=s, r=0.105)
    return pts

y0 = 0.30
s1 = group_snapshot(axB, -1.45, y0, [INF, SUS, SUS, SUS])
s2 = group_snapshot(axB, 0.10, y0, [INF, INF, SUS, SUS])
s3 = group_snapshot(axB, 1.65, y0, [INF, INF, INF, SUS])
arrow(axB, (-0.83, y0), (-0.52, y0), color=SEC, lw=1.1)
arrow(axB, (0.72, y0), (1.03, y0), color=SEC, lw=1.1)
for x, lab in ((-1.45, "one seed"), (0.10, "it infects a\nsecond member"),
               (1.65, "which extends the\ngroup's active period")):
    axB.text(x, y0 - 0.80, lab, ha="center", va="top", fontsize=8.3, color=SEC,
             linespacing=1.35)

# the active-clock bar
axB.plot([-1.45, 1.95], [-1.28, -1.28], "-", color=INF, lw=3.4, alpha=0.30,
         solid_capstyle="round")
axB.plot([-1.45, 0.62], [-1.28, -1.28], "-", color=INF, lw=3.4,
         solid_capstyle="round")
axB.text(-1.45, -1.50, "group active while any member is infectious",
         fontsize=8.3, color=SEC, ha="left", va="top")
axB.annotate("", xy=(1.95, -1.28), xytext=(0.62, -1.28),
             arrowprops=dict(arrowstyle="-|>", color=INF, lw=1.4, alpha=0.75))
axB.text(1.30, -1.10, "extended", fontsize=8.0, color=INF, ha="center")

axB.text(0.10, 1.98,
         r"$C(m,\lambda,\theta) \;>\; (m-1)\,T$",
         ha="center", fontsize=11.5, color=INK)
axB.text(0.10, 1.62,
         "each intra-group infection lengthens the clock, so members\n"
         "keep being exposed after the seed has recovered",
         ha="center", va="top", fontsize=8.3, color=SEC, linespacing=1.4)
axB.set_xlim(-2.15, 2.35); axB.set_ylim(-1.95, 2.40)
axB.set_title("intra-group cascade", fontsize=10.5, pad=6, color=INK)
tag(axB, "b")

# =============================================================== (c) branching
v = np.array([0.0, 0.25])
inc = [v, np.array([-1.35, 0.95]), np.array([-1.35, -0.30])]     # arriving group
blob(axC, inc, MUTED, pad=26, alpha=0.13)
for p in inc[1:]:
    node(axC, p, state=REC)
arrow(axC, (-1.30, 0.35), (-0.22, 0.28), color=INF, lw=1.5)
axC.text(-1.40, -0.62, "the group it was\ninfected through", fontsize=8.3,
         color=SEC, ha="center", va="top", linespacing=1.35)
axC.text(-0.72, 0.62, r"$-\,\delta_{ab}$", fontsize=10.5, color=INK, ha="center")
axC.text(-0.72, 0.30, "not counted\nagain", fontsize=7.8, color=MUTED,
         ha="center", va="top", linespacing=1.3)

out1 = [v, np.array([1.38, 0.86]), np.array([0.80, 1.24])]        # other layer-1
out2 = [v, np.array([1.45, -0.28]), np.array([1.05, -1.05]),
        np.array([0.20, -1.20])]                                  # layer-2
blob(axC, out1, L1); blob(axC, out2, L2)
for p in out1[1:] + out2[1:]:
    node(axC, p)
node(axC, v, state=INF, r=0.145, lw=1.5, ec=INK)

axC.text(0.05, 2.24, r"$C(m_b,\lambda_b,\theta_b)\;\times\;"
                     r"\left[\dfrac{\langle k^{(a)}k^{(b)}\rangle}"
                     r"{\langle k^{(a)}\rangle}-\delta_{ab}\right]$",
         fontsize=11, color=INK, ha="center", va="top")
axC.text(0.05, -1.74, "each remaining group is ignited and yields $C$ new\n"
         "infections: the branching process lives on groups, not nodes",
         fontsize=8.5, color=SEC, ha="center", va="top", linespacing=1.4)
axC.set_xlim(-2.15, 2.15); axC.set_ylim(-2.35, 2.40)
axC.set_title("group-level branching", fontsize=10.5, pad=6, color=INK)
tag(axC, "c")

fig.savefig("figure1_mechanism.pdf", bbox_inches="tight")
fig.savefig("figure1_mechanism.png", bbox_inches="tight", dpi=210)
print("saved figure1_mechanism.pdf / .png")
