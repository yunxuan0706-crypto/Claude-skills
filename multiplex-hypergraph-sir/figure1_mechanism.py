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
  (c) the group-level branching step and the excess subtraction. A node reached
      via a layer-a group re-ignites its other groups; with X_ab = <k^a k^b>/<k^a>
      it belongs to X_ab groups of layer b, but the arrival group -- itself a
      layer-a group -- is already spent, so N_ab = C(m_b)[X_ab - delta_ab]. The
      subtraction is diagonal: a layer-a successor carries X_aa-1, a layer-b
      successor the full X_ab.
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
    # same weight/size as the panel tags of figures 2-5, so the five figures
    # carry one labelling convention
    ax.text(0.0, 1.0, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")
    ax.text(0.048, 0.997, title, transform=ax.transAxes, fontsize=10,
            va="top", ha="left")


fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(12.6, 4.35))
fig.subplots_adjust(left=0.012, right=0.988, bottom=0.03, top=0.92, wspace=0.10)
for ax in (axA, axB, axC):
    ax.set_aspect("equal"); ax.axis("off")
for x in (0.3455, 0.6725):                      # thin rules between panels
    fig.add_artist(Line2D([x, x], [0.05, 0.93], color=LGY, lw=0.7))

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
axA.set_xlim(-1.55, 2.55); axA.set_ylim(-2.50, 1.95)
panel_label(axA, "a", "multiplex hypergraph")

# ================================================================ (b)
rb = 0.42
y_nv, y_cs = 1.16, -0.44                     # naive row / cascade row


def grp(ax, cx, cy, n_inf, r=rb):
    pts = group(ax, (cx, cy), r, 4, np.deg2rad(90))
    for j, p in enumerate(pts):
        node(ax, p, fill=F_I if j < n_inf else F_S)
    return pts

# --- row 1: the naive count -- only the seed transmits -------------------
q = grp(axB, -1.05, y_nv, 1)
for j in (1, 2, 3):
    arrow(axB, q[0], q[j], color=RED, lw=0.85, ms=6.5)
axB.text(-1.68, y_nv, "naive", fontsize=8.8, color=GRY, ha="right", va="center")
axB.text(0.02, y_nv, r"$\Rightarrow\;(m-1)\,T$", fontsize=10, va="center",
         ha="left")

# --- row 2: the cascade --------------------------------------------------
p1 = grp(axB, -1.05, y_cs, 1)
p2 = grp(axB, 0.12, y_cs, 2)
p3 = grp(axB, 1.29, y_cs, 3)
arrow(axB, p1[0], p1[1], color=RED, lw=0.9, ms=7)
arrow(axB, p2[1], p2[2], color=RED, lw=0.9, ms=7)
axB.text(-1.68, y_cs, "cascade", fontsize=8.8, color=BLK, ha="right", va="center")
for cx, tri in ((-1.05, "(1,3)"), (0.12, "(2,2)"), (1.29, "(3,1)")):
    axB.text(cx, y_cs - 0.62, rf"$(i,s)={tri}$", ha="center", va="top",
             fontsize=8.4)
axB.text(1.96, y_cs, r"$\Rightarrow\;C$", fontsize=10, va="center", ha="left")

# --- the two active periods, aligned under the rows ----------------------
ya = -1.72
axB.plot([-1.05, 0.20], [ya + 0.24, ya + 0.24], "-", color=BLK, lw=1.2,
         solid_capstyle="butt")
axB.plot([-1.05, 0.20], [ya, ya], "-", color=BLK, lw=1.2, solid_capstyle="butt")
axB.plot([0.20, 1.55], [ya, ya], "-", color=RED, lw=1.2, solid_capstyle="butt")
for x, y in ((-1.05, ya + 0.24), (0.20, ya + 0.24), (-1.05, ya), (0.20, ya),
             (1.55, ya)):
    axB.plot([x, x], [y - 0.065, y + 0.065], "-", color=BLK, lw=0.8)
axB.text(-1.18, ya + 0.24, "naive", fontsize=8.4, color=GRY, ha="right",
         va="center")
axB.text(-1.18, ya, "cascade", fontsize=8.4, color=BLK, ha="right", va="center")
axB.text(0.88, ya + 0.15, "extension", fontsize=8.4, color=RED, ha="center")
axB.text(0.25, ya - 0.30, "active period of the group", fontsize=8.4,
         color=GRY, ha="center", va="top")

axB.text(0.25, ya - 0.72, r"$C(m,\lambda,\theta)\;>\;(m-1)\,T$", ha="center",
         va="top", fontsize=10.5)
axB.set_xlim(-2.45, 2.55); axB.set_ylim(-2.95, 1.95)
panel_label(axB, "b", "intra-group cascade")

# ================================================================ (c)
# A node reached VIA a layer-a group (type a) re-ignites its other groups.
# Among its layer-a groups the arrival one is already spent, so it is removed
# -> the -delta subtraction. Because the arrival is a layer-a group, that
# removal only touches the same-layer (a=b) count: a layer-a successor carries
# X_aa-1, a layer-b successor carries the full X_ab.
rc = 0.52
u = np.zeros(2)
cIn = np.array([-1.24, -0.02])                       # arrival group (layer a)
ca = np.array([0.74, 0.86])                          # layer-a successor
cb = np.array([0.98, -0.70])                         # layer-b successor

gIn = group(axC, cIn, rc, 3, np.deg2rad(-18), color=GRY)   # solid = layer a
ga = group(axC, ca, rc, 3, np.deg2rad(205))               # solid = layer a
gb = group(axC, cb, rc, 4, np.deg2rad(150), dashed=True)  # dashed = layer b
for p in gIn[1:]:
    node(axC, p, fill=F_R, ec=GRY, lw=0.75)               # arrival: spent members
for p in ga[1:] + gb[1:]:
    node(axC, p)
arrow(axC, gIn[1], u, color=RED, lw=1.05)                 # arrival -> u
arrow(axC, u, ga[0], color=RED, lw=0.95, rad=-0.12)       # u re-ignites layer a
arrow(axC, u, gb[0], color=RED, lw=0.95, rad=0.12)        # u re-ignites layer b
node(axC, u, fill=F_I, r=0.112, lw=1.15)

# cross the arrival group out: it is the -delta term (removed from the count)
axC.plot([cIn[0] - 0.1, cIn[0] + 0.1], [cIn[1] - 0.1, cIn[1] + 0.1], "-",
         color=BLK, lw=1.3, zorder=8)
axC.plot([cIn[0] - 0.1, cIn[0] + 0.1], [cIn[1] + 0.1, cIn[1] - 0.1], "-",
         color=BLK, lw=1.3, zorder=8)

axC.text(cIn[0], cIn[1] + 0.80, r"layer $a$", fontsize=8.6, color=GRY, ha="center")
axC.text(cIn[0], cIn[1] - 0.86, r"$-\,\delta_{ab}$  ($a\!=\!b$ only)",
         fontsize=9, color=BLK, ha="center", va="center")
axC.text(ca[0] + 0.72, ca[1] + 0.10, r"layer $a$", fontsize=8.6, color=GRY, ha="left")
axC.text(ca[0] + 0.72, ca[1] - 0.24, r"$C\,(X_{aa}\!-\!1)$", fontsize=9,
         color=BLK, ha="left", va="center")
axC.text(cb[0] + 0.72, cb[1] + 0.10, r"layer $b$", fontsize=8.6, color=GRY, ha="left")
axC.text(cb[0] + 0.72, cb[1] - 0.24, r"$C\,X_{ab}$", fontsize=9,
         color=BLK, ha="left", va="center")
axC.text(0.30, -2.02, r"$N_{ab} = C\,(X_{ab} - \delta_{ab})$", ha="center",
         va="center", fontsize=10.5)
axC.set_xlim(-1.95, 2.55); axC.set_ylim(-2.50, 1.95)
panel_label(axC, "c", "group-level branching")

fig.savefig("figure1_mechanism.pdf", bbox_inches="tight")
fig.savefig("figure1_mechanism.png", bbox_inches="tight", dpi=210)
print("saved figure1_mechanism.pdf / .png")
