"""Fig. 1, PowerPoint-reproducible variant.

Every element maps to one shape from PowerPoint's own gallery:
  hyperedge / node -> Oval          arrow -> Line (with arrowhead)
  K box            -> Rounded Rectangle    label -> Text Box
No freeform curves and no Bezier arrows, so the figure can be rebuilt by hand.
"""
import sys, os, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from style import *

W, H = 190.0, 145.0
OUT_MM = float(os.environ.get("FIG_WIDTH_MM", 190.0))
fig = plt.figure(figsize=(OUT_MM / 25.4, OUT_MM / W * H / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
rule(ax, 95, 14, 95, 141); rule(ax, 5, 76, 186, 76)

# ═════════════════════════════  (a)  three circles meeting at u
panel_head(ax, 5, 140, "a", "Multiplex hypergraph")
u = (36, 113)
R1, R2 = 10.0, 11.0

def through_u(th, R, m, layer):
    """Circle of radius R whose centre lies at angle th from u; u is a member."""
    c = on_ell(u, R, R, 0, th)
    ts = [(th + 180 + i * 360 / m) % 360 for i in range(m)]
    return hyperedge(ax, c, R, R, 0, ts, layer)

e1 = through_u(155, R1, 3, L1)
e2 = through_u(25,  R1, 3, L1)
e3 = through_u(270, R2, 4, L2)
for p in e1[1:] + e2[1:] + e3[1:]: node(ax, p, "S")
node(ax, u, "U", r=2.0)
txt(ax, 39, 110.2, "$u$", size=FS_MATH, style="italic")
txt(ax, 21, 123, "$e_1$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 51, 123, "$e_2$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 36, 97, "$e_3$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 66, 121, "$k^{(a)}(u)=2$", size=FS_MATH)
txt(ax, 66, 115, "$k^{(b)}(u)=1$", size=FS_MATH)
txt(ax, 66, 106, "$P(\\mathbf{k})\\rightarrow\\Psi,\\,\\psi_a$", size=FS_MATH)

# ═════════════════════════════  (b)  one circle + straight arrows
panel_head(ax, 99, 140, "b", "Edge-based closure")
ce = (126, 118)
ut, nI, nS, nR = hyperedge(ax, ce, 11, 11, 0, [180, 90, 0, 270], L1)
node(ax, nI, "I"); node(ax, nS, "S"); node(ax, nR, "R"); node(ax, ut, "U", r=2.0)
txt(ax, 116.5, 114.5, "$u$", size=FS_MATH, style="italic")
arrow(ax, nI, ut, color=C_TRAN, sb=2.7)
txt(ax, 122.5, 121.5, "$\\beta_a$", size=FS_SM, color=C_TRAN)
arrow(ax, ut, (104, 112), color=C_TRAN, sa=2.7, sb=0.3, alpha=0.5)
cross(ax, (109, 115), s=1.5)
txt(ax, 130, 112, "$x^{(a)}_{sir}$", size=FS_MATH, ha="center")
cc = (168, 124)
cl = hyperedge(ax, cc, 9, 9, 0, [45, 135, 225, 315],
               dict(ec="#B8BCC0", fc="#B8BCC0", alpha_f=0.18, ls="-"), lw=0.7)
for p in cl: node(ax, p, "S", r=1.1, lw=0.55, alpha=0.5)
txt(ax, 168, 112.5, "rest of network", size=FS_SM, ha="center", color="#6A6A6A")
arrow(ax, (159.5, 121), nS, color=C_COUP, sa=1.5, sb=2.5)
txt(ax, 150, 124, "$h_a$", size=FS_MATH, ha="center", color=C_COUP)

cyc, lab = [105, 124, 141, 159, 178], ["$x^{(a)}_{sir}$", "$A_a$", "$\\Phi_a$",
                                       "$\\psi_a(\\Phi)$", "$h_a$"]
YC = 97.0
for x, s in zip(cyc, lab): txt(ax, x, YC, s, size=FS_MATH, ha="center")
for i in range(4):
    arrow(ax, (cyc[i] + 5.8, YC), (cyc[i + 1] - 5.8, YC),
          color="#4A4A4A", lw=0.8, sa=0, sb=0, ms=5)
# feedback drawn as three straight segments -- a hand-drawable elbow
FB, GY = 88.0, "#4A4A4A"
rule(ax, 178, YC - 3.2, 178, FB, GY, 0.8, z=4)
rule(ax, 178, FB, 105, FB, GY, 0.8, z=4)
arrow(ax, (105, FB), (105, YC - 3.2), color=GY, lw=0.8, sa=0, sb=0, ms=5)
txt(ax, 141, 85.5, "closure", size=FS_SM, ha="center", color="#7A7A7A")
txt(ax, 99, 79.5, "$\\Rightarrow\;S(t),\\,I(t),\\,R(t),\;R(\\infty)$", size=FS_MATH)

# ═════════════════════════════  (c)  circles, straight arrows
panel_head(ax, 5, 73, "c", "Within-hyperedge factor $C$")
txt(ax, 5, 58, "direct", size=FS_LBL)
seed, s1, s2, s3 = hyperedge(ax, (38, 58), 9, 9, 0, [135, 45, 315, 225], L1)
for p in (s1, s2, s3): node(ax, p, "S")
node(ax, seed, "I")
for p in (s1, s2, s3): arrow(ax, seed, p, color=C_TRAN, lw=0.8)
txt(ax, 56, 58, "$\\Rightarrow\;(m-1)\\,T$", size=FS_MATH)

txt(ax, 5, 34, "secondary", size=FS_LBL)
for cx, states, arr in [(30, "ISSS", (0, 1)), (52, "IISS", (1, 3)), (74, "IISI", None)]:
    pos = hyperedge(ax, (cx, 34), 7, 7, 0, [135, 45, 225, 315], L1)
    for p, st in zip(pos, states): node(ax, p, st, r=1.45)
    if arr: arrow(ax, pos[arr[0]], pos[arr[1]], color=C_TRAN, lw=0.8, sa=1.9, sb=2.4)
for cx, s in zip([30, 52, 74], ["$(i,s)=(1,3)$", "$(2,2)$", "$(3,1)$"]):
    txt(ax, cx, 22.0, s, size=FS_SM, ha="center")
for x in (41, 63):
    arrow(ax, (x - 2, 34), (x + 2, 34), color="#4A4A4A", lw=0.8, sa=0, sb=0, ms=4.5)
txt(ax, 84, 34, "$\\Rightarrow\;C$", size=FS_MATH)
txt(ax, 46, 16.0, "$C\;\\geq\;(m-1)\\,T$", size=FS_MATH, ha="center")

# ═════════════════════════════  (d)  three ellipses meeting at F
panel_head(ax, 99, 73, "d", "Excess hyperdegree $B$")
F = (134, 44)
ar = hyperedge(ax, (122, 44), 12, 6.5, 0, [0, 120, 240], L1, alpha_scale=0.42)
node(ax, ar[1], "I", alpha=0.40); node(ax, ar[2], "R", alpha=0.40)
cross(ax, (120, 44), s=2.1, color="#606060", lw=1.4)
txt(ax, 110, 32, "$-\\,\\delta_{ab}$", size=FS_MATH, color="#6A6A6A")
ta = hyperedge(ax, on_ell(F, 13, 13, 0, 44), 13, 5.5, 44, [180, 60, 300], L1)
tb = hyperedge(ax, on_ell(F, 14, 14, 0, -46), 14, 6.0, -46, [180, 55, 0, 305], L2)
for p in ta[1:] + tb[1:]: node(ax, p, "S")
node(ax, F, "I", r=2.0)
txt(ax, 163, 66, "layer $a$", size=FS_SM); txt(ax, 163, 61.5, "$B_{aa}$", size=FS_MATH)
txt(ax, 165, 30, "layer $b$", size=FS_SM); txt(ax, 165, 25.5, "$B_{ab}$", size=FS_MATH)
roundbox(ax, 99, 15.0, 44, 7.2)
txt(ax, 102, 18.6, "$K_{ab}=B_{ab}\\,C_b$", size=FS_MATH + 0.8, weight="bold")

# ═════════════════════════════  legend
rule(ax, 5, 11.0, 186, 11.0)
for x, st, s in [(6, "S", "susceptible"), (30, "I", "infected"),
                 (52, "R", "recovered"), (75, "U", "test node $u$")]:
    node(ax, (x, 6.0), st, r=1.55); txt(ax, x + 3.0, 6.0, s, size=FS_SM)
rule(ax, 100, 6.0, 106, 6.0, C_INK, 0.9); txt(ax, 108, 6.0, "layer $a$", size=FS_SM)
rule(ax, 126, 6.0, 132, 6.0, C_INK, 0.9, dashed=True)
txt(ax, 134, 6.0, "layer $b$", size=FS_SM)

if os.environ.get("EXPORT_SPEC"):
    export_spec(os.environ["EXPORT_SPEC"], fig, ax, W, H)
tag = "%dmm" % round(OUT_MM)
fig.savefig("fig1_ppt_%s.pdf" % tag); fig.savefig("fig1_ppt_%s.png" % tag, dpi=400)
print("wrote fig1_ppt_%s.pdf  (%.1f x %.1f mm)" % (tag, OUT_MM, OUT_MM / W * H))
