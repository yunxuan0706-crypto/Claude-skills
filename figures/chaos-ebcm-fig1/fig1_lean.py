"""Fig. 1, lean variant: the figure keeps only labels that name something drawn
in it. Every sentence and every equation number lives in the caption instead."""
import sys, pathlib, os
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from style import *

W, H = 190.0, 138.0
OUT_MM = float(os.environ.get("FIG_WIDTH_MM", 190.0))
fig = plt.figure(figsize=(OUT_MM / 25.4, OUT_MM / W * H / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
rule(ax, 95, 14, 95, 134)
rule(ax, 5, 72, 186, 72)

# ───────────────────────────────  (a)
panel_head(ax, 5, 133, "a", "Multiplex hypergraph")
u = (36, 108)
e1 = [u, (22, 121), (14, 107)]
e2 = [u, (50, 121), (58, 107)]
e3 = [u, (22, 94), (50, 94), (36, 85)]
blob(ax, e3, L2); blob(ax, e1, L1); blob(ax, e2, L1)
for p in e1[1:] + e2[1:] + e3[1:]: node(ax, p, "S")
node(ax, u, "U", r=2.0)
txt(ax, 41.5, 105.6, "$u$", size=FS_MATH, style="italic")
txt(ax, 15.5, 115, "$e_1$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 56.5, 115, "$e_2$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 36, 89.5, "$e_3$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 68, 116, "$k^{(a)}(u)=2$", size=FS_MATH)
txt(ax, 68, 110, "$k^{(b)}(u)=1$", size=FS_MATH)
txt(ax, 68, 101, "$P(\\mathbf{k})\\rightarrow\\Psi,\\,\\psi_a$", size=FS_MATH)

# ───────────────────────────────  (b)
panel_head(ax, 99, 133, "b", "Edge-based closure")
ut, n1, n2, n3 = (112, 114), (127, 124), (139, 113), (126, 103)
blob(ax, [ut, n1, n2, n3], L1)
node(ax, n1, "I"); node(ax, n2, "S"); node(ax, n3, "R")
node(ax, ut, "U", r=2.0)
txt(ax, 115.5, 111.4, "$u$", size=FS_MATH, style="italic")
arrow(ax, n1, ut, color=C_TRAN, sb=2.7)
txt(ax, 124.0, 117.5, "$\\beta_a$", size=FS_SM, color=C_TRAN)
arrow(ax, ut, (101, 108), color=C_TRAN, sa=2.7, sb=0.3, alpha=0.5)
cross(ax, (105.5, 110.6), s=1.5)
txt(ax, 131, 107.5, "$x^{(a)}_{sir}$", size=FS_MATH, ha="center")
cl = [(160, 120), (170, 123), (175, 115), (164, 111)]
blob(ax, cl, dict(ec="#B8BCC0", fc="#B8BCC0", alpha_f=0.18, ls="-"), pad=3.4, lw=0.7)
for p in cl: node(ax, p, "S", r=1.1, lw=0.55, alpha=0.5)
txt(ax, 167, 105.5, "rest of network", size=FS_SM, ha="center", color="#6A6A6A")
arrow(ax, (160, 118), n2, color=C_COUP, sa=2.2, sb=2.5, rad=0.20)
txt(ax, 151, 119.5, "$h_a$", size=FS_MATH, ha="center", color=C_COUP)

cyc, lab = [105, 124, 141, 159, 178], ["$x^{(a)}_{sir}$", "$A_a$", "$\\Phi_a$",
                                       "$\\psi_a(\\Phi)$", "$h_a$"]
YC = 89.0
for x, s in zip(cyc, lab): txt(ax, x, YC, s, size=FS_MATH, ha="center")
for i in range(4):
    arrow(ax, (cyc[i] + 5.8, YC), (cyc[i + 1] - 5.8, YC),
          color="#4A4A4A", lw=0.8, sa=0, sb=0, ms=5)
curve(ax, (178, YC - 2.6), (105, YC - 2.6), -0.22, "#4A4A4A", 0.8)
txt(ax, 99, 75.5, "$\\Rightarrow\;S(t),\\,I(t),\\,R(t),\;R(\\infty)$", size=FS_MATH)

# ───────────────────────────────  (c)
panel_head(ax, 5, 69, "c", "Within-hyperedge factor $C$")
txt(ax, 5, 54, "direct", size=FS_LBL)
dc, dr_ = (38, 54), 6.5
seed = (dc[0], dc[1] + dr_)
tg = [(dc[0] - dr_, dc[1]), (dc[0] + dr_, dc[1]), (dc[0], dc[1] - dr_)]
blob(ax, [seed] + tg, L1, pad=3.8)
for p in tg: node(ax, p, "S")
node(ax, seed, "I")
for p in tg: arrow(ax, seed, p, color=C_TRAN, lw=0.8)
txt(ax, 55, 54, "$\\Rightarrow\;(m-1)\\,T$", size=FS_MATH)

txt(ax, 5, 31, "secondary", size=FS_LBL)
for (cx, cy), states, arr in [((30, 31), "ISSS", (0, 3)), ((52, 31), "IISS", (1, 2)),
                              ((74, 31), "IIIS", None)]:
    r = 5.2
    pos = [(cx, cy + r), (cx - r, cy), (cx, cy - r), (cx + r, cy)]
    blob(ax, pos, L1, pad=3.0)
    for p, s in zip(pos, states): node(ax, p, s, r=1.45)
    if arr: arrow(ax, pos[arr[0]], pos[arr[1]], color=C_TRAN, lw=0.8, sa=1.9, sb=2.4)
for cx, s in zip([30, 52, 74], ["$(i,s)=(1,3)$", "$(2,2)$", "$(3,1)$"]):
    txt(ax, cx, 20.0, s, size=FS_SM, ha="center")
for x in (41, 63):
    arrow(ax, (x - 1.8, 31), (x + 1.8, 31), color="#4A4A4A", lw=0.8, sa=0, sb=0, ms=4.5)
txt(ax, 84, 31, "$\\Rightarrow\;C$", size=FS_MATH)
txt(ax, 46, 15.0, "$C\;\\geq\;(m-1)\\,T$", size=FS_MATH, ha="center")

# ───────────────────────────────  (d)
panel_head(ax, 99, 69, "d", "Excess hyperdegree $B$")
F = (131, 44)
arr_e = [F, (112, 51), (110, 38)]
blob(ax, arr_e, L1, alpha_scale=0.42)
node(ax, arr_e[1], "I", alpha=0.40); node(ax, arr_e[2], "R", alpha=0.40)
cross(ax, (116, 43), s=2.1, color="#606060", lw=1.4)
txt(ax, 112, 30.5, "$-\\,\\delta_{ab}$", size=FS_MATH, color="#6A6A6A")
ta = [F, (148, 59), (161, 54)]
tb = [F, (147, 29), (159, 25), (164, 34)]
blob(ax, tb, L2); blob(ax, ta, L1)
for p in ta[1:] + tb[1:]: node(ax, p, "S")
node(ax, F, "I", r=2.0)
txt(ax, 168, 60, "layer $a$", size=FS_SM); txt(ax, 168, 55.5, "$B_{aa}$", size=FS_MATH)
txt(ax, 172, 33, "layer $b$", size=FS_SM); txt(ax, 172, 28.5, "$B_{ab}$", size=FS_MATH)
roundbox(ax, 99, 14.6, 44, 7.2)
txt(ax, 102, 18.2, "$K_{ab}=B_{ab}\\,C_b$", size=FS_MATH + 0.8, weight="bold")

# ───────────────────────────────  legend
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
fig.savefig("fig1_lean_%s.pdf" % tag); fig.savefig("fig1_lean_%s.png" % tag, dpi=400)
print("wrote fig1_lean_%s.pdf  (%.1f x %.1f mm)" % (tag, OUT_MM, OUT_MM / W * H))
