import sys, numpy as np
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from style import *
from matplotlib.patches import FancyBboxPatch

W, H = 190.0, 158.0                       # layout canvas, in mm
# Physical output width. 190 mm = Elsevier/CSF double column;
# 177.8 mm = REVTeX 4 two-column \textwidth (7 in) for AIP Chaos.
OUT_MM = float(__import__("os").environ.get("FIG_WIDTH_MM", 190.0))
fig = plt.figure(figsize=(OUT_MM / 25.4, OUT_MM / W * H / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
ax.plot([95, 95], [15, 154], color="#DCDCDC", lw=0.6, zorder=0)
ax.plot([5, 186], [81, 81], color="#DCDCDC", lw=0.6, zorder=0)

# ═══════════════════════════════════════════════  (a) structure
panel_head(ax, 5, 155, "a", "Multiplex hypergraph and hyperdegrees")
u = (36, 124)
e1 = [u, (22, 137), (14, 123)]
e2 = [u, (50, 137), (58, 123)]
e3 = [u, (22, 110), (50, 110), (36, 101)]
blob(ax, e3, L2); blob(ax, e1, L1); blob(ax, e2, L1)
for p in e1[1:] + e2[1:] + e3[1:]: node(ax, p, "S")
node(ax, u, "U", r=2.0)
txt(ax, 41.5, 121.8, "$u$", size=FS_MATH, style="italic")
txt(ax, 15.5, 131, "$e_1$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 56.5, 131, "$e_2$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 36, 105.5, "$e_3$", size=FS_SM, ha="center", color="#5A5A5A")
txt(ax, 68, 133, "$k^{(a)}(u)=2$", size=FS_MATH)
txt(ax, 68, 127, "$k^{(b)}(u)=1$", size=FS_MATH)
txt(ax, 68, 119, "layer $a$:  $m_a=3$", size=FS_SM, color="#5A5A5A")
txt(ax, 68, 114, "layer $b$:  $m_b=4$", size=FS_SM, color="#5A5A5A")
txt(ax, 5, 94.5, "joint hyperdegree distribution  $P(\\mathbf{k})$,   $\\mathbf{k}=(k^{(1)},\\ldots,k^{(M)})$", size=FS_SM)
txt(ax, 5, 90.0, "$\\Psi(\\mathbf{x})$ — uniformly sampled node                                   [Eq. (3)]", size=FS_SM)
txt(ax, 5, 85.5, "$\\psi_a(\\mathbf{x})$ — node reached via layer $a$, arrival edge excluded    [Eq. (4)]", size=FS_SM)

# ═══════════════════════════════════════  (b) edge-based closure, finite prevalence
panel_head(ax, 99, 155, "b", "Edge-based closure at finite prevalence")
ut = (127, 136)
n1, n2, n3 = (139, 144), (149, 135), (138, 127)
blob(ax, [ut, n1, n2, n3], L1)
node(ax, n1, "I"); node(ax, n2, "S"); node(ax, n3, "R")
node(ax, ut, "U", r=2.0)
txt(ax, 130.8, 133.6, "$u$", size=FS_MATH, style="italic")
txt(ax, 99, 142.0, "layer-$a$ hyperedge $e$", size=FS_SM, color="#5A5A5A")
arrow(ax, n1, ut, color=C_TRAN, sb=2.7)
txt(ax, 130.5, 141.5, "$\\beta_a$", size=FS_SM, color=C_TRAN)
arrow(ax, ut, (113, 136), color=C_TRAN, sa=2.7, sb=0.3, alpha=0.5)
cross(ax, (118.5, 136.0), s=1.5)
txt(ax, 99, 130.5, "$u$ cannot transmit", size=FS_SM, color=C_TRAN)
cl = [(165, 143), (175, 146), (180, 138), (169, 134)]
blob(ax, cl, dict(ec="#B8BCC0", fc="#B8BCC0", alpha_f=0.18, ls="-"), pad=3.4, lw=0.7)
for p in cl: node(ax, p, "S", r=1.1, lw=0.55, alpha=0.5)
txt(ax, 172, 129.0, "rest of network", size=FS_SM, ha="center", color="#6A6A6A")
arrow(ax, (165, 140), n2, color=C_COUP, sa=2.2, sb=2.5, rad=0.20)
txt(ax, 158, 145.0, "$h_a$", size=FS_MATH, ha="center", color=C_COUP)
txt(ax, 99, 116.0, "$x^{(a)}_{sir}$ — joint SIR composition of the $m_a\\!-\\!1$ other members  [Eq. (9)]", size=FS_SM)
txt(ax, 99, 111.5, "$\\Phi_a$ — $e$ has not infected $u$;   $A_a$ — $e$ is transmitting  [Eqs. (10)–(12)]", size=FS_SM)

cyc = [105, 124, 141, 159, 178]
lab = ["$x^{(a)}_{sir}$", "$A_a$", "$\\Phi_a$", "$\\psi_a(\\Phi)$", "$h_a$"]
eqs = ["(11)", "(12)", "(4)", "(13)"]
YC = 103.0
for x, s in zip(cyc, lab): txt(ax, x, YC, s, size=FS_MATH, ha="center")
for i in range(4):
    x0, x1 = cyc[i] + 5.8, cyc[i + 1] - 5.8
    arrow(ax, (x0, YC), (x1, YC), color="#4A4A4A", lw=0.8, sa=0, sb=0, ms=5)
    if eqs[i]: txt(ax, (x0 + x1) / 2, YC + 3.4, eqs[i], size=FS_SM, ha="center", color="#7A7A7A")
ax.add_patch(FancyArrowPatch((178, YC - 2.6), (105, YC - 2.6), arrowstyle="-|>",
                             mutation_scale=5, lw=0.8, color="#4A4A4A", zorder=4,
                             connectionstyle="arc3,rad=-0.22"))
txt(ax, 141, 91.6, "closure   [Eq. (14)]", size=FS_SM, ha="center", color="#7A7A7A",
    bbox=dict(facecolor="white", edgecolor="none", pad=1.4))
txt(ax, 99, 85.5, "$\\Rightarrow$   $S(t),\\,I(t),\\,R(t)$   and   $R(\\infty)$   [Eqs. (8), (17), (18)]", size=FS_MATH)

# ═══════════════════════════════  (c) rare-infection limit I: within-hyperedge C
panel_head(ax, 5, 78, "c", "Rare-infection limit I: within-hyperedge factor $C_b$")
txt(ax, 5, 63.5, "direct\nonly", size=FS_LBL, va="center", linespacing=1.35)
dc, dr_ = (40, 63.5), 7.2
seed = (dc[0], dc[1] + dr_)
tg = [(dc[0] - dr_, dc[1]), (dc[0] + dr_, dc[1]), (dc[0], dc[1] - dr_)]
blob(ax, [seed] + tg, L1, pad=3.8)
for p in tg: node(ax, p, "S")
node(ax, seed, "I")
for p in tg: arrow(ax, seed, p, color=C_TRAN, lw=0.8)
txt(ax, 58, 63.5, "$\\Rightarrow\;(m-1)\\,T$", size=FS_MATH)

txt(ax, 5, 45, "with\nsecondary", size=FS_LBL, va="center", linespacing=1.35)
for (cx, cy), states, arr in [((30, 45), "ISSS", (0, 3)), ((52, 45), "IISS", (1, 2)),
                              ((74, 45), "IIIS", None)]:
    r = 5.2
    pos = [(cx, cy + r), (cx - r, cy), (cx, cy - r), (cx + r, cy)]
    blob(ax, pos, L1, pad=3.0)
    for p, s in zip(pos, states): node(ax, p, s, r=1.45)
    if arr: arrow(ax, pos[arr[0]], pos[arr[1]], color=C_TRAN, lw=0.8, sa=1.9, sb=2.4)
for cx, s in zip([30, 52, 74], ["$(i,s)=(1,3)$", "$(i,s)=(2,2)$", "$(i,s)=(3,1)$"]):
    txt(ax, cx, 35.5, s, size=FS_SM, ha="center")
for x in (41, 63):
    arrow(ax, (x - 1.8, 45), (x + 1.8, 45), color="#4A4A4A", lw=0.8, sa=0, sb=0, ms=4.5)
txt(ax, 84, 45, "$\\Rightarrow\;C$", size=FS_MATH)

bx = 26.0
for y, segs in [(29.5, [(0, 20, "#333333")]), (25.0, [(0, 20, "#333333"), (20, 36, C_TRAN)])]:
    for a, b, c in segs:
        ax.plot([bx + a, bx + b], [y, y], color=c, lw=1.9, solid_capstyle="butt", zorder=4)
        ax.plot([bx + b, bx + b], [y - 1.1, y + 1.1], color=c, lw=0.7, zorder=4)
ax.plot([bx, bx], [23.9, 30.6], color="#777777", lw=0.7, zorder=3)
txt(ax, 5, 29.5, "direct", size=FS_SM)
txt(ax, 5, 25.0, "secondary", size=FS_SM)
txt(ax, 64, 25.0, "extension", size=FS_SM, color=C_TRAN)
txt(ax, 26, 20.5, "time with at least one infectious node in $e$", size=FS_SM, color="#6A6A6A")
txt(ax, 46, 15.5, "$C(m,\\lambda,\\theta)\;\\geq\;(m-1)\\,T$      (equality iff $m=2$)",
    size=FS_MATH, ha="center")

# ═══════════════════  (d) rare-infection limit II: excess hyperdegree, K = B C
panel_head(ax, 99, 78, "d", "Rare-infection limit II: excess hyperdegree $B_{ab}$")
F = (131, 57)
arr_e = [F, (112, 64), (110, 51)]
blob(ax, arr_e, L1, alpha_scale=0.42)
node(ax, arr_e[1], "I", alpha=0.40); node(ax, arr_e[2], "R", alpha=0.40)
cross(ax, (116.0, 56.0), s=2.1, color="#606060", lw=1.4)
txt(ax, 99, 43.5, "arrival hyperedge", size=FS_SM, color="#6A6A6A")
txt(ax, 99, 39.0, "excluded  ($-\\delta_{ab}$)", size=FS_SM, color="#6A6A6A")
ta = [F, (147, 69), (160, 65)]
tb = [F, (145, 40), (156, 36), (163, 43)]
blob(ax, tb, L2); blob(ax, ta, L1)
for p in ta[1:] + tb[1:]: node(ax, p, "S")
node(ax, F, "I", r=2.0)
ax.annotate("newly infected", xy=F, xytext=(126, 72.5), fontsize=FS_SM, color="#5A5A5A",
            zorder=9, va="bottom", ha="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#9A9A9A", shrinkA=2, shrinkB=6))
txt(ax, 167, 67.5, "layer $a$", size=FS_SM)
txt(ax, 171, 41.5, "layer $b$", size=FS_SM)
txt(ax, 99, 29.0, "a type-$a$ node reaches $B_{ab}$ remaining layer-$b$ hyperedges", size=FS_SM)
txt(ax, 99, 24.0, "$B_{ab}=\\langle k^{(a)}k^{(b)}\\rangle/\\langle k^{(a)}\\rangle-\\delta_{ab}$   [Eq. (23)]", size=FS_SM)
ax.add_patch(FancyBboxPatch((99, 13.4), 44, 7.2, boxstyle="round,pad=0,rounding_size=1.2",
                            fc="#F2F4F6", ec="#B9C0C7", lw=0.7, zorder=3))
txt(ax, 102, 17.0, "$K_{ab}=B_{ab}\\,C_b$", size=FS_MATH + 0.8, weight="bold")
txt(ax, 147, 17.0, "$\\rho(K)=1\;\\Rightarrow\;\\lambda_c$   [Eqs. (24), (25)]", size=FS_SM)

# ═══════════════════════════════════════════════  legend
ax.plot([5, 186], [11.6, 11.6], color="#DCDCDC", lw=0.6, zorder=0)
for x, st, s in [(6, "S", "susceptible"), (27, "I", "infected"),
                 (45, "R", "recovered"), (64, "U", "test node $u$")]:
    node(ax, (x, 6.2), st, r=1.55); txt(ax, x + 3.0, 6.2, s, size=FS_SM)
ax.plot([86, 92], [6.2, 6.2], color=C_INK, lw=0.9); txt(ax, 94, 6.2, "layer $a$", size=FS_SM)
ax.plot([108, 114], [6.2, 6.2], color=C_INK, lw=0.9, ls=(0, (3, 2)))
txt(ax, 116, 6.2, "layer $b$", size=FS_SM)
arrow(ax, (130, 6.2), (136, 6.2), color=C_TRAN, lw=0.9, sa=0, sb=0, ms=5)
txt(ax, 138, 6.2, "transmission", size=FS_SM)
arrow(ax, (158, 6.2), (164, 6.2), color=C_COUP, lw=0.9, sa=0, sb=0, ms=5)
txt(ax, 166, 6.2, "coupling", size=FS_SM)

tag = "%dmm" % round(OUT_MM)
fig.savefig("fig1_csf_%s.pdf" % tag); fig.savefig("fig1_csf_%s.png" % tag, dpi=400)
print("wrote fig1_csf_%s.pdf  (%.1f x %.1f mm)" % (tag, OUT_MM, OUT_MM / W * H))
