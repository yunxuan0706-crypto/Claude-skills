"""Figure 2: the within-group cascade C(m,lambda,theta), quantified and checked.

C is the quantity that separates this model from a pairwise one: a group stays
active for the whole cascade, so a single seed infects more than the naive
(m-1)T that independent pairwise transmission would give. The figure is laid
out as a 2x3 grid -- (a-d) the quantity itself, (e-f) its verification against
a direct single-group Gillespie simulation that shares no code with the
recursion.

  (a) C(m,lambda) against the naive (m-1)T, m = 2..6;
  (b) the excess C/[(m-1)T] - 1 in percent, which is what the closure would
      lose if the cascade were dropped;
  (c) C over the (lambda, m) plane on a perceptually uniform map, stretched
      by a power norm so the low end opens up; the C=1 locus -- where one
      group alone already replaces its seed -- is the one heavy contour;
  (d) the recursion's internal state lattice u(i,s) for m=6, lambda=1, values
      printed on the cells and the one cell the closure reads off marked;

Both maps are colour-vision-safe and monotone in lightness, so the panels
reproduce in greyscale; overlays are cased rather than the maps altered.
  (e) recursion vs Gillespie on twelve (m,lambda,theta) settings;
  (f) the same comparison in units of the Monte-Carlo standard error.

MC results are cached in figure_cascade_data.json; delete it to re-measure.
"""
import json
import os
from functools import lru_cache

import numpy as np
import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator

from theory import cascade_C

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE, PLUM = "#1C9B8E", "#E76F51", "#2F5FD0", "#7D6B9E"
BAND2, BAND1 = "#eef2f6", "#dbe3ec"
CACHE = "figure_cascade_data.json"

# Overlays on a perceptually uniform map have to survive the whole range of
# that map, so every line and label drawn on one is cased -- the standard
# treatment, and the reason the colormaps themselves need no alteration.
def cased(lw, fg="#17171a"):
    return [pe.withStroke(linewidth=lw, foreground=fg)]

# the twelve settings the Monte-Carlo check runs on
MC_CASES = [(2, 0.35, 1), (2, 1.0, 1), (3, 0.5, 1), (3, 1.0, 1), (3, 2.0, 1),
            (4, 0.8, 1), (4, 1.0, 1), (4, 2.0, 1), (5, 1.0, 1), (5, 0.5, 1),
            (6, 1.0, 1), (6, 0.6, 1)]
NRUNS = 400_000


def u_lattice(m, lam, theta=1):
    """u(i,s) on the whole reachable lattice, same recursion as cascade_C."""
    @lru_cache(maxsize=None)
    def u(i, s):
        if i == 0 or s == 0:
            return 0.0
        a = lam * s if i >= theta else 0.0
        return (a / (a + i)) * (1.0 + u(i + 1, s - 1)) + (i / (a + i)) * u(i - 1, s)

    Z = np.full((m + 1, m + 1), np.nan)
    for i in range(0, m + 1):
        for s in range(0, m + 1):
            if 0 < i + s <= m:
                Z[i, s] = u(i, s)
    return Z


def cascade_C_mc(m, lam, theta, nruns, rng):
    """Vectorised Gillespie on one isolated m-group: seed one infected member,
    count how many further members are ever infected. No shared code with the
    recursion -- only the model definition is common."""
    i = np.ones(nruns, dtype=np.int64)
    s = np.full(nruns, m - 1, dtype=np.int64)
    alive = np.ones(nruns, dtype=bool)
    while alive.any():
        idx = np.flatnonzero(alive)
        ii, ss = i[idx], s[idx]
        r_inf = lam * ss * (ii >= theta)
        r_rec = ii.astype(float)
        p_inf = r_inf / (r_inf + r_rec)
        infect = rng.random(idx.size) < p_inf
        i[idx] = np.where(infect, ii + 1, ii - 1)
        s[idx] = np.where(infect, ss - 1, ss)
        alive[idx] = i[idx] > 0
    add = (m - s) - 1.0
    return add.mean(), add.std(ddof=1) / np.sqrt(nruns)


def mc_data():
    if os.path.exists(CACHE):
        return json.load(open(CACHE))
    rng = np.random.default_rng(20260828)
    th, mc, se = [], [], []
    for m, lam, t in MC_CASES:
        c_th = cascade_C(m, lam, t)
        c_mc, c_se = cascade_C_mc(m, lam, t, NRUNS, rng)
        th.append(c_th); mc.append(c_mc); se.append(c_se)
        print(f"  m={m} lam={lam} th={t}: rec={c_th:.5f} MC={c_mc:.5f}+-{c_se:.5f} "
              f"({abs(c_mc-c_th)/c_se:.2f}s)")
    d = {"cases": MC_CASES, "th": th, "mc": mc, "se": se, "nruns": NRUNS}
    json.dump(d, open(CACHE, "w"), indent=1)
    return d


def main():
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "axes.labelpad": 3.0, "text.color": INK,
        "xtick.color": SEC, "ytick.color": SEC,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "axes.grid": True, "grid.color": "#ecebe5", "grid.linewidth": 0.55,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.solid_capstyle": "round",
    })
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.2))
    (axA, axB, axC), (axD, axE, axF) = axes
    fig.subplots_adjust(left=0.055, right=0.975, bottom=0.085, top=0.945,
                        wspace=0.30, hspace=0.34)

    def tag(ax, s):
        ax.text(-0.155, 1.045, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minor(ax):
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which="both", top=False, right=False)

    MS = [2, 3, 4, 5, 6]
    cmap_m = plt.get_cmap("viridis")
    col = {m: cmap_m(0.08 + 0.78 * j / (len(MS) - 1)) for j, m in enumerate(MS)}

    # ------------------------------------------------ (a) C vs the naive (m-1)T
    lam = np.linspace(0.01, 2.5, 300)
    for m in MS:
        Cv = np.array([cascade_C(m, L, 1) for L in lam])
        naive = (m - 1) * lam / (1 + lam)
        axA.plot(lam, Cv, "-", color=col[m], lw=1.8, zorder=3, label=f"$m={m}$")
        axA.plot(lam, naive, "--", color=col[m], lw=1.0, alpha=0.55,
                 dashes=(3.4, 2.2), zorder=2)
    axA.axhline(1.0, color=MUTED, lw=0.8, ls=(0, (5, 3)), zorder=1)
    axA.set_xlim(0, 2.5); axA.set_ylim(0, 5.2)
    axA.set_xlabel(r"$\lambda$")
    axA.set_ylabel(r"cascade  $C(m,\lambda,\theta{=}1)$")
    axA.legend(loc="upper left", fontsize=8.2, handlelength=1.5,
               labelspacing=0.32, borderaxespad=0.4)
    axA.text(0.97, 0.06, "dashed: naive $(m-1)T$", transform=axA.transAxes,
             ha="right", va="bottom", fontsize=8.4, color=SEC, style="italic")
    minor(axA); tag(axA, "a")

    # ------------------------------------------------------ (b) the excess in %
    for m in MS[1:]:
        Cv = np.array([cascade_C(m, L, 1) for L in lam])
        naive = (m - 1) * lam / (1 + lam)
        axB.plot(lam, 100 * (Cv / naive - 1), "-", color=col[m], lw=1.8,
                 label=f"$m={m}$")
    for m, y in ((3, 11.1), (4, 21.5), (5, 31.0)):
        axB.plot([1.0], [y], "o", ms=4.6, mfc="white", mec=col[m], mew=1.5, zorder=5)
    axB.axvline(1.0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    for m, y, dy in ((3, 11.1, -8), (4, 21.5, -1), (5, 31.0, 6)):
        axB.annotate(f"{y}%", xy=(1.0, y), xytext=(10, dy), textcoords="offset points",
                     fontsize=8.0, color=col[m], va="center", ha="left",
                     bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.6))
    axB.set_xlim(0, 2.5); axB.set_ylim(0, None)
    axB.set_xlabel(r"$\lambda$")
    axB.set_ylabel(r"excess  $C/[(m{-}1)T]-1$   (%)")
    axB.legend(loc="upper left", fontsize=8.2, handlelength=1.5,
               labelspacing=0.32, borderaxespad=0.4)
    minor(axB); tag(axB, "b")

    # -------------------------------------------- (c) C over the (lambda, m) plane
    mg = np.arange(2, 9)
    lg = np.linspace(0.02, 2.0, 180)
    Z = np.array([[cascade_C(int(m), L, 1) for L in lg] for m in mg])
    edges_m = np.arange(1.5, 9.5, 1.0)
    dl = lg[1] - lg[0]
    edges_l = np.concatenate([lg - dl / 2, [lg[-1] + dl / 2]])
    # The field keeps the perceptually uniform, colour-vision-safe map the
    # rest of the figure is sampled from; the result is carried by the
    # contours, not by a hue trick, so the panel stays monotone in lightness
    # and reproduces in greyscale. Thin lines at C=2,4,6 make the field
    # quantitative; the C=1 locus -- where one group alone already replaces
    # its seed -- is the only heavy line, and every overlay is cased.
    # C spans a decade over this plane and the region that matters is the
    # low end, so the map is stretched by a power norm -- the colourbar
    # carries its own ticks, so the scale stays readable while C=1 lands
    # near the middle of the ramp instead of in the darkest corner.
    zmax = float(np.ceil(Z.max()))
    im = axC.pcolormesh(edges_l, edges_m, Z, cmap="viridis",
                        norm=mpl.colors.PowerNorm(gamma=0.55, vmin=0.0,
                                                  vmax=zmax),
                        rasterized=True)
    thin = [v for v in (2.0, 4.0) if v < zmax]
    cs_t = axC.contour(lg, mg, Z, levels=thin, colors="white", linewidths=1.0)
    cs_t.set_path_effects(cased(2.1))
    for t in axC.clabel(cs_t, fmt="%g", fontsize=7.0, colors="white"):
        t.set_path_effects(cased(1.5))
    cs = axC.contour(lg, mg, Z, levels=[1.0], colors="white", linewidths=2.2)
    cs.set_path_effects(cased(3.6))
    for t in axC.clabel(cs, fmt={1.0: r"$C=1$"}, fontsize=8.8, colors="white",
                        manual=[(1.20, 2.80)]):
        t.set_path_effects(cased(1.6))
    # m=2 stays below the line over the whole row: a pairwise-sized group
    # never replaces its own seed, however large lambda gets
    for (x, y), s in (((1.72, 2.0), r"$C<1$"), ((1.62, 7.5), r"$C>1$")):
        axC.text(x, y, s, fontsize=9, color="white", ha="center", va="center",
                 path_effects=cased(1.5))
    axC.set_yticks(mg)
    axC.set_xlabel(r"$\lambda$"); axC.set_ylabel(r"group size  $m$")
    axC.grid(False)
    cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.03,
                      ticks=np.arange(0, zmax + 0.5, 1.0))
    cb.set_label(r"$C$", fontsize=9); cb.ax.tick_params(labelsize=7.5)
    cb.ax.axhline(1.0, color="white", lw=1.5, path_effects=cased(2.8))
    cb.outline.set_linewidth(0.7); cb.outline.set_edgecolor("#c6c5bf")
    axC.tick_params(which="both", top=False, right=False)
    tag(axC, "c")

    # ------------------------------------- (d) the recursion lattice u(i,s), m=8
    # An annotated lattice is a different object from the continuous field of
    # (c), so it takes a different scale -- a single-hue sequential ramp, the
    # conventional choice for a table one reads numbers off, and likewise
    # monotone in lightness. Cell text follows the cell's own luminance, so
    # the values read at both ends of the ramp without touching the map.
    m0, lam0 = 6, 1.0
    U = u_lattice(m0, lam0, 1)
    cmap_d = plt.get_cmap("Blues")
    umax = float(np.nanmax(U))
    axD.set_facecolor("white")                 # unreachable cells stay blank
    axD.imshow(U, origin="lower", cmap=cmap_d, aspect="equal",
               interpolation="nearest", vmin=0.0, vmax=umax)
    for i in range(m0 + 1):
        for s in range(m0 + 1):
            if np.isnan(U[i, s]):
                continue
            r, g, b, _ = cmap_d(U[i, s] / umax)
            dark = 0.299 * r + 0.587 * g + 0.114 * b < 0.55
            # i=0 (died out) and s=0 (group exhausted) are absorbing, u=0
            # there by definition; keep them as structure, not as data
            absorbing = i == 0 or s == 0
            axD.text(s, i, "0" if absorbing else f"{U[i, s]:.2f}",
                     ha="center", va="center", zorder=3,
                     fontsize=6.4 if absorbing else 6.9,
                     alpha=0.55 if absorbing else 1.0,
                     color=("white" if dark else INK))
    axD.add_patch(Rectangle((m0 - 1.5, 0.5), 1.0, 1.0, fill=False, ec=CORAL,
                            lw=2.0, zorder=6, path_effects=cased(3.6, "white")))
    axD.annotate(rf"$C=u(1,m{{-}}1)={U[1, m0 - 1]:.3f}$",
                 xy=(m0 - 1, 1.58), xytext=(m0 - 1.6, 5.25),
                 fontsize=8.4, color=CORAL, ha="center",
                 path_effects=cased(2.4, "white"),
                 arrowprops=dict(arrowstyle="->", color=CORAL, lw=1.1,
                                 connectionstyle="arc3,rad=-0.16",
                                 shrinkA=3, shrinkB=3,
                                 path_effects=cased(2.8, "white")))
    axD.set_xlabel(r"susceptible  $s$"); axD.set_ylabel(r"infected  $i$")
    axD.set_xticks(range(0, m0 + 1)); axD.set_yticks(range(0, m0 + 1))
    axD.set_title(rf"$u(i,s)$,  $m={m0}$,  $\lambda={lam0:g}$", fontsize=9.2, pad=4)
    axD.grid(False)
    axD.tick_params(which="both", top=False, right=False)
    tag(axD, "d")

    # ------------------------------------------ (e) recursion vs Gillespie MC
    d = mc_data()
    th = np.array(d["th"]); mc = np.array(d["mc"]); se = np.array(d["se"])
    lo, hi = 0, max(th.max(), mc.max()) * 1.09
    axE.plot([lo, hi], [lo, hi], "-", color=MUTED, lw=0.9, zorder=1)
    axE.errorbar(th, mc, yerr=se, fmt="D", ms=4.6, color=BLUE, mfc=BLUE,
                 ecolor=BLUE, elinewidth=1.0, capsize=2.2, zorder=3)
    axE.set_xlim(lo, hi); axE.set_ylim(lo, hi); axE.set_aspect("equal")
    axE.set_xlabel(r"$C$  from the recursion (2.5)–(2.6)")
    axE.set_ylabel(r"$C$  from single-group Gillespie")
    axE.text(0.97, 0.10, f"{len(th)} settings\n$4\\times10^{{5}}$ runs each",
             transform=axE.transAxes, ha="right", va="bottom", fontsize=8.4,
             color=SEC)
    minor(axE); tag(axE, "e")

    # ------------------------------------------------------------- (f) MC pulls
    pull = (mc - th) / se
    xs = np.arange(1, len(pull) + 1)
    axF.axhspan(-2, 2, color=BAND2, lw=0, zorder=0)
    axF.axhspan(-1, 1, color=BAND1, lw=0, zorder=0)
    axF.axhline(0, color=MUTED, lw=0.9, zorder=1)
    axF.plot(xs, pull, "o", ms=5.0, color=BLUE, zorder=3)
    axF.set_xlim(0.4, len(pull) + 0.6); axF.set_ylim(-3.4, 3.4)
    axF.set_xticks(xs)
    axF.set_xticklabels([f"{m},{l:g}" for m, l, _ in d["cases"]],
                        rotation=45, fontsize=6.8, ha="right",
                        rotation_mode="anchor")
    axF.set_xlabel(r"setting  $(m,\lambda)$", labelpad=1.5)
    axF.set_ylabel(r"$(C_{\rm MC}-C_{\rm rec})/\sigma$")
    axF.text(0.5, 0.955, r"$\pm1\sigma,\ \pm2\sigma$", transform=axF.transAxes,
             ha="center", va="top", fontsize=8.4, color=SEC)
    axF.yaxis.set_minor_locator(AutoMinorLocator(2))
    axF.tick_params(which="both", top=False, right=False)
    axF.grid(axis="x", visible=False)
    tag(axF, "f")

    fig.savefig("figure_cascade.pdf", bbox_inches="tight")
    fig.savefig("figure_cascade.png", bbox_inches="tight", dpi=250)
    print("saved figure_cascade.png/.pdf")
    print(f"  worst MC pull: {np.abs(pull).max():.2f} sigma")
    for m in (3, 4, 5):
        c = cascade_C(m, 1.0, 1); n = (m - 1) * 0.5
        print(f"  m={m}, lam=1: C={c:.4f}  naive={n:.4f}  excess={100*(c/n-1):.1f}%")


if __name__ == "__main__":
    main()
