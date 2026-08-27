"""Figure 2: the group closure is correct, two independent ways.

  (a) full transient: I(t) from the closure vs an exact Gillespie simulation
      (N=6000, 400 runs), the two indistinguishable;
  (b) the subcritical determination of the threshold: the reciprocal outbreak
      size eps/R(inf) = 1/chi falls linearly to zero, its x-intercept is lambda_c
      (three example configurations, extrapolation dashed to the zero);
  (c) the eight configurations' subcritical lambda_c against the spectral value
      rho(N)=1, on the identity line;
  (d) the deviation in units of the regression standard error.

All data cached: figure2_data.json (a) and figure_data.json (b-d).
"""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE = "#1C9B8E", "#E76F51", "#2F5FD0"
PT = "#2F5FD0"
BAND2, BAND1, BAND = "#eef2f6", "#dbe3ec", "#e4e3db"


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
    fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(9.3, 7.0))
    fig.subplots_adjust(left=0.088, right=0.975, bottom=0.085, top=0.955,
                        wspace=0.26, hspace=0.30)

    def tag(ax, s):
        ax.text(-0.155, 1.04, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minor(ax, nx=2, ny=2):
        ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
        ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
        ax.tick_params(which="both", top=False, right=False)

    # ------------------------------------------------ (a) I(t): closure vs sim
    d2 = json.load(open("figure2_data.json"))
    TG = np.array(d2["t"]); Im = np.array(d2["Im"]); Ici = np.array(d2["Ici"])
    Ic = np.array(d2["Ic"])
    axA.fill_between(TG, Im - Ici, Im + Ici, color=BAND, lw=0, zorder=1)
    axA.plot(TG, Im, "-", color=TEAL, lw=1.6, zorder=3, label="simulation")
    axA.plot(TG, Ic, "--", color=CORAL, lw=1.7, dashes=(4, 2.4), zorder=4,
             label="group closure")
    axA.set_xlim(0, 20); axA.set_ylim(bottom=0)
    axA.set_xlabel(r"time $t$"); axA.set_ylabel(r"infected fraction $I(t)$")
    minor(axA)
    axA.legend(loc="upper right", fontsize=8.6, handlelength=1.8,
               handletextpad=0.6, labelspacing=0.5)
    tag(axA, "a")

    # ------------------------------------------- (b) subcritical extrapolation
    data = json.load(open("figure_data.json"))
    by = {r["name"]: r for r in data}
    show = [("2-layer m=(3,4)", TEAL, r"2-layer $m=(3,4)$"),
            ("6-reg. pairwise", CORAL, "6-reg. pairwise"),
            ("m=3 deg{2,3}", BLUE, r"1-layer $m=3$")]
    axB.axhline(0, color=MUTED, lw=0.8, zorder=1)
    handles = []
    for name, col, lab in show:
        r = by[name]; lcr = r["lc_rho"]
        xn = np.array(r["lams"]) / lcr; ys = np.array(r["ys"]); a, b = r["fit"]
        xext = r["lc_extrap"] / lcr
        axB.plot(np.linspace(xn[0], xext, 60), a + b * np.linspace(xn[0], xext, 60) * lcr,
                 "--", color=col, lw=1.0, alpha=0.75, dashes=(4, 2.4), zorder=2)
        ln, = axB.plot(np.linspace(xn[0], xn[-1], 40),
                       a + b * np.linspace(xn[0], xn[-1], 40) * lcr, "-",
                       color=col, lw=1.7, zorder=3, label=lab)
        axB.plot(xn, ys, "o", ms=4.0, mfc=col, mec="white", mew=0.7, zorder=4)
        handles.append(ln)
    axB.plot([1.0], [0.0], "o", ms=6, mfc="white", mec=SEC, mew=1.0, zorder=5)
    axB.set_xlim(0.876, 1.026); axB.set_ylim(-0.009, 0.098)
    axB.set_xlabel(r"$\lambda/\lambda_c$"); axB.set_ylabel(r"$\varepsilon/R(\infty)=1/\chi$")
    axB.xaxis.set_major_locator(MultipleLocator(0.05)); minor(axB, nx=5)
    axB.legend(handles=handles, loc="upper right", fontsize=8.0,
               handlelength=1.5, handletextpad=0.55, labelspacing=0.45)
    tag(axB, "b")

    # ------------------------------------------- (c) extrap vs rho(N)=1
    lcr = np.array([r["lc_rho"] for r in data])
    lce = np.array([r["lc_extrap"] for r in data])
    sig = np.array([r["sigma"] for r in data])
    lo, hi = 0.06, 0.42
    axC.plot([lo, hi], [lo, hi], "-", color=MUTED, lw=1.0, zorder=1)
    axC.errorbar(lcr, lce, yerr=sig, fmt="D", ms=4.5, mfc=PT, mec="white", mew=0.7,
                 ecolor=PT, elinewidth=1.0, capsize=2.2, capthick=0.8, zorder=3)
    axC.set_xlim(lo, hi); axC.set_ylim(lo, hi); axC.set_aspect("equal")
    axC.set_xlabel(r"$\lambda_c$  from $\rho(N)=1$")
    axC.set_ylabel(r"$\lambda_c$  from extrapolation")
    axC.xaxis.set_major_locator(MultipleLocator(0.1))
    axC.yaxis.set_major_locator(MultipleLocator(0.1)); minor(axC)
    maxrel = np.max(np.abs(lce - lcr) / lcr)
    axC.text(0.96, 0.06, "max rel. dev.\n" + r"$6.5\times10^{-4}$", transform=axC.transAxes,
             fontsize=8.4, color=SEC, ha="right", va="bottom", linespacing=1.4)
    tag(axC, "c")

    # ------------------------------------------- (d) pull
    pull = (lce - lcr) / sig
    order = np.argsort(lcr)
    axD.axhspan(-2, 2, color=BAND2, lw=0, zorder=0)
    axD.axhspan(-1, 1, color=BAND1, lw=0, zorder=0)
    axD.axhline(0, color=MUTED, lw=0.9, zorder=1)
    axD.plot(lcr[order], pull[order], "o", ms=5.0, mfc=PT, mec="white", mew=0.8, zorder=3)
    axD.set_ylim(-3, 3); axD.set_xlim(0.06, 0.42)
    axD.set_xlabel(r"$\lambda_c$"); axD.set_ylabel(r"$(\lambda_c^{\mathrm{extr}}-\lambda_c^{\rho})/\sigma$")
    axD.xaxis.set_major_locator(MultipleLocator(0.1)); axD.set_yticks([-2, -1, 0, 1, 2])
    axD.xaxis.set_minor_locator(AutoMinorLocator(2)); axD.tick_params(which="both", top=False, right=False)
    axD.text(0.5, 0.955, r"$\pm1\sigma,\ \pm2\sigma$", transform=axD.transAxes,
             fontsize=8.4, color=SEC, ha="center", va="top")
    tag(axD, "d")

    fig.savefig("figure2_validation.pdf", bbox_inches="tight")
    fig.savefig("figure2_validation.png", bbox_inches="tight", dpi=200)
    print("saved figure2_validation.png/.pdf ; max pull=%.2f, max rel dev=%.2e"
          % (np.max(np.abs(pull)), maxrel))


if __name__ == "__main__":
    main()
