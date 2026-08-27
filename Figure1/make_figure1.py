#!/usr/bin/env python3
"""Figure 1 -- design space, single-chain dimensions, single-component phase behavior.
Run: python make_figure1.py"""
import os
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(HERE, "figure1_data.csv")).sort_values("SCD").reset_index(drop=True)
df["rg_sc"] = df["sc_Rg_nm"];       df["rg_sc_e"] = df["sc_Rg_sem_nm"]
df["rg_in"] = df["incond_Rg_nm"];   df["rg_in_e"] = df["incond_Rg_sem_nm"]
df["D"]     = df["incond_D"];       df["D_e"]     = df["incond_D_sem"]
df["dense"] = df["dense_used_mM"];  df["dense_e"] = df["dense_used_sem"]
df["csat"]  = df["dilute_KDE_uM"];  df["csat_e"]  = df["dilute_KDE_sem"]

cmap = plt.cm.RdYlBu
norm = mpl.colors.Normalize(vmin=-22, vmax=2)
col = lambda v: cmap(norm(v))
is_wt = (df["key"] == "LAF-1_IDR_WT").values

def sctr(ax, x, y, xerr=None, yerr=None, s=60):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="none", ecolor="0.55",
                elinewidth=0.7, capsize=2, zorder=2)
    ax.scatter(x[~is_wt], y[~is_wt], c=[col(v) for v in df["SCD"].values[~is_wt]],
               s=s, edgecolors="k", linewidths=0.6, zorder=3)
    ax.scatter(x[is_wt], y[is_wt], c="k", s=s, edgecolors="k", linewidths=0.6, zorder=4)

def shade_nollps(ax, dense_cut=12.0, csat_cut=300.0):
    y0, y1 = ax.get_ylim()
    if ax.get_yscale() == "log":
        ax.axhspan(csat_cut, y1, color="0.85", alpha=0.55, zorder=0, lw=0); ty = 0.90
    else:
        ax.axhspan(y0, dense_cut, color="0.85", alpha=0.55, zorder=0, lw=0); ty = 0.12
    ax.set_ylim(y0, y1)
    ax.text(0.03, ty, "No LLPS\n($T_c$ < 298 K)", transform=ax.transAxes,
            ha="left", va="center", fontsize=8, color="0.35")

def add_yx(ax):
    lo = min(ax.get_xlim()[0], ax.get_ylim()[0]); hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.6, zorder=1, label="y = x")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)

def pear(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    return np.corrcoef(x[m], y[m])[0, 1]

mpl.rcParams.update({"font.size": 11, "font.family": "sans-serif", "mathtext.default": "regular"})
fig = plt.figure(figsize=(15, 9.5))
gs = GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.36, left=0.06, right=0.97, top=0.94, bottom=0.07)

axA = fig.add_subplot(gs[0, 0])
sctr(axA, df["SCD"], df["SHD"])
for xi, yi, l in zip(df["SCD"], df["SHD"], df["label"]):
    axA.annotate(l, (xi, yi), fontsize=6.5, ha="center", va="bottom", xytext=(0, 3), textcoords="offset points")
axA.set_xlabel("SCD (charge patterning)"); axA.set_ylabel("SHD (hydropathy patterning)")
axA.set_title("A   Sequence design space", loc="left", fontweight="bold"); axA.grid(alpha=0.25)

axB = fig.add_subplot(gs[0, 1])
sctr(axB, df["SCD"], df["rg_sc"], yerr=df["rg_sc_e"])
axB.text(0.05, 0.93, f"r = {pear(df['SCD'], df['rg_sc']):.2f}", transform=axB.transAxes, va="top")
axB.set_xlabel("SCD"); axB.set_ylabel(r"single-chain $\langle R_g \rangle$ (nm)")
axB.set_title("B   Single-chain size vs SCD", loc="left", fontweight="bold"); axB.grid(alpha=0.25)

axC = fig.add_subplot(gs[0, 2])
sctr(axC, df["SCD"], df["dense"], yerr=df["dense_e"])
shade_nollps(axC)
axC.set_xlabel("SCD"); axC.set_ylabel("Dense-phase conc. (mM)")
axC.set_title("C   Dense-phase concentration vs SCD", loc="left", fontweight="bold"); axC.grid(alpha=0.25)

axD = fig.add_subplot(gs[1, 0])
axD.errorbar(df["SCD"], df["csat"], yerr=df["csat_e"], fmt="none", ecolor="0.55", elinewidth=0.7, capsize=2, zorder=2)
axD.scatter(df["SCD"][~is_wt], df["csat"][~is_wt], c=[col(v) for v in df["SCD"].values[~is_wt]], s=60, edgecolors="k", linewidths=0.6, zorder=3)
axD.scatter(df["SCD"][is_wt], df["csat"][is_wt], c="k", s=60, edgecolors="k", zorder=4)
axD.set_yscale("log"); shade_nollps(axD)
axD.set_xlabel("SCD"); axD.set_ylabel(r"Dilute conc. $c_{sat}$ ($\mu$M)")
axD.set_title("D   Saturation concentration vs SCD", loc="left", fontweight="bold"); axD.grid(alpha=0.25, which="both")

axE = fig.add_subplot(gs[1, 1])
sctr(axE, df["rg_sc"], df["rg_in"], xerr=df["rg_sc_e"], yerr=df["rg_in_e"])
add_yx(axE)
axE.set_xlabel(r"single-chain $\langle R_g \rangle$ (nm)"); axE.set_ylabel(r"Multi-chain $\langle R_g \rangle$ (nm)")
axE.set_title("E   Chains swell to a common size", loc="left", fontweight="bold")
axE.legend(loc="lower right", fontsize=9, frameon=False); axE.grid(alpha=0.25)

axF = fig.add_subplot(gs[1, 2])
sctr(axF, df["SCD"], df["D"], yerr=df["D_e"])
axF.text(0.05, 0.93, f"r = {pear(df['SCD'], df['D']):.2f}", transform=axF.transAxes, va="top")
axF.set_xlabel("SCD"); axF.set_ylabel(r"Multi-chain diffusivity (cm$^2$/s)")
axF.set_title("F   Dynamics vs SCD", loc="left", fontweight="bold"); axF.grid(alpha=0.25)

sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=[axA, axB, axC, axD, axE, axF], pad=0.01, aspect=40, shrink=0.6)
cb.set_label("SCD")

fig.savefig(os.path.join(HERE, "Figure1.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "Figure1.png"), dpi=300, bbox_inches="tight")
print("saved Figure1.png")
