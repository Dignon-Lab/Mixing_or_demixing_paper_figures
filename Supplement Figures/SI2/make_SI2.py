#!/usr/bin/env python3
"""SI Figure 2 -- single-chain and single-component characterisation (extra panels for Fig 1).
Run: python make_SI2.py"""
import os
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(HERE, "figure1_data.csv")).sort_values("SCD").reset_index(drop=True)
df["rg_sc"] = df["sc_Rg_nm"];      df["rg_sc_e"] = df["sc_Rg_sem_nm"]
df["nu_e"]  = df["nu_sem"]
df["tau"]   = df["incond_tau_ns"]; df["tau_e"] = df["incond_tau_sem_ns"]
df["D"]     = df["incond_D"];      df["D_e"]   = df["incond_D_sem"]
df["dense"] = df["dense_used_mM"]; df["dense_e"] = df["dense_used_sem"]

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

def pear(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    return np.corrcoef(x[m], y[m])[0, 1]

mpl.rcParams.update({"font.size": 11, "font.family": "sans-serif", "mathtext.default": "regular"})
fig = plt.figure(figsize=(11, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32, left=0.09, right=0.93, top=0.93, bottom=0.08)

axA = fig.add_subplot(gs[0, 0])
sctr(axA, df["SHD"], df["rg_sc"], yerr=df["rg_sc_e"])
axA.text(0.05, 0.93, f"r = {pear(df['SHD'], df['rg_sc']):.2f}", transform=axA.transAxes, va="top")
axA.set_xlabel("SHD"); axA.set_ylabel(r"single-chain $\langle R_g \rangle$ (nm)")
axA.set_title("A   Chain size vs SHD", loc="left", fontweight="bold"); axA.grid(alpha=0.25)

axB = fig.add_subplot(gs[0, 1])
sctr(axB, df["SCD"], df["nu"], yerr=df["nu_e"])
axB.axhline(0.5, ls=":", color="0.5", lw=1.0)
axB.text(0.03, 0.06, r"$\theta$ ($\nu$=0.5)", transform=axB.transAxes, fontsize=8, color="0.4")
axB.set_xlabel("SCD"); axB.set_ylabel(r"Multi-chain scaling exponent $\nu$")
axB.set_title("B   Scaling exponent vs SCD", loc="left", fontweight="bold"); axB.grid(alpha=0.25)

axC = fig.add_subplot(gs[1, 0])
sctr(axC, df["SCD"], df["tau"], yerr=df["tau_e"])
axC.text(0.95, 0.93, f"r = {pear(df['SCD'], df['tau']):.2f}", transform=axC.transAxes, va="top", ha="right")
axC.set_xlabel("SCD"); axC.set_ylabel(r"Multi-chain $\tau_{R_g}$ (ns)")
axC.set_title("C   Chain relaxation time vs SCD", loc="left", fontweight="bold"); axC.grid(alpha=0.25)

axD = fig.add_subplot(gs[1, 1])
sctr(axD, df["dense"], df["D"], xerr=df["dense_e"], yerr=df["D_e"])
axD.text(0.95, 0.93, f"r = {pear(df['dense'], df['D']):.2f}", transform=axD.transAxes, va="top", ha="right")
axD.set_xlabel("dense-phase conc. (mM)"); axD.set_ylabel(r"Multi-chain diffusivity (cm$^2$/s)")
axD.set_title("D   Dynamics vs condensate density", loc="left", fontweight="bold"); axD.grid(alpha=0.25)

sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=[axA, axB, axC, axD], pad=0.02, aspect=40, shrink=0.6)
cb.set_label("SCD")

fig.savefig(os.path.join(HERE, "FigureS2.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS2.png"), dpi=300, bbox_inches="tight")
print("saved FigureS2.png")
