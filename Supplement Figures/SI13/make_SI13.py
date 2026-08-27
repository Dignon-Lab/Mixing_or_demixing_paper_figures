#!/usr/bin/env python3
"""SI Figure 13 -- LAF-1 V2 + LAF-1 V3: Rg relaxation time, dilute-phase
concentration, and demixing index vs composition.
Run: python make_SI13.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
CA, CB = "#2166ac", "#b2182b"
name1, name2 = "LAF-1 V2", "LAF-1 V3"
PCT = np.array([5, 25, 50, 75, 95])

TAU_V2 = np.array([5.009, 4.883, 5.011, 4.686, 4.726]); TAU_V2_SEM = np.array([1.290, 1.047, 1.182, 1.230, 2.417])
TAU_V3 = np.array([6.339, 6.321, 6.655, 6.211, 7.263]); TAU_V3_SEM = np.array([1.784, 1.770, 1.898, 1.659, 3.228])
V2_DIL = np.array([0.2073, 4.0291, 2.8179, 3.3362, 28.2524])
V3_DIL = np.array([1.5964, 3.7051, 0.4851, 0.2378, 0.0178])

z = np.load(os.path.join(HERE, "phase_voxel_LAF_1_V2__LAF_1_V3.npz"), allow_pickle=True)
order = np.argsort(np.asarray(z["pct"], int))
PCT_D = np.asarray(z["pct"], int)[order]; DEMIX_D = np.asarray(z["D"], float)[order]

mpl.rcParams.update({"font.size": 11, "font.family": "sans-serif", "mathtext.default": "regular"})
fig = plt.figure(figsize=(15, 4.6))
gs = GridSpec(1, 3, figure=fig, wspace=0.34, left=0.06, right=0.97, top=0.85, bottom=0.15)

axG = fig.add_subplot(gs[0, 0])
axG.errorbar(PCT, TAU_V2, yerr=TAU_V2_SEM, fmt="o-", color=CA, lw=1.9, capsize=3, label="A = " + name1)
axG.errorbar(PCT, TAU_V3, yerr=TAU_V3_SEM, fmt="s-", color=CB, lw=1.9, capsize=3, label="B = " + name2)
axG.set_xlabel("% V2 (component A)"); axG.set_ylabel(r"$\tau_{R_g}$ (ns)")
axG.set_title("G   Rg relaxation time", loc="left", fontweight="bold")
axG.axvline(50, color="0.8", ls=":", lw=1.0); axG.legend(fontsize=8); axG.grid(alpha=0.25)

axH = fig.add_subplot(gs[0, 1])
axH.plot(PCT, V2_DIL, "o-", color=CA, lw=1.9, label="A = " + name1)
axH.plot(PCT, V3_DIL, "s-", color=CB, lw=1.9, label="B = " + name2)
axH.set_xlabel("% V2 (component A)"); axH.set_ylabel("dilute-phase conc. ($\\mu$M)")
axH.set_title("H   Saturation concentration", loc="left", fontweight="bold")
axH.axvline(50, color="0.8", ls=":", lw=1.0); axH.legend(fontsize=8); axH.grid(alpha=0.25)

axI = fig.add_subplot(gs[0, 2])
axI.axhspan(0.0, 0.2, color="green", alpha=0.08, label="mixed (D < 0.2)")
axI.plot(PCT_D, DEMIX_D, "o-", color="#333333", lw=1.9, zorder=5, label="demixing index D")
for x, y in zip(PCT_D, DEMIX_D):
    axI.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
axI.set_xlabel("% V2 (component A)"); axI.set_ylabel("demixing index  D")
axI.set_ylim(-0.02, 1.02)
axI.set_title("Demixing index vs composition", loc="left", fontweight="bold")
axI.axvline(50, color="0.85", ls=":", lw=1.0); axI.legend(fontsize=8); axI.grid(alpha=0.25)

fig.suptitle(f"SI Figure 13 -- {name1} + {name2}: dynamics and demixing", fontsize=13, fontweight="bold", y=0.99)
fig.savefig(os.path.join(HERE, "FigureS13.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS13.png"), dpi=300, bbox_inches="tight")
print("saved FigureS13.png")
