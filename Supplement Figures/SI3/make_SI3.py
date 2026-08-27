#!/usr/bin/env python3
"""SI Figure 3 -- Ideal-mixing characterisation for LAF-1 V2 + protein9
(contact maps + chain properties vs composition).
Run: python make_SI3.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, "contact_maps_LAF1V2__protein9.npz")
CA, CB = "#2166ac", "#1b7837"

PCT        = np.array([5, 25, 50, 75, 95])
RG_V2      = np.array([36.893, 36.841, 36.606, 36.752, 36.674])
RG_V2_SEM  = np.array([0.3577, 0.2272, 0.0983, 0.1117, 0.0860])
RG_P9      = np.array([36.920, 36.915, 37.115, 36.789, 36.962])
RG_P9_SEM  = np.array([0.0843, 0.0883, 0.1328, 0.1596, 0.3965])
NU_V2      = np.array([0.5798, 0.5818, 0.5794, 0.5803, 0.5795])
NU_V2_SEM  = np.array([0.0030, 0.0016, 0.0009, 0.0010, 0.0008])
NU_P9      = np.array([0.5802, 0.5803, 0.5820, 0.5797, 0.5798])
NU_P9_SEM  = np.array([0.0006, 0.0006, 0.0011, 0.0013, 0.0032])
D_V2       = np.array([3.943e-6, 2.884e-6, 2.953e-6, 3.279e-6, 3.239e-6])
D_V2_SEM   = np.array([np.nan,   1.725e-7, 1.256e-7, 1.202e-7, 1.189e-7])
D_P9       = np.array([3.083e-6, 2.899e-6, 2.785e-6, 3.335e-6, 3.296e-6])
D_P9_SEM   = np.array([1.440e-7, 3.367e-8, 9.905e-8, 3.840e-7, np.nan])

d = np.load(NPZ, allow_pickle=True)
AA, AB, BB = d["AA"], d["AB"], d["BB"]
name1, name2 = str(d["name1"]), str(d["name2"])

mpl.rcParams.update({"font.size": 10, "axes.titlesize": 11})
fig = plt.figure(figsize=(13.5, 8.4))
gs = GridSpec(2, 3, figure=fig, hspace=0.44, wspace=0.42, left=0.06, right=0.9, top=0.90, bottom=0.09)

vmax = np.percentile(np.concatenate([AA.ravel(), AB.ravel(), BB.ravel()]), 99)
maps = [(AA, f"A   A-A   {name1}-{name1}"), (AB, f"B   A-B   {name1}-{name2}"), (BB, f"C   B-B   {name2}-{name2}")]
ims = []
for k, (M, t) in enumerate(maps):
    ax = fig.add_subplot(gs[0, k])
    im = ax.imshow(M, origin="lower", cmap="magma", vmin=0, vmax=vmax, aspect="equal", interpolation="nearest")
    ims.append(im)
    ax.set_title(t, loc="left", fontweight="bold", fontsize=10.5)
    ax.set_xlabel("residue index  j")
    if k == 0:
        ax.set_ylabel("residue index  i")
    ax.text(0.96, 0.05, f"max = {M.max():.1e}", transform=ax.transAxes, ha="right", va="bottom", color="w", fontsize=8)
cax = fig.add_axes([0.915, 0.565, 0.010, 0.30])
cb = fig.colorbar(ims[0], cax=cax)
cb.set_label("contact frequency\n(per chain-pair per frame)", fontsize=8.5)

def tw(ax, ya, ya_e, yb, yb_e, ylabel, title, letter, loga=False):
    ea = np.where(np.isfinite(ya_e), ya_e, 0.0)
    eb = np.where(np.isfinite(yb_e), yb_e, 0.0)
    ax.errorbar(PCT, ya, yerr=ea, fmt="o-", color=CA, lw=1.9, capsize=3, label=f"A = {name1}")
    ax.errorbar(PCT, yb, yerr=eb, fmt="D-", color=CB, lw=1.9, capsize=3, label=f"B = {name2}")
    ax.set_xlabel("% V2 (component A)"); ax.set_ylabel(ylabel)
    ax.set_title(f"{letter}   {title}", loc="left", fontweight="bold")
    if loga: ax.set_yscale("log")
    ax.axvline(50, color="0.8", ls=":", lw=1.0)
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

tw(fig.add_subplot(gs[1, 0]), RG_V2, RG_V2_SEM, RG_P9, RG_P9_SEM, r"$\langle R_g\rangle$ (A)", "Radius of gyration", "D")
tw(fig.add_subplot(gs[1, 1]), NU_V2, NU_V2_SEM, NU_P9, NU_P9_SEM, r"scaling exponent $\nu$", "Internal scaling", "E")
tw(fig.add_subplot(gs[1, 2]), D_V2, D_V2_SEM, D_P9, D_P9_SEM, r"D (cm$^2$/s)", "Translational diffusivity", "F", loga=True)

fig.suptitle(f"SI Figure 3 -- Ideal-mixing characterisation:  {name1} + {name2}", fontsize=13, fontweight="bold", y=0.975)
fig.savefig(os.path.join(HERE, "FigureS3.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS3.png"), dpi=300, bbox_inches="tight")
print("saved FigureS3.png")
