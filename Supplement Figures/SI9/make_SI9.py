#!/usr/bin/env python3
"""SI Figure 9 -- Near-zero self-associator characterisation for LAF-1 V2 + protein13
(contact maps + chain properties vs composition).
Run: python make_SI5.py"""
import os
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, "contact_maps_LAF1V2__protein13.npz")
MASTER = os.path.join(HERE, "two_comp_master.csv")
CA, CB = "#2166ac", "#b2182b"

m = pd.read_csv(MASTER); m = m[m.partner == "protein13"].sort_values("pctV2")
PCT = m.pctV2.values

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
    ax.errorbar(PCT, ya, yerr=ya_e, fmt="o-", color=CA, lw=1.9, capsize=3, label=f"A = {name1}")
    ax.errorbar(PCT, yb, yerr=yb_e, fmt="D-", color=CB, lw=1.9, capsize=3, label=f"B = {name2}")
    ax.set_xlabel("% V2 (component A)"); ax.set_ylabel(ylabel)
    ax.set_title(f"{letter}   {title}", loc="left", fontweight="bold")
    if loga: ax.set_yscale("log")
    ax.axvline(50, color="0.8", ls=":", lw=1.0)
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

tw(fig.add_subplot(gs[1, 0]), m.Rg_V2.values, m.Rg_V2_sem.values, m.Rg_B.values, m.Rg_B_sem.values,
   r"$\langle R_g\rangle$ (A)", "Radius of gyration", "D")
tw(fig.add_subplot(gs[1, 1]), m.nu_V2.values, m.nu_V2_sem.values, m.nu_B.values, m.nu_B_sem.values,
   r"scaling exponent $\nu$", "Internal scaling", "E")
tw(fig.add_subplot(gs[1, 2]), m.D_V2.values, m.D_V2_sem.values, m.D_B.values, m.D_B_sem.values,
   r"D (cm$^2$/s)", "Translational diffusivity", "F", loga=True)

fig.suptitle(f"SI Figure 9 -- Near-zero self-associator characterisation:  {name1} + {name2}", fontsize=13, fontweight="bold", y=0.975)
fig.savefig(os.path.join(HERE, "FigureS9.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS9.png"), dpi=300, bbox_inches="tight")
print("saved FigureS9.png")
