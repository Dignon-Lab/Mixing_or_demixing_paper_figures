#!/usr/bin/env python3
"""SI Figure 12 -- LAF-1 V2 + LAF-1 V3 characterisation: contact maps and
chain properties vs composition.
Run: python make_SI12.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
CA, CB = "#2166ac", "#b2182b"
name1, name2 = "LAF-1 V2", "LAF-1 V3"
PCT = np.array([5, 25, 50, 75, 95])

RG_V2 = np.array([36.20, 36.62, 36.76, 36.83, 36.47]); RG_V2_SEM = np.array([1.11, 0.33, 0.26, 0.22, 0.10])
RG_V3 = np.array([36.36, 36.41, 36.56, 36.50, 36.36]); RG_V3_SEM = np.array([0.13, 0.08, 0.39, 0.47, 1.61])
NU_V2 = np.array([0.5479, 0.5502, 0.5508, 0.5514, 0.5503]); NU_V2_SEM = np.array([0.0017]*5)
NU_V3 = np.array([0.5489, 0.5492, 0.5502, 0.5506, 0.5491]); NU_V3_SEM = np.array([0.0018]*5)
D_V2 = np.array([3.433e-6, 3.664e-6, 3.762e-6, 4.018e-6, 4.420e-6])
D_V3 = np.array([3.123e-6, 3.279e-6, 3.174e-6, 3.691e-6, 3.285e-6])
def _e(a): return np.where(np.isfinite(np.asarray(a, float)), a, 0.0)

cm = np.load(os.path.join(HERE, "contact_maps_LAF1V2__LAF1V3.npz"), allow_pickle=True)
AA, AB, BB = cm["AA"], cm["AB"], cm["BB"]

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
    ax.errorbar(PCT, ya, yerr=_e(ya_e), fmt="o-", color=CA, lw=1.9, capsize=3, label=f"A = {name1}")
    ax.errorbar(PCT, yb, yerr=_e(yb_e), fmt="s-", color=CB, lw=1.9, capsize=3, label=f"B = {name2}")
    ax.set_xlabel("% V2 (component A)"); ax.set_ylabel(ylabel)
    ax.set_title(f"{letter}   {title}", loc="left", fontweight="bold")
    if loga: ax.set_yscale("log")
    ax.axvline(50, color="0.8", ls=":", lw=1.0); ax.legend(fontsize=8); ax.grid(alpha=0.25)

tw(fig.add_subplot(gs[1, 0]), RG_V2, RG_V2_SEM, RG_V3, RG_V3_SEM, r"$\langle R_g\rangle$ (A)", "Radius of gyration", "D")
tw(fig.add_subplot(gs[1, 1]), NU_V2, NU_V2_SEM, NU_V3, NU_V3_SEM, r"scaling exponent $\nu$", "Internal scaling", "E")
tw(fig.add_subplot(gs[1, 2]), D_V2, None, D_V3, None, r"D (cm$^2$/s)", "Translational diffusivity", "F", loga=True)

fig.suptitle(f"SI Figure 12 -- {name1} + {name2} characterisation", fontsize=13, fontweight="bold", y=0.975)
fig.savefig(os.path.join(HERE, "FigureS12.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS12.png"), dpi=300, bbox_inches="tight")
print("saved FigureS12.png")
