#!/usr/bin/env python3
"""SI Figure 4 -- LAF-1 V2 + protein9 occupy the same dense slab and stay mixed
at every composition (298 K): demixing index D vs composition.
Note: the raw-trajectory density-profile panel is not reproduced here (needs
the production trajectory, not distributed in this repository).
Run: python make_SI4.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, "phase_voxel_LAF_1_V2__protein9.npz")

z = np.load(NPZ, allow_pickle=True)
pct = np.asarray(z["pct"], dtype=int)
D = np.asarray(z["D"], float)
order = np.argsort(pct)

fig, ax = plt.subplots(figsize=(6.0, 4.6))
ax.plot(pct[order], D[order], "o-", color="#7b3294", lw=1.9, ms=8)
ax.axhline(0.2, color="0.6", ls=":", lw=1.0)
ax.text(0.03, 0.90, "well-mixed threshold  D < 0.2", transform=ax.transAxes, fontsize=8.5, color="0.4")
ax.set_xlabel("% V2 (component A)")
ax.set_ylabel("demixing index  D")
ax.set_ylim(0, 1)
ax.set_title("SI Figure 4B -- Demixing index vs composition\n(V2 + protein9)", loc="left", fontweight="bold")
ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "FigureS4.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS4.png"), dpi=300, bbox_inches="tight")
print("saved FigureS4.png")
