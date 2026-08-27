#!/usr/bin/env python3
"""SI Figure 6 -- Density profiles and demixing index for LAF-1 V2 + protein2
across composition (298 K). Note: the raw-trajectory density-profile panel is
not reproduced here (needs the production trajectory, not distributed in this
repository); the demixing index panel is fully reproducible from the DDM voxel
export.
Run: python make_SI6.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, "phase_voxel_LAF_1_V2__protein2.npz")

z = np.load(NPZ, allow_pickle=True)
pct = np.asarray(z["pct"], dtype=int)
D = np.asarray(z["D"], float)
order = np.argsort(pct)

fig, ax = plt.subplots(figsize=(6.0, 4.6))
ax.plot(pct[order], D[order], "o-", color="#018571", lw=1.9, ms=8)
ax.axhline(0.2, color="0.6", ls=":", lw=1.0)
ax.text(0.03, 0.90, "well-mixed threshold  D < 0.2", transform=ax.transAxes, fontsize=8.5, color="0.4")
ax.set_xlabel("% V2 (component A)")
ax.set_ylabel("demixing index  D")
ax.set_ylim(0, 1)
ax.set_title("SI Figure 6B -- Demixing index vs composition\n(V2 + protein2)", loc="left", fontweight="bold")
ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "FigureS6.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS6.png"), dpi=300, bbox_inches="tight")
print("saved FigureS6.png")
