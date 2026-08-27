#!/usr/bin/env python3
"""Figure 7 -- LAF-1 V2 (A) + LAF-1 V3 (B), 298 K.
Run: python make_figure7.py"""
import os
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
F1  = os.path.join(HERE, "figure1_data.csv")
NPZ = os.path.join(HERE, "phase_voxel_LAF_1_V2__LAF_1_V3.npz")

CA, CB, CAB = "#2166ac", "#b2182b", "#7b3294"
name1, name2 = "LAF-1 V2", "LAF-1 V3"
PCT = np.array([5, 25, 50, 75, 95])

V2_DENSE = np.array([1.3123, 5.9693, 11.6858, 16.7984, 20.9949])
V3_DENSE = np.array([24.3319, 19.0379, 12.3210, 5.7454, 1.1162])

cm = np.load(os.path.join(HERE, "contact_maps_LAF1V2__LAF1V3.npz"), allow_pickle=True)
sAA, sAB, sBB = float(cm["sAA"]), float(cm["sAB"]), float(cm["sBB"])
cvc = np.load(os.path.join(HERE, "contacts_vs_composition_LAF1V2__LAF1V3.npz"), allow_pickle=True)
pctA = cvc["pctA"]

f1 = pd.read_csv(F1)
scd_v2 = float(f1.loc[f1.key == "LAF-1_V2", "SCD"].iloc[0]); shd_v2 = float(f1.loc[f1.key == "LAF-1_V2", "SHD"].iloc[0])
scd_v3 = float(f1.loc[f1.key == "LAF-1_V3", "SCD"].iloc[0]); shd_v3 = float(f1.loc[f1.key == "LAF-1_V3", "SHD"].iloc[0])

def draw_phase_diagram(ax, npz_path, add_inset=True, legend_fs=8.0):
    z = np.load(npz_path, allow_pickle=True)
    pct = np.asarray(z["pct"], int)
    dil1, den1 = np.asarray(z["dil1"], float), np.asarray(z["den1"], float)
    dil2, den2 = np.asarray(z["dil2"], float), np.asarray(z["den2"], float)
    n1, n2 = str(z["name1"]), str(z["name2"])
    den1e = np.asarray(z["den1_err"], float) if "den1_err" in z else np.zeros_like(den1)
    den2e = np.asarray(z["den2_err"], float) if "den2_err" in z else np.zeros_like(den2)
    order = np.argsort(pct)
    colors = plt.cm.plasma(np.linspace(0.08, 0.85, len(order)))
    tie, all_tot = [], []
    for rank, k in enumerate(order):
        p = pct[k]
        d1, d2 = np.asarray(z[f"d1_{p}"], float), np.asarray(z[f"d2_{p}"], float)
        all_tot.append(d1 + d2)
        tie.append(dict(p=p, d1=d1, d2=d2, color=colors[rank],
                        dil=(dil1[k], dil2[k]), den=(den1[k], den2[k]), den_err=(den1e[k], den2e[k])))
    all_tot = np.concatenate(all_tot)
    lim = float(np.percentile(all_tot, 99.9)) * 1.05
    dens = np.concatenate([den1, den2]); dens = dens[np.isfinite(dens)]
    if dens.size: lim = max(lim, float(dens.max()) * 1.10)
    axins = ax.inset_axes([0.55, 0.55, 0.42, 0.42]) if add_inset else None
    for t in tie:
        ax.scatter(t["d1"], t["d2"], s=12, alpha=0.16, color=t["color"], rasterized=True, zorder=1)
        dil, den = t["dil"], t["den"]
        hd, hn = np.all(np.isfinite(dil)), np.all(np.isfinite(den))
        if hd and hn:
            ax.plot([dil[0], den[0]], [dil[1], den[1]], color=t["color"], lw=2.2, alpha=0.9, zorder=3)
        if hn:
            de = t["den_err"]
            if np.all(np.isfinite(de)) and np.any(np.asarray(de) > 0):
                ax.errorbar(den[0], den[1], xerr=de[0], yerr=de[1], fmt="none",
                            ecolor="black", elinewidth=1.1, capsize=3, capthick=1.1, zorder=4)
            ax.scatter(den[0], den[1], color=t["color"], s=80, edgecolors="black", lw=1.4, zorder=5,
                       label=f"{t['p']}:{100 - t['p']}")
        if add_inset and hd:
            axins.scatter(dil[0], dil[1], color=t["color"], s=55, edgecolors="black", lw=1.1, zorder=5)
            if hn:
                axins.plot([dil[0], den[0]], [dil[1], den[1]], color=t["color"], lw=1.8, alpha=0.9, zorder=4)
    ax.plot([0, lim], [lim, 0], color="gray", ls="--", alpha=0.6, lw=1.4, zorder=0)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    ax.set_xlabel(f"{n1}  (Res./Voxel)"); ax.set_ylabel(f"{n2}  (Res./Voxel)")
    ax.legend(title="dense-phase points   (% V2 : % V3)", fontsize=legend_fs,
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=5,
              framealpha=0.9, title_fontsize=legend_fs, columnspacing=1.0, handletextpad=0.35)
    if add_inset:
        axins.set_xlim(0, 2); axins.set_ylim(0, 2); axins.set_aspect("equal")
        axins.grid(True, alpha=0.3); axins.set_title("dilute end", fontsize=8.5)
        axins.tick_params(labelsize=6.5)

mpl.rcParams.update({"font.size": 11, "font.family": "sans-serif", "mathtext.default": "regular"})
fig = plt.figure(figsize=(15.5, 8.6))
gs = GridSpec(2, 3, figure=fig, width_ratios=[1.0, 1.0, 1.28],
              hspace=0.40, wspace=0.42, left=0.06, right=0.955, top=0.90, bottom=0.09)

axA = fig.add_subplot(gs[0, 0])
axA.scatter(f1.SCD, f1.SHD, s=52, c="0.80", edgecolors="0.5", linewidths=0.5, zorder=2)
axA.scatter([scd_v2], [shd_v2], s=52, c=CA, edgecolors="k", lw=1.0, zorder=4, label=name1 + " (A)")
axA.scatter([scd_v3], [shd_v3], s=52, c=CB, edgecolors="k", lw=1.0, zorder=5, marker="s", label=name2 + " (B)")
axA.set_xlabel("SCD"); axA.set_ylabel("SHD")
axA.set_title("A   The mixed pair on the SCD-SHD map", loc="left", fontweight="bold")
axA.legend(fontsize=8.5, loc="lower left", frameon=True); axA.grid(alpha=0.25)

axB = fig.add_subplot(gs[0, 1])
bars = axB.bar(["A-A\n(V2-V2)", "A-B\n(V2-V3)", "B-B\n(V3-V3)"], [sAA, sAB, sBB],
               color=[CA, CAB, CB], edgecolor="k", lw=0.6)
for b, v in zip(bars, [sAA, sAB, sBB]):
    axB.text(b.get_x() + b.get_width() / 2, v + 0.15, f"{v:.1f}", ha="center", fontsize=10, fontweight="bold")
axB.set_ylabel("inter-chain contacts\nper chain-pair per frame")
axB.set_title("B   Which interactions dominate", loc="left", fontweight="bold")
axB.set_ylim(0, max(sAA, sAB, sBB) * 1.20); axB.grid(alpha=0.3, axis="y")

axC = fig.add_subplot(gs[1, 0])
axC.plot(pctA, cvc["sAA"], "o-", color=CA, lw=1.9, label="A-A (V2-V2)")
axC.plot(pctA, cvc["sAB"], "s-", color=CAB, lw=1.9, label="A-B (V2-V3)")
axC.plot(pctA, cvc["sBB"], "D-", color=CB, lw=1.9, label="B-B (V3-V3)")
axC.set_xlabel("% V2 (component A)"); axC.set_ylabel("contacts per chain-pair\nper frame")
axC.set_title("C   Contact strength vs composition", loc="left", fontweight="bold")
axC.set_ylim(0, max(cvc["sBB"].max(), cvc["sAB"].max(), cvc["sAA"].max()) * 1.15)
axC.axvline(50, color="0.8", ls=":", lw=1.0); axC.legend(fontsize=8.2, loc="lower right", frameon=True)
axC.grid(alpha=0.25)

axD = fig.add_subplot(gs[1, 1])
TOTAL_DENSE = V2_DENSE + V3_DENSE
axD.plot(PCT, TOTAL_DENSE, "s--", color="0.5", lw=1.6, zorder=2, label="total")
axD.plot(PCT, V2_DENSE, "o-", color=CA, lw=1.9, zorder=3, label="V2")
axD.plot(PCT, V3_DENSE, "s-", color=CB, lw=1.9, zorder=3, label="V3")
axD.set_xlabel("% V2 (component A)"); axD.set_ylabel("dense-phase conc. (mM)")
axD.set_title("D   Phase behaviour vs composition", loc="left", fontweight="bold")
axD.axvline(50, color="0.8", ls=":", lw=1.0)
axD.legend(fontsize=8, loc="upper center", ncol=3, frameon=True); axD.grid(alpha=0.25)

axE = fig.add_subplot(gs[:, 2])
axE.set_title("E   Phase diagram (DDM tie-lines)", loc="left", fontweight="bold")
draw_phase_diagram(axE, NPZ, add_inset=True)

fig.suptitle(f"Figure 7 -- {name1} + {name2}  (298 K)", fontsize=14, fontweight="bold", y=0.975)
fig.savefig(os.path.join(HERE, "Figure7.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "Figure7.png"), dpi=300, bbox_inches="tight")
print("saved Figure7.png")
