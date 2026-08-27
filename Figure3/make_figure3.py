#!/usr/bin/env python3
"""Figure 3 -- Strong self-associator: LAF-1 V2 (A) + protein2 (B), 298 K.
Run: python make_figure3.py"""
import os, glob
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
F1  = os.path.join(HERE, "figure1_data.csv")
NPZ = os.path.join(HERE, "contact_maps_LAF1V2__protein2.npz")
NPZ_COMP = os.path.join(HERE, "contacts_vs_composition_LAF1V2__protein2.npz")
MASTER = os.path.join(HERE, "two_comp_master.csv")

CA, CB, CAB = "#2166ac", "#018571", "#7b3294"

def find_phase_voxel():
    hits = sorted(glob.glob(os.path.join(HERE, "phase_voxel_*.npz")))
    return hits[0] if hits else None

def draw_phase_diagram(ax, npz_path, add_inset=True, legend_fs=8.0, minor="B"):
    z = np.load(npz_path, allow_pickle=True)
    pct = np.asarray(z["pct"], int)
    dil1, den1 = np.asarray(z["dil1"], float), np.asarray(z["den1"], float)
    dil2, den2 = np.asarray(z["dil2"], float), np.asarray(z["den2"], float)
    name1, name2 = str(z["name1"]), str(z["name2"])
    order  = np.argsort(pct)
    colors = plt.cm.plasma(np.linspace(0.08, 0.85, len(order)))
    tie, all_tot = [], []
    for rank, k in enumerate(order):
        p = pct[k]
        d1, d2 = np.asarray(z[f"d1_{p}"], float), np.asarray(z[f"d2_{p}"], float)
        all_tot.append(d1 + d2)
        tie.append(dict(p=p, d1=d1, d2=d2, color=colors[rank],
                        dil=(dil1[k], dil2[k]), den=(den1[k], den2[k])))
    all_tot = np.concatenate(all_tot)
    lim = float(np.percentile(all_tot, 99.9)) * 1.05
    dens = np.concatenate([den1, den2]); dens = dens[np.isfinite(dens)]
    if dens.size:
        lim = max(lim, float(dens.max()) * 1.10)
    dils = np.concatenate([dil1, dil2]); dils = dils[np.isfinite(dils)]
    lim_in = max(float(dils.max()) * 3.0, lim * 0.02) if dils.size else lim * 0.05
    axins = ax.inset_axes([0.55, 0.55, 0.42, 0.42]) if add_inset else None
    for t in tie:
        ax.scatter(t["d1"], t["d2"], s=12, alpha=0.16, color=t["color"], rasterized=True, zorder=1)
        dil, den = t["dil"], t["den"]
        hd, hn = np.all(np.isfinite(dil)), np.all(np.isfinite(den))
        if hd and hn:
            ax.plot([dil[0], den[0]], [dil[1], den[1]], color=t["color"], lw=2.2, alpha=0.9, zorder=3)
        if hn:
            ax.scatter(den[0], den[1], color=t["color"], s=80, edgecolors="black", lw=1.4, zorder=5,
                       label=f"{t['p']}:{100 - t['p']}")
        if add_inset and hd:
            axins.scatter(dil[0], dil[1], color=t["color"], s=55, edgecolors="black", lw=1.1, zorder=5)
            if hn:
                axins.plot([dil[0], den[0]], [dil[1], den[1]], color=t["color"], lw=1.8, alpha=0.9, zorder=4)
    ax.plot([0, lim], [lim, 0], color="gray", ls="--", alpha=0.6, lw=1.4, zorder=0)
    for c in [t for t in ax.get_xticks() if 0 <= t <= lim]:
        ax.plot([0, c], [c, 0], color="gray", ls="--", alpha=0.22, lw=0.7, zorder=0)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    ax.set_xlabel(f"{name1}  (Res./Voxel)"); ax.set_ylabel(f"{name2}  (Res./Voxel)")
    ax.legend(title=f"dense-phase points   (% V2 : % {minor})", fontsize=legend_fs,
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=5,
              framealpha=0.9, title_fontsize=legend_fs, columnspacing=1.0, handletextpad=0.35)
    if add_inset:
        axins.set_xlim(0, 0.5); axins.set_ylim(0, 0.5); axins.set_aspect("equal")
        axins.grid(True, alpha=0.3); axins.set_title("dilute end", fontsize=8.5)
        axins.tick_params(labelsize=6.5)

def main():
    d = np.load(NPZ, allow_pickle=True)
    sAA, sAB, sBB = float(d["sAA"]), float(d["sAB"]), float(d["sBB"])
    name1, name2 = str(d["name1"]), str(d["name2"])

    dc = np.load(NPZ_COMP, allow_pickle=True)
    pct = dc["pctA"]

    m = pd.read_csv(MASTER); m = m[m.partner == "protein2"].sort_values("pctV2")
    PCT = m.pctV2.values
    V2_DENSE, V2_DENSE_SEM = m.dense_V2.values, m.dense_V2_sem.values
    B_DENSE, B_DENSE_SEM = m.dense_B.values, m.dense_B_sem.values

    f1 = pd.read_csv(F1)
    scd_v2 = f1.loc[f1.key == "LAF-1_V2", "SCD"].values[0]
    shd_v2 = f1.loc[f1.key == "LAF-1_V2", "SHD"].values[0]
    scd_b = f1.loc[f1.key == "protein2", "SCD"].values[0]
    shd_b = f1.loc[f1.key == "protein2", "SHD"].values[0]

    mpl.rcParams.update({
        "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "font.family": "sans-serif", "mathtext.default": "regular",
    })
    fig = plt.figure(figsize=(15.5, 8.6))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.0, 1.0, 1.28],
                  hspace=0.40, wspace=0.42, left=0.06, right=0.955, top=0.90, bottom=0.09)

    axA = fig.add_subplot(gs[0, 0])
    axA.scatter(f1.SCD, f1.SHD, s=52, c="0.80", edgecolors="0.5", linewidths=0.5, zorder=2)
    axA.scatter([scd_v2], [shd_v2], s=175, c=CA, edgecolors="k", lw=1.0, zorder=4, label=f"{name1} (A)")
    axA.scatter([scd_b], [shd_b], s=105, c=CB, edgecolors="k", lw=1.0, zorder=5, marker="D", label=f"{name2} (B)")
    axA.set_xlabel("SCD"); axA.set_ylabel("SHD")
    axA.set_title("A   The mixed pair (strong self-associator)", loc="left", fontweight="bold")
    axA.legend(fontsize=8.5, loc="lower left", frameon=True); axA.grid(alpha=0.25)

    axB = fig.add_subplot(gs[0, 1])
    bars = axB.bar(["A-A\n(V2-V2)", "A-B\n(V2-P2)", "B-B\n(P2-P2)"], [sAA, sAB, sBB],
                   color=[CA, CAB, CB], edgecolor="k", lw=0.6)
    for b, v in zip(bars, [sAA, sAB, sBB]):
        axB.text(b.get_x() + b.get_width() / 2, v + 0.12, f"{v:.1f}", ha="center", fontsize=10, fontweight="bold")
    axB.set_ylabel("inter-chain contacts\nper chain-pair per frame")
    axB.set_title("B   Which interactions dominate", loc="left", fontweight="bold")
    axB.set_ylim(0, max(sAA, sAB, sBB) * 1.20); axB.grid(alpha=0.3, axis="y")

    axC = fig.add_subplot(gs[1, 0])
    axC.plot(pct, dc["sAA"], "o-", color=CA, lw=1.9, label="A-A (V2-V2)")
    axC.plot(pct, dc["sAB"], "s-", color=CAB, lw=1.9, label="A-B (V2-P2)")
    axC.plot(pct, dc["sBB"], "D-", color=CB, lw=1.9, label="B-B (P2-P2)")
    axC.set_xlabel("% V2 (component A)"); axC.set_ylabel("contacts per chain-pair\nper frame")
    axC.set_title("C   Contact strength vs composition", loc="left", fontweight="bold")
    axC.set_ylim(0, max(dc["sBB"].max(), dc["sAB"].max(), dc["sAA"].max()) * 1.15)
    axC.axvline(50, color="0.8", ls=":", lw=1.0, zorder=0)
    axC.legend(fontsize=8.2, loc="lower right", frameon=True); axC.grid(alpha=0.25)

    axD = fig.add_subplot(gs[1, 1])
    axD.errorbar(PCT, V2_DENSE, yerr=V2_DENSE_SEM, fmt="o-", color=CA, lw=1.9, capsize=3, label="V2 dense")
    axD.errorbar(PCT, B_DENSE, yerr=B_DENSE_SEM, fmt="D-", color=CB, lw=1.9, capsize=3, label="protein2 dense")
    axD.set_xlabel("% V2 (component A)"); axD.set_ylabel("dense-phase conc. (mM)")
    axD.set_title("D   Phase behaviour vs composition", loc="left", fontweight="bold")
    axD.axvline(50, color="0.8", ls=":", lw=1.0, zorder=0)
    axD.legend(fontsize=8.5, loc="upper center", frameon=True); axD.grid(alpha=0.25)

    axE = fig.add_subplot(gs[:, 2])
    axE.set_title("E   Phase diagram (DDM tie-lines)", loc="left", fontweight="bold")
    pv = find_phase_voxel()
    if pv:
        draw_phase_diagram(axE, pv, add_inset=True, minor="P2")
    else:
        axE.text(0.5, 0.5, "phase_voxel_*.npz not found", ha="center", va="center", fontsize=10, color="0.4",
                 transform=axE.transAxes)
        axE.set_xticks([]); axE.set_yticks([])

    fig.suptitle(f"Figure 3 -- A strong self-associator mixes with the reference:  {name1} + {name2}  (298 K)",
                 fontsize=13.5, fontweight="bold", y=0.975)
    fig.savefig(os.path.join(HERE, "Figure3.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "Figure3.png"), dpi=300, bbox_inches="tight")
    print("saved Figure3.png")

if __name__ == "__main__":
    main()
