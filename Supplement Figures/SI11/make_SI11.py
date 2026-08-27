#!/usr/bin/env python3
"""SI Figure 11 -- Tc vs D and D vs SCD at 50:50, split by whether the partner
phase-separates on its own at 298 K (screen of 18 partners + LAF-1 V2).
Run: python make_SI11.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import linregress
from screen_io import load

HERE = os.path.dirname(os.path.abspath(__file__))
NON_PS = ['WT', 'P13', 'P14', 'P15']
TC = {'P1': 339.08, 'P2': 328.67, 'P3': 331.64, 'P4': 339.85, 'P5': 323.79,
      'P6': 312.42, 'P7': 328.67, 'P8': 317.78, 'P9': 312.49, 'P10': 311.84,
      'P11': 302.85, 'P12': 312.05, 'P13': 277.46, 'P14': 274.37, 'P15': 283.52,
      'WT': 282.60, 'V1': 302.60, 'V3': 319.90}
D_YLIM = (-0.03, 0.80); SCD_XLIM = (-22, 1.5)

M = load()
M['Tc'] = M['label'].map(TC)
partners = sorted(M.label.unique(), key=lambda k: M[M.label == k].SCD.iloc[0])
scd = {k: M[M.label == k].SCD.iloc[0] for k in partners}
PS = [p for p in partners if p not in NON_PS]
scd_norm = mpl.colors.Normalize(vmin=min(scd.values()), vmax=max(scd.values()))
scd_cmap = plt.cm.viridis

def add_fit(ax, x, y, color='k'):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2 or np.ptp(x) == 0:
        return np.nan
    lr = linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, lr.intercept + lr.slope * xs, '--', color=color, lw=1.4, zorder=2)
    return lr.rvalue

def scatter(ax, s, xcol):
    sc = ax.scatter(s[xcol], s.D, s=80, c=s.SCD, cmap=scd_cmap, norm=scd_norm,
                     edgecolors='k', linewidths=0.5, zorder=3)
    for row in s.itertuples():
        ax.annotate(row.label, (getattr(row, xcol), row.D), fontsize=8, fontweight='bold',
                    ha='center', va='bottom', xytext=(0, 7), textcoords='offset points')
    return sc

mpl.rcParams.update({"font.size": 11, "font.family": "sans-serif", "mathtext.default": "regular"})
fig = plt.figure(figsize=(11.5, 9.5))
gs = GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.34, left=0.08, right=0.92, top=0.93, bottom=0.08)
s50 = M[M.pct == 50].copy()

axA = fig.add_subplot(gs[0, 0])
s = s50[s50.label.isin(NON_PS)]
scatter(axA, s, 'Tc')
r = add_fit(axA, s.Tc, s.D)
axA.text(0.97, 0.95, f"r = {r:.2f}", transform=axA.transAxes, va="top", ha="right", fontsize=11)
axA.set_xlim(270, 290); axA.set_ylim(*D_YLIM)
axA.set_xlabel("Single-component $T_c$ (K)"); axA.set_ylabel("Demixing index  D  (50:50)")
axA.set_title("A   Tc vs D -- do NOT self-phase-separate (4)", loc="left", fontweight="bold"); axA.grid(alpha=0.25)

axB = fig.add_subplot(gs[0, 1])
s = s50[s50.label.isin(PS)]
scatter(axB, s, 'Tc')
r = add_fit(axB, s.Tc, s.D)
axB.text(0.97, 0.95, f"r = {r:.2f}", transform=axB.transAxes, va="top", ha="right", fontsize=11)
axB.set_xlim(300, 350); axB.set_ylim(0, 0.2)
axB.set_xlabel("Single-component $T_c$ (K)"); axB.set_ylabel("Demixing index  D  (50:50)")
axB.set_title("B   Tc vs D -- DO self-phase-separate (14)", loc="left", fontweight="bold"); axB.grid(alpha=0.25)

axC = fig.add_subplot(gs[1, 0])
s = s50[s50.label.isin(NON_PS)]
scatter(axC, s, 'SCD')
axC.set_xlim(*SCD_XLIM); axC.set_ylim(*D_YLIM)
axC.set_xlabel("partner charge patterning  SCD"); axC.set_ylabel("demixing index  D  (50:50)")
axC.set_title("C   D vs SCD -- do NOT self-phase-separate (4)", loc="left", fontweight="bold"); axC.grid(alpha=0.25)

axD = fig.add_subplot(gs[1, 1])
s = s50[s50.label.isin(PS)]
scatter(axD, s, 'SCD')
axD.set_xlim(*SCD_XLIM); axD.set_ylim(*D_YLIM)
axD.set_xlabel("partner charge patterning  SCD"); axD.set_ylabel("demixing index  D  (50:50)")
axD.set_title("D   D vs SCD -- DO self-phase-separate (14)", loc="left", fontweight="bold"); axD.grid(alpha=0.25)

sm = mpl.cm.ScalarMappable(norm=scd_norm, cmap=scd_cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=[axA, axB, axC, axD], pad=0.02, aspect=40, shrink=0.6)
cb.set_label("partner charge patterning  SCD")

fig.suptitle("SI Figure 11 -- Screen supporting panels: Tc and SCD vs demixing index", fontsize=13.5, fontweight="bold", y=0.97)
fig.savefig(os.path.join(HERE, "FigureS11.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "FigureS11.png"), dpi=300, bbox_inches="tight")
print("saved FigureS11.png")
