#!/usr/bin/env python3
"""Figure 6 -- Two-component mixing landscape: LAF-1 V2 + all 18 partners, 298 K.
Run: python make_figure6.py"""
import os
import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import linregress
from screen_io import load

HERE = os.path.dirname(os.path.abspath(__file__))
COMPS = [5, 25, 50, 75, 95]
NON_PS = ['WT', 'P13', 'P14', 'P15']   # do not phase-separate on their own at 298 K
DEMIXERS = ['P14', 'WT', 'P13', 'P15']
TC = {'P1': 339.08, 'P2': 328.67, 'P3': 331.64, 'P4': 339.85, 'P5': 323.79,
      'P6': 312.42, 'P7': 328.67, 'P8': 317.78, 'P9': 312.49, 'P10': 311.84,
      'P11': 302.85, 'P12': 312.05, 'P13': 277.46, 'P14': 274.37, 'P15': 283.52,
      'WT': 282.60, 'V1': 302.60, 'V3': 319.90}

M = load()
M['Tc'] = M['label'].map(TC)
partners = sorted(M.label.unique(), key=lambda k: M[M.label == k].SCD.iloc[0])
scd = {k: M[M.label == k].SCD.iloc[0] for k in partners}
PS = [p for p in partners if p not in NON_PS]
scd_norm = mpl.colors.Normalize(vmin=min(scd.values()), vmax=max(scd.values()))
scd_cmap = plt.cm.viridis
def pcol(k): return scd_cmap(scd_norm(scd[k]))

def add_fit(ax, x, y, color='k', ls='--', lw=1.7):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2 or np.ptp(x) == 0:
        return np.nan
    lr = linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, lr.intercept + lr.slope * xs, ls, color=color, lw=lw, zorder=2)
    return lr.rvalue

mpl.rcParams.update({"font.size": 11.5, "font.family": "sans-serif", "mathtext.default": "regular"})
fig = plt.figure(figsize=(15, 6.2))
gs = GridSpec(1, 2, figure=fig, wspace=0.30, left=0.06, right=0.93, top=0.88, bottom=0.12)

axA = fig.add_subplot(gs[0, 0])
for k in partners:
    s = M[M.label == k].sort_values('pct')
    axA.plot(s.pct, s.D, '-o', ms=4.5, lw=1.5, color=pcol(k), alpha=0.85)
for k in DEMIXERS:
    s = M[(M.label == k) & (M.pct == 50)]
    axA.annotate(k, (50, s.D.values[0]), fontsize=8.5, fontweight='bold', xytext=(4, 2), textcoords='offset points')
axA.set_xlabel('% V2 (component A)'); axA.set_ylabel('demixing index  D')
axA.set_title('A   Demixing vs composition (18 partners)', loc='left', fontweight='bold')
axA.set_xticks(COMPS); axA.grid(alpha=0.25)

axB = fig.add_subplot(gs[0, 1])
s50 = M[M.pct == 50].copy()
specs = [('Super-critical (Tc<298K, no self-LLPS)', NON_PS, 's', '#c1272d'),
         ('Sub-critical (Tc>298K, self-LLPS)', PS, 'o', '#1b1f5c')]
for name, grp, mk, fc in specs:
    d = s50[s50.label.isin(grp)]
    axB.scatter(d.Tc, d.D, s=85, marker=mk, c=d.SCD, cmap=scd_cmap, norm=scd_norm,
                edgecolors='k', linewidths=0.5, zorder=3)
    add_fit(axB, d.Tc, d.D, color=fc)
for row in s50.itertuples():
    axB.annotate(row.label, (row.Tc, row.D), fontsize=7, fontweight='bold', ha='center', va='bottom',
                 xytext=(0, 6), textcoords='offset points')
axB.axvline(298.0, color='0.5', ls=':', lw=1.2)
axB.set_xlim(270, 345); axB.set_ylim(-0.03, 0.80)
axB.set_xlabel(r'Single-component critical temperature  $T_c$  (K)')
axB.set_ylabel('Demixing index  D  (50:50)')
axB.set_title('B   D at 50:50 vs partner $T_c$', loc='left', fontweight='bold')
axB.grid(alpha=0.25)

sm = mpl.cm.ScalarMappable(norm=scd_norm, cmap=scd_cmap); sm.set_array([])
cax = fig.add_axes([0.945, 0.12, 0.012, 0.76])
cb = fig.colorbar(sm, cax=cax); cb.set_label('partner charge patterning  SCD', fontsize=10)

fig.suptitle('Figure 6 -- Two-component mixing landscape: 18 partners + LAF-1 V2 (298 K)',
             fontsize=14.5, fontweight='bold', y=0.975)
fig.savefig(os.path.join(HERE, 'Figure6.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(HERE, 'Figure6.png'), dpi=250, bbox_inches='tight')
print('saved Figure6.png')
