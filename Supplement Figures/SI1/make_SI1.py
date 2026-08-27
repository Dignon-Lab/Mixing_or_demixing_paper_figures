#!/usr/bin/env python3
"""SI Figure 1 -- single-component critical temperature Tc vs charge patterning.
Run: python make_SI1.py"""
import os
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# label, SCD, SHD, SAD, dG, rhoc, Tc(K)
data = [
    ('P1',  -20.9683, 5.3534, 0.02178, -8.4560, 179.61, 339.08),
    ('P2',  -19.7489, 5.3851, 0.06990, -8.4805, 187.80, 328.67),
    ('P3',  -19.5157, 5.3140, 0.01690, -8.4250, 184.47, 331.64),
    ('P4',  -14.9950, 5.2777, 0.05917, -9.1627, 184.43, 339.85),
    ('P5',  -15.0000, 5.3300, 0.03210, -8.9231, 176.27, 323.79),
    ('P6',  -15.0000, 5.3816, 0.03940, -8.9920, 177.14, 312.42),
    ('P7',   -9.9812, 5.3706, 0.04850, -8.6177, 166.73, 328.67),
    ('P8',   -9.8382, 5.2798, 0.04280, -8.7060, 170.41, 317.78),
    ('P9',  -10.0008, 5.3300, 0.02420, -8.7439, 168.68, 312.49),
    ('P10',  -5.0022, 5.3300, 0.02340, -8.4160, 165.43, 311.84),
    ('P11',  -4.9891, 5.3857, 0.04430, -8.2059, 170.91, 302.85),
    ('P12',  -5.0320, 5.2579, 0.06960, -8.3050, 180.99, 312.05),
    ('P13',   0.0423, 5.2672, 0.04880, -4.0120, 155.97, 277.46),
    ('P14',   0.0000, 5.3300, 0.02655, -4.4784, 162.22, 274.37),
    ('P15',  -0.2352, 5.3954, 0.06620, -5.2371, 168.23, 283.52),
    ('WT',    0.0200, 5.3390, 0.01900, np.nan,  np.nan, 282.60),
    ('V1',   -5.0000, 5.3270, 0.03100, np.nan,  np.nan, 302.60),
    ('V2',   -9.9970, 5.3280, 0.02000, np.nan,  np.nan, 307.70),
    ('V3',  -15.0020, 5.3490, 0.02800, np.nan,  np.nan, 319.90),
]
df = pd.DataFrame(data, columns=['label', 'SCD', 'SHD', 'SAD', 'dG', 'rhoc', 'Tc'])
df = df.sort_values('SCD').reset_index(drop=True)

fig, ax = plt.subplots(figsize=(9.0, 4.6))
norm = mpl.colors.Normalize(vmin=df.SCD.min(), vmax=df.SCD.max())
cmap = plt.cm.viridis
x = np.arange(len(df))

sc = ax.scatter(x, df.Tc, c=df.SCD, cmap=cmap, norm=norm, s=110,
                edgecolors='k', linewidths=0.5, zorder=3)
for xi, row in zip(x, df.itertuples()):
    ax.annotate(row.label, (xi, row.Tc), fontsize=8, fontweight='bold',
                ha='center', va='bottom', xytext=(0, 8), textcoords='offset points')
ax.margins(y=0.10)
ax.set_xticks(x); ax.set_xticklabels(df.label, rotation=90)
ax.set_xlabel('protein (ordered by SCD: blocky -> charge-mixed)')
ax.set_ylabel(r'critical temperature  $T_c$  (K)')
ax.set_title('SI Figure 1 -- single-component $T_c$ vs charge patterning', loc='left', fontweight='bold')
ax.grid(alpha=0.25, axis='y'); ax.margins(x=0.02)
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label('charge patterning  SCD')

fig.tight_layout()
fig.savefig(os.path.join(HERE, 'FigureS1.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(HERE, 'FigureS1.png'), dpi=300, bbox_inches='tight')
print('saved FigureS1.png')
