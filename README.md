# Designed IDPs phase separate and mix or demix according to sequence design parameters

Arjun Singh and Gregory L. Dignon, Rutgers University.

This repository contains the data and plotting scripts needed to reproduce every
figure in the main text and Supporting Information of the paper. Each figure has
its own folder with one Python script and the data it reads; running the script
regenerates that figure's PDF/PNG from scratch.

## Layout

```
Figure1/ ... Figure7/   one folder per main-text figure
SI/SI1/  ... SI/SI13/   one folder per supporting-information figure
```

Each folder contains `make_figureN.py` (or `make_SIN.py`), the small CSV/NPZ
data files it reads, and produces `FigureN.pdf`/`.png` (or `FigureSN.pdf`/`.png`)
in the same folder when run.

## Requirements

Python 3.9+ with:

```
pip install -r requirements.txt
```

(numpy, pandas, matplotlib, scipy)

## Running a figure

```
cd Figure2
python make_figure2.py
```

Each script is self-contained and only reads files in its own folder.

## Data sources

- `figure1_data.csv` -- per-sequence SCD/SHD/Rg/nu/dense/dilute/diffusivity table
  for the 19-member LAF-1 RGG library (single-chain and single-component slabs).
- `two_comp_master.csv` -- per-composition Rg/nu/dense/dilute/tauRg/diffusivity
  for LAF-1 V2 mixed with each of protein2, protein9, protein13, protein15,
  parsed from the make-whole-corrected two-component analysis run.
- `contact_maps_*.npz` / `contacts_vs_composition_*.npz` -- residue-residue
  intermolecular contact maps and their composition dependence, from a 1.0 nm
  C-alpha cutoff contact analysis.
- `phase_voxel_*.npz` -- per-voxel domain-decomposition (DDM) exports (35 A
  voxels for the individual pairs, Figures 2-5 and 7; the same method at 20 A
  underlies the 18-partner screen in Figure 6) used for the phase-diagram
  panels and the demixing index D.
- `analysis_37377.out` -- the full make-whole-corrected two-component analysis
  report (all 18 LAF-1 V2 + partner datasets, 5 compositions each), used
  directly by the Figure 6 / SI Figure 11 screen scripts.

## Note on density-profile panels

A few SI figures (SI4, SI6, SI8, SI10) pair a demixing-index panel with a
slab density-profile panel in the published figure. The demixing-index panel
is fully reproduced here from the DDM voxel export. The density-profile panel
needs the raw production trajectory, which is not distributed in this
repository (see the paper's Data and code availability statement); those
scripts reproduce the demixing-index panel only and say so in a comment.

## Citation

If you use this code or data, please cite the paper.
