# Figures

One directory per figure of the article, containing its LaTeX source
(`fig_*.tex`), the data it plots, and the images it includes. Field snapshots
(density, polarity, velocity) are included as PNG images; the raw simulation
output they were made from is not distributed, as it amounts to several
terabytes.

## Compiling

Each directory has a `main.tex` that compiles the figure on its own, at the
size it has in the article:

```bash
./build_all.sh          # all figures
./build_all.sh Fig6     # one figure
```

This produces `FigN/main.pdf`. It needs `pdflatex` with the usual `pgfplots`
and TikZ packages.

## Contents

| Figure | Shows | Data from |
|---|---|---|
| 1 | schematic of the model; density snapshots at $\rho_0 = 0.45, 0.6, 0.75$, $L = 50$ | [`params/DF_1.csv`](../params/DF_1.csv) |
| 2 | defect pairs: number of defects over time, density and velocity around a pair, critical annihilation distance | [`params/DF_2.csv`](../params/DF_2.csv) for (a); two-defect simulations for (b–d), see [`params/`](../params) |
| 3 | steady state of a single defect | the 1D notebook, [`1D/`](../1D) |
| 4 | defect density in the $(\rho_0, \zeta_\rho)$ plane for six renewal times | [`params/DF_4.csv`](../params/DF_4.csv); white dots from the linear stability analysis |
| 5 | linear stability of the homogeneous states | linear stability analysis; data in `Fig5/` |
| 6 | asymptotic states for different target densities, $L = 10$ | [`params/DF_6.csv`](../params/DF_6.csv) |
| 7 | classification of the states in the $(\tau, \rho_0)$ plane | [`params/DF_7.csv`](../params/DF_7.csv) |
| 8 | square and hexagonal defect lattices | [`params/DF_8.csv`](../params/DF_8.csv) |
| 9 | shape order against target density and defect density | [`params/DF_9.csv`](../params/DF_9.csv) |
| 10 | density correlation function and its exponential fits | [`params/DF_10.csv`](../params/DF_10.csv) |
| 11 | density snapshots for $L = 50$ | [`params/DF_11.csv`](../params/DF_11.csv) |
| 12 | the lattices of figure 8: Voronoi tessellation and shape order over time | [`params/DF_12.csv`](../params/DF_12.csv) |
| 13 | the lattices of figure 8 from five different initial conditions | [`params/DF_13.csv`](../params/DF_13.csv) |

The quantities shown — defect density, structural persistence time, areas of
low-density regions, shape order — are defined in the appendices of the
article. Defect densities are counted per unit area of the computational
domain, $(1008 \Delta x)^2 = 101.6$ for $L = 10$.

The data of figure 3 can be regenerated with the notebook in [`1D/`](../1D).
The linear stability analysis in [`LSA/`](../LSA) is described in the main
[README](../README.md#linear-stability-analysis).
