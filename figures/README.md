# Paper figures

One directory per figure, holding its LaTeX/TikZ source, the data it plots and
the rendered panels. Each `main.tex` is a standalone wrapper that compiles that
figure on its own; `abrv.tex` and `some_command.tex` are the shared macro files
the sources expect.

Panels showing fields (density, polarity, velocity) are committed as **PNG
snapshots**. The raw `.jld` output they were rendered from is not published —
at 5008² in `Float64` a single snapshot is about 1 GB, and a full run set runs
to terabytes.

Only images a figure actually `\includegraphics` are kept. The working
directory they came from held some 300 MB of candidates and superseded
variants; what the paper compiles is a fraction of that.

**Images belong to their figure's directory, never to a shared pool.** The
names collide: `density.png` exists as three different images across the
lattice and triple-panel figures, and `voronoi.png` as two. Flattening them
would silently overwrite one with another.

## Building a figure on its own

```bash
./build_all.sh          # all of them
./build_all.sh Fig6     # just one
```

Every figure has the same layout: `main.tex` is the wrapper, and the picture
itself lives in a `fig_*.tex` beside it, with its data and images. Building a
wrapper produces `main.pdf`.

**All ten reproduce the figure PDFs that were submitted, exactly** — verified by
comparing bounding boxes against `Dumoulin_etal_FigN.pdf` in the resubmission.

The wrappers are *not* interchangeable, because the submitted figures were not
all produced the same way:

| | figures | how the wrapper is set up |
|---|---|---|
| Built **inside the article**, externalised by TikZ | 2, 4, 6, 8, 9 | must recreate the article's text block and font |
| Built as a **standalone project** and included as a PDF | 1, 3, 5, 7, 10 | keeps that project's own preamble |

For the first group the figure is written in terms of `\linewidth` and inherits
the article's font, so the wrapper reproduces both:

- **Text width set explicitly** — `\textwidth`, `\columnwidth`, `\linewidth`
  and `\hsize` to 246.0 pt for a one-column `figure` (2, 4, 8, 9) or 510.0 pt
  for a two-column `figure*` (6). The standalone class's `varwidth` option does
  not do this: with the `tikz` option the class crops to the picture and leaves
  `\linewidth` at its own default of 345.0 pt.
- **9 pt Computer Modern on a 10.5 pt baseline** — what revtex4-2 uses *inside*
  a figure, not the 10 pt document size, and not Latin Modern.
- **`abrv.tex` loaded**, since the figures use its shorthands (`\zr`, `\Dmu`,
  …). Without it those labels silently disappear.

The second group never went through the article's text block at all, so none of
that applies: each keeps the preamble it was built with, and changing it changes
the figure. Fig1 renders at `varwidth=180mm` with Latin Modern at 10 pt, and
Fig5 sets `\def\LW{246.0pt}` and sizes everything from that without ever
referring to `\linewidth`. **Do not regenerate these wrappers from the
one-column template** — the output stops matching the article.

## Where each figure comes from

One script per figure, named for it. Each takes the run set through `DATA_DIR`
and writes into `FIG_DIR`:

```bash
DATA_DIR=/path/to/fig9-runs/ FIG_DIR=$PWD/figures/Fig9/ \
    julia --project=. Analysis/make_fig9.jl
```

| Figure | Script | Established by |
|---|---|---|
| **1** | [`make_fig1.jl`](../Analysis/make_fig1.jl) | `make_plot_L50` writes `phases_L50_3.png`; rows 2, 5, 8 of the L=50 table are ρ₀ = 0.45, 0.60, 0.75, the caption's three densities |
| **2** | [`make_fig2.jl`](../Analysis/make_fig2.jl) | the two-defect routines; the figure ships a pre-built PDF rather than data |
| **3** | [`1D/models.pluto.jl`](../1D/models.pluto.jl) | its export cells; data byte-identical to [`1D/figdata/`](../1D/figdata) |
| **4** | [`make_fig4.jl`](../Analysis/make_fig4.jl) | `make_heatmap`'s six `tkd` values give 12×14×6 = 1008 rows, exactly `Ndef.csv` |
| **5** | [`make_fig5.jl`](../Analysis/make_fig5.jl) | regenerates `tau1.csv` and `tau5.csv` byte for byte |
| **6** | [`make_fig6.jl`](../Analysis/make_fig6.jl) | `make_full_heatmap_idx` / `make_zoom_heatmap_idx` |
| **7** | [`make_fig7.jl`](../Analysis/make_fig7.jl) | `df_tikz_phase_diagram` writes `DF_tikz_norm_adjusted.csv` |
| **8** | [`make_fig8.jl`](../Analysis/make_fig8.jl) | `plot_dens_and_angle` writes `density.png` and `angle.png` |
| **9** | [`make_fig9.jl`](../Analysis/make_fig9.jl) | `gamma_as_Ndefpm` writes `g34_rho-pm.csv`; see below for panel (b) |
| **10** | [`make_fig10.jl`](../Analysis/make_fig10.jl) | `make_csv_exp_plot` writes `DF_tikz_exp_<t>.csv` |
| **11** | [`make_fig11.jl`](../Analysis/make_fig11.jl) | the same routine as figure 1, over more of the L=50 series |
| **12** | [`make_fig12.jl`](../Analysis/make_fig12.jl) | `figure_triple` writes the density, Voronoi and shape-order panels |
| **13** | — | **not identified**; see below |

Two things the mapping turned up.

**Figure 9(b) comes from commented-out code.** `gamma_as_Ndefpm` writes panel
(a) live, and its final block — building `x = Ndefects` — is commented out. The
committed `g34_densdef.csv` is that block's output with `x` divided by the area
of the *padded* domain, (1008 × 0.01)² = 101.6064: the ratio between
`g34_Ndef-pm.csv` and `g34_densdef.csv` is exactly that on every row. The live
filter `Ndefects > 150` is the caption's "defect density larger than 1.5".

**Three scripts draw figures the revision removed.** The earlier version had
flow-alignment, anisotropic-stress and saturation figures; the submitted one
does not. They are parked in
[`Analysis/not_in_paper.jl`](../Analysis/not_in_paper.jl) rather than mixed in
with the figure scripts.

## Still to do

**Figure 13 has no identified generator.** Nothing in `Analysis/` obviously
draws the ten-lattice panel, and its standalone project holds only the built
PDF.

**Five data files have no generator in the repository.** `df_zrc.csv` and
`zrc_plus.csv` (Fig4's critical-activity dots, which the caption obtains
analytically), `omegarho.csv` and `zrc_eq{,2}.csv` (Fig5), and `DF_tikz.csv`
(Fig7, where the script writes `DF_tikz_norm_adjusted.csv` instead). The first
four are analytic rather than measured, so they likely come from the same kind
of calculation as [`LSA.jl`](../Analysis/LSA.jl).

**The scripts have not been run.** No simulation output is published with this
repository, so the drivers are wired from the evidence above rather than
executed. Treat the argument lists as a starting point.

**Snapshot counts are not fixed.** Runs were configured to `t_fin` but some
stop earlier (see [`params/README.md`](../params/README.md)), so anything that
indexes snapshots positionally will not read the same time point in a
reproduction as it did originally.

## Running the analysis

Each figure script takes its run set from `DATA_DIR` and writes into `FIG_DIR`:

```bash
DATA_DIR=/path/to/runset/ FIG_DIR=$PWD/figures/Fig9/ \
    julia --project=. Analysis/make_fig9.jl
```

Figures 5, 7 and 10 need no simulation output: figure 5 is computed from the
linearised equations, and the phase-diagram tables ship in
[`Analysis/PDpaper/`](../Analysis/PDpaper), so `DATA_DIR` is only needed there
to rebuild them from your own snapshots.

[`Analysis/MakePlots.jl`](../Analysis/MakePlots.jl) stands apart from the
figure scripts: it renders the standard density, angle, velocity and order
panels for any run set and reproduces no particular figure.

See the caveats in the top-level [README](../README.md#analysis).
