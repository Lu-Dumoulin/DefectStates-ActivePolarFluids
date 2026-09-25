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

The panels were drawn from simulation output that is not published here, so
this is a record of which routine produced which figure rather than a pipeline
you can run. The routines are in
[`Analysis/figure_routines.jl`](../Analysis/figure_routines.jl); the
measurement primitives they build on are in
[`Analysis/core.jl`](../Analysis/core.jl).

| Figure | Routine | Established by |
|---|---|---|
| **1** | `make_plot_L50` | writes `phases_L50_3.png`; rows 2, 5, 8 of the L=50 table are ρ₀ = 0.45, 0.60, 0.75, the caption's three densities |
| **2** | `dist_2defects`, `plot_2defects_zoom` | the two-defect routines; the figure ships a pre-built PDF rather than data |
| **3** | [`1D/models.pluto.jl`](../1D/models.pluto.jl) | **reproducible from this repository**; its data is byte-identical to [`1D/figdata/`](../1D/figdata) |
| **4** | `make_heatmap` | its six `tkd` values give 12×14×6 = 1008 rows, exactly `Ndef.csv` |
| **5** | [`Analysis/linear_stability.jl`](../Analysis/linear_stability.jl) | **reproducible**; regenerates `tau1.csv` and `tau5.csv` byte for byte, with no simulation output |
| **6** | `make_full_heatmap_idx`, `make_zoom_heatmap_idx` | the snapshot grid and the zoomed panels |
| **7** | `df_tikz_phase_diagram` in [`PhaseDiagram.jl`](../Analysis/PhaseDiagram.jl) | writes `DF_tikz_norm_adjusted.csv` |
| **8** | `plot_dens_and_angle` in [`MakePlots.jl`](../Analysis/MakePlots.jl) | writes `density.png` and `angle.png` |
| **9** | `gamma_as_Ndefpm` | writes `g34_rho-pm.csv`; see below for panel (b) |
| **10** | `make_csv_exp_plot` in [`PhaseDiagram.jl`](../Analysis/PhaseDiagram.jl) | writes `DF_tikz_exp_<t>.csv` |
| **11** | `make_plot_L50`, over more of the L=50 series | the caption notes ρ₀ = 0.45 and 0.55 appear here but not in figure 6 |
| **12** | `figure_triple` | the density, Voronoi and shape-order panels |
| **13** | — | **not identified** |

Two things the mapping turned up.

**Figure 9(b) comes from commented-out code.** `gamma_as_Ndefpm` writes panel
(a) live, and its final block — building `x = Ndefects` — is commented out. The
committed `g34_densdef.csv` is that block's output with `x` divided by the area
of the *padded* domain, (1008 × 0.01)² = 101.6064: the ratio between
`g34_Ndef-pm.csv` and `g34_densdef.csv` is exactly that on every row. The live
filter `Ndefects > 150` is the caption's "defect density larger than 1.5".

**Three routines drew figures the revision removed.** The earlier version had
flow-alignment, anisotropic-stress and saturation figures; the submitted one
does not. They are in
[`Analysis/not_in_paper.jl`](../Analysis/not_in_paper.jl).

## Still to do

**Figure 13 has no identified generator.** Nothing in `Analysis/` obviously
draws the ten-lattice panel, and its standalone project holds only the built
PDF.

**The critical-activity curves come from the linear stability analysis** —
Fig4's white dots (`zrc_plus.csv`, and the superseded `df_zrc.csv`) and Fig5's
`zrc_eq{,2}.csv` are the same quantity, ζ_ρ^c, solved self-consistently with
a = 4ζ_ρ^c/3. The routine that writes them is not in this repository, and it
could not be reconstructed: reading `zrc_plus.csv` in the panel coordinates its
`fig_ndef.tex` uses, no arrangement of the equation reproduces the committed
numbers — the closest of eight variants is ~7 % out.

**That search turned up a discrepancy worth checking.** The article's
Eq. (zetaRhoDeltaMuC) reads

    zeta_c = A(chi=0)/(3 rho0) + ( sqrt(A(chi=0) gamma / (3 rho0^3)) + sqrt(2/(3 rho0^3 tau)) )^2

while `ζc_simp` in [`LSA.jl`](../Analysis/LSA.jl), which reproduces Fig5's
`tau1.csv` and `tau5.csv` byte for byte, differs from it in three places: it
keeps the χ term in `A` instead of setting χ = 0, divides by `3 rho0^2` rather
than `3 rho0^3`, and **subtracts** the second root rather than adding it. One
of the two is not what the other describes.

**`omegarho.csv` and `DF_tikz.csv` also have no generator here** — the latter
because `make_fig7.jl` writes `DF_tikz_norm_adjusted.csv` instead.

**The scripts have not been run.** No simulation output is published with this
repository, so the drivers are wired from the evidence above rather than
executed. Treat the argument lists as a starting point.

**Snapshot counts are not fixed.** Runs were configured to `t_fin` but some
stop earlier (see [`params/README.md`](../params/README.md)), so anything that
indexes snapshots positionally will not read the same time point in a
reproduction as it did originally.

## What can actually be run

Two things, both without any simulation output:

```bash
julia --project=. Analysis/linear_stability.jl     # figure 5's data
```

and the figure 3 notebook in [`1D/`](../1D). Everything else in `Analysis/` is
reference material — see the note at the top of
[`figure_routines.jl`](../Analysis/figure_routines.jl).

The simulations themselves are reproducible: [`2D/`](../2D) with the parameter
tables in [`params/`](../params).
