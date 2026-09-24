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

| Figure | Article label | Drawn by |
|---|---|---|
| **Fig1** | `fig:schema` | [`make_fig_schema.jl`](../Analysis/make_fig_schema.jl) — `make_plot_L50()` renders `phases_L50_3.png` |
| **Fig2** | `fig:stab_defects` | [`make_fig_2defects.jl`](../Analysis/make_fig_2defects.jl) |
| **Fig3** | `fig:singleDefect` | [`1D/models.pluto.jl`](../1D/models.pluto.jl) — **fully reproducible here**; the data is byte-identical to [`1D/figdata/`](../1D/figdata) |
| **Fig4** | `fig:Ndefects` | [`make_fig_ndef.jl`](../Analysis/make_fig_ndef.jl) |
| **Fig5** | `fig:linstab` | [`make_fig_LSA.jl`](../Analysis/make_fig_LSA.jl) — `tau1.csv` and `tau5.csv` **regenerate byte for byte** |
| **Fig6** | `fig:states` | [`make_fig_phases.jl`](../Analysis/make_fig_phases.jl) |
| **Fig7** | `fig:phasediagram` | [`PhaseDiagram.jl`](../Analysis/PhaseDiagram.jl) — `df_tikz_phase_diagram()`; tables ship in [`Analysis/PDpaper/`](../Analysis/PDpaper) |
| **Fig8** | `fig:lattice1` | [`make_fig_lattice.jl`](../Analysis/make_fig_lattice.jl) |
| **Fig9** | `fig:gamma_rho` | [`make_fig_gamma_rho.jl`](../Analysis/make_fig_gamma_rho.jl) |
| **Fig10** | `fig:tauc` | [`PhaseDiagram.jl`](../Analysis/PhaseDiagram.jl) — `make_csv_exp_plot(t)` |

Most of these need simulation output that is not published here; see the
caveats in the top-level [README](../README.md#analysis).

## Still to do

**Two of Fig5's data files are unaccounted for.** `make_fig_LSA.jl` regenerates
`tau1.csv` and `tau5.csv` byte for byte, but nothing here produces
`omegarho.csv` (six growth-rate curves against ρ₀ over 0.3–0.8) or
`zrc_eq.csv` / `zrc_eq2.csv` (critical-activity curves in `yp`/`zp`/`wp` and
`ynop`/`znop`/`wnop` variants, apparently with and without polarity).
[`LSA.jl`](../Analysis/LSA.jl) defines the 2×2 coefficients they would come
from but computes no eigenvalues, produces no variant curves, and writes no
files.

**Figures 11, 12 and 13 are not here** — the phase panels at L=50, the extended
lattice figure (submitted as two files, `Fig12ac` and `Fig12df`) and the
ten-lattice figure. Their standalone projects exist alongside the others.

## Running the analysis

Each script takes its run set from the environment:

```bash
DATA_DIR=/path/to/runset/ FIG_DIR=/path/to/figures/ julia --project=. Analysis/make_fig_schema.jl
```

The per-figure scripts available so far, and the figures they belong to:

| Script | Draws |
|---|---|
| `make_fig_schema.jl` | the schematic / L50 phase panels — **Fig1** |
| `make_fig_phases.jl` | the phase overview heatmaps |
| `make_fig_lattice.jl` | the square and hexagonal lattice triples |
| `make_fig_ndef.jl` | defect counts vs density and activity |
| `make_fig_gamma_rho.jl` | the shape order Γ against ρ₀ |
| `make_fig_nu.jl` | the flow-alignment (ν₁) panels |
| `make_fig_aniso.jl` | the anisotropic-stress panels |
| `make_fig_2defects.jl` | the two-defect figures |
| `make_fig_saturation.jl` | the saturation figures |
| `PhaseDiagram.jl` | the phase diagram and exponent fits — **Fig7**, **Fig10** |
| `make_fig_LSA.jl` | the linear-stability curves — **Fig5** |

These are named by content, not by figure number, because the revised
numbering is not settled; only Fig1, Fig3, Fig5, Fig7 and Fig10 are known. They
will be renamed `make_figN.jl` once it is, to line up with the per-figure
`DF_N.csv` parameter tables.

The phase-diagram tables ship in [`Analysis/PDpaper/`](../Analysis/PDpaper),
so `PhaseDiagram.jl` finds them without any setup; `DATA_DIR` overrides that
when you want to rebuild them from your own snapshots.

`Analysis/MakePlots.jl` is separate from all of these: it renders the standard
panels for any run set and reproduces no particular figure. See the caveats in
the top-level [README](../README.md#analysis).
