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

## What each figure shows, and where its data came from

The analysis that produced these panels ran against simulation output that is
not published here — a single snapshot is about a gigabyte and a run set
reaches terabytes — so this is a description of how each figure was obtained,
not a pipeline. The methods are described in the article's appendices.

| Figure | Shows | Obtained from |
|---|---|---|
| **1** | schematic, and density snapshots at ρ₀ = 0.45, 0.6, 0.75 | three L=50 runs, [`params/DF_1.csv`](../params/DF_1.csv) |
| **2** | defect pair: count against time, density and velocity maps, critical annihilation distance | two-defect runs |
| **3** | steady state of a single defect | **reproducible**: [`1D/models.pluto.jl`](../1D/models.pluto.jl); its data is byte-identical to [`1D/figdata/`](../1D/figdata) |
| **4** | defect density over the (ρ₀, ζ_ρ) plane at six renewal times | [`params/DF_4.csv`](../params/DF_4.csv), 1008 runs; the white dots are the critical activity from the stability analysis |
| **5** | linear stability of the homogeneous state | **reproducible**: [`LSA/linear_stability.jl`](../LSA/linear_stability.jl) regenerates `tau1.csv` and `tau5.csv` byte for byte, with no simulation output |
| **6** | asymptotic states at different target densities | [`params/DF_6.csv`](../params/DF_6.csv), eight L=10 runs at τ=5, ζ_ρ=4 |
| **7** | state classification over τ and ρ₀, with defect density, persistence time and low-density areas | the phase-diagram runs |
| **8** | square and hexagonal defect lattices | [`params/DF_8.csv`](../params/DF_8.csv), two runs to t = 2·10⁶ |
| **9** | shape function against target density and against defect density | [`params/DF_9.csv`](../params/DF_9.csv), 1680 runs |
| **10** | density correlation against time, with exponential fits | the phase-diagram runs |
| **11** | density snapshots at L=50 across the density range | [`params/DF_11.csv`](../params/DF_11.csv) |
| **12** | the lattices of figure 8 with Voronoi tessellation and shape order | [`params/DF_12.csv`](../params/DF_12.csv) |
| **13** | the same two parameter sets from different initial conditions | two run sets, five initial conditions each |

Defect density is counted per unit area of the padded domain, (1008 × 0.01)² =
101.6064, not of the nominal L = 10. The figure 9 threshold "defect density
larger than 1.5" is a count of 150.

## Still to do

**Figure 13** and several data files have no generator here: `omegarho.csv` and
`zrc_eq{,2}.csv` (figure 5), `zrc_plus.csv` (figure 4's critical-activity dots,
which come from the same stability analysis), and `DF_tikz.csv` (figure 7).

**The article's Eq. for the critical activity and `LSA.jl` do not agree.** The
article reads

    zeta_c = A(chi=0)/(3 rho0) + ( sqrt(A(chi=0) gamma / (3 rho0^3)) + sqrt(2/(3 rho0^3 tau)) )^2

while `ζc_simp` in [`LSA/LSA.jl`](../LSA/LSA.jl), which reproduces figure 5's
`tau1.csv` and `tau5.csv` byte for byte, keeps the χ term in `A` instead of
setting χ = 0, divides by `3 rho0^2` rather than `3 rho0^3`, and **subtracts**
the second root rather than adding it. One of the two is not what the other
describes — worth resolving before publication.

## What can be run

```bash
julia --project=. LSA/linear_stability.jl     # figure 5's data
```

and the figure 3 notebook in [`1D/`](../1D). The simulations themselves are
reproducible from [`2D/`](../2D) with the tables in [`params/`](../params).
