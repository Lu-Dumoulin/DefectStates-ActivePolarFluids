# Defect states in compressible active polar fluids with turnover

CUDA simulation code for

> L. Dumoulin, C. Blanch-Mercader, K. Kruse,
> *Defect states in compressible active polar fluids with turnover*,
> [arXiv:2506.03795](https://arxiv.org/abs/2506.03795) (2025).

This repository is an archival snapshot of the code that produced the published
simulations. It holds two independent pieces:

- **[`FFT_2D_P_L50/`](FFT_2D_P_L50/)** — the 2D GPU run set behind the main
  results: a pseudo-spectral solver written directly against
  [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) and run as a Slurm array job on
  the [Baobab](https://doc.eresearch.unige.ch/hpc/start) cluster at the
  University of Geneva.
- **[`1D/`](1D/)** — reduced 1D models of a single defect core and of a defect
  pair, as an interactive Pluto notebook, together with the data and sources for
  the figures they feed. These accompany the revised version of the paper.

A more general, better-documented and actively maintained implementation of the
same physics — CPU/CUDA/Metal via ParallelStencil.jl, with a Pluto front end and
support for nematic and nematopolar order — lives in
[HydraFluids](https://github.com/Lu-Dumoulin/HydraFluids). **If you want to build
on this model, start there.** Use this repository if you want the exact code
behind the figures.

---

## The model

A compressible polar fluid on a 2D periodic square domain, with density $\rho$,
polarity $\mathbf{P}$ and velocity $\mathbf{v}$.

**Density** — advection, diffusion through the chemical potential $\mu$, and
turnover at rate $k_d$ towards the homeostatic density $\rho_0$:

$$\partial_t \rho = -\nabla\cdot(\rho\mathbf{v}) + M\,\Delta\mu + k_d(\rho_0-\rho)$$

**Polarity** — advection, co-rotation, flow alignment ($\nu_1$), and relaxation
along the molecular field $\mathbf{h}$ at rate $\gamma$:

$$\partial_t \mathbf{P} = -(\mathbf{v}\cdot\nabla)\mathbf{P} - \boldsymbol{\omega}\cdot\mathbf{P} - \nu_1 \mathbf{P}\cdot\mathbf{u} + \gamma\mathbf{h}$$

**Momentum** — Stokes flow with substrate friction $\xi$, solved in Fourier
space, with the stress carrying the Ericksen and flow-alignment contributions
plus the active terms $\zeta_\rho \Delta\mu\, \rho^3$ and
$\zeta_p \Delta\mu\, \rho P_i P_j$.

The exact expressions for $\mu$, $\mathbf{h}$ and $\sigma_{ij}$ are in
[`FFT_2D_P_L50/kernels.jl`](FFT_2D_P_L50/kernels.jl); see the paper for their
derivation.

Spatial derivatives and the Stokes solve are spectral (`rfft` along $x$, `fft`
along $z$); time stepping is explicit Euler with an adaptive step capped by a
CFL condition on $\mathbf{v}$ and by `dtmin`.

---

## Layout

```
.
├── FFT_2D_P_L50/          # the simulation, exactly as run
│   ├── AllInputParam.jl   # defines the sweep, writes DF.csv
│   ├── DF.csv             # the 13 parameter sets used for the paper
│   ├── InputParameters.jl # reads one DF.csv row, sets up grid and FFT plans
│   ├── kernels.jl         # the CUDA kernels
│   ├── 2D.jl              # entry point: fields, time loop, snapshot writing
│   └── MakePlots.jl       # figures from the .jld snapshots
├── 1D/                    # reduced 1D models (revised version of the paper)
│   ├── models.pluto.jl    # Pluto notebook: radial and two-defect solvers
│   ├── figdata/           # CSVs it writes, read directly by the figures
│   ├── fig_active_polar.tex        # pgfplots figure built from figdata/
│   ├── fig_active_polar_loglog.tex # log-log version of the same
│   ├── fig_core_analytic.tex       # core profiles vs the analytic small-r form
│   └── standalone_fig.tex          # wrapper to compile a figure on its own
├── Utilities/             # lab-internal helpers the scripts include
│   ├── using.jl           # using_pkg / using_mod (install-on-demand)
│   ├── JulUtils.jl        # generate_dataframe, file helpers
│   └── PictUtils.jl       # png/gif helpers used by MakePlots.jl
├── slurm/submit.sh        # the array job
└── Project.toml
```

The directory layout is preserved from the cluster so that the `include("../Utilities/...")`
paths in the simulation files work unchanged.

---

## The parameter set

`DF.csv` holds one row per simulation. The sweep is over the homeostatic density
$\rho_0 \in [0.4, 1.0]$ in steps of $0.05$ (13 runs), everything else fixed:

| Column | Symbol | Value | Meaning |
|---|---|---|---|
| `L` | $L$ | 50 | box side |
| `dx`, `dz` | $\Delta x$ | 0.01 | grid spacing → 5000 points/side, padded to 5008 (multiple of 16) |
| `dtmin` | | 0.01 | upper bound on the adaptive time step |
| `rho0` | $\rho_0$ | 0.40 … 1.00 | **swept** — homeostatic density |
| `kd` | $k_d$ | 0.2 | turnover rate |
| `ap` | $a_p$ | 0.1 | polar ordering coefficient |
| `kp` | $k_p$ | 1e-4 | Frank constant |
| `gamma` | $\gamma$ | 1.0 | rotational mobility |
| `nu1` | $\nu_1$ | 0.0 | flow alignment |
| `zetarho` | $\zeta_\rho$ | 4.0 | isotropic active stress |
| `zetap`, `zetap2` | $\zeta_p$ | 0.0 | anisotropic active stress |
| `ar` | $a_r$ | 16/3 | set to $\lvert\zeta_\rho\rvert\cdot 4/3$ by `AllInputParam.jl` |
| `xi` | $\xi$ | 1.0 | substrate friction |
| `M` | $M$ | 1e-4 | mobility (divided by $a_r$ at load time) |
| `seed` | | 1 | RNG seed for the initial noise |

Initial condition: $\rho = \rho_0(1 + 0.002\,\eta)$ and
$\mathbf{P} = 0.002\,\eta'$, with $\eta,\eta'$ uniform in $[-0.5,0.5]$ from
`Random.seed!(seed)` — so runs are reproducible bit-for-bit on the same GPU and
CUDA version.

`t_fin`, `t_prin` and `t_check` (integration horizon, snapshot interval and
adaptive-step interval) are set at the bottom of `InputParameters.jl`, not in
`DF.csv`. The paper runs used `t_fin = 150000`, `t_prin = 1000` against a 12 h
wall-clock limit.

To regenerate `DF.csv` from the sweep definition:

```bash
julia --project=. FFT_2D_P_L50/AllInputParam.jl
```

The committed `DF.csv` is the one used for the paper; regenerating overwrites it.

---

## Running

Requires Julia ≥ 1.10 and an NVIDIA GPU. The paper runs used A100 40 GB / 80 GB
cards with Julia supplied by Baobab's `module load Julia`. Note the memory
footprint: at 5008² in `Float64`, the field and Fourier arrays come to several
GB of device memory, so the grid as configured needs a 40 GB-class card.

**One simulation, by hand:**

```bash
export DATA_DIR=/path/to/output/    # trailing slash
export SIM_IDX=7                    # row of DF.csv
julia --project=. --optimize=3 FFT_2D_P_L50/2D.jl
```

**The full array on Slurm:**

```bash
sbatch --chdir=FFT_2D_P_L50 slurm/submit.sh
```

`submit.sh` is reconstructed from the lab's `Code2Cluster.jl` job generator with
the settings used for these runs (13-task array, 20 concurrent, 12 h, one A100).
Adjust the partition and constraint for your site.

**Output.** Simulation `idx` writes `<DATA_DIR>/<idx>/Data/data<t>.jld` every
`t_prin` time units, where `<t>` is the integration time zero-padded to 10
digits. Each file holds `rho` (Nx×Nz), `v` and `P` (Nx×Nz×2) as `Float64`. The
run aborts with status 1 if `rho` goes NaN.

**Figures:**

```bash
export DATA_DIR=/path/to/output/    # must also contain DF.csv
export FIG_DIR=/path/to/figures/
julia --project=. FFT_2D_P_L50/MakePlots.jl
```

---

## The 1D models

[`1D/models.pluto.jl`](1D/models.pluto.jl) is a self-contained
[Pluto](https://plutojl.org/) notebook holding reduced models built from the
same free energy as the 2D solver, with the symbols renamed:

| 1D | 2D | |
|---|---|---|
| `A` | `ar` | compressibility |
| `χ` | `ap` | polar ordering |
| `κ` | `kp` | Frank constant |
| `τ` | `1/kd` | turnover time |

$$f = \tfrac{A}{4}\rho^4 - \tfrac{\chi}{2\rho_0}\rho^3 p^2 + \tfrac{\chi}{4}\rho^2 p^4 + \tfrac{\kappa}{2}\rho^2\Big[(\partial_r p)^2 + \tfrac{p^2}{r^2}\Big]$$

Three solvers are selectable in the notebook:

- **One defect (polar, radial)** — the 1D radial reduction around a single
  defect core, on a staggered grid (ρ at cell centres, `p` and `v` at faces)
  with finite volumes and IMEX implicit diffusion. This is the one the committed
  figure data comes from.
- **Two defects (explicit)** and **(implicit)** — a Cartesian pseudo-spectral
  solver for a defect pair, closer in spirit to the 2D code.

It also sweeps the turnover time τ and records the velocity field `v(r, τ)`.

### Running it

```julia
using Pkg; Pkg.add("Pluto"); import Pluto; Pluto.run()
```

then open `1D/models.pluto.jl`. The notebook pins its own package versions in
the embedded `PLUTO_PROJECT_TOML`/`PLUTO_MANIFEST_TOML` blocks, so it needs
nothing from this repository's `Project.toml` and installs its own environment
on first run.

Parameters are set with the sliders, and the run starts when you flip the
**ready** switch and press *Confirm*. Two checkboxes near the bottom write the
CSVs into `figdata/`, overwriting what is committed here.

### Figure data

The committed `figdata/` is the data the figures use, exactly as the notebook
wrote it. `p` and `v` live on the N+1 cell faces and ρ on the N centres, which
is why the row counts differ by one.

| File | Columns | |
|---|---|---|
| `polar_fields_rho.csv` | `r, rho` | density profile through the core |
| `polar_fields_Pv.csv` | `r, P, v` | polarity and velocity profiles |
| `polar_balance.csv` | `r, transport, diffusive, turnover` | the three terms of the density balance |
| `polar_vsweep_summary.csv` | `tau, vmax, r_vmax, rdiv, v_rdiv` | one row per τ |
| `polar_vsweep_heatmap.csv` | `r, tau, v` | full v(r, τ), 301 radii × 50 τ values |

### Figures

The `.tex` files are pgfplots figures that read those CSVs directly (comma
separated, via `\addplot table`):

- `fig_active_polar.tex` and `fig_active_polar_loglog.tex` — the profiles, the
  balance and the τ-sweep summary
- `fig_core_analytic.tex` — the core profiles against the analytic small-r form
  ρ ≈ ρ̄(1 + r²/λ_ρ²), p ≈ r/λ_p, v ≈ r/τ_v

`standalone_fig.tex` wraps one of them for compiling on its own, but it also
`\input`s a `some_command.tex` of shared macros that lives in the paper's LaTeX
project and is **not** included here — supply your own, or strip that line.

---

## Relation to the code as it was run

The simulation files are byte-identical to the cluster copies, with three
exceptions, each marked in place:

1. **`InputParameters.jl`** — the three hard-coded cluster paths (`dir`,
   `localpath`, `idx` from `SLURM_ARRAY_TASK_ID`) now read from environment
   variables with defaults. The original values are kept as comments.
2. **`InputParameters.jl`** — added `mkpath(joinpath(file, "Data"))`. The
   original code never created this directory; on the cluster it already existed
   from earlier runs, so a fresh checkout would have failed on the first
   snapshot write.
3. **`MakePlots.jl`** — the two hard-coded figure paths likewise read from the
   environment.

`AllInputParam.jl` and `InputParameters.jl` are otherwise only commented and
stripped of commented-out dead code; both were checked by parsing the original
and the annotated version and comparing syntax trees. `AllInputParam.jl` still
regenerates the committed `DF.csv` byte-for-byte.

`kernels.jl` and `2D.jl` have additionally been tidied and de-allocated. The
intent is that results stay bit-for-bit identical, so every change was chosen to
leave the floating-point arithmetic alone:

- **`kernels.jl`** — all 38 arithmetic expressions are carried over verbatim
  (verified by extracting and diffing them). What changed is `@inbounds` on the
  three kernels that lacked it, a shared inlined `thread_indices()` in place of
  the index arithmetic repeated in all four, and docstrings.
- **`2D.jl`** — dropped the unused `ρn` and `F` allocations (~600 MB of device
  memory); replaced `minimum([a,b,c])` with `min(a,b,c)`, dropped an `abs.` on
  a sum of squares and a multiply by an array of ones, and fused three
  broadcasts that each materialised a full-grid temporary. Each of these four
  was checked bitwise on CPU over random data and the signed-zero / tie / Inf /
  NaN edge cases. Snapshot writing moved into `write_snapshot`, and the noise
  amplitude, step-growth factor and CFL fraction became named constants holding
  the same literals.

**None of this has been run on a GPU** — it was prepared on a machine without
one. Before trusting it, confirm on a GPU node:

```bash
tools/verify_bitwise.sh HEAD~1 1 20 10
```

which runs one simulation for a couple of snapshots under both the previous and
the current revision and compares every dataset bit for bit.

Not done: the time loop still allocates roughly 9 GB of full-grid temporaries
per step, because each `W * x` and `Wi * (factor .* x)` materialises a fresh
array. Reusing preallocated buffers via `mul!` would remove nearly all of it,
but the inverse plan may clobber its input and the destinations are views, so
that change needs a GPU to validate rather than reasoning alone.

`Project.toml` and `slurm/submit.sh` are new; the cluster runs installed packages
on demand into the default environment via `Utilities/using.jl`, and the job
script was generated on the fly. No `Manifest.toml` is shipped, because one
resolved off-cluster would pin the wrong CUDA artifacts — this does mean package
versions are not pinned to those used for the paper.

Everything under `1D/` is copied verbatim and is not annotated: Pluto manages
the cell markers and cell-order block in `models.pluto.jl`, and the notebook
already carries its equations in markdown cells.

`Utilities/` is vendored verbatim from a lab-internal directory. `JulUtils.jl`
pulls in Makie for a `screensize()` helper that nothing here calls; it is kept so
the file matches the original.

---

## Citing

```bibtex
@article{dumoulin2025defect,
  title  = {Defect states in compressible active polar fluids with turnover},
  author = {Dumoulin, Ludovic and Blanch-Mercader, Carles and Kruse, Karsten},
  journal = {arXiv preprint arXiv:2506.03795},
  year   = {2025}
}
```

## License

MIT — see [LICENSE](LICENSE).
