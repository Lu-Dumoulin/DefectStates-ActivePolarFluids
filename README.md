# Defect states in compressible active polar fluids with turnover

CUDA simulation code for

> L. Dumoulin, C. Blanch-Mercader, K. Kruse,
> *Defect states in compressible active polar fluids with turnover*,
> [arXiv:2506.03795](https://arxiv.org/abs/2506.03795) (2025).

This repository is an archival snapshot of the code that produced the published
simulations. It is the `FFT_2D_P_L50` run set: a pseudo-spectral solver for a
compressible active polar fluid with turnover, written directly against
[CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) and run as a Slurm array job on
the [Baobab](https://doc.eresearch.unige.ch/hpc/start) cluster at the University
of Geneva.

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

`Project.toml` and `slurm/submit.sh` are new; the cluster runs installed packages
on demand into the default environment via `Utilities/using.jl`, and the job
script was generated on the fly. No `Manifest.toml` is shipped, because one
resolved off-cluster would pin the wrong CUDA artifacts — this does mean package
versions are not pinned to those used for the paper.

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
