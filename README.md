# Defect states in compressible active polar fluids with renewal

Simulation code for

> L. Dumoulin, C. Blanch-Mercader, K. Kruse,
> *Defect states in compressible active polar fluids with renewal*,
> PRX Life (2026), [doi:10.1103/bnx5-sryc](https://doi.org/10.1103/bnx5-sryc);
> [arXiv:2506.03795](https://arxiv.org/abs/2506.03795).

> **To explore the model, start with [HydraFluids](https://github.com/Lu-Dumoulin/HydraFluids).**
> It implements the same equations with an interactive
> [Pluto](https://plutojl.org/) interface, runs on a laptop CPU as well as on
> a GPU, and needs no cluster. It is the easiest way to play with the system
> and recovers qualitatively the same states as in the article — defect
> lattices, foams, waves, homogeneous polar states. It also supports nematic
> and nematopolar order, and is actively maintained.

This repository contains the code used to produce the results of the article
itself, and is the place to look for the exact simulations behind each figure:

- **[`2D/`](2D/)** — the 2D solver, a pseudo-spectral GPU code written with
  [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl);
- **[`params/`](params/)** — the parameters of every simulation shown in the
  figures, one table per figure;
- **[`1D/`](1D/)** — the reduced 1D model of a single defect (figure 3), as a
  [Pluto](https://plutojl.org/) notebook;
- **[`LSA/`](LSA/)** — the linear stability analysis of the homogeneous states;
- **[`figures/`](figures/)** — the LaTeX sources and data of the figures.

The 2D solver needs an NVIDIA GPU, and the large simulations of the article need
a GPU cluster.

---

## Contents

```
.
├── 2D/
│   ├── 2D.jl              entry point: initial condition, time loop, output
│   ├── kernels.jl         GPU kernels: free energy, stress, force balance, polarity
│   ├── InputParameters.jl reads one row of a parameter table, sets up grid and FFTs
│   ├── AllInputParam.jl   builds DF.csv, the L=50 density series
│   └── DF.csv             the L=50 density series (figures 1 and 11)
├── params/                DF_N.csv, the simulations behind figure N
├── 1D/                    1D model notebook and the data of figure 3
├── LSA/                   linear stability analysis of the homogeneous states
├── figures/               one directory per figure: LaTeX source and data
├── slurm/submit.sh        Slurm array job
├── test/                  test suite
├── Utilities/             helper functions used by the simulation scripts
└── Project.toml
```

---

## The model

A compressible polar fluid on a 2D periodic square domain, described by the
density $\rho$, the polarity $\mathbf{P}$ and the velocity $\mathbf{v}$.

**Density** — advection, diffusion driven by the chemical potential $\mu$, and
renewal at rate $k_d = 1/\tau$ towards the target density $\rho_0$:

$$\partial_t \rho = -\nabla\cdot(\rho\mathbf{v}) + \gamma \Delta\mu + k_d(\rho_0-\rho)$$

**Polarity** — flow alignment and relaxation along the molecular field
$\mathbf{h}$, with $D/Dt$ the co-rotational derivative (advection and
rotation by the flow) and $\mathbf{u}$ the strain rate:

$$\frac{D\mathbf{P}}{Dt} = -\nu \mathbf{P}\cdot\mathbf{u} + \frac{1}{\Gamma}\mathbf{h}$$

**Force balance** — viscous flow with substrate friction $\xi$, driven by the
Ericksen stress and the active stresses $\zeta_\rho \Delta\mu \rho^3$ and
$\zeta_p \Delta\mu \rho P_i P_j$.

The full expressions are derived in the article and implemented in
[`2D/kernels.jl`](2D/kernels.jl). Spatial derivatives and the force balance are
computed in Fourier space; time integration is explicit Euler with an adaptive
time step.

---

## Requirements

- Julia ≥ 1.10
- an NVIDIA GPU for the 2D solver. The L=10 simulations (1008² grid points)
  run on any recent card; the L=50 ones (5008²) need about 40 GB of GPU memory.

Install the dependencies once, from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The 1D notebook and the linear stability analysis do not need a GPU.

---

## Running a simulation

Every simulation is one row of a parameter table. To run row 1 of the table
behind figure 9:

```bash
cd 2D
DF_FILE=../params/DF_9.csv DATA_DIR=/path/to/output/ SIM_IDX=1 \
    julia --project=.. --optimize=3 2D.jl
```

| Variable | Meaning | Default |
|---|---|---|
| `DF_FILE` | parameter table, relative to `2D/` | `DF.csv` |
| `SIM_IDX` | row of the table to run | `SLURM_ARRAY_TASK_ID`, else 1 |
| `DATA_DIR` | output directory | `data/L50/` at the repository root |

### On a Slurm cluster

To run a whole table, one task per row:

```bash
export DF_FILE=../params/DF_9.csv
export DATA_DIR=/path/to/output/
sbatch --chdir=2D --array=1-1680%20 slurm/submit.sh
```

`--array` must match the number of rows of the table. Adapt the partition,
GPU constraint and time limit in [`slurm/submit.sh`](slurm/submit.sh) to your
cluster.

### Output

Simulation `n` writes `DATA_DIR/n/Data/data<t>.jld` every `t_prin` time units,
where `<t>` is the time, zero-padded to 10 digits. Each file contains the
fields `rho`, `v` and `P` as `Float64` arrays. The run stops with status 1 if
the density becomes NaN.

### Parameter tables

Each row of a table gives, in the notation of the article (Table I):

| Column | Symbol | Meaning |
|---|---|---|
| `L` | $L$ | system size |
| `rho0` | $\rho_0$ | target density |
| `kd` | $1/\tau$ | renewal rate — **the article uses the renewal time $\tau$** |
| `zetarho` | $\zeta_\rho$ | isotropic active stress |
| `zetap`, `zetap2` | $\zeta_p$ | anisotropic active stress |
| `ar` | $a$ | $4\zeta_\rho/3$, or $4/3$ when $\zeta_\rho = 0$ |
| `ap` | $\chi$ | polar ordering coefficient |
| `kp` | $\kappa$ | Frank constant |
| `M` | $\gamma$ | density mobility |
| `gamma` | $\Gamma$ | rotational mobility |
| `nu1` | $\nu$ | flow alignment |
| `xi` | $\xi$ | substrate friction |
| `seed` | | random seed of the initial condition |
| `dx`, `dz` | $\Delta x$ | grid spacing |
| `dtmin` | $\Delta t_\text{max}$ | largest allowed time step |
| `t_fin` | | final time |
| `t_prin` | | interval between saved snapshots |
| `t_check` | | interval between time-step updates |

The initial condition is the uniform isotropic state with a small random
perturbation: $\rho = \rho_0(1 + 0.002 \eta)$, $\mathbf{P} = 0.002 \eta'$,
with $\eta, \eta'$ uniform in $[-0.5, 0.5]$. For a given seed, a simulation is
reproducible bit for bit on the same GPU model and CUDA version.

**Two-defect simulations.** Adding a column `D` to a table starts the
simulation instead from a pair of oppositely charged defects, a distance `D`
apart (in units of the system size). The simulation then stops as soon as the
two defects either annihilate or move apart. With `D = 0`, the code searches
for the critical distance: it increases `D` from 0.01 in steps of 0.01 and
reports the first value at which the two defects do not annihilate. This is
used for figure 2(c).

See [`params/README.md`](params/README.md) for the table behind each figure.

---

## Reproducing the figures

| Figure | How |
|---|---|
| 1, 2(a), 4, 6–13 | run the simulations of [`params/DF_N.csv`](params/) |
| 2(b–d) | two-defect simulations, see above |
| 3 | [`1D/models.pluto.jl`](1D/models.pluto.jl) |
| 5 | linear stability analysis, see below; the data plotted is in [`figures/Fig5/`](figures/Fig5/) |

The analysis of the simulation output is described in the appendices of the
article. [`figures/README.md`](figures/README.md) lists what each figure shows
and which simulations it uses, and explains how to compile the figures.

### Figure 3: the 1D model

```julia
using Pkg; Pkg.add("Pluto"); import Pluto; Pluto.run()
```

then open [`1D/models.pluto.jl`](1D/models.pluto.jl). The notebook installs its
own dependencies on first use. Set the parameters, turn on the **ready** switch
and press **Confirm**; two checkboxes at the bottom write the data of figure 3
into `1D/figdata/`.

The notebook contains a radial model of a single defect, used for figure 3, and
a Cartesian model of a defect pair. Its symbols relate to those of the 2D code
as `A` = `ar`, `χ` = `ap`, `κ` = `kp`, `τ` = `1/kd`.

### Linear stability analysis

```bash
julia --project=. LSA/linear_stability.jl
```

linearises the dynamic equations about the homogeneous polar state and computes
the critical activity $\zeta_\rho^c$ as a function of $\rho_0$, for $\tau = 1$
and $\tau = 5$. It writes `tau1.csv` and `tau5.csv` into `figures/Fig5/`, each
holding a closed-form estimate (`y`) and the threshold found by integrating the
linearised equations (`yy`). It needs no simulation output and runs in under a
minute.

The data plotted in figure 5 is provided in [`figures/Fig5/`](figures/Fig5/):
`zrc_eq2.csv` for panel (a) and `omegarho.csv` for panel (b). The closed-form
expression of the critical activity used in figure 4 is given in the article.

---

## Tests

```bash
julia --project=. test/runtests.jl
```

checks the parameter tables, the solver sources, the linear stability analysis
and the compilation of every figure. It does not run the GPU solver. The figure
tests are skipped if `pdflatex` is not installed.

---

## Notes

- The simulation code is the one used for the article, adapted to run outside
  the cluster it was written for: paths and run parameters are read from the
  environment and the parameter table, and the arithmetic is unchanged.
- No `Manifest.toml` is provided, since the simulations ran on Linux with CUDA;
  package versions are therefore not pinned.
- `Utilities/` holds helper functions shared by several of our projects, copied
  unchanged.

## Citing

```bibtex
@article{dumoulin2026defect,
  title         = {Defect states in compressible active polar fluids with renewal},
  author        = {Dumoulin, Ludovic and Blanch-Mercader, Carles and Kruse, Karsten},
  journal       = {PRX Life},
  year          = {2026},
  doi           = {10.1103/bnx5-sryc},
  eprint        = {2506.03795},
  archivePrefix = {arXiv}
}
```

## License

MIT — see [LICENSE](LICENSE).
