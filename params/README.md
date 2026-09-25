# Input parameters, per figure

`DF_<N>.csv` holds the simulations behind figure N — one row each, in the
schema [`2D/InputParameters.jl`](../2D/InputParameters.jl) reads. Regenerate
them all with:

```bash
julia --project=. params/make_dataframes.jl
```

To run a figure's set, point the solver at its table:

```bash
DF_FILE=../params/DF_9.csv DATA_DIR=/scratch/fig9/ SIM_IDX=1 \
    julia --project=. --optimize=3 2D/2D.jl
```

`DF_FILE` is relative to `2D/`; unset, the solver reads `2D/DF.csv` as before.

## Where the values come from

Table I of the article fixes everything common to every run: `ap` (χ) = 0.1,
`nu1` (ν) = 0, `gamma` (Γ) = 1, `kp` (κ) = 10⁻⁴, `M` (γ) = 10⁻⁴, `zetap` = 0,
`xi` = 1, and `ar` = 4ζ_ρ/3, falling back to 4/3 when ζ_ρ = 0. The numerical
appendix fixes `dx` = `dz` = 10⁻² and `dtmin` (Δt_max) = 10⁻². Each caption
then supplies what that figure varies.

**The article writes the renewal time τ; the solver takes the rate `kd` = 1/τ.**

| File | Rows | Figure varies |
|---|---|---|
| `DF_1.csv` | 3 | L=50, ρ₀ ∈ {0.45, 0.6, 0.75}, τ=5, ζ_ρ=4 |
| `DF_2.csv` | 150 | L=10, ρ₀=1, ζ_ρ=4, τ ∈ {10, 1, 0.2}, 50 initial conditions each — **panel (a) only** |
| `DF_4.csv` | 1008 | L=10, 12 ρ₀ × 14 ζ_ρ × 6 τ |
| `DF_6.csv` | 8 | L=10, ρ₀ ∈ {0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 1.2}, τ=5, ζ_ρ=4 |
| `DF_8.csv` | 2 | L=10, (ρ₀=0.7, τ=1, ζ_ρ=10) and (ρ₀=1.3, τ=0.2, ζ_ρ=1) |
| `DF_7.csv` | 120 | L=10, 12 ρ₀ × 10 τ, ζ_ρ=4 |
| `DF_9.csv` | 1680 | L=10, 12 ρ₀ × 14 ζ_ρ × 10 τ |
| `DF_10.csv` | 3 | L=10, ζ_ρ=4, (ρ₀=0.7, τ=1), (1.2, 0.5), (0.7, 5) |
| `DF_11.csv` | 13 | L=50, ρ₀ = 0.40 : 0.05 : 1.00, τ=5, ζ_ρ=4 |
| `DF_12.csv` | 2 | the same two solutions as Fig8 |
| `DF_13.csv` | 10 | Fig8's two solutions from five initial conditions each |

`DF_11.csv` reproduces [`2D/DF.csv`](../2D/DF.csv) exactly. `DF_9.csv` was
checked against the table the sweep actually used and matched all 1680
parameter sets; that table is not published here, so the test suite pins the
three grids instead. `DF_4.csv` has 1008 rows, exactly the row count of
[`figures/Fig4/Ndef.csv`](../figures/Fig4/Ndef.csv).

## Not generated, and why

### Figure 2, panels (b) to (d)

Those panels are two-defect runs: a pair of oppositely charged defects placed a
set distance apart, rather than the noisy uniform state everything else starts
from. Reproducing them needs two things this repository does not have.

The initial condition is one. `2D/2D.jl` starts every run from a small random
perturbation; seeding a defect pair needs a different initialiser, which is not
part of the published solver.

The separation is the other. Those runs sweep an initial defect distance `D`
alongside ρ₀ and ζ_ρ, and `D` is not a column of this schema. Panel (b,d) is a
single solution at ρ₀ = 0.6, τ = 1, ζ_ρ = 12; panel (c) sweeps ρ₀ ∈ {0.6, 0.8,
1.0, 1.2} against ζ_ρ ∈ {1, 4, 8, 12} at τ = 1, over twelve separations.

### Not a simulation table

- **Fig3** is the 1D model, whose inputs are notebook fields rather than a row
  of `DF.csv`: χ=1, γ=10⁻², τ=0.2, κ=10⁻³, ζ_ρ=0. Note these differ from
  Table I, which is correct — it is a different reduction — but the τ range for
  the inset sweep is not given.
- **Fig5** is linear stability analysis, computed from the equations rather than
  from simulation output. See [`Analysis/make_fig_LSA.jl`](../Analysis/make_fig_LSA.jl).

## The `ar` column

`ar` is the compressibility coefficient a', which Table I ties to the activity:
a' = 4ζ_ρ'/3, or 4/3 when ζ_ρ = 0. The solver reads it from the table.

It used not to. `2D/InputParameters.jl` recomputed it and ignored the column,
and the two were associated differently — `AllInputParam.jl` scaled a
precomputed 4/3, giving `zr*(4/3)`, while the solver computed `(zr*4)/3`. Those
differ in the last bit for ζ_ρ = 10, 14 and 20, so the column never held the
value the run used, and editing it had no effect.

All three now associate it the same way, as `(zr*4)/3`: the sweep writes it,
these tables carry it, and the solver reads it. The change moved no run — the
value every simulation used is what the solver computed, which is exactly what
the column now holds. The archival `2D/DF.csv` is unaffected and still
regenerates byte for byte, because it sweeps only ζ_ρ = 4, where the two
associations agree.

## The integration horizon

`t_fin`, `t_prin` and `t_check` are columns of these tables, so a figure's run
is fully described by its row. `2D/InputParameters.jl` reads them when present
and otherwise falls back to 150000 / 1000 / 1, the values the paper's L=50 runs
used — so the archival [`2D/DF.csv`](../2D/DF.csv), which has no such columns,
behaves exactly as before.

`t_check` is 1 everywhere: Δt is retuned once per unit of simulated time.

`t_fin` and `t_prin` are taken from the configuration of the run set each
figure came from:

| Figure | t_fin | t_prin | run set |
|---|---|---|---|
| 1, 11 | 150 000 | 1 000 | the L=50 series |
| 4, 9 | 200 000 | 1 000 | the 1680-row sweep |
| 6 | 200 000 | 1 000 | the L=10 detailed-density runs |
| 8, 12, 13 | 2 000 000 | 10 000 | the seven-day lattice runs |
| 2 | 150 000 | 1 000 | the saturation runs |
| 7, 10 | 300 000 | 1 000 | the phase-diagram runs; Fig10's t_f = 287·10³ fits inside |

Note `t_prin` is 10 000 rather than 1 000 for the lattice figures: the long runs
snapshot ten times less often.

**`t_fin` is the configured cap, not the length of every run, and it does not
contradict the article.** The sweeps were set to 200 000, while the article
quotes a typical simulated time of 10⁵. Both are right: the cluster sometimes
went down mid-run, and runs that had already reached a steady state were not
restarted, so a number of them stop somewhere between 10⁵ and the cap. 10⁵ is
what every run is guaranteed to have reached; 200 000 is what to configure to
reproduce them. Expect a reproduction to run longer than some of the original
solutions did, and do not treat the difference as an error in either place.

Figure 2's two-defect runs, which have no table here, used `t_fin` = 30 000
with `t_prin` = 500, alongside a much shorter set at 1 000 / 20 for the
critical-distance measurement.
