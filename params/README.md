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
| `DF_4.csv` | 1008 | L=10, 12 ρ₀ × 14 ζ_ρ × 6 τ |
| `DF_8.csv` | 2 | L=10, (ρ₀=0.7, τ=1, ζ_ρ=10) and (ρ₀=1.3, τ=0.2, ζ_ρ=1) |
| `DF_9.csv` | 1680 | L=10, 12 ρ₀ × 14 ζ_ρ × 10 τ |
| `DF_11.csv` | 13 | L=50, ρ₀ = 0.40 : 0.05 : 1.00, τ=5, ζ_ρ=4 |
| `DF_12.csv` | 2 | the same two solutions as Fig8 |

Two of these are checkable against tables that were actually used, and both
match exactly: `DF_9.csv` reproduces the 1680 parameter sets of
[`Analysis/PDpaper/DF.csv`](../Analysis/PDpaper/DF.csv), and `DF_11.csv`
reproduces [`2D/DF.csv`](../2D/DF.csv). `DF_4.csv` has 1008 rows, exactly the
row count of [`figures/Fig4/Ndef.csv`](../figures/Fig4/Ndef.csv), with the six
renewal times taken from `make_heatmap`'s `tkd = [10, 5, 1, 0.2, 0.1, 0.01]`.

## Not generated, and why

### Missing from the article

These would need a sentence in the manuscript before the table can be written.

- **Fig2** — panel (a) is "for different values of τ" without saying which, and
  panel (c) sweeps ρ₀ and ζ_ρ without giving either range. Panels (b,d) are
  fully specified (ρ₀=0.6, τ=1, ζ_ρ=12, L=10), and (a) states ρ₀=1, ζ_ρ=4,
  L=10, 50 runs per τ — so only the τ list and the (c) grid are missing.
- **Fig6** — "different target densities ρ₀", but not which. τ=5, ζ_ρ=4, L=10
  are given. The Fig11 caption says ρ₀ = 0.45 and 0.55 are *not* among them,
  which bounds but does not determine the set.
- **Fig7** — neither the τ and ρ₀ grids nor ζ_ρ are stated. The committed
  `DF_tikz.csv` has 120 rows on a 10 × 12 index grid, consistent with the ten τ
  and twelve ρ₀ used elsewhere; ζ_ρ is unrecoverable from the figure data.
- **Fig10** — the three solutions give ρ₀ and τ but not ζ_ρ, and no system size.
- **Fig13** — "two different sets of parameters", unnamed. They are plausibly
  Fig8's two, but the caption does not say so, and the five initial conditions
  are not identified by seed.

### Not a simulation table

- **Fig3** is the 1D model, whose inputs are notebook fields rather than a row
  of `DF.csv`: χ=1, γ=10⁻², τ=0.2, κ=10⁻³, ζ_ρ=0. Note these differ from
  Table I, which is correct — it is a different reduction — but the τ range for
  the inset sweep is not given.
- **Fig5** is linear stability analysis, computed from the equations rather than
  from simulation output. See [`Analysis/make_fig_LSA.jl`](../Analysis/make_fig_LSA.jl).

## The integration horizon

`t_fin`, `t_prin` and `t_check` are columns of these tables, so a figure's run
is fully described by its row. `2D/InputParameters.jl` reads them when present
and otherwise falls back to 150000 / 1000 / 1, the values the paper's L=50 runs
used — so the archival [`2D/DF.csv`](../2D/DF.csv), which has no such columns,
behaves exactly as before.

`t_check` is 1 everywhere: Δt is retuned once per unit of simulated time.

`t_fin` and `t_prin` are read from the `InputParameters.jl` of the run set each
figure came from, on the cluster under `Code/`:

| Figure | t_fin | t_prin | run set |
|---|---|---|---|
| 1, 11 | 150 000 | 1 000 | `FFT_2D_P_L50` |
| 4, 9 | 200 000 | 1 000 | `FFT_2D_P_AdptDt` — its `DF.csv` is the same 1680-row sweep |
| 8, 12 | 2 000 000 | 10 000 | `FFT_2D_P_7D`, the seven-day lattice runs |

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

Horizons for figures whose tables are not yet generated, from the same source:
Fig2 ran at `t_fin` = 30 000 with `t_prin` = 500 (`FFT_2D_P_2defects`), with a
second, much shorter set at 1 000 / 20 for the critical-distance measurement
(`FFT_2D_P_2defects2`); Fig7 and Fig10 ran at 300 000 / 1 000
(`FFT_2D_P_phasediag2`), which is consistent with Fig10's stated t_f = 287·10³.
