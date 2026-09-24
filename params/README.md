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

**`t_fin` and `t_prin` are provisional.** They are what the article says about
run lengths, not what the runs were configured with, and are set in one place —
the `HORIZON` table at the top of `make_dataframes.jl` — so they can be
corrected in one edit followed by a regeneration.

| Figure | t_fin | from |
|---|---|---|
| 1, 11 | 150 000 | the L=50 cluster runs |
| 4, 9 | 100 000 | "a total simulated time of 10⁵" |
| 8, 12 | 2 000 000 | caption: snapshot at t = 2·10⁶ |

Fig9 needs checking: the article reads Γ_p at t = 1.4·10⁵, which is past the
10⁵ horizon assumed here. Fig10's t_f = 287·10³ is stated but the rest of its
row is not, so no table is generated for it.
