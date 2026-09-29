# Parameters of the simulations, per figure

`DF_N.csv` lists the simulations used in figure N of the article, one per row.
To run row 1 of the table behind figure 9, from `2D/`:

```bash
DF_FILE=../params/DF_9.csv DATA_DIR=/path/to/output/ SIM_IDX=1 \
    julia --project=.. --optimize=3 2D.jl
```

and see the main [README](../README.md#on-a-slurm-cluster) to run a whole table
on a cluster. The columns are described in the main README.

All tables use the common values of Table I of the article — $\chi = 0.1$,
$\kappa = 10^{-4}$, $\gamma = 10^{-4}$, $\Gamma = 1$, $\nu = 0$, $\zeta_p = 0$,
$\xi = 1$, $a = 4\zeta_\rho/3$ — and the discretisation of the numerical
appendix, $\Delta x = 10^{-2}$ and $\Delta t_\text{max} = 10^{-2}$.

**The article gives the renewal time $\tau$; the tables give the renewal rate
`kd` $= 1/\tau$.**

## The tables

| File | Figure | Runs | Parameters | `t_fin` | `t_prin` |
|---|---|---|---|---|---|
| `DF_1.csv` | 1 | 3 | $L = 50$, $\rho_0 \in \{0.45, 0.6, 0.75\}$, $\tau = 5$, $\zeta_\rho = 4$ | 150 000 | 1 000 |
| `DF_2.csv` | 2(a) | 150 | $L = 10$, $\rho_0 = 1$, $\zeta_\rho = 4$, $\tau \in \{0.2, 1, 10\}$, 50 seeds each | 150 000 | 1 000 |
| `DF_4.csv` | 4 | 1008 | $L = 10$, $\rho_0 = 0.4, 0.5, \ldots, 1.5$, $\zeta_\rho \in \{0, 1, 2, 4, 6, \ldots, 24\}$, $\tau \in \{0.1, 0.2, 1, 5, 10, 100\}$ | 200 000 | 1 000 |
| `DF_6.csv` | 6 | 8 | $L = 10$, $\rho_0 \in \{0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 1.2\}$, $\tau = 5$, $\zeta_\rho = 4$ | 200 000 | 1 000 |
| `DF_7.csv` | 7 | 120 | $L = 10$, $\rho_0 = 0.4, 0.5, \ldots, 1.5$, 10 values of $\tau$ from 0.1 to 100, $\zeta_\rho = 4$ | 300 000 | 1 000 |
| `DF_8.csv` | 8 | 2 | $L = 10$; $(\rho_0, \tau, \zeta_\rho) = (0.7, 1, 10)$ and $(1.3, 0.2, 1)$ | 2 000 000 | 10 000 |
| `DF_9.csv` | 9 | 1680 | $L = 10$, the grid of `DF_4.csv` with all 10 values of $\tau$ | 200 000 | 1 000 |
| `DF_10.csv` | 10 | 3 | $L = 10$, $\zeta_\rho = 4$; $(\rho_0, \tau) = (0.7, 1)$, $(1.2, 0.5)$, $(0.7, 5)$ | 300 000 | 1 000 |
| `DF_11.csv` | 11 | 13 | $L = 50$, $\rho_0 = 0.40, 0.45, \ldots, 1.00$, $\tau = 5$, $\zeta_\rho = 4$ | 150 000 | 1 000 |
| `DF_12.csv` | 12 | 2 | the two simulations of figure 8 | 2 000 000 | 10 000 |
| `DF_13.csv` | 13 | 10 | the two simulations of figure 8, with seeds 1 to 5 | 2 000 000 | 10 000 |

The ten renewal times are $\tau$ = 0.1, 0.2, 0.5, 1, 1.25, 1.67, 2.5, 5, 10 and
100.

`t_fin` is the time the simulations were set to reach. Some of the original
simulations stopped earlier, after reaching a steady state; every one reached
at least $t = 10^5$.

## Figures without a table

- **Figure 2(b–d)** uses two-defect simulations, which start from a pair of
  defects rather than from a random perturbation; see the main
  [README](../README.md#parameter-tables). Panels (b) and (d) show one
  simulation at $\rho_0 = 0.6$, $\tau = 1$, $\zeta_\rho = 12$, $L = 10$, with
  `t_fin` = 30 000 and `t_prin` = 500. Panel (c) gives the critical distance
  for $\rho_0 \in \{0.6, 0.8, 1.0, 1.2\}$ and $\zeta_\rho \in \{1, 4, 8, 12\}$
  at $\tau = 1$, $L = 10$, obtained with `D = 0`, `t_fin` = 1 000 and
  `t_prin` = 20.
- **Figure 3** comes from the 1D notebook in [`1D/`](../1D), with $\chi = 1$,
  $\gamma = 10^{-2}$, $\tau = 0.2$, $\kappa = 10^{-3}$, $\zeta_\rho = 0$.
- **Figure 5** comes from the linear stability analysis in [`LSA/`](../LSA)
  and needs no simulation.

## Regenerating the tables

```bash
julia --project=. params/make_dataframes.jl
```

rebuilds every `DF_N.csv` from the values above.
