# =============================================================================
#  Builds DF_<N>.csv, the input parameters for the simulations behind figure N.
#
#      julia --project=. params/make_dataframes.jl
#
#  Each row is one simulation, in the schema 2D/InputParameters.jl reads. To run
#  a figure's set, point the solver at its table:
#
#      DF_FILE=../params/DF_9.csv DATA_DIR=/scratch/fig9/ SIM_IDX=1 \
#          julia --project=. --optimize=3 2D/2D.jl
#
#  Sources, in order of authority:
#    - Table I of the article for the parameters common to every run
#    - each figure's caption for what that figure varies
#    - App. "Numerical methods" for the discretisation
#    - Analysis/PDpaper/DF.csv, the table the phase-diagram runs actually used
#
#  The article writes the renewal time tau; the solver takes the rate kd = 1/tau.
#
#  Figures whose parameters the article does not pin down are listed in
#  README.md rather than guessed at here.
# =============================================================================

include("../Utilities/using.jl")
using_pkg("CSV, DataFrames")

const OUT = @__DIR__

# --- Integration horizon, per figure ----------------------------------------
#  t_check is 1 everywhere: Delta t is retuned once per unit of simulated time.
#  t_fin and t_prin are per figure. PROVISIONAL below - taken from what the
#  article states about run lengths - and to be replaced with the values the
#  runs actually used.
const T_CHECK = 1
const HORIZON = Dict(          # figure => (t_fin, t_prin)
     1 => (150_000, 1_000),    # L=50 series, as the cluster runs were set
     4 => (100_000, 1_000),    # article: "a total simulated time of 1e5"
     8 => (2_000_000, 1_000),  # caption: snapshot at t = 2e6
     9 => (100_000, 1_000),    # Gamma_p is read at t = 1.4e5 - CHECK
    11 => (150_000, 1_000),    # L=50 series
    12 => (2_000_000, 1_000),  # same solutions as Fig8
)

# --- Table I: common to every run -------------------------------------------
const COMMON = (
    rhocr  = 0.0,
    ap     = 0.1,      # chi'
    nu1    = 0.0,      # nu'
    gamma  = 1.0,      # Gamma'
    kp     = 1.0e-4,   # kappa'
    zetap  = 0.0,      # zeta_p'
    zetap2 = 0.0,
    xi     = 1.0,
    M      = 1.0e-4,   # gamma'
    dx     = 0.01,     # App.: Delta x = 1e-2
    dz     = 0.01,
    dtmin  = 0.01,     # App.: Delta t_max = 1e-2
)

# a' = 4 zeta_rho'/3, and 4/3 when zeta_rho = 0  (Table I footnote)
ar_of(zr) = zr == 0 ? 4/3 : abs(zr)*4/3

"""
    sweep(; L, rho0, tau, zetarho, seed=[1], t_fin, t_prin)

One row per combination, in the column order InputParameters.jl expects.
"""
function sweep(; L, rho0, tau, zetarho, seed=[1], t_fin, t_prin)
    rows = NamedTuple[]
    for s in seed, t in tau, zr in zetarho, r in rho0
        push!(rows, (fn=0, L=float(L), rho0=float(r), kd=1/float(t),
                     rhocr=COMMON.rhocr, ap=COMMON.ap, nu1=COMMON.nu1,
                     gamma=COMMON.gamma, kp=COMMON.kp, zetap=COMMON.zetap,
                     zetap2=COMMON.zetap2, zetarho=float(zr), ar=ar_of(zr),
                     xi=COMMON.xi, M=COMMON.M, seed=float(s),
                     dx=COMMON.dx, dz=COMMON.dz, dtmin=COMMON.dtmin,
                     t_fin=t_fin, t_prin=t_prin, t_check=T_CHECK))
    end
    df = DataFrame(rows)
    df.fn = 1:nrow(df)
    return df
end

write_df(n, df, note) = begin
    CSV.write(joinpath(OUT, "DF_$(n).csv"), df)
    println(rpad("  DF_$(n).csv", 16), rpad(nrow(df), 6), note)
end

println("figure        rows  parameters")

# Fig1 - fig:schema. rho0 = 0.45, 0.6, 0.75; tau = 5; zeta_rho = 4; L = 50.
write_df(1, sweep(L=50, rho0=[0.45,0.6,0.75], tau=[5], zetarho=[4],
                  t_fin=HORIZON[1][1], t_prin=HORIZON[1][2]),
         "L=50, rho0 in {0.45,0.6,0.75}, tau=5, zeta=4")

# Fig4 - fig:Ndefects. (rho0, zeta_rho) plane at six renewal times.
# The six are make_heatmap's tkd = [10,5,1,0.2,0.1,0.01]; 12*14*6 = 1008 rows,
# which is exactly the row count of figures/Fig4/Ndef.csv.
write_df(4, sweep(L=10, rho0=0.4:0.1:1.5, tau=[0.1,0.2,1,5,10,100],
                  zetarho=[0,1,2,4,6,8,10,12,14,16,18,20,22,24],
                  t_fin=HORIZON[4][1], t_prin=HORIZON[4][2]),
         "L=10, 12 rho0 x 14 zeta x 6 tau")

# Fig8 / Fig12 - fig:lattice1 and fig:lattice1extended, same two solutions.
lattice = sweep(L=10, rho0=[0.7], tau=[1], zetarho=[10],
                t_fin=HORIZON[8][1], t_prin=HORIZON[8][2])
append!(lattice, sweep(L=10, rho0=[1.3], tau=[0.2], zetarho=[1],
                       t_fin=HORIZON[8][1], t_prin=HORIZON[8][2]))
lattice.fn = 1:nrow(lattice)
write_df(8,  copy(lattice), "L=10, (0.7,tau=1,zeta=10) and (1.3,tau=0.2,zeta=1)")
write_df(12, copy(lattice), "same two solutions as Fig8")

# Fig9 - fig:gamma_rho. The full sweep; matches Analysis/PDpaper/DF.csv.
write_df(9, sweep(L=10, rho0=0.4:0.1:1.5,
                  tau=[0.1,0.2,0.5,1,1.25,1/0.6,2.5,5,10,100],
                  zetarho=[0,1,2,4,6,8,10,12,14,16,18,20,22,24],
                  t_fin=HORIZON[9][1], t_prin=HORIZON[9][2]),
         "L=10, 12 rho0 x 14 zeta x 10 tau (1680)")

# Fig11 - fig:statesL50. The L=50 series; this is 2D/DF.csv.
write_df(11, sweep(L=50, rho0=0.4:0.05:1.0, tau=[5], zetarho=[4],
                  t_fin=HORIZON[11][1], t_prin=HORIZON[11][2]),
         "L=50, rho0 = 0.40:0.05:1.00, tau=5, zeta=4")

println("\nNot generated - see README.md:")
println("  Fig2, Fig6, Fig7, Fig10, Fig13  (parameters not fully given in the article)")
println("  Fig3  (1D model, parameters are notebook inputs, not a DF row)")
println("  Fig5  (linear stability analysis; no simulation)")
