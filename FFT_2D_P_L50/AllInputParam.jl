# =============================================================================
#  Builds the parameter table DF.csv: one row per simulation.
#
#  Run once, before submitting the jobs:
#      julia --project=. FFT_2D_P_L50/AllInputParam.jl
#
#  Each `t...` vector below lists the values to sweep for one parameter;
#  `generate_dataframe` (Utilities/JulUtils.jl) takes their full factorial
#  product, so the number of simulations is the product of the lengths.
#
#  For this run set only rho0 is swept — 13 values from 0.4 to 1.0 — at fixed
#  turnover kd = 0.2 and isotropic activity zetarho = 4. The same script drives
#  the other run sets of the study with different vectors here.
#
#  NOTE: this overwrites DF.csv. The committed DF.csv is the table used for the
#  paper; running this reproduces it byte-for-byte.
# =============================================================================

include("../Utilities/using.jl")
using_pkg("DelimitedFiles, CSV, DataFrames, JLD, Random")
using_mod(".JulUtils")

dir = @__DIR__

# --- Geometry and time step --------------------------------------------------
# System size (square): L=50 at dx=1e-2 gives a 5000² grid, padded to 5008².
tL = [50.0]
tdx = [1e-2]
tdz = [1e-2]
tdt = [1e-2]                # upper bound on the adaptive step (`dtmin`)

# --- Density and turnover ----------------------------------------------------
tρ0 = Array(0.4:0.05:1.0)   # SWEPT: homeostatic density, 13 values
tρcr = [0.0]                # unused in this run set
tkd = [0.2]                 # turnover rate

## Polar
tap = [0.1]                 # polar ordering coefficient
# P - v, h, ΔμP
tν1 = [0]                   # flow alignment: off for the published runs
tγ = [1.0]                  # rotational mobility
tkp = [1.0e-4]              # Frank constant

## Stress
# Active Polar/Nematic ζpP_zP_z, ζqQ
tζp = [0]                   # anisotropic active stress: off
tζp2 = [0.0]
tζρ = [4]                   # isotropic active stress ∝ ρ³

## Rho
tar = [4.0/3.0]             # rescaled below by |ζρ|
## Friction
txi = [1.0]

## Diffusion
tM = [1e-4] #, 1e-3] (D = M*ar/ar = M)

# Seed for noise
tseed = Array(1:1)          # one realisation per parameter set

listname = ["L", "rho0", "kd", "rhocr", "ap", "nu1", "gamma", "kp", "zetap", "zetap2", "zetarho", "ar", "xi", "M", "seed", "dx", "dz", "dtmin"]
listtab = [tL, tρ0, tkd, tρcr, tap, tν1, tγ, tkp, tζp, tζp2, tζρ, tar, txi, tM, tseed, tdx, tdz, tdt]

df = generate_dataframe(listname, listtab; fn="NO");

# Post-process the factorial grid:
#  - tie the passive pressure coefficient `ar` to the activity, ar = |ζρ|·4/3,
#    so the two stay comparable across a sweep in ζρ (InputParameters.jl
#    recomputes the same value, so the column is informational);
#  - drop extensile activity when there is no turnover, a combination with no
#    steady state. Irrelevant here (kd = 0.2 > 0) but kept as it ran.
for i=1:nrow(df)
    kd = df[i, :kd]
    zr = df[i, :zetarho]
    df[i, :ar] = zr == 0.0 ? 4.0/3.0 : abs(zr)*df[i, :ar]
    df[i, :zetarho] = (zr>0) & (df[i,:kd]==0) ? 0.0 : zr
end
unique!(df)
# `fn` is the simulation index: row n is run by SLURM_ARRAY_TASK_ID = n.
insertcols!(df, 1, :fn => 1:nrow(df))
CSV.write(joinpath(dir,"DF.csv"), df)

println("Number of sims: ", nrow(df))
