# =============================================================================
#  Per-simulation setup: reads one row of DF.csv and defines every constant
#  the kernels and the time loop use as a global.
#
#  Included by kernels.jl, which is included by 2D.jl. Nothing here is a
#  function: the file executes top to bottom and leaves its bindings in Main.
#  Most are `const` so the kernels specialise on concrete Float64 types.
# =============================================================================

# ==============================================================================
#  CONFIGURATION  --  the only block that differs from the code run for the paper
# ==============================================================================
#  On Baobab (UNIGE) these three lines were hard-coded to:
#
#      dir       = "/srv/beegfs/scratch/users/d/dumoulil/Data/P-series/L50/"
#      localpath = "Z:/L50/"
#      idx       = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
#
#  They are now read from the environment so the code runs anywhere. Everything
#  below this block, and every other simulation file, is byte-identical to what
#  produced the published results.
#
#      DATA_DIR  output root; simulation `idx` writes to <DATA_DIR>/<idx>/Data/
#                (default: ./data/L50/ next to this repository)
#      SIM_IDX   which row of DF.csv to run; defaults to SLURM_ARRAY_TASK_ID
#                under Slurm, or 1 when running a single job by hand
# ==============================================================================
dir = get(ENV, "DATA_DIR", abspath(joinpath(@__DIR__, "..", "data", "L50")) * "/")
localpath = get(ENV, "LOCAL_DATA_DIR", "Z:/L50/")
idx = Base.parse(Int, get(ENV, "SIM_IDX", get(ENV, "SLURM_ARRAY_TASK_ID", "1")))
@show fn = "$idx/"
file = joinpath(dir, fn)
mkpath(file)
# Added for this release: on the cluster the <idx>/Data/ folders already existed
# from earlier runs, so the original code never created them. 2D.jl saves into
# string(file, "Data/data<t>.jld"), which fails on a fresh checkout without this.
mkpath(joinpath(file, "Data"))
println("path_c = ", dir)
println("path_l = ", localpath)
println(idx)

# `using_pkg` installs anything missing on first use — the cluster runs had no
# Project.toml. With the Project.toml shipped here it just resolves normally.
include("../Utilities/using.jl")
using_pkg("FFTW, Distributions, DelimitedFiles, CSV, DataFrames, Dates, Printf, JLD, CUDA, Random")

# One row of the parameter table = one simulation. See README for the columns.
# DF_FILE selects a different table - the per-figure ones in params/ - and is
# resolved relative to this directory when given as a relative path.
dir_df = @__DIR__
df_file = get(ENV, "DF_FILE", "DF.csv")
df_path = isabspath(df_file) ? df_file : joinpath(dir_df, df_file)
isfile(df_path) || error("No parameter table at $df_path (set DF_FILE)")
df = CSV.read(df_path, DataFrame)[idx,:]

Tf = Float64

# --- Grid and time step ------------------------------------------------------
const Δx::Float64 = Tf(df[:dx]); const Δz::Float64 = Tf(df[:dz])
const Δx2::Float64 = Δx*Δx; const Δz2::Float64 = Δz*Δz
# Despite the name, Δtmin is the *upper* bound on the adaptive step set in
# 2D.jl; Δt starts there and is then limited by the CFL condition.
Δtmin::Float64 = Tf(df[:dtmin])
global Δt::Float64 = Δtmin
@show Δx2/Δt*0.5
const L = Tf(df[:L])
N = L / Δx
Nline = N

# The grid is padded up to a multiple of the 16x16 block size so that the
# real-space kernels need no bounds check. For the paper runs L=50, Δx=0.01
# gives N=5000 -> Bx=313 blocks -> Nx=Nz=5008.
WrapsT = 16
Bx = ceil(Int, N/WrapsT)
Bz = ceil(Int, Nline/WrapsT)
block_dim = (WrapsT, WrapsT)
grid_dim = (Bx, Bz)
# Fourier-space launches only need the rfft half-grid; this over-covers it
# slightly, which is why those kernels guard on `i > Lkx`.
gridFFT_dim = (div(Bx,2)+1, Bz)

const Nx::Int = WrapsT * Bx
const Nz::Int = WrapsT * Bz

# --- FFT plans and wavevectors ----------------------------------------------
# Real-to-complex along x (hence Lkx = Nx/2+1 modes), complex along z.
# Plans are created once and reused for every transform in the time loop.
const W = plan_rfft(CUDA.ones(Float64, Nx, Nz))
const Wi = inv(W)
const kx = CuArray{Float64}(2*pi*rfftfreq(Nx, 1/Δx))
const kz = CuArray{Float64}(2*pi*fftfreq(Nz, 1/Δz))
const kx2 = kx.*kx
const kz2 = kz.*kz
const Lkx::Int = length(kx)

# --- Physical parameters (non-dimensionalised) -------------------------------
const ρ0::Float64 = Tf(df[:rho0])   # homeostatic density; polar order sets in above it
ρcr = Tf(df[:rhocr])                # unused in this run set
const kd::Float64 = Tf(df[:kd])     # turnover rate
kpoly = kd*ρ0                       # unused in this run set

## Polar
const ap::Float64 = Tf(df[:ap])     # polar ordering coefficient
const ν1::Float64 = Tf(df[:nu1]);   # flow alignment (0 for the published runs)
const ν2::Float64 = Tf(0.0);        # alignment with compression
const γ::Float64 = Tf(df[:gamma])   # rotational mobility, ∂ₜP ⊃ γh
const ε::Float64 = Tf(0.0);         # active polarity source, ∂ₜP ⊃ εΔμρP
const kp::Float64 = Tf(df[:kp])     # Frank constant (one-constant approximation)

## Stress
const ζp::Float64 = Tf(df[:zetap])      # anisotropic active stress ∝ PᵢPⱼ
const ζp2::Float64 = Tf(df[:zetap2])    # anisotropic active stress, second form
const Δμ::Float64 = Tf(1.0)             # chemical drive; sets the activity scale
const ζρ::Float64 = Tf(df[:zetarho])    # isotropic active stress ∝ ρ³
const ξ::Float64 = Tf(df[:xi])          # substrate friction

## Rho
# ar is recomputed here rather than read from DF.csv (AllInputParam.jl writes
# the same value into the `ar` column); tying it to |ζρ| keeps the passive
# pressure comparable to the active stress across the sweep.
const ar::Float64 = ζρ==0 ? 4/3 : abs(ζρ)*4/3
const M::Float64 = Tf(df[:M])/ar        # mobility; D = M·ar so this fixes D
const sd = Int(df[:seed])               # RNG seed for the initial noise

# --- Integration horizon -----------------------------------------------------
# Taken from the parameter table when it provides them, since the horizon
# differs by figure (see params/README.md). Tables without these columns - the
# archival 2D/DF.csv among them - fall back to the values used for the paper's
# L=50 runs, sized for the 12 h wall-clock limit of the Slurm job.
t_fin   = hasproperty(df, :t_fin)   ? df[:t_fin]   : 150000
t_prin  = hasproperty(df, :t_prin)  ? df[:t_prin]  : 1000
t_check = hasproperty(df, :t_check) ? df[:t_check] : 1

# Flush denormals to zero: they are worthless here and slow on GPU.
set_zero_subnormals(true)
