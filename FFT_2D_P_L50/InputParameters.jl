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
# dir = "/srv/beegfs/scratch/users/d/dumoulil/Data/P-series/Dt/" 
# idx = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])


include("../Utilities/using.jl")
using_pkg("FFTW, Distributions, DelimitedFiles, CSV, DataFrames, Dates, Printf, JLD, CUDA, Random")

dir_df = @__DIR__
df = CSV.read(joinpath(dir_df,"DF.csv"), DataFrame)[idx,:]

# tfiles_names = isdir(string(file,"Data/")) ? readdir(string(file,"Data/")) : []
# num_files = length(tfiles_names)
# num_files < 4 ? nothing : return 0

Tf = Float64

# Δx = 1e-2; Δz = 1e-2
const Δx::Float64 = Tf(df[:dx]); const Δz::Float64 = Tf(df[:dz])
const Δx2::Float64 = Δx*Δx; const Δz2::Float64 = Δz*Δz
Δtmin::Float64 = Tf(df[:dtmin])
global Δt::Float64 = Δtmin
@show Δx2/Δt*0.5
const L = Tf(df[:L])
N = L / Δx
Nline = N 

WrapsT = 16
Bx = ceil(Int, N/WrapsT)
Bz = ceil(Int, Nline/WrapsT)
block_dim = (WrapsT, WrapsT)
grid_dim = (Bx, Bz)
gridFFT_dim = (div(Bx,2)+1, Bz)

const Nx::Int = WrapsT * Bx
const Nz::Int = WrapsT * Bz

### FFT
const W = plan_rfft(CUDA.ones(Float64, Nx, Nz))
const Wi = inv(W)
const kx = CuArray{Float64}(2*pi*rfftfreq(Nx, 1/Δx))
const kz = CuArray{Float64}(2*pi*fftfreq(Nz, 1/Δz))
const kx2 = kx.*kx
const kz2 = kz.*kz
const Lkx::Int = length(kx)

# Non-dimensionalized parameters
const ρ0::Float64 = Tf(df[:rho0])
ρcr = Tf(df[:rhocr])
const kd::Float64 = Tf(df[:kd])
kpoly = kd*ρ0

## Polar
const ap::Float64 = Tf(df[:ap])
# P - v, h, ΔμP
const ν1::Float64 = Tf(df[:nu1]); 
const ν2::Float64 = Tf(0.0); 
const γ::Float64 = Tf(df[:gamma])
const ε::Float64 = Tf(0.0);
# h - ap\rhoP-apP^3+k_pΔP, β1 normal + β2 tangential anchoring 
const kp::Float64 = Tf(df[:kp])

## Stress
# Active Polar/Nematic ζP_zP_z, ζ1Q
const ζp::Float64 = Tf(df[:zetap])
const ζp2::Float64 = Tf(df[:zetap2])
const Δμ::Float64 = Tf(1.0)
# Active rho
const ζρ::Float64 = Tf(df[:zetarho])
# Friction
const ξ::Float64 = Tf(df[:xi])#1.0)
## Rho
const ar::Float64 = ζρ==0 ? 4/3 : abs(ζρ)*4/3 #Tf(df[:ar])
const M::Float64 = Tf(df[:M])/ar#*10.0
# seed
const sd = Int(df[:seed])

# d = Tf(df[:D])

## Debug
# t_fin = 40000
# t_prin = 2000
# t_check = 1
## 2 defects ?
# t_fin = 30000
# t_prin = 500
# t_check = 1
# # 7 days
# t_fin = 2000000
# t_prin = 10000
# t_check = 1
# 12h
t_fin = 150000
t_prin = 1000
t_check = 1
set_zero_subnormals(true)
