include("../Utilities/using.jl")
using_pkg("DelimitedFiles, CSV, DataFrames, JLD, Random")
using_mod(".JulUtils")

dir = @__DIR__

# System size (square)
tL = [50.0]#reverse([0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0])# [10.0]
tdx = [1e-2]
tdz = [1e-2]
tdt = [1e-2]#, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4]
# Non-dimensionalized parameters
tρ0 = Array(0.4:0.05:1.0)#Array(0.4:0.1:1.0)#[0.6, 0.7, 0.8]#, 0.7]#[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5] #0.4:0.2:1.2 #[0.8]#, 1.2]
# tρ0 = [0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 1.0, 1.2]#, 0.7]#[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5] #0.4:0.2:1.2 #[0.8]#, 1.2]
tρcr = [0.0]
# tkd = reverse([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 5.0])#reverse([0.0, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 10.0]) #[0.01, 0.05, 0.10, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]#, 1.0]
tkd = [0.2]#0.25, 0.5, 1.0] #reverse([0.0, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 10.0]) #[0.01, 0.05, 0.10, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]#, 1.0]

## Polar
tap = [0.1]
# P - v, h, ΔμP
tν1 = [0]#[-0.5, 0, 0.5]#, -0.5, 0.5, -1.1, 1.1]#, -0.01, 0.01, -5.0, 5.0] # [-1.1, 0, 1.1]#, 0.6, -0.6]
tγ = [1.0] # [1.0, 5.0] #[0.1, 1.0]#, 0.5, 1.0, 2.0]
tkp = [1.0e-4] #, 0.01]#, 1.0e-2, 1.0e-1, 1.0]

## Stress
# Active Polar/Nematic ζpP_zP_z, ζqQ
tζp = [0]#[-1, -0.5, 0.0, 0.5, 1.0]#Array(-0.015:0.001:0.015) #[-0.2, 0.0, 0.2]# [-1.0, -0.5, 0.0, 0.5, 1.0]#, 1.0]
tζp2 = [0.0] #[-1.0, 0.0, 1.0]#[-0.5, 0.0, 0.5]#, 1.0]
tζρ = [4]#[-8, -4, -1, 0, 1, 4.0, 8.0]#, 6.0, 8.0, 10.0] #[0.0,1.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0] #, 4.0, 8.0, 16.0, 32.0] #

## Rho
tar = [4.0/3.0]#, 8.0/3.0] # times a #.*4.0/3.0
## Friction
txi = [1.0] #, 1.0]#, 10.0, 100.0]

## Diffusion
tM = [1e-4] #, 1e-3] (D = M*ar/ar = M)

# Distance betwween defects
# tD = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5].*0.5#Array(0.02:0.02:0.2)

# Seed for noise
tseed = Array(1:1)#50)

listname = ["L", "rho0", "kd", "rhocr", "ap", "nu1", "gamma", "kp", "zetap", "zetap2", "zetarho", "ar", "xi", "M", "seed", "dx", "dz", "dtmin"]#, "D"]
listtab = [tL, tρ0, tkd, tρcr, tap, tν1, tγ, tkp, tζp, tζp2, tζρ, tar, txi, tM, tseed, tdx, tdz, tdt]#, tD];

df = generate_dataframe(listname, listtab; fn="NO");
# df = DataAPI.allcombinations(DataFrame,listname, listtab) #generate_dataframe(listname, listtab; fn="NO");
for i=1:nrow(df)
    kd = df[i, :kd]
    zr = df[i, :zetarho]#*kd
    # df[i, :zetap] = df[i, :zetap] * zr
    df[i, :ar] = zr == 0.0 ? 4.0/3.0 : abs(zr)*df[i, :ar]
    df[i, :zetarho] = (zr>0) & (df[i,:kd]==0) ? 0.0 : zr
end
unique!(df)
insertcols!(df, 1, :fn => 1:nrow(df))
CSV.write(joinpath(dir,"DF.csv"), df)

println("Number of sims: ", nrow(df))
