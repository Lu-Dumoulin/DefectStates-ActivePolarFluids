# =============================================================================
#  Shared configuration for the analysis scripts.
#
#  Included (directly or indirectly) by every make_fig_*.jl. Set:
#
#      DATA_DIR   directory holding DF.csv and one <idx>/Data/ folder per run
#      FIG_DIR    where figures are written (defaults to DATA_DIR)
#
#  On the original machine these were hard-coded Windows drive paths, one per
#  run set (F:/ZetaP_v2/, Z:/2defects/, D:/PDpaper/, ...); pass the run set you
#  want through DATA_DIR instead.
# =============================================================================

# Resolve the run set first: loading the plotting and image stack takes a while,
# and there is no point paying for it only to fail on a missing path.
dir_df = get(ENV, "DATA_DIR", "")
isempty(dir_df) && error("""
    DATA_DIR is not set.

    Point it at the run set you want to analyse - the directory holding DF.csv
    and one <idx>/Data/ folder per simulation:

        DATA_DIR=/path/to/runset/ julia --project=. Analysis/<script>.jl

    Optionally set FIG_DIR for the output (it defaults to DATA_DIR).
    """)
dir_fig = get(ENV, "FIG_DIR", dir_df)
endswith(dir_df, "/") || (dir_df *= "/")
endswith(dir_fig, "/") || (dir_fig *= "/")
isfile(dir_df*"DF.csv") || error("No DF.csv in $dir_df - is DATA_DIR pointing at a run set?")

include("../Utilities/using.jl")
using_pkg("CairoMakie, JLD, Printf, LaTeXStrings, FixedPointNumbers, DelimitedFiles, CSV, DataFrames, FileIO, Base64, Colors, LinearAlgebra, Statistics, Images, ImageSegmentation, FFTW, DelaunayTriangulation")
using_mod(".PictUtils")
CairoMakie.activate!(type = "png")

empty_theme = Theme(
    Axis = (
        backgroundcolor = :transparent,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xticklabelsvisible = false, 
        yticklabelsvisible = false,
        xgridcolor = :transparent,
        ygridcolor = :transparent,
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0.0,0.0),
        yautolimitmargin = (0.0,0.0),
    )
)

df = CSV.read(dir_df*"DF.csv", DataFrame)
N_sim = nrow(df)

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
stepp = 15

G_avg = [0.475, 0.472, 0.461, 0.498, 0.522]
