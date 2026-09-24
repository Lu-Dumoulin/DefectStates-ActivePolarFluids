# =============================================================================
#  Figure 11 - fig:statesL50
#
#  Density snapshots at L=50 over the full range of target densities, the
#  larger-system counterpart of figure 6.
#
#      DATA_DIR=/path/to/fig11-runs/ FIG_DIR=<repo>/figures/Fig11/ \
#          julia --project=. Analysis/make_fig11.jl
#
#  Inputs : params/DF_11.csv - the L=50 series, rho0 = 0.40:0.05:1.00 at
#           tau = 5 and zeta_rho = 4. This is the same table as 2D/DF.csv.
#  Outputs: the density panel grid.
#
#  Drawn by the same routine as figure 1, over more of the series: the caption
#  notes that rho0 = 0.45 and 0.55 appear here but not in figure 6. The grid
#  shape below is a starting point - set Lz and Lx to the layout you want.
# =============================================================================

include("make_fig1.jl")

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    make_plot_L50(Array(1:9); Lz = 3, Lx = 3)
end
