# =============================================================================
#  Figure 8 - fig:lattice1
#
#  Density and polarity-angle maps for a square and a hexagonal defect lattice.
#
#      DATA_DIR=/path/to/fig8-runs/ FIG_DIR=<repo>/figures/Fig8/ \
#          julia --project=. Analysis/make_fig8.jl
#
#  Inputs : params/DF_8.csv - two runs, (rho0=0.7, tau=1, zeta=10) and
#           (rho0=1.3, tau=0.2, zeta=1), both to t_fin = 2e6 snapshotting
#           every 10000.
#  Outputs: density.png / angle.png for the first solution and
#           density-2.png / angle-2.png for the second.
#
#  plot_dens_and_angle writes into its own subdirectory per index; the figure
#  expects the four files flat, so they are collected here.
# =============================================================================

include("MakePlots.jl")

"""
    collect_panels(idx, suffix)

Run plot_dens_and_angle for one simulation and move its density/angle panels
to the flat names the figure reads.
"""
function collect_panels(idx, suffix)
    plot_dens_and_angle(idx)
    src = joinpath(dir_fig, "density_and_angle_tri_$(idx)/")
    for name in ("density", "angle")
        from = joinpath(src, name*".png")
        isfile(from) || (@warn "not produced: $from"; continue)
        mv(from, joinpath(dir_fig, name*suffix*".png"); force=true)
    end
end

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    collect_panels(1, "")     # square lattice
    collect_panels(2, "-2")   # hexagonal lattice
end
