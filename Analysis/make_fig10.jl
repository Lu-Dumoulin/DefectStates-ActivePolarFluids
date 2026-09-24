# =============================================================================
#  Figure 10 - fig:tauc
#
#  Density-density correlation against time for three solutions, with the
#  exponential fits that define the structural persistence time T_p.
#
#      FIG_DIR=<repo>/figures/Fig10/ julia --project=. Analysis/make_fig10.jl
#
#  Inputs : Analysis/PDpaper/, as for figure 7; runs without simulation output.
#  Outputs: DF_tikz_exp_<t>.csv and DF_tikz_exp_3_<t>.csv
#
#  The figure ships four of these, at t = 20 and t = 50, from both routines.
#  The caption's t_f = 287e3 sits inside the 300000 horizon those runs used.
#  No DF_10.csv yet: the caption gives rho0 and tau for the three solutions but
#  neither zeta_rho nor the system size - see params/README.md.
# =============================================================================

include("PhaseDiagram.jl")

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    for t in (20, 50)
        make_csv_exp_plot(t)
        make_csv_exp_plot_3(t)
    end
end
