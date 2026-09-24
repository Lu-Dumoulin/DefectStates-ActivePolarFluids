# =============================================================================
#  Figure 7 - fig:phasediagram
#
#  Classification of states over the renewal parameters tau and rho0, coloured
#  by defect density, structural persistence time, and the mean and standard
#  deviation of low-density areas.
#
#      FIG_DIR=<repo>/figures/Fig7/ julia --project=. Analysis/make_fig7.jl
#
#  Inputs : Analysis/PDpaper/, the phase-diagram tables, which ship with the
#           repository - so this runs without simulation output. DATA_DIR
#           overrides them to rebuild from your own snapshots.
#  Outputs: DF_tikz_norm_adjusted.csv, the table the figure reads.
#
#  The runs behind those tables used t_fin = 300000, t_prin = 1000
#  (Code/FFT_2D_P_phasediag2 on the cluster). No DF_7.csv yet: the article
#  gives neither the tau and rho0 grids nor zeta_rho - see params/README.md.
# =============================================================================

include("PhaseDiagram.jl")

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    df_tikz_phase_diagram()
end
