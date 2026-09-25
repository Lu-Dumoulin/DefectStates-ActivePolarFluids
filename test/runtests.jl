# Test suite.
#
#     julia --project=. test/runtests.jl              # everything available
#     julia --project=. test/runtests.jl params       # one group
#
# What can be checked here is bounded by what this repository ships. There is
# no simulation output and the solver needs a GPU, so the kernels are not
# exercised; what is checked is that the parameter tables agree with the ones
# the runs used, that every source parses, that the figures build and are wired
# to the data they plot, and that the one analysis that needs no simulation
# output - the linear stability analysis - reproduces its committed data
# exactly.
#
# The hygiene group checks that nothing published here depends on, or
# discloses, the machine the runs were done on.
#
# Optional groups skip rather than fail when their tools are missing:
#   figures     needs pdflatex
#   submitted   needs pdflatex, ghostscript and SUBMITTED_DIR

using Test

include(joinpath(@__DIR__, "support", "helpers.jl"))

const GROUPS = ["hygiene", "params", "solver", "analysis", "figures", "submitted"]
selected = isempty(ARGS) ? GROUPS : ARGS

for g in selected
    g in GROUPS || error("unknown group $g; choose from $(join(GROUPS, ", "))")
end

@testset "DefectStates-ActivePolarFluids" begin
    for g in selected
        include(joinpath(@__DIR__, "test_$(g).jl"))
    end
end
