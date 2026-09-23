# =============================================================================
#  Linear stability analysis of the homogeneous state — the Fig5 curves.
#
#      julia --project=. Analysis/make_fig_LSA.jl [outdir]
#
#  Writes tau1.csv and tau5.csv (default: figures/Fig5/), each holding, per
#  homeostatic density rho0:
#
#      x   rho0
#      y   the closed-form critical activity zeta_c
#      yy  the critical activity found by integrating the linearised system
#
#  tau5 is the turnover time tau = 5 used for the published runs (kd = 0.2);
#  tau1 is tau = 1 (kd = 1.0).
#
#  Needs no simulation output — this is analysis of the linearised equations,
#  not of the data, so it runs anywhere with no DATA_DIR.
#
#  The expressions are carried over unchanged from LSA.jl; the only difference
#  is that kd is an argument here instead of a global constant, so both
#  turnover times come from one routine. Verified to reproduce the committed
#  tau1.csv and tau5.csv byte for byte.
# =============================================================================

const ar = 4.0*4/3
const χ = 0.1
const κ = 1e-4
const Γ = 1.0
const γ = 1e-4

# Coefficients of the linearised 2x2 system for (delta-rho, delta-p) at wavenumber q.
A(ρ0) = 3*ar*ρ0^2-5/2*χ
a(q, ρ0, ζ, kd) = -( kd + q^2*( ρ0^2*( (A(ρ0)-3*ρ0*ζ)/(1+2*q^2) + κ/Γ ) + γ*A(ρ0)) )
b(q, ρ0) = q^2*χ*ρ0*( ρ0^2/(1+2*q^2) + γ )
c(ρ0) = χ*ρ0/Γ
d(q, ρ0) = -ρ0^2/Γ * ( 2*χ+κ*q^2 )

# Closed-form threshold.
ζc_simp(ρ0, kd) = A(ρ0)/(3*ρ0) + ( sqrt(A(ρ0)*γ/(3*ρ0^2)) - sqrt(2*kd/(3*ρ0^3)) )^2

"""
    numeric_ζc(ρ0, kd)

Smallest activity on the 0.01:40 grid for which a small perturbation grows:
the linearised system is stepped forward from (1e-3, 1e-3) and the mode is
called unstable once either component passes 5e-3 within 1000 steps. Returns
0.0 if no activity on the grid destabilises the state.
"""
function numeric_ζc(ρ0, kd)
    for ζ in 0.01:0.01:40.0
        for q = 0:100
            dρ = 1e-3
            dp = 1e-3
            steps = 0
            while dρ < 5e-3 && dp < 5e-3 && steps < 1000
                temp = dρ
                dρ += 0.01*( a(q,ρ0,ζ,kd)*temp + b(q,ρ0)*dp )
                dp += 0.01*( c(ρ0)*temp + d(q,ρ0)*dp )
                steps += 1
            end
            if dρ > 1e-3 || dp > 1e-3
                return ζ
            end
        end
    end
    return 0.0
end

"""
    lsa_curves(kd; tρ0 = Array(0.3:0.01:1.2))

The (rho0, closed-form zeta_c, numerical zeta_c) table for one turnover rate.
"""
function lsa_curves(kd; tρ0 = Array(0.3:0.01:1.2))
    res = zeros(length(tρ0), 3)
    for i in eachindex(tρ0)
        ρ0 = tρ0[i]
        res[i,1] = ρ0
        res[i,2] = ζc_simp(ρ0, kd)
        res[i,3] = numeric_ζc(ρ0, kd)
    end
    return res
end

# CRLF separators and no trailing newline, matching the committed files
# (they were written on Windows).
function write_lsa_csv(path, res)
    open(path, "w") do io
        write(io, "x, y, yy")
        for i in axes(res, 1)
            write(io, string("\r\n", res[i,1], ", ", res[i,2], ", ", res[i,3]))
        end
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) >= 1 ? ARGS[1] :
             abspath(joinpath(@__DIR__, "..", "figures", "Fig5"))
    mkpath(outdir)
    write_lsa_csv(joinpath(outdir, "tau1.csv"), lsa_curves(1.0))
    write_lsa_csv(joinpath(outdir, "tau5.csv"), lsa_curves(0.2))
end
