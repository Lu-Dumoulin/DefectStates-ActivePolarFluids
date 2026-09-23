# =============================================================================
#  CUDA kernels for the compressible active polar fluid
#
#  Three kernels are launched per time step from 2D.jl:
#
#    kernel_comp_diff_factor   once, at setup: the spectral derivative factors
#    kernel_comp_μhσ!          chemical potential μ, molecular field h, stress σ
#    kernel_solve_kspace!      force balance -> velocity, in Fourier space
#    kernel_comp_P!            explicit Euler update of the polarity P
#
#  Indexing convention throughout: x is the first array axis, z the second.
#  Each thread owns one grid point (i,j); the launch geometry is set in
#  InputParameters.jl (16x16 blocks). Real-space kernels cover the full Nx*Nz
#  grid, which is padded up to a multiple of 16, so no bounds check is needed.
#  Fourier-space kernels run on the rfft half-grid and must guard i > Lkx.
#
#  See arXiv:2506.03795 for the derivation of these expressions.
# =============================================================================

include("InputParameters.jl")

# Precompute the spectral multipliers used for every derivative:
#   ∂x -> i*kx      ∂z -> i*kz      Δ -> -(kx² + kz²)
# Storing them once avoids recomputing k-vectors inside the time loop.
function kernel_comp_diff_factor(factordx, factordz, factorΔ, kx, kz, kx2, kz2, Lkx)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if i>Lkx
        return nothing
    end
    factordx[i,j] = im*kx[i]
    factordz[i,j] = im*kz[j]
    factorΔ[i,j] = -kx2[i]-kz2[j]
    return nothing
end

# -----------------------------------------------------------------------------
# Chemical potential μ = δF/δρ, molecular field h = -δF/δP, and stress σ.
#
# The free-energy density is
#     f = (1/4) ar ρ⁴ - (1/2) ap ρ³ P²/ρ0 + (1/4) ap ρ² (P²)² + (1/2) kp ρ² |∇P|²
# i.e. a compressible (ρ⁴) part, a density-controlled polar ordering part whose
# sign flips at ρ = ρ0, and a one-constant Frank energy of stiffness kp.
#
# Inputs Dρ and DP hold the spectral derivatives computed in 2D.jl:
#     Dρ[:,:,1:2]  = ∂xρ, ∂zρ
#     DP[:,:,1:3]  = ∂xPx, ∂zPx, ΔPx
#     DP[:,:,4:6]  = ∂xPz, ∂zPz, ΔPz
#
# Output σ is stored as the four components (xx, xz, zx, zz) in slices 1..4;
# note xz and zx differ because the Ericksen stress has an antisymmetric part.
#
# A variant of the Frank energy that distinguishes splay from bend (extra
# 2∂xPx∂zPz - 2∂xPz∂zPx cross terms in μ, h, f and σ) was explored during the
# study; the published runs use the one-constant form implemented here.
# -----------------------------------------------------------------------------
function kernel_comp_μhσ!(σ, μ, h, ρ_, Dρ, P, DP, ar, ap, kp, ζρ, ζp, ζp2, Δμ, ν1, ν2, ρ0)
    # Indexing (i,j)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    @inbounds begin
        ρ = ρ_[i,j]
        Px = P[i,j,1]; Pz = P[i,j,2]
        ∂xρ = Dρ[i,j,1]; ∂zρ = Dρ[i,j,2]
        ∂xPx = DP[i,j,1]; ∂zPx = DP[i,j,2]; ΔPx = DP[i,j,3]
        ∂xPz = DP[i,j,4]; ∂zPz = DP[i,j,5]; ΔPz = DP[i,j,6]

        P2 = Px*Px + Pz*Pz

        # μ = δf/δρ
        μ[i,j] = ar*ρ^3 - 3*ρ^2*0.5*ap*P2/ρ0 + 0.5*ρ*ap*P2*P2 + ρ*kp*(∂xPx^2+∂zPx^2+∂xPz^2+∂zPz^2)

        # h = -δf/δP: ordering term (changes sign at ρ = ρ0), saturation, Frank
        # relaxation, plus the ∇ρ·∇P terms from the ρ-dependent Frank stiffness.
        hx = ρ^2*(ap*ρ*Px/ρ0 - ap*Px*P2 + kp*ΔPx) + 2*kp*ρ*(∂xρ*∂xPx+∂zρ*∂zPx); h[i,j,1] = hx
        hz = ρ^2*(ap*ρ*Pz/ρ0 - ap*Pz*P2 + kp*ΔPz) + 2*kp*ρ*(∂xρ*∂xPz+∂zρ*∂zPz); h[i,j,2] = hz

        # f - ρμ: the isotropic (pressure-like) part common to σxx and σzz
        f_ρμ = 0.25*ar*ρ^4 + ρ^2*( -0.5*ap*ρ/ρ0*P2 + 0.25*ap*P2*P2 + 0.5*kp*(∂xPx^2+∂zPx^2+∂xPz^2+∂zPz^2) ) - ρ*μ[i,j]

        # Symmetric stress: Ericksen + flow alignment (ν1, ν2) + active terms.
        # ζρ Δμ ρ³ is the isotropic active stress, ζp/ζp2 the anisotropic ones.
        σ[i,j,1] = f_ρμ - ρ^2*kp*(∂xPx^2+∂xPz^2) + ν1*Px*hx + ν2*(Px*hx+Pz*hz) + ζp*Δμ*Px*Px*ρ + ζp2*ρ*Δμ*(Px*Px-Pz*Pz) + ζρ*Δμ*ρ^3
        σ[i,j,4] = f_ρμ - ρ^2*kp*(∂zPx^2+∂zPz^2) + ν1*Pz*hz + ν2*(Px*hx+Pz*hz) + ζp*Δμ*Pz*Pz*ρ - ζp2*ρ*Δμ*(Px*Px-Pz*Pz) + ζρ*Δμ*ρ^3
        σ[i,j,2] = - ρ^2*kp*(∂xPx*∂zPx+∂xPz*∂zPz) + 0.5*ν1*(Px*hz+Pz*hx) + ζp*Δμ*Px*Pz*ρ

        # Antisymmetric (torque) part of the Ericksen stress: σxz = s + a, σzx = s - a
        σ_anti = 0.5*(Px*hz-hx*Pz)
        σ[i,j,3] = σ[i,j,2] - σ_anti
        σ[i,j,2] += σ_anti

    end
    return nothing
end

# -----------------------------------------------------------------------------
# Force balance, solved pointwise in Fourier space.
#
# Momentum conservation for a viscous fluid on a substrate with friction ξ:
#
#     η(2∂ₓ² + ∂_z²)vx + η ∂ₓ∂_z vz - ξ vx + ∂ⱼσ_xⱼ = 0        (and x <-> z)
#
# The viscous stress is not carried in σ; it is applied here analytically with
# the viscosity nondimensionalised to η = 1. In Fourier space this is a 2x2
# linear system per mode,
#
#     (-2kx² - kz² - ξ) ṽx  -  kx kz ṽz  =  -i(kx σ̃xx + kz σ̃xz)
#     (-2kz² - kx² - ξ) ṽz  -  kx kz ṽx  =  -i(kx σ̃zx + kz σ̃zz)
#
# solved below by eliminating ṽx (hence `fact`, the inverse of the xx operator)
# and back-substituting. Being local in k, this replaces a global Stokes solve
# with one flop-cheap kernel per step — the reason the scheme is spectral.
# -----------------------------------------------------------------------------
function kernel_solve_kspace!(fv, fσ, kx_, kz_, Lkx, ξ)  # Indexing (i,j)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if i>Lkx
        return nothing
    end
    kx = kx_[i];     kx2 = kx*kx
    kz = kz_[j];     kz2 = kz*kz
    kxz = kx*kz
    fσxx = fσ[i,j,1]; fσzz = fσ[i,j,4]
    fσxz = fσ[i,j,2]; fσzx = fσ[i,j,3]

    fact = 1.0/(-2*kx2-kz2-ξ)
    fvzij = im / (2*kz2 + kxz*kxz*fact + kx2 + ξ) * (kz*fσzz + kx*fσzx + kxz*(kx*fσxx+kz*fσxz)*fact)
    fv[i,j,1] = (-im*(kx*fσxx+kz*fσxz)+kxz*fvzij)*fact

    fv[i,j,2] = fvzij
    return nothing
end

# -----------------------------------------------------------------------------
# Explicit Euler step for the polarity:
#
#     ∂ₜP = -(v·∇)P - ω·P - ν1 P·u - ν2 P(∇·v) + γ h + ε Δμ ρ P
#
# with u the symmetric and ω the antisymmetric part of ∇v. The published runs
# use ν1 = ν2 = ε = 0, so only advection, co-rotation and relaxation along h
# are active; the terms are kept for completeness.
#
# ∂v holds the velocity gradients assembled in 2D.jl:
#     ∂v[:,:,1] = ∂xvx   ∂v[:,:,2] = ∂xvz   ∂v[:,:,3] = ∂zvx   ∂v[:,:,4] = ∂zvz
# -----------------------------------------------------------------------------
function kernel_comp_P!(P, h, DP, v, ∂v, ρ_, ν1, ν2, γ, ε, Δμ, Δt)
   # Indexing (i,j)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    vx = v[i,j,1]; vz = v[i,j,2]
    vxx = ∂v[i,j,1]                     # ∂xvx
    vxz = 0.5*(∂v[i,j,2]+∂v[i,j,3])     # symmetric shear rate u_xz
    vzz = ∂v[i,j,4]                     # ∂zvz
    ωxz = 0.5*(∂v[i,j,2]-∂v[i,j,3])     # vorticity
    Px = P[i,j,1]; Pz = P[i,j,2]
    hx = h[i,j,1]; hz = h[i,j,2]

    ρ = ρ_[i,j]

    P[i,j,1] = Px + Δt*(-vx*DP[i,j,1] - vz*DP[i,j,2] - ωxz*Pz - ν1*(Px*vxx+Pz*vxz) - ν2*Px*(vxx+vzz) + γ*hx + ε*Δμ*ρ*Px)
    P[i,j,2] = Pz + Δt*(-vx*DP[i,j,4] - vz*DP[i,j,5] + ωxz*Px - ν1*(Px*vxz+Pz*vzz) - ν2*Pz*(vxx+vzz) + γ*hz + ε*Δμ*ρ*Pz)

    return nothing
end
