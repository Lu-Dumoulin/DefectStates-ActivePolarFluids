include("InputParameters.jl")

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

# function kernel_ini_P!(P, Nx, Ny, d)
#     # Indexing (i,j)
#     i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
#     @inbounds begin
#         # ic = floor(Int, Ny*0.5)i>floor(Int, Nx*0.5) ? floor(Int, Nx*0.75) : floor(Int, Nx*0.25)
#         c_def = i>floor(Int, Nx*0.5) ? -1 : 1
#         ic = c_def == -1 ? floor(Int, Nx*(0.5+d/2)) : floor(Int, Nx*(0.5-d/2))
#         jc = floor(Int, Ny*0.5)
#         θ = atan(j-jc,i-ic)
#         P[i,j,1] =  c_def == 1 ? cos(θ) : -cos(θ)
#         P[i,j,2] = sin(θ)
#     end
#     return nothing
# end

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
        
        μ[i,j] = ar*ρ^3 - 3*ρ^2*0.5*ap*P2/ρ0 + 0.5*ρ*ap*P2*P2 + ρ*kp*(∂xPx^2+∂zPx^2+∂xPz^2+∂zPz^2)
        # μ[i,j] = ar*ρ^3 - 3*ρ^2*0.5*ap*P2/ρ0 + 0.5*ρ*ap*P2*P2 + ρ*kp*(∂xPx^2+∂zPx^2+∂xPz^2+∂zPz^2+2*∂xPx*∂zPz-2*∂xPz*∂zPx)
        hx = ρ^2*(ap*ρ*Px/ρ0 - ap*Px*P2 + kp*ΔPx) + 2*kp*ρ*(∂xρ*∂xPx+∂zρ*∂zPx); h[i,j,1] = hx 
        # hx = ρ^2*(ap*ρ*Px/ρ0 - ap*Px*P2 + kp*ΔPx) + 2*kp*ρ*(∂xρ*(∂xPx+∂zPz)+∂zρ*(∂zPx-∂xPz)); h[i,j,1] = hx 
        hz = ρ^2*(ap*ρ*Pz/ρ0 - ap*Pz*P2 + kp*ΔPz) + 2*kp*ρ*(∂xρ*∂xPz+∂zρ*∂zPz); h[i,j,2] = hz 
        # hz = ρ^2*(ap*ρ*Pz/ρ0 - ap*Pz*P2 + kp*ΔPz) + 2*kp*ρ*(∂xρ*(∂xPz-∂zPx)+∂zρ*(∂zPz+∂xPx)); h[i,j,2] = hz 
        # f - ρμ
        f_ρμ = 0.25*ar*ρ^4 + ρ^2*( -0.5*ap*ρ/ρ0*P2 + 0.25*ap*P2*P2 + 0.5*kp*(∂xPx^2+∂zPx^2+∂xPz^2+∂zPz^2) ) - ρ*μ[i,j]
        # f_ρμ = 0.25*ar*ρ^4 + ρ^2*( -0.5*ap*ρ/ρ0*P2 + 0.25*ap*P2*P2 + 0.5*kp*(∂xPx^2+∂zPx^2+∂xPz^2+∂zPz^2+2*∂xPx*∂zPz-2*∂xPz*∂zPx) ) - ρ*μ[i,j]
        
        # σ = sym Ericksen + flow aligment + actif + hydro pressur in bρ^4
        σ[i,j,1] = f_ρμ - ρ^2*kp*(∂xPx^2+∂xPz^2) + ν1*Px*hx + ν2*(Px*hx+Pz*hz) + ζp*Δμ*Px*Px*ρ + ζp2*ρ*Δμ*(Px*Px-Pz*Pz) + ζρ*Δμ*ρ^3
        # σ[i,j,1] = f_ρμ - ρ^2*kp*(∂xPx^2+∂xPz^2+2*∂xPx*∂zPz-2*∂xPz*∂zPx) + ν1*Px*hx + ν2*(Px*hx+Pz*hz) + ζp*Δμ*Px*Px*ρ + ζp2*ρ*Δμ*(Px*Px-Pz*Pz) + ζρ*Δμ*ρ^3
        σ[i,j,4] = f_ρμ - ρ^2*kp*(∂zPx^2+∂zPz^2) + ν1*Pz*hz + ν2*(Px*hx+Pz*hz) + ζp*Δμ*Pz*Pz*ρ - ζp2*ρ*Δμ*(Px*Px-Pz*Pz) + ζρ*Δμ*ρ^3
        # σ[i,j,4] = f_ρμ - ρ^2*kp*(∂zPx^2+∂zPz^2+2*∂xPx*∂zPz-2*∂xPz*∂zPx) + ν1*Pz*hz + ν2*(Px*hx+Pz*hz) + ζp*Δμ*Pz*Pz*ρ - ζp2*ρ*Δμ*(Px*Px-Pz*Pz) + ζρ*Δμ*ρ^3
        σ[i,j,2] = - ρ^2*kp*(∂xPx*∂zPx+∂xPz*∂zPz) + 0.5*ν1*(Px*hz+Pz*hx) + ζp*Δμ*Px*Pz*ρ
        
        # σ += anti-sym part of Ericksen
        σ_anti = 0.5*(Px*hz-hx*Pz)
        σ[i,j,3] = σ[i,j,2] - σ_anti
        σ[i,j,2] += σ_anti
        
    end
    return nothing
end

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

function kernel_comp_P!(P, h, DP, v, ∂v, ρ_, ν1, ν2, γ, ε, Δμ, Δt)
   # Indexing (i,j)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    vx = v[i,j,1]; vz = v[i,j,2]
    vxx = ∂v[i,j,1]
    vxz = 0.5*(∂v[i,j,2]+∂v[i,j,3])
    vzz = ∂v[i,j,4]
    ωxz = 0.5*(∂v[i,j,2]-∂v[i,j,3])
    Px = P[i,j,1]; Pz = P[i,j,2]
    hx = h[i,j,1]; hz = h[i,j,2]

    ρ = ρ_[i,j]

    P[i,j,1] = Px + Δt*(-vx*DP[i,j,1] - vz*DP[i,j,2] - ωxz*Pz - ν1*(Px*vxx+Pz*vxz) - ν2*Px*(vxx+vzz) + γ*hx + ε*Δμ*ρ*Px)
    P[i,j,2] = Pz + Δt*(-vx*DP[i,j,4] - vz*DP[i,j,5] + ωxz*Px - ν1*(Px*vxz+Pz*vzz) - ν2*Pz*(vxx+vzz) + γ*hz + ε*Δμ*ρ*Pz)

    return nothing
end