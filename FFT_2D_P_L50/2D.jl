# =============================================================================
#  Entry point for one simulation.
#
#  Run as:  SIM_IDX=<row of DF.csv> DATA_DIR=<output root> julia --optimize=3 2D.jl
#  Under Slurm, SIM_IDX defaults to SLURM_ARRAY_TASK_ID (see slurm/submit.sh).
#
#  Includes kernels.jl, which in turn includes InputParameters.jl, so every
#  physical parameter and grid constant is already in scope as a global here.
#
#  Scheme (pseudo-spectral, explicit Euler in time). Per step:
#    1. FFT ρ and P, multiply by the spectral factors -> ∇ρ, ∇P, ΔP
#    2. kernel_comp_μhσ!    -> μ, h, σ in real space
#    3. FFT σ, kernel_solve_kspace! -> ṽ, then inverse FFT -> v and ∇v
#    4. kernel_comp_P!      -> advance P
#    5. advance ρ by the continuity equation with diffusion and turnover
#
#  All spatial derivatives are spectral, so the periodic box is exact to
#  round-off; only the time integration is approximate.
#
#  L. Dumoulin, C. Blanch-Mercader, K. Kruse, arXiv:2506.03795
# =============================================================================

include("kernels.jl")

function main()
    @inbounds begin
        # Initial noise is drawn on the CPU from a seeded RNG and copied to the
        # device, so a run is reproducible from `seed` in DF.csv regardless of
        # the GPU's own RNG stream.
        Random.seed!(sd)
        rho_noise_cpu = rand(Float64, Nx, Nz)
        rho_noise_gpu = CUDA.ones(Float64, Nx, Nz)
        copyto!(rho_noise_gpu, rho_noise_cpu)
        P_noise_cpu = rand(Float64, Nx, Nz, 2)
        P_noise_gpu = CUDA.ones(Float64, Nx, Nz, 2)
        copyto!(P_noise_gpu, P_noise_cpu)

        # Fields. Start from the homogeneous isotropic state perturbed by 0.1%:
        # ρ = ρ0(1 + 0.002η) and P = 0.002η', with η, η' uniform in [-0.5, 0.5].
        # Defects are not seeded — they emerge from this noise.
        ρ = (CUDA.ones(Float64, Nx, Nz) .+ 0.002*( rho_noise_gpu .-0.5))*ρ0
        ρn = CUDA.zeros(Float64, Nx, Nz)          # allocated but unused
        v = CUDA.zeros(Float64, Nx, Nz, 2);
        P = CUDA.ones(Float64, Nx, Nz, 2) .* ( P_noise_gpu .-0.5) * 0.002
        h = CUDA.zeros(Float64, Nx, Nz, 2)        # molecular field
        μ = CUDA.zeros(Float64, Nx, Nz)           # chemical potential
        σ = CUDA.zeros(Float64, Nx, Nz, 4)        # stress (xx, xz, zx, zz)
        F = CUDA.zeros(Float64, Nx, Nz, 2)        # allocated but unused

        # Host-side staging buffers for snapshot writing
        ρ_cpu = zeros(Float64, Nx, Nz)
        v_cpu = zeros(Float64, Nx, Nz, 2)
        P_cpu = zeros(Float64, Nx, Nz, 2)

        # Derivatives. DP packs, per component i: ∂xPi, ∂zPi, ΔPi
        DP = CUDA.zeros(Float64, Nx, Nz, 6)
        Dρ = CUDA.zeros(Float64, Nx, Nz, 2)       # ∂xρ, ∂zρ
        ∂v = CUDA.zeros(Float64, Nx, Nz, 4)       # ∂xvx, ∂xvz, ∂zvx, ∂zvz

        # Fourier-space buffers, on the rfft half-grid (Lkx x Nz)
        fv = CUDA.zeros(ComplexF64, Lkx, Nz, 2)
        fρ = CUDA.zeros(ComplexF64, Lkx, Nz)
        fP = CUDA.zeros(ComplexF64, Lkx, Nz, 2)
        fσ = CUDA.zeros(ComplexF64, Lkx, Nz, 4)

        # Spectral multipliers: i*kx, i*kz and -(kx²+kz²)
        factordx = CUDA.zeros(ComplexF64, Lkx, Nz)
        factordz = CUDA.zeros(ComplexF64, Lkx, Nz)
        factorΔ = CUDA.zeros(Float64, Lkx, Nz)

        @cuda threads = block_dim blocks = gridFFT_dim kernel_comp_diff_factor(factordx, factordz, factorΔ, kx, kz, kx2, kz2, Lkx)

        # Compile each kernel once and reuse the handle, so the launch
        # configuration is not recomputed on every one of millions of steps.
        comp_μhσ! = @cuda launch=false kernel_comp_μhσ!(σ, μ, h, ρ, Dρ, P, DP, ar, ap, kp, ζρ, ζp, ζp2, Δμ, ν1, ν2, ρ0)
        solve_kspace! = @cuda launch=false kernel_solve_kspace!(fv, fσ, kx, kz, Lkx, ξ)
        comp_P! = @cuda launch=false kernel_comp_P!(P, h, DP, v, ∂v, ρ, ν1, ν2, γ, ε, Δμ, Δt)

        # t is physical time, not a step counter: the step size adapts, so
        # `prin` and `check` are the next times to write a snapshot / retune Δt.
        t = 0; prin = 0; check = 0;
        CUDA.@time while t < t_fin + t_prin
           # --- snapshot -------------------------------------------------
           if t >= prin || t ≈ t_fin
                # Bail out rather than burn GPU hours on a diverged run.
                any(isnan, ρ) && return 1
                numm = @sprintf("%010d", t)  #print time iteration step, number format
                global nm = string(file,"Data/data",numm,".jld")
                copyto!(ρ_cpu, ρ)
                copyto!(v_cpu, v)
                copyto!(P_cpu, P)
                save(nm, "rho", ρ_cpu, "v", v_cpu, "P", P_cpu)
                println("idx = ", idx, ", t = ", floor(Int,t), ", dt = ", Δt)
                prin += t_prin
            end
            # --- adaptive time step ---------------------------------------
            # Grow Δt by 25% per check, but cap it by `Δtmin` (an upper bound
            # despite the name, set from DF.csv) and by a CFL condition that
            # keeps a fluid element from crossing 5% of a cell per step.
            if t >= check
                global @views Δt = minimum([Δt*1.25, Δtmin, 0.05*Δx/sqrt(maximum(abs.(v[:,:,1].^2+v[:,:,2].^2)))])
                check += t_check
            end

            # --- spectral derivatives of ρ and P --------------------------
            fρ .= W * ρ

            for i=1:2
               @views fP[:,:,i] .= W * P[:,:,i]
               @views DP[:,:,3*i-2] .= Wi * (factordx .* fP[:,:,i])
               @views DP[:,:,3*i-1] .= Wi * (factordz .* fP[:,:,i])
               @views DP[:,:,3*i] .= Wi * (factorΔ .* fP[:,:,i])
            end
            @views Dρ[:,:,1] .= Wi * (factordx .* fρ)
            @views Dρ[:,:,2] .= Wi * (factordz .* fρ)

            # --- constitutive relations: μ, h, σ --------------------------
            comp_μhσ!(σ, μ, h, ρ, Dρ, P, DP, ar, ap, kp, ζρ, ζp, ζp2, Δμ, ν1, ν2, ρ0; threads = block_dim, blocks = grid_dim)
            for i=1:4
                @views fσ[:,:,i] .= W * σ[:,:,i]
            end

            # --- force balance -> velocity --------------------------------
            solve_kspace!(fv, fσ, kx, kz, Lkx, ξ; threads = block_dim, blocks = gridFFT_dim)

            for i=1:2
                @views ∂v[:,:,i] .= Wi * (factordx .* fv[:,:,i])
                @views ∂v[:,:,i+2] .= Wi * (factordz .* fv[:,:,i])
                @views v[:,:,i] .= Wi * fv[:,:,i]
            end

            # --- advance P ------------------------------------------------
            comp_P!(P, h, DP, v, ∂v, ρ, ν1, ν2, γ, ε, Δμ, Δt; threads = block_dim, blocks = grid_dim)

            # --- advance ρ: -∇·(ρv) + M Δμ + turnover ---------------------
            # The flux divergence is taken spectrally on the product ρv, so the
            # advection is conservative up to the round-off of the transforms.
            @views ∂xρvx = Wi * (factordx .* (W * (ρ.*v[:,:,1]) ) )
            @views ∂zρvz = Wi * (factordz .* (W * (ρ.*v[:,:,2]) ) )

            ∇M∇μ = M * Wi * (factorΔ .* (W * μ))
            @. ρ += Δt*(-∂xρvx - ∂zρvz + ∇M∇μ + kd*(ρ0-ρ))# + kpoly
            t += Δt
        end
    end
end

main()
