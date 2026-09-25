# =============================================================================
#  Entry point for one simulation.
#
#  Run as:  SIM_IDX=<row of DF.csv> DATA_DIR=<output root> julia --optimize=3 2D.jl
#  Under Slurm, SIM_IDX defaults to SLURM_ARRAY_TASK_ID (see slurm/submit.sh).
#
#  Includes kernels.jl, which includes InputParameters.jl, so every physical
#  parameter and grid constant is already in scope as a global here.
#
#  Scheme (pseudo-spectral, explicit Euler in time). Per step:
#    1. FFT ρ and P, multiply by the spectral factors -> ∇ρ, ∇P, ΔP
#    2. kernel_comp_μhσ!            -> μ, h, σ in real space
#    3. FFT σ, kernel_solve_kspace! -> ṽ, then inverse FFT -> v and ∇v
#    4. kernel_comp_P!              -> advance P
#    5. advance ρ by the continuity equation with diffusion and turnover
#
#  All spatial derivatives are spectral, so the periodic box is exact to
#  round-off; only the time integration is approximate.
#
#  L. Dumoulin, C. Blanch-Mercader, K. Kruse, arXiv:2506.03795
# =============================================================================

include("kernels.jl")

# --- Constants of the scheme, at the values used for the paper ---------------
const NOISE_AMPLITUDE = 0.002   # relative size of the initial perturbation
const DT_GROWTH       = 1.25    # factor Δt may grow by at each check
const CFL_FRACTION    = 0.05    # cell fraction a flow may cross in one step

# Two-defect runs only. The separations the critical-distance scan walks
# through, and the two tests that decide a pair has annihilated: no polarity
# core left anywhere, or the cores closing back in on each other.
const SEPARATION_SCAN = 0.01:0.01:0.50
const CORE_GONE       = 0.995   # min |P| above this means no core in the domain
const SEPARATED       = 5       # cells the pair must gain to count as surviving

"""
    write_snapshot(file, t, idx, Δt, ρ, v, P, ρ_cpu, v_cpu, P_cpu)

Copy the fields back to the host and write `<file>Data/data<t>.jld`, with `t`
zero-padded to 10 digits. Holds `rho` (Nx×Nz) plus `v` and `P` (Nx×Nz×2).
The staging buffers are passed in so nothing is allocated per snapshot.
"""
function write_snapshot(file, t, idx, Δt, ρ, v, P, ρ_cpu, v_cpu, P_cpu)
    numm = @sprintf("%010d", t)
    nm = string(file, "Data/data", numm, ".jld")
    copyto!(ρ_cpu, ρ)
    copyto!(v_cpu, v)
    copyto!(P_cpu, P)
    save(nm, "rho", ρ_cpu, "v", v_cpu, "P", P_cpu)
    println("idx = ", idx, ", t = ", floor(Int, t), ", dt = ", Δt)
    return nothing
end

"""
    main(d = D) -> Int

Integrate one parameter set to `t_fin`, writing a snapshot every `t_prin`.
Returns 1 if the density went NaN.

For two-defect runs `d` is the initial separation, and the return value says
what became of the pair: `0` if it annihilated, `d` if it survived. That is
what lets the caller scan for the critical separation.
"""
function main(d = D)
    @inbounds begin
        # Initial noise is drawn on the host from a seeded RNG and copied to the
        # device, so a run is reproducible from `seed` in DF.csv regardless of
        # the GPU's own RNG stream. The draw order (ρ then P) is part of that.
        Random.seed!(sd)
        ρ_noise = CuArray(rand(Float64, Nx, Nz))

        # Start from the homogeneous isotropic state perturbed by 0.1%:
        # ρ = ρ0(1 + 0.002η), η uniform in [-0.5, 0.5].
        ρ = (1.0 .+ NOISE_AMPLITUDE .* (ρ_noise .- 0.5)) .* ρ0

        # The polarity either carries the same perturbation, with defects left
        # to emerge from it - which is what every published large-domain run
        # did - or is seeded with a defect pair a distance D apart, for the
        # two-defect runs. The draw order is unchanged in the noise case.
        P = if seed_defect_pair
            Pd = CUDA.ones(Float64, Nx, Nz, 2)
            @cuda threads = block_dim blocks = grid_dim kernel_ini_P!(Pd, Nx, Nz, d)
            Pd
        else
            P_noise = CuArray(rand(Float64, Nx, Nz, 2))
            (P_noise .- 0.5) .* NOISE_AMPLITUDE
        end

        # Working fields
        v = CUDA.zeros(Float64, Nx, Nz, 2)
        h = CUDA.zeros(Float64, Nx, Nz, 2)        # molecular field
        μ = CUDA.zeros(Float64, Nx, Nz)           # chemical potential
        σ = CUDA.zeros(Float64, Nx, Nz, 4)        # stress (xx, xz, zx, zz)

        # Host-side staging buffers for snapshot writing
        ρ_cpu = zeros(Float64, Nx, Nz)
        v_cpu = zeros(Float64, Nx, Nz, 2)
        P_cpu = zeros(Float64, Nx, Nz, 2)

        # Derivatives. DP packs, per component i: ∂xPi, ∂zPi, ΔPi
        DP = CUDA.zeros(Float64, Nx, Nz, 6)
        Dρ = CUDA.zeros(Float64, Nx, Nz, 2)       # ∂xρ, ∂zρ
        ∂v = CUDA.zeros(Float64, Nx, Nz, 4)       # ∂xvx, ∂xvz, ∂zvx, ∂zvz

        # Fourier-space buffers, on the rfft half-grid (Lkx × Nz)
        fv = CUDA.zeros(ComplexF64, Lkx, Nz, 2)
        fρ = CUDA.zeros(ComplexF64, Lkx, Nz)
        fP = CUDA.zeros(ComplexF64, Lkx, Nz, 2)
        fσ = CUDA.zeros(ComplexF64, Lkx, Nz, 4)

        # Spectral multipliers: i·kx, i·kz and -(kx²+kz²)
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
        # `prin` and `check` are the next times to snapshot / retune Δt.
        t = 0; prin = 0; check = 0
        t_end = t_fin + t_prin

        CUDA.@time while t < t_end
            # --- snapshot ---------------------------------------------------
            if t >= prin || t ≈ t_fin
                # Bail out rather than burn GPU hours on a diverged run.
                any(isnan, ρ) && return 1
                write_snapshot(file, t, idx, Δt, ρ, v, P, ρ_cpu, v_cpu, P_cpu)
                prin += t_prin

                # --- two-defect runs only --------------------------------
                # Follow the pair and stop as soon as its fate is decided,
                # rather than integrating to the horizon either way.
                if seed_defect_pair && t > 1
                    # no polarity core anywhere: the pair has annihilated
                    @views sqrt(minimum(P[:,:,1].^2 .+ P[:,:,2].^2)) > CORE_GONE && return 0

                    # the cores are the |P| minima on the mid-row, one per
                    # half; `dist` is how far the pair has moved apart since
                    # it was seeded, in units of the domain width.
                    @views dist = abs(-findmin(P[1:div(Nx,2),div(Nz,2),1].^2 .+ P[1:div(Nx,2),div(Nz,2),2].^2)[2] + div(Nx,2)-1+findmin(P[div(Nx,2):end,div(Nz,2),1].^2 .+ P[div(Nx,2):end,div(Nz,2),2].^2)[2])/Nx - d
                    dist < 0 && return 0                      # closing in
                    t > 20 && dist > SEPARATED/Nx && return d # holding apart
                end
            end

            # --- adaptive time step -----------------------------------------
            # Grow Δt by 25% per check, capped by `Δtmin` (an upper bound
            # despite the name) and by a CFL condition that keeps a fluid
            # element from crossing 5% of a cell per step.
            if t >= check
                global @views Δt = min(Δt*DT_GROWTH, Δtmin,
                                       CFL_FRACTION*Δx/sqrt(maximum(v[:,:,1].^2 .+ v[:,:,2].^2)))
                check += t_check
            end

            # --- spectral derivatives of ρ and P ----------------------------
            fρ .= W * ρ

            for i = 1:2
                @views fP[:,:,i] .= W * P[:,:,i]
                @views DP[:,:,3*i-2] .= Wi * (factordx .* fP[:,:,i])
                @views DP[:,:,3*i-1] .= Wi * (factordz .* fP[:,:,i])
                @views DP[:,:,3*i]   .= Wi * (factorΔ  .* fP[:,:,i])
            end
            @views Dρ[:,:,1] .= Wi * (factordx .* fρ)
            @views Dρ[:,:,2] .= Wi * (factordz .* fρ)

            # --- constitutive relations: μ, h, σ ----------------------------
            comp_μhσ!(σ, μ, h, ρ, Dρ, P, DP, ar, ap, kp, ζρ, ζp, ζp2, Δμ, ν1, ν2, ρ0; threads = block_dim, blocks = grid_dim)
            for i = 1:4
                @views fσ[:,:,i] .= W * σ[:,:,i]
            end

            # --- force balance -> velocity ----------------------------------
            solve_kspace!(fv, fσ, kx, kz, Lkx, ξ; threads = block_dim, blocks = gridFFT_dim)

            for i = 1:2
                @views ∂v[:,:,i]   .= Wi * (factordx .* fv[:,:,i])
                @views ∂v[:,:,i+2] .= Wi * (factordz .* fv[:,:,i])
                @views v[:,:,i]    .= Wi * fv[:,:,i]
            end

            # --- advance P ---------------------------------------------------
            comp_P!(P, h, DP, v, ∂v, ρ, ν1, ν2, γ, ε, Δμ, Δt; threads = block_dim, blocks = grid_dim)

            # --- advance ρ: -∇·(ρv) + M Δμ + turnover ------------------------
            # The flux divergence is taken spectrally on the product ρv, so the
            # advection is conservative up to the round-off of the transforms.
            @views ∂xρvx = Wi * (factordx .* (W * (ρ .* v[:,:,1])))
            @views ∂zρvz = Wi * (factordz .* (W * (ρ .* v[:,:,2])))

            ∇M∇μ = M * Wi * (factorΔ .* (W * μ))
            @. ρ += Δt*(-∂xρvx - ∂zρvz + ∇M∇μ + kd*(ρ0-ρ))
            t += Δt
        end
    end
end

if scan_separation
    # The critical separation: the smallest seeding distance at which the pair
    # survives instead of annihilating. main returns d only when it survives.
    for d in SEPARATION_SCAN
        if d == main(d)
            println("critical separation: idx = ", idx, ", d = ", d)
            break
        end
    end
else
    main()
end
