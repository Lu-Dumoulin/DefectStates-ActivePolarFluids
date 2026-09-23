include("kernels.jl")

function main()
    @inbounds begin
        Random.seed!(sd)
        rho_noise_cpu = rand(Float64, Nx, Nz)
        rho_noise_gpu = CUDA.ones(Float64, Nx, Nz)
        copyto!(rho_noise_gpu, rho_noise_cpu)
        P_noise_cpu = rand(Float64, Nx, Nz, 2)
        P_noise_gpu = CUDA.ones(Float64, Nx, Nz, 2)
        copyto!(P_noise_gpu, P_noise_cpu)
        #Fields
        # ρ = (CUDA.ones(Float64, Nx, Nz) .+ 0.002*(CUDA.rand(Float64, Nx, Nz) .-0.5))*ρ0
        ρ = (CUDA.ones(Float64, Nx, Nz) .+ 0.002*( rho_noise_gpu .-0.5))*ρ0
        ρn = CUDA.zeros(Float64, Nx, Nz)
        v = CUDA.zeros(Float64, Nx, Nz, 2); 
        # P = CUDA.ones(Float64, Nx, Nz, 2) .* ( CUDA.rand(Float64, Nx, Nz, 2) .-0.5) * 0.002
        P = CUDA.ones(Float64, Nx, Nz, 2) .* ( P_noise_gpu .-0.5) * 0.002
        
        # P = CUDA.zeros(Float64, Nx, Nz, 2) .+ ( CUDA.rand(Float64, Nx, Nz, 2) .-0.5) * 0.002
        # @cuda threads = block_dim blocks = grid_dim kernel_ini_P!(P, Nx, Nz, d) ##### SPECIFIC
        h = CUDA.zeros(Float64, Nx, Nz, 2)
        μ = CUDA.zeros(Float64, Nx, Nz)
        σ = CUDA.zeros(Float64, Nx, Nz, 4)
        F = CUDA.zeros(Float64, Nx, Nz, 2)
        # CPU Fields
        ρ_cpu = zeros(Float64, Nx, Nz)
        v_cpu = zeros(Float64, Nx, Nz, 2)
        P_cpu = zeros(Float64, Nx, Nz, 2)

        # Derivatives
        DP = CUDA.zeros(Float64, Nx, Nz, 6)
        Dρ = CUDA.zeros(Float64, Nx, Nz, 2)
        ∂v = CUDA.zeros(Float64, Nx, Nz, 4)

        # FFT:
        fv = CUDA.zeros(ComplexF64, Lkx, Nz, 2)
        fρ = CUDA.zeros(ComplexF64, Lkx, Nz)
        fP = CUDA.zeros(ComplexF64, Lkx, Nz, 2)
        fσ = CUDA.zeros(ComplexF64, Lkx, Nz, 4)

        # Factor for FFT derivatives
        factordx = CUDA.zeros(ComplexF64, Lkx, Nz)
        factordz = CUDA.zeros(ComplexF64, Lkx, Nz)
        factorΔ = CUDA.zeros(Float64, Lkx, Nz)

        @cuda threads = block_dim blocks = gridFFT_dim kernel_comp_diff_factor(factordx, factordz, factorΔ, kx, kz, kx2, kz2, Lkx)
        
        comp_μhσ! = @cuda launch=false kernel_comp_μhσ!(σ, μ, h, ρ, Dρ, P, DP, ar, ap, kp, ζρ, ζp, ζp2, Δμ, ν1, ν2, ρ0)
        solve_kspace! = @cuda launch=false kernel_solve_kspace!(fv, fσ, kx, kz, Lkx, ξ)
        comp_P! = @cuda launch=false kernel_comp_P!(P, h, DP, v, ∂v, ρ, ν1, ν2, γ, ε, Δμ, Δt)
        
        t = 0; prin = 0; check = 0;# incrkp = 200000
        # CUDA.@time for t in 0:Nt
        CUDA.@time while t < t_fin + t_prin
           # if (t%prin==0)
           if t >= prin || t ≈ t_fin
                any(isnan, ρ) && return 1
                numm = @sprintf("%010d", t)  #print time iteration step, number format
                global nm = string(file,"Data/data",numm,".jld")
                copyto!(ρ_cpu, ρ)
                copyto!(v_cpu, v)
                copyto!(P_cpu, P)
                save(nm, "rho", ρ_cpu, "v", v_cpu, "P", P_cpu)
                println("idx = ", idx, ", t = ", floor(Int,t), ", dt = ", Δt)
                prin += t_prin
                # (t > 1 && sqrt(minimum(P[:,:,1].^2 .+ P[:,:,2].^2)) > 0.995) && return 1# : nothing ##### SPECIFIC 
            end
            if t >= check
                global @views Δt = minimum([Δt*1.25, Δtmin, 0.05*Δx/sqrt(maximum(abs.(v[:,:,1].^2+v[:,:,2].^2)))])
                check += t_check
            end
            # if t >= incrkp
            #     global kp *= 2
            #     incrkp += 200000
            # end
            
            fρ .= W * ρ
            
            for i=1:2
               @views fP[:,:,i] .= W * P[:,:,i]
               @views DP[:,:,3*i-2] .= Wi * (factordx .* fP[:,:,i])
               @views DP[:,:,3*i-1] .= Wi * (factordz .* fP[:,:,i])
               @views DP[:,:,3*i] .= Wi * (factorΔ .* fP[:,:,i])
            end
            @views Dρ[:,:,1] .= Wi * (factordx .* fρ)
            @views Dρ[:,:,2] .= Wi * (factordz .* fρ)
            comp_μhσ!(σ, μ, h, ρ, Dρ, P, DP, ar, ap, kp, ζρ, ζp, ζp2, Δμ, ν1, ν2, ρ0; threads = block_dim, blocks = grid_dim)
            for i=1:4
                @views fσ[:,:,i] .= W * σ[:,:,i]
            end

            solve_kspace!(fv, fσ, kx, kz, Lkx, ξ; threads = block_dim, blocks = gridFFT_dim) 

            for i=1:2
                @views ∂v[:,:,i] .= Wi * (factordx .* fv[:,:,i])
                @views ∂v[:,:,i+2] .= Wi * (factordz .* fv[:,:,i])
                @views v[:,:,i] .= Wi * fv[:,:,i]
            end
            comp_P!(P, h, DP, v, ∂v, ρ, ν1, ν2, γ, ε, Δμ, Δt; threads = block_dim, blocks = grid_dim)
            
            @views ∂xρvx = Wi * (factordx .* (W * (ρ.*v[:,:,1]) ) )
            @views ∂zρvz = Wi * (factordz .* (W * (ρ.*v[:,:,2]) ) )
            
            ∇M∇μ = M * Wi * (factorΔ .* (W * μ))
            @. ρ += Δt*(-∂xρvx - ∂zρvz + ∇M∇μ + kd*(ρ0-ρ))# + kpoly
            t += Δt
        end
    end
end

main()
