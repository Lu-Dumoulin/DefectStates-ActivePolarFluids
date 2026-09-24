include("../Utilities/using.jl")
using_pkg("CairoMakie, JLD, Printf, LaTeXStrings, FixedPointNumbers, DelimitedFiles, CSV, DataFrames, FileIO, Base64, Colors, LinearAlgebra, Statistics, Images, ImageSegmentation, FFTW, DelaunayTriangulation, ImageView, Random, LsqFit, GaussianMixtures, ShiftedArrays")
using_mod(".PictUtils")
CairoMakie.activate!(type = "png")

# The phase-diagram tables ship with the repository, so the figure functions
# run out of the box. Override with DATA_DIR to point at your own run set; the
# passes that rebuild DF_analyse.csv from .jld snapshots need one, the figure
# functions only need these tables. Originally hard-coded to "D:/PDpaper/".
dir_df = get(ENV, "DATA_DIR", abspath(joinpath(@__DIR__, "PDpaper")) * "/")
endswith(dir_df, "/") || (dir_df *= "/")

df1_ = CSV.read(dir_df*"DF.csv", DataFrame)
df1b_ = CSV.read(dir_df*"DF2.csv", DataFrame)
df2 = CSV.read(dir_df*"DF_ph2.csv", DataFrame)
df2b = CSV.read(dir_df*"DF2_ph2.csv", DataFrame)

df1 = df1_[df1_.zetarho .== 4.0, :]
df1b = df1b_[df1_.zetarho .== 4.0, :]

trho0 = sort(unique(df1[:,:rho0]))
tkd = sort(unique(df1[:,:kd]))

Nsim = length(trho0)*length(tkd)

# for x in eachrow(df2[:, [:rho0, :kd]])
#     if x ∈ eachrow(df1[:, [:rho0, :kd]])
#        println(df1[(df1.rho0 .== x.rho0) .& (df1.kd .== x.kd), :fn])
#     end
# end

# Two further run sets this file reads snapshots from, beyond the tables in
# PDpaper/. Defaults keep the original names; override to point at your own.
dird1 = get(ENV, "ADPTDT_DIR",    joinpath(dir_df, "AdptDt") * "/")
dird2 = get(ENV, "PHASEDIAG_DIR", joinpath(dir_df, "phasediag2") * "/")

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function expend_periodic_array(ar, ex)
    Nx, Ny = size(ar)
    nNx = Nx + 2*ex
    nNy = Ny + 2*ex
    nar = zeros(Float64, nNx, nNy)
    # Mid
    nar[1+ex:Nx+ex, 1+ex:Ny+ex] .= ar
    # Left
    nar[1+ex:Nx+ex, 1:ex] .= ar[:,Ny-ex+1:Ny]
    # Right
    nar[1+ex:Nx+ex, nNy-ex+1:nNy] .= ar[:,1:ex]
    # Top
    nar[1:ex, 1+ex:Ny+ex] .= ar[Nx-ex+1:Nx,:]
    # Bot
    nar[Nx+1+ex:nNx, 1+ex:Ny+ex] .= ar[1:ex,:]
    # Top Left
    nar[1:ex,1:ex] .= ar[Nx-ex+1:Nx,Ny-ex+1:Ny]
    # Bot Right
    nar[Nx+ex+1:nNx,Ny+ex+1:nNy] .= ar[1:ex,1:ex]
    # Top Right
    nar[1:ex,Ny+ex+1:nNy] .= ar[Nx-ex+1:Nx,1:ex]
    # Bot Left
    nar[Nx+ex+1:nNx,1:ex] .= ar[1:ex,Ny-ex+1:Ny]
    return nar
end

function shift_array(ar, x, y)
    s = ShiftedArrays.circshift(ar, (x, y))
    return copy(s)
end

df = isfile(joinpath(dir_df,"DF_analyse.csv")) ? CSV.read(dir_df*"DF_analyse.csv", DataFrame) : DataFrame(rho0 = zeros(Float64, Nsim), kd = zeros(Float64, Nsim), tauc = zeros(Float64, Nsim), pNseg=zeros(Int, Nsim), pmean=zeros(Float32, Nsim), pstd=zeros(Float32, Nsim), pArea=zeros(Float32, Nsim), rNseg=zeros(Int, Nsim), rmean=zeros(Float32, Nsim), rstd=zeros(Float32, Nsim), rArea=zeros(Float32, Nsim), ndefects=zeros(Int, Nsim))


    tab_i_tau = sort(unique(df[:, :kd]), rev=true)
    tab_i_rho0 = sort(unique(df[:, :rho0]))

function seg_last_t(dir_data)

    dird = dir_data*"Data/"
    files_names = readdir(dird)
    num_files = length(files_names)
    println("Data in $dird")

    sttopp = num_files # == 21 ? num_files : num_files-1

    if sttopp > 1
        # dist = []
        # for i=1:sttopp
            density = load(string(dird,files_names[sttopp]), "rho");

            ex = 200
            tex = round(Int, ex*2)
            rho = expend_periodic_array(density, ex);
            @show meanrho = mean(rho); @show stdrho = std(rho)
            @show maximum(rho); @show minimum(rho)
            # bwrho = (rho .< meanrho) .* (maximum(rho)-minimum(rho) > 1e-4)
            # bwrho = (rho .< meanrho-stdrho) .& (maximum(rho)-minimum(rho) > 1e-4)
            bwrho = .!((rho .< (meanrho*0.9 + minimum(rho)*0.1 )) .& (maximum(rho)-minimum(rho) > 1e-5))
            # bwrho = (rho .< meanrho*0.99)#-0.3*(meanrho - minimum(rho)))
            @show sum(bwrho)

            if sum(bwrho) > 0
                distrho = 1 .- distance_transform(feature_transform(bwrho));
                th_dist = mean(distrho) > -2 ? -1 : mean(distrho)*0.5
                segrho = felzenszwalb(distrho .< th_dist, 5, 50);
                Nsegrho = length(segrho.segment_labels)
                arrrho = [Float64(segrho.segment_pixel_count[k]) for k = 1:Nsegrho]
                length(arrrho) > 0 ? deleteat!(arrrho, findmax(arrrho)[2]) : nothing
                if length(arrrho) > 0
                    arrrho *= 1e-4
                    # @show Nsegrho
                    # @show mean(arrrho)
                    # @show std(arrrho)
                    # @show area = sum(arrrho)
                end
                return Nsegrho, arrrho
                # push!(dist, arrrho)
            end

        # end
        # maxlen = maximum(length.(dist))
        # df_dist = DataFrame()
        # for k = 1:sttopp
        #     df_dist[!, "$k"] = vcat(dist[k], fill(missing, maxlen - length(dist[k])))
        # end
        # # return df_dist
        # CSV.write(dir_save_df, df_dist)
    end
end

function seg_all_t(idx, dir_data)
    dir_save_df = dir_df*"Area_ditributions/"
    mkpath(dir_save_df)
    path = string("$idx/")
    dir_save_df *= string("$(idx).csv")
    dird = dir_data*"Data/"
    files_names = readdir(dird)
    num_files = length(files_names)
    println("Data in $dird")

    sttopp = num_files # == 21 ? num_files : num_files-1

    if sttopp > 1
        dist = []
        for i=1:sttopp
            density = load(string(dird,files_names[i]), "rho");

            ex = 200
            tex = round(Int, ex*2)
            rho = expend_periodic_array(density, ex);
            @show meanrho = mean(rho); @show stdrho = std(rho)
            @show maximum(rho); @show minimum(rho)
            # bwrho = (rho .< meanrho) .* (maximum(rho)-minimum(rho) > 1e-4)
            # bwrho = (rho .< meanrho-stdrho) .& (maximum(rho)-minimum(rho) > 1e-4)
            bwrho = .!((rho .< (meanrho*0.9 + minimum(rho)*0.1 )) .& (maximum(rho)-minimum(rho) > 1e-5))
            # bwrho = (rho .< meanrho*0.99)#-0.3*(meanrho - minimum(rho)))
            @show sum(bwrho)

            if sum(bwrho) > 0
                distrho = 1 .- distance_transform(feature_transform(bwrho));
                th_dist = mean(distrho) > -2 ? -1 : mean(distrho)*0.5
                segrho = felzenszwalb(distrho .< th_dist, 5, 50);
                Nsegrho = length(segrho.segment_labels)
                arrrho = [Float64(segrho.segment_pixel_count[k]) for k = 1:Nsegrho]
                length(arrrho) > 0 ? deleteat!(arrrho, findmax(arrrho)[2]) : nothing
                if length(arrrho) > 0
                    arrrho = arrrho * 1e-4
                    @show Nsegrho
                    @show mean(arrrho)
                    @show std(arrrho)
                    @show area = sum(arrrho)
                end
                push!(dist, arrrho)
            end

        end
        maxlen = maximum(length.(dist))
        df_dist = DataFrame()
        for k = 1:sttopp
            df_dist[!, "$k"] = vcat(dist[k], fill(missing, maxlen - length(dist[k])))
        end
        # return df_dist
        CSV.write(dir_save_df, df_dist)
    end
end

function seg_all_t_forall_rho_kd()
    for r ∈ trho0, kd ∈ tkd
        dir_data = ""
        idx = 0
        tmp = df2[(df2.rho0 .== r) .& (df2.kd .== kd), :fn]
        if isempty(tmp)
            idx += df1[(df1.rho0 .== r) .& (df1.kd .== kd), :fn][1]
            dir_data *= dird1*string("$idx/")
        else
            idx += tmp[1]
            dir_data *= dird2*string("$idx/")
        end
        println(idx, " ", dir_data)
        seg_all_t(idx, dir_data)
    end
end


@inline exp_decay(t, p) = p[1] .* exp.(-t ./ p[2])
# p = [A0, tau_c]

function ev_autocor(dir_data, t=10)
    dird = dir_data*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    if length(files_names) > t+5
        ρ_end = load(string(dird,files_names[end]), "rho");
        avg_ρ_end = mean(ρ_end)
        δρ_end = ρ_end .- avg_ρ_end;
        res = zeros(t)
        for i=1:t
            ρ_past = load(string(dird,files_names[end-i]), "rho");
            avg_ρ_past = mean(ρ_past)
            δρ_past = ρ_past .- avg_ρ_past;
            res[i] = mean( δρ_end .* δρ_past ) / sqrt( mean(δρ_end .^ 2) * mean(δρ_past .^ 2) )
        end 
        p0 = [1.0, 2.0]
        fit = curve_fit(exp_decay, 1:t, res, p0)
        p = coef(fit)
        A0, τc = p
        return τc
    else 
        return -1
    end
end

function ev_autocor_shift(dir_data, t=10)
    dird = dir_data*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    if length(files_names) > t+5
        ρ_end = load(string(dird,files_names[end]), "rho");
        avg_ρ_end = mean(ρ_end)
        δρ_end = ρ_end .- avg_ρ_end;
        res = zeros(t)
        shift = zeros(t+1,2)
        for i=1:t
            ρ_past = load(string(dird,files_names[end-i]), "rho");
            avg_ρ_past = mean(ρ_past)
            tmp = zeros(7,7)
            for x=-3:3, y=-3:3
                ρ_shift = shift_array(ρ_past, x+shift[i,1], y+shift[i,2])
                δρ_past = ρ_shift .- avg_ρ_past;
                tmp[x+4,y+4] = mean( δρ_end .* δρ_past ) / sqrt( mean(δρ_end .^ 2) * mean(δρ_past .^ 2) )
            end
            r, id = findmax(tmp)
            res[i] = r
            shift[i+1,1] = id[1] - 4 + shift[i,1]; shift[i+1,2] = id[2] - 4 + shift[i,2]
        end 
        p0 = [1.0, 2.0]
        fit = curve_fit(exp_decay, 1:t, res, p0)
        p = coef(fit)
        A0, τc = p
        return τc, shift
    else 
        return -1
    end
end

function ev_autocor_shift_std(dir_data, t=10)
    dird = dir_data*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    if length(files_names) > t+5
        ρ_end = load(string(dird,files_names[end]), "rho");
        avg_ρ_end = mean(ρ_end)
        δρ_end = ρ_end .- avg_ρ_end;
        @show crit_end = maximum(abs.(δρ_end)) / avg_ρ_end
        res = zeros(t)
        shift = zeros(t+1,2)
        for i=1:t
            ρ_past = load(string(dird,files_names[end-i]), "rho");
            avg_ρ_past = mean(ρ_past)
            tmp = zeros(7,7)
            for x=-3:3, y=-3:3
                ρ_shift = shift_array(ρ_past, x+shift[i,1], y+shift[i,2])
                δρ_past = ρ_shift .- avg_ρ_past;
                crit_past = maximum(abs.(δρ_past)) / avg_ρ_past
                tmp[x+4,y+4] = (crit_end < 1e-6 && crit_past < 1e-6) ? 1.0 : mean( δρ_end .* δρ_past ) / sqrt( mean(δρ_end .^ 2) * mean(δρ_past .^ 2) )
            end
            r, id = findmax(tmp)
            res[i] = r
            shift[i+1,1] = id[1] - 4 + shift[i,1]; shift[i+1,2] = id[2] - 4 + shift[i,2]
        end 
        p0 = [1.0, 2.0]
        fit = curve_fit(exp_decay, 1:t, res, p0)
        p = coef(fit)
        A0, τc = p
        return τc, shift, res
    else 
        return -1
    end
end

function make_df_corr_time(t)
    dfc = DataFrame(rho0 = Float64[], kd = Float64[], tauc = Float64[])
    for r ∈ trho0, kd ∈ tkd
        dir_data = ""
        idx = 0
        tmp = df2[(df2.rho0 .== r) .& (df2.kd .== kd), :fn]
        if isempty(tmp)
            idx += df1[(df1.rho0 .== r) .& (df1.kd .== kd), :fn][1]
            dir_data *= dird1*string("$idx/")
        else
            idx += tmp[1]
            dir_data *= dird2*string("$idx/")
        end
        # println(idx, " ", dir_data)
        τc = ev_autocor(dir_data, t)
        push!(dfc, [tho0, kd, τc])
    end
    CSV.write(dir_df*"DF_tauc_$t.csv", dfc)
end

function analyse(t)
    row = 0
    for r ∈ trho0, kd ∈ tkd
        row += 1
        df[row, :rho0] = r; df[row, :kd] = kd
        dir_data = ""
        idx = 0
        tmp = df2[(df2.rho0 .== r) .& (df2.kd .== kd), :fn]
        if isempty(tmp)
            idx += df1[(df1.rho0 .== r) .& (df1.kd .== kd), :fn][1]
            dir_data *= dird1*string("$idx/")
            df[row, :ndefects] = df1b_[idx, :Ndefects]
        else
            idx += tmp[1]
            dir_data *= dird2*string("$idx/")
            df[row, :ndefects] = df2b[idx, :Ndefects]
        end
        τc = ev_autocor(dir_data, t)
        df[row, :tauc] = τc 

        nseg, arrho = seg_last_t(dir_data)
        df[row, :rNseg] = nseg
        df[row, :rmean] = mean(arrho)
        df[row, :rstd] = std(arrho)
        df[row, :rArea] = sum(arrho)
        println(df[row, :])
    end
    CSV.write(dir_df*"DF_analyse.csv", df)
end

function analyse_shift(t)
    row = 0
    for r ∈ trho0, kd ∈ tkd
        row += 1
        df[row, :rho0] = r; df[row, :kd] = kd
        dir_data = ""
        idx = 0
        tmp = df2[(df2.rho0 .== r) .& (df2.kd .== kd), :fn]
        if isempty(tmp)
            idx += df1[(df1.rho0 .== r) .& (df1.kd .== kd), :fn][1]
            dir_data *= dird1*string("$idx/")
            df[row, :ndefects] = df1b_[idx, :Ndefects]
        else
            idx += tmp[1]
            dir_data *= dird2*string("$idx/")
            df[row, :ndefects] = df2b[idx, :Ndefects]
        end
        τc, shift = ev_autocor_shift(dir_data, t)
        df[row, :tauc] = τc 

        nseg, arrho = seg_last_t(dir_data)
        df[row, :rNseg] = nseg
        df[row, :rmean] = mean(arrho)
        df[row, :rstd] = std(arrho)
        df[row, :rArea] = sum(arrho)
        println(df[row, :])
        println(shift)
    end
    CSV.write(dir_df*"DF_analyse2b.csv", df)
end

function analyse_shift_std(t)
    df = DataFrame(rho0 = zeros(Float64, Nsim), kd = zeros(Float64, Nsim), tauc = zeros(Float64, Nsim), pNseg=zeros(Int, Nsim), pmean=zeros(Float32, Nsim), pstd=zeros(Float32, Nsim), pArea=zeros(Float32, Nsim), rNseg=zeros(Int, Nsim), rmean=zeros(Float32, Nsim), rstd=zeros(Float32, Nsim), rArea=zeros(Float32, Nsim), ndefects=zeros(Int, Nsim))
    df_tikz = DataFrame(y = zeros(Float64, Nsim), x = zeros(Float64, Nsim), q = zeros(Float64, Nsim), s=zeros(Float32, Nsim), w=zeros(Float32, Nsim), a=zeros(Float32, Nsim))
    row = 0
    for r ∈ trho0, kd ∈ tkd
        row += 1
        df[row, :rho0] = r; df[row, :kd] = kd
        tau_i = findfirst(isequal(kd), tab_i_tau)
        rho0_i = findfirst(isequal(r), tab_i_rho0)
        df_tikz[row, :y] = rho0_i; df_tikz[row, :x] = tau_i
        dir_data = ""
        idx = 0
        tmp = df2[(df2.rho0 .== r) .& (df2.kd .== kd), :fn]
        if isempty(tmp)
            idx += df1[(df1.rho0 .== r) .& (df1.kd .== kd), :fn][1]
            dir_data *= dird1*string("$idx/")
            df[row, :ndefects] = df1b_[idx, :Ndefects]
            df_tikz[row, :a] = df1b_[idx, :Ndefects]/10/10
        else
            idx += tmp[1]
            dir_data *= dird2*string("$idx/")
            df[row, :ndefects] = df2b[idx, :Ndefects]
            df_tikz[row, :a] = df2b[idx, :Ndefects]/10/10
        end
        τc, shift, _ = ev_autocor_shift_std(dir_data, t)
        df[row, :tauc] = τc 
        df_tikz[row, :q] = τc 

        nseg, arrho = seg_last_t(dir_data)
        df[row, :rNseg] = nseg
        df[row, :rmean] = mean(arrho)
        df_tikz[row, :s] = isnan(mean(arrho)) ? 0 : mean(arrho)
        df[row, :rstd] = std(arrho)
        df_tikz[row, :w] = isnan(std(arrho)) ? 0 : std(arrho)
        df[row, :rArea] = sum(arrho)
        println(df_tikz[row, :])
        println(shift)
    end
    CSV.write(dir_df*"DF_analyse4_$t.csv", df)
    CSV.write(dir_df*"DF_analyse_tikz_rho_$t.csv", df_tikz)
end


function phase_diagram()
    df = CSV.read(dir_df*"DF_analyse2.csv", DataFrame)
    N_points = Nsim

    tab_i_tau = sort(unique(df[:, :kd]), rev=true)
    tab_i_rho0 = sort(unique(df[:, :rho0]))

    x1 = []; y1 = [] # Homogeneous no defects
    x2 = []; y2 = [] # Foams (thin or large)
    x3 = []; y3 = [] # Transition from foams to waves
    x4 = []; y4 = [] # Waves 
    x5 = []; y5 = [] # Vortex
    x6 = []; y6 = [] # Sea of defects
    x7 = []; y7 = [] # ?

    for i=1:N_points
        tau = findfirst(isequal(df[i, :kd]), tab_i_tau)
        rho0 = findfirst(isequal(df[i,:rho0]), tab_i_rho0)
        idx = findfirst(x -> (x.rho0 == tab_i_rho0[rho0] && x.kd == tab_i_tau[tau]), eachrow(df))

        ravg = df[idx, :rmean]; rstd = df[idx, :rstd]; rarea = df[idx,:rArea]/(14*14); Ndef = df[idx, :ndefects]/10/10; τc = df[idx, :tauc]

        ndef__ = Ndef < 0.1 ? 0 : 1 # (Ndef < 0.5 ? 1 : 2)
        # ravg__ = ravg == 0 ? 0 : (ravg < 0.5 ? 1 : 2)
        # rarea__ = rarea <= 0.01 ? 0 : (rarea < 0.25 ? 1 : 2)
        runi__ = ravg == 0 ? -1 : (rstd <= 0.9*ravg ? 0 : 1)
        τc__ = τc < 25 ? 0 : 1# (τc <= 100 ? 1 : 2)

        if ndef__ == 0 && τc__ == 1
            push!(x1, tau); push!(y1, rho0)
        elseif ndef__ == 0 && τc__ == 0                                     # Waves
            push!(x4, tau); push!(y4, rho0)
        elseif ndef__ == 1 && τc__ == 0 && runi__ == 0                      # Foams
            push!(x2, tau); push!(y2, rho0)
        elseif ndef__ == 1 && τc__ == 1 && runi__ == 0                      # Sea of defect
            push!(x6, tau); push!(y6, rho0)
        elseif ndef__ == 1 && τc__ == 1 && runi__ == 1                      # Vortex Glass
            push!(x5, tau); push!(y5, rho0)
        elseif ndef__ == 1 && τc__ == 0 && runi__ == 1                      # Transition
            push!(x3, tau); push!(y3, rho0)
        else
            push!(x7, tau); push!(y7, rho0)
        end

    end

    fig = Figure(size=(500,400))
    ax = Makie.Axis(fig[1,1], title=latexstring(L"\zeta_\rho = ",4), xticks = (1:length(tab_i_tau), [@sprintf("%2.2f", 1/i) for i in tab_i_tau]), yticks = (1:length(tab_i_rho0), ["$i" for i in tab_i_rho0]), xlabel=L"\tau", ylabel = L"\rho_0", xticklabelrotation = pi/4)
    !isempty(x1) ? scatter!(ax, x1, y1, marker = :xcross, color = 10*ones(size(x1)), colormap = :RdPu, colorrange=(2,10), label = "HUPS") : nothing
    !isempty(x2) ? scatter!(ax, x2, y2, marker = :circle, color = 5*ones(size(x2)), colormap = :RdPu, colorrange=(2,10), label = "Foams") : nothing
    !isempty(x3) ? scatter!(ax, x3, y3, marker = :star5, color = 6*ones(size(x3)), colormap = :RdPu, colorrange=(2,10), label = "F + W") : nothing
    !isempty(x4) ? scatter!(ax, x4, y4, marker = :cross, color = 7*ones(size(x4)), colormap = :RdPu, colorrange=(2,10), label = "Waves") : nothing
    !isempty(x5) ? scatter!(ax, x5, y5, marker = :diamond, color = 8*ones(size(x5)), colormap = :RdPu, colorrange=(2,10), label = "Vortex") : nothing
    !isempty(x6) ? scatter!(ax, x6, y6, marker = :rect, color = 9*ones(size(x6)), colormap = :RdPu, colorrange=(2,10), label = "Sea") : nothing
    !isempty(x7) ? scatter!(ax, x7, y7, marker = '?', color = :black, colormap = :RdPu, colorrange=(2,10), label = "?") : nothing
    fig[1,2] = Legend(fig, ax, framevisible = false)
    return fig
end

function plot_all_zr_diag_3(t)
    dfa = CSV.read(dir_df*"DF2.csv", DataFrame)
    all_zr = sort(unique(dfa[:, :zetarho]))
    for i in eachindex(all_zr)
        fig = phase_diagram_zr_3(all_zr[i], t)
        save(joinpath(dir_df, "zr_phase_diag$(i)_$(t)_3.png"), fig)
    end
end

r0 = 0.7; k=0.2; println(df1[(df1.rho0 .== r0) .& (df1.kd .== k), :fn], " ", df2[(df2.rho0 .== r0) .& (df2.kd .== k), :fn]); println(df[(df.rho0 .== r0) .& (df.kd .== k), :])


function df_tikz_phase_diagram()
    df = CSV.read(dir_df*"DF_analyse2.csv", DataFrame)
    df_tikz = DataFrame(x=zeros(Nsim), y=zeros(Nsim), q=zeros(Nsim), w=zeros(Nsim), a=zeros(Nsim), s=zeros(Nsim))

    tab_i_tau = sort(unique(df[:, :kd]), rev=true)
    tab_i_rho0 = sort(unique(df[:, :rho0]))

    for i=1:Nsim
        tau = findfirst(isequal(df[i, :kd]), tab_i_tau)
        rho0 = findfirst(isequal(df[i,:rho0]), tab_i_rho0)
        idx = findfirst(x -> (x.rho0 == tab_i_rho0[rho0] && x.kd == tab_i_tau[tau]), eachrow(df))

        ravg = df[idx, :rmean]; rstd = df[idx, :rstd]; rarea = df[idx,:rArea]/(14*14); Ndef = df[idx, :ndefects]/10/10; τc = df[idx, :tauc]
        df_tikz[i,:x] = tau
        df_tikz[i,:y] = rho0
        df_tikz[i,:q] = τc
        df_tikz[i,:w] = isnan(rstd) ? 0.0 : rstd
        df_tikz[i,:a] = Ndef
        df_tikz[i,:s] = isnan(ravg) ? 0.0 : rstd

    end
    @show maximum(df_tikz[:,:q])
    @show maximum(df_tikz[:,:w])
    @show maximum(df_tikz[:,:s])
    @show maximum(df_tikz[:,:a])
    df_tikz[!,:q] ./= maximum([maximum(df_tikz[:,:q]), 100])
    df_tikz[!,:w] ./= maximum([maximum(df_tikz[:,:w]), 2])
    df_tikz[!,:a] ./= maximum(df_tikz[:,:a])
    df_tikz[!,:s] ./= maximum([maximum(df_tikz[:,:s]), 2])
    CSV.write(dir_df*"DF_tikz_norm_adjusted.csv", df_tikz)
    return df_tikz
end

function df_tikz_phase_diagram_cut()
    df = CSV.read(dir_df*"DF_analyse2.csv", DataFrame)
    df_tikz = DataFrame(x=zeros(Nsim), y=zeros(Nsim), q=zeros(Nsim), w=zeros(Nsim), a=zeros(Nsim), s=zeros(Nsim))

    tab_i_tau = sort(unique(df[:, :kd]), rev=true)
    tab_i_rho0 = sort(unique(df[:, :rho0]))

    for i=1:Nsim
        tau = findfirst(isequal(df[i, :kd]), tab_i_tau)
        rho0 = findfirst(isequal(df[i,:rho0]), tab_i_rho0)
        idx = findfirst(x -> (x.rho0 == tab_i_rho0[rho0] && x.kd == tab_i_tau[tau]), eachrow(df))

        ravg = df[idx, :rmean]; rstd = df[idx, :rstd]; rarea = df[idx,:rArea]/(14*14); Ndef = df[idx, :ndefects]/10/10; τc = df[idx, :tauc]
        df_tikz[i,:x] = maximum([50, tau])
        df_tikz[i,:y] = rho0
        df_tikz[i,:q] = τc
        df_tikz[i,:w] = isnan(rstd) ? 0.0 : rstd
        df_tikz[i,:a] = Ndef
        df_tikz[i,:s] = isnan(ravg) ? 0.0 : rstd

    end
    # @show maximum(df_tikz[:,:q])
    # @show maximum(df_tikz[:,:w])
    # @show maximum(df_tikz[:,:s])
    # @show maximum(df_tikz[:,:a])
    df_tikz[!,:q] ./= maximum([maximum(df_tikz[:,:q]), 100])
    df_tikz[!,:w] ./= maximum([maximum(df_tikz[:,:w]), 2])
    df_tikz[!,:a] ./= maximum(df_tikz[:,:a])
    df_tikz[!,:s] ./= maximum([maximum(df_tikz[:,:s]), 2])
    CSV.write(dir_df*"DF_tikz_norm_adjusted.csv", df_tikz)
    return df_tikz
end

function make_csv_exp_plot(t=20)
    df_exp_curv = DataFrame(x=0:t, y1=ones(t+1), y2=ones(t+1), y3=ones(t+1))
    t1, _, res = ev_autocor_shift_std(string(dird2, "69/"), t)
    df_exp_curv[!,:y1][2:end] = res
    t2, _, res = ev_autocor_shift_std(string(dird1, "376/"), t)
    df_exp_curv[!,:y2][2:end] = res
    t3, _, res = ev_autocor_shift_std(string(dird2, "28/"), t)
    df_exp_curv[!,:y3][2:end] = res
    @show t1, t2, t3
    CSV.write(dir_df*"DF_tikz_exp_$t.csv", df_exp_curv)
end

function ev_autocor_shift_std_2(dir_data, t=10)
    dird = dir_data*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    if length(files_names) > t+5
        ρ_end = load(string(dird,files_names[end-t]), "rho");
        avg_ρ_end = mean(ρ_end)
        δρ_end = ρ_end .- avg_ρ_end;
        @show crit_end = maximum(abs.(δρ_end)) / avg_ρ_end
        res = zeros(t)
        shift = zeros(t+1,2)
        for i=1:t
            ρ_past = load(string(dird,files_names[end-t+i]), "rho");
            avg_ρ_past = mean(ρ_past)
            tmp = zeros(7,7)
            for x=-3:3, y=-3:3
                ρ_shift = shift_array(ρ_past, x+shift[i,1], y+shift[i,2])
                δρ_past = ρ_shift .- avg_ρ_past;
                crit_past = maximum(abs.(δρ_past)) / avg_ρ_past
                tmp[x+4,y+4] = (crit_end < 1e-6 && crit_past < 1e-6) ? 1.0 : mean( δρ_end .* δρ_past ) / sqrt( mean(δρ_end .^ 2) * mean(δρ_past .^ 2) )
            end
            r, id = findmax(tmp)
            res[i] = r
            shift[i+1,1] = id[1] - 4 + shift[i,1]; shift[i+1,2] = id[2] - 4 + shift[i,2]
        end 
        p0 = [1.0, 2.0]
        fit = curve_fit(exp_decay, 1:t, res, p0)
        p = coef(fit)
        A0, τc = p
        return τc, shift, res
    else 
        return -1
    end
end

function ev_autocor_shift_std_3(dir_data, t=10)
    dird = dir_data*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    if length(files_names) > t+5
        ρ_end = load(string(dird,files_names[end]), "rho");
        avg_ρ_end = mean(ρ_end)
        δρ_end = ρ_end .- avg_ρ_end;
        @show crit_end = maximum(abs.(δρ_end)) / avg_ρ_end
        res = ones(t+1)
        shift = zeros(t+1,2)
        for i=1:t
            ρ_past = load(string(dird,files_names[end-i]), "rho");
            avg_ρ_past = mean(ρ_past)
            tmp = zeros(7,7)
            for x=-3:3, y=-3:3
                ρ_shift = shift_array(ρ_past, x+shift[i,1], y+shift[i,2])
                δρ_past = ρ_shift .- avg_ρ_past;
                crit_past = maximum(abs.(δρ_past)) / avg_ρ_past
                tmp[x+4,y+4] = (crit_end < 1e-6 && crit_past < 1e-6) ? 1.0 : mean( δρ_end .* δρ_past ) / sqrt( mean(δρ_end .^ 2) * mean(δρ_past .^ 2) )
            end
            r, id = findmax(tmp)
            res[i+1] = r
            shift[i+1,1] = id[1] - 4 + shift[i,1]; shift[i+1,2] = id[2] - 4 + shift[i,2]
        end 
        p0 = [1.0, 2.0]
        fit = curve_fit(exp_decay, 0:t, res, p0)
        p = coef(fit)
        A0, τc = p
        return τc, shift, res
    else 
        return -1
    end
end

function make_csv_exp_plot_3(t=20)
    df_exp_curv = DataFrame(x=0:t, y1=ones(t+1), y2=ones(t+1), y3=ones(t+1))
    t1, _, res = ev_autocor_shift_std_3(string(dird2, "69/"), t)
    df_exp_curv[!,:y1][1:end] = res
    t2, _, res = ev_autocor_shift_std_3(string(dird1, "376/"), t)
    df_exp_curv[!,:y2][1:end] = res
    t3, _, res = ev_autocor_shift_std_3(string(dird2, "28/"), t)
    df_exp_curv[!,:y3][1:end] = res
    @show t1, t2, t3
    CSV.write(dir_df*"DF_tikz_exp_3_$t.csv", df_exp_curv)
end