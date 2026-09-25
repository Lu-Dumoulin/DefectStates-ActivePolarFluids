# =============================================================================
#  The routines that produced the figures, for reference.
#
#  This file is not runnable and is not meant to be. The figures were drawn
#  from simulation output that is not published with this repository - a single
#  snapshot is about a gigabyte and a run set reaches terabytes - so nothing
#  here can be executed against anything you can obtain from this repository.
#
#  It is kept because the article's appendix describes how the measured
#  quantities were obtained, and this is that code: defect counting, the shape
#  function, the low-density area statistics and so on. The measurement
#  primitives themselves are in core.jl; what is collected here is the
#  per-figure work built on top of them.
#
#  figures/README.md records which routine drew which panel.
#
#  What IS reproducible from this repository:
#    - the simulations themselves, from 2D/ and the tables in params/
#    - figure 3, from the notebook in 1D/
#    - figure 5's data, from linear_stability.jl - it needs no simulation output
# =============================================================================

include("MakePlots.jl")

# ----------------------------------------------------------------------------
#  Figure 1 - fig:schema
# ----------------------------------------------------------------------------

function make_plot_L50(tidx=Array(1:9); part = 1.0, Lz = 3, Lx = 3)
    # set_theme!(empty_theme)
    # df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    # tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    
    for i = 1:length(tidx)
        idx = tidx[i]
        dird = dir_df*"$idx/"*"Data/"
        rho = load(string(dird,readdir(dird)[end]), "rho")
        # mrho = maximum(rho)
        # rhomax = mrho > rhomax ? mrho : rhomax
        # mrho = minimum(rho)
        # rhomin = mrho < rhomin ? mrho : rhomin
        rhomax[i] = maximum(rho)
        rhomin[i] = minimum(rho)
    end
    
    fig = Figure(size=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lx for j in 1:Lz]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        # zr = df2[idx,:zetarho]
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        last_file = readdir(dird)[end]
        rho = load(string(dird,last_file), "rho")
        # Nx, Nz = size(rho)
        # pNx = round(Int, Nx*part)
        # pNz = round(Int, Nz*part)
        # x = 1:steppp:pNx
        # y = 1:steppp:pNz
        # P_cpu = load(string(dird,last_file), "P")
        # Px_cpu = P_cpu[x,y,1]
        # Pz_cpu = P_cpu[x,y,2]
        # arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
        heatmap!(axs[i+1], rho, colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        # arrows!(axs[i], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
        # Box(gls[i][1,1], color = (:white, 0.75), strokecolor = :black, strokewidth = 1)
        # Label(gls[i][1,1], L"\zeta_\rho = %$zr", padding = (10, 10, 10, 10), color = :black, font = :bold, fontsize = 30)
        println(i, ", ", rhomin[i], ", ", rhomax[i])
    end
    tightlimits!.(axs)
    hideydecorations!.(axs)
    hidexdecorations!.(axs)
    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 1)
    Colorbar(fig[:, end+1], colormap=:RdPu_8, colorrange=(rhomin[1],rhomax[1]), ticksvisible=false, size = 25, ticklabelsize=0) #, label= L"\rho", labelsize = 30, ticklabelsize=24)
    mkpath(dir_fig)
    save(joinpath(dir_fig, "phases_L50_3.png"), fig)
    return fig
end

function make_plot_L50_3zoom(tidx=[2,5,8], zoom=3; part = 1.0, Lz = 2, Lx = 3)
    # set_theme!(empty_theme)
    # df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    # tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    
    for i = 1:length(tidx)
        idx = tidx[i]
        dird = dir_df*"$idx/"*"Data/"
        rho = load(string(dird,readdir(dird)[end]), "rho")
        # mrho = maximum(rho)
        # rhomax = mrho > rhomax ? mrho : rhomax
        # mrho = minimum(rho)
        # rhomin = mrho < rhomin ? mrho : rhomin
        rhomax[i] = maximum(rho)
        rhomin[i] = minimum(rho)
    end
    
    fig = Figure(size=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    # gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lx for j in 1:Lz]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        # zr = df2[idx,:zetarho]
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        last_file = readdir(dird)[end]
        rho = load(string(dird,last_file), "rho")
        # Nx, Nz = size(rho)
        # pNx = round(Int, Nx*part)
        # pNz = round(Int, Nz*part)
        # x = 1:steppp:pNx
        # y = 1:steppp:pNz
        # P_cpu = load(string(dird,last_file), "P")
        # Px_cpu = P_cpu[x,y,1]
        # Pz_cpu = P_cpu[x,y,2]
        # arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
        heatmap!(axs[i+3], rho, colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        # arrows!(axs[i], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
        # Box(gls[i][1,1], color = (:white, 0.75), strokecolor = :black, strokewidth = 1)
        # Label(gls[i][1,1], L"\zeta_\rho = %$zr", padding = (10, 10, 10, 10), color = :black, font = :bold, fontsize = 30)
        println(i, ", ", rhomin[i], ", ", rhomax[i])
    end
    idx = tidx[zoom]
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    last_file = readdir(dird)[end]
    rho = load(string(dird,last_file), "rho")
    heatmap!(axs[3], rho[1:1001,end-1000:end], colormap=:RdPu_8, colorrange=(rhomax[zoom]-rhomin[zoom]<0.01 ? (0, rhomax[zoom]) : (rhomin[zoom],rhomax[zoom])) )
    hidespines!(axs[1], :r)
    hidespines!(axs[2], :l)
    tightlimits!.(axs)
    hideydecorations!.(axs)
    hidexdecorations!.(axs)
    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 1)
    Colorbar(fig[:, end+1], colormap=:RdPu_8, colorrange=(rhomin[1],rhomax[1]), ticksvisible=false, size = 25, ticklabelsize=0) #, label= L"\rho", labelsize = 30, ticklabelsize=24)
    mkpath(dir_fig)
    save(joinpath(dir_fig, "phases_L50_3zoom.png"), fig)
    return fig
end

# ----------------------------------------------------------------------------
#  Figure 2 - fig:stab_defects
# ----------------------------------------------------------------------------

function dist_2defects()
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    "dist" in names(df2) ? nothing : df2[!,:dist] .= 0
    for idx=1:nrow(df2)
        if df2[idx, :Ndefects] > 0
            print(idx, ", ")
            path = string("$idx/")
            dird5 = dir_fig*path*"Defects_pos/"
            last_file = readdir(dird5)[end]
            pts = load(string(dird5,last_file), "pos_unique")
            Ndef = load(string(dird5,last_file), "N")

            x = []
            for i=1:Ndef
                p = pts[i,:]
                p[2] >500 && p[2]<1000 ? push!(x, p[1]) : nothing
            end
            df2[idx, :dist] = length(x) > 0 ? abs(x[2]-x[1]) : 0
        end
    end
    CSV.write(joinpath(dir_df,"DF2.csv"), df2)
end

function heatmap_dist_2defects(kk)
    df3 = CSV.read(dir_df*"dist_min.csv", DataFrame)
    # kk = 2
    D = df3[df3.k .== kk, :Dmin]*1008*1e-2
    noNaND = filter(!isnan, D)
    f,a,h = heatmap(df3[df3.k .== kk, :r], df3[df3.k .== kk, :zr], D, colorrange=(minimum(noNaND),maximum(noNaND)), nan_color=(:grey,0.75))
    a.xlabel = L"\rho_0"
    a.ylabel = L"\zeta_\rho"
    a.xlabelsize = 30
    a.ylabelsize = 30
    a.xticklabelsize = 24
    a.yticklabelsize = 24
    a.title = kk == 1 ? L"\tau = 1" : L"\tau = 0.5"
    a.titlesize = 30
    Colorbar(f[:,2], label= L"d", labelsize = 30, ticklabelsize=24, limits=(minimum(noNaND),maximum(noNaND)))
    save(joinpath(dir_fig, "dist_kd=$kk.png"), f)
    
    dens = map(x->isnan(x) ? 0 : 1/(x*x), D/10.08)
    md = 500#maximum(dens)
    f,a,h = heatmap(df3[df3.k .== kk, :r], df3[df3.k .== kk, :zr], dens, colorrange=(0,md))#))#, nan_color=(:grey,0.75))
    a.xlabel = L"\rho_0"
    a.ylabel = L"\zeta_\rho"
    a.xlabelsize = 30
    a.ylabelsize = 30
    a.xticklabelsize = 24
    a.yticklabelsize = 24
    a.title = kk == 1 ? L"\tau = 1" : L"\tau = 0.5"
    a.titlesize = 30
    Colorbar(f[:,2], label= L"N^u", labelsize = 30, ticklabelsize=24, limits=(0,md))#maximum(dens)))
    save(joinpath(dir_fig, "upper_kd=$kk.png"), f)
end

function heatmap_dist_2defects_b(rr)
    df3 = CSV.read(dir_df*"dist_min.csv", DataFrame)
    D = df3[df3.r .== rr, :Dmin]*1008*1e-2
    k = df3[df3.r .== rr, :k]
    z = df3[df3.r .== rr, :zr]
    noNaND = filter(!isnan, D)
    f,a,h = heatmap(k, z, D, colorrange=(minimum(noNaND),maximum(noNaND)), nan_color=(:grey,0.75))
    a.xlabel = L"k_d"
    a.ylabel = L"\zeta_\rho"
    a.xlabelsize = 30
    a.ylabelsize = 30
    a.xticklabelsize = 24
    a.yticklabelsize = 24
    a.title = L"\rho_0 = %$rr"
    a.titlesize = 30
    a.xscale = Makie.pseudolog10
    Colorbar(f[:,2], label= L"d", labelsize = 30, ticklabelsize=24, limits=(minimum(noNaND),maximum(noNaND)))
    save(joinpath(dir_fig, "dist_r0=$rr.png"), f)
    
    dens = map(x->isnan(x) ? 0 : 1/(x*x), D/10.08)
    md = maximum(dens)
    # f,a,h = heatmap(log.(df3[df3.r .== rr, :k]), df3[df3.r .== rr, :zr], dens, colorrange=(0,md))#))#, nan_color=(:grey,0.75))
    f,a,h = heatmap(k, z, dens, colorrange=(0,md), colorscale = Makie.pseudolog10)#))#, nan_color=(:grey,0.75))
    a.xlabel = L"k_d"
    a.ylabel = L"\zeta_\rho"
    a.xlabelsize = 30
    a.ylabelsize = 30
    a.xticklabelsize = 24
    a.yticklabelsize = 24
    a.title = L"\rho_0 = %$rr"
    a.titlesize = 30
    a.xscale = Makie.pseudolog10
    for i in 1:length(dens)
        text!(a, "$(round(Int, dens[i]))", position = (k[i], z[i]), color = :black, align = (:center, :center))
    end
    Colorbar(f[:,2], label= L"N^u", labelsize = 30, ticklabelsize=24, limits=(0,md), scale = Makie.pseudolog10, ticks=[0,50,100,200,400,800,1600])#maximum(dens)))
    save(joinpath(dir_fig, "upper_r0=$rr.png"), f)
end

function velocity_defects(idx)
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    trajp = []
    trajn = []
    if df2[idx, :Ndefects] > 0
        path = string("$idx/")
        dird5 = dir_fig*path*"Defects_pos/"
        for i in readdir(dird5)
            pts = load(string(dird5,i), "pos_unique")
            if length(pts) > 0
                push!(trajp, pts[1,1])
                push!(trajn, pts[2,1])
            end
        end
    end
    return trajp, trajn
end

function plot_2defects_zoom(idx; dt = 1, sc=1)
    # set_theme!(theme_black())
    # IJulia.clear_output(true)
    # was pinned to "Z:/2defects/" and "D:/2defects_white/"
    path = string("$idx/")
    dirf = dir_fig*path*"Density_zoom_b/"
    dirf2 = dir_fig*path*"Angle_zoom/"
    # dirf3 = dir_fig*path*"Velocity_max/"
    dirf3 = dir_fig*path*"Velocity_zoom_b/"
    dirf4 = dir_fig*path*"Order_zoom/"
    # dirf5 = dir_fig*path*"Defects/"
    
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []
    mkpath(dirf)
    mkpath(dirf2)
    mkpath(dirf3)
    mkpath(dirf4)
    # mkpath(dirf5)
    num_files = length(files_names)
    num_pict = length(filter!(endswith(".png"), readdir(dirf)));
    println("Data in $dird")
    println("Figure will be saved in $dirf")
    println(num_files-num_pict, " / ", num_files) 
    
    if true #num_files > 0 && num_files > num_pict
        # try
        stepp = 5
        Nx, Nz = size(load(string(dird,files_names[1]), "rho"))
        sx = 200
        sz = 100
        x = 0:stepp:sx
        y = 0:stepp:sz
        hasP = haskey(load(string(dird,files_names[1])), "P")
        hasQ = haskey(load(string(dird,files_names[1])), "Q")
        # catch
        # end
        # x, y = meshgrid(1:stepp:Nz,1:stepp:Nx)
        sttopp = num_files # == 21 ? num_files : num_files-1
        # sttopp = (sttopp==0) ? sttopp+1 : sttopp
        for i=1:5#sttopp
            if true#!isfile(string(dirf,@sprintf("%04i",i),".png")) #|| hasQ
                # try
                t = floor(Int,parse(Int, files_names[i][5:end-4])*dt)
                    density = load(string(dird,files_names[i]), "rho")[400:600,450:550]

                    v_cpu = load(string(dird,files_names[i]), "v")[400:600,450:550,:]
                    max_v = sqrt(maximum(v_cpu[:,:,1].^2 .+ v_cpu[:,:,2].^2))
                    v_cpu ./= max_v
                    vx_cpu = v_cpu[1:stepp:end,1:stepp:end,1]
                    vz_cpu = v_cpu[1:stepp:end,1:stepp:end,2]


                    P_cpu = hasP ? load(string(dird,files_names[i]), "P")[400:600,450:550,:] : nothing
                    Px_cpu = hasP ? P_cpu[1:stepp:end,1:stepp:end,1] : nothing
                    Pz_cpu = hasP ? P_cpu[1:stepp:end,1:stepp:end,2] : nothing

                    Q_cpu = hasQ ? load(string(dird,files_names[i]), "Q")[400:600,450:550,:] : nothing
                    Q1_cpu = hasQ ? Q_cpu[1:stepp:end,1:stepp:end,1] : nothing
                    Q2_cpu = hasQ ? Q_cpu[1:stepp:end,1:stepp:end,2] : nothing
                    S = hasQ ? sqrt.(Q1_cpu.^2 .+ Q2_cpu.^2) : nothing
                    theta = hasQ ?  0.5*atan.(Q2_cpu.+1e-6, Q1_cpu) : (hasP ? atan.(Pz_cpu.+1e-6,Px_cpu) : 0.0*similar(vx_cpu) )
                    Qx = hasQ ? S.*cos.(theta) : nothing
                    Qz = hasQ ? S.*sin.(theta) : nothing


                    fig1 = Figure(size = (500*sc, 250*sc))
                        ax1 = Makie.Axis(fig1[1,1], aspect = DataAspect())
                        # ax1.title = string("Density and orientation " * string(t))
                        # h1 = heatmap!(ax1, 1:4:Nx, 1:4:Nz, density[1:4:end, 1:4:end], colormap=:RdPu)
                        h1 = heatmap!(ax1, density, colormap=:RdPu_8)
                        # try
                        #     Colorbar(fig1[1, 1], h1, ticklabelsize=15, size=20, vertical = false)
                        # catch
                        # end
                        if hasQ
                            arrow_size_Q = vec(S)
                            arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 10.0, linewidth=2.0)
                            # arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 8.0, linewidth=2.0)
                        end
                        if hasP
                            arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
                            arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*12.0, linecolor=:black, lengthscale = 4.0, linewidth=3.0, arrowhead = '▲', align = :center)
                            # arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P, linecolor=:black, lengthscale = 8.0, linewidth=2.0, arrowhead = '▲')
                        end
                        tightlimits!(ax1)
                        hideydecorations!(ax1)
                        hidexdecorations!(ax1)
                        # colgap!(fig1.layout, 4)
                        # rowgap!(fig1.layout, 0)
                    save(string(dirf,@sprintf("%04i",i),".png"), fig1)

                    # fig2 = Figure(size = (1000, 1025))
                    #     ax1 = Makie.Axis(fig2[2,1], aspect = DataAspect())
                    #     subs = hasQ ? " of Q " : " of P "
                    #     # ax1.title = string("Angle"* subs * string(t))
                    #     theta2 = hasQ ?  0.5*atan.(Q_cpu[:,:,2].+1e-6, Q_cpu[:,:,1]) : (hasP ? atan.(P_cpu[:,:,2].+1e-6,P_cpu[:,:,1]) : 0.0*similar(vx_cpu))
                    #     h1 = heatmap!(ax1, theta2, colormap=:hsv, colorrange = (-π, π))
                    #     try
                    #         Colorbar(fig2[1, 1], h1, ticklabelsize=15, size=20, vertical = false, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
                    #     catch
                    #     end
                    #     tightlimits!(ax1)
                    #     hideydecorations!(ax1)
                    #     hidexdecorations!(ax1)
                    #     colgap!(fig1.layout, 4)
                    #     rowgap!(fig1.layout, 0)
                    # save(string(dirf2,@sprintf("%04i",i),".png"), fig2)

                    fig3 = Figure(size = (500*sc, 250*sc))
                        ax1 = Makie.Axis(fig3[1,1], aspect = DataAspect())
                        # ax1.title = string("Velocity field " * string(t) * ", max v = " * string(max_v))
                        h1 = heatmap!(ax1, density, colormap=:RdPu_8)
                        arrow_size_v = vec(norm.(Vec2f.(vx_cpu,vz_cpu)))
                        arrows!(ax1, x, y, vx_cpu, vz_cpu, arrowsize = arrow_size_v*12.0, color=:white, lengthscale = 5.0, linewidth=3.0, arrowhead = '▲', align = :center) # :blue3
                        # try
                        #     Colorbar(fig3[1, 1], h1, ticklabelsize=15, size=20, vertical = false)
                        # catch
                        # end
                        tightlimits!(ax1)
                        hideydecorations!(ax1)
                        hidexdecorations!(ax1)
                        # colgap!(fig1.layout, 4)
                        # rowgap!(fig1.layout, 0)
                    save(string(dirf3,@sprintf("%04i",i),".png"), fig3)
                
                    # fig4 = Figure(size = (1000, 1025))
                    #     ax1 = Makie.Axis(fig4[2,1], aspect = DataAspect())
                    #     subs = hasQ ? " of Q " : " of P "
                    #     ax1.title = string("Order"* subs * string(t))
                    #     ordre = hasQ ? sqrt.(Q_cpu[:,:,1].^2 .+ Q_cpu[:,:,2].^2) : (hasP ? sqrt.(P_cpu[:,:,1].^2 .+ P_cpu[:,:,2].^2) : 0.0*similar(vx_cpu))
                    #     h1 = heatmap!(ax1, ordre, colorrange = (0, 2))
                    #     try
                    #         Colorbar(fig4[1, 1], h1, ticklabelsize=15, size=20, vertical = false) #, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
                    #     catch
                    #     end
                    #     if hasQ
                    #         arrow_size_Q = vec(S)
                    #         arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 10.0, linewidth=2.0)
                    #         # arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 8.0, linewidth=2.0)
                    #     end
                    #     if hasP
                    #         arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
                    #         arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
                    #     end
                    # save(string(dirf4,@sprintf("%04i",i),".png"), fig4)
                
                    # fig5 = Figure(size = (1000, 1000))
                    #     ax1 = Axis(fig5[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                    #     defect_idx = count_defect(P_cpu[:,:,1], P_cpu[:,:,2], 10)
                    #     s1 = scatter!(ax1, defect_idx)
                    #     Ndefect = length(defect_idx)
                    #     ax1.title = string("Defect position (number of defect = $Ndefect) " * string(i-1))
                    # save(string(dirf5,@sprintf("%04i",i),".png"), fig5)
                    # defect_idx = nothing
                
                # catch
                # end
                
            end
        end
        # pngstogif(dirf, dirf, string(idx), 1)
        # pngstogif(dirf2, dirf2, string(idx), 1)

    end
    if true #num_files > 0 && ( num_files > num_pict || !isfile(dirf*string(idx)*".gif") )
        # maxv = zeros(num_files)
        # for i=1:sttopp
        #     v_cpu = load(string(dird,files_names[i]), "v")
        #     vx_cpu = v_cpu[:,:,1]
        #     vz_cpu = v_cpu[:,:,2]
        #     maxv[i] = sqrt(maximum(vx_cpu.^2+vz_cpu.^2))
        # end
        # fig3 = Figure(size = (1000, 700))
        #     ax3 = Axis(fig3[1,1], xlabel = "x", ylabel = "max v")#, aspect = DataAspect())
        #     lines!(ax3, maxv)
        #     save(string(dirf3,string(idx),".png"), fig3)
        pngstogif(dirf, dirf, string(idx), 8)
        # pngstogif(dirf2, dirf2, string(idx), 8)
        length(readdir(dirf3)) > 0 ? pngstogif(dirf3, dirf3, string(idx), 8) : nothing
        # length(readdir(dirf4)) > 0 ? pngstogif(dirf4, dirf4, string(idx), 8) : nothing
        # length(readdir(dirf5)) > 0 ? pngstogif(dirf5, dirf5, string(idx), 4) : nothing
    end
    # hasP = haskey(load(string(dird,files_names[end])), "P")
    # hasQ = haskey(load(string(dird,files_names[end])), "Q")
    # hasP ? println("max|P| = ",maximum(abs.(load(string(dird,files_names[end]), "P")))) : nothing
    # hasQ ? println("max|Q| = ",maximum(abs.(load(string(dird,files_names[end]), "Q")))) : nothing
    # v_cpu = load(string(dird,files_names[end]), "v")
    # vx_cpu = v_cpu[:,:,1]
    # vz_cpu = v_cpu[:,:,2]
    # println("mav||v|| =", sqrt(maximum(vx_cpu.^2+vz_cpu.^2)) )
    # pngstogif(dirf, dir_fig*"/", string(1), 10, Nimg=101)
end

# ----------------------------------------------------------------------------
#  Figure 4 - fig:Ndefects
# ----------------------------------------------------------------------------

function make_heatmap(tkd = [10,5,1.0,0.2,0.1,0.01])
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    df3 = df2[in.(df2.kd, Ref(tkd)), :]
    # Nrho = length(unique(df2[:,:rho0]))
    # Nkd = length(unique(df2[:,:kd]))
    x1 = sort(unique(df3[:,:rho0])); Lx1 = length(unique(df3[:,:rho0]))
    x2 = sort(1 ./ unique(df3[:,:kd])); Lx2 = length(unique(df3[:,:kd]))
    y = sort(unique(df3[:,:zetarho]))
    Ny = length(unique(df3[:,:zetarho]))
    heatN = zeros(Lx1, Lx2, Ny)# .* NaN
    heatG4 = zeros(Lx1, Lx2, Ny)# .* NaN
    heatG3 = zeros(Lx1, Lx2, Ny)# .* NaN
    heatt = zeros(Lx1, Lx2, Ny)# .* NaN
    
    xcsv = []; ycsv = []; zcsv = []
    for i=1:Lx1, j=1:Ny, k=1:Lx2
        push!(xcsv, i)
        push!(ycsv, j)
        push!(zcsv, k)
    end
    Ncsv = length(xcsv)
    ndef = DataFrame(i = 1:Ncsv, x = xcsv, y = ycsv, z = zcsv, N = rand(Ncsv), G3 = rand(Ncsv), G4 = rand(Ncsv))
    
    for idx = 1:Ncsv#N_sim
        # for row = 1:Ncsv#N_sim
        # idx = df3[row, :fn]
        zr = df3[idx,:zetarho]
        r0 = df3[idx,:rho0]
        kd = df3[idx,:kd]
        i = findfirst(isequal(r0),x1) 
        k = findfirst(isequal(zr),y)
        j = findfirst(isequal(1/kd),x2)
        heatN[i, j, k] = df3[idx,:Ndefects]
        heatG4[i, j, k] = df3[idx,:Ndefects] > 50 ? df3[idx,:g4] : NaN
        heatG3[i, j, k] = df3[idx,:Ndefects] > 50 ? df3[idx,:g3] : NaN
        heatt[i, j, k] = df3[idx,:Ndefects] > 50 ? df3[idx,:theta4] : NaN
        
        l = ndef[(ndef.x .== i) .& (ndef.y .== k) .& (ndef.z .== j), :i][1]
        ndef[l, :N] = isnan(heatN[i, j, k]) ? -1 : heatN[i, j, k]
        ndef[l, :G3] = isnan(heatG3[i, j, k]) ? -1 : heatG3[i, j, k]
        ndef[l, :G4] = isnan(heatG4[i, j, k]) ? -1 : heatG4[i, j, k]
    end
    
    #### FIGURE NUMBER OF DEFECTS ####
#     fig = Figure(resolution=(2000,990))
#         axs = [Makie.Axis(fig[i,j]; aspect = DataAspect(), xticks = (1:2:12, [ string(x1[k]) for k in 1:2:Lx1]), yticks = (1:Ny, [ string(k) for k in y]), xaxisposition = i==1 ? :top : :bottom, yaxisposition = j==5 ? :right : :left, xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30) for i in 1:2 for j in 1:5] #, limits = (0, Lx1+0.5, -0, Ny+0.5)

#         for idx in 1:10
#             heatmap!(axs[idx], heatN[:,idx,:], colormap=:hot, nan_color = :black, colorrange=(0,400))
#             text!(axs[idx], 9, 13, text = L"\tau = %$(round(x2[idx], digits=1))", color = :white, font = :bold, fontsize = 24)#, strokewidth=2, strokecolor= :black)
#             # contour!(axs[idx], heatN[:,idx,:], levels=50:1:51, color=:white, linewidth=2)
#     end
#     axs[3].xlabel = L"\rho_0"
#     axs[8].xlabel = L"\rho_0"
#     axs[1].ylabel = L"\zeta_\rho"
#     axs[6].ylabel = L"\zeta_\rho"
#     hideydecorations!.(axs[2:5])
#     hideydecorations!.(axs[7:end])
#     colgap!(fig.layout, 1)
#     rowgap!(fig.layout, 1)
#     Colorbar(fig[:, end+1], colormap=:hot, limits=(0,400), label= "Number of defects", labelsize = 30, ticklabelsize=24)# L"N_\text{defects}")
    
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "Ndefect.png"), fig)
    
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","Ndefb.csv"), ndef)
    
#     #### FIGURE GAMMA 4 ####
#     fig2 = Figure(resolution=(2000,990))
#     axs2 = [Makie.Axis(fig2[i,j]; aspect = DataAspect(), xticks = (1:2:12, [ string(x1[k]) for k in 1:2:Lx1]), yticks = (1:Ny, [ string(k) for k in y]), xaxisposition = i==1 ? :top : :bottom, yaxisposition = j==5 ? :right : :left, xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30) for i in 1:2 for j in 1:5]

#     for idx in 1:10
#         heatmap!(axs2[idx], replace(heatG4[:,idx,:], NaN=>0.0), colormap=:vik, nan_color = :black, colorrange=(0,1))
#         text!(axs2[idx], 9, 13, text = L"\tau = %$(round(x2[idx], digits=1))", color = :white, font = :bold, fontsize = 24)#, strokewidth=2, strokecolor= :black)
#     end
#     axs2[3].xlabel = L"\rho_0"
#     axs2[8].xlabel = L"\rho_0"
#     axs2[1].ylabel = L"\zeta_\rho"
#     axs2[6].ylabel = L"\zeta_\rho"
#     hideydecorations!.(axs2[2:5])
#     hideydecorations!.(axs2[7:end])
#     colgap!(fig2.layout, 1)
#     rowgap!(fig2.layout, 1)
#     Colorbar(fig2[:, end+1], colormap=:vik, limits=(0,1), label= L"\langle|\gamma_4|\rangle", labelsize = 30, ticklabelsize=24)# L"N_\text{defects}")
    
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "Gamma4.png"), fig2)
    
#     #### FIGURE GAMMA 3 ####
#     fig3 = Figure(resolution=(2000,990))
#     axs3 = [Makie.Axis(fig3[i,j]; aspect = DataAspect(), xticks = (1:2:12, [ string(x1[k]) for k in 1:2:Lx1]), yticks = (1:Ny, [ string(k) for k in y]), xaxisposition = i==1 ? :top : :bottom, yaxisposition = j==5 ? :right : :left, xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30) for i in 1:2 for j in 1:5]

#     for idx in 1:10
#         heatmap!(axs3[idx], replace(heatG3[:,idx,:], NaN=>0.0), colormap=:vik, nan_color = :black, colorrange=(0,1))
#         text!(axs3[idx], 9, 13, text = L"\tau = %$(round(x2[idx], digits=1))", color = :white, font = :bold, fontsize = 24)#, strokewidth=2, strokecolor= :black)
#     end
#     axs3[3].xlabel = L"\rho_0"
#     axs3[8].xlabel = L"\rho_0"
#     axs3[1].ylabel = L"\zeta_\rho"
#     axs3[6].ylabel = L"\zeta_\rho"
#     hideydecorations!.(axs3[2:5])
#     hideydecorations!.(axs3[7:end])
#     colgap!(fig3.layout, 1)
#     rowgap!(fig3.layout, 1)
#     Colorbar(fig3[:, end+1], colormap=:vik, limits=(0,1), label= L"\langle|\gamma_3|\rangle", labelsize = 30, ticklabelsize=24)# L"N_\text{defects}")
    
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "Gamma3.png"), fig3)
#     return fig, fig2, fig3
end


function Ndefects_as_act(rho0=[0.6, 0.9, 1.0, 1.1], kd=0.6)
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    x2 = sort(unique(rho0)); Lx2 = length(x2)
    x1 = sort(unique(df2[:,:zetarho])); Lx1 = length(x1)
    yn = zeros(Lx1, Lx2)# .* NaN
    yG4 = zeros(Lx1, Lx2)# .* NaN
    yG3 = zeros(Lx1, Lx2)# .* NaN
    for idx = 1:N_sim
        zr = df2[idx,:zetarho]
        r0 = df2[idx,:rho0]
        kdi = df2[idx,:kd]
        if r0 in rho0 && kdi == kd
            j = findfirst(isequal(r0),x2) 
            i = findfirst(isequal(zr),x1)
            yn[i, j] = df2[idx,:Ndefects]
            yG4[i, j] = df2[idx,:Ndefects] > 25 ? df2[idx,:g4] : 0.0
            yG3[i, j] = df2[idx,:Ndefects] > 25 ? df2[idx,:g3] : 0.0
        end
    end
    
    colors = Makie.wong_colors()[1:Lx2]
    
    #### FIGURE NUMBER OF DEFECTS ####
    f1 = Figure(resolution=(1000,800))
    a1 = Makie.Axis(f1[1,1], xlabel = L"\zeta_\rho", ylabel="Number of defects", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    for i = 1:Lx2
        scatter!(a1, x1, yn[:,i], markersize=15, color=colors[i], label=L"\rho_0=%$(x2[i])")
        lines!(a1, x1, yn[:,i], linewidth=4, color=colors[i], linestyle=:dash)
    end
    axislegend(; position = :rt, labelsize = 24)
    f1
    mkpath(dir_fig)
    save(joinpath(dir_fig, "As_zetarho_Ndefects.png"), f1)
    
    #### FIGURE GAMMA 3 & 4 ####
    my = maximum(yn)
    f2 = Figure(resolution=(1000,800))
    a2 = Makie.Axis(f2[1,1], xlabel = L"\zeta_\rho", ylabel=L"\langle|\gamma_p|\rangle", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    for i = 1:Lx2
        scatter!(a2, x1, yG4[:,i], linewidth=4, color=colors[i], marker = :rect, markersize = 15)
        lines!(a2, x1, yG4[:,i], linewidth=4, color=yn[:,i]/my, colormap=[(colors[i], 0.0), (colors[i], 1.0)], linestyle = :dot, colorrange=(0,1))
        scatter!(a2, x1, yG3[:,i], linewidth=4, color=colors[i], marker = :utriangle, markersize = 15)
        lines!(a2, x1, yG3[:,i], linewidth=4, color=yn[:,i]/my, colormap=[(colors[i], 0.0), (colors[i], 1.0)], linestyle = :dash, colorrange=(0,1))#, 
    end
    markertype = [:rect, :utriangle]
    linstyle = [:dot, :dash]
    group_type = [[LineElement(color = :black, linestyle = linstyle[mt]), MarkerElement(marker = markertype[mt], color = :black, strokecolor = :transparent,  markersize = 15)] for mt in 1:length(markertype)]
    group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors]
    Legend(f2[1,1], [group_type, group_color], [[L"p=4", L"p=3"], [L"\rho_0=%$(i)" for i in x2]], labelsize = 24, nbanks = 2, orientation = :horizontal, tellheight = false, tellwidth=false, [L"p", L"\rho_0"], titleposition = :top, titlesize=24, halign = :right, valign = :top, margin = (10, 10, 10, 10))#, ) 
    Colorbar(f2[1, 2], limits = (0,1), colormap = [(:black, 0.0), (:black, 1.0)], label= "Normalized number of defects", labelsize = 24, ticklabelsize=20)
    f2
    mkpath(dir_fig)
    save(joinpath(dir_fig, "As_zetarho_Gamma34.png"), f2)
    return f1, f2
end

function Ndefects_as_rho0(zetarho=[0.0, 1.0, 4.0, 8.0, 12.0], kd=0.6, name_file_addon = "")
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    x2 = sort(unique(zetarho)); Lx2 = length(x2)
    x1 = sort(unique(df2[:,:rho0])); Lx1 = length(x1)
    yn = zeros(Lx1, Lx2)# .* NaN
    yG4 = zeros(Lx1, Lx2)# .* NaN
    yG3 = zeros(Lx1, Lx2)# .* NaN
    for idx = 1:N_sim
        zr = df2[idx,:zetarho]
        r0 = df2[idx,:rho0]
        kdi = df2[idx,:kd]
        if zr in zetarho && kdi == kd
            j = findfirst(isequal(zr),x2) 
            i = findfirst(isequal(r0),x1)
            yn[i, j] = df2[idx,:Ndefects]
            yG4[i, j] = df2[idx,:Ndefects] > 25 ? df2[idx,:g4] : 0.0
            yG3[i, j] = df2[idx,:Ndefects] > 25 ? df2[idx,:g3] : 0.0
        end
    end
    
    colors = Makie.wong_colors()[1:Lx2]
    
    #### FIGURE NUMBER OF DEFECTS ####
    f1 = Figure(resolution=(1000,800))
    a1 = Makie.Axis(f1[1,1], xlabel = L"\rho_0", ylabel="Number of defects", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    for i = 1:Lx2
        scatter!(a1, x1, yn[:,i], markersize=15, color=colors[i], label=L"\zeta_\rho=%$(x2[i])")
        lines!(a1, x1, yn[:,i], linewidth=4, color=colors[i], linestyle=:dash)
    end
    axislegend(; position = :rt, labelsize = 24)
    mkpath(dir_fig)
    save(joinpath(dir_fig, string("As_rho0_Ndefects", name_file_addon,".png"), f1)
    
    #### FIGURE GAMMA 3 & 4 ####
    my = maximum(yn)
    f2 = Figure(resolution=(1000,800))
    a2 = Makie.Axis(f2[1,1], xlabel = L"\rho_0", ylabel=L"\langle|\gamma_p|\rangle", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    for i = 1:Lx2
        scatter!(a2, x1, yG4[:,i], linewidth=4, color=[(colors[i], 1.0 - Int(yG4[k,i]==0.0)) for k = 1:Lx1], marker = :rect, markersize = 15)#yn[:,i]/maximum(yn[:,i])*20)#15)
        # lines!(a2, yG4[:,i], linewidth=4, color=colors[i], linestyle = :dot)#, alpha = yn[:,i]/maximum(yn[:,i]))#, linestyle = :rect
        lines!(a2, x1, yG4[:,i], linewidth=4, color=yn[:,i]/my, colormap=[(colors[i], 0.0), (colors[i], 1.0)], linestyle = :dot, colorrange=(0,1))#, alpha = yn[:,i]/maximum(yn[:,i]))#, linestyle = :rect
        scatter!(a2, x1, yG3[:,i], linewidth=4, color=[(colors[i], 1.0 - Int(yG4[k,i]==0.0)) for k = 1:Lx1], marker = :utriangle, markersize = 15)#yn[:,i]/maximum(yn[:,i])*20)#(1 .-yn[:,i]/maximum(yn[:,i]))*15)#15)
        # lines!(a2, yG3[:,i], linewidth=4, color=colors[i], linestyle = :dash)#, 
        lines!(a2, x1, yG3[:,i], linewidth=4, color=yn[:,i]/my, colormap=[(colors[i], 0.0), (colors[i], 1.0)], linestyle = :dash, colorrange=(0,1))#, 
    end
    markertype = [:rect, :utriangle]
    linstyle = [:dot, :dash]
    group_type = [[LineElement(color = :black, linestyle = linstyle[mt]), MarkerElement(marker = markertype[mt], color = :black, strokecolor = :transparent,  markersize = 15)] for mt in 1:length(markertype)]
    group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors]
    Legend(f2[1,1], [group_type, group_color], [[L"p=4", L"p=3"], [L"\zeta_\rho=%$(i)" for i in x2]], labelsize = 24, nbanks = 2, orientation = :horizontal, tellheight = false, tellwidth=false, [L"p", L"\zeta_\rho"], titleposition = :top, titlesize=24, halign = :right, valign = :bottom, margin = (10, 10, 10, 10))#, ) 
    Colorbar(f2[1, 2], limits = (0,1), colormap = [(:black, 0.0), (:black, 1.0)], label= "Normalized number of defects", labelsize = 24, ticklabelsize=20)
    mkpath(dir_fig)
    save(joinpath(dir_fig, string("As_rho0_Gamma34",name_file_addon,".png"), f2)
    
    #### FIGURE GAMMA 3 & 4 ####
    # my = maximum(yn)
    mg = maximum([maximum(yG4)-G_avg[3],maximum(yG3)-G_avg[2]])
    mg += 0.05
    f2b = Figure(resolution=(1000,800))
    a2b = Makie.Axis(f2b[1,1], xlabel = L"\rho_0", ylabel=L"\langle|\gamma_p|\rangle-\gamma_P^0", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    for i = 1:Lx2
        scatter!(a2b, x1, yG4[:,i].-G_avg[3], linewidth=4, color=[(colors[i], 1.0 - Int(yG4[k,i]==0.0)) for k = 1:Lx1], marker = :rect, markersize = 15)#yn[:,i]/maximum(yn[:,i])*20)#15)
        # lines!(a2, yG4[:,i], linewidth=4, color=colors[i], linestyle = :dot)#, alpha = yn[:,i]/maximum(yn[:,i]))#, linestyle = :rect
        lines!(a2b, x1, yG4[:,i].-G_avg[3], linewidth=4, color=yn[:,i]/my, colormap=[(colors[i], 0.0), (colors[i], 1.0)], linestyle = :dot, colorrange=(0,1))#, alpha = yn[:,i]/maximum(yn[:,i]))#, linestyle = :rect
        scatter!(a2b, x1, yG3[:,i].-G_avg[2], linewidth=4, color=[(colors[i], 1.0 - Int(yG4[k,i]==0.0)) for k = 1:Lx1], marker = :utriangle, markersize = 15)#yn[:,i]/maximum(yn[:,i])*20)#(1 .-yn[:,i]/maximum(yn[:,i]))*15)#15)
        # lines!(a2, yG3[:,i], linewidth=4, color=colors[i], linestyle = :dash)#, 
        lines!(a2b, x1, yG3[:,i].-G_avg[2], linewidth=4, color=yn[:,i]/my, colormap=[(colors[i], 0.0), (colors[i], 1.0)], linestyle = :dash, colorrange=(0,1))#, 
    end
    markertype = [:rect, :utriangle]
    linstyle = [:dot, :dash]
    group_type = [[LineElement(color = :black, linestyle = linstyle[mt]), MarkerElement(marker = markertype[mt], color = :black, strokecolor = :transparent,  markersize = 15)] for mt in 1:length(markertype)]
    group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors]
    Legend(f2b[1,1], [group_type, group_color], [[L"p=4", L"p=3"], [L"\zeta_\rho=%$(i)" for i in x2]], labelsize = 24, nbanks = 2, orientation = :horizontal, tellheight = false, tellwidth=false, [L"p", L"\zeta_\rho"], titleposition = :top, titlesize=24, halign = :right, valign = :bottom, margin = (10, 10, 10, 10))#, ) 
    Colorbar(f2b[1, 2], limits = (0,1), colormap = [(:black, 0.0), (:black, 1.0)], label= "Normalized number of defects", labelsize = 24, ticklabelsize=20)
    ylims!(a2b, -mg, mg)
    mkpath(dir_fig)
    save(joinpath(dir_fig, string("As_rho0_Gamma34_mean",name_file_addon,".png"), f2b)
    
    return f1, f2, f2b
end

# ----------------------------------------------------------------------------
#  Figure 6 - fig:states
# ----------------------------------------------------------------------------

function make_full_heatmap_idx(tkd=[0.2], tr0=[0.4,0.5,0.6,0.65,0.7,0.75,0.8,1.2], tzr=[4]; part = 1.0, Lz = 3, Lx = 3)
    # set_theme!(empty_theme)
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    
    for i = 1:length(tidx)
        idx = tidx[i]
        dird = dir_df*"$idx/"*"Data/"
        rho = load(string(dird,readdir(dird)[end]), "rho")
        # mrho = maximum(rho)
        # rhomax = mrho > rhomax ? mrho : rhomax
        # mrho = minimum(rho)
        # rhomin = mrho < rhomin ? mrho : rhomin
        rhomax[i] = maximum(rho)
        rhomin[i] = minimum(rho)
    end
    
    fig = Figure(resolution=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lx for j in 1:Lz]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        zr = df2[idx,:zetarho]
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        last_file = readdir(dird)[end]
        rho = load(string(dird,last_file), "rho")
        Nx, Nz = size(rho)
        pNx = round(Int, Nx*part)
        pNz = round(Int, Nz*part)
        x = 1:steppp:pNx
        y = 1:steppp:pNz
        P_cpu = load(string(dird,last_file), "P")
        Px_cpu = P_cpu[x,y,1]
        Pz_cpu = P_cpu[x,y,2]
        arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))

        heatmap!(axs[i], rho[1:pNx,1:pNz], colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        arrows!(axs[i], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
        # Box(gls[i][1,1], color = (:white, 0.75), strokecolor = :black, strokewidth = 1)
        # Label(gls[i][1,1], L"\zeta_\rho = %$zr", padding = (10, 10, 10, 10), color = :black, font = :bold, fontsize = 30)
    end
    tightlimits!.(axs)
    hideydecorations!.(axs)
    hidexdecorations!.(axs)
    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 1)
    Colorbar(fig[:, end+1], colormap=:RdPu_8, colorrange=(rhomin[1],rhomax[1]), ticksvisible=false, size = 25, ticklabelsize=0) #, label= L"\rho", labelsize = 30, ticklabelsize=24)
    mkpath(dir_fig)
    save(joinpath(dir_fig, "phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_v2.png"), fig)
    return fig
end

function make_zoom_heatmap_idx(tkd=[0.2], tr0=[0.4,0.7, 0.75,0.8], tzr=[4]; part = 0.25, Lz = 2, Lx = 2, posx=400, posz=750)
    # set_theme!(empty_theme)
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    for i = 1:length(tidx)
        idx = tidx[i]
        dird = dir_df*"$idx/"*"Data/"
        rho = load(string(dird,readdir(dird)[end]), "rho")
        # mrho = maximum(rho)
        # rhomax = mrho > rhomax ? mrho : rhomax
        # mrho = minimum(rho)
        # rhomin = mrho < rhomin ? mrho : rhomin
        rhomax[i] = maximum(rho)
        rhomin[i] = minimum(rho)
    end
    
    fig = Figure(resolution=(800,800))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lx for j in 1:Lz]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        zr = df2[idx,:zetarho]
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        last_file = readdir(dird)[end]
        rho = load(string(dird,last_file), "rho")
        Nx, Nz = size(rho)
        pNx = round(Int, Nx*part)
        pNz = round(Int, Nz*part)
        x = posx+10:stepp:posx+pNx-10
        y = posz+10:stepp:posz+pNz-10
        xx = 10:stepp:pNx-10
        yy = 10:stepp:pNz-10
        P_cpu = load(string(dird,last_file), "P")
        Px_cpu = P_cpu[x,y,1]
        Pz_cpu = P_cpu[x,y,2]
        arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))

        heatmap!(axs[i], rho[posx:posx+pNx,posz:posz+pNz], colormap=:RdPu_8, colorrange=(rhomin[i],rhomax[i]))
        arrows!(axs[i], xx, yy, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
        # Box(gls[i][1,1], color = (:white, 0.75), strokecolor = :black, strokewidth = 1)
        # Label(gls[i][1,1], L"\zeta_\rho = %$zr", padding = (10, 10, 10, 10), color = :black, font = :bold, fontsize = 30)
    end
    tightlimits!.(axs)
    hideydecorations!.(axs)
    hidexdecorations!.(axs)
    # Colorbar(fig[:, end+1], colormap=:viridis, colorrange=(rhomin[1],rhomax[1]), ticksvisible=true, size = 25, ticklabelsize=0) #, label= L"\rho", labelsize = 30, ticklabelsize=24)
    colgap!(fig.layout, 1)
    rowgap!(fig.layout, 0.1)
    mkpath(dir_fig)
    save(joinpath(dir_fig, "phases_zoom_$(tkd[1])_$(tr0[1])_$(tzr[1])_v2.png"), fig)
    return fig
end

# ----------------------------------------------------------------------------
#  Figure 9 - fig:gamma_rho
# ----------------------------------------------------------------------------

function gamma_as_Ndef(stepp=30)
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    Ndef=Array(25:stepp:400)
    Gammap = zeros(length(Ndef), 5, 2)
    for i=1:length(Ndef)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g2]
        Gammap[i,1,1] = mean(t)
        Gammap[i,1,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g3]
        Gammap[i,2,1] = mean(t)
        Gammap[i,2,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g4]
        Gammap[i,3,1] = mean(t)
        Gammap[i,3,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g5]
        Gammap[i,4,1] = mean(t)
        Gammap[i,4,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g6]
        Gammap[i,5,1] = mean(t)
        Gammap[i,5,2] = var(t)
    end
    
    fig = Figure(size = (1000, 800))
    ax = Makie.Axis(fig[1,1], xlabel="Number of defects", ylabel=L"\langle|\gamma_p|\rangle", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    errorbars!(ax, Ndef, Gammap[:,1,1], Gammap[:,1,2], whiskerwidth = 12, color=:violet)
    scatter!(ax, Ndef, Gammap[:,1,1], color=:violet, makersize=8, label=L"p=2")
    errorbars!(ax, Ndef, Gammap[:,2,1], Gammap[:,2,2], whiskerwidth = 12, color=:red)
    scatter!(ax, Ndef, Gammap[:,2,1], color=:red, makersize=8, label=L"p=3")
    errorbars!(ax, Ndef, Gammap[:,3,1], Gammap[:,3,2], whiskerwidth = 12, color=:black)
    scatter!(ax, Ndef, Gammap[:,3,1], color=:black, makersize=8, label=L"p=4")
    errorbars!(ax, Ndef, Gammap[:,4,1], Gammap[:,4,2], whiskerwidth = 12, color=:green)
    scatter!(ax, Ndef, Gammap[:,4,1], color=:green, makersize=8, label=L"p=5")
    errorbars!(ax, Ndef, Gammap[:,5,1], Gammap[:,5,2], whiskerwidth = 12, color=:blue)
    scatter!(ax, Ndef, Gammap[:,5,1], color=:blue, makersize=8, label=L"p=6")
    axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
    mkpath(dir_fig)
    save(joinpath(dir_fig, "GammapvsNdef.png"), fig)
    fig2 = Figure(size = (1000, 800))
    ax2 = Makie.Axis(fig2[1,1], xlabel="Number of defects", ylabel=L"\langle|\gamma_p|\rangle-\gamma_P^0", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    errorbars!(ax2, Ndef, Gammap[:,1,1].-G_avg[1], Gammap[:,1,2], whiskerwidth = 12, color=:violet)
    scatter!(ax2, Ndef, Gammap[:,1,1].-G_avg[1], color=:violet, makersize=8, label=L"p=2")
    errorbars!(ax2, Ndef, Gammap[:,2,1].-G_avg[2], Gammap[:,2,2], whiskerwidth = 12, color=:red)
    scatter!(ax2, Ndef, Gammap[:,2,1].-G_avg[2], color=:red, makersize=8, label=L"p=3")
    errorbars!(ax2, Ndef, Gammap[:,3,1].-G_avg[3], Gammap[:,3,2], whiskerwidth = 12, color=:black)
    scatter!(ax2, Ndef, Gammap[:,3,1].-G_avg[3], color=:black, makersize=8, label=L"p=4")
    errorbars!(ax2, Ndef, Gammap[:,4,1].-G_avg[4], Gammap[:,4,2], whiskerwidth = 12, color=:green)
    scatter!(ax2, Ndef, Gammap[:,4,1].-G_avg[4], color=:green, makersize=8, label=L"p=5")
    errorbars!(ax2, Ndef, Gammap[:,5,1].-G_avg[5], Gammap[:,5,2], whiskerwidth = 12, color=:blue)
    scatter!(ax2, Ndef, Gammap[:,5,1].-G_avg[5], color=:blue, makersize=8, label=L"p=6")
    axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
    mkpath(dir_fig)
    save(joinpath(dir_fig, "GammapvsNdef2.png"), fig2)
    
    gamma_Ndef = DataFrame(x=Ndef)
    for p=1:5
        gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
        gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
    end
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","gamma_Ndef.csv"), gamma_Ndef)
    
    Ndefects = df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])
    g2 = df2[:,:g2] .* .!isnan.(df2[:,:Ndefects]) 
    g3 = df2[:,:g3] .* .!isnan.(df2[:,:Ndefects]) 
    g4 = df2[:,:g4] .* .!isnan.(df2[:,:Ndefects]) 
    g5 = df2[:,:g5] .* .!isnan.(df2[:,:Ndefects]) 
    g6 = df2[:,:g6] .* .!isnan.(df2[:,:Ndefects]) 
    fig3 = Figure(size = (1000, 800))
    ax3 = Makie.Axis(fig3[1,1], xlabel="Number of defects", ylabel=L"\langle|\gamma_p|\rangle-\gamma_P^0", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    xlims!(ax3, 5, 380)
    # scatter!(ax3, Ndefects, g2.-G_avg[1], color=:violet, makersize=8, label=L"p=2")
    scatter!(ax3, Ndefects, g3.-G_avg[2], color=:red, makersize=8, label=L"p=3")
    scatter!(ax3, Ndefects, g4.-G_avg[3], color=:black, makersize=8, label=L"p=4")
    # scatter!(ax3, Ndefects, g5.-G_avg[4], color=:green, makersize=8, label=L"p=5")
    # scatter!(ax3, Ndefects, g6.-G_avg[5], color=:blue, makersize=8, label=L"p=6")
    axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
    mkpath(dir_fig)
    save(joinpath(dir_fig, "GammapvsNdef3.png"), fig3)
    
    g34_Ndef = DataFrame(x = Ndefects, y1=g3.-G_avg[2], y2=g4.-G_avg[3])
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","g34_Ndef.csv"), g34_Ndef)
return fig, fig2, fig3
end

function gamma_as_Ndefpm(stepp=30)
    # Reads the run set configured through DATA_DIR. This used to be pinned to
    # "E:/AdptDt/", which silently ignored it.
    df2 = CSV.read(dir_df*"DF2-pm.csv", DataFrame)
    r = df2[:,:rho0]
    g4 = df2[:,:g4pm]
    g6 = df2[:,:g6pm]
    
    ru = unique(r)
    res = zeros(length(ru), 4)
    for i in eachindex(ru)
        idx_range = findall( x-> (df2[x,:rho0] == ru[i]) & (df2[x,:Ndefects] > 150), 1:nrow(df2));
        res[i,1] = mean(df2[idx_range, :g4pm])-G_avg[3]
        res[i,2] = std(df2[idx_range, :g4pm])
        res[i,3] = mean(df2[idx_range, :g6pm])-G_avg[5]
        res[i,4] = std(df2[idx_range, :g6pm])
    end
    g34_rho = DataFrame(x = ru, gq = res[:,1], gq_sd = res[:,2], gh = res[:,3], gh_sd = res[:,4])
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","g34_rho-pm.csv"), g34_rho)

    f, a, l = scatter(ru, res[:,1] .- G_avg[3]);
    scatter!(a, ru, res[:,2] .- G_avg[5])
    f

#     gamma_Ndef = DataFrame(x=Ndef)
#     for p=1:5
#         gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
#         gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
#     end
#     CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","gamma_Ndef-pm.csv"), gamma_Ndef)
    
#     Ndefects = df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])
#     g2 = df2[:,:g2pm] .* .!isnan.(df2[:,:Ndefects]) 
#     g3 = df2[:,:g3pm] .* .!isnan.(df2[:,:Ndefects]) 
#     g4 = df2[:,:g4pm] .* .!isnan.(df2[:,:Ndefects]) 
#     g5 = df2[:,:g5pm] .* .!isnan.(df2[:,:Ndefects]) 
#     g6 = df2[:,:g6pm] .* .!isnan.(df2[:,:Ndefects]) 
#     fig3 = Figure(size = (1000, 800))
#     ax3 = Makie.Axis(fig3[1,1], xlabel="Number of defects", ylabel=L"\langle|\gamma_p|\rangle-\gamma_P^0", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
#     xlims!(ax3, 5, 380)
#     # scatter!(ax3, Ndefects, g2.-G_avg[1], color=:violet, makersize=8, label=L"p=2")
#     scatter!(ax3, Ndefects, g4.-G_avg[3], color=:black, makersize=8, label=L"p=4")
#     scatter!(ax3, Ndefects, g6.-G_avg[5], color=:blue, makersize=8, label=L"p=6")
#     # scatter!(ax3, Ndefects, g5.-G_avg[4], color=:green, makersize=8, label=L"p=5")
#     # scatter!(ax3, Ndefects, g6.-G_avg[5], color=:blue, makersize=8, label=L"p=6")
#     axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "GammapvsNdef3-pm.png"), fig3)
    
#     g34_Ndef = DataFrame(x = Ndefects, y1=g4.-G_avg[3], y2=g6.-G_avg[5])
#     CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","g34_Ndef-pm.csv"), g34_Ndef)
# return fig, fig2, fig3
end

function gamma_as_ddef(lidx = N_sim)
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)[1:lidx,:] 
    ddef = zeros(lidx)
    dird = dir_df*"1/Data/"
    Nx, Ny = size(load(string(dird,readdir(dird)[end]), "rho"))
    Ndef =  df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])
    r =  df2[:,:Rad] .* .!isnan.(df2[:,:Rad])
    @. ddef = Ndef*pi*r^2/(Nx*Ny)
    g2 = df2[:,:g2] .* .!isnan.(df2[:,:Ndefects]) 
    g3 = df2[:,:g3] .* .!isnan.(df2[:,:Ndefects]) 
    g4 = df2[:,:g4] .* .!isnan.(df2[:,:Ndefects]) 
    g5 = df2[:,:g5] .* .!isnan.(df2[:,:Ndefects]) 
    g6 = df2[:,:g6] .* .!isnan.(df2[:,:Ndefects]) 
    
    fig = Figure(size = (1000, 800))
    ax = Makie.Axis(fig[1,1], xlabel="Defects <density> ", ylabel=L"\langle|\gamma_p|\rangle", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    # errorbars!(ax, Ndef, Gammap[:,1,1], Gammap[:,1,2], whiskerwidth = 12, color=:violet)
    scatter!(ax, ddef, g2, color=:violet, makersize=8, label=L"p=2")
    # errorbars!(ax, Ndef, Gammap[:,2,1], Gammap[:,2,2], whiskerwidth = 12, color=:red)
    scatter!(ax, ddef, g3, color=:red, makersize=8, label=L"p=3")
    # errorbars!(ax, Ndef, Gammap[:,3,1], Gammap[:,3,2], whiskerwidth = 12, color=:black)
    scatter!(ax, ddef, g4, color=:black, makersize=8, label=L"p=4")
    # errorbars!(ax, Ndef, Gammap[:,4,1], Gammap[:,4,2], whiskerwidth = 12, color=:green)
    scatter!(ax, ddef, g5, color=:green, makersize=8, label=L"p=5")
    # errorbars!(ax, Ndef, Gammap[:,5,1], Gammap[:,5,2], whiskerwidth = 12, color=:blue)
    scatter!(ax, ddef, g6, color=:blue, makersize=8, label=L"p=6")
    axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
    # mkpath(dir_fig)
    # save(joinpath(dir_fig, "GammapvsNdef.png"), fig)
    fig2 = Figure(size = (1000, 800))
    ax2 = Makie.Axis(fig2[1,1], xlabel="Defects <density> ", ylabel=L"\langle|\gamma_p|\rangle", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    xlims!(ax2, 0,0.3)
    # errorbars!(ax, Ndef, Gammap[:,1,1], Gammap[:,1,2], whiskerwidth = 12, color=:violet)
    scatter!(ax2, ddef, g2.-G_avg[1], color=:violet, makersize=8, label=L"p=2")
    # errorbars!(ax, Ndef, Gammap[:,2,1], Gammap[:,2,2], whiskerwidth = 12, color=:red)
    scatter!(ax2, ddef, g3.-G_avg[2], color=:red, makersize=8, label=L"p=3")
    # errorbars!(ax, Ndef, Gammap[:,3,1], Gammap[:,3,2], whiskerwidth = 12, color=:black)
    scatter!(ax2, ddef, g4.-G_avg[3], color=:black, makersize=8, label=L"p=4")
    # errorbars!(ax, Ndef, Gammap[:,4,1], Gammap[:,4,2], whiskerwidth = 12, color=:green)
    scatter!(ax2, ddef, g5.-G_avg[4], color=:green, makersize=8, label=L"p=5")
    # errorbars!(ax, Ndef, Gammap[:,5,1], Gammap[:,5,2], whiskerwidth = 12, color=:blue)
    scatter!(ax2, ddef, g6.-G_avg[5], color=:blue, makersize=8, label=L"p=6")
    axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
    # mkpath(dir_fig)
    # save(joinpath(dir_fig, "GammapvsNdef.png"), fig)
return fig, fig2
end

# ----------------------------------------------------------------------------
#  Figure 12 - fig:lattice1extended
# ----------------------------------------------------------------------------

function figure_triple(idx, k=2)
    mkpath(joinpath(dir_fig, "triple_$idx_$k/"))
    ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    lfile_path = dird*lfile_name
    num_files = length(files_names)
       
    dird5 = dir_fig*path*"Defects_pos/"
    files_names5 = isdir(dird5) ? readdir(dird5) : []
    lfile_name5 = last(files_names5)
    lfile_path5 = dird5*lfile_name5
    
      
    density = load(string(dird,lfile_name), "rho")

    Nx, Nz = size(density)
    x = 1:stepp:Nx
    y = 1:stepp:Nz
    
    P_cpu = load(string(dird,lfile_name), "P")
    Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    
    fig1 = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density)
    Colorbar(fig1[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
    arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1.layout, 4)
    rowgap!(fig1.layout, 0)
    save(joinpath(dir_fig, "triple_$idx_$k/density.png"), fig1)
    # fig1
    
    
    pts = load(lfile_path5, "pos")
    color = load(lfile_path5, "colors")
    tri = triangulate(pts')
    vorn = voronoi(tri, true)
    map_points = vorn.generators
    fig2 = Figure(size = (1000, 1025))
    ax = Makie.Axis(fig2[2,1], limits=(ex, ex+Nx, ex, ex+Nz), aspect = DataAspect())
    sc = voronoiplot!(ax, vorn, strokewidth=4, strokecolor=:black)
    s1 = scatter!(ax, pts', color=color, markersize = 15)
    # elem_1 = MarkerElement(color = :blue, marker = :circle, markersize = 25)
    # elem_2 = MarkerElement(color = :red, marker = :circle, markersize = 25)
    elem_1 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    elem_2 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
    hideydecorations!(ax)
    hidexdecorations!(ax)
    colgap!(fig2.layout, 0)
    rowgap!(fig2.layout, -5)
    tightlimits!(ax)
    save(joinpath(dir_fig, "triple_$idx_$k/voronoi.png"), fig2)
    # fig2
    
    mean_gamma_t = zeros(num_files,5)
    mean_theta_t = zeros(num_files,5)
    for file in 1:num_files
        pts = Array{Float64,2}
        Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
        pts = load(string(dird5,@sprintf("%04i",file),".jld"), "pos")

        if length(pts) > 50
            tri = triangulate(pts')
            vorn = voronoi(tri, true)
            map_points = vorn.generators
            polys = vorn.polygons
            boundary = vorn.boundary_polygons
            Npoly = length(polys)

            rθ = zeros(Npoly, 30, 2)
            for i=1:Npoly
                if !(i in boundary)
                    c = DelaunayTriangulation.get_centroid(vorn, i)
                    if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
                        po = get_polygon(vorn, i)
                        for l in 1:length(po)-1
                            a = get_polygon_point(vorn, po[l])
                            v1 = [a[1]-c[1], a[2]-c[2]]
                            v2 = [1,0]
                            ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
                            rθ[i,l,1] = ac
                            rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
                        end
                    end
                end 
            end

            for p=2:6
                aγcomp = zeros(ComplexF64, Npoly)
                aγ = zeros(Npoly)
                aθ = zeros(Npoly)
                for i=1:Npoly
                    if !(i in boundary)
                    c = DelaunayTriangulation.get_centroid(vorn, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

                            Δp = sum(rθ[i,:,1].^p)
                            γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
                            # aγ[i] = γ
                            aγcomp[i] = γ
                            aγ[i] = abs.(γ)
                            aθ[i] = angle.(γ)/p
                        end
                    end
                end
                mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
                mgamma = mean(filter!(x -> x != 0, aγ))
                mtheta = mean(filter!(x -> x != 0, aθ))
                mean_gamma_t[file, p-1] = mgamma
                mean_theta_t[file, p-1] = mtheta
            end
        end
    end
    gam = DataFrame(x=(1:(num_files-1)))
    for i=2:6
        gam[!,"g$i"] = mean_gamma_t[2:num_files, i-1]
    end
    CSV.write(joinpath(joinpath(dir_fig, "triple_$idx_$k/"),"gamma.csv"), gam)
    return fig1, fig2, gam
    # return fig1, fig2
end

function figure_triple_pm(idx, k)
    dir = joinpath(dir_fig, "triple_pm_$(idx)_$k/")
    mkpath(dir)
    ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    lfile_path = dird*lfile_name
    num_files = length(files_names)
       
    # dird5 = dir_fig*path*"Defects_pos/"
    # files_names5 = isdir(dird5) ? readdir(dird5) : []
    # lfile_name5 = last(files_names5)
    # lfile_path5 = dird5*lfile_name5
    
      
    density = load(string(dird,lfile_name), "rho")

    Nx, Nz = size(density)
    x = 1:stepp:Nx
    y = 1:stepp:Nz
    
    P_cpu = load(string(dird,lfile_name), "P")
    Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    
    fig1 = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density, colormap=:RdPu_8)
    Colorbar(fig1[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
    arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1.layout, 4)
    rowgap!(fig1.layout, 0)
    save(string(dir,"density.png"), fig1)
    # fig1
    
    
    # pts = load(lfile_path5, "pos")
    # color = load(lfile_path5, "colors")
    # si = load(lfile_path5, "charges")
    # it = findall(x-> x<-0.2, -si)
    # ptsm = zeros(length(it), 2)
    # for i in 1:length(it)
    #     ptsm[i,:] .= Int.(pts[it[i],:])
    # end
        
    # it = findall(x->x>0.2, -si)
    # ptsp = zeros(length(it), 2)
    # for i in 1:length(it)
    #     ptsp[i,:] .= Int.(pts[it[i],:])
    # end

    # # All defects
    # tri = triangulate(pts')
    # vorn = voronoi(tri, true)
    # map_points = vorn.generators
    # # Minus defects
    # trim = triangulate(ptsm')
    # vornm = voronoi(trim, true)
    # map_pointsm = vornm.generators
    # polysm = vornm.polygons
    # boundarym = vornm.boundary_polygons
    # Npolym = length(polysm)
    # # Positive defects
    # trip = triangulate(ptsp')
    # vornp = voronoi(trip, true)
    # map_pointsp = vornp.generators
    # polysp = vornp.polygons
    # boundaryp = vornp.boundary_polygons
    # Npolyp = length(polysp)
    

    # fig2 = Figure(size = (1000, 1025))
    # ax = Makie.Axis(fig2[2,1], limits=(ex, ex+Nx, ex, ex+Nz), aspect = DataAspect())
    # sc = voronoiplot!(ax, vornm, strokewidth=4, strokecolor=:blue, strokestyle=:dashed)
    # sc = voronoiplot!(ax, vornp, strokewidth=4, strokecolor=:red, strokestyle=:dashed)
    # s1 = scatter!(ax, pts', color=color, markersize = 15)
    # elem_1 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    # elem_2 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    # Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
    # hideydecorations!(ax)
    # hidexdecorations!(ax)
    # colgap!(fig2.layout, 0)
    # rowgap!(fig2.layout, -5)
    # tightlimits!(ax)
    # save(joinpath(dir_fig, "triple_pm_$k/voronoi.png"), fig2)
    # # fig2
    
#     mean_gamma_t = zeros(num_files,5)
#     mean_theta_t = zeros(num_files,5)
        
        
        
#     for file in 1:num_files
#         pts = Array{Float64,2}
#         Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
#         pts = load(string(dird5,@sprintf("%04i",file),".jld"), "pos")
#         color = load(string(dird5,@sprintf("%04i",file),".jld"), "colors")
#         si = load(string(dird5,@sprintf("%04i",file),".jld"), "charges")
          
#         if length(pts) > 50
#             it = findall(x-> x<-0.2, -si)
#             ptsm = zeros(length(it), 2)
#             for i in 1:length(it)
#                 ptsm[i,:] .= Int.(pts[it[i],:])
#             end
            
#             it = findall(x->x>0.2, -si)
#             ptsp = zeros(length(it), 2)
#             for i in 1:length(it)
#                 ptsp[i,:] .= Int.(pts[it[i],:])
#             end  
            
#             # Negative defects
#             trim = triangulate(ptsm')
#             vornm = voronoi(trim, true)
#             map_pointsm = vornm.generators
#             polysm = vornm.polygons
#             boundarym = vornm.boundary_polygons
#             Npolym = length(polysm)
#             # Positive defects
#             trip = triangulate(ptsp')
#             vornp = voronoi(trip, true)
#             map_pointsp = vornp.generators
#             polysp = vornp.polygons
#             boundaryp = vornp.boundary_polygons
#             Npolyp = length(polysp)

#             rθ = zeros(Npolym+Npolyp, 30, 2)
#             for i=1:Npolym
#                 if !(i in boundarym)
#                     c = DelaunayTriangulation.get_centroid(vornm, i)
#                     if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
#                         po = get_polygon(vornm, i)
#                         for l in 1:length(po)-1
#                             a = get_polygon_point(vornm, po[l])
#                             v1 = [a[1]-c[1], a[2]-c[2]]
#                             v2 = [1,0]
#                             ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
#                             rθ[i,l,1] = ac
#                             rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
#                         end
#                     end
#                 end 
#             end
#             for i=1:Npolyp
#                 if !(i in boundaryp)
#                     c = DelaunayTriangulation.get_centroid(vornp, i)
#                     if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
#                         po = get_polygon(vornp, i)
#                         for l in 1:length(po)-1
#                             a = get_polygon_point(vornp, po[l])
#                             v1 = [a[1]-c[1], a[2]-c[2]]
#                             v2 = [1,0]
#                             ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
#                             rθ[i+Npolym,l,1] = ac
#                             rθ[i+Npolym,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
#                         end
#                     end
#                 end 
#             end

#             for p=2:6
#                 aγcomp = zeros(ComplexF64, Npolym+Npolyp)
#                 aγ = zeros(Npolym+Npolyp)
#                 aθ = zeros(Npolym+Npolyp)
#                 for i=1:Npolym
#                     if !(i in boundarym)
#                     c = DelaunayTriangulation.get_centroid(vornm, i)
#                         if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

#                             Δp = sum(rθ[i,:,1].^p)
#                             γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
#                             # aγ[i] = γ
#                             aγcomp[i] = γ
#                             aγ[i] = abs.(γ)
#                             aθ[i] = angle.(γ)/p
#                         end
#                     end
#                 end
#                 for i=1:Npolyp
#                     if !(i in boundaryp)
#                     c = DelaunayTriangulation.get_centroid(vornp, i)
#                         if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

#                             Δp = sum(rθ[Npolym+i,:,1].^p)
#                             γ = 1/Δp*sum(rθ[Npolym+i,:,1].^p .* exp.(im*p*rθ[Npolym+i,:,2]) )
#                             # aγ[i] = γ
#                             aγcomp[Npolym+i] = γ
#                             aγ[Npolym+i] = abs.(γ)
#                             aθ[Npolym+i] = angle.(γ)/p
#                         end
#                     end
#                 end
#                 mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
#                 mgamma = mean(filter!(x -> x != 0, aγ))
#                 mtheta = mean(filter!(x -> x != 0, aθ))
#                 mean_gamma_t[file, p-1] = mgamma
#                 mean_theta_t[file, p-1] = mtheta
#             end
#         end
#     end
#     gam = DataFrame(x=(1:(num_files-1)))
#     for i=2:6
#         gam[!,"g$i"] = mean_gamma_t[2:num_files, i-1]
#     end
#     CSV.write(joinpath(joinpath(dir_fig, "triple_pm_$k/"),"gamma.csv"), gam)
#     return fig1, fig2, gam
    return fig1#, fig2
end

function figure_triple_zoom(idx, k=2)
    mkpath(joinpath(dir_fig, "triple_zoom_$k/"))
    ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    lfile_path = dird*lfile_name
    num_files = length(files_names)
       
    dird5 = dir_fig*path*"Defects_pos/"
    files_names5 = isdir(dird5) ? readdir(dird5) : []
    lfile_name5 = last(files_names5)
    lfile_path5 = dird5*lfile_name5
    
      
    density = load(string(dird,lfile_name), "rho")

    Nx, Nz = size(density)
    ti = floor(Int, Nx/3)
    tN = floor(Int, Nx-Nx/3+1)
    x = ti:10:tN
    y = ti:10:tN
    @show length(y)
    xx = 1+10:10:ti+10
    yy = 1+10:10:ti+10
    P_cpu = load(string(dird,lfile_name), "P")
    # Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    # Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    Px_cpu = P_cpu[x,y,1]
    Pz_cpu = P_cpu[x,y,2]
    
    
    fig1 = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density[ti-10:tN,ti-10:tN])
    Colorbar(fig1[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
    arrows!(ax1, xx, yy, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*15.0, linecolor=:black, lengthscale = 8.0, linewidth=4.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1.layout, 4)
    rowgap!(fig1.layout, 0)
    save(joinpath(dir_fig, "triple_zoom_$k/density.png"), fig1)
    # fig1
    
    
    # pts = load(lfile_path5, "pos")
    # color = load(lfile_path5, "colors")
    # tri = triangulate(pts')
    # vorn = voronoi(tri, true)
    # map_points = vorn.generators
    # fig2 = Figure(size = (1000, 1025))
    # ax = Makie.Axis(fig2[2,1], limits=(ex, ex+Nx, ex, ex+Nz), aspect = DataAspect())
    # sc = voronoiplot!(ax, vorn, strokewidth=4, strokecolor=:black)
    # s1 = scatter!(ax, pts', color=color, markersize = 15)
    # # elem_1 = MarkerElement(color = :blue, marker = :circle, markersize = 25)
    # # elem_2 = MarkerElement(color = :red, marker = :circle, markersize = 25)
    # elem_1 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    # elem_2 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    # Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
    # hideydecorations!(ax)
    # hidexdecorations!(ax)
    # colgap!(fig2.layout, 0)
    # rowgap!(fig2.layout, -5)
    # tightlimits!(ax)
    # save(joinpath(dir_fig, "triple_zoom_$k/voronoi.png"), fig2)
    # # fig2
    
#     mean_gamma_t = zeros(num_files,5)
#     mean_theta_t = zeros(num_files,5)
#     for file in 1:num_files
#         pts = Array{Float64,2}
#         Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
#         pts = load(string(dird5,@sprintf("%04i",file),".jld"), "pos")

#         if length(pts) > 50
#             tri = triangulate(pts')
#             vorn = voronoi(tri, true)
#             map_points = vorn.generators
#             polys = vorn.polygons
#             boundary = vorn.boundary_polygons
#             Npoly = length(polys)

#             rθ = zeros(Npoly, 30, 2)
#             for i=1:Npoly
#                 if !(i in boundary)
#                     c = DelaunayTriangulation.get_centroid(vorn, i)
#                     if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
#                         po = get_polygon(vorn, i)
#                         for l in 1:length(po)-1
#                             a = get_polygon_point(vorn, po[l])
#                             v1 = [a[1]-c[1], a[2]-c[2]]
#                             v2 = [1,0]
#                             ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
#                             rθ[i,l,1] = ac
#                             rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
#                         end
#                     end
#                 end 
#             end

#             for p=2:6
#                 aγcomp = zeros(ComplexF64, Npoly)
#                 aγ = zeros(Npoly)
#                 aθ = zeros(Npoly)
#                 for i=1:Npoly
#                     if !(i in boundary)
#                     c = DelaunayTriangulation.get_centroid(vorn, i)
#                         if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

#                             Δp = sum(rθ[i,:,1].^p)
#                             γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
#                             # aγ[i] = γ
#                             aγcomp[i] = γ
#                             aγ[i] = abs.(γ)
#                             aθ[i] = angle.(γ)/p
#                         end
#                     end
#                 end
#                 mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
#                 mgamma = mean(filter!(x -> x != 0, aγ))
#                 mtheta = mean(filter!(x -> x != 0, aθ))
#                 mean_gamma_t[file, p-1] = mgamma
#                 mean_theta_t[file, p-1] = mtheta
#             end
#         end
#     end
#     gam = DataFrame(x=(1:(num_files-1)))
#     for i=2:6
#         gam[!,"g$i"] = mean_gamma_t[2:num_files, i-1]
#     end
#     CSV.write(joinpath(joinpath(dir_fig, "triple_zoom_$k/"),"gamma.csv"), gam)
    # return fig1, fig2, gam
    return fig1#, fig2
end

function figure_triple_zoom_pm(idx, k)
    mkpath(joinpath(dir_fig, "triple_zoom_pm_$k/"))
    ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    lfile_path = dird*lfile_name
    num_files = length(files_names)
       
    dird5 = dir_fig*path*"Defects_pos/"
    files_names5 = isdir(dird5) ? readdir(dird5) : []
    lfile_name5 = last(files_names5)
    lfile_path5 = dird5*lfile_name5
    
      
    density = load(string(dird,lfile_name), "rho")

    Nx, Nz = size(density)
    ti = floor(Int, Nx/3)
    tN = floor(Int, Nx-Nx/3+1)
    x = ti:10:tN
    y = ti:10:tN
    @show length(y)
    xx = 1+10:10:ti+10
    yy = 1+10:10:ti+10
    P_cpu = load(string(dird,lfile_name), "P")
    # Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    # Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    Px_cpu = P_cpu[x,y,1]
    Pz_cpu = P_cpu[x,y,2]
    
    
    fig1 = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density[ti-10:tN,ti-10:tN])
    Colorbar(fig1[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
    arrows!(ax1, xx, yy, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*15.0, linecolor=:black, lengthscale = 8.0, linewidth=4.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1.layout, 4)
    rowgap!(fig1.layout, 0)
    save(joinpath(dir_fig, "triple_zoom_pm_$k/density.png"), fig1)
    # fig1
    
    
    pts = load(lfile_path5, "pos")
    color = load(lfile_path5, "colors")
    si = load(lfile_path5, "charges")
    it = findall(x-> x<-0.2, -si)
    ptsm = zeros(length(it), 2)
    for i in 1:length(it)
        ptsm[i,:] .= Int.(pts[it[i],:])
    end
        
    it = findall(x->x>0.2, -si)
    ptsp = zeros(length(it), 2)
    for i in 1:length(it)
        ptsp[i,:] .= Int.(pts[it[i],:])
    end

    # All defects
    tri = triangulate(pts')
    vorn = voronoi(tri, true)
    map_points = vorn.generators
    # Minus defects
    trim = triangulate(ptsm')
    vornm = voronoi(trim, true)
    map_pointsm = vornm.generators
    polysm = vornm.polygons
    boundarym = vornm.boundary_polygons
    Npolym = length(polysm)
    # Positive defects
    trip = triangulate(ptsp')
    vornp = voronoi(trip, true)
    map_pointsp = vornp.generators
    polysp = vornp.polygons
    boundaryp = vornp.boundary_polygons
    Npolyp = length(polysp)
    

    fig2 = Figure(size = (1000, 1025))
    ax = Makie.Axis(fig2[2,1], limits=(ex, ex+Nx, ex, ex+Nz), aspect = DataAspect())
    sc = voronoiplot!(ax, vornm, strokewidth=4, strokecolor=(:blue, 0.25) )
    sc = voronoiplot!(ax, vornp, strokewidth=4, strokecolor=(:red, 0.25) )
    s1 = scatter!(ax, pts', color=color, markersize = 15)
    elem_1 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    elem_2 = MarkerElement(color = :white, marker = :circle, markersize = 25)
    Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
    hideydecorations!(ax)
    hidexdecorations!(ax)
    colgap!(fig2.layout, 0)
    rowgap!(fig2.layout, -5)
    tightlimits!(ax)
    save(joinpath(dir_fig, "triple_zoom_pm_$k/voronoi.png"), fig2)
    # fig2
    
    mean_gamma_t = zeros(num_files,5)
    mean_theta_t = zeros(num_files,5)
        
        
        
    for file in 1:num_files
        pts = Array{Float64,2}
        Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
        pts = load(string(dird5,@sprintf("%04i",file),".jld"), "pos")
        color = load(string(dird5,@sprintf("%04i",file),".jld"), "colors")
        si = load(string(dird5,@sprintf("%04i",file),".jld"), "charges")
          
        if length(pts) > 50
            it = findall(x-> x<-0.2, -si)
            ptsm = zeros(length(it), 2)
            for i in 1:length(it)
                ptsm[i,:] .= Int.(pts[it[i],:])
            end
            
            it = findall(x->x>0.2, -si)
            ptsp = zeros(length(it), 2)
            for i in 1:length(it)
                ptsp[i,:] .= Int.(pts[it[i],:])
            end  
            
            # Negative defects
            trim = triangulate(ptsm')
            vornm = voronoi(trim, true)
            map_pointsm = vornm.generators
            polysm = vornm.polygons
            boundarym = vornm.boundary_polygons
            Npolym = length(polysm)
            # Positive defects
            trip = triangulate(ptsp')
            vornp = voronoi(trip, true)
            map_pointsp = vornp.generators
            polysp = vornp.polygons
            boundaryp = vornp.boundary_polygons
            Npolyp = length(polysp)

            rθ = zeros(Npolym+Npolyp, 30, 2)
            for i=1:Npolym
                if !(i in boundarym)
                    c = DelaunayTriangulation.get_centroid(vornm, i)
                    if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
                        po = get_polygon(vornm, i)
                        for l in 1:length(po)-1
                            a = get_polygon_point(vornm, po[l])
                            v1 = [a[1]-c[1], a[2]-c[2]]
                            v2 = [1,0]
                            ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
                            rθ[i,l,1] = ac
                            rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
                        end
                    end
                end 
            end
            for i=1:Npolyp
                if !(i in boundaryp)
                    c = DelaunayTriangulation.get_centroid(vornp, i)
                    if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
                        po = get_polygon(vornp, i)
                        for l in 1:length(po)-1
                            a = get_polygon_point(vornp, po[l])
                            v1 = [a[1]-c[1], a[2]-c[2]]
                            v2 = [1,0]
                            ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
                            rθ[i+Npolym,l,1] = ac
                            rθ[i+Npolym,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
                        end
                    end
                end 
            end

            for p=2:6
                aγcomp = zeros(ComplexF64, Npolym+Npolyp)
                aγ = zeros(Npolym+Npolyp)
                aθ = zeros(Npolym+Npolyp)
                for i=1:Npolym
                    if !(i in boundarym)
                    c = DelaunayTriangulation.get_centroid(vornm, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

                            Δp = sum(rθ[i,:,1].^p)
                            γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
                            # aγ[i] = γ
                            aγcomp[i] = γ
                            aγ[i] = abs.(γ)
                            aθ[i] = angle.(γ)/p
                        end
                    end
                end
                for i=1:Npolyp
                    if !(i in boundaryp)
                    c = DelaunayTriangulation.get_centroid(vornp, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

                            Δp = sum(rθ[Npolym+i,:,1].^p)
                            γ = 1/Δp*sum(rθ[Npolym+i,:,1].^p .* exp.(im*p*rθ[Npolym+i,:,2]) )
                            # aγ[i] = γ
                            aγcomp[Npolym+i] = γ
                            aγ[Npolym+i] = abs.(γ)
                            aθ[Npolym+i] = angle.(γ)/p
                        end
                    end
                end
                mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
                mgamma = mean(filter!(x -> x != 0, aγ))
                mtheta = mean(filter!(x -> x != 0, aθ))
                mean_gamma_t[file, p-1] = mgamma
                mean_theta_t[file, p-1] = mtheta
            end
        end
    end
    gam = DataFrame(x=(1:(num_files-1)))
    for i=2:6
        gam[!,"g$i"] = mean_gamma_t[2:num_files, i-1]
    end
    CSV.write(joinpath(joinpath(dir_fig, "triple_zoom_pm_$k/"),"gamma.csv"), gam)
    return fig1, fig2, gam
    # return fig1, fig2
end
