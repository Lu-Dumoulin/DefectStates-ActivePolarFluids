# =============================================================================
#  Figures that did not survive the revision.
#
#  The earlier version of the article had figures for flow alignment (nu_1),
#  anisotropic active stress and saturation. The submitted version does not:
#  its thirteen figures are schema, stab_defects, singleDefect, Ndefects,
#  linstab, states, phasediagram, lattice1, gamma_rho, tauc, statesL50,
#  lattice1extended and 10lattices.
#
#  Kept so the analysis behind the dropped panels stays inspectable. Nothing in
#  figures/ comes from here.
# =============================================================================

include("MakePlots.jl")

# ----- flow alignment (nu_1) -----
function make_small_heatmap_nu(tkd=[0.2], tr0=[0.8], tzr=[4], tnu=[-1.1, 1.1]; part = 1.0, Lz = 2, Lx = 1)
    # set_theme!(empty_theme)
    df2 = CSV.read(dir_df*"DF.csv", DataFrame)
    tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& in.(df2.nu1, Ref(tnu)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    pos=zeros(Int,length(tidx))
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
        @show pos[i] = Int(length(tnu)*(findfirst(isequal(df2[idx,:rho0]), tr0)-1) + findfirst(isequal(df2[idx,:nu1]), tnu))
    end
    
    fig = Figure(resolution=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lz for j in 1:Lx]
    
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

        heatmap!(axs[pos[i]], rho[1:pNx,1:pNz], colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        # arrows!(axs[pos[i]], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
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
    save(joinpath(dir_fig, "phases_small_$(tkd[1])_$(tr0[1])_$(tzr[1])_nu_v2.png"), fig)
    return fig
end

function make_full_heatmap_nu(tkd=[0.2], tr0=[0.75,0.8], tzr=[4], tnu=[-1.1, -0.5, 0.5, 1.1]; part = 1.0, Lz = 2, Lx = 4)
    # set_theme!(empty_theme)
    df2 = CSV.read(dir_df*"DF.csv", DataFrame)
    tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& in.(df2.nu1, Ref(tnu)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    pos=zeros(Int,length(tidx))
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
        @show pos[i] = Int(length(tnu)*(findfirst(isequal(df2[idx,:rho0]), tr0)-1) + findfirst(isequal(df2[idx,:nu1]), tnu))
    end
    
    fig = Figure(resolution=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lz for j in 1:Lx]
    
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

        heatmap!(axs[pos[i]], rho[1:pNx,1:pNz], colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        # arrows!(axs[pos[i]], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
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
    save(joinpath(dir_fig, "phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_nu_v2.png"), fig)
    return fig
end

function make_full_heatmap_nu2(tkd=[0.2], tr0=[0.65,0.75,0.8], tzr=[4], tnu=[-1.1, -0.5, 0.5, 1.1]; part = 1.0, Lz = 3, Lx = 4)
    # set_theme!(empty_theme)
    df2 = CSV.read(dir_df*"DF.csv", DataFrame)
    tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& in.(df2.nu1, Ref(tnu)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    pos=zeros(Int,length(tidx))
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
        @show pos[i] = Int(length(tnu)*(findfirst(isequal(df2[idx,:rho0]), tr0)-1) + findfirst(isequal(df2[idx,:nu1]), tnu))
    end
    
    fig = Figure(resolution=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lz for j in 1:Lx]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        zr = df2[idx,:zetarho]
        r0 = df2[idx,:rho0]
        nu = df2[idx,:nu1]
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

        if r0 != 0.65 || nu != -1.1
            heatmap!(axs[pos[i]], rho[1:pNx,1:pNz], colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
            # arrows!(axs[pos[i]], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
        end
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
    save(joinpath(dir_fig, "phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_nu2_v2.png"), fig)
    return fig
end

# ----- anisotropic active stress -----
function make_plot_aniso(tidx=[1598, 1828, 2096, 2171]; part = 1.0, Lz = 2, Lx = 2)
    # set_theme!(empty_theme)
    dird_l = "Z:/zetanu/"
    dir_df_l = "Z:/zetanu/"
    df2 = CSV.read("Z:/zetanu/DF.csv", DataFrame)
    # tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& in.(df2.nu1, Ref(tnu)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    pos=zeros(Int,length(tidx))
    for i = 1:length(tidx)
        idx = tidx[i]
        dird_l = dir_df_l*"$idx/"*"Data/"
        rho = load(string(dird_l,readdir(dird_l)[end]), "rho")
        # mrho = maximum(rho)
        # rhomax = mrho > rhomax ? mrho : rhomax
        # mrho = minimum(rho)
        # rhomin = mrho < rhomin ? mrho : rhomin
        rhomax[i] = maximum(rho)
        rhomin[i] = minimum(rho)
        # @show pos[i] = Int(length(tnu)*(findfirst(isequal(df2[idx,:rho0]), tr0)-1) + findfirst(isequal(df2[idx,:nu1]), tnu))
    end
    
    fig = Figure(resolution=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lz for j in 1:Lx]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        zr = df2[idx,:zetarho]
        path = string("$idx/")
        dird_l = dir_df_l*path*"Data/"
        last_file = readdir(dird_l)[end]
        rho = load(string(dird_l,last_file), "rho")
        Nx, Nz = size(rho)
        pNx = round(Int, Nx*part)
        pNz = round(Int, Nz*part)
        x = 1:steppp:pNx
        y = 1:steppp:pNz
        P_cpu = load(string(dird_l,last_file), "P")
        Px_cpu = P_cpu[x,y,1]
        Pz_cpu = P_cpu[x,y,2]
        arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))

        heatmap!(axs[i], rho[1:pNx,1:pNz], colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        # arrows!(axs[i], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
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
    save(joinpath(dir_fig, "zetanu_v2.png"), fig)
    return fig
end

function make_plot_aniso2(tidx=[9,10,12,13]; part = 1.0, Lz = 2, Lx = 2)
    # set_theme!(empty_theme)
    dird_l = "Z:/Aniso2/"
    dir_df_l = "Z:/Aniso2/"
    df2 = CSV.read("Z:/Aniso2/DF.csv", DataFrame)
    # tidx = df2[in.(df2.kd, Ref(tkd)) .& in.(df2.rho0, Ref(tr0)) .& in.(df2.nu1, Ref(tnu)) .& (in.(df2.zetarho, Ref(tzr))), :fn]
    steppp = 15
    rhomax = zeros(length(tidx))
    rhomin = ones(length(tidx))
    pos=zeros(Int,length(tidx))
    for i = 1:length(tidx)
        idx = tidx[i]
        dird_l = dir_df_l*"$idx/"*"Data/"
        rho = load(string(dird_l,readdir(dird_l)[end]), "rho")
        # mrho = maximum(rho)
        # rhomax = mrho > rhomax ? mrho : rhomax
        # mrho = minimum(rho)
        # rhomin = mrho < rhomin ? mrho : rhomin
        rhomax[i] = maximum(rho)
        rhomin[i] = minimum(rho)
        # @show pos[i] = Int(length(tnu)*(findfirst(isequal(df2[idx,:rho0]), tr0)-1) + findfirst(isequal(df2[idx,:nu1]), tnu))
    end
    
    fig = Figure(resolution=(800*Lx+25,800*Lz))
    axs = [Makie.Axis(fig[i,j]; aspect = DataAspect()) for i in 1:Lz for j in 1:Lx]
    gls = [GridLayout(fig[i,j]; tellwidth = false, tellheight = false, halign = :right, valign = :top) for i in 1:Lz for j in 1:Lx]
    
    for i = 1:length(tidx)
        idx = tidx[i]
        zr = df2[idx,:zetarho]
        path = string("$idx/")
        dird_l = dir_df_l*path*"Data/"
        last_file = readdir(dird_l)[end]
        rho = load(string(dird_l,last_file), "rho")
        Nx, Nz = size(rho)
        pNx = round(Int, Nx*part)
        pNz = round(Int, Nz*part)
        x = 1:steppp:pNx
        y = 1:steppp:pNz
        P_cpu = load(string(dird_l,last_file), "P")
        Px_cpu = P_cpu[x,y,1]
        Pz_cpu = P_cpu[x,y,2]
        arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))

        heatmap!(axs[i], rho[1:pNx,1:pNz], colormap=:RdPu_8, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        # arrows!(axs[i], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
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
    save(joinpath(dir_fig, "aniso2_v2.png"), fig)
    return fig
end

# ----- saturation -----
function plot_saturation(tidx = [9, 29, 119])
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    Nidx = length(tidx)
    colors = Makie.wong_colors()[1:Nidx]
    tt = zeros(Nidx)
    for i=1:Nidx
        idx = tidx[i]
        path = string("$idx/")
        dird5 = dir_fig*path*"Defects_pos/"
        tt[i] = length(readdir(dird5))
    end
    Nt = Int(minimum(tt))
    Nd = zeros(Nidx, Nt)
    for i=1:Nidx
        idx = tidx[i]
        path = string("$idx/")
        dird5 = dir_fig*path*"Defects_pos/"
        for j=1:Nt
            Nd[i,j] = load(string(dird5,@sprintf("%04i",j),".jld"), "N")
        end
    end
    f, a, ls = lines(Nd[1,2:end], linewidth=4, color = colors[1], label=L"\tau=0.2")
    a.xlabel = "Time"
    a.ylabel = "Number of defects"
    a.xlabelsize = 30
    a.ylabelsize = 30
    a.xticklabelsize = 24
    a.yticklabelsize = 24
    lines!(a, Nd[2,2:end], linewidth=4, color = colors[2], label=L"\tau=1")
    lines!(a, Nd[3,2:end], linewidth=4, color = colors[3], label=L"\tau=10")
    # Legend(f[1,1], labelsize = 24, tellheight = false, tellwidth=false, halign = :right, valign = :bottom, margin = (10, 10, 10, 10))#, )
    axislegend(; position = :rc, labelsize = 24)
    save(joinpath(dir_fig, "saturation.png"), f)
    sat = DataFrame(x=(1:(Nt-1))*1000)
    for i=1:Nidx
        sat[!,"y$i"] = Nd[i,2:Nt]
    end
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","saturation.csv"), sat)
end

function plot_saturation2()
    df = CSV.read(dir_df*"DF.csv", DataFrame)
    kd_unique = unique(df[:, :kd])
    Nunique = length(kd_unique)
    tt = zeros(N_sim)
    
    for idx=1:N_sim
        path = string("$idx/")
        dird5 = dir_fig*path*"Defects_pos/"
        tt[idx] = length(readdir(dird5))
        # @show idx, tt[idx]
    end
    Nt = Int(minimum(tt[1:N_sim]))
    @show findmin(tt[1:N_sim])
    # Nd = zeros(N_sim, Nt)
    
    sat = DataFrame(x=1:(Nt-1))
    Nd_kd = zeros(Nunique, round(Int, N_sim/Nunique), Nt)
    for idx=1:N_sim
        path = string("$idx/")
        dird5 = dir_fig*path*"Defects_pos/"
        i = findfirst(isequal(df[idx,:kd]), kd_unique)
        j = mod(idx,1:round(Int, N_sim/Nunique))
        sat[!,"t$i$j"] .= 0.0
        sat[!,"c$i$j"] .= 0.0
        for t=2:Nt
            Nd_kd[i, j, t] = load(string(dird5,@sprintf("%04i",t),".jld"), "N")
            sat[t-1,"t$i$j"] = Nd_kd[i,j,t]
            sat[t-1,"c$i$j"] = Nd_kd[i,j,t]/3
        end
    end
    
    # for i=1:Nunique
    #     sat[!,"m$i"] .= 0.0
    #     sat[!,"ps$i"] .= 0.0
    #     sat[!,"ms$i"] .= 0.0
    #     sat[!,"min$i"] .= 0.0
    #     sat[!,"max$i"] .= 0.0
    #     sat[!,"max$i"] .= 0.0
    #     for t=1:Nt-1
    #         sat[t,"m$i"] = mean(Nd_kd[i,:,t+1])
    #         sat[t,"ps$i"] = mean(Nd_kd[i,:,t+1]) + std(Nd_kd[i,:,t+1])
    #         sat[t,"ms$i"] = mean(Nd_kd[i,:,t+1]) - std(Nd_kd[i,:,t+1])
    #         sat[t,"min$i"] = minimum(Nd_kd[i,:,t-1])
    #         sat[t,"max$i"] = maximum(Nd_kd[i,:,t-1])
    #     end
    # end
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","saturation3.csv"), sat)
end
