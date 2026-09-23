# =============================================================================
#  Analysis and plotting for the "nu" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_nu.jl
# =============================================================================

include("MakePlots.jl")

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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_small_$(tkd[1])_$(tr0[1])_$(tzr[1])_nu_v2.png", fig)
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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_nu_v2.png", fig)
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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_nu2_v2.png", fig)
    return fig
end

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    make_small_heatmap_nu()
    make_full_heatmap_nu()
    make_full_heatmap_nu2()
end
