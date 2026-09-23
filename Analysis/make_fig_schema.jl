# =============================================================================
#  Analysis and plotting for the "schema" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_schema.jl
# =============================================================================

include("MakePlots.jl")

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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_L50_3.png", fig)
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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_L50_3zoom.png", fig)
    return fig
end

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    make_plot_L50()
    make_plot_L50_3zoom()
end
