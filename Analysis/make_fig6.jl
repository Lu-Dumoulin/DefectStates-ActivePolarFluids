# =============================================================================
#  Figure 6 - fig:states
#
#  Asymptotic states at different target densities.
#
#      DATA_DIR=/path/to/fig6-runs/ FIG_DIR=<repo>/figures/Fig6/ \
#          julia --project=. Analysis/make_fig6.jl
#
#  Inputs : the L=10 runs at tau=5, zeta_rho=4
#  Outputs: the density/polarity snapshot grid and the zoomed panels
#
#  No DF_6.csv yet - the article does not say which rho0 values are shown.
# =============================================================================

include("MakePlots.jl")

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

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(joinpath(dir_fig, "Data_Fig_Tikz"))
    make_full_heatmap_idx()
    make_zoom_heatmap_idx()
end
