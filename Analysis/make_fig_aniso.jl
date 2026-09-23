# =============================================================================
#  Analysis and plotting for the "aniso" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_aniso.jl
# =============================================================================

include("MakePlots.jl")

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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/zetanu_v2.png", fig)
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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/aniso2_v2.png", fig)
    return fig
end

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    make_plot_aniso()
    make_plot_aniso2()
end
