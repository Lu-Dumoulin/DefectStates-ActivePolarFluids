# =============================================================================
#  Slides for a talk (soutenance), NOT paper figures. Kept for reference only;
#  nothing in figures/ is produced here.
# =============================================================================

include("MakePlots.jl")

function plot_soutenance(idx; dt = 1)
    set_theme!(theme_black())
    dir_figg="Z:/fig_soutenance_b/"

    path = string("$idx/")
    dirf = dir_figg*path*"Density/"
    dirf2 = dir_figg*path*"Angle/"
    dirf3 = dir_figg*path*"Velocity/"
    dirf4 = dir_figg*path*"Order/"
    
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    mkpath(dirf)
    mkpath(dirf2)
    mkpath(dirf3)
    mkpath(dirf4)
    num_files = length(files_names)
    num_pict = length(filter!(endswith(".png"), readdir(dirf)));
    println("Data in $dird")
    println("Figure will be saved in $dirf")
    println(num_files-num_pict, " / ", num_files) 
    
    if true#num_files > 0 && num_files > num_pict
        # try
        Nx, Nz = size(load(string(dird,files_names[1]), "rho"))
        x = 1:stepp:Nx
        y = 1:stepp:Nz
        hasP = haskey(load(string(dird,files_names[1])), "P")
        hasQ = haskey(load(string(dird,files_names[1])), "Q")
        sttopp = num_files # == 21 ? num_files : num_files-1
        for i=2:4:sttopp
            if true#!isfile(string(dirf,@sprintf("%04i",i),".png")) #|| hasQ
                # try
                t = floor(Int,parse(Int, files_names[i][5:end-4])*dt)
                    density = load(string(dird,files_names[i]), "rho")

                    v_cpu = load(string(dird,files_names[i]), "v")
                    max_v = sqrt(maximum(v_cpu[:,:,1].^2 .+ v_cpu[:,:,2].^2))
                    v_cpu ./= max_v
                    vx_cpu = v_cpu[1:stepp:end,1:stepp:end,1]
                    vz_cpu = v_cpu[1:stepp:end,1:stepp:end,2]


                    P_cpu = hasP ? load(string(dird,files_names[i]), "P") : nothing
                    Px_cpu = hasP ? P_cpu[1:stepp:end,1:stepp:end,1] : nothing
                    Pz_cpu = hasP ? P_cpu[1:stepp:end,1:stepp:end,2] : nothing

                    Q_cpu = hasQ ? load(string(dird,files_names[i]), "Q") : nothing
                    Q1_cpu = hasQ ? Q_cpu[1:stepp:end,1:stepp:end,1] : nothing
                    Q2_cpu = hasQ ? Q_cpu[1:stepp:end,1:stepp:end,2] : nothing
                    S = hasQ ? sqrt.(Q1_cpu.^2 .+ Q2_cpu.^2) : nothing
                    theta = hasQ ?  0.5*atan.(Q2_cpu.+1e-6, Q1_cpu) : (hasP ? atan.(Pz_cpu.+1e-6,Px_cpu) : 0.0*similar(vx_cpu) )
                    Qx = hasQ ? S.*cos.(theta) : nothing
                    Qz = hasQ ? S.*sin.(theta) : nothing


                    fig1 = Figure(size = (1000, 1000))
                        ax1 = Makie.Axis(fig1[1,1], aspect = DataAspect())
                        # ax1.title = string("Density and orientation " * string(t))
                        # h1 = heatmap!(ax1, 1:4:Nx, 1:4:Nz, density[1:4:end, 1:4:end], colormap=:RdPu)
                        h1 = heatmap!(ax1, density)
                        try
                            # Colorbar(fig1[1, 1], h1, ticklabelsize=15, size=20, vertical = false)
                        catch
                        end
                        if hasQ
                            arrow_size_Q = vec(S)
                            arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 10.0, linewidth=2.0)
                            # arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 8.0, linewidth=2.0)
                        end
                        if hasP
                            arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
                            arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
                            # arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P, linecolor=:black, lengthscale = 8.0, linewidth=2.0, arrowhead = '▲')
                        end
                        tightlimits!(ax1)
                        hideydecorations!(ax1)
                        hidexdecorations!(ax1)
                        # colgap!(fig1.layout, 4)
                        # rowgap!(fig1.layout, 0)
                    save(string(dirf,@sprintf("%04i",i),".png"), fig1)

                    # fig2 = Figure(size = (1000, 1000))
                    #     ax1 = Makie.Axis(fig2[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                    #     subs = hasQ ? " of Q " : " of P "
                    #     ax1.title = string("Angle"* subs * string(t))
                    #     theta2 = hasQ ?  0.5*atan.(Q_cpu[:,:,2].+1e-6, Q_cpu[:,:,1]) : (hasP ? atan.(P_cpu[:,:,2].+1e-6,P_cpu[:,:,1]) : 0.0*similar(vx_cpu))
                    #     h1 = heatmap!(ax1, theta2, colormap=:hsv, colorrange = (-π, π))
                    #     try
                    #         Colorbar(fig2[1, 2], h1, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
                    #     catch
                    #     end
                    # save(string(dirf2,@sprintf("%04i",i),".png"), fig2)

                    # fig3 = Figure(size = (1000, 1000))
                    #     ax1 = Makie.Axis(fig3[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                    #     ax1.title = string("Velocity field " * string(t) * ", max v = " * string(max_v))
                    #     h1 = heatmap!(ax1, density)
                    #     arrow_size_v = vec(norm.(Vec2f.(vx_cpu,vz_cpu)))
                    #     arrows!(ax1, x, y, vx_cpu, vz_cpu, arrowsize = arrow_size_v*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
                    #     try
                    #         Colorbar(fig3[1, 2], h1)
                    #     catch
                    #     end
                    # save(string(dirf3,@sprintf("%04i",i),".png"), fig3)
                
                    # fig4 = Figure(size = (1000, 1000))
                    #     ax1 = Makie.Axis(fig4[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                    #     subs = hasQ ? " of Q " : " of P "
                    #     ax1.title = string("Order"* subs * string(t))
                    #     ordre = hasQ ? sqrt.(Q_cpu[:,:,1].^2 .+ Q_cpu[:,:,2].^2) : (hasP ? sqrt.(P_cpu[:,:,1].^2 .+ P_cpu[:,:,2].^2) : 0.0*similar(vx_cpu))
                    #     h1 = heatmap!(ax1, ordre, colorrange = (0, 2))
                    #     try
                    #         Colorbar(fig4[1, 2], h1) #, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
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
    if true#num_files > 0 && ( num_files > num_pict || !isfile(dirf*string(idx)*".gif") )
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
        # length(readdir(dirf3)) > 0 ? pngstogif(dirf3, dirf3, string(idx), 8) : nothing
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

function make_full_heatmap_idx_sout(tkd=[0.2], tr0=[0.4,0.5,0.6,0.65,0.7,0.75,0.8,1.2], tzr=[4]; part = 1.0, Lz = 3, Lx = 3)
    set_theme!(theme_black())
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

        heatmap!(axs[i], rho[1:pNx,1:pNz], colormap=:viridis, colorrange=(rhomax[i]-rhomin[i]<0.01 ? (0, rhomax[i]) : (rhomin[i],rhomax[i])) )
        arrows!(axs[i], x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=i<=0 ? :white : :black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
        # Box(gls[i][1,1], color = (:white, 0.75), strokecolor = :black, strokewidth = 1)
        # Label(gls[i][1,1], L"\zeta_\rho = %$zr", padding = (10, 10, 10, 10), color = :black, font = :bold, fontsize = 30)
    end
    tightlimits!.(axs)
    hideydecorations!.(axs)
    hidexdecorations!.(axs)
    # colgap!(fig.layout, 1)
    # rowgap!(fig.layout, 1)
    Colorbar(fig[:, end+1], colormap=:viridis, colorrange=(rhomin[1],rhomax[1]), ticksvisible=false, size = 25, ticklabelsize=0) #, label= L"\rho", labelsize = 30, ticklabelsize=24)
    isdir("Z:/fig_soutenance_b/") ? nothing : mkpath("Z:/fig_soutenance_b/")
    save("Z:/fig_soutenance_b/phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_v2.png", fig)
    return fig
end

function make_zoom_heatmap_idx_sout(tkd=[0.2], tr0=[0.4,0.7, 0.75,0.8], tzr=[4]; part = 0.25, Lz = 2, Lx = 2, posx=400, posz=750)
    set_theme!(theme_black())
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

        heatmap!(axs[i], rho[posx:posx+pNx,posz:posz+pNz], colormap=:viridis, colorrange=(rhomin[i],rhomax[i]))
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
    isdir("Z:/fig_soutenance_b/") ? nothing : mkpath("Z:/fig_soutenance_b/")
    save("Z:/fig_soutenance_b/phases_zoom_$(tkd[1])_$(tr0[1])_$(tzr[1])_v2.png", fig)
    return fig
end

function figure_triple_zoom_pm_sout(idx, k)
    mkpath("Z:/fig_soutenance_b/triple_zoom_pm_$k/")
    set_theme!(theme_black())

    ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    lfile_path = dird*lfile_name
    num_files = length(files_names)
       
    dird5 = dir_df*path*"Defects_pos/"
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
    save("Z:/fig_soutenance_b/triple_zoom_pm_$k/density.png", fig1)
    # fig1

    x = 1:stepp:Nx
    y = 1:stepp:Nz
    Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    fig1b = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1b[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density)
    Colorbar(fig1b[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
    arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
    # arrows!(ax1, P_cpu[:,:,1], P_cpu[:,:,2], arrowsize = arrow_size_P*15.0, linecolor=:black, lengthscale = 8.0, linewidth=4.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1b.layout, 4)
    rowgap!(fig1b.layout, 0)
    save("Z:/fig_soutenance_b/triple_zoom_pm_$k/densityb.png", fig1b)
    
    
    pts = load(lfile_path5, "pos")
    color = load(lfile_path5, "colors")
    color = [i== :blue ? :cyan : :magenta for i in color]
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
    sc = voronoiplot!(ax, vornm, strokewidth=4, strokecolor=(:cyan, 0.5), color=(:black,0) )
    sc = voronoiplot!(ax, vornp, strokewidth=4, strokecolor=(:magenta, 0.5), color=(:black,0) )
    s1 = scatter!(ax, pts', color=color, markersize = 15)
    elem_1 = MarkerElement(color = :black, marker = :circle, markersize = 25)
    elem_2 = MarkerElement(color = :black, marker = :circle, markersize = 25)
    Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
    hideydecorations!(ax)
    hidexdecorations!(ax)
    colgap!(fig2.layout, 0)
    rowgap!(fig2.layout, -5)
    tightlimits!(ax)
    save("Z:/fig_soutenance_b/triple_zoom_pm_$k/voronoi.png", fig2)
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
    CSV.write(joinpath("Z:/fig_soutenance_b/triple_zoom_pm_$k/","gamma.csv"), gam)
    return fig1, fig2, gam
    # return fig1, fig2
end

function figure_triple_zoom_sout(idx, k=2)
    mkpath("Z:/fig_soutenance_b/triple_zoom_$k/")
    ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    lfile_path = dird*lfile_name
    num_files = length(files_names)
       
    dird5 = dir_df*path*"Defects_pos/"
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
    save("Z:/fig_soutenance_b/triple_zoom_$k/density.png", fig1)
    # fig1

    x = 1:stepp:Nx
    y = 1:stepp:Nz
    Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    fig1b = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1b[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density)
    Colorbar(fig1b[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
    arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
    # arrows!(ax1, P_cpu[:,:,1], P_cpu[:,:,2], arrowsize = arrow_size_P*15.0, linecolor=:black, lengthscale = 8.0, linewidth=4.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1b.layout, 4)
    rowgap!(fig1b.layout, 0)
    save("Z:/fig_soutenance_b/triple_zoom_$k/densityb.png", fig1b)
    
    pts = load(lfile_path5, "pos")
    color = load(lfile_path5, "colors")
    color = [i== :blue ? :cyan : :magenta for i in color]
    tri = triangulate(pts')
    vorn = voronoi(tri, true)
    map_points = vorn.generators
    fig2 = Figure(size = (1000, 1025))
    ax = Makie.Axis(fig2[2,1], limits=(ex, ex+Nx, ex, ex+Nz), aspect = DataAspect())
    sc = voronoiplot!(ax, vorn, strokewidth=4, strokecolor=:white, color=(:black,0.0))
    s1 = scatter!(ax, pts', color=color, markersize = 15)
    # elem_1 = MarkerElement(color = :blue, marker = :circle, markersize = 25)
    # elem_2 = MarkerElement(color = :red, marker = :circle, markersize = 25)
    elem_1 = MarkerElement(color = :black, marker = :circle, markersize = 25)
    elem_2 = MarkerElement(color = :black, marker = :circle, markersize = 25)
    Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
    hideydecorations!(ax)
    hidexdecorations!(ax)
    colgap!(fig2.layout, 0)
    rowgap!(fig2.layout, -5)
    tightlimits!(ax)
    save("Z:/fig_soutenance_b/triple_zoom_$k/voronoi.png", fig2)
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
    CSV.write(joinpath("Z:/fig_soutenance_b/triple_zoom_$k/","gamma.csv"), gam)
    # return fig1, fig2, gam
    # return fig1#, fig2
end

function get_Sk_sout(idx, W, kx, kz, hN)
    set_theme!(theme_black())
    # dir_figg="Z:/fig_soutenance_b/"*dir_fig[4:end-1]*"b/"
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []

    num_files = length(files_names)
       
    hasP = haskey(load(string(dird,files_names[1])), "P")
    hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
    # dirf6a = dir_figg*path*"Skrho/"
    # dirf6b = dir_figg*path*"SkP/"
    # mkpath(dirf6a)
    # mkpath(dirf6b)
    # lto = zeros(num_files)
    # ltr = zeros(num_files)
    for file in num_files:num_files
        P_cpu = hasP ? load(string(dird,files_names[file]), "P") : nothing
        rho = load(string(dird,files_names[file]), "rho")
        # rho = maximum(rho) .- rho
        Px = P_cpu[:,:,1]
        Py = P_cpu[:,:,2]
        Nx, Ny = size(Px)
        sko = zeros(hN, hN); skr = zeros(hN, hN)
        ksko = zeros(hN, hN); kskr = zeros(hN, hN)
        order = sqrt.(Px.^2 .+ Py.^2)
        fo = (W * order)
        fr = (W * rho)
        # for i=3:hN-2, j=3:hN-2
            # n_shell = 0
            # for ii=-2:2, jj=-2:2
            #     # sko[i,j] += sqrt(ii*ii+jj*jj) < 2 ? fo[i+ii,j+jj]*fo[Nx-(i-1)-ii,Ny-(j-1)-jj] : 0
            #     # skr[i,j] += sqrt(ii*ii+jj*jj) < 2 ? fr[i+ii,j+jj]*fr[Nx-(i-1)-ii,Ny-(j-1)-jj] : 0
            #     sko[i,j] += sqrt(ii*ii+jj*jj) < 2 ? abs(fo[i+ii,j+jj])^2 : 0#*fo[Nx-(i-1)-ii,Ny-(j-1)-jj] : 0
            #     skr[i,j] += sqrt(ii*ii+jj*jj) < 2 ? abs(fr[i+ii,j+jj])^2 : 0#*fr[Nx-(i-1)-ii,Ny-(j-1)-jj] : 0
            # end
            # ksko[i,j] = sqrt(kx[i]*kx[i]+kz[j]*kz[j])*sko[i,j]/13
            # kskr[i,j] = sqrt(kx[i]*kx[i]+kz[j]*kz[j])*skr[i,j]/13
        # end
        # sko ./= 13
        # skr ./= 13
        # for i=2:hN-1, j=2:hN-1
        #     sko[i,j] = abs(fo[i,j])^2
        #     skr[i,j] = abs(fr[i,j])^2
        # end
        for i=2:hN-1, j=2:hN-1
            ii = Nx-i+2
            jj = Ny-j+2
            if kx[i] != -kx[ii] || kz[j] != -kz[jj]
                println(i, " ", k[i], " ", k[ii])
                return nothing
            end
            sko[i,j] = abs(fo[i,j]*fo[ii,jj])
            skr[i,j] = abs(fr[i,j]*fr[ii,jj])
        end

        # fig = Figure(size = (1000, 1000))
        # ax1 = Makie.Axis(fig[1,1], xlabel = L"k_x", ylabel = L"k_y", aspect = DataAspect())
        # # ax1.title = string("Structure factor Order")
        # heatmap!(ax1, sko[1:100, 1:100])
        # # hideydecorations!(ax1)
        # # hidexdecorations!(ax1)
        # # tightlimits!(ax1)
        # save(string(dirf6a,@sprintf("order_%04i",file),".png"), fig)

        # fig = Figure(size = (1000, 1000))
        # ax1 = Makie.Axis(fig[1,1], xlabel = L"k_x", ylabel = L"k_y", aspect = DataAspect())
        # # ax1.title = string("Structure factor Density")
        # heatmap!(ax1, skr[1:100, 1:100])
        # # hideydecorations!(ax1)
        # # hidexdecorations!(ax1)
        # # tightlimits!(ax1)
        # save(string(dirf6b,@sprintf("density_%04i",file),".png"), fig)
        @show findmax(sko[:,:])
        @show kxx = findmax(sko[:,:])[2][1]
        @show kzz = findmax(sko[:,:])[2][2]
        @show sqrt(kx[kxx]^2 + kz[kzz]^2)
        # lto[file] = 1/sqrt(kx[kxx]^2 + kz[kzz]^2)
        # ltr[file] = 1/sqrt(kx[kxx]^2 + kz[kzz]^2)
    end
    # fig = Figure(size = (1000, 1000))
    #     ax1 = Makie.Axis(fig[1,1], xlabel = L"t", ylabel = L"length")#, aspect = DataAspect())
    #     # ax1.title = string("l(t) Order")
    #     lines!(ax1, lto, linewidth=5)
    #     # hideydecorations!(ax1)
    #     # hidexdecorations!(ax1)
    #     # tightlimits!(ax1)
    #     save(string(dirf6a,"order_lt",".png"), fig)
    # fig = Figure(size = (1000, 1000))
    #     ax1 = Makie.Axis(fig[1,1], xlabel = L"t", ylabel = L"length")#, aspect = DataAspect())
    #     # ax1.title = string("l(t) Density")
    #     lines!(ax1, ltr, linewidth=5)
    #     # hideydecorations!(ax1)
    #     # hidexdecorations!(ax1)
    #     # tightlimits!(ax1)
    #     save(string(dirf6b,"density_lt",".png"), fig)
    # return lto, ltr
end

function get_spectra_sout(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    rho = load(string(dird,files_names[1]), "rho")
    Nx, Ny = size(rho)
    hN = Int(Nx/2)
    
    W = plan_fft(ones(Float64, Nx, Ny))
    # Wi = inv(W)
    kx = 2*pi*fftfreq(Nx, 1/1e-2)
    kz = 2*pi*fftfreq(Ny, 1/1e-2)
    # tab_l = zeros(N_sim,2,length(files_names))
    for idx = 1:N_sim
        # lP, lrho = get_spectrum(idx, W, kx, kz)
        # lP, lrho = get_Sk_sout(idx, W, kx, kz, hN)
        get_Sk_sout(idx, W, kx, kz, hN)
        # println(idx, " ", lP, " ", lrho)
        # tab_l[idx,1,1:length(lP)] .= lP
        # tab_l[idx,2,1:length(lrho)] .= lrho
    end
    # return tab_l
end

function convert_gif_to_mp4(fn)
    str = string("""-i $fn.gif -f lavfi -i anullsrc -vf "scale='trunc(in_w/2)*2':'trunc(in_h/2)*2',format=yuv420p,fps=8.33" -movflags +faststart -shortest $fn.mp4""")
    # run(`ffmpeg $str`)
    println("ffmpeg $str ")
end
