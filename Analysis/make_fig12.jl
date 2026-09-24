# =============================================================================
#  Figure 12 - fig:lattice1extended
#
#  Square and hexagonal lattices: density, Voronoi tessellation and the
#  time evolution of the shape order.
#
#      DATA_DIR=/path/to/fig12-runs/ FIG_DIR=<repo>/figures/Fig12/ \
#          julia --project=. Analysis/make_fig12.jl
#
#  Inputs : params/DF_12.csv (2 runs, t_fin 2e6)
#  Outputs: density, Voronoi and shape-order panels for both lattices
# =============================================================================

include("MakePlots.jl")

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

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(joinpath(dir_fig, "Data_Fig_Tikz"))
    figure_triple(1)
    figure_triple(2)
end
