# =============================================================================
#  Generic plotting from simulation output.
#
#  Renders density / angle / velocity / order heatmaps, defect overlays and
#  Voronoi tessellations for any run set, straight from the .jld snapshots.
#  Nothing here reproduces a specific paper figure — for that see the
#  make_fig_*.jl scripts, which include this file.
#
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/MakePlots.jl
# =============================================================================

include("core.jl")

function Plot_per(idx; dt = 1, sc=1)
    # set_theme!(theme_black())
    # IJulia.clear_output(true)
    # dir_df = "F:/Prot_FFT_2D_per/"*pa
    # dir_fig = "F:/Prot_FFT_2D_per/"*pa
    path = string("$idx/")
    dirf = dir_fig*path*"Density_karsten_ba/"
    # dirfb = dir_fig*path*"Density_only/"
    # dirf2 = dir_fig*path*"Angle/"
    # dirf3 = dir_fig*path*"Velocity_max/"
    # dirf3 = dir_fig*path*"Velocity/"
    # dirf4 = dir_fig*path*"Order/"
    # dirf5 = dir_fig*path*"Defects/"
    
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []
    mkpath(dirf)
    # mkpath(dirfb)
    # mkpath(dirf2)
    # mkpath(dirf3)
    # mkpath(dirf4)
    # mkpath(dirf5)
    num_files = length(files_names)
    num_pict = length(filter!(endswith(".png"), readdir(dirf)));
    println("Data in $dird")
    println("Figure will be saved in $dirf")
    println(num_files-num_pict, " / ", num_files) 
    
    if num_files > 0 && num_files > num_pict
        # try
        Nx, Nz = size(load(string(dird,files_names[1]), "rho"))
        x = 1:stepp:Nx
        y = 1:stepp:Nz
        hasP = haskey(load(string(dird,files_names[1])), "P")
        hasQ = haskey(load(string(dird,files_names[1])), "Q")
        # catch
        # end
        # x, y = meshgrid(1:stepp:Nz,1:stepp:Nx)
        sttopp = num_files # == 21 ? num_files : num_files-1
        # sttopp = (sttopp==0) ? sttopp+1 : sttopp
        rmin = 10
        rmax = 0
        for i=1:sttopp 
            density = load(string(dird,files_names[i]), "rho")
            rmin = minimum([rmin, minimum(density)])
            rmax = maximum([rmax, maximum(density)])
        end
        for i=1:sttopp
            if !isfile(string(dirf,@sprintf("%04i",i),".png")) #|| hasQ
                # try
                t = floor(Int,parse(Int, files_names[i][5:end-4])*dt)
                    density = load(string(dird,files_names[i]), "rho")

                    # v_cpu = load(string(dird,files_names[i]), "v")
                    # max_v = sqrt(maximum(v_cpu[:,:,1].^2 .+ v_cpu[:,:,2].^2))
                    # v_cpu ./= max_v
                    # vx_cpu = v_cpu[1:stepp:end,1:stepp:end,1]
                    # vz_cpu = v_cpu[1:stepp:end,1:stepp:end,2]


                    P_cpu = hasP ? load(string(dird,files_names[i]), "P") : nothing
                    Px_cpu = hasP ? P_cpu[1:stepp:end,1:stepp:end,1] : nothing
                    Pz_cpu = hasP ? P_cpu[1:stepp:end,1:stepp:end,2] : nothing

                    # Q_cpu = hasQ ? load(string(dird,files_names[i]), "Q") : nothing
                    # Q1_cpu = hasQ ? Q_cpu[1:stepp:end,1:stepp:end,1] : nothing
                    # Q2_cpu = hasQ ? Q_cpu[1:stepp:end,1:stepp:end,2] : nothing
                    # S = hasQ ? sqrt.(Q1_cpu.^2 .+ Q2_cpu.^2) : nothing
                    # theta = hasQ ?  0.5*atan.(Q2_cpu.+1e-6, Q1_cpu) : (hasP ? atan.(Pz_cpu.+1e-6,Px_cpu) : 0.0*similar(vx_cpu) )
                    # Qx = hasQ ? S.*cos.(theta) : nothing
                    # Qz = hasQ ? S.*sin.(theta) : nothing


                    # fig1 = Figure(size = (500*sc, 500*sc))
                    #     ax1 = Makie.Axis(fig1[1,1], aspect = DataAspect())
                    # #     # ax1.title = string("Density and orientation " * string(t))
                    #     # h1 = heatmap!(ax1, 1:4:Nx, 1:4:Nz, density[1:4:end, 1:4:end], colormap=:RdPu)
                    #     # h1 = heatmap!(ax1, density, colormap=:RdPu, colorrange=(rmin*0.9,rmax*1.1))
                    #     h1 = heatmap!(ax1, density, colormap=:RdPu_8, colorrange=(rmin,rmax))#
                    #     # h1 = heatmap!(ax1, density, colormap=:RdPu, colorrange=(rmin-0.1*rmin,rmax+0.1*rmax))
                    #     # try
                    #     #     Colorbar(fig1[1, 1], h1, ticklabelsize=15, size=20, vertical = false)
                    #     # catch
                    #     # end
                    #     if hasP
                    #         arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
                    #         arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*6.0, color=:black, lengthscale = 10.0, linewidth=1.0, arrowhead = '▲', align = :center)
                    #         # arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, color=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
                    #     end
                    #     if hasQ
                    #         arrow_size_Q = vec(S)
                    #         arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 10.0, linewidth=1.0)
                    #         # arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 8.0, linewidth=2.0)
                    #     end
                    #     tightlimits!(ax1)
                    #     hideydecorations!(ax1)
                    #     hidexdecorations!(ax1)
                    #     # colgap!(fig1.layout, 4)
                    #     # rowgap!(fig1.layout, 0)
                    # save(string(dirf,@sprintf("%04i",i),".png"), fig1)

                    fig1 = Figure(size = (1000, 1000))
                        ax1 = Makie.Axis(fig1[1,1], aspect = DataAspect())
                        # ax1.title = string("Density and orientation " * string(t))
                        # h1 = heatmap!(ax1, 1:4:Nx, 1:4:Nz, density[1:4:end, 1:4:end], colormap=:RdPu)
                        h1 = heatmap!(ax1, density, colorrange=(rmin,rmax), colormap=:RdPu_8)
                        Colorbar(fig1[1, 2], h1)
                        arrows2d!(ax1, x, y, Px_cpu, Pz_cpu, lengthscale = 10, align = :center) 
                        tightlimits!(ax1)
                        hideydecorations!(ax1)
                        hidexdecorations!(ax1)
                        # colgap!(fig1.layout, 4)
                        # rowgap!(fig1.layout, 0)
                    save(string(dirf,@sprintf("%04i",i),".png"), fig1)

                    # fig2 = Figure(size = (1000*sc, 1000*sc))
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

                    # fig3 = Figure(size = (1000*sc, 1000*sc))
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
                
                    # fig4 = Figure(size = (1000*sc, 1000*sc))
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
    if num_files > 0 && ( num_files > num_pict || !isfile(dirf*string(idx)*".gif") )
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
        # pngstogif(dirfb, dirfb, string(idx), 8)
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

function plot_dens_angle(idx)
    dir = joinpath(dir_fig, "density_angle_$(idx)/")
    mkpath(dir)
    # ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    # lfile_path = dird*lfile_name
    # num_files = length(files_names)
        
    density = load(string(dird,lfile_name), "rho")

    Nx, Nz = size(density)
    x = 1:stepp:Nx
    y = 1:stepp:Nz
    
    P_cpu = load(string(dird,lfile_name), "P")
    Px_cpu = P_cpu[1:stepp:end,1:stepp:end,1]
    Pz_cpu = P_cpu[1:stepp:end,1:stepp:end,2]
    
    fig1 = Figure(size = (1000, 1025))
    ax1 = Makie.Axis(fig1[2,1], aspect = DataAspect())
    h1 = heatmap!(ax1, density[1:div(Nx,2),1:end], colormap=:RdPu_8)
    Colorbar(fig1[1, 1], h1, ticklabelsize=0, size=25, vertical = false)
    arrow_size_P = vec(norm.(Vec2f.(Px_cpu[1:div(length(x),2),1:end],Pz_cpu[1:div(length(x),2),1:end])))
    arrows!(ax1, x[1:div(length(x),2)], y, Px_cpu[1:div(length(x),2),1:end], Pz_cpu[1:div(length(x),2),1:end], arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    ax2 = Makie.Axis(fig1[2,2], aspect = DataAspect())
    theta = atan.(P_cpu[:,:,2].+1e-6,P_cpu[:,:,1])
    h2 = heatmap!(ax2, theta[div(Nx,2):end,1:end], colormap=:hsv, colorrange = (-π, π))
    try
        Colorbar(fig1[1, 2], h2, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]), ticklabelsize=0, size=25, vertical = false)
    catch
    end
    tightlimits!(ax2)
    hideydecorations!(ax2)
    hidexdecorations!(ax2)
    colgap!(fig1.layout, 0)
    rowgap!(fig1.layout, 0)
    save(string(dir,"density_angle.png"), fig1)
end

function plot_dens_and_angle(idx)
    dir = joinpath(dir_fig, "density_and_angle_tri_$(idx)/")
    mkpath(dir)
    # ex = 200
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    lfile_name = last(files_names)
    # lfile_path = dird*lfile_name
    # num_files = length(files_names)
        
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
    tightlimits!(ax1)
    hideydecorations!(ax1)
    hidexdecorations!(ax1)
    colgap!(fig1.layout, 0)
    rowgap!(fig1.layout, 0)
    save(string(dir,"density.png"), fig1)

    fig2 = Figure(size = (1000, 1025))
    ax2 = Makie.Axis(fig2[2,1], aspect = DataAspect())
    theta = atan.(P_cpu[:,:,2].+1e-6,P_cpu[:,:,1])
    h2 = heatmap!(ax2, theta, colormap=:hsv, colorrange = (-π, π))
    Colorbar(fig2[1, 1], h2, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]), ticklabelsize=0, size=25, vertical = false)
    tightlimits!(ax2)
    hideydecorations!(ax2)
    hidexdecorations!(ax2)
    colgap!(fig2.layout, 0)
    rowgap!(fig2.layout, 0)
    save(string(dir,"angle.png"), fig2)
end

function plt_defect(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []

    num_files = length(files_names)
       
    hasP = haskey(load(string(dird,files_names[1])), "P")
    hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
    dirf5 = dir_fig*path*"Defects/"
    dird5 = dir_fig*path*"Defects_pos/"
    mkpath(dirf5)
    mkpath(dird5)
    
    Npict = length(readdir(dirf5))
    
    for i=1:num_files
        if !isfile(string(dirf5,@sprintf("%04i",i),".png")) || !isfile(string(dird5,@sprintf("%04i",i),".jld"))
            t = parse(Int, files_names[i][5:end-4])
            P_cpu = hasP ? load(string(dird,files_names[i]), "P") : nothing
            Nx, Ny = size(P_cpu[:,:,1])
            # Px = P_cpu[:,:,1]
            # Py = P_cpu[:,:,2]
            ex = 200
            tex = round(Int, ex*2)
            Px = expend_periodic_array(P_cpu[:,:,1],ex)
            Py = expend_periodic_array(P_cpu[:,:,2],ex)
            # Nx, Ny = size(Px)
            theta2 = atan.(Py.+1e-6,Px)
            P_cpu = nothing
            # Px = load(string(dird,files_names[i]), "P")[:,:,1]
            # Py = load(string(dird,files_names[i]), "P")[:,:,2]
            defects_pos = count_defect(Px, Py, 10)
            Ndef = length(defects_pos[:,1])
            ux = ones(Ndef)*Nx*2
            uy = ones(Ndef)*Ny*2
            for ii=1:Ndef
                xi = defects_pos[ii,1]
                yi = defects_pos[ii,2]
                if (xi<Nx+ex+1) & (xi>ex) & (yi < Ny+ex+1) & (yi > ex)
                    ux[ii] = defects_pos[ii,1]
                    uy[ii] = defects_pos[ii,2]
                    for j=ii+1:Ndef
                        xj = defects_pos[j,1]; 
                        yj = defects_pos[j,2]; 
                        if (xj<Nx+ex+1) & (xj>ex) & (yj < Ny+ex+1) & (yj > ex)
                            xj = xi-xj < -Nx*0.8 ? xj - Nx : xj
                            xj = xi-xj > Nx*0.8 ? xj + Nx : xj
                            yj = yi-yj < -Ny*0.8 ? yj - Ny : yj
                            yj = yi-yj > Ny*0.8 ? yj + Ny : yj
                            dij = sqrt((xi-xj)^2+(yi-yj)^2)
                            if dij < 5
                                ux[ii] = Nx*2
                                uy[ii] = Ny*2
                            end
                        end
                    end
                end
            end
            filter!(x->x!=Nx*2, ux)
            filter!(x->x!=Ny*2, uy)
            defects_pos_unique = hcat(ux,uy)
            
            if isnothing(defects_pos)
                # fig5 = Figure(size = (1100, 1000))
                # ax1 = Makie.Axis(fig5[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect(), limits = (-5, Nx+tex+5, -5, Ny+tex+5))
                # ax1.title = string("Defect position (number of defect = NAN) " * string(t))
                
                # elem_1 = MarkerElement(color = :blue, marker = :circle, markersize = 15)
                # elem_2 = MarkerElement(color = :red, marker = :circle, markersize = 15)
                # Legend(fig5[1, 2], [elem_1, elem_2], ["-1", "+1"], "Defect charge", framevisible = false)
                # save(string(dirf5,@sprintf("%04i",i),".png"), fig5)
            else
                # defects_charge = cara_defects(defects_pos_unique, theta2)
                defects_charge = cara_defects(defects_pos, theta2)
                color = [ abs(i)>0.2 ? (i > 0 ? :red : :blue) : :green for i in defects_charge ]
                Px = nothing; Py = nothing
                Ndefect = length(defects_pos_unique[:,1]) #sum([Int((defects_pos[i,1]> 50) && (defects_pos[i,1]<=(Nx+50)) && (defects_pos[i,2]>50) && (defects_pos[i,2]<=(Ny+50))) for i in 1:length(defects_pos[:,1])])
                # fig5 = Figure(size = (1100, 1000))
                    # ax1 = Makie.Axis(fig5[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect(), limits = (-5, Nx+5, -5, Ny+5))
                    # ax1 = Makie.Axis(fig5[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect(), limits = (-5, Nx+tex+5, -5, Ny+tex+5))
                    # s1 = scatter!(ax1, defects_pos_unique', color=color)
                    # s1 = scatter!(ax1, defects_pos', color=color)
                    # ax1.title = string("Defect position (number of defect = $Ndefect) " * string(t))
                    # lines!(ax1, [-10,Nx+tex+10], [ex,ex], color = :green)
                    # lines!(ax1, [ex,ex], [-10,Ny+tex+10], color = :green)
                    # lines!(ax1, [-10,Nx+tex+10], [Ny+ex,Ny+ex], color = :green)
                    # lines!(ax1, [Nx+ex,Nx+ex], [-10,Ny+tex+10], color = :green)
                # elem_1 = MarkerElement(color = :blue, marker = :circle, markersize = 15)
                # elem_2 = MarkerElement(color = :red, marker = :circle, markersize = 15)
                # Legend(fig5[1, 2], [elem_1, elem_2], ["-1", "+1"], "Defect charge", framevisible = false)
                # save(string(dirf5,@sprintf("%04i",i),".png"), fig5)
                save(string(dird5,@sprintf("%04i",i),".jld"), "N", Ndefect, "pos", defects_pos, "pos_unique", defects_pos_unique, "charges", defects_charge, "colors", color)
                defects_pos = nothing
            end
        end
    end
    # length(readdir(dirf5)) > 0 && (Npict < num_files || !isfile(dirf5*string(idx)*".gif") ) ? pngstogif(dirf5, dirf5, string(idx), 8) : nothing
    # pngstogif(dirf5, dirf5, string(idx), 8)# : nothing
end

function voronoi_idx(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []

    num_files = length(files_names)
       
    hasP = haskey(load(string(dird,files_names[1])), "P")
    hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
    dirf7 = dir_fig*path*"Gamma/"
    dirf8 = dir_fig*path*"Voronoi/"
    mkpath(dirf7)
    mkpath(dirf8)
    
    mean_gamma_comp_t = zeros(num_files, 5)
    mean_gamma_t = zeros(num_files, 5)
    mean_theta_t = zeros(num_files, 5)
    
    for file in 1:num_files
        P_cpu = hasP ? load(string(dird,files_names[file]), "P") : nothing
        # rho = load(string(dird,files_names[file]), "rho")
        Px = P_cpu[:,:,1]
        Py = P_cpu[:,:,2]
        Nx, Ny = size(Px)
        hN = floor(Int, Nx/2); margin = 100
        pts = count_defect(Px, Py, 10)
        # pts = zeros(111,2)
        # for i=0:10, j=0:10
        #     pts[i*10+j+1,1] = i*100
        #     pts[i*10+j+1,2] = j*100
        # end
        Ndefects = size(pts)[1]
        if Ndefects > 50
            for i=1:Ndefects
                xi = pts[i,1]
                yi = pts[i,2]
                if abs(xi-hN) > hN-margin || abs(yi-hN) > hN-margin
                    xn = abs(xi-hN) > hN-margin ? xi-sign(xi-hN)*Nx : xi        
                    yn = abs(yi-hN) > hN-margin ? yi-sign(yi-hN)*Ny : yi
                    pts = vcat(pts, [xn yn])
                end
            end
            tri = triangulate(pts')
            vorn = voronoi(tri, true)
            # vorn = centroidal_smooth(voronoi(tri, true))
            bpoly = vorn.boundary_polygons
            poly_colors = [ poly in bpoly ? :red : :blue for poly in vorn.polygons]
            # fig, ax, sc = voronoiplot(vorn, strokecolor=poly_colors, markersize=9) #:red, markersize=9)
            fig = Figure(size = (1000, 1000))
            ax = Makie.Axis(fig[1,1], aspect = DataAspect())
            sc = voronoiplot!(ax, vorn, strokecolor=:red, markersize=9)
            triplot!(ax, tri, strokewidth=1, strokecolor=(:black, 0.4))
            limits!(ax, 0,1008,0,1008)
            ax.title = string("Voronoi Tessellations " * string(file-1))
            save(string(dirf8,@sprintf("voronoi_%04i",file),".png"), fig)
            polys = vorn.polygons
            map_points = vorn.generators
            boundary = vorn.boundary_polygons
            Npoly = length(polys)
            rθ = zeros(Npoly, 20, 2)
            for i=1:Npoly
                if !(i in boundary)
                    c = DelaunayTriangulation.get_centroid(vorn, i)
                    po = get_polygon(vorn, i)
                    for l in 1:length(po)-1
                        a = get_polygon_point(vorn, po[l])
                        v1 = [a[1]-c[1], a[2]-c[2]]
                        v2 = [1,0]
                        # v2 = [0,1]
                        ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
                        rθ[i,l,1] = ac
                        # rθ[i,l,2] = mod(acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2]), 2π)
                        rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
                    end
                end 
            end
            
            fig = Figure(size = (1600, 800))
            ax = Makie.Axis(fig[1, 1], xlabel = L"p", ylabel = "Distribution", title = string("Distribution of magnitude of shape function " * string(file-1)) )
            for p=2:6
                aγcomp = zeros(ComplexF64, Npoly)
                aγ = zeros(Npoly)
                aθ = zeros(Npoly)
                for i=1:Npoly
                    if !(i in boundary)
                        Δp = sum(rθ[i,:,1].^p)
                        γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
                        # aγ[i] = γ
                        aγcomp[i] = γ
                        aγ[i] = abs.(γ)
                        aθ[i] = angle.(γ)/p
                    end
                end
                mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
                mgamma = mean(filter!(x -> x != 0, aγ))
                mtheta = mean(filter!(x -> x != 0, aθ))
                hist!(ax, aγ, normalization = :pdf, scale_to=0.6, offset=p, direction=:x, bins = 20)
                ylims!(ax, 0, 1)
                # hist!(ax, aγ, normalization = :pdf, scale_to=0.6, offset=p, direction=:x, bar_labels = :values, label_offset=p, label_formatter=x-> round(x, digits=2), bins = 20)
                # println(p, ", ", mgamma, ", ", mtheta)
                mean_gamma_comp_t[file, p-1] = abs(mgammacomp)
                mean_gamma_t[file, p-1] = mgamma
                mean_theta_t[file, p-1] = mtheta
            end
            save(string(dirf7,@sprintf("gamma_%04i",file),".png"), fig)
        end
    end
    isfile(string(dirf7,"mgt",".png")) ? rm(string(dirf7,"mgt",".png")) : nothing
    isfile(string(dirf7,"mgt_comp",".png")) ? rm(string(dirf7,"mgt_comp",".png")) : nothing
    length(readdir(dirf7)) > 0 ? pngstogif(dirf7, dirf7, string(idx), 8) : nothing
    length(readdir(dirf8)) > 0 ? pngstogif(dirf8, dirf8, string(idx), 4) : nothing
    f, a, l = lines(mean_gamma_t[:,1], label="2")
    for i=2:5
        lines!(a, mean_gamma_t[:,i], label="$(i+1)")
    end
    f[1, 2] = Legend(f, a)
    f2, a, l = lines(mean_gamma_comp_t[:,1], label="2")
    for i=2:5
        lines!(a, mean_gamma_comp_t[:,i], label="$(i+1)")
    end
    f2[1, 2] = Legend(f2, a)
    save(string(dirf7,"mgt",".png"), f)
    save(string(dirf7,"mgt_comp",".png"), f2)
    return mean_gamma_t, mean_theta_t
end

function voronoi2_idx(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []

    num_files = length(files_names)
       
    hasP = haskey(load(string(dird,files_names[1])), "P")
    hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
    dird5 = dir_fig*path*"Defects_pos/"
    
    dirf7 = dir_fig*path*"Gamma/"
    dirf8 = dir_fig*path*"Voronoi/"
    mkpath(dirf7)
    mkpath(dirf8)
    
    mean_gamma_comp_t = zeros(num_files, 5)
    mean_gamma_t = zeros(num_files, 5)
    mean_theta_t = zeros(num_files, 5)
    ex = 200
    
    for file in 1:num_files
        pts = Array{Float64,2}
        Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
        t = parse(Int, files_names[file][5:end-4])
        
        # if !isfile(string(dirf8,@sprintf("%04i",file),".png"))
        if !isfile(string(dird5,@sprintf("%04i",file),".jld"))
            # println("compute")
            P_cpu = hasP ? load(string(dird,files_names[file]), "P") : nothing
            Nx, Ny = size(P_cpu[:,:,1])
            Px = expend_periodic_array(P_cpu[:,:,1],ex)
            Py = expend_periodic_array(P_cpu[:,:,2],ex)

            pts = count_defect(Px, Py, 10)
        else
            # println("load")
            pts = load(string(dird5,@sprintf("%04i",file),".jld"), "pos")
            color = load(string(dird5,@sprintf("%04i",file),".jld"), "colors")
        end
        if length(pts) > 50
            tri = triangulate(pts')
            vorn = voronoi(tri, true)
            # smooth_vorn = centroidal_smooth(vorn; maxiters=2500)
            map_points = vorn.generators

            fig = Figure(size = (1000, 1000))
                ax = Makie.Axis(fig[1,1], aspect = DataAspect())
                # sc = voronoiplot!(ax, vorn, strokecolor=:red, markersize=9)
                sc = voronoiplot!(ax, vorn, strokecolor=:black, color=(:white,0))
                s1 = scatter!(ax, pts', color=color)
                # triplot!(ax, tri, strokewidth=1, strokecolor=(:black, 0.4))
                # scatter!(map_pointsap_points[i][1] for i=1:length(pts[:,1])], [map_points[i][2] for i=1:length(pts[:,1])], markersize=9, color=:blue)
                # scatter!(ax, [DelaunayTriangulation.get_centroid(vorn, i) for i in 1:length(vorn.polygons)], markersize=9, color=:green)
                limits!(ax, 0,Nx+2*ex+5,0,Ny+2*ex+5)
                lines!(ax, [-10,Nx+2*ex+10], [ex,ex], color = :green, linewidth=5)
                lines!(ax, [ex,ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                lines!(ax, [-10,Nx+2*ex+10], [Ny+ex,Ny+ex], color = :green, linewidth=5)
                lines!(ax, [Nx+ex,Nx+ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                ax.title = string("Voronoi Tessellations " * string(t))
                save(string(dirf8,@sprintf("voronoi_%04i",file),".png"), fig)

            polys = vorn.polygons
            # map_points = vorn.generators
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
                            # v2 = [0,1]
                            ac = sqrt((a[1]-c[1])^2+(a[2]-c[2])^2)
                            rθ[i,l,1] = ac
                            # rθ[i,l,2] = mod(acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2]), 2π)
                            rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
                        end
                    end
                end 
            end

            fig = Figure(size = (1600, 800))
            ax = Makie.Axis(fig[1, 1], xlabel = L"p", ylabel = "Distribution", title = string("Distribution of magnitude of shape function " * string(t)) )
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
                hist!(ax, aγ, normalization = :pdf, scale_to=0.6, offset=p, direction=:x, bins = 20)
                ylims!(ax, 0, 1)
                ## println(p, ", ", mgamma, ", ", mtheta)
                mean_gamma_comp_t[file, p-1] = abs(mgammacomp)
                mean_gamma_t[file, p-1] = mgamma
                mean_theta_t[file, p-1] = mtheta
            end
            save(string(dirf7,@sprintf("gamma_%04i",file),".png"), fig)
        end
        # end
    end
    isfile(string(dirf7,"mgt",".png")) ? rm(string(dirf7,"mgt",".png")) : nothing
    isfile(string(dirf7,"mgt_comp",".png")) ? rm(string(dirf7,"mgt_comp",".png")) : nothing
    length(readdir(dirf7)) > 0 ? pngstogif(dirf7, dirf7, string(idx), 8) : nothing
    length(readdir(dirf8)) > 0 ? pngstogif(dirf8, dirf8, string(idx), 8) : nothing
    wcolor= Makie.wong_colors()
    f= Figure(size = (1000, 800))
    a = Makie.Axis(f[1, 1], xlabel = "Time", ylabel = L"\gamma_p")
    l = lines!(a, mean_gamma_t[:,1], label=L"p=2", linewidth = 4)
    for i=2:5
        lines!(a, mean_gamma_t[:,i], label=L"p=%$(i+1)", linewidth = 4)
    end
    f[1, 2] = Legend(f, a)
    scatter!(a, num_files, mean_gamma_t[end,3], color = wcolor[3], makersize=8)
    text!(a, num_files*0.85, mean_gamma_t[end,3]*1.02+0.02, text = L"\gamma_4 = %$(round(mean_gamma_t[end,3], digits=2))", color = wcolor[3], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_t[end,2], color = wcolor[2], makersize=8)
    text!(a, num_files*0.85, mean_gamma_t[end,2]*1.02+0.02, text = L"\gamma_3 = %$(round(mean_gamma_t[end,2], digits=2))", color = wcolor[2], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_t[end,5], color = wcolor[5], makersize=8)
    text!(a, num_files*0.85, mean_gamma_t[end,5]*1.02+0.02, text = L"\gamma_6 = %$(round(mean_gamma_t[end,5], digits=2))", color = wcolor[5], font = :bold, fontsize = 24)
    limits!(a, 0, num_files+1.05, 0, 1.0)
    
    f2, a, l = lines(mean_gamma_comp_t[:,1], label=L"p=2")
    for i=2:5
        lines!(a, mean_gamma_comp_t[:,i], label=L"p=%$(i+1)", linewidth = 4)
    end
    f2[1, 2] = Legend(f2, a)
    save(string(dirf7,"mgt",".png"), f)
    # save(string(dirf7,"mgt_comp",".png"), f2)
    return mean_gamma_t, mean_theta_t
end

function voronoi3_idx(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []

    num_files = length(files_names)
       
    hasP = haskey(load(string(dird,files_names[1])), "P")
    hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
    dird5 = dir_fig*path*"Defects_pos/"
    
    dirf7 = dir_fig*path*"Gamma_2/"
    dirf8 = dir_fig*path*"Voronoi_2/"
    mkpath(dirf7)
    mkpath(dirf8)
    
    dirf7m = dir_fig*path*"Gamma_2m/"
    dirf8m = dir_fig*path*"Voronoi_2m/"
    mkpath(dirf7m)
    mkpath(dirf8m)
    
    dirf7p = dir_fig*path*"Gamma_2p/"
    dirf8p = dir_fig*path*"Voronoi_2p/"
    mkpath(dirf7p)
    mkpath(dirf8p)
    
    mean_gamma_comp_t = zeros(num_files, 5)
    mean_gamma_t = zeros(num_files, 5)
    mean_theta_t = zeros(num_files, 5)
    mean_gamma_tm = zeros(num_files, 5)
    mean_theta_tm = zeros(num_files, 5)
    mean_gamma_tp = zeros(num_files, 5)
    mean_theta_tp = zeros(num_files, 5)
    ex = 200
    
    for file in 1:num_files
        pts = Array{Float64,2}
        ptsp = Array{Float64,2}
        Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
        t = parse(Int, files_names[file][5:end-4])
        
        if !isfile(string(dird5,@sprintf("%04i",file),".jld"))
            println("Call plt_defect(idx) first")
            return 0
        end
        pts = load(string(dird5,@sprintf("%04i",file),".jld"), "pos")
        si = load(string(dird5,@sprintf("%04i",file),".jld"), "charges")
        color = load(string(dird5,@sprintf("%04i",file),".jld"), "colors")
        it = findall(x-> x<-0.2, si)
        ptsm = zeros(length(it), 2)
        for i in 1:length(it)
            ptsm[i,:] .= Int.(pts[it[i],:])
        end
        # ptsm .= pts[it,:]
        
        it = findall(x->x>0.2, si)
        ptsp = zeros(length(it), 2)
        for i in 1:length(it)
            ptsp[i,:] .= Int.(pts[it[i],:])
        end
        # ptsp .= pts[it,:]
    
        
        if length(pts) > 50
            tri = triangulate(pts')
            vorn = voronoi(tri, true)
            map_points = vorn.generators
            polys = vorn.polygons
            boundary = vorn.boundary_polygons
            Npoly = length(polys)
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
            
            
            fig = Figure(size = (1000, 1000))
                ax = Makie.Axis(fig[1,1], aspect = DataAspect())
                scp = voronoiplot!(ax, vornp, strokecolor=:blue, )
                scm = voronoiplot!(ax, vornm, strokecolor=:red, )
                s1 = scatter!(ax, pts', color=color)
                limits!(ax, 0,Nx+2*ex+5,0,Ny+2*ex+5)
                lines!(ax, [-10,Nx+2*ex+10], [ex,ex], color = :green, linewidth=5)
                lines!(ax, [ex,ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                lines!(ax, [-10,Nx+2*ex+10], [Ny+ex,Ny+ex], color = :green, linewidth=5)
                lines!(ax, [Nx+ex,Nx+ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                ax.title = string("Voronoi Tessellations " * string(t))
                save(string(dirf8,@sprintf("voronoi_%04i",file),".png"), fig)

            
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

            fig = Figure(size = (1600, 800))
            ax = Makie.Axis(fig[1, 1], xlabel = L"p", ylabel = "Distribution", title = string("Distribution of magnitude of shape function " * string(t)) )
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
                            aγcomp[i] = γ
                            aγ[i] = abs.(γ)
                            aθ[i] = angle.(γ)/p
                        end
                    end
                end
                mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
                mgamma = mean(filter!(x -> x != 0, aγ))
                mtheta = mean(filter!(x -> x != 0, aθ))
                hist!(ax, aγ, normalization = :pdf, scale_to=0.6, offset=p, direction=:x, bins = 20)
                ylims!(ax, 0, 1)
                mean_gamma_comp_t[file, p-1] = abs(mgammacomp)
                mean_gamma_t[file, p-1] = mgamma
                mean_theta_t[file, p-1] = mtheta
            end
            save(string(dirf7,@sprintf("gamma_%04i",file),".png"), fig)
            
           
            fig = Figure(size = (1000, 1000))
                ax = Makie.Axis(fig[1,1], aspect = DataAspect())
                sc = voronoiplot!(ax, vornm, strokecolor=:red, markersize=9)
                limits!(ax, 0,Nx+2*ex+5,0,Ny+2*ex+5)
                lines!(ax, [-10,Nx+2*ex+10], [ex,ex], color = :green, linewidth=5)
                lines!(ax, [ex,ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                lines!(ax, [-10,Nx+2*ex+10], [Ny+ex,Ny+ex], color = :green, linewidth=5)
                lines!(ax, [Nx+ex,Nx+ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                ax.title = string("Voronoi Tessellations " * string(t))
                save(string(dirf8m,@sprintf("voronoi_%04i",file),".png"), fig)



            rθ = zeros(Npolym, 30, 2)
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

            fig = Figure(size = (1600, 800))
            ax = Makie.Axis(fig[1, 1], xlabel = L"p", ylabel = "Distribution", title = string("Distribution of magnitude of shape function " * string(t)) )
            for p=2:6
                aγcomp = zeros(ComplexF64, Npoly)
                aγ = zeros(Npoly)
                aθ = zeros(Npoly)
                for i=1:Npolym
                    if !(i in boundarym)
                    c = DelaunayTriangulation.get_centroid(vornm, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

                            Δp = sum(rθ[i,:,1].^p)
                            γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
                            aγcomp[i] = γ
                            aγ[i] = abs.(γ)
                            aθ[i] = angle.(γ)/p
                        end
                    end
                end
                mgammacomp = mean(filter!(x -> real(x) != 0, aγcomp))
                mgamma = mean(filter!(x -> x != 0, aγ))
                mtheta = mean(filter!(x -> x != 0, aθ))
                hist!(ax, aγ, normalization = :pdf, scale_to=0.6, offset=p, direction=:x, bins = 20)
                ylims!(ax, 0, 1)
                mean_gamma_comp_t[file, p-1] = abs(mgammacomp)
                mean_gamma_tm[file, p-1] = mgamma
                mean_theta_tm[file, p-1] = mtheta
            end
            save(string(dirf7m,@sprintf("gamma_%04i",file),".png"), fig)
            
           

            fig = Figure(size = (1000, 1000))
                ax = Makie.Axis(fig[1,1], aspect = DataAspect())
                sc = voronoiplot!(ax, vornp, strokecolor=:red, markersize=9)
                limits!(ax, 0,Nx+2*ex+5,0,Ny+2*ex+5)
                lines!(ax, [-10,Nx+2*ex+10], [ex,ex], color = :green, linewidth=5)
                lines!(ax, [ex,ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                lines!(ax, [-10,Nx+2*ex+10], [Ny+ex,Ny+ex], color = :green, linewidth=5)
                lines!(ax, [Nx+ex,Nx+ex], [-10,Ny+2*ex+10], color = :green, linewidth=5)
                ax.title = string("Voronoi Tessellations " * string(t))
                save(string(dirf8p,@sprintf("voronoi_%04i",file),".png"), fig)



            rθ = zeros(Npolyp, 30, 2)
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
                            rθ[i,l,1] = ac
                            rθ[i,l,2] = acos(dot(v2/norm(v2),v1/norm(v1)))*sign(a[2]-c[2])
                        end
                    end
                end 
            end

            fig = Figure(size = (1600, 800))
            ax = Makie.Axis(fig[1, 1], xlabel = L"p", ylabel = "Distribution", title = string("Distribution of magnitude of shape function " * string(t)) )
            for p=2:6
                aγcomp = zeros(ComplexF64, Npoly)
                aγ = zeros(Npoly)
                aθ = zeros(Npoly)
                for i=1:Npolyp
                    if !(i in boundaryp)
                    c = DelaunayTriangulation.get_centroid(vornp, i)
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
                hist!(ax, aγ, normalization = :pdf, scale_to=0.6, offset=p, direction=:x, bins = 20)
                ylims!(ax, 0, 1)
                mean_gamma_comp_t[file, p-1] = abs(mgammacomp)
                mean_gamma_tp[file, p-1] = mgamma
                mean_theta_tp[file, p-1] = mtheta
            end
            save(string(dirf7p,@sprintf("gamma_%04i",file),".png"), fig)
        end
        # end
    end
    isfile(string(dirf7,"mgt",".png")) ? rm(string(dirf7,"mgt",".png")) : nothing
    isfile(string(dirf7,"mgt_comp",".png")) ? rm(string(dirf7,"mgt_comp",".png")) : nothing
    length(readdir(dirf7)) > 0 ? pngstogif(dirf7, dirf7, string(idx), 8) : nothing
    length(readdir(dirf8)) > 0 ? pngstogif(dirf8, dirf8, string(idx), 8) : nothing
    wcolor= Makie.wong_colors()
    f= Figure(size = (1000, 800))
    a = Makie.Axis(f[1, 1], xlabel = "Time", ylabel = L"\gamma_p")
    l = lines!(a, mean_gamma_t[:,1], label=L"p=2", linewidth = 4)
    for i=2:5
        lines!(a, mean_gamma_t[:,i], label=L"p=%$(i+1)", linewidth = 4)
    end
    f[1, 2] = Legend(f, a)
    scatter!(a, num_files, mean_gamma_t[end,3], color = wcolor[3], makersize=8)
    text!(a, num_files*0.85, mean_gamma_t[end,3]*1.02+0.02, text = L"\gamma_4 = %$(round(mean_gamma_t[end,3], digits=2))", color = wcolor[3], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_t[end,2], color = wcolor[2], makersize=8)
    text!(a, num_files*0.85, mean_gamma_t[end,2]*1.02+0.02, text = L"\gamma_3 = %$(round(mean_gamma_t[end,2], digits=2))", color = wcolor[2], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_t[end,5], color = wcolor[5], makersize=8)
    text!(a, num_files*0.85, mean_gamma_t[end,5]*1.02+0.02, text = L"\gamma_6 = %$(round(mean_gamma_t[end,5], digits=2))", color = wcolor[5], font = :bold, fontsize = 24)
    limits!(a, 0, num_files+1.05, 0, 1.0)
    
    f2, a, l = lines(mean_gamma_comp_t[:,1], label=L"p=2")
    for i=2:5
        lines!(a, mean_gamma_comp_t[:,i], label=L"p=%$(i+1)", linewidth = 4)
    end
    f2[1, 2] = Legend(f2, a)
    save(string(dirf7,"mgt",".png"), f)
    
    # Neg defects
    f= Figure(size = (1000, 800))
    a = Makie.Axis(f[1, 1], xlabel = "Time", ylabel = L"\gamma_p")
    l = lines!(a, mean_gamma_tm[:,1], label=L"p=2", linewidth = 4)
    for i=2:5
        lines!(a, mean_gamma_tm[:,i], label=L"p=%$(i+1)", linewidth = 4)
    end
    f[1, 2] = Legend(f, a)
    scatter!(a, num_files, mean_gamma_tm[end,3], color = wcolor[3], makersize=8)
    text!(a, num_files*0.85, mean_gamma_tm[end,3]*1.02+0.02, text = L"\gamma_4 = %$(round(mean_gamma_tm[end,3], digits=2))", color = wcolor[3], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_tm[end,2], color = wcolor[2], makersize=8)
    text!(a, num_files*0.85, mean_gamma_tm[end,2]*1.02+0.02, text = L"\gamma_3 = %$(round(mean_gamma_tm[end,2], digits=2))", color = wcolor[2], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_tm[end,5], color = wcolor[5], makersize=8)
    text!(a, num_files*0.85, mean_gamma_tm[end,5]*1.02+0.02, text = L"\gamma_6 = %$(round(mean_gamma_tm[end,5], digits=2))", color = wcolor[5], font = :bold, fontsize = 24)
    limits!(a, 0, num_files+1.05, 0, 1.0)
    
    save(string(dirf7m,"mgt",".png"), f)

    # Pos defects
    f= Figure(size = (1000, 800))
    a = Makie.Axis(f[1, 1], xlabel = "Time", ylabel = L"\gamma_p")
    l = lines!(a, mean_gamma_tp[:,1], label=L"p=2", linewidth = 4)
    for i=2:5
        lines!(a, mean_gamma_tp[:,i], label=L"p=%$(i+1)", linewidth = 4)
    end
    f[1, 2] = Legend(f, a)
    scatter!(a, num_files, mean_gamma_tp[end,3], color = wcolor[3], makersize=8)
    text!(a, num_files*0.85, mean_gamma_tp[end,3]*1.02+0.02, text = L"\gamma_4 = %$(round(mean_gamma_tp[end,3], digits=2))", color = wcolor[3], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_tp[end,2], color = wcolor[2], makersize=8)
    text!(a, num_files*0.85, mean_gamma_tp[end,2]*1.02+0.02, text = L"\gamma_3 = %$(round(mean_gamma_tp[end,2], digits=2))", color = wcolor[2], font = :bold, fontsize = 24)
    scatter!(a, num_files, mean_gamma_tp[end,5], color = wcolor[5], makersize=8)
    text!(a, num_files*0.85, mean_gamma_tp[end,5]*1.02+0.02, text = L"\gamma_6 = %$(round(mean_gamma_tp[end,5], digits=2))", color = wcolor[5], font = :bold, fontsize = 24)
    limits!(a, 0, num_files+1.05, 0, 1.0)
    
    save(string(dirf7p,"mgt",".png"), f)
    
    return mean_gamma_t, mean_theta_t
end

function plot_NandR(N_sim=N_sim)
    fig_N_act = Figure(size = (1000, 1000))
    ax_N_act = Makie.Axis(fig_N_act[1,1], xlabel = L"\frac{k_d}{a_\rho}", ylabel = L"\rho_0")#, aspect = DataAspect())
    ax_N_act.title = latexstring("Nomber of defects for \$\\zeta_\\rho = a_\\rho\$")
    fig_N_pas = Figure(size = (1000, 1000))
    ax_N_pas = Makie.Axis(fig_N_pas[1,1], xlabel = L"\frac{k_d}{a_\rho}", ylabel = L"\rho_0")#, aspect = DataAspect())
    ax_N_pas.title = latexstring("Nomber of defects for \$\\zeta_\\rho = 0\$")
    fig_R_act = Figure(size = (1000, 1000))
    ax_R_act = Makie.Axis(fig_R_act[1,1], xlabel = L"\frac{k_d}{a_\rho}", ylabel = L"\rho_0")#, aspect = DataAspect())
    ax_R_act.title = latexstring("Radius of defects for \$\\zeta_\\rho = a_\\rho\$")
    fig_R_pas = Figure(size = (1000, 1000))
    ax_R_pas = Makie.Axis(fig_R_pas[1,1], xlabel = L"\frac{k_d}{a_\rho}", ylabel = L"\rho_0")#, aspect = DataAspect())
    ax_R_pas.title = latexstring("Radius of defects for \$\\zeta_\\rho = 0\$")
    fig_E_act = Figure(size = (1000, 1000))
    ax_E_act = Makie.Axis(fig_E_act[1,1], xlabel = L"\frac{k_d}{\rho_0^2} \quad (\log_{10})", ylabel = L"a_\rho\left(\rho_0-\frac{3}{4}\right)")#, aspect = DataAspect())
    ax_E_act.title = latexstring("\$\\varepsilon\$ with \$\\zeta_\\rho\\neq0\$")#string(L"\frac{\epsilon}{\rho_0} \text{ for } \zeta_\rho = a_\rho")
    fig_E_pas = Figure(size = (1000, 1000))
    ax_E_pas = Makie.Axis(fig_E_pas[1,1], xlabel = L"\frac{k_d}{\rho_0^2} \quad (\log_{10})", ylabel = L"a_\rho\rho_0")#, aspect = DataAspect())
    ax_E_pas.title = latexstring("\$\\varepsilon\$ with \$\\zeta_\\rho=0\$")
    pointx_NR_act = [];   pointy_NR_act = []
    pointx_KAcr_act = []; pointx_KAcr_pas = [];
    pointy_KAcr_act = [];
    color_N_act = [];    color_N_pas = []
    pointx_NR_pas = [];   pointy_NR_pas = []
    color_R_act = [];    color_R_pas = []
    pointx_E_act = [];   pointy_E_act = []
    pointx_E_pas = [];   pointy_E_pas = []
    color_E_act = [];    color_E_pas = []
    color_E_act_th = [];    color_E_pas_th = []
    include(string(dir_df,"Ndefects.jl"))
    mkpath(string(dir_df,"/1/Others/"))
    for idx=1:N_sim
        zp = df[idx, :zetap]
        g = df[idx, :gamma]
        nu = df[idx, :nu1]
        if zp == 0 && g == 1 && nu == 0 
            kd = df[idx, :kd]/1e-3
            ar = df[idx, :ar]
            a = df[idx, :a]
            ρ0 = df[idx, :rho0]
            if ρ0 > 0.3 && !( (a > 0) & (kd == 0) ) && !( (ar > 2) & (kd == 0) )
                fract = a>0 ? kd/ar : -kd/ar
                ridx = Ndefects[idx] < 3 ? 0 : Radius[idx]
                if fract == 0
                    push!(pointx_NR_pas, fract)
                    push!(pointy_NR_pas, ρ0)
                    push!(color_N_pas, Ndefects[idx])
                    # Ndefects[idx] < 3 ? push!(color2a, 0) : push!(color2a, Radius[idx])
                    Ndefects[idx] < 3 ? push!(color_R_pas, 0) : push!(color_R_pas, log10(Radius[idx]))
                    
                    push!(pointx_E_pas, kd)
                    push!(pointy_E_pas, ar*ρ0)
                    push!(color_E_pas, Epsi[idx])
                    push!(color_E_pas_th, 0.1/(4*ar*ρ0))
                    
                    push!(pointx_KAcr_pas, log10(kdoar(ρ0,0,ridx)))
                    
                elseif fract > 0
                    push!(pointx_NR_act, log10(abs(fract)) )
                    push!(pointy_NR_act, ρ0)
                    push!(color_N_act, Ndefects[idx])
                    # Ndefects[idx] < 3 ? push!(color2a, 0) : push!(color2a, Radius[idx])
                    Ndefects[idx] < 3 ? push!(color_R_act, 0) : push!(color_R_act, log10(Radius[idx]))                    
                    
                    if (ρ0 > 0.75) && (ar < 20) #|| kd > 4                    
                        push!(pointx_E_act, log10(kd/ρ0^2))
                        push!(pointy_E_act, ar*ρ0-ar*0.75)
                        push!(color_E_act, Epsi[idx])
                        push!(color_E_act_th, 0.1/(4*ar*(ρ0-0.75)+2*kd*ridx^2/3/ρ0^2))
                        push!(pointx_KAcr_act, log10(kdoar(ρ0,1,ridx)))
                        push!(pointy_KAcr_act, ρ0)
                    end
                    
                else
                    push!(pointx_NR_pas, log10(abs(fract)) )
                    push!(pointy_NR_pas, ρ0)
                    push!(color_N_pas, Ndefects[idx])
                    # Ndefects[idx] < 3 ? push!(color2b, 0) : push!(color2b, Radius[idx])
                    Ndefects[idx] < 3 ? push!(color_R_pas, 0) : push!(color_R_pas, log10(Radius[idx]))
                    
                    push!(pointx_E_pas, log10(kd/ρ0^2))
                    push!(pointy_E_pas, ar*ρ0)
                    push!(color_E_pas, Epsi[idx])
                    push!(color_E_pas_th, 0.1/(4*ar*ρ0+2*kd*ridx^2/3/ρ0^2))
                    
                    push!(pointx_KAcr_pas, log10(kdoar(ρ0,0,ridx)))
                end
                # push!(pointy, ρ0)
                # push!(color1, Ndefects[idx])
                # Ndefects[idx] < 3 ? push!(color2, 0) : push!(color2, log10(Radius[idx]))
            end
        end
    end

    s_N_act = scatter!(ax_N_act, Vector{Float32}(pointx_NR_act), Vector{Float32}(pointy_NR_act), color = Vector{Float32}(color_N_act), colormap = :thermal, colorrange=(minimum(color_N_act), maximum(color_N_act)) )
    l_N_act = lines!(ax_N_act, Vector{Float32}(pointx_KAcr_act), Vector{Float32}(pointy_KAcr_act))
    s_R_act = scatter!(ax_R_act, Vector{Float32}(pointx_NR_act), Vector{Float32}(pointy_NR_act), color = Vector{Float32}(color_R_act), colormap = :thermal)
    Colorbar(fig_N_act[1, 2], s_N_act, label="Number of defect")
    Colorbar(fig_R_act[1, 2], s_R_act, label="Radius of defect")
    save(string(dir_df,"/1/Others/","ndefects_active.png"), fig_N_act)
    save(string(dir_df,"/1/Others/","radius_active.png"), fig_R_act)
    
    s_N_pas = scatter!(ax_N_pas, Vector{Float32}(pointx_NR_pas), Vector{Float32}(pointy_NR_pas), color = Vector{Float32}(color_N_pas), colormap = :thermal, colorrange=(minimum(color_N_pas), maximum(color_N_pas)) )
    l_N_pas = lines!(ax_N_pas, Vector{Float32}(pointx_KAcr_pas), Vector{Float32}(pointy_NR_pas))
    s_R_pas = scatter!(ax_R_pas, Vector{Float32}(pointx_NR_pas), Vector{Float32}(pointy_NR_pas), color = Vector{Float32}(color_R_pas), colormap = :thermal)
    Colorbar(fig_N_pas[1, 2], s_N_pas, label="Number of defect")
    Colorbar(fig_R_pas[1, 2], s_R_pas, label="Radius of defect")
    save(string(dir_df,"/1/Others/","ndefects_passive.png"), fig_N_pas)
    save(string(dir_df,"/1/Others/","radius_passive.png"), fig_R_pas)
    
    s_E_act_th = scatter!(ax_E_act, Vector{Float32}(pointx_E_act), Vector{Float32}(pointy_E_act), color = Vector{Float32}(color_E_act_th), colormap = :thermal, marker = Circle, alpha = 0.5, markersize = 24) #, colorrange=(minimum(color_E_pas), maximum(color_E_pas))
    s_E_act = scatter!(ax_E_act, Vector{Float32}(pointx_E_act), Vector{Float32}(pointy_E_act), color = Vector{Float32}(color_E_act), colormap = :thermal, marker=Circle)#, colorrange=(minimum(color_E_act_th), maximum(color_E_act_th)) )#, colorrange=(minimum(color_E_act_th), maximum(color_E_act_th)) )
    s_E_pas_th = scatter!(ax_E_pas, Vector{Float32}(pointx_E_pas), Vector{Float32}(pointy_E_pas), color = Vector{Float32}(color_E_pas_th), colormap = :thermal, marker = Circle, alpha = 0.5, markersize = 24) #, colorrange=(minimum(color_E_pas), maximum(color_E_pas))
    s_E_pas = scatter!(ax_E_pas, Vector{Float32}(pointx_E_pas), Vector{Float32}(pointy_E_pas), color = Vector{Float32}(color_E_pas), colormap = :thermal, marker=Circle) #, colorrange=(minimum(color_E_pas), maximum(color_E_pas))
    Colorbar(fig_E_act[1, 2], s_E_act, label="Simulation")
    Colorbar(fig_E_act[1, 3], s_E_act_th, label="Theory")
    Colorbar(fig_E_pas[1, 2], s_E_pas, label="Simulation")
    Colorbar(fig_E_pas[1, 3], s_E_pas_th, label="Theory")
    save(string(dir_df,"/1/Others/","epsi_active.png"), fig_E_act)
    save(string(dir_df,"/1/Others/","epsi_passive.png"), fig_E_pas)
    
    println(mean(color_E_pas_th./color_E_pas), ", ", mean(color_E_pas./color_E_pas_th))
    println(std(color_E_pas_th./color_E_pas), ", ", std(color_E_pas./color_E_pas_th))
    @show pointx_KAcr_act
end

# --- entry point -------------------------------------------------------------
# Renders the standard panels for every simulation in DATA_DIR.
if abspath(PROGRAM_FILE) == @__FILE__
    for idx in 1:N_sim
        Plot_per(idx)
    end
end
