include("../Utilities/using.jl")
using_pkg("CairoMakie, JLD, Printf, LaTeXStrings, FixedPointNumbers, DelimitedFiles, CSV, DataFrames, FileIO, Base64, Colors, LinearAlgebra, Statistics")
using_mod(".PictUtils")
CairoMakie.activate!(type = "png")

# ==============================================================================
#  CONFIGURATION  --  the only block that differs from the code run for the paper
# ==============================================================================
#  Figures for the paper were produced on a local machine after downloading the
#  cluster output, with these two lines hard-coded to:
#
#      dir_df  = homedir()*"/Data/P_noturnover/"
#      dir_fig = homedir()*"/Data/P_noturnover/"
#
#  DATA_DIR must contain DF.csv and one <idx>/Data/ folder per simulation, i.e.
#  the same layout InputParameters.jl writes. FIG_DIR receives the figures.
# ==============================================================================
dir_df = get(ENV, "DATA_DIR", abspath(joinpath(@__DIR__, "..", "data", "L50")) * "/")
dir_fig = get(ENV, "FIG_DIR", dir_df)

df = CSV.read(dir_df*"DF.csv", DataFrame)

N_sim = nrow(df)
# rm(dir_fig*"Figure_per/", recursive=true)

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
steppx = 10
steppz = 10
stepp = 15

# # function PlotFFT1D_per(idx; dir_df=dir_df, dir_fig=dir_fig)
# function Plot_per(idx)
#     # IJulia.clear_output(true)
#     # dir_df = "F:/Prot_FFT_2D_per/"*pa
#     # dir_fig = "F:/Prot_FFT_2D_per/"*pa
#     path = string("$(idx)/")
#     dirf = dir_fig*path*"Density/"
#     dirf2 = dir_fig*path*"Angle/"
#     # dirf3 = dir_fig*path*"Velocity_max/"
#     dirf3 = dir_fig*path*"Velocity/"
    
#     dird = dir_df*path*"Data/"
#     # dt = df[idx, :prin]*df[idx, :dt];
#     files_names = isdir(dird) ? readdir(dird) : []
#     mkpath(dirf)
#     mkpath(dirf2)
#     mkpath(dirf3)
#     num_files = length(files_names)
#     num_pict = length(filter!(endswith(".png"), readdir(dirf)));
#     println("Data in $dird")
#     println("Figure will be saved in $dirf")
#     println(num_files-num_pict, " / ", num_files) 
    
#     if num_files > 0 && num_files > num_pict
#         # try
#         Nx, Nz = size(load(string(dird,files_names[1]), "rho"))
#         x = 1:steppx:Nx
#         y = 1:steppz:Nz
#         hasP = haskey(load(string(dird,files_names[1])), "P")
#         hasQ = haskey(load(string(dird,files_names[1])), "Q")
#         # catch
#         # end
#         # x, y = meshgrid(1:stepp:Nz,1:stepp:Nx)
#         sttopp = num_files
#         # sttopp = (sttopp==0) ? sttopp+1 : sttopp
#         for i=1:sttopp
#             P_cpu = hasP ? load(string(dird,files_names[i]), "P") : nothing
#             P2_min = sqrt(minimum(P_cpu[:,:,1].^2 .+P_cpu[:,:,2].^2))
#             println(P2_min)
#             if !isfile(string(dirf,@sprintf("%04i",i),".png")) #|| hasQ
#                 # try
#                     density = load(string(dird,files_names[i]), "rho")

#                     v_cpu = load(string(dird,files_names[i]), "v")
#                     max_v = sqrt(maximum(v_cpu[:,:,1].^2 .+ v_cpu[:,:,2].^2))
#                     v_cpu ./= max_v
#                     mean_vx = mean(density.*v_cpu[:,:,1])
#                 vx_cpu = v_cpu[1:steppx:end,1:steppz:end,1]
#                 vz_cpu = v_cpu[1:steppx:end,1:steppz:end,2]


#                     P_cpu = hasP ? load(string(dird,files_names[i]), "P") : nothing
#                 Px_cpu = hasP ? P_cpu[1:steppx:end,1:steppz:end,1] : nothing
#                 Pz_cpu = hasP ? P_cpu[1:steppx:end,1:steppz:end,2] : nothing

#                     Q_cpu = hasQ ? load(string(dird,files_names[i]), "Q") : nothing
#                 Q1_cpu = hasQ ? Q_cpu[1:steppx:end,1:steppz:end,1] : nothing
#                 Q2_cpu = hasQ ? Q_cpu[1:steppx:end,1:steppz:end,2] : nothing
#                     S = hasQ ? sqrt.(Q1_cpu.^2 .+ Q2_cpu.^2) : nothing
#                     theta = hasQ ?  0.5*atan.(Q2_cpu.+1e-6, Q1_cpu) : atan.(Pz_cpu.+1e-6,Px_cpu)
#                     Qx = hasQ ? S.*cos.(theta) : nothing
#                     Qz = hasQ ? S.*sin.(theta) : nothing


#                     fig = Figure(size = (2000, 2000))
#                         ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
#                         ax1.title = string("Density and orientation " * string(i-1))
#                         h1 = heatmap!(ax1, density)
#                         try
#                             Colorbar(fig[1, 2], h1)
#                         catch
#                         end
#                         if hasQ
#                             arrow_size_Q = vec(S)
#                             arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 10.0, linewidth=2.0)
#                             # arrows!(ax1, x, y, Qx, Qz, arrowsize = arrow_size_Q, linecolor=:white, arrowhead=' ', lengthscale = 8.0, linewidth=2.0)
#                         end
#                         if hasP
#                             arrow_size_P = vec(norm.(Vec2f.(Px_cpu,Pz_cpu)))
#                             println(maximum(arrow_size_P))
#                             arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
#                             # arrows!(ax1, x, y, Px_cpu, Pz_cpu, arrowsize = arrow_size_P, linecolor=:black, lengthscale = 8.0, linewidth=2.0, arrowhead = '▲')
#                         end
#                     save(string(dirf,@sprintf("%04i",i),".png"), fig)

#                     fig2 = Figure(size = (1000, 1000))
#                         ax1 = Axis(fig2[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
#                         subs = hasQ ? " of Q " : " of P "
#                         ax1.title = string("Angle"* subs * string(i-1))
#                         theta2 = hasQ ?  0.5*atan.(Q_cpu[:,:,2].+1e-6, Q_cpu[:,:,1]) : atan.(P_cpu[:,:,2].+1e-6,P_cpu[:,:,1])
#                         h1 = heatmap!(ax1, theta2, colormap=:hsv)
#                         try
#                             Colorbar(fig2[1, 2], h1)
#                         catch
#                         end
#                     save(string(dirf2,@sprintf("%04i",i),".png"), fig2)

#                     fig3 = Figure(size = (1000, 1000))
#                         ax1 = Axis(fig3[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
#                         ax1.title = string("Velocity field " * string(i-1) * ", max v = " * string(max_v) * ", mean rho*vx = " * string(mean_vx))
#                         h1 = heatmap!(ax1, density)
#                         arrow_size_v = vec(norm.(Vec2f.(vx_cpu,vz_cpu)))
#                         arrows!(ax1, x, y, vx_cpu, vz_cpu, arrowsize = arrow_size_v*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
#                         try
#                             Colorbar(fig3[1, 2], h1)
#                         catch
#                         end
#                     save(string(dirf3,@sprintf("%04i",i),".png"), fig3)
#                 # catch
#                 # end
                
#             end
#         end
#         # pngstogif(dirf, dirf, string(idx), 1)
#         # pngstogif(dirf2, dirf2, string(idx), 1)

#     end
#     if num_files > 0 && ( num_files > num_pict || !isfile(dirf*string(idx)*".gif") )
#         pngstogif(dirf, dirf, string(idx), 8)
#         pngstogif(dirf2, dirf2, string(idx), 8)
#         length(readdir(dirf3)) > 0 ? pngstogif(dirf3, dirf3, string(idx), 8) : nothing
#     end
# end

function Plot_per(idx; dt = 1)
    # IJulia.clear_output(true)
    # dir_df = "F:/Prot_FFT_2D_per/"*pa
    # dir_fig = "F:/Prot_FFT_2D_per/"*pa
    path = string("$idx/")
    dirf = dir_fig*path*"Density/"
    dirf2 = dir_fig*path*"Angle/"
    # dirf3 = dir_fig*path*"Velocity_max/"
    dirf3 = dir_fig*path*"Velocity/"
    dirf4 = dir_fig*path*"Order/"
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
        for i=1:sttopp
            if !isfile(string(dirf,@sprintf("%04i",i),".png")) #|| hasQ
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


                    fig1 = Figure(size = (1000, 1025))
                        ax1 = Makie.Axis(fig1[2,1], aspect = DataAspect())
                        # ax1.title = string("Density and orientation " * string(t))
                        # h1 = heatmap!(ax1, 1:4:Nx, 1:4:Nz, density[1:4:end, 1:4:end], colormap=:RdPu)
                        h1 = heatmap!(ax1, density)
                        try
                            Colorbar(fig1[1, 1], h1, ticklabelsize=15, size=20, vertical = false)
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
                        colgap!(fig1.layout, 4)
                        rowgap!(fig1.layout, 0)
                    save(string(dirf,@sprintf("%04i",i),".png"), fig1)

                    fig2 = Figure(size = (1000, 1000))
                        ax1 = Makie.Axis(fig2[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                        subs = hasQ ? " of Q " : " of P "
                        ax1.title = string("Angle"* subs * string(t))
                        theta2 = hasQ ?  0.5*atan.(Q_cpu[:,:,2].+1e-6, Q_cpu[:,:,1]) : (hasP ? atan.(P_cpu[:,:,2].+1e-6,P_cpu[:,:,1]) : 0.0*similar(vx_cpu))
                        h1 = heatmap!(ax1, theta2, colormap=:hsv, colorrange = (-π, π))
                        try
                            Colorbar(fig2[1, 2], h1, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
                        catch
                        end
                    save(string(dirf2,@sprintf("%04i",i),".png"), fig2)

                    fig3 = Figure(size = (1000, 1000))
                        ax1 = Makie.Axis(fig3[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                        ax1.title = string("Velocity field " * string(t) * ", max v = " * string(max_v))
                        h1 = heatmap!(ax1, density)
                        arrow_size_v = vec(norm.(Vec2f.(vx_cpu,vz_cpu)))
                        arrows!(ax1, x, y, vx_cpu, vz_cpu, arrowsize = arrow_size_v*10.0, linecolor=:black, lengthscale = 10.0, linewidth=2.0, arrowhead = '▲', align = :center)
                        try
                            Colorbar(fig3[1, 2], h1)
                        catch
                        end
                    save(string(dirf3,@sprintf("%04i",i),".png"), fig3)
                
                    fig4 = Figure(size = (1000, 1000))
                        ax1 = Makie.Axis(fig4[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
                        subs = hasQ ? " of Q " : " of P "
                        ax1.title = string("Order"* subs * string(t))
                        ordre = hasQ ? sqrt.(Q_cpu[:,:,1].^2 .+ Q_cpu[:,:,2].^2) : (hasP ? sqrt.(P_cpu[:,:,1].^2 .+ P_cpu[:,:,2].^2) : 0.0*similar(vx_cpu))
                        h1 = heatmap!(ax1, ordre, colorrange = (0, 2))
                        try
                            Colorbar(fig4[1, 2], h1) #, ticks = ([-π, -π / 2, 0, π / 2, π], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]))
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
                        end
                    save(string(dirf4,@sprintf("%04i",i),".png"), fig4)
                
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
        pngstogif(dirf2, dirf2, string(idx), 8)
        length(readdir(dirf3)) > 0 ? pngstogif(dirf3, dirf3, string(idx), 8) : nothing
        length(readdir(dirf4)) > 0 ? pngstogif(dirf4, dirf4, string(idx), 8) : nothing
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


function plot_soutenance(idx; dt = 1)
    set_theme!(theme_black())
    dir_figg=dir_fig[1:end-1]*"b/"#"Z:/fig_soutenance_b/"

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
        stopp = num_files # == 21 ? num_files : num_files-1
        md=zeros(num_files,2)
        for i=1:stopp
            density = load(string(dird,files_names[i]), "rho")
            md[i,1] = minimum(density)
            md[i,2] = maximum(density)
        end
        minrho = minimum(md[:,1])
        maxrho = maximum(md[:,2])
        for i=2:4:stopp
            if !isfile(string(dirf,@sprintf("%04i",i),".png")) #|| hasQ
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
                        h1 = heatmap!(ax1, density, colorrange = (minrho, maxrho))
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
