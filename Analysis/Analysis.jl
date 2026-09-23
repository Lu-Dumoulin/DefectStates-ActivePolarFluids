include("../Utilities/using.jl")
using_pkg("CairoMakie, JLD, Printf, LaTeXStrings, FixedPointNumbers, DelimitedFiles, CSV, DataFrames, FileIO, Base64, Colors, LinearAlgebra, Statistics, Images, ImageSegmentation, FFTW, DelaunayTriangulation")
using_mod(".PictUtils")
CairoMakie.activate!(type = "png")
# CairoMakie.activate!(type = "svg")

# dir_df = "Z:/fts/"
# dir_fig = "Z:/fts/"

# dir_df = "Z:/FlowLattice/"
# dir_fig = "Z:/FlowLattice/"

# dir_df = "Z:/Det2/"
# dir_fig = "Z:/Det2/"

# dir_df = "Z:/LatticeAct/"
# dir_fig = "Z:/LatticeAct/"

# dir_df = "Z:/Aniso2/"
# dir_fig = "Z:/Aniso2/"

# dir_df = "E:/AdptDt/"
# dir_fig = "E:/AdptDt/"

# dir_df = "E:/AdptDt7D/"
# dir_fig = "E:/AdptDt7D/"

# dir_df = "Z:/2defects/"
# dir_fig = "Z:/2defects/"

# dir_df = "E:/tri7D/"
# dir_fig = "E:/tri7D/"

# dir_df = "E:/Coex7D/"
# dir_fig = "E:/Coex7D/"

# dir_df = "Z:/sat/"
# dir_fig = "Z:/sat/"

# dir_df = "Z:/hexatosquare/"
# dir_fig = "Z:/hexatosquare/"

# dir_df = "Z:/zetanu/"
# dir_fig = "D:/zetanu/"

# dir_df = "Z:/Dt/"
# dir_fig = "Z:/Dt/"

# dir_df = "D:/L/"
# dir_fig = "D:/L/"

# dir_df = "D:/small/"
# dir_fig = "D:/small/"

# dir_df = "D:/dx5/"
# dir_fig = "D:/dx5/"

# dir_df = "Z:/Qscan/"
# dir_fig = "D:/Qscan/"

# dir_df = "Z:/zetalamb/"
# dir_fig = "D:/zetalamb/"

# dir_df = "Z:/PQscan/"
# dir_fig = "D:/PQscan/"

# dir_df = "D:/nu_paper2/"
# dir_fig = "D:/nu_paper2/"

# dir_df = "D:/Anisov2/"
# dir_fig = "D:/Anisov2/"

dir_df = "F:/ZetaP_v2/"
dir_fig = "F:/ZetaP_v2/"

# dir_df = "D:/phasediag2/"
# dir_fig = "D:/phasediag2/"

# dir_df = "F:/sq5/"
# dir_fig = "F:/sq5/"

# dir_df = "F:/hx5/"
# dir_fig = "F:/hx5/"

# dir_df = "F:/L50/"
# dir_fig = "F:/L50/"

empty_theme = Theme(
    Axis = (
        backgroundcolor = :transparent,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xticklabelsvisible = false, 
        yticklabelsvisible = false,
        xgridcolor = :transparent,
        ygridcolor = :transparent,
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0.0,0.0),
        yautolimitmargin = (0.0,0.0),
    )
)

df = CSV.read(dir_df*"DF.csv", DataFrame)

N_sim = nrow(df)
# rm(dir_fig*"Figure_per/", recursive=true)

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
stepp = 15

G_avg = [0.475, 0.472, 0.461, 0.498, 0.522]
# function PlotFFT1D_per(idx; dir_df=dir_df, dir_fig=dir_fig)
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

# function PlotFFT2D_per_nov(pa=pa; gif=true)
#     dir_df = "F:/Prot_FFT_2D_per/"*pa
#     dir_fig = "F:/Prot_FFT_2D_per/"*pa
#     # path = df[idx, :fn]
#     # dirf = dir_fig*path*"Figure_per/"
#     dirf = dir_fig*"Figure_nov2/"
#     # dird = dir_df*path*"Data/"
#     dird = dir_df*"Data/"
#     # dt = df[idx, :prin]*df[idx, :dt];
#     # kd = df[idx, :kd];        
#     files_names = readdir(dird)
#     mkpath(dirf)
#     num_files = length(files_names)
#     num_pict = length(readdir(dirf));
#     println("Data in $dird")
#     println("Figure will be saved in $dirf")
#     println(num_files-num_pict)
#     if num_files > num_pict
#         Nx, Nz = size(load(string(dird,files_names[1]), "rho"))
#         for i=1:num_files
#             namepng = string(files_names[i][1:end-4],".png")
#             # if !isfile(string(dirf,@sprintf("%08i",i),".png"))
#             if !isfile(string(dirf,namepng))
#                 fig, axs = plt.subplots(2, 3, figsize=(25,12))
                
#                 ax1 = axs[1,1]
#                 ax5 = axs[1,2]
#                 ax6 = axs[1,3]
#                 ax2 = axs[2,1]
#                 ax7 = axs[2,2]
#                 ax8 = axs[2,3]
                
#                 density = load(string(dird,files_names[i]), "rho")
#                 v_cpu = load(string(dird,files_names[i]), "v")
#                 vx_cpu = v_cpu[:,:,1]
#                 vz_cpu = v_cpu[:,:,2]
#                 P_cpu = load(string(dird,files_names[i]), "P")
#                 Px_cpu = P_cpu[:,:,1]
#                 Pz_cpu = P_cpu[:,:,2]
#                 Q_cpu = load(string(dird,files_names[i]), "Q")
#                 Q1_cpu = Q_cpu[:,:,1]
#                 Q2_cpu = Q_cpu[:,:,2]
#                 # v = load(string(dird,files_names[i]), "velocity")
#                 PyPlot.title(string(i-1))
#                 ax1.set_title("Density")
#                 h1 = ax1.pcolormesh(density')
#                 colorbar(h1, ax=ax1)
                
#                 ax5.set_title("Px")
#                 h5 = ax5.pcolormesh(Px_cpu')
#                 colorbar(h5, ax=ax5)
                
#                 ax6.set_title("Pz")
#                 h6 = ax6.pcolormesh(Pz_cpu')
#                 colorbar(h6, ax=ax6)
                
#                 ax7.set_title("Q1")
#                 h7 = ax7.pcolormesh(Q1_cpu')
#                 colorbar(h7, ax=ax7)
                
#                 ax8.set_title("Q2")
#                 h8 = ax8.pcolormesh(Q2_cpu')
#                 colorbar(h8, ax=ax8)
                
#                 # savefig(string(dirf,@sprintf("%04i",i),".png"))
#                 savefig(string(dirf,namepng))
#                 close("all")
#             end
#         end
#     end
#     # pngstogif(dirf, dir_fig*path*"/", string(1), 10)
#     if gif 
#         pngstogif(dirf, dir_fig*"/", string("nov"), 10, Nimg=101)
#     end
# end

# include("MakePlots.jl")
# PlotFFT2D_per_nov("PQ-clu-fft-big-23/", gif=false)
# PlotFFT2D_per_nov("PQ-clu-fft-big-23m/", gif=false)
# PlotFFT2D_per_nov("PQ-clu-fft-big-f10/", gif=false)

function count_defect(Px, Py, size_defect)
    Nx, Ny = size(Px)
    order = sqrt.(Px.^2 .+ Py.^2)
    test = order .* (order .< 0.25*maximum(order)) 
    nzer = findall(!iszero, test)
    test = nothing
    lnzer = length(nzer)
    if lnzer > 0#1000#0
        return use_segmentation(order, Nx, Ny)
    end
    dist_matrix = zeros(lnzer,lnzer)
    idxnzer = hcat(getindex.(nzer, 1), getindex.(nzer,2))
    for i=1:lnzer, j=1:lnzer
        x1 = idxnzer[i,1]; x2 = idxnzer[j,1]
        y1 = idxnzer[i,2]; y2 = idxnzer[j,2]
        x1 = abs(x1-x2) > Nx*0.5 ? x1 = x1 + sign(x2-x1)*Nx : x1
        y1 = abs(y1-y2) > Ny*0.5 ? y1 = y1 + sign(y2-y1)*Ny : y1
        dist_matrix[i,j] = sqrt((x1-x2)^2+(y1-y2)^2)
    end
    alclust= Vector{Vector{Any}}()#[[]] 
    for i=1:lnzer
        idxci = 0
        for k=1:length(alclust)
            if i ∈ alclust[k]
                idxci += k
            end
        end
        for j=i+1:lnzer
            if dist_matrix[i,j] < size_defect
                idxcj = 0
                for k=1:length(alclust)
                    if j ∈ alclust[k]
                        idxcj += k
                    end
                end
                idxcj == idxci == 0 ? push!(alclust, [i, j]) : nothing
                idxci == idxcj != 0 ? break : nothing
                idxcj != idxci == 0 ? push!(alclust[idxcj], i) : nothing
                idxcj != idxci == 0 ? idxci = idxcj : nothing
                idxci != idxcj == 0 ? push!(alclust[idxci], j) : nothing
                if idxci != idxcj != 0 
                    push!(alclust, vcat(alclust[idxci],alclust[idxcj]))
                    idxci > idxcj ? deleteat!(alclust, idxci) : deleteat!(alclust, idxcj)
                    idxci > idxcj ? deleteat!(alclust, idxcj) : deleteat!(alclust, idxci)
                    idxci = length(alclust)
                end
                idxcj == idxci == 0 ? idxci = length(alclust) : nothing
            end
        end
    end
    Ndefect = length(alclust)
    defect_idx = zeros(Ndefect,2)
    for i=1:Ndefect
        tmp = [order[idxnzer[j,1],idxnzer[j,2]] for j in alclust[i] ]
        defect_idx[i, :] = idxnzer[alclust[i][findmin(tmp)[2]],:]
    end
    order = nothing
    test = nothing
    nzer = nothing
    dist_matrix = nothing
    idxnzer = nothing
    alclust = nothing
    
    return defect_idx
    
#     # Clear if copy of one
#     dist_matrix = zeros(Ndefect, Ndefect)
#     for i=1:Ndefect, j=i+1:Ndefect
#         x1 = defect_idx[i,1]; x2 = defect_idx[j,1]
#         y1 = defect_idx[i,2]; y2 = defect_idx[j,2]
#         x1 = abs(x1-x2) > Nx*0.5 ? x1 = x1 + sign(x2-x1)*Nx : x1
#         y1 = abs(y1-y2) > Ny*0.5 ? y1 = y1 + sign(y2-y1)*Ny : y1
#         dist_matrix[i,j] = sqrt((x1-x2)^2+(y1-y2)^2)
#     end
#     idx_of_unique = []
#     for i=1:Ndefect
#         doublon = findall(x-> ((x<size_defect+5) & (x!=0)), dist_matrix[i,:])
#         if length(doublon) >= 1
#             for j in doublon
#                 i < j ? push!(idx_of_unique, i) : push!(idx_of_unique, j)
#             end
#         else
#             push!(idx_of_unique, i)
#         end
#     end
            
#     return defect_idx[unique!(idx_of_unique), :]
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

function cara_defects(defects_pos, theta2)
    Nx, Ny = size(theta2)
    ndefects = length(defects_pos[:,1])
    box_size = 5
    length_ = box_size-1
    half_size = round(Int, length_*0.5)
    defect_charge=zeros(ndefects)
    for d=1:ndefects
        cx, cy = Int.(defects_pos[d,:])
        cycle = zeros(4*length_)
        for i=0:length_-1
            xi = mod(cx-half_size+i,1:Nx)
            yi = mod(cy-half_size, 1:Ny)
            cycle[i+1]= theta2[xi,yi]
            xi = mod(cx+half_size,1:Nx)
            yi = mod(cy-half_size+i,1:Ny)
            cycle[i+1+length_] = theta2[xi,yi]
            xi = mod(cx+half_size-i,1:Nx)
            yi = mod(cy+half_size, 1:Ny)
            cycle[i+1+2*length_]= theta2[xi,yi]
            xi = mod(cx-half_size,1:Nx)
            yi = mod(cy+half_size-i,1:Ny)
            cycle[i+1+3*length_] = theta2[xi,yi]
        end
        defect_charge[d] = sum(diff(sort!(cycle, rev = sign(sum(sign.(diff(cycle)))) >0 )))/(2*π)
    end
    return defect_charge
end
    
# cd("../FFT_2D_Clean/")
# include("../Utilities/using.jl")
# using_pkg("CairoMakie, JLD, Printf, LaTeXStrings, FixedPointNumbers, DelimitedFiles, CSV, DataFrames, FileIO, Base64, Colors, LinearAlgebra")
# using_mod(".PictUtils")
# CairoMakie.activate!(type = "png")
# # CairoMakie.activate!(type = "svg")

# # dir_df = "F:/2D_P_Q_PQ_kplonger/"
# # dir_fig = "F:/2D_P_Q_PQ_kplonger/"

# dir_df = "F:/rho2long/"
# dir_fig = "F:/rho2long/"

# df = CSV.read(dir_df*"DF.csv", DataFrame)

# N_sim = nrow(df)
# # rm(dir_fig*"Figure_per/", recursive=true)

# meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
# stepp = 15
# idx = 306
# path = string("$idx/")
#     dird = dir_df*path*"Data/"
#     # dt = df[idx, :prin]*df[idx, :dt];
#     files_names = isdir(dird) ? readdir(dird) : []

#     num_files = length(files_names)
# Nx, Nz = size(load(string(dird,files_names[1]), "rho"))
#         x = 1:stepp:Nx
#         y = 1:stepp:Nz
#         hasP = haskey(load(string(dird,files_names[1])), "P")
#         hasQ = haskey(load(string(dird,files_names[1])), "Q")
#  P_cpu = hasP ? load(string(dird,files_names[61]), "P") : nothing
# Px = P_cpu[:,:,1]
# Py = P_cpu[:,:,2]

function all_N_and_R(N_sim=N_sim)
    Ndefects = zeros(N_sim)
    Radius = zeros(N_sim)
    Epsi = zeros(N_sim)
    for idx=1:N_sim
        # println(idx)
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        
        files_names = isdir(dird) ? readdir(dird) : []

        num_files = length(files_names)
       
        hasP = haskey(load(string(dird,files_names[end])), "P")
        hasQ = haskey(load(string(dird,files_names[end])), "Q")

    
        i=num_files
        P_cpu = hasP ? load(string(dird,files_names[i]), "P") : nothing
        Px = P_cpu[:,:,1]
        Py = P_cpu[:,:,2]
        theta2 = atan.(Py.+1e-6,Px)
        P_cpu = nothing
        Nx, Ny = size(Px)
        S = Nx*Ny*1e-2*1e-2
        defects_pos = count_defect(Px, Py, 10)
        if !isnothing(defects_pos)
            Ndefects[idx] = length(defects_pos[:,1])
            rho = load(string(dird,files_names[end]), "rho")
            rho_0 = Statistics.mean(rho)
            rho_in = minimum(rho)
            rho_out = maximum(rho)
            Radius[idx] = sqrt(S*(rho_0-rho_out)/(rho_in-rho_out)/Ndefects[idx]/π)
            Epsi[idx] = rho_0-rho_in
        end
        println(idx, " ", Ndefects[idx]," ", Radius[idx], " ", Epsi[idx])
    end
    return Ndefects, Radius, Epsi
end

kdoar(ρ0, act, r) = 6/r^2*(ρ0^3-0.75*act*ρ0^2)/(3*0.1*ρ0^2/1e-4-1)

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

function use_segmentation(order, Nx, Nz)
    # Make binary mask
    bw = (order .< 0.25*maximum(order)) 
   
    # Segment defects
    dist = 1 .- distance_transform(feature_transform(bw));
    seg = felzenszwalb(dist .< -8, 5, 200)

    # Change label borders
    # lmap = labels_map(seg)
    # for i=1:Nx
    #     if (lmap[i,1] != 1) & (lmap[i,Nz] != 1) && lmap[i,1] != lmap[i,Nz] 
    #         lmap[findall(x->x==lmap[i,Nz], lmap)] .= lmap[i,1]
    #     end
    # end
    # for j=1:Nz
    #     if (lmap[1,j] != 1) & (lmap[Nx,j] != 1) && lmap[1,j] != lmap[Nx,j]
    #         lmap[findall(x->x==lmap[Nx,j], lmap)] .= lmap[1,j]
    #     end
    # end
    
    # Remove one seg at border
    # lmap = labels_map(seg)
    #################################################################################################################################
    # for i=1:Nx
    #     if (labels_map(seg)[i,1] != 1) & (labels_map(seg)[i,Nz] != 1) && labels_map(seg)[i,1] != labels_map(seg)[i,Nz]
    #         spc1 = segment_pixel_count(seg, labels_map(seg)[i,1])
    #         spc2 = segment_pixel_count(seg, labels_map(seg)[i,Nz])
    #         if spc1 >= spc2 
    #             rem_segment!(seg, labels_map(seg)[i,Nz], ret1)
    #         else
    #             rem_segment!(seg, labels_map(seg)[i,1], ret1)
    #         end
    #     end
    # end
    # for j=1:Nz
    #     if (labels_map(seg)[1,j] != 1) & (labels_map(seg)[Nx,j] != 1) && labels_map(seg)[1,j] != labels_map(seg)[Nx,j]
    #         spc1 = segment_pixel_count(seg, labels_map(seg)[1,j])
    #         spc2 = segment_pixel_count(seg, labels_map(seg)[Nx,j])
    #         if spc1 >= spc2 
    #             rem_segment!(seg, labels_map(seg)[Nx,j], ret1)
    #         else
    #             rem_segment!(seg, labels_map(seg)[1,j], ret1)
    #         end
    #     end
    # end
    #################################################################################################################################
    
    nseg = length(segment_labels(seg))-1
    defect_idx = zeros(nseg, 2)
    # Find position of defects (center of seg)
    for i=2:length(segment_labels(seg))
        points = findall(x->x==segment_labels(seg)[i], labels_map(seg))
        xc = mean([p[1] for p in points])
        yc = mean([p[2] for p in points])
        defect_idx[i-1,1] = round(Int,xc)
        defect_idx[i-1,2] = round(Int,yc)
    end
    # imshow(map(i->get_random_color(i), labels_map(seg)) .* (1 .-bw))
    return defect_idx
end

ret1(x,y) = 1.0

function get_spectra(N_sim=N_sim)
    path = string("1/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    rho = load(string(dird,files_names[1]), "rho")
    Nx, Ny = size(rho)
    hN = Int(Nx/2)
    
    W = plan_fft(ones(Float64, Nx, Ny))
    # Wi = inv(W)
    kx = 2*pi*fftfreq(Nx, 1/1e-2)
    kz = 2*pi*fftfreq(Ny, 1/1e-2)
    tab_l = zeros(N_sim,2,length(files_names))
    for idx = 1:N_sim
        # lP, lrho = get_spectrum(idx, W, kx, kz)
        lP, lrho = get_Sk(idx, W, kx, kz, hN)
        println(idx, " ", lP, " ", lrho)
        tab_l[idx,1,1:length(lP)] .= lP
        tab_l[idx,2,1:length(lrho)] .= lrho
    end
    return tab_l
end

function get_Sk(idx, W, kx, kz, hN)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    # dt = df[idx, :prin]*df[idx, :dt];
    files_names = isdir(dird) ? readdir(dird) : []

    num_files = length(files_names)
       
    hasP = haskey(load(string(dird,files_names[1])), "P")
    hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
    dirf6 = dir_fig*path*"Others/"
    mkpath(dirf6)
    lto = zeros(num_files)
    ltr = zeros(num_files)
    for file in 1:num_files
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
        for i=2:hN-1, j=2:hN-1
            sko[i,j] = abs(fo[i,j])^2
            skr[i,j] = abs(fr[i,j])^2
        end

        fig = Figure(size = (1000, 1000))
        ax1 = Makie.Axis(fig[1,1], xlabel = L"k_x", ylabel = L"k_y")#, aspect = DataAspect())
        ax1.title = string("Structure factor Order")
        heatmap!(ax1, sko[1:100, 1:100])
        save(string(dirf6,@sprintf("order_%04i",file),".png"), fig)

        fig = Figure(size = (1000, 1000))
        ax1 = Makie.Axis(fig[1,1], xlabel = L"k_x", ylabel = L"k_y")#, aspect = DataAspect())
        ax1.title = string("Structure factor Density")
        heatmap!(ax1, skr[1:100, 1:100])
        save(string(dirf6,@sprintf("density_%04i",file),".png"), fig)
    
        lto[file] = 1/sqrt(kx[Int(findmax(sko[2,:])[2])]^2 + kz[2]^2)
        ltr[file] = 1/sqrt(kx[Int(findmax(skr[2,:])[2])]^2 + kz[2]^2)
    end
    fig = Figure(size = (1000, 1000))
        ax1 = Makie.Axis(fig[1,1], xlabel = L"t", ylabel = L"length")#, aspect = DataAspect())
        ax1.title = string("l(t) Order")
        lines!(ax1, lto)
        save(string(dirf6,"order_lt",".png"), fig)
    fig = Figure(size = (1000, 1000))
        ax1 = Makie.Axis(fig[1,1], xlabel = L"t", ylabel = L"length")#, aspect = DataAspect())
        ax1.title = string("l(t) Density")
        lines!(ax1, ltr)
        save(string(dirf6,"density_lt",".png"), fig)
    return lto, ltr
end

# function get_spectrum(idx, W, kx, kz)
#     path = string("$idx/")
#     dird = dir_df*path*"Data/"
#     # dt = df[idx, :prin]*df[idx, :dt];
#     files_names = isdir(dird) ? readdir(dird) : []

#     num_files = length(files_names)
       
#     hasP = haskey(load(string(dird,files_names[1])), "P")
#     hasQ = haskey(load(string(dird,files_names[1])), "Q")
    
#     dirf6 = dir_fig*path*"Others/"
#     mkpath(dirf6)
    
#     P_cpu = hasP ? load(string(dird,files_names[num_files]), "P") : nothing
#     rho = load(string(dird,files_names[num_files]), "rho")
#     Px = P_cpu[:,:,1]
#     Py = P_cpu[:,:,2]
#     Nx, Ny = size(Px)
#     order = sqrt.(Px.^2 .+ Py.^2)
#     fo = (W * order)
#     # afo = abs.((W * order))
#     afo[1,1] = 0.0
#     # ii = findmax(afo[1,1:100])
#     fig = Figure(size = (1000, 1000))
#     ax1 = Makie.Axis(fig[1,1], xlabel = L"k_x", ylabel = L"k_y")#, aspect = DataAspect())
#     ax1.title = string("Wave number Order")
#     heatmap!(ax1, afo[:,1:floor(Int,Ny/2)])
#     save(string(dirf6,@sprintf("%04i",num_files),".png"), fig)
#     ii = Int(findmax(afo[1,1:floor(Int,Ny/2)])[2])
    
#     afr = abs.((W * rho))
#     afr[1,1] = 0
#     fig = Figure(size = (1000, 1000))
#     ax1 = Makie.Axis(fig[1,1], xlabel = L"k_x", ylabel = L"k_y")#, aspect = DataAspect())
#     ax1.title = string("Wave number Density")
#     heatmap!(ax1, afr[:,1:floor(Int,Ny/2)])
#     save(string(dirf6,@sprintf("density_%04i",num_files),".png"), fig)
#     iii = Int(findmax(afr[1,1:floor(Int,Ny/2)])[2])
#     return 1/kx[ii], 1/kx[iii]
# end
    

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

# function make_passive_fig(list_idx)
#     for idx=1:N_sim
#         zp = df[idx, :zetap]
#         g = df[idx, :gamma]
#         nu = df[idx, :nu1]
#         if zp == 0 && g == 1 && nu == 0 
#             kd = df[idx, :kd]/1e-3
#             ar = df[idx, :ar]
#             a = df[idx, :a]
#             ρ0 = df[idx, :rho0]
#             if ρ0 > 0.3 && !( (a > 0) & (kd == 0) ) && !( (ar > 2) & (kd == 0) )
#     fig = Figure(resolution=(1200,1000))
#     axs = [Axis(fig[i,j]; aspect = DataAspect()) for i in 1:6 for j in 1:5]
#     for i in 1:30
#         idx = list_idx[i]
#         path = string("$idx/"); 
#         dird = dir_df*path*"Data/"; 
#         files_names = isdir(dird) ? readdir(dird) : return 0; 
#         ss_namefile = string(dird,files_names[end])
#         density = load(ss_namefile, "rho")
#         heatmap!(axs[i], density)
#     end
#     hidedecorations!.(axs)
#     colgap!(fig.layout, 1)
#     rowgap!(fig.layout,1)
#     hidespines!.(axs)
#     return fig
# end

function expend_periodic_array(ar, ex)
    Nx, Ny = size(ar)
    nNx = Nx + 2*ex
    nNy = Ny + 2*ex
    nar = zeros(Float64, nNx, nNy)
    # Mid
    nar[1+ex:Nx+ex, 1+ex:Ny+ex] .= ar
    # Left
    nar[1+ex:Nx+ex, 1:ex] .= ar[:,Ny-ex+1:Ny]
    # Right
    nar[1+ex:Nx+ex, nNy-ex+1:nNy] .= ar[:,1:ex]
    # Top
    nar[1:ex, 1+ex:Ny+ex] .= ar[Nx-ex+1:Nx,:]
    # Bot
    nar[Nx+1+ex:nNx, 1+ex:Ny+ex] .= ar[1:ex,:]
    # Top Left
    nar[1:ex,1:ex] .= ar[Nx-ex+1:Nx,Ny-ex+1:Ny]
    # Bot Right
    nar[Nx+ex+1:nNx,Ny+ex+1:nNy] .= ar[1:ex,1:ex]
    # Top Right
    nar[1:ex,Ny+ex+1:nNy] .= ar[Nx-ex+1:Nx,1:ex]
    # Bot Left
    nar[Nx+ex+1:nNx,1:ex] .= ar[1:ex,Ny-ex+1:Ny]
return nar
end

function expend_periodic_array2(ar, ex)
    Nx, Nz = size(ar)
    nNx = Nx + 2*ex
    nNz = Nz + 2*ex
    inX = ex+1:Nx+ex; leftn = 1:ex; rightn = Nx+ex+1:nNx
    inZ = ex+1:Nz+ex; botn = 1:ex; topn = Nz+ex+1:nNz
    left = 1:ex; right = Nx-ex+1:Nx
    bot = 1:ex; top = Nz-ex+1:Nz
    nar = zeros(Float64, nNx, nNz)
    # Mid
    nar[inX, inZ] .= ar
    # Left In
    nar[leftn, inZ] .= ar[right,:]
    # nar[1+ex:Nx+ex, 1:ex] .= ar[:,Ny-ex+1:Ny]
    # Right
    nar[rightn, inZ] .= ar[left,:]
    # nar[1+ex:Nx+ex, nNy-ex+1:nNy] .= ar[:,1:ex]
    # Top
    nar[inX, topn] .= ar[:, bot]
    # nar[1:ex, 1+ex:Ny+ex] .= ar[Nx-ex+1:Nx,:]
    # Bot
    nar[inX, botn] .= ar[:, top]
    # nar[Nx+1+ex:nNx, 1+ex:Ny+ex] .= ar[1:ex,:]
    # Top Left
    nar[leftn, topn] .= ar[right, bot]
    # nar[1:ex,1:ex] .= ar[Nx-ex+1:Nx,Ny-ex+1:Ny]
    # Bot Right
    nar[rightn, botn] .= ar[left, top]
    # nar[Nx+ex+1:nNx,Ny+ex+1:nNy] .= ar[1:ex,1:ex]
    # Top Right
    nar[rightn, topn] .= ar[left, bot]
    # nar[1:ex,Ny+ex+1:nNy] .= ar[Nx-ex+1:Nx,1:ex]
    # Bot Left
    nar[leftn, botn] .= ar[right, top]
    # nar[Nx+ex+1:nNx,1:ex] .= ar[1:ex,Ny-ex+1:Ny]
return nar
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

function all_G6pm(N = N_sim, N_ini = 1)
    df2 = DataFrame()
    if isfile(joinpath(dir_df,"DF2-pm.csv")) 
        df2 = CSV.read(dir_df*"DF2-pm.csv", DataFrame)
    elseif isfile(joinpath(dir_df,"DF2.csv")) 
        df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
        df2[!, :g2pm] .= NaN
        df2[!, :g3pm] .= NaN
        df2[!, :g4pm] .= NaN
        df2[!, :g5pm] .= NaN
        df2[!, :g6pm] .= NaN
        df2[!, :theta2pm] .= NaN
        df2[!, :theta3pm] .= NaN
        df2[!, :theta4pm] .= NaN
        df2[!, :theta5pm] .= NaN
        df2[!, :theta6pm] .= NaN
    end

    for idx=N_ini:N
        print(idx, ", ")
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        files_names = isdir(dird) ? readdir(dird) : []
        num_files = length(files_names)
        dird5 = dir_fig*path*"Defects_pos/"
        ex = 200
        file = num_files
        ptsp = Array{Float64,2}
        Nx, Ny = size(load(string(dird,files_names[file]), "rho"))
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
            vornm = voronoi(trim, clip=true)
            map_pointsm = vornm.generators
            polysm = vornm.polygons
            boundarym = vornm.boundary_polygons
            Npolym = length(polysm)
            # Positive defects
            trip = triangulate(ptsp')
            vornp = voronoi(trip, clip=true)
            map_pointsp = vornp.generators
            polysp = vornp.polygons
            boundaryp = vornp.boundary_polygons
            Npolyp = length(polysp)
            
            aγ = zeros(Npolym+Npolyp,5)
            aθ = zeros(Npolym+Npolyp,5)

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
                print(p, ", ")
                for i=1:Npolym
                    if !(i in boundarym)
                    c = DelaunayTriangulation.get_centroid(vornm, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

                            Δp = sum(rθ[i,:,1].^p)
                            γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
                            aγ[i,p-1] = abs.(γ)
                            aθ[i,p-1] = angle.(γ)/p
                        end
                    end
                end
                for i=1:Npolyp
                    if !(i in boundaryp)
                    c = DelaunayTriangulation.get_centroid(vornp, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1

                            Δp = sum(rθ[Npolym+i,:,1].^p)
                            γ = 1/Δp*sum(rθ[Npolym+i,:,1].^p .* exp.(im*p*rθ[Npolym+i,:,2]) )
                            aγ[Npolym+i,p-1] = abs.(γ)
                            aθ[Npolym+i,p-1] = angle.(γ)/p
                        end
                    end
                end
            end
            println()
            df2[idx, :g2pm] = mean(filter!(x -> x != 0, aγ[:,1]))
            df2[idx, :g3pm] = mean(filter!(x -> x != 0, aγ[:,2]))
            df2[idx, :g4pm] = mean(filter!(x -> x != 0, aγ[:,3]))
            df2[idx, :g5pm] = mean(filter!(x -> x != 0, aγ[:,4]))
            df2[idx, :g6pm] = mean(filter!(x -> x != 0, aγ[:,5]))
            df2[idx, :theta2pm] = mean(filter!(x -> x != 0, aθ[:,1]))
            df2[idx, :theta3pm] = mean(filter!(x -> x != 0, aθ[:,2]))
            df2[idx, :theta4pm] = mean(filter!(x -> x != 0, aθ[:,3]))
            df2[idx, :theta5pm] = mean(filter!(x -> x != 0, aθ[:,4]))
            df2[idx, :theta6pm] = mean(filter!(x -> x != 0, aθ[:,5]))
            CSV.write(joinpath(dir_df,"DF2-pm.csv"), df2)
        end
    end
    # CSV.write(joinpath(dir_df,"DF2-pm.csv"), df2)
end

function all_N_R_G(N = N_sim, N_ini = 1, ex=200)
    df2 = DataFrame()
    if isfile(joinpath(dir_df,"DF2.csv")) 
        df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    else
        df2 = DataFrame(df)
        df2[!, :Ndefects] .= NaN
        df2[!, :epsi] .= NaN
        df2[!, :Rad] .= NaN
        df2[!, :g2] .= NaN
        df2[!, :g3] .= NaN
        df2[!, :g4] .= NaN
        df2[!, :g5] .= NaN
        df2[!, :g6] .= NaN
        df2[!, :theta2] .= NaN
        df2[!, :theta3] .= NaN
        df2[!, :theta4] .= NaN
        df2[!, :theta5] .= NaN
        df2[!, :theta6] .= NaN
    end

    for idx=N_ini:N
        path = string("$idx/")
        dird = dir_df*path*"Data/"
        dird5 = dir_fig*path*"Defects_pos/"
        files_names = isdir(dird) ? readdir(dird) : []
        files_name5 = isdir(dird5) ? last(readdir(dird5)) : []
        num_files = length(files_names)
        if num_files > 5
            pts = Array{Float64,2}
            ptsp = Array{Float64,2}
            Nx, Ny = size(load(string(dird,files_names[num_files]), "rho"))

            if !isfile(string(dird5,files_name5))
                println("Call plt_defect($idx) first")
                return 0
            end
            defects_pos = load(string(dird5,files_name5), "pos")
#             si = load(string(dird5,files_name5,".jld"), "charges")
#             it = findall(x-> x<-0.2, si)
#             ptsm = zeros(length(it), 2)
#             for i in 1:length(it)
#                 ptsm[i,:] .= Int.(pts[it[i],:])
#             end

#             it = findall(x->x>0.2, si)
#             ptsp = zeros(length(it), 2)
#             for i in 1:length(it)
#                 ptsp[i,:] .= Int.(pts[it[i],:])
#             end
            
            defects_pos_unique = load(string(dird5,files_name5), "pos_unique")
            Ndef_u = length(defects_pos_unique[:,1])
            rho = load(string(dird,files_names[num_files]), "rho")
            rho_0 = Statistics.mean(rho)
            rho_in = minimum(rho)
            rho_out = maximum(rho)
            df2[idx, :epsi] = rho_0-rho_in
            if length(defects_pos) > 2#500000
                tri = triangulate(defects_pos')
                vorn = voronoi(tri, clip=true)
                map_points = vorn.generators
                polys = vorn.polygons
                boundary = vorn.boundary_polygons
                Npoly = length(polys)

                rθ = zeros(Npoly, 20, 2)
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

                aγ = zeros(Npoly,5)
                aθ = zeros(Npoly,5)
                for i=1:Npoly
                    if !(i in boundary)
                    c = DelaunayTriangulation.get_centroid(vorn, i)
                        if c[1] > ex && c[2]>ex && c[1]<Nx+ex+1 && c[2]<Ny+ex+1
                            for p=2:6
                                Δp = sum(rθ[i,:,1].^p)
                                γ = 1/Δp*sum(rθ[i,:,1].^p .* exp.(im*p*rθ[i,:,2]) )
                                aγ[i,p-1] = abs.(γ)
                                aθ[i,p-1] = angle.(γ)/p
                            end
                        end
                    end
                end
                df2[idx, :g2] = mean(filter!(x -> x != 0, aγ[:,1]))
                df2[idx, :g3] = mean(filter!(x -> x != 0, aγ[:,2]))
                df2[idx, :g4] = mean(filter!(x -> x != 0, aγ[:,3]))
                df2[idx, :g5] = mean(filter!(x -> x != 0, aγ[:,4]))
                df2[idx, :g6] = mean(filter!(x -> x != 0, aγ[:,5]))
                df2[idx, :theta2] = mean(filter!(x -> x != 0, aθ[:,1]))
                df2[idx, :theta3] = mean(filter!(x -> x != 0, aθ[:,2]))
                df2[idx, :theta4] = mean(filter!(x -> x != 0, aθ[:,3]))
                df2[idx, :theta5] = mean(filter!(x -> x != 0, aθ[:,4]))
                df2[idx, :theta6] = mean(filter!(x -> x != 0, aθ[:,5]))
                df2[idx, :Rad] = sqrt(Nx*Ny*(rho_0-rho_out)/(rho_in-rho_out)/Ndef_u/π)
            end
            df2[idx, :Ndefects] = Ndef_u
            # df2[idx, :g4] = mgamma
            # df2[idx, :theta4] = mangle
            println(idx, ", ", Ndef_u)
        end 
    end
    CSV.write(joinpath(dir_df,"DF2.csv"), df2)
    # return df2
end

   
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
    
#     isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
#     save("Z:/Fig_paper/Ndefect.png", fig)
    
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","Ndefb.csv"), ndef)
    
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
    
#     isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
#     save("Z:/Fig_paper/Gamma4.png", fig2)
    
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
    
#     isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
#     save("Z:/Fig_paper/Gamma3.png", fig3)
#     return fig, fig2, fig3
end
    
function comp_correlation()
    df2 = CSV.read(dir_df*"DF2.csv", DataFrame)
    dd = df2[:,:Ndefects] .> 150
    n = sum(dd)
    Nd = zeros(n)
    G4 = zeros(n)
    i=0
    for idx=1:nrow(df2)
        if dd[idx] == 1
            Nd[i] = df2[idx,:Ndefects]
            G4[i] = df2[idx,:g4]
            i+=1
        end
    end
    return cor(Nd, G4)
end
    
# hist(filter(!isequal(0), df2[:,:theta4] .* .!isnan.(df2[:,:theta4]) ))
# hist(filter(!isequal(0), df2[:,:g4] .* .!isnan.(df2[:,:g4]) ))
# hist(filter(!isequal(0), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects]) ))

# j = 1
# tbl = Array(50:25:350)
# tbu = Array(400:-25:100)
# ltbl = length(tbl)
# ltbu = length(tbu)
# nn = ltbl*ltbu
# tab = zeros(ltbl,ltbu,5)
# for l=1:ltbl, u=1:ltbu
#     bu = tbu[u]; bl = tbl[l]
#     if bl < bu
#         dd = (df2[:,:Ndefects] .< bu) .* (df2[:,:Ndefects] .> bl)
#         n = sum(dd)
#         Nd = zeros(n)
#         gg = zeros(n,5)
#         i=1
#         for idx=1:nrow(df2)
#             if dd[idx] == 1
#                 Nd[i] = df2[idx,:Ndefects]
#                 gg[i,5] = df2[idx,:g6]
#                 gg[i,4] = df2[idx,:g5]
#                 gg[i,3] = df2[idx,:g4]
#                 gg[i,2] = df2[idx,:g3]
#                 gg[i,1] = df2[idx,:g2]
#                 i+=1
#             end
#         end
#         tab[l,u,1] = cor(Nd, gg[:,1])
#         tab[l,u,2] = cor(Nd, gg[:,2])
#         tab[l,u,3] = cor(Nd, gg[:,3])
#         tab[l,u,4] = cor(Nd, gg[:,4])
#         tab[l,u,5] = cor(Nd, gg[:,5])
#         println(l, " ", u, " ", cor(Nd, gg[:,2]))
#     end
# end
# f = Figure(resolution=(1800,1000))
#     axs = [j*i==6 ? nothing : Makie.Axis(f[i,j+1]; aspect = DataAspect(), xaxisposition = i==1 ? :top : :bottom, yaxisposition = j==3 ? :right : :left, xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30) for i in 1:2 for j in 1:3] 
# for idx in 1:5
#     heatmap!(axs[idx], tbl, tbu, tab[:,:,idx], colormap = :vik, colorrange=(-1,1))
#     text!(axs[idx], 50, 375, text = L"p = %$(idx+1)", color = :white, font = :bold, fontsize = 24)
# end
# Colorbar(f[:,1], colormap=:vik, limits=(-1,1), label= "Correlation", labelsize = 30, ticklabelsize=24)
# f


# # @show mean(df2[findall(x->x>0.65, df2[:,:g2] .* .!isnan.(df2[:,:g2])), :Ndefects])
# @show mean(df2[findall(x->x>0.65, df2[:,:g3] .* .!isnan.(df2[:,:g3])), :Ndefects])
# @show mean(df2[findall(x->x>0.65, df2[:,:g4] .* .!isnan.(df2[:,:g4])), :Ndefects])
# # @show mean(df2[findall(x->x>0.65, df2[:,:g5] .* .!isnan.(df2[:,:g5])), :Ndefects])
# @show mean(df2[findall(x->x>0.65, df2[:,:g6] .* .!isnan.(df2[:,:g6])), :Ndefects])
# @show mean(df2[findall(x->x>300, df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g2])
# @show mean(df2[findall(x->x>300, df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g3])
# @show mean(df2[findall(x->x>300, df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g4])
# @show mean(df2[findall(x->x>300, df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g5])
# @show mean(df2[findall(x->x>300, df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g6])

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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/GammapvsNdef.png", fig)
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/GammapvsNdef2.png", fig2)
    
    gamma_Ndef = DataFrame(x=Ndef)
    for p=1:5
        gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
        gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
    end
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","gamma_Ndef.csv"), gamma_Ndef)
    
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/GammapvsNdef3.png", fig3)
    
    g34_Ndef = DataFrame(x = Ndefects, y1=g3.-G_avg[2], y2=g4.-G_avg[3])
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","g34_Ndef.csv"), g34_Ndef)
return fig, fig2, fig3
end

function gamma_as_Ndefpm(stepp=30)
    df2 = CSV.read(dir_df*"DF2-pm.csv", DataFrame)
    Ndef=Array(25:stepp:400)
    Gammap = zeros(length(Ndef), 5, 2)
    for i=1:length(Ndef)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g2pm]
        Gammap[i,1,1] = mean(t)
        Gammap[i,1,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g3pm]
        Gammap[i,2,1] = mean(t)
        Gammap[i,2,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g4pm]
        Gammap[i,3,1] = mean(t)
        Gammap[i,3,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g5pm]
        Gammap[i,4,1] = mean(t)
        Gammap[i,4,2] = var(t)
        t = df2[findall(x-> (x>=(i-1)*stepp && x < i*stepp), df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])), :g6pm]
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/GammapvsNdef-pm.png", fig)
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/GammapvsNdef2-pm.png", fig2)
    
    gamma_Ndef = DataFrame(x=Ndef)
    for p=1:5
        gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
        gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
    end
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","gamma_Ndef-pm.csv"), gamma_Ndef)
    
    Ndefects = df2[:,:Ndefects] .* .!isnan.(df2[:,:Ndefects])
    g2 = df2[:,:g2pm] .* .!isnan.(df2[:,:Ndefects]) 
    g3 = df2[:,:g3pm] .* .!isnan.(df2[:,:Ndefects]) 
    g4 = df2[:,:g4pm] .* .!isnan.(df2[:,:Ndefects]) 
    g5 = df2[:,:g5pm] .* .!isnan.(df2[:,:Ndefects]) 
    g6 = df2[:,:g6pm] .* .!isnan.(df2[:,:Ndefects]) 
    fig3 = Figure(size = (1000, 800))
    ax3 = Makie.Axis(fig3[1,1], xlabel="Number of defects", ylabel=L"\langle|\gamma_p|\rangle-\gamma_P^0", xticklabelsize=24, yticklabelsize=24, xlabelsize=30, ylabelsize=30)
    xlims!(ax3, 5, 380)
    # scatter!(ax3, Ndefects, g2.-G_avg[1], color=:violet, makersize=8, label=L"p=2")
    scatter!(ax3, Ndefects, g4.-G_avg[3], color=:black, makersize=8, label=L"p=4")
    scatter!(ax3, Ndefects, g6.-G_avg[5], color=:blue, makersize=8, label=L"p=6")
    # scatter!(ax3, Ndefects, g5.-G_avg[4], color=:green, makersize=8, label=L"p=5")
    # scatter!(ax3, Ndefects, g6.-G_avg[5], color=:blue, makersize=8, label=L"p=6")
    axislegend(; position = :rt, labelsize = 24)#, markersize = 5);
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/GammapvsNdef3-pm.png", fig3)
    
    g34_Ndef = DataFrame(x = Ndefects, y1=g4.-G_avg[3], y2=g6.-G_avg[5])
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","g34_Ndef-pm.csv"), g34_Ndef)
return fig, fig2, fig3
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
    # isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    # save("Z:/Fig_paper/GammapvsNdef.png", fig)
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
    # isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    # save("Z:/Fig_paper/GammapvsNdef.png", fig)
return fig, fig2
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/As_zetarho_Ndefects.png", f1)
    
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save("Z:/Fig_paper/As_zetarho_Gamma34.png", f2)
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save(string("Z:/Fig_paper/As_rho0_Ndefects", name_file_addon,".png"), f1)
    
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save(string("Z:/Fig_paper/As_rho0_Gamma34",name_file_addon,".png"), f2)
    
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
    isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
    save(string("Z:/Fig_paper/As_rho0_Gamma34_mean",name_file_addon,".png"), f2b)
    
    return f1, f2, f2b
end

# make_full_heatmap_idx([0.2], [0.4,0.5,0.6,0.65,0.7,0.75,0.8,1.2], [4]; part = 1.0, Lz = 3, Lx = 3)
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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_$(tkd[1])_$(tr0[1])_$(tzr[1])_v2.png", fig)
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
    isdir("D:/Fig_paper/") ? nothing : mkpath("D:/Fig_paper/")
    save("D:/Fig_paper/phases_zoom_$(tkd[1])_$(tr0[1])_$(tzr[1])_v2.png", fig)
    return fig
end

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
    save("Z:/Fig_paper/dist_kd=$kk.png", f)
    
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
    save("Z:/Fig_paper/upper_kd=$kk.png", f)
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
    save("Z:/Fig_paper/dist_r0=$rr.png", f)
    
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
    save("Z:/Fig_paper/upper_r0=$rr.png", f)
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
# using Images, ImageView
# segments = felzenszwalb(1 .- bw, 10, 400)  # removes segments with fewer than 400 pixels
# imshow(map(i->segment_mean(segments,i), labels_map(segments)))
# labels = labels_map(segments)
# colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))

# bw = (rhoex .> 1.05*0.6) 
# # using Images, ImageView
# segments = felzenszwalb(1 .- bw, 10, 400)  # removes segments with fewer than 100 pixels
# imshow(map(i->segment_mean(segments,i), labels_map(segments)))
# labels = labels_map(segments)
# colored_labels = IndirectArray(labels, distinguishable_colors(maximum(labels)))

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
    save("Z:/Fig_paper/saturation.png", f)
    sat = DataFrame(x=(1:(Nt-1))*1000)
    for i=1:Nidx
        sat[!,"y$i"] = Nd[i,2:Nt]
    end
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","saturation.csv"), sat)
end

function figure_triple(idx, k=2)
    mkpath("Z:/Fig_paper/triple_$idx_$k/")
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
    save("Z:/Fig_paper/triple_$idx_$k/density.png", fig1)
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
    save("Z:/Fig_paper/triple_$idx_$k/voronoi.png", fig2)
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
    CSV.write(joinpath("Z:/Fig_paper/triple_$idx_$k/","gamma.csv"), gam)
    return fig1, fig2, gam
    # return fig1, fig2
end

function figure_triple_pm(idx, k)
    dir = "D:/Fig_paper/triple_pm_$(idx)_$k/"
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
    # save("D:/Fig_paper/triple_pm_$k/voronoi.png", fig2)
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
#     CSV.write(joinpath("Z:/Fig_paper/triple_pm_$k/","gamma.csv"), gam)
#     return fig1, fig2, gam
    return fig1#, fig2
end

function figure_triple_zoom(idx, k=2)
    mkpath("Z:/Fig_paper/triple_zoom_$k/")
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
    save("Z:/Fig_paper/triple_zoom_$k/density.png", fig1)
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
    # save("Z:/Fig_paper/triple_zoom_$k/voronoi.png", fig2)
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
#     CSV.write(joinpath("Z:/Fig_paper/triple_zoom_$k/","gamma.csv"), gam)
    # return fig1, fig2, gam
    return fig1#, fig2
end

function figure_triple_zoom_pm(idx, k)
    mkpath("Z:/Fig_paper/triple_zoom_pm_$k/")
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
    save("Z:/Fig_paper/triple_zoom_pm_$k/density.png", fig1)
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
    save("Z:/Fig_paper/triple_zoom_pm_$k/voronoi.png", fig2)
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
    CSV.write(joinpath("Z:/Fig_paper/triple_zoom_pm_$k/","gamma.csv"), gam)
    return fig1, fig2, gam
    # return fig1, fig2
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
    CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","saturation3.csv"), sat)
end

function profile(idx)
    dir_ = "Z:/2defects/"
    df_ = CSV.read(dir_*"DF.csv", DataFrame)
    
    path = string("$idx/")
    dird = dir_*path*"Data/"
    dirdp = dir_*path*"Defects_pos/"
    file_names = readdir(dird)
    Nt = length(file_names)
    marg = 25
    x = -marg:marg
    for t=1:Nt
        pts = load(string(dirdp,@sprintf("%04i",t),".jld"), "pos")[3:4,:]
        si = load(string(dirdp,@sprintf("%04i",t),".jld"), "charges")[3:4,:]
        zz = Int(pts[1,2]-200)
        x1 = Int(pts[1,1]-200)
        x2 = Int(pts[2,1]-200)
        den = load(string(dird,file_names[t]), "rho")
        Nx, Nz = size(den)
        dm_line = den[x2-marg:x2+marg, zz]
        dm_line .-= minimum(dm_line)
        dm_line ./= maximum(dm_line)
        dp_line = den[x1-marg:x1+marg, zz]
        dp_line .-= minimum(dp_line)
        dp_line ./= maximum(dp_line)
        pm_line = (load(string(dird,file_names[t]), "P")[x2-marg:x2+marg, zz,1]).^2
        pp_line = (load(string(dird,file_names[t]), "P")[x1-marg:x1+marg, zz,1]).^2
        f,a,l = lines(x,dm_line)
        lines!(a,x,dp_line, linestyle = :dash)
        lines!(a,x,pm_line)
        lines!(a,x,pp_line, linestyle = :dash)
        f
    end
end

function maxrho(idx)
    path = string("$idx/")
       
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []
    
    num_files = length(files_names)
    
    rho = load(string(dird,files_names[end]), "rho")
    println(minimum(rho), ", ", maximum(rho))
end

function plot_2defects_zoom(idx; dt = 1, sc=1)
    # set_theme!(theme_black())
    # IJulia.clear_output(true)
    dir_df = "Z:/2defects/"
    dir_fig = "D:/2defects_white/"
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

#make_full_heatmap_idx([0.2], [0.4,0.5,0.6,0.65,0.7,0.75,0.8,1.2], [4]; part = 1.0, Lz = 3, Lx = 3)
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

# function plot_ΔP(idx)
#     path = string("$idx/")
#     dird = dir_df*path*"Data/"
#     file_name = isdir(dird) ? readdir(dird)[end] : []
#     P_cpu = load(string(dird,file_name), "P")
#     # rho = load(string(dird,file_name), "rho")
#     # Px = P_cpu[:,:,1]
#     # Pz = P_cpu[:,:,2]
#     ex = 200
#     # tex = round(Int, ex*2)
#     Px = expend_periodic_array2(P_cpu[:,:,1], ex)
#     Pz = expend_periodic_array2(P_cpu[:,:,2], ex)
#     @show Nx, Nz = size(Px)
#     DP = zeros(Nx, Nz)

#     Δx = 0.01;   Δz = 0.01
#     Δx2 = Δx*Δx; Δz2 = Δz*Δz

#     for i=1:Nx, j=1:Nz
#         i_ = mod(i-1,1:Nx); ip = mod(i+1,1:Nx)
#         jm = mod(j-1,1:Nz); jp = mod(j+1,1:Nz)
#         ΔPx = (Px[i_,j] - 2*Px[i,j] + Px[ip,j])/Δx2 + (Px[i,jm] - 2*Px[i,j] + Px[i,jp])/Δz2 # ΔPx
#         ΔPz = (Pz[i_,j] - 2*Pz[i,j] + Pz[ip,j])/Δx2 + (Pz[i,jm] - 2*Pz[i,j] + Pz[i,jp])/Δz2 # ΔPz
#         DP[i,j] = ΔPx + ΔPz
#     end

#     # heatmap(DP[3:end-3, 3:end-3])
#     # @show maximum(DP)
#     # @show minimum(DP)
#     # @show findmin(DP)
#     fig1 = Figure(size = (1000, 1000))
#     ax = Makie.Axis(fig1[1,1], aspect = DataAspect())
#     heatmap!(ax, DP[3:end-3, 3:end-3])

#     # elem_1 = MarkerElement(color = :blue, marker = :circle, markersize = 25)
#     # elem_2 = MarkerElement(color = :red, marker = :circle, markersize = 25)
#     # elem_1 = MarkerElement(color = :black, marker = :circle, markersize = 25)
#     # elem_2 = MarkerElement(color = :black, marker = :circle, markersize = 25)
#     # Legend(fig2[1, 1], [elem_1, elem_2], ["                    ", "                    "], framevisible = false, orientation = :horizontal)
#     hideydecorations!(ax)
#     hidexdecorations!(ax)
#     # colgap!(fig2.layout, 0)
#     # rowgap!(fig2.layout, -5)
#     tightlimits!(ax)
#     save("Z:/Temp/DP.png", fig1)
#     using Image, ImageView, ImageSegmentation, Random 
#     img = load("Z:/Temp/DP.png")
#     imgg = green.(img) .* (green.(img) .> 0.3) .* (green.(img) .< 0.6 )
#     imshow(imgg)
#     segments = felzenszwalb(Gray.(imgg), 10, 800)
#     imshow(map(i->get_random_color(i), labels_map(segments)))
# end
# function get_random_color(seed)
#     Random.seed!(seed)
#     rand(RGB{N0f8})
# end


function convert_gif_to_mp4(fn)
    str = string("""-i $fn.gif -f lavfi -i anullsrc -vf "scale='trunc(in_w/2)*2':'trunc(in_h/2)*2',format=yuv420p,fps=8.33" -movflags +faststart -shortest $fn.mp4""")
    # run(`ffmpeg $str`)
    println("ffmpeg $str ")
end

function plot_dens_angle(idx)
    dir = "D:/Fig_paper/density_angle_$(idx)/"
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
    dir = "D:/Fig_paper/density_and_angle_tri_$(idx)/"
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

# A(x, z) = z>0 ? 3*4/3*z*x^2 - 5/2*0.1 : 3*4/3*x^2 - 5/2*0.1
# f(x, z, t) = A(x,z)/(3*x) + (sqrt(A(x,z)*1e-4/(3*x^2))-sqrt(2/(3*x^3*t)))^2
# B(q,x,z) = q^2*(x^2*(A(x,z)-3*x*z)/(2*q^2+1)+1e-4*A(x,z))
# C(q, x, z, t) = (0.1-1e-4*q^2)*x^2*q^2/(2*q^2+1)-x*(0.1+1e-4*q^2)*1e-4*q^2
# tq = 1 ./ Array(0.01:0.01:10);
# trho = Array(0.3:0.01:1.2);
# ttau = [1, 5];
# tzr = Array(1:0.01:30)
# Lx = length(trho)
# Ly = length(ttau)
# Lz = length(tzr)
# res = zeros(Lx,Ly);
# for i = 1:Lx, j = 1:Ly
#     res[i,j] = f(trho[i], 4, ttau[j])
# end
# ff, a, l = lines(trho, res[:,1]);
# for i in 2:Ly
#     lines!(a, trho, res[:,i])
# end
# ff

# tzr = vcat(0,Array(1:0.01:99));
# Lz = length(tzr);
# Lq = length(tq);
# resb = 100*ones(Lx, Ly);
# for i=1:Lx, j=1:Ly
#     @show rho0 = trho[i]
#     @show tau = ttau[j]
#     # global zrc = 100
#     # zrcq = 100*ones(Lq)
#     ρ = 1e-3*ones(Lq)
#     p = 1e-3*ones(Lq)
#     ρt = zeros(Lq)
#     pt = zeros(Lq)
#     for k = 1:Lz
#         zr = tzr[k]
#         # global zrc = B(1/10,trho[i],tzr[k],ttau[j]) < -1 ? minimum([zrc, tzr[k]]) : zrc
#         # for q in 1:Lq
#         for t=1:1000
#             ρt .= ρ .+ 0.01*( -(1/tau .+ B.(tq,rho0,zr)) .* ρ)# .+ C.(tq,rho0,zr,tau) .*p )
#             # pt .= p .+ 0.01*( 0.1*rho0*ρ .- rho0^2*(2*0.1 .+ tq.*tq*1e-4).*p )
#             copyto!(ρ, ρt)
#             # copyto!(p, pt)
#             # zrcq[q] = B(tq[q],trho[i],tzr[k],ttau[j]) < -1 ? minimum([zrc, tzr[k]]) : zrc
#         end
#         if any(ρ .> 1e-3)
#             @show resb[i,j] = zr
#             break
#         end
#         # end
#         # zrc = minimum([zrc, minimum(zrcq)])
#     end
#     # resb[i,j] = zrc
# end
# ff2, a2, l2 = lines(trho, resb[:,1]);
# for i in 2:Ly
#     lines!(a2, trho, resb[:,i])
# end
# ylims!(a2,(0,5))
# ff2
# dfzrc = DataFrame(x=trho, y=res[:,1], z=res[:,2], yy=resb[:,1], zz=resb[:,2]);
# mkpath("D:/Fig_paper/Data_Fig_Tikz/");
# CSV.write(joinpath("D:/Fig_paper/Data_Fig_Tikz/","zrc.csv"), dfzrc)

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


function gamma_as_Ndefpm(stepp=30)
    dir_df = "E:/AdptDt/"
# dir_fig = "E:/AdptDt/"
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
    CSV.write(joinpath("D:/Fig_paper/Data_Fig_Tikz/","g34_rho-pm.csv"), g34_rho)

    f, a, l = scatter(ru, res[:,1] .- G_avg[3]);
    scatter!(a, ru, res[:,2] .- G_avg[5])
    f

#     gamma_Ndef = DataFrame(x=Ndef)
#     for p=1:5
#         gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
#         gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
#     end
#     CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","gamma_Ndef-pm.csv"), gamma_Ndef)
    
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
#     isdir("Z:/Fig_paper/") ? nothing : mkpath("Z:/Fig_paper/")
#     save("Z:/Fig_paper/GammapvsNdef3-pm.png", fig3)
    
#     g34_Ndef = DataFrame(x = Ndefects, y1=g4.-G_avg[3], y2=g6.-G_avg[5])
#     CSV.write(joinpath("Z:/Fig_paper/Data_Fig_Tikz/","g34_Ndef-pm.csv"), g34_Ndef)
# return fig, fig2, fig3
end



function maxV(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    # num_files = length(files_names)
    v_cpu = load(string(dird,files_names[end]), "v")
    @show max_v = sqrt(maximum(v_cpu[:,:,1].^2 .+ v_cpu[:,:,2].^2))
    return max_v
end