# =============================================================================
#  Analysis and plotting for the "2defects" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_2defects.jl
# =============================================================================

include("MakePlots.jl")

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

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    dist_2defects()

    # These need arguments (the simulation index, and which panel):
    # heatmap_dist_2defects(kk)
    # heatmap_dist_2defects_b(rr)
    # velocity_defects(idx)
    # plot_2defects_zoom(idx)
end
