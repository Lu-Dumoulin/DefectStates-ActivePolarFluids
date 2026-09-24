# =============================================================================
#  Figure 4 - fig:Ndefects
#
#  Defect density in the (rho0, zeta_rho) plane at six renewal times.
#
#      DATA_DIR=/path/to/fig4-runs/ FIG_DIR=.../figures/Fig4/ \
#          julia --project=. Analysis/make_fig4.jl
#
#  Inputs : params/DF_4.csv  (12 rho0 x 14 zeta_rho x 6 tau = 1008 runs)
#  Outputs: Ndef.csv, the heatmap data the figure reads
#
#  make_heatmap's tkd is the six renewal rates, i.e. tau = 0.1, 0.2, 1, 5, 10
#  and 100; 12 x 14 x 6 = 1008 is exactly the row count of figures/Fig4/Ndef.csv.
#
#  Not reproduced here: zrc_plus.csv (and the superseded df_zrc.csv), the
#  critical-activity curves drawn as white dots. They come from the linear
#  stability analysis - the same quantity as figure 5(a) - solving
#  Eq. (zetaRhoDeltaMuC) self-consistently with a = 4 zeta_c/3. See
#  make_fig5.jl and LSA.jl; the exact routine is not in this repository, and
#  figures/README.md records what is known about it.
# =============================================================================

include("MakePlots.jl")
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
    
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "Ndefect.png"), fig)
    
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","Ndefb.csv"), ndef)
    
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
    
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "Gamma4.png"), fig2)
    
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
    
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "Gamma3.png"), fig3)
#     return fig, fig2, fig3
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, "As_zetarho_Ndefects.png"), f1)
    
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, "As_zetarho_Gamma34.png"), f2)
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, string("As_rho0_Ndefects", name_file_addon,".png"), f1)
    
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, string("As_rho0_Gamma34",name_file_addon,".png"), f2)
    
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, string("As_rho0_Gamma34_mean",name_file_addon,".png"), f2b)
    
    return f1, f2, f2b
end


# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(joinpath(dir_fig, "Data_Fig_Tikz"))
    make_heatmap([10, 5, 1.0, 0.2, 0.1, 0.01])
end
