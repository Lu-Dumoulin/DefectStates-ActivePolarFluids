# =============================================================================
#  Analysis and plotting for the "ndef" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_ndef.jl
# =============================================================================

include("MakePlots.jl")

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

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    Ndefects_as_act()
    Ndefects_as_rho0()
end
