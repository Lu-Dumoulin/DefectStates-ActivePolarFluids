# =============================================================================
#  Analysis and plotting for the "gamma rho" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_gamma_rho.jl
# =============================================================================

include("MakePlots.jl")

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

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    gamma_as_Ndef()
    gamma_as_Ndefpm()
    gamma_as_ddef()
    make_heatmap()
end
