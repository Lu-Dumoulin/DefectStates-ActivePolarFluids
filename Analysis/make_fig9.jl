# =============================================================================
#  Figure 9 - fig:gamma_rho
#
#  Shape function against target density and against defect density.
#
#      DATA_DIR=/path/to/fig9-runs/ FIG_DIR=<repo>/figures/Fig9/ \
#          julia --project=. Analysis/make_fig9.jl
#
#  Inputs : params/DF_9.csv (1680 runs) analysed into DF2-pm.csv
#  Outputs: g34_rho-pm.csv (panel a) and g34_densdef.csv (panel b)
#
#  Panel (b) comes from the commented-out block at the end of gamma_as_Ndefpm,
#  which builds x = Ndefects. The committed g34_densdef.csv is that same data
#  with x divided by the area of the padded domain, (1008*0.01)^2 = 101.6064 -
#  the ratio between the two committed files is exactly that, on every row.
#  The live code keeps only Ndefects > 150, i.e. defect density > 1.5, which is
#  the threshold the caption states.
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, "GammapvsNdef.png"), fig)
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, "GammapvsNdef2.png"), fig2)
    
    gamma_Ndef = DataFrame(x=Ndef)
    for p=1:5
        gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
        gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
    end
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","gamma_Ndef.csv"), gamma_Ndef)
    
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
    mkpath(dir_fig)
    save(joinpath(dir_fig, "GammapvsNdef3.png"), fig3)
    
    g34_Ndef = DataFrame(x = Ndefects, y1=g3.-G_avg[2], y2=g4.-G_avg[3])
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","g34_Ndef.csv"), g34_Ndef)
return fig, fig2, fig3
end

function gamma_as_Ndefpm(stepp=30)
    # Reads the run set configured through DATA_DIR. This used to be pinned to
    # "E:/AdptDt/", which silently ignored it.
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
    CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","g34_rho-pm.csv"), g34_rho)

    f, a, l = scatter(ru, res[:,1] .- G_avg[3]);
    scatter!(a, ru, res[:,2] .- G_avg[5])
    f

#     gamma_Ndef = DataFrame(x=Ndef)
#     for p=1:5
#         gamma_Ndef[!,"y$p"] = Gammap[:,p,1].-G_avg[p]
#         gamma_Ndef[!,"error$p"] = Gammap[:,p,2]
#     end
#     CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","gamma_Ndef-pm.csv"), gamma_Ndef)
    
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
#     mkpath(dir_fig)
#     save(joinpath(dir_fig, "GammapvsNdef3-pm.png"), fig3)
    
#     g34_Ndef = DataFrame(x = Ndefects, y1=g4.-G_avg[3], y2=g6.-G_avg[5])
#     CSV.write(joinpath(dir_fig, "Data_Fig_Tikz","g34_Ndef-pm.csv"), g34_Ndef)
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
    # mkpath(dir_fig)
    # save(joinpath(dir_fig, "GammapvsNdef.png"), fig)
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
    # mkpath(dir_fig)
    # save(joinpath(dir_fig, "GammapvsNdef.png"), fig)
return fig, fig2
end

# --- entry point -------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(joinpath(dir_fig, "Data_Fig_Tikz"))
    gamma_as_Ndefpm()
end
