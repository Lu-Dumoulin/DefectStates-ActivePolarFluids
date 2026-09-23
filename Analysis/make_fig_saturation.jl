# =============================================================================
#  Analysis and plotting for the "saturation" figure.
#
#  Run against the run set that figure uses:
#      DATA_DIR=/path/to/runset/ julia --project=. Analysis/make_fig_saturation.jl
# =============================================================================

include("MakePlots.jl")

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

# --- entry point -------------------------------------------------------------
# Runs when this file is executed directly. The calls below use the default
# arguments the functions were written with; check them against the run set in
# DATA_DIR before trusting the output.
if abspath(PROGRAM_FILE) == @__FILE__
    plot_saturation()
    plot_saturation2()
end
