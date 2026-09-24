# =============================================================================
#  Shared analysis machinery: defect detection and characterisation, image
#  segmentation, structure factors, correlations and batch passes over a run
#  set. No figure is assembled here.
# =============================================================================

include("config.jl")

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

function profile(idx)
    dir_ = dir_df                      # was pinned to "Z:/2defects/"
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

function maxV(idx)
    path = string("$idx/")
    dird = dir_df*path*"Data/"
    files_names = isdir(dird) ? readdir(dird) : []

    # num_files = length(files_names)
    v_cpu = load(string(dird,files_names[end]), "v")
    @show max_v = sqrt(maximum(v_cpu[:,:,1].^2 .+ v_cpu[:,:,2].^2))
    return max_v
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

# helpers defined mid-file in the original
kdoar(ρ0, act, r) = 6/r^2*(ρ0^3-0.75*act*ρ0^2)/(3*0.1*ρ0^2/1e-4-1)
ret1(x,y) = 1.0
