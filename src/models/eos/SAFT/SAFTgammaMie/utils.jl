function gc_to_comp_sites(sites::SiteParam,groups::GroupParam)
    #given some groups and some sites calculated over those groups
    #calculates "flattened" sites
    comps = groups.components
    sitenames = deepcopy(sites.sites)
    gcnames = groups.flattenedgroups
    for i in 1:length(sitenames)
        gci = gcnames[i]
        sitei = sitenames[i]
        for a in 1:length(sitei)
            sitei[a] = gci * '{' * sitei[a] * '}'
        end
    end

    flattened_comp_sitenames = collect(Iterators.flatten(sitenames))
    flattened_comp_nsites = collect(Iterators.flatten(sites.n_sites))
    siteidx = deepcopy(sites.i_sites)
    for i in 1:length(siteidx)
        sitei = siteidx[i]
        for a in 1:length(sitei)
            sitei[a] = i
        end
    end 
    flattened_gcidx = collect(Iterators.flatten(siteidx))
    flattened_comp_isites = collect(Iterators.flatten(sites.i_sites))
    idxdict = Dict((flattened_gcidx[i],flattened_comp_isites[i]) => i for i in 1:length(flattened_gcidx) )


    n_flattenedsites = Vector{Vector{Int}}(undef,length(comps))
    flat_groups = groups.n_flattenedgroups
    for i in 1:length(n_flattenedsites) 
        n_flattenedsites[i] = flat_groups[i][flattened_gcidx] .* flattened_comp_nsites
    end
    k = length(flattened_comp_sitenames)
    pairs = [(comps[i],[flattened_comp_sitenames[j] =>n_flattenedsites[i][j] for j in 1:k]) for i in 1:length(comps)]
    n_flattenedsites

    for pair in pairs
        pairname,pairvec = pair
        filter!(x->!iszero(last(x)),pairvec)
    end


    return SiteParam(pairs),idxdict 
end


#returns a compressed assoc matrix corresponding to the indices
#of the old gc values, arranged by component instead.
#that is, 
function gc_to_comp_assoc_idx(param::AssocParam,sites::SiteParam,idxdict)
    pvals = param.values
    vals,outer,inner = pvals.values, pvals.outer_indices,pvals.inner_indices
    ngc = length(vals)
    nsites = length(sites.flattenedsites)
    comps = sites.components
    ncomps = length(comps)
    site1 = Vector{Int64}(undef,ngc)
    site2 = Vector{Int64}(undef,ngc)
    assoc_idxdict = Dict{Tuple{Int,Int},Int}()

    for ii in 1:ngc
        i,j = outer[ii]
        a,b = inner[ii]
        s1 =  idxdict[(i,a)]
        s2 = idxdict[(j,b)]
        site1[ii] = s1
        site2[ii] = s2
        assoc_idxdict[(s1,s2)] = ii
        assoc_idxdict[(s2,s1)] = ii
    end

    #return site1,site2
    x = Matrix{Matrix{Int}}(undef,ncomps,ncomps)
    for i = 1:ncomps
        xi = zeros(Int,nsites,nsites)
        for (a,b) in Iterators.product(sites.i_sites[i],sites.i_sites[i])
            abval = get(assoc_idxdict,(a,b),0)
            
            xi[a,b] = abval
            xi[b,a] = abval
            
        end
        x[i,i] = xi
        for j = 1:i-1
            xi = zeros(Int,nsites,nsites)
            for (a,b) in Iterators.product(sites.i_sites[i],sites.i_sites[j])
                abval = get(assoc_idxdict,(a,b),0)
                
                xi[a,b] = abval
                xi[b,a] = abval
            end
            x[i,j] = xi
            x[j,i] = xi
        end
    end
    val =  Compressed4DMatrix(x)
    new_inneridx = val.inner_indices
    for (idx,(i,j),(a,b)) ∈ indices(val)
        _a = sites.i_flattenedsites[i][a]
        _b = sites.i_flattenedsites[j][b]
        new_inneridx[idx] = (_a,_b)
    end

    n = findall(!=((0,0)),new_inneridx)
    values = val.values[n]
    outer_indices = val.outer_indices[n]
    inner_indices = new_inneridx[n]
    outer_size = val.outer_size
    inner_size = val.inner_size
    newvals = Compressed4DMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
end
    
  
