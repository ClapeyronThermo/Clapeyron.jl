function gc_to_comp_sites(sites::SiteParam,groups::GroupParam)
    #given some groups and some sites calculated over those groups
    #calculates "flattened" sites

    comps = groups.components #component names

    #each GC/GCsite pair encodes a tuple in GC space
    #we use this traslator later to convert any assoc param into it's component equivalent
    site_translator = Vector{Vector{Tuple{Int,Int}}}(undef,length(comps))

    #shortcut for non-assoc case
    if length(sites.n_sites.v) == 0
        new_sites = SiteParam(groups.components)
        return new_sites,site_translator
    end

    sitenames = deepcopy(sites.sites) #group sites
    gcnames = groups.flattenedgroups #group names
    gc_groups = groups.groups
    #with this, each site now has an unique name
    for i in 1:length(sitenames)
        gci = gcnames[i]
        sitei = sitenames[i]
        for a in 1:length(sitei)
            sitei[a] = gci * '/' * sitei[a]
        end
    end

    #now we fill our new component sites,
    gc_n_sites = sites.n_sites.v #should be of the same size as flattened_comp_sitenames
    comp_n_sites = Vector{Vector{Int}}(undef,length(comps))
    comp_sites = Vector{Vector{String}}(undef,length(comps))


    for i in 1:length(comps)
        n_sites_i = Int[]
        comp_n_sites[i] = n_sites_i

        sites_i = String[]
        comp_sites[i] = sites_i
        site_translator_i = NTuple{2,Int}[]
        site_translator[i] = site_translator_i
        gc_name_i = gc_groups[i]
        for k in 1:length(gc_name_i)
            for (w,comp_gcname) in enumerate(Iterators.flatten(sitenames))
                gcname_ik = gc_name_i[k]
                lookup_cgname = gcname_ik * '/'
                if startswith(comp_gcname,lookup_cgname)
                    #fill sites, n_sites
                    push!(sites_i,comp_gcname)
                    push!(n_sites_i,gc_n_sites[w])

                    #fill translation between gc_gcsite combination and original indices for assoc
                    gc_ik = findfirst(isequal(gcname_ik),groups.flattenedgroups)
                    push!(site_translator_i,(gc_ik,length(n_sites_i)))
                end
            end
        end
    end
    new_sites = SiteParam(comps,comp_sites,comp_n_sites,sites.sourcecsvs)

    return new_sites,site_translator
end


function gc_to_comp_sites(param::AssocParam,sites::SiteParam,site_translator)

    #shortcut for non-assoc case
    if length(sites.n_sites.v) == 0
        new_val = Compressed4DMatrix{eltype(param)}()
        return AssocParam(param.name,sites.components,new_val,sites.sites,param.sourcecsvs,param.sources)
    end
    new_val = assoc_similar(sites,eltype(param))
    for i in 1:length(sites.components)
        site_translator_i = site_translator[i]
        for j in 1:i
            ij_pair = new_val[i,j]
            #display(TextDisplay(stdout),MIME"text/plain"(),ij_pair)
            site_translator_j = site_translator[j]
            aa,bb = length(site_translator_i),length(site_translator_j)
            for a in 1:length(site_translator_i)
                i_gc,a_gc = site_translator_i[a]
                for b in 1:length(site_translator_j)
                    j_gc,b_gc = site_translator_j[b]
                    #absolute index, relative to the Compressed4DMatrix
                    idx = validindex(ij_pair,a,b)
                    if idx != 0 #if the index is valid
                        ijab_val = param[i_gc,j_gc][a_gc,b_gc]
                        if !_iszero(ijab_val) #if the value is not zero
                            ij_pair.vec.values[idx] =ijab_val
                        end
                    end
                end
            end
        end
    end
    dropzeros!(new_val) #clean all zero values
    return AssocParam(param.name,sites.components,new_val,sites.sites,param.sourcecsvs,param.sources)
end




