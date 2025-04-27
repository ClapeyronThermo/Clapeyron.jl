"""
    SiteParam
Struct holding site parameters.
Is built by parsing all association parameters in the input CSV files.
It has the following fields:
* `components`: a list of all components (or groups in Group Contribution models)
* `sites`: a list containing a list of all sites corresponding to each component (or group) in the components field
* `n_sites`: a list of the site multiplicities corresponding to each site in `flattenedsites`
* `flattenedsites`: a list of all unique sites
* `i_sites`: an iterator that goes through the indices corresponding  to each site in `flattenedsites`
* `n_flattenedsites`: the site multiplicities corresponding to each site in `flattenedsites`
* `i_flattenedsites`: an iterator that goes through the indices for each flattened site
Let's explore the sites in a 3-component `SAFTGammaMie` model:
```julia
julia> model3 = SAFTgammaMie([    
                "ethanol",
                ("nonadecanol", ["CH3"=>1, "CH2"=>18, "OH"=>1]),     
                ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])
                               ])
SAFTgammaMie{BasicIdeal} with 3 components:
 "ethanol"
 "nonadecanol"
 "ibuprofen"
Contains parameters: segment, shapefactor, lambda_a, lambda_r, sigma, epsilon, epsilon_assoc, bondvol 
julia> model3.sites
SiteParam with 8 sites:
 "CH2OH": "H" => 1, "e1" => 2     
 "CH3": (no sites)
 "CH2": (no sites)
 "OH": "H" => 1, "e1" => 2        
 "COOH": "e2" => 2, "H" => 1, "e1" => 2
 "aCCH": (no sites)
 "aCCH2": (no sites)
 "aCH": (no sites)
julia> model3.sites.flattenedsites
3-element Vector{String}:
 "H"
 "e1"
 "e2"
julia> model3.sites.i_sites       
8-element Vector{Vector{Int64}}:
 [1, 2]
 []
 []
 [1, 2]
 [1, 2, 3]
 []
 []
 []
julia> model3.sites.n_sites       
8-element Vector{Vector{Int64}}:
 [1, 2]
 []
 []
 [1, 2]
 [2, 1, 2]
 []
 []
 []
```
"""
struct SiteParam <: ClapeyronParam
    components::Array{String,1}
    sites::Array{Array{String,1},1}
    n_sites::PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int}, Vector{Int}, SubArray{Int, 1, Vector{Int}, Tuple{UnitRange{Int64}}, true}}
    i_sites::Array{Array{Int,1},1}
    flattenedsites::Array{String,1}
    n_flattenedsites::Array{Array{Int,1},1}
    i_flattenedsites::Array{Array{Int,1},1}
    sourcecsvs::Array{String,1}
    site_translator::Union{Nothing,Vector{Vector{NTuple{2,Int}}}}
end

function SiteParam(components::Vector{String},sites::Array{Array{String,1},1},n_sites::Vector{Vector{Int}},sourcecsvs = String[],site_translator = nothing)
    n_sites = PackedVectorsOfVectors.pack(n_sites)
    param = SiteParam(components, 
    sites, 
    n_sites,
    Vector{Vector{Int}}(undef,0), 
    String[],
    Vector{Vector{Int}}(undef,0),
    Vector{Vector{Int}}(undef,0),
    sourcecsvs,
    site_translator)
        
    recombine!(param)
    return param
end

function Base.show(io::IO, mime::MIME"text/plain", param::SiteParam)
    print(io,"SiteParam ")
    len = length(param.components)
    println(io,"with ", len, " component", ifelse(len==1, ":", "s:"))
    site_print_i(ii,val) = __show_group_i(ii,val,"(no sites)")
    show_pairs(io,param.components,zip(param.sites,param.n_sites),": ",site_print_i)
end

function Base.show(io::IO, param::SiteParam)
    print(io,"SiteParam[")
    function wrap_print(io,val)
        print(io,'[')
        __show_group_i(io,val)
        print(io,']')
    end
    show_pairs(io,param.components,zip(param.sites,param.n_sites)," => ",wrap_print,pair_separator = ", ")
    print(io,"]")
end


function SiteParam(pairs::Dict{String,SingleParam{Int}},allcomponentsites)
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    sites = allcomponentsites
    ncomps = length(components)
    sourcecsvs = String[]
    for x in values(pairs)
        append!(sourcecsvs,x.sourcecsvs)  
    end
    unique!(sourcecsvs)

    n_sites = [[pairs[sites[i][j]].values[i] for j ∈ 1:length(sites[i])] for i ∈ 1:ncomps]  # or groupsites
    return SiteParam(components,sites,n_sites,sourcecsvs)
end

function SiteParam(input::PARSED_GROUP_VECTOR_TYPE,sourcecsvs::Vector{String}=String[])
    components = [first(i) for i ∈ input]
    raw_sites = [last(i) for i ∈ input]
    sites = [first.(sitepairs) for sitepairs ∈ raw_sites]
    n_sites = [last.(sitepairs) for sitepairs ∈ raw_sites]
    return SiteParam(components,sites,n_sites,sourcecsvs)
end

function SiteParam(components::Vector{String})
    n = length(components)
    return SiteParam(
    components,
    [String[] for _ ∈ 1:n],
    PackedVectorsOfVectors.pack([Int[] for _ ∈ 1:n]),
    [Int[] for _ ∈ 1:n],
    String[],
    [Int[] for _ ∈ 1:n],
    [Int[] for _ ∈ 1:n],
    String[],
    nothing)
end

function recombine!(param::SiteParam)
    components = param.components
    sites = param.sites
    n_sites = param.n_sites
    ℂ = length(components)

    #initialization of flattenedsites
    #flattenedsites = unique!(reduce(vcat,sites))
    
    flattenedsites = param.flattenedsites
    resize!(flattenedsites,0)
    for site ∈ sites
        append!(flattenedsites,site)
    end
    unique!(flattenedsites)
    #initialization of i_sites
    #i_sites = [[findfirst(isequal(site), flattenedsites) for site ∈ componentsites] for componentsites ∈ sites]

    i_sites = param.i_sites
    resize!(i_sites,ℂ)
    for i in 1:ℂ
        site = sites[i]
        if !isassigned(i_sites,i)
            i_sites[i] = Int[]
        end
        i_site = i_sites[i]
        resize!(i_site,length(site))
        for j in eachindex(i_site)
            i_site[j] = findfirst(isequal(site[j]), flattenedsites)::Int
        end
    end

    #initialization of n_flattenedsites, i_flattenedsites

    flat𝕊 = length(flattenedsites)
    n_flattenedsites = param.n_flattenedsites
    i_flattenedsites = param.i_flattenedsites
    
    resize!(n_flattenedsites,ℂ)
    resize!(i_flattenedsites,ℂ)
    
    for i in 1:ℂ
        if !isassigned(n_flattenedsites,i)
            n_flattenedsites[i] = Int[]
        end

        if !isassigned(i_flattenedsites,i)
            i_flattenedsites[i] = Int[]
        end
        n_flatsite = n_flattenedsites[i]
        i_flatsite = i_flattenedsites[i]

        resize!(n_flatsite,flat𝕊)
        resize!(i_flatsite,flat𝕊)

        n_flatsite .= 0
        i_flatsite .= 0

        setindex!(n_flatsite,n_sites[i],i_sites[i])
        setindex!(i_flatsite,1:length(i_sites[i]),i_sites[i])
    end
    #TODO: do something about site_translator
    #but if anyone tries to change the components of a SiteParam in a GC-comp mixed context, they already are
    #generating site_translator again.
    return param
end
"""
    assoc_similar(param::SiteParam)
    assoc_similar(param::SiteParam,::Type{𝕋}) where 𝕋 <:Number)

returns a `Clapeyron.Compressed4DMatrix` with the smae number of components as the input `AssocParam`, with the same element type as `𝕋`.
All site combinations are filled.
"""
function assoc_similar(param::SiteParam,::Type{𝕋}) where 𝕋 <:Number
    comps = 1:length(param.components)
    indices = NTuple{4,Int}[]
    values = 𝕋[]
    _0 = zero(𝕋)
    __set_idx_4d!(param.sites,values,indices)
    Compressed4DMatrix(values,indices)
end

function assoc_similar!(x::Compressed4DMatrix,param::SiteParam)
    resize!(x.values,0)
    _4dindices = NTuple{4,Int}[]
    __set_idx_4d!(param.sites,values,_4dindices)
    resize!(x.outer_indices,length(values))
    resize!(x.inner_indices,length(values))
    for i in 1:length(values)
        i,j,a,b = _4dindices[i]
        x.outer_indices[i] = (i,j)
        x.inner_indices[i] = (a,b)
    end
    return Compressed4DMatrix(values,x.outer_indices,x.inner_indices,x.outer_size,x.inner_size)
end

function Compressed4DMatrix(param::SiteParam)
    return assoc_similar(param,Float64)
end

#build dense assocparam from sites.
function AssocParam{T}(name,sites::SiteParam) where T <: Number
    values = assoc_similar(sites,T)
    return AssocParam(name,sites.components,values,sites.sites)
end

function gc_to_comp_sites(sites::SiteParam,groups::GroupParameter)
    #given some groups and some sites calculated over those groups
    #calculates "flattened" sites

    comps = groups.components #component names

    #each GC/GCsite pair encodes a tuple in GC space
    #we use this traslator later to convert any assoc param into it's component equivalent
    site_translator = Vector{Vector{Tuple{Int,Int}}}(undef,length(comps))

    #shortcut for non-assoc case
    if length(sites.n_sites.v) == 0
        new_sites = SiteParam(groups.components)
        return new_sites
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
    gc_n_groups = groups.n_groups
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
            #=
            each gc group has it's own amount of sites.
            to correctly translate those, we need an accumulator for each group
            
            =#
            gc_site_count = 0
            for (w,comp_gcname) in enumerate(Iterators.flatten(sitenames))
                gcname_ik = gc_name_i[k]
                lookup_cgname = gcname_ik * '/'
                if startswith(comp_gcname,lookup_cgname)
                    #fill sites, n_sites
                    push!(sites_i,comp_gcname)
                    push!(n_sites_i,gc_n_sites[w]*gc_n_groups[i][k])
                    #fill translation between gc_gcsite combination and original indices for assoc
                    #we add one to the length of gc_site_count. this accounts for adding the next site.
                    gc_site_count += 1
                    gc_ik = findfirst(isequal(gcname_ik),groups.flattenedgroups)::Int
                    push!(site_translator_i,(gc_ik,gc_site_count))
                end
            end
        end
    end

    new_sites = SiteParam(comps,comp_sites,comp_n_sites,sites.sourcecsvs,site_translator)

    return new_sites
end

#=
Utilities to create "group-component" sites

if we re-index the association sites that are present in a GC-based approach, the result is the same.
=#
function gc_to_comp_sites(param::AssocParam,sites::SiteParam)

    #shortcut for non-assoc case
    if length(sites.n_sites.v) == 0
        new_val = Compressed4DMatrix{eltype(param)}()
        return AssocParam(param.name,sites.components,new_val,sites.sites,param.sourcecsvs,param.sources)
    end
    site_translator = sites.site_translator
    new_val = assoc_similar(sites,eltype(param))
    for i in 1:length(sites.components)
        site_translator_i = site_translator[i]
        for j in 1:length(sites.components)
            ij_pair = new_val[i,j]
            #display(TextDisplay(stdout),MIME"text/plain"(),ij_pair)
            site_translator_j = site_translator[j]
            aa,bb = length(site_translator_i),length(site_translator_j)
            for a in 1:length(site_translator_i)
                i_gc,a_gc = site_translator_i[a]
                for b in 1:length(site_translator_j)
                    j_gc,b_gc = site_translator_j[b]
                    #absolute index, relative to the Compressed4DMatrix
                    idx = validindex(ij_pair,a,b,false)
                    if idx != 0
                        ijab1 = param[i_gc,j_gc]
                        idx_ijab = validindex(ijab1,a_gc,b_gc,false)
                        if idx_ijab != 0
                            ijab_val = ijab1.values.values[idx_ijab]
                            if !_iszero(ijab_val) #if the value is not zero
                                ij_pair.values.values[idx] =ijab_val
                            end
                        end
                    end
                end
            end
        end
    end
    dropzeros!(new_val) #clean all zero values
    return AssocParam(param.name,sites.components,new_val,sites.sites,param.sourcecsvs,param.sources)
end

function get_group_idx(model::EoSModel,i,j,a,b)
    return get_group_idx(model.sites,i,j,a,b)
end

function get_group_idx(param::SiteParam,i,j,a,b)
    site_translator::Vector{Vector{NTuple{2,Int}}} = param.site_translator
    k,_ = site_translator[i][a]
    l,_ = site_translator[j][b]
  return k,l
end

