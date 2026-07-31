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
    n_sites::PackedVector{Int}
    i_sites::Array{Array{Int,1},1}
    flattenedsites::Array{String,1}
    n_flattenedsites::Array{Array{Int,1},1}
    i_flattenedsites::Array{Array{Int,1},1}
    sourcecsvs::Array{String,1}
    site_translator::Union{Nothing,Vector{Vector{NTuple{2,Int}}}}
end

function SiteParam(components,sites,n_sites,sourcecsvs = String[],site_translator = nothing)
    if n_sites isa PackedVector{Int}
        _n_sites = n_sites
    else
        _n_sites = PackedVectorsOfVectors.pack(n_sites)
    end

    param = SiteParam(components,
    sites,
    _n_sites,
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

function SiteParam(input::Vector{Tuple{String, Vector{Pair{String,T}}}},sourcecsvs=String[]) where {T<:Number}
    components = [first(i) for i ∈ input]
    raw_sites = [last(i) for i ∈ input]
    sites = [first.(sitepairs) for sitepairs ∈ raw_sites]
    n_sites = [last.(sitepairs) for sitepairs ∈ raw_sites]
    return SiteParam(components,sites,n_sites,sourcecsvs)
end

function SiteParam(components::Vector{String},sourcecsvs::Vector{String},::Val{SITE_TRANSLATOR}) where SITE_TRANSLATOR
    n = length(components)
    sites = [String[] for _ ∈ 1:n]
    n_sites = packed_zeros!(Int,zeros(Int,n))
    sourcecsvs = String[]
    if SITE_TRANSLATOR
        site_translator = nothing
    else
        site_translator = Vector{NTuple{2,Int}}[]
    end
    SiteParam(components,sites,n_sites,sourcecsvs,site_translator)
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

Returns a `Clapeyron.Compressed4DMatrix` with the same number of components as the input `AssocParam`, with the same element type as `𝕋`.
All site combinations are filled.
"""
function assoc_similar(param::SiteParam,::Type{𝕋}) where 𝕋 <:Number
    Compressed4DMatrices.c4d_from_site_offsets(𝕋,copy(param.n_sites.p))
end

assoc_similar(m) = assoc_similar(m,paramtype(m))

function Compressed4DMatrix(param::SiteParam)
    return assoc_similar(param,Float64)
end

#build dense assocparam from sites.
function AssocParam{T}(name,sites::SiteParam) where T <: Number
    values = assoc_similar(sites,T)
    return AssocParam(name,sites.components,values,sites.sites)
end

#=
Utilities to create "group-component" sites

if we re-index the association sites that are present in a GC-based approach, the result is the same.
=#
function gc_to_comp_sites!(out::SiteParam,sites::SiteParam,groups::GroupParameter)
    #translates association sites defined per-group (sites indexed by (g,s))
    #into association sites defined per-component (indexed by (c,s)), by
    #concatenating, for each component, the site lists of its constituent groups.
    #site_translator[c][s] = (g,s_g) lets an AssocParam computed in GC-space be
    #looked up in component-space later.
    comps = groups.components
    nc = length(comps)
    
    if out.components !== comps
        resize!(out.components,nc)
        out.components .= comps
    end
    
    site_translator = out.site_translator::Vector{Vector{NTuple{2,Int}}}
    resize!(site_translator,nc)

    out_sites = out.sites
    resize!(out_sites,nc)

    n_sites_p = out.n_sites.p
    n_sites_v = out.n_sites.v
    resize!(n_sites_p,nc + 1)
    resize!(n_sites_v,0)
    n_sites_p[1] = 1
    offset = 1

    is_empty = iszero(sites.n_sites.v)
    if length(sites.sourcecsvs) != 0 && sites.sourcecsvs !== out.sourcecsvs
        resize!(out.sourcecsvs,length(sites.sourcecsvs))
        out.sourcecsvs .= sites.sourcecsvs
    end

    gc_sitenames = sites.sites #site names, indexed by flattened group
    gc_names = groups.flattenedgroups
    i_groups_per_comp = groups.i_groups #i_groups_per_comp[i][k] = flattened-group-index of the k-th group in component i
    n_groups = groups.n_groups #gc_n_groups[i][k] = multiplicity of that group in component i

    for i in 1:nc
        if isassigned(out_sites,i)
            sites_i = out_sites[i]
            resize!(sites_i,0)
        else
            sites_i = String[]
            out_sites[i] = sites_i
        end

        if isassigned(site_translator,i)
            translator_i = site_translator[i]
            resize!(translator_i,0)
        else
            translator_i = Tuple{Int,Int}[]
            site_translator[i] = translator_i
        end

        if !is_empty
            i_groups_i = i_groups_per_comp[i]
            n_groups_i = n_groups[i]

            for k in 1:length(i_groups_i)
                g = i_groups_i[k] #flattened-group-index
                n_gc_k = n_groups_i[k] #amount of groups k in component i
                gname = gc_names[g]
                g_sitenames = gc_sitenames[g]
                gc_n_sites = sites.n_sites[g]

                for s in 1:length(g_sitenames)
                    push!(sites_i,gname * '/' * g_sitenames[s])
                    push!(n_sites_i, gc_n_sites[s]*n_gc_k) #amount of sites equal to n_gc * n_sites(gc)
                    push!(translator_i,(g,s))
                    offset += 1
                end
            end
        end

        n_sites_p[i+1] = offset
    end

    recombine!(out)
    return out
end

function _gc_to_comp_sites!(m::Compressed4DMatrix{T1},m_gc::Compressed4DMatrix{T2},sites::SiteParam) where {T1,T2}

    if length(sites.n_sites.v) == 0
        resize!(m.values,0)
        resize!(m.indices,0)
        resize!(m.site_offsets,1)
        m.site_offsets[1] = 1
        return m
    end

    site_offsets = sites.n_sites.p
    resize!(m.site_offsets,length(site_offsets))
    m.site_offsets .= site_offsets
    length(m_gc.values) == 0 && return m
    site_translator = sites.site_translator
    Compressed4DMatrices.extend!(m)
    for (idx,(i,j),(a,b)) in indices(m)
        igc,jgc,agc,bgc = get_group_ijab(sites,i,j,a,b)
        idx_gc = validindex(m_gc,(igc,jgc,agc,bgc))
        if !iszero(idx_gc)
            m[idx] = convert(T1,m_gc[idx_gc])
        end
    end
    dropzeros!(m)
    return m
end

_gc_to_comp_sites!(m::AssocParam,m_gc::AssocParam,sites::SiteParam) = _gc_to_comp_sites!(m.values,m_gc.values,sites)
_gc_to_comp_sites!(m::Compressed4DMatrix,m_gc::AssocParam,sites::SiteParam) = _gc_to_comp_sites!(m,m_gc.values,sites)
_gc_to_comp_sites!(m::AssocParam,m_gc::Compressed4DMatrix,sites::SiteParam) = _gc_to_comp_sites!(m.values,m_gc,sites)

function gc_to_comp_sites(sites::SiteParam,groups::GroupParameter)
    out = SiteParam(groups.components,sites.sourcecsvs,Val(true))
    return gc_to_comp_sites!(out,sites,groups)
end

gc_to_comp_sites!(out,param::Compressed4DMatrix,sites::SiteParam) = _gc_to_comp_sites!(out,param,sites)

function gc_to_comp_sites!(out,param::AssocParam,sites::SiteParam)
    m = _gc_to_comp_sites!(out,param,sites)
    return AssocParam(param.name,sites.components,m,sites.sites,param.sourcecsvs,param.sources)
end

function gc_to_comp_sites(param::AssocParam{T},sites::SiteParam) where T
    m = _gc_to_comp_sites!(Compressed4DMatrix{T}(),param,sites)
    return AssocParam(param.name,sites.components,m,sites.sites,param.sourcecsvs,param.sources)
end

gc_to_comp_sites(param::Compressed4DMatrix{T},sites::SiteParam) where T = _gc_to_comp_sites!(Compressed4DMatrix{T}(),param,sites)

#=
site_translator utilities
=#
function get_group_idx(model::EoSModel,i,j,a,b)
    return get_group_idx(model.sites,i,j,a,b)
end

function get_group_idx(param::SiteParam,i,j,a,b)
    site_translator::Vector{Vector{NTuple{2,Int}}} = param.site_translator
    k,_ = site_translator[i][a]
    l,_ = site_translator[j][b]
  return k,l
end

function get_group_ijab(param::SiteParam,i,j,a,b)
    site_translator::Vector{Vector{NTuple{2,Int}}} = param.site_translator
    k,s1 = site_translator[i][a]
    l,s2 = site_translator[j][b]
    return k,l,s1,s2
end
