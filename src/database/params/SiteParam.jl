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
end




function SiteParam(components::Vector{String},sites::Array{Array{String,1},1},n_sites::Vector{Vector{Int}},sourcecsvs = String[])
    n_sites = PackedVectorsOfVectors.pack(n_sites)
    param = SiteParam(components, 
    sites, 
    n_sites,
    Vector{Vector{Int}}(undef,0), 
    String[],
    Vector{Vector{Int}}(undef,0),
    Vector{Vector{Int}}(undef,0),
    sourcecsvs)
        
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
    raw_sites =  [last(i) for i ∈ input]
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
    String[])
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
            i_site[j] = findfirst(isequal(site[j]), flattenedsites)
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

function Compressed4DMatrix(param::SiteParam)
    return assoc_similar(param,Float64)
end