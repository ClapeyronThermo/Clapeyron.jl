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


function Base.show(io::IO, mime::MIME"text/plain", param::SiteParam)
    print(io,"SiteParam ")
    len = length(param.components)
    println(io,"with ", len, " component", ifelse(len==1, ":", "s:"))
    
    for i in 1:length(param.components)
        
        print(io, " \"", param.components[i], "\": ")
        firstloop = true
        if length(param.n_sites[i]) == 0
            print(io,"(no sites)")
        end
        for j in 1:length(param.n_sites[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.sites[i][j], "\" => ", param.n_sites[i][j])
            firstloop = false
        end
        i != length(param.components) && println(io)
    end 
end

function Base.show(io::IO, param::SiteParam)
    print(io,"SiteParam[")
    len = length(param.components)
    
    for i in 1:length(param.components)
        
        print(io, "\"", param.components[i], "\" => [")
        firstloop = true
    
        for j in 1:length(param.n_sites[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.sites[i][j], "\" => ", param.n_sites[i][j])
            firstloop = false
        end
        print(io,']')
        i != length(param.components) && print(io,", ")
    end
    print(io,"]")
end


function SiteParam(pairs::Dict{String,SingleParam{Int}},allcomponentsites)
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    sites = allcomponentsites
    ncomps = length(components)
    sourcecsvs = String[]
    for x in values(pairs)
        vcat(sourcecsvs,x.sourcecsvs)  
    end
    if length(sourcecsvs) >0
        unique!(sourcecsvs)
    end
    n_sites = [[pairs[sites[i][j]].values[i] for j ∈ 1:length(sites[i])] for i ∈ 1:ncomps]  # or groupsites
    length_sites = [length(componentsites) for componentsites ∈ sites]
    i_sites = [1:length_sites[i] for i ∈ 1:ncomps]
    flattenedsites = unique!(reduce(vcat,sites,init = String[]))
    len_flattenedsites = length(flattenedsites)
    #i_flattenedsites = 1:len_flattenedsites
    n_flattenedsites = [zeros(Int,len_flattenedsites) for _ ∈ 1:ncomps]
    i_flattenedsites = [zeros(Int,len_flattenedsites) for _ ∈ 1:ncomps]
    for i in 1:ncomps
        setindex!.(n_flattenedsites,n_sites,i_sites)
        flatsites_i = i_flattenedsites[i]
        for (aidx,a) in Base.pairs(i_sites[i])
            flatsites_i[a] = aidx
        end
    end
    return SiteParam(components, 
    sites, 
    PackedVectorsOfVectors.pack(n_sites),
    i_sites, 
    flattenedsites,
    n_flattenedsites, 
    i_flattenedsites,
    sourcecsvs)
end

function SiteParam(input::PARSED_GROUP_VECTOR_TYPE,sourcecsvs::Vector{String}=String[])
    components = [first(i) for i ∈ input]
    raw_sites =  [last(i) for i ∈ input]
    sites = [first.(sitepairs) for sitepairs ∈ raw_sites]
    n_sites = [last.(sitepairs) for sitepairs ∈ raw_sites]
    flattenedsites = unique!(reduce(vcat,sites,init = String[]))
    i_sites = [[findfirst(isequal(site), flattenedsites) for site ∈ componentsites] for componentsites ∈ sites]
    len_flattenedsites = length(flattenedsites)
    i_flattenedsites = [zeros(Int,len_flattenedsites) for _ ∈ 1:length(input)]
    n_flattenedsites = [zeros(Int,len_flattenedsites) for _ ∈ 1:length(input)]
    for i in 1:length(input)
        setindex!.(n_flattenedsites,n_sites,i_sites)
        flatsites_i = i_flattenedsites[i]
        for (aidx,a) in Base.pairs(i_sites[i])
            flatsites_i[a] = aidx
        end
    end

    return SiteParam(components, 
    sites, 
    PackedVectorsOfVectors.pack(n_sites),
    i_sites, 
    flattenedsites,
    n_flattenedsites, 
    i_flattenedsites,
    sourcecsvs)
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