abstract type ClapeyronParam end
abstract type EoSParam end
export EoSParam

struct SingleParameter{T,V<:AbstractVector{T}} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::V
    ismissingvalues::Array{Bool,1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

"""
    SingleParam{T}

Struct designed to contain single parameters. Basically a vector with some extra info.

## Creation:
```julia-repl
julia> mw = SingleParam("molecular weight",["water","ammonia"],[18.01,17.03])
SingleParam{Float64}("molecular weight") with 2 components:
 "water" => 18.01
 "ammonia" => 17.03

julia> mw.values
2-element Vector{Float64}:
 18.01
 17.03

julia> mw.components
2-element Vector{String}:
 "water"
 "ammonia"

julia> mw2 = SingleParam(mw,"new name")
SingleParam{Float64}("new name") with 2 components:
 "water" => 18.01
 "ammonia" => 17.03

julia> has_oxigen = [true,false]; has_o = SingleParam(mw2,has_oxigen)
SingleParam{Bool}("new name") with 2 components:
 "water" => true
 "ammonia" => false

```

## Example usage in models:

```
function molecular_weight(model,molar_frac)
    mw = model.params.mw.values
    res = zero(eltype(molarfrac))
    for i in @comps #iterating through all components
        res += molar_frac[i]*mw[i]
    end
    return res
end
```
"""
const SingleParam{T} = SingleParameter{T,Vector{T}} where T

SingleParam(name,components,values,missingvals,src,sourcecsv) = SingleParameter(name,components,values,missingvals,src,sourcecsv)
function Base.convert(::Type{SingleParam{String}},param::SingleParam{<:AbstractString})::SingleParam{String}
    values = String.(param.values)
    return (param.name,param.components,values,param.missingvals,param.src,param.sourcecsv)
end
function Base.show(io::IO, param::SingleParameter)
    print(io, typeof(param), "(\"", param.name, "\")[")
    for component in param.components
        component != first(param.components) && print(io, ",")
        print(io, "\"", component, "\"")
    end
    print(io, "]")
end

function Base.show(io::IO, ::MIME"text/plain", param::SingleParameter)
    len = length(param.values)
    print(io, "SingleParam{",eltype(param.values), "}(\"", param.name)
    println(io, "\") with ", len, " component", ifelse(len==1, ":", "s:"))
    i = 0
    for (name, val, miss) in zip(param.components, param.values, param.ismissingvalues)
        i += 1
        if i > 1
            println(io)
        end
        if miss == false
            if typeof(val) <: AbstractString
                print(io, " \"", name, "\" => \"", val, "\"")
            else
                print(io, " \"", name, "\" => ", val)
            end
        else
            print(io, " \"", name, " => -")
        end
    end
end

function SingleParam(x::SingleParameter,name=x.name)
    return SingleParam(name, x.components,deepcopy(x.values), deepcopy(x.ismissingvalues), x.sourcecsvs, x.sources)
end

#a barebones constructor, in case we dont build from csv
function SingleParam(
    name::String,
    components::Vector{String},
    values::Vector{T},
    sourcecsvs = String[],
    sources = String[]
    ;default = _zero(T)) where T
    _values,_ismissingvalues = defaultmissing(values,default)
    TT = eltype(_values)
    return  SingleParam{TT}(name,components, _values, _ismissingvalues, sourcecsvs, sources)
end


function SingleParam(x::SingleParameter, v::Vector)
    _values,_ismissingvalues = defaultmissing(v)
    return SingleParam(x.name, x.components,_values, _ismissingvalues , x.sourcecsvs, x.sources)
end


struct PairParameter{T,V<:AbstractMatrix{T},D} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::V
    diagvalues::D
    ismissingvalues::Array{Bool,2}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end
"""
    PairParam{T}

Struct designed to contain pair data. used a matrix as underlying data storage.

## Creation:
```julia-repl
julia> kij = PairParam("interaction params",["water","ammonia"],[0.1 0.0;0.1 0.0])
PairParam{Float64}["water", "ammonia"]) with values:
2×2 Matrix{Float64}:
 0.1  0.0
 0.1  0.0

julia> kij.values
2×2 Matrix{Float64}:
 0.1  0.0
 0.1  0.0

julia> kij.diagvalues
2-element view(::Vector{Float64}, 
1:3:4) with eltype Float64:
 0.1
 0.0
```

## Example usage in models:

```julia
#lets compute ∑xᵢxⱼkᵢⱼ
function alpha(model,x)
    kij = model.params.kij.values
    ki = model.params.kij.diagvalues
    res = zero(eltype(molarfrac))
    for i in @comps 
        @show ki[i] #diagonal values
        for j in @comps 
            res += x[i]*x[j]*kij[i,j]
        end
    end
    return res
end
```
"""
const PairParam{T} = PairParameter{T,Matrix{T},SubArray{T, 1, Vector{T}, Tuple{StepRange{Int64, Int64}}, true}} where T

PairParam(name,components,values,diagvals, missingvals,src,sourcecsv) = PairParameter(name,components,values,diagvals,missingvals,src,sourcecsv)

#unsafe constructor
function PairParam(name,components,values)
    missingvals = fill(false,size(values))
    diagvals = view(values, diagind(values))
    src = String[]
    sourcecsv = String[]
    return PairParam(name,components,values,diagvals, missingvals,src,sourcecsv)
end
function PairParam(name::String,
                    components::Array{String,1},
                    values::Array{T,2},
                    ismissingvalues = fill(false,length(components),length(components)),
                    sourcecsvs::Array{String,1} = String[], 
                    sources::Array{String,1} = String[]) where T
    
    _values,_ismissingvalues = defaultmissing(values)
    diagvalues = view(_values, diagind(_values))
    if !all(ismissingvalues)
        _ismissingvalues = ismissingvalues
    end
    return PairParam(name, components,_values, diagvalues, _ismissingvalues, sourcecsvs, sources)
end

function PairParam(x::PairParameter,name::String=x.name)
    values = deepcopy(x.values)
    diagvalues = view(values,diagind(values))
    return PairParam(name, x.components,values ,diagvalues, deepcopy(x.ismissingvalues), x.sourcecsvs, x.sources)
end

function PairParam(x::SingleParameter,name::String=x.name)
    pairvalues = singletopair(x.values,missing)
    for i in 1:length(x.values)
        if x.ismissingvalues[i]
            pairvalues[i,i] = missing
        end
    end
    _values,_ismissingvalues = defaultmissing(pairvalues)
    diagvalues = view(_values, diagind(_values))
    return PairParam(name, x.components, _values,diagvalues,_ismissingvalues,x.sourcecsvs, x.sources)
end
#=
function PairParam(x::PairParam, v::Matrix,name::String=x.name)
    return PairParam(name, x.components,deepcopy(v), x.sourcecsvs, x.sources)
end
function PairParam(x::SingleParam, v::Vector,name::String=x.name)
    pairvalues = singletopair(v,missing)
    return PairParam(x.name, x.components, pairvalues,x.sourcecsvs, x.sources)
end
function PairParam(x::SingleParam, v::Matrix,name::String=x.name)
    return PairParam(x.name, x.components, deepcopy(v),x.sourcecsvs, x.sources)
end
=#
#barebones constructor by list of pairs.


function Base.show(io::IO,mime::MIME"text/plain",param::PairParameter) 
    print(io,"PairParam{",eltype(param.values),"}")
    show(io,param.components)
    println(io,") with values:")
    show(io,mime,param.values)
end

function Base.show(io::IO,param::PairParameter)
    print(io, "PairParam{",eltype(param.values),"}", "(\"", param.name, "\")[")
    print(io,Base.summary(param.values))
    print(io,"]")
end
"""
    AssocParam{T}

Struct holding association parameters.
"""
struct AssocParam{T} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::Compressed4DMatrix{T,Vector{T}}
    sites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function AssocParam(name::String,components::Vector{String},values::MatrixofMatrices,allcomponentsites,sourcecsvs,sources) where T
    _values = Compressed4DMatrix(values)
    return AssocParam(name,components,_values,allcomponentsites,sourcecsvs,sources)
end

#=
function AssocParam(x::AssocParam{T}) where T
    return AssocParam{T}(x.name,x.components, deepcopy(x.values), x.sites, x.sourcecsvs, x.sources)
end

function AssocParam{T}(x::AssocParam, v::Matrix{Matrix{T}}) where T
    return AssocParam{T}(x.name, x.components,Compressed4DMatrix(v), x.sites, x.sourcecsvs, x.sources)
end

function AssocParam{T}(name::String,components::Vector{String}) where T
    n = length(components)
    return AssocParam{T}(name, 
    components,
    Compressed4DMatrix{T}(),
    [String[] for _ ∈ 1:n], 
    String[],
    String[])
end
=#
function Base.show(io::IO,mime::MIME"text/plain",param::AssocParam{T}) where T
    print(io,"AssocParam{",string(T),"}")
    print(io,param.components)
    println(io,") with values:")
    comps = param.components
    vals = param.values
    sitenames = param.sites
    for (idx,(i,j),(a,b)) in indices(vals)
        try
        s1 = sitenames[i][a]
        s2 = sitenames[j][b]
        print(io,"(\"",comps[i],"\", \"",s1,"\")")
        print(io," >=< ")
        print(io,"(\"",comps[j],"\", \"",s2,"\")")
        print(io,": ")
        println(io,vals.values[idx])
        catch
        println("error at i = $i, j = $j a = $a, b = $b")
        end
    end
end

function Base.show(io::IO,param::AssocParam)
    print(io, typeof(param), "(\"", param.name, "\")")
    print(io,param.values.values)
end
const PARSED_GROUP_VECTOR_TYPE =  Vector{Tuple{String, Vector{Pair{String, Int64}}}}

"""
    GroupParam

Struct holding group parameters.contains:
* `components`: a list of all components
* `groups`: a list of groups names for each component
* `i_groups`: a list containing the number of groups for each component
* `n_groups`: a list of the group multiplicity of each group corresponding to each group in `i_groups`
* `flattenedgroups`: a list of all unique groups--the parameters correspond to this list
* `n_flattenedgroups`: the group multiplicities corresponding to each group in `flattenedgroups`
* `i_flattenedgroups`: an iterator that goes through the indices for each flattened group

You can create a group param by passing a `Vector{Tuple{String, Vector{Pair{String, Int64}}}}.
For example:
```julia-repl
julia> grouplist = [
           ("ethanol", ["CH3"=>1, "CH2"=>1, "OH"=>1]), 
           ("nonadecanol", ["CH3"=>1, "CH2"=>18, "OH"=>1]),
           ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])];

julia> groups = GroupParam(grouplist)
GroupParam with 3 components:
 "ethanol": "CH3" => 1, "CH2" => 1, "OH" => 1
 "nonadecanol": "CH3" => 1, "CH2" => 18, "OH" => 1    
 "ibuprofen": "CH3" => 3, "COOH" => 1, "aCCH" => 1, "aCCH2" => 1, "aCH" => 4

julia> groups.flattenedgroups
7-element Vector{String}:
 "CH3"
 "CH2"
 "OH"
 "COOH"
 "aCCH"
 "aCCH2"
 "aCH"

julia> groups.i_groups
3-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [1, 2, 3]
 [1, 4, 5, 6, 7]

julia> groups.n_groups
3-element Vector{Vector{Int64}}:
 [1, 1, 1]
 [1, 18, 1]
 [3, 1, 1, 1, 4]

julia> groups.n_flattenedgroups
 3-element Vector{Vector{Int64}}:
 [1, 1, 1, 0, 0, 0, 0]
 [1, 18, 1, 0, 0, 0, 0]
 [3, 0, 0, 1, 1, 1, 4]
```

if you have CSV with group data, you can also pass those, to automatically query the missing groups in your input vector:

```julia-repl
julia> grouplist = [
           "ethanol", 
           ("nonadecanol", ["CH3"=>1, "CH2"=>18, "OH"=>1]),
           ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])];

           julia> groups = GroupParam(grouplist, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"])
           GroupParam with 3 components:
            "ethanol": "CH2OH" => 1, "CH3" => 1
            "nonadecanol": "CH3" => 1, "CH2" => 18, "OH" => 1    
            "ibuprofen": "CH3" => 3, "COOH" => 1, "aCCH" => 1, "aCCH2" => 1, "aCH" => 4
```
In this case, `SAFTGammaMie` files support the second order group `CH2OH`.

"""
struct GroupParam <: ClapeyronParam
    components::Array{String,1}
    groups::Array{Array{String,1},1}
    n_groups::Array{Array{Int,1},1}
    i_groups::Array{Array{Int,1},1}
    flattenedgroups::Array{String,1}
    n_flattenedgroups::Array{Array{Int,1},1}
    i_flattenedgroups::UnitRange{Int}
    sourcecsvs::Array{String,1}
end

function GroupParam(input::PARSED_GROUP_VECTOR_TYPE,sourcecsvs::Vector{String}=String[])
    components = [first(i) for i ∈ input]
    raw_groups =  [last(i) for i ∈ input]
    groups = [first.(grouppairs) for grouppairs ∈ raw_groups]
    n_groups = [last.(grouppairs) for grouppairs ∈ raw_groups]
    flattenedgroups = unique!(reduce(vcat,groups))
    i_groups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ groups]
    len_flattenedgroups = length(flattenedgroups)
    i_flattenedgroups = 1:len_flattenedgroups
    n_flattenedgroups = [zeros(Int,len_flattenedgroups) for _ ∈ 1:length(input)]
    for i in length(input)
        setindex!.(n_flattenedgroups,n_groups,i_groups)
    end

    return GroupParam(components, 
    groups, 
    n_groups,
    i_groups, 
    flattenedgroups,
    n_flattenedgroups, 
    i_flattenedgroups,
    sourcecsvs)
end

function Base.show(io::IO, mime::MIME"text/plain", param::GroupParam)
    print(io,"GroupParam ")
    len = length(param.components)
    println(io,"with ", len, " component", ifelse(len==1, ":", "s:"))
    
    for i in 1:length(param.components)
        
        print(io, " \"", param.components[i], "\": ")
        firstloop = true
        for j in 1:length(param.n_groups[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.groups[i][j], "\" => ", param.n_groups[i][j])
            firstloop = false
        end
        i != length(param.components) && println(io)
    end 
end

function Base.show(io::IO, param::GroupParam)
    print(io,"GroupParam[")
    len = length(param.components)
    
    for i in 1:length(param.components)
        
        print(io, "\"", param.components[i], "\" => [")
        firstloop = true
        for j in 1:length(param.n_groups[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.groups[i][j], "\" => ", param.n_groups[i][j])
            firstloop = false
        end
        print(io,']')
        i != length(param.components) && print(io,", ")
    end
    print(io,"]")
end

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

"""
    AssocOptions(;rtol = 1e-12,atol = 1e-12,max_iters = 1000,dampingfactor = 0.5)

Struct containing iteration parameters for the solver of association sites.

"""
@Base.kwdef struct AssocOptions <: ClapeyronParam
    rtol::Float64 = 1e-12
    atol::Float64 = 1e-12
    max_iters::Int = 1000
    dampingfactor::Float64 = 0.5
end

is_splittable(::AssocOptions) = false

function pack_vectors(x::AbstractVector{<:AbstractVector})
    return PackedVectorsOfVectors.pack(x)
end

function pack_vectors(x::SparseMatrixCSC{<:AbstractVector})
    return SparsePackedMofV(x)
end

function pack_vectors(param::SingleParameter{<:AbstractVector})
    name,components,vals,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    vals = pack_vectors(vals)
    return SingleParam(name,components,vals,missingvals,srccsv,src)
end

function pack_vectors(param::PairParameter{<:AbstractVector})
    name,components,vals,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    vals = pack_vectors(vals)
    return PairParam(name,components,vals,nothing,missingvals,srccsv,src)
end

function pack_vectors(params::Vararg{SingleParameter{T},N}) where {T<:Number,N}
    param = first(params)
    name,components,_,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    len = length(params)
    vals = [zeros(len) for _ in params]

    for i in 1:length(vals)
        vali = vals[i]
        for (k,par) in pairs(params)
            vali[k] = par.values[i]
        end
    end
    vals = PackedVectorsOfVectors.pack(vals)
    SingleParam(name,components,vals,missingvals,srccsv,src)
end

const PackedVectorSingleParam{T} = Clapeyron.SingleParameter{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}, PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{T}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}}

const PackedSparsePairParam{T} = Clapeyron.PairParameter{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}, SparsePackedMofV{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, 
true}, PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{T}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}}, Nothing} where T
export PackedSparsePairParam
export PackedVectorSingleParam