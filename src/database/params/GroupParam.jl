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
    n_groups_cache::PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}}
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
    n_groups_cache = PackedVectorsOfVectors.packed_fill(0.0,fill(len_flattenedgroups,length(input)))
    for i in 1:length(input)
        setindex!(n_flattenedgroups[i],n_groups[i],i_groups[i])
        setindex!(n_groups_cache[i],n_groups[i],i_groups[i])
    end
    
    return GroupParam(components, 
    groups, 
    n_groups,
    i_groups, 
    flattenedgroups,
    n_flattenedgroups,
    n_groups_cache,
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
