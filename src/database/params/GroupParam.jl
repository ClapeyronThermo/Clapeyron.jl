"""
    GroupParam
Struct holding group parameters.contains:
* `components`: a list of all components
* `groups`: a list of groups names for each component
* `grouptype`: used to differenciate between different group models.
* `i_groups`: a list containing the number of groups for each component
* `n_groups`: a list of the group multiplicity of each group corresponding to each group in `i_groups`
* `flattenedgroups`: a list of all unique groups--the parameters correspond to this list
* `n_flattenedgroups`: the group multiplicities corresponding to each group in `flattenedgroups`
You can create a group param by passing a `Vector{Tuple{String, Vector{Pair{String, Int64}}}}.
For example:
```julia-repl
julia> grouplist = [
           ("ethanol", ["CH3"=>1, "CH2"=>1, "OH"=>1]),
           ("nonadecanol", ["CH3"=>1, "CH2"=>18, "OH"=>1]),
           ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])];
julia> groups = GroupParam(grouplist)
GroupParam(:unkwown) with 3 components:
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
    grouptype::Symbol
    n_groups::Array{Array{Int,1},1}
    i_groups::Array{Array{Int,1},1}
    flattenedgroups::Array{String,1}
    n_flattenedgroups::Array{Array{Int,1},1}
    n_groups_cache::PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{Float64}, SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}}
    sourcecsvs::Array{String,1}
end

function GroupParam(imput::PARSED_GROUP_VECTOR_TYPE)
    return GroupParam(imput,:unknown,String[])
end


#given components, groups, n_groups, reconstitute GroupParam
function recombine!(param::GroupParam)
    components = param.components
    groups = param.groups
    n_groups = param.n_groups
    ℂ = length(components)
    𝔾 = length(groups)
    #initialization of flattenedgroups
    #flattenedgroups = unique!(reduce(vcat,groups))
    
    flattenedgroups = param.flattenedgroups
    resize!(flattenedgroups,0)
    for group ∈ groups
        append!(flattenedgroups,group)
    end
    unique!(flattenedgroups)
    #initialization of i_groups
    #i_groups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ groups]
    i_groups = param.i_groups
    resize!(i_groups,ℂ)
    for i in 1:ℂ
        group = groups[i]
        if !isassigned(i_groups,i)
            i_groups[i] = Int[]
        end
        i_group = i_groups[i]
        resize!(i_group,length(group))
        for j in eachindex(i_group)
            i_group[j] = findfirst(isequal(group[j]), flattenedgroups)
        end
    end

    #initialization of n_flattenedgroups, n_groups_cache
    #n_flattenedgroups = [zeros(Int,len_flattenedgroups) for _ ∈ 1:length(input)]
    #=for i in 1:length(input)
        setindex!(n_flattenedgroups[i],n_groups[i],i_groups[i])
        setindex!(n_groups_cache[i],n_groups[i],i_groups[i])
    end =#
    flat𝔾 = length(flattenedgroups)
    n_flattenedgroups = param.n_flattenedgroups
    
    #resizing of Packed Vector of Vectors
    n_groups_cache = param.n_groups_cache
    _p,_v = n_groups_cache.p,n_groups_cache.v
    resize!(_v,flat𝔾*ℂ)
    resize!(_p,ℂ + 1)
    resize!(n_flattenedgroups,ℂ)
    fill!(_v,0.0)

    for i in eachindex(_p)
        _p[i] = 1 + (i-1)*flat𝔾
    end

    for i in 1:ℂ
        if !isassigned(n_flattenedgroups,i)
            n_flattenedgroups[i] = Float64[]
        end
        n_flatgroup = n_flattenedgroups[i]
        resize!(n_flatgroup,flat𝔾)
        n_flatgroup .= 0.0
        setindex!(n_flatgroup,n_groups[i],i_groups[i])
        setindex!(n_groups_cache[i],n_groups[i],i_groups[i])
    end

    return param
end

function GroupParam(input::PARSED_GROUP_VECTOR_TYPE,grouptype::Symbol,sourcecsvs::Vector{String})
    components = [first(i) for i ∈ input]
    raw_groups =  [last(i) for i ∈ input]
    groups = [first.(grouppairs) for grouppairs ∈ raw_groups]
    n_groups = [last.(grouppairs) for grouppairs ∈ raw_groups]
    flattenedgroups = String[]
    i_groups = Vector{Vector{Int}}(undef,0)
    n_flattenedgroups = Vector{Vector{Int}}(undef,0)
    n_groups_cache = PackedVofV(Int[],Float64[])

    param =  GroupParam(components,
    groups,
    grouptype,
    n_groups,
    i_groups,
    flattenedgroups,
    n_flattenedgroups,
    n_groups_cache,
    sourcecsvs)
    #do the rest of the work here
    recombine!(param)
    return param
end



function Base.show(io::IO, mime::MIME"text/plain", param::GroupParam)
    print(io,"GroupParam(:",param.grouptype,") ")
    len = length(param.components)
    println(io,"with ", len, " component", ifelse(len==1, ":", "s:"))
    show_groups(io,param)
end

function Base.show(io::IO, param::GroupParam)
    print(io,"GroupParam[")
    function wrap_print(io,val)
        print(io,'[')
        __show_group_i(io,val)
        print(io,']')
    end
    show_pairs(io,param.components,zip(param.groups,param.n_groups)," => ",wrap_print,pair_separator = ", ")
    print(io,"]")
end