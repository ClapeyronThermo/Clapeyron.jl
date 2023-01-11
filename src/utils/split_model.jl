function each_split_model(param::AbstractVector,I)
    val = param[I]
    eltype(param) <: AbstractArray && return deepcopy(val)
    return val
end

function each_split_model(param::AbstractMatrix,I)
    val =  param[I,I]
    eltype(param) <: AbstractArray && return deepcopy(val)
    return val
end

function each_split_model(y::SparseMatrixCSC{<:AbstractVector},I)
    x = y[I,I]
    m,n,colptr,rowval,nzval = x.m,x.n,x.colptr,x.rowval,x.nzval
    return SparseMatrixCSC(m,n,colptr,rowval,nzval)
end

function each_split_model(y::SparsePackedMofV,I)
    idx = y.idx[I,I]
    if iszero(length(y.storage))
        return SparsePackedMofV(y.storage,idx)
    end

    if iszero(nnz(idx))
        st = y.storage
        storage = PackedVofV([1],zeros(eltype(st.v),0))
        return SparsePackedMofV(storage,idx)
    else
        str = y.storage[nnz(idx)]
        storage = PackedVectorsOfVectors.pack(str)
        return SparsePackedMofV(storage,idx)
    end
end

function each_split_model(param::PackedVofV,I)
    val =  PackedVectorsOfVectors.pack(param[I])
    return val
end

function each_split_model(assoc::Compressed4DMatrix{T},I) where T
    len = length(assoc.values)
    iszero(len) && return Compressed4DMatrix{T}()
    old_idx = assoc.outer_indices
    idx_bool = findall(x -> (first(x) ∈ I)&(last(x) ∈ I),old_idx)
    iszero(length(idx_bool)) && return Compressed4DMatrix{T}()
    values = assoc.values[idx_bool]
    outer_indices = assoc.outer_indices[idx_bool]
    inner_indices = assoc.inner_indices[idx_bool]
    out_val = length(I)
    outer_size = (out_val,out_val)
    inner_size = assoc.inner_size
    len2 = length(outer_indices)
    for i ∈ 1:len2
        i1,j1 = outer_indices[i]
        i2,j2 = findfirst(==(i1),I),findfirst(==(j1),I)
        outer_indices[i] = (i2,j2)
    end
    return Compressed4DMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
end

function each_split_model(param::SingleParameter,I)
    return SingleParameter(
        param.name,
        param.components[I],
        each_split_model(param.values,I),
        param.ismissingvalues[I],
        param.sourcecsvs,
        param.sources
    )
end

function each_split_model(param::PairParameter,I)
    value = each_split_model(param.values,I)
    ismissingvalues = param.ismissingvalues[I,I]
    components = param.components[I]
    res = PairParameter(
            param.name,
            components,
            value,
            ismissingvalues,
            param.sourcecsvs,
            param.sources
            )
    return res
end

function each_split_model(param::AssocParam,I)
    _value  = each_split_model(param.values,I)
    return AssocParam(
            param.name,
            param.components[I],
            _value,
            param.sites[I],
            param.sourcecsvs,
            param.sources
            )
end

function each_split_model(param::GroupParam,I)
    grouptype = param.grouptype
    components = param.components[I]
    groups = param.groups[I]
    n_groups = param.n_groups[I]
    sourcecsvs = param.sourcecsvs

    #unique, but without allocating sets.
    idx = zeros(Int,length(param.flattenedgroups))
    for i in I
        group_i = param.groups[i]
        for k in 1:length(group_i)
            j = findfirst(==(group_i[k]),param.flattenedgroups)
            idx[j] = j
        end
    end
    zero_idx = iszero.(idx)
    nonzero_idx = @. !zero_idx
    _idx = view(idx,nonzero_idx)

    len_groups = length(_idx)

    flattenedgroups = param.flattenedgroups[_idx]
    i_groups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ groups]
    n_flattenedgroups = Vector{Vector{Int64}}(undef,length(I))
    for (k,i) in pairs(I)
        pii = param.n_flattenedgroups[i]
        n_flattenedgroups[k] = pii[_idx]
    end
    n_groups_cache  = PackedVectorsOfVectors.packed_fill(0.0,(length(ni) for ni in n_flattenedgroups))

    for (k,i) in pairs(I)
        pii = param.n_groups_cache[i]
        true_n = @view(pii[_idx])
        n_groups_cache[k] .= true_n
    end

    return GroupParam(
        components,
        groups,
        grouptype,
        n_groups,
        i_groups,
        flattenedgroups,
        n_flattenedgroups,
        n_groups_cache,
        sourcecsvs)
end
"""
    split_model(model::EoSModel)
Takes in a model for a multi-component system and returns a vector of model for each pure system.
## Example:
```julia-repl
julia> gerg2 = GERG2008(["propane","pentane"])
GERG008 model with 2 components:
"propane"
"pentane"
julia> split_model(gerg2)
2-element Vector{GERG2008}:
 GERG2008("propane")
 GERG2008("pentane")
```
"""
function split_model end

"""
    is_splittable(model)::Bool
Trait to determine if a `EoSModel` should be splitted by itself or can be simply filled into a vector.
This is useful in the case of models without any parameters, as those models are impossible by definition to split, because they don't have any underlying data.
The Default is `is_splittable(model) = true`.
"""
is_splittable(model) = true
is_splittable(null::Union{Nothing,Missing}) = false
is_splittable(::Number) = false
is_splittable(::String) = false

function split_model(param::SingleParameter,
    splitter =split_model(1:length(param.components)))
    return [each_split_model(param,i) for i ∈ splitter]
end

#this conversion is lossy, as interaction between two or more components are lost.

function split_model(param::PairParameter,
    splitter = split_model(1:length(param.components)))
    return [each_split_model(param,i) for i ∈ splitter]
end

function split_model(param::AbstractVector,splitter = ([i] for i ∈ 1:length(param)))
    return [each_split_model(param,i) for i ∈ splitter]
end

function split_model(param::UnitRange{Int},splitter = ([i] for i ∈ 1:length(param)))
    return [1:length(i) for i ∈ splitter]
end

#this conversion is lossy, as interaction between two or more components are lost.
#also, this conversion stores the site values for other components. (those are not used)
function split_model(param::AssocParam{T},
    splitter = split_model(1:length(param.components))) where T
    return [each_split_model(param,i) for i ∈ splitter]
end

#this param has a defined split form
function split_model(groups::GroupParam,
    splitter = split_model(collect(1:length(groups.components))))
    return [each_split_model(groups,i) for i in splitter]
end

function split_model(param::SiteParam,
    splitter = split_model(param.components))
    function generator(I)
        return SiteParam(
            param.components[I],
            param.sites[I],
            each_split_model(param.n_sites,I),
            param.i_sites[I],
            param.flattenedsites,
            param.n_flattenedsites[I],
            param.i_flattenedsites[I],
        param.sourcecsvs)
    end
    return [generator(i) for i ∈ splitter]
end

function split_model(Base.@nospecialize(params::EoSParam),splitter)
    T = typeof(params)
    split_paramsvals = (split_model(getfield(params,i),splitter) for i  ∈ fieldnames(T))
    return T.(split_paramsvals...)
end

function group_splitter(group,splitted_groups)
    flattenedgroups = group.flattenedgroups
    res = Vector{Vector{Int64}}(undef,length(splitted_groups))
    for (i,groupi) in pairs(splitted_groups)
        res[i] = indexin(groupi.flattenedgroups,flattenedgroups)
    end
    return res
end

split_model(model::EoSModel,subset=nothing) = auto_split_model(model,subset)

function auto_split_model(Base.@nospecialize(model::EoSModel),subset=nothing)
        
    try
        allfields = Dict{Symbol,Any}()
        #if has_groups(typeof(model))
        #    raw_splitter = model.groups.i_groups
        #    subset !== nothing && throw("using subsets is not supported with Group Contribution models")
        #else
        #    raw_splitter = split_model(Vector(1:length(model.components)))
        #end

        raw_splitter = split_model(Vector(1:length(model.components)))

        if subset === nothing
            splitter = raw_splitter
        elseif eltype(subset) <: Integer
            splitter = raw_splitter[subset]
        elseif eltype(subset) <: AbstractVector
            splitter = subset
        else
            throw("invalid type of subset.")
        end

        if hasfield(typeof(model),:groups)
            gc_split = split_model(model.groups,splitter)
            allfields[:groups] = gc_split
            allfields[:components] = split_model(model.groups.components,splitter)
            splitter = group_splitter(model.groups,gc_split)
        end

        len = length(splitter)
        M = typeof(model)
        allfieldnames = fieldnames(M)

        #add here any special keys, that behave as non_splittable values
        for modelkey in [:references]
            if modelkey in allfieldnames
                if !haskey(allfields,modelkey)
                    allfields[modelkey] = fill(getproperty(model,modelkey),len)
                end
            end
        end

        for modelkey ∈ allfieldnames
            if !haskey(allfields,modelkey)
                modelx = getproperty(model,modelkey)
                if is_splittable(modelx)
                    if modelx isa EoSModel
                        allfields[modelkey]= split_model(modelx,subset)
                    else
                        allfields[modelkey]= split_model(modelx,splitter)
                    end
                else
                    allfields[modelkey] = fill(modelx,len)
                end
            end
        end

        return [M((allfields[k][i] for k ∈ fieldnames(M))...) for i ∈ 1:len]::Vector{M}
    catch e
        M = typeof(model)
        @error "$M cannot be splitted"
        rethrow(e)
    end
end

##fallback,around 50 times slower if there is any need to read csvs

function simple_split_model(Base.@nospecialize(model::EoSModel),subset = nothing)
    MODEL = typeof(model)
    pure = Vector{MODEL}(undef,0)
    if subset === nothing
        comps = model.components
    else
        comps = model.components[subset]
    end
    for comp ∈ comps
        push!(pure,MODEL([comp]))
    end
    return pure
end

export split_model

#=
function _split_model(groups::GroupParam)
    len = length(groups.components)
    function generator(i)
        return GroupParam(
        [groups.components[i]],
        [groups.groups[i]],
        [groups.n_groups[i]],
        [collect(1:length(groups.n_groups[i]))],
        groups.flattenedgroups[groups.i_groups[i]],
        [groups.n_groups[i]],
        each_split_model(groups.n_groups_cache,groups.i_groups[i]),
        1:length(groups.n_groups[i]),
        groups.sourcecsvs)
    end
    [generator(i) for i ∈ 1:len]
end
=#