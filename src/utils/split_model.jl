#==
There are a lot of splitting mechanisms here. this is a sort of summary:

- each_split_model: given an index, will return a param "slice" - oblivious of GC distinction
- core_split_model: (core_split_model(x,i) for i in selected_splitter). differenciates between comp and GC models
  on EoSModels defaults to split_model.
- split_model: exposed function, only used by EoSModels. defaults to auto_split_model
- auto_split_model: recursive split_model for EoSModels it uses core_split_model
- simple_split_model: old option, reconstructs a model by calling the constructor each time, very slow.

on recursive EoSModels there is the following recursion:
    split_model
    -> auto_split_model
       -> core_split_model
          -> split_model

to allow splitting on an arbitrary struct, it is recommendable to overload core_split_model(struct::Struct,splitter::ModelSplitter)
each_split_model is used on simple structs that allow a simple splitting, it is a good abstraction but it is not recomendable to overload.
==#


#=each_split_model=#

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
    if isnothing(param.diagvalues)
        diagvalue = nothing
    else
        diagvalue = view(value,diagind(value))
    end
    ismissingvalues = param.ismissingvalues[I,I]
    components = param.components[I]
    res = PairParameter(
            param.name,
            components,
            value,
            diagvalue,
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
    i_flattenedgroups = 1:len_groups
    
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
        n_groups,
        i_groups,
        flattenedgroups,
        n_flattenedgroups,
        n_groups_cache,
        i_flattenedgroups,
        sourcecsvs)
end

function each_split_model(param::SiteParam,I)
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

#=split_model=#

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

function group_splitter(group,splitted_groups)
    flattenedgroups = group.flattenedgroups
    res = Vector{Vector{Int64}}(undef,length(splitted_groups))
    for (i,groupi) in pairs(splitted_groups)
        res[i] = indexin(groupi.flattenedgroups,flattenedgroups)
    end    
    return res
end

#=ModelSplitter=#

struct ModelSplitter
    group::Vector{String}
    comp::Vector{String}
    group_splitter::Vector{Vector{Int}}
    comp_splitter::Vector{Vector{Int}}
end

Base.length(model::ModelSplitter) = length(model.comp)

function init_splitter(x,::Nothing)
    raw_split = [[i] for i ∈ 1:length(x)]
    return _init_splitter(x,raw_split)
end

init_splitter(x,i) = _init_splitter(x,i)

function init_splitter(x::GroupParam,::Nothing)
    raw_split = [[i] for i ∈ 1:length(x.components)]
    return _init_splitter(x,raw_split)
end

function init_splitter(x,split::AbstractVector{<:Integer})
    raw_split = [[i] for i ∈ 1:length(x)]
    return _init_splitter(x,raw_split[split])
end

function _init_splitter(x::Vector{String},splitter)
    group = comp = x
    group_splitter = comp_splitter = splitter
    all_fields = Dict{Symbol,Any}()
    model_splitter = ModelSplitter(group,comp,group_splitter,comp_splitter)
    return model_splitter,all_fields
end

function _init_splitter(groups::GroupParam,subset)
    all_fields = Dict{Symbol,Any}()
    gc_split = [each_split_model(groups,i) for i in subset]
    group = groups.flattenedgroups
    comp = groups.components
    all_fields[:groups] = gc_split
    all_fields[:components] = [each_split_model(comp,i) for i in subset]
    _group_splitter = group_splitter(groups,gc_split)
    comp_splitter = subset
    model_splitter = ModelSplitter(group,comp,_group_splitter,comp_splitter)
    return model_splitter,all_fields
end

#=core_split_model=#

function splitter_index(val,splitter)
    if has_groups(val)
        return splitter.group_splitter
    else
        comps = components(val)
        isnothing(comps) && return splitter.comp_splitter
        if comps == splitter.group
            return splitter.group_splitter
        else
            return splitter.comp_splitter
        end
    end
end
function core_split_model(param::ClapeyronParam,splitter::ModelSplitter)
    idx = splitter_index(param,splitter)
    return [each_split_model(param,i) for i in idx]
end

function core_split_model(param::EoSModel,splitter::ModelSplitter)
    idx = splitter_index(param,splitter)
    return split_model(param,idx)
end

function core_split_model(param::AbstractVector,splitter::ModelSplitter)
    return [each_split_model(param,i) for i in splitter.comp_splitter]
end

function core_split_model(Base.@nospecialize(params),splitter::ModelSplitter)
    T = typeof(params)
    split_paramsvals = (core_split_model(getfield(params,i),splitter) for i  ∈ fieldnames(T))
    return T.(split_paramsvals...)
end

split_model(model::EoSModel,subset=nothing) = auto_split_model(model,subset)

function auto_split_model(Base.@nospecialize(model::EoSModel),subset=nothing)
    try
        if has_groups(model)
            init = groups(model)
        else
            init = components(model)
        end
        #==
        model_splitter contains:
        -comp: component keys
        -group: group keys
        -group_splitter: group indices
        -comp_splitter: component indices
        in the case of component models, group_splitter == comp_splitter
        ==#

        model_splitter,allfields = init_splitter(init,subset)
        len = length(model_splitter.comp_splitter)
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
                    allfields[modelkey]= core_split_model(modelx,model_splitter)
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
