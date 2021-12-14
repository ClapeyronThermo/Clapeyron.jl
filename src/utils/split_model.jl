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

function split_model(param::SingleParam{T},
    splitter =split_model(1:length(param.components))) where T    
    function generator(I)
    return SingleParam(
    param.name,
    param.components[I],
    deepcopy(param.values[I]),#for SingleParam{Vector{X}}
    param.ismissingvalues[I],
    param.sourcecsvs,
    param.sources
    )
    end
    return [generator(i) for i ∈ splitter]
end

function split_model(param::ClapeyronParam,groups::GroupParam)
    return split_model(param,groups.i_groups)
end
#this conversion is lossy, as interaction between two or more components are lost.

function split_model(param::PairParam{T},
    splitter = split_model(1:length(param.components))) where T
    function generator(I) 
        value = param.values[I,I]
        diagvalue = view(value,diagind(value))
        ismissingvalues = param.ismissingvalues[I,I]
        components = param.components[I]
        return PairParam{T}(
                param.name,
                components,
                value,
                diagvalue,
                ismissingvalues,
                param.sourcecsvs,
                param.sources
                )
    end 
    return [generator(I) for I ∈ splitter]
end

function split_model(param::AbstractVector,splitter = ([i] for i ∈ 1:length(param)))
    return [param[i] for i ∈ splitter]
end

function split_model(param::UnitRange{Int},splitter = ([i] for i ∈ 1:length(param)))
    return [1:length(i) for i ∈ splitter]
end

function _split_model(assoc::CompressedAssocMatrix{T},I) where T
    len = length(assoc.values)
    iszero(len) && return CompressedAssocMatrix{T}()
    old_idx = assoc.outer_indices 
    idx_bool = findall(x -> (first(x) ∈ I)&(last(x) ∈ I),old_idx)
    iszero(length(idx_bool)) && return CompressedAssocMatrix{T}()
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
    return CompressedAssocMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
end

#this conversion is lossy, as interaction between two or more components are lost.
#also, this conversion stores the site values for other components. (those are not used)
function split_model(param::AssocParam{T},
    splitter = split_model(1:length(param.components))) where T
    function generator(I)     
        _value  = _split_model(param.values,I)
     
        return AssocParam{T}(
                param.name,
                param.components[I],
                _value,
                param.sites[I],
                param.sourcecsvs,
                param.sources
                )
        end
    return [generator(I) for I ∈ splitter]
end

#this param has a defined split form
function split_model(groups::GroupParam)
    len = length(groups.components)
    function generator(i)
        return GroupParam(
        [groups.components[i]],
        [groups.groups[i]],
        [groups.n_groups[i]],
        [collect(1:length(groups.n_groups[i]))],
        groups.flattenedgroups[groups.i_groups[i]],
        [groups.n_groups[i]],
        1:length(groups.n_groups[i]),
        groups.sourcecsvs)
    end
    [generator(i) for i ∈ 1:len]
end

function split_model(param::SiteParam,
    splitter = split_model(param.components))
    function generator(I)
        return SiteParam(
            param.components[I],
            param.sites[I],
            param.n_sites[I],
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
    split_paramsvals = [split_model(getfield(params,i),splitter) for i  ∈ fieldnames(T)]
    return T.(split_paramsvals...)
end

export SingleParam, SiteParam, PairParam, AssocParam, GroupParam
#

split_model(model::EoSModel,subset=nothing) = auto_split_model(model,subset)

function auto_split_model(Base.@nospecialize(model::EoSModel),subset=nothing)
    try
        allfields = Dict{Symbol,Any}()
        if has_groups(typeof(model))
            raw_splitter = model.groups.i_groups
            subset !== nothing && throw("using subsets is not supported with Group Contribution models")
        else
            raw_splitter = split_model(Vector(1:length(model.components)))
        end
        if subset === nothing
            splitter = raw_splitter
        elseif eltype(subset) <: Integer
            splitter = raw_splitter[subset]
        elseif eltype(subset) <: AbstractVector
            splitter = subset
        else
            throw("invalid type of subset.")
        end
        
        len = length(splitter)
        M = typeof(model)
        allfieldnames = fieldnames(M)
        if hasfield(typeof(model),:groups) #TODO implement a splitter that accepts a subset
            allfields[:groups] = split_model(model.groups)
            allfields[:components] = getfield.(allfields[:groups],:components)
            for modelkey in allfieldnames
                if !haskey(allfields,modelkey)
                    modelx = getproperty(model,modelkey)
                    if modelx isa Vector{<:EoSModel}
                        allfields[modelkey] = split_model(modelx) #this is an api problem
                    end
                end
            end
        end
        
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

        return [M((allfields[k][i] for k ∈ fieldnames(M))...) for i ∈ 1:len]
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