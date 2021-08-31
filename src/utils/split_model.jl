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

function split_model(param::SingleParam{T},
    splitter =split_model(1:length(param.components))) where T    
    function generator(I)
    return SingleParam(
    param.name,
    param.components[I],
    param.values[I],
    param.ismissingvalues[I],
    param.sourcecsvs,
    param.sources
    )
    end
    return [generator(i) for i in splitter]
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
    return [generator(I) for I in splitter]
end

function split_model(param::AbstractVector)
    return [[xi] for xi in param]
end 

#this conversion is lossy, as interaction between two or more components are lost.
#also, this conversion stores the site values for other components. (those are not used)
function split_model(param::AssocParam{T},
    splitter = split_model(1:length(param.components))) where T
    function generator(I)     
        _value  = param.values[I,I]
        _ismissingvalues = param.ismissingvalues[I,I]
        return AssocParam{T}(
                param.name,
                param.components[I],
                _value,
                _ismissingvalues,
                param.allcomponentsites[I],
                param.sourcecsvs,
                param.sources
                )
        end
    return [generator(I) for I in splitter]
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
    [generator(i) for i in 1:len]
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
            1:length(param.flattenedsites),
        param.sourcecsvs)
    end
    return [generator(i) for i in splitter]
end

function split_model(Base.@nospecialize(params::EoSParam),splitter)
    T = typeof(params)
    split_paramsvals = [split_model(getfield(params,i),splitter) for i  in fieldnames(T)]
    return T.(split_paramsvals...)
end

export SingleParam, SiteParam, PairParam, AssocParam, GroupParam
#

split_model(model::EoSModel) = auto_split_model(model)

function auto_split_model(Base.@nospecialize(model::EoSModel))
    try
        allfields = Dict{Symbol,Any}()
        
        if has_groups(typeof(model))
            splitter = model.groups.i_groups
        else
            splitter = split_model(1:length(model.components))
        end
        len = length(splitter)

        if hasfield(typeof(model),:groups)
            allfields[:groups] = split_model(model.groups)
        end
        M = typeof(model)

        len_comps = length(model.components)
        allfields[:components] = split_model(model.components)
        
        if hasfield(typeof(model),:icomponents)
            allfields[:icomponents] = [1:1 for _ in 1:len_comps]
        end
        
        #process all model fields
        modelfields = filter(x->getproperty(model,x) isa EoSModel,fieldnames(M))
        for modelkey in modelfields
            modelx = getproperty(model,modelkey)
            if is_splittable(modelx)
                allfields[modelkey]= split_model(modelx)
            else
                allfields[modelkey] = fill(modelx,len_comps)
            end
        end

        modelfields = filter(x->getproperty(model,x) isa Vector{<:EoSModel},fieldnames(M))
        for modelkey in modelfields
            modelx = getproperty(model,modelkey)
            allfields[modelkey] = [[tup] for tup in modelx]
        end

        #process all empty (Missing,Nothing) fields
        emptyfields = filter(x->getproperty(model,x) isa Union{Nothing,Missing},fieldnames(M))

        for emptykey in emptyfields
             allfields[emptykey] = fill(model.emptykey,len_comps)
        end
        
        
        if hasfield(typeof(model),:params)
            allfields[:params] = split_model(model.params,splitter)
        end

        if hasfield(typeof(model),:references)
            allfields[:references] = fill(model.references,len_comps)
        end

        if hasfield(typeof(model),:absolutetolerance)
            allfields[:absolutetolerance] = fill(model.absolutetolerance,len_comps)
        end
        if hasfield(typeof(model),:sites)
            allfields[:sites] = split_model(model.sites,splitter)
        end


        return [M((allfields[k][i] for k in fieldnames(M))...) for i in 1:len]
    catch e
        rethrow(e)
        @show model
        return simple_split_model(model)
    end
end

##fallback,around 50 times slower if there is any need to read csvs

function simple_split_model(Base.@nospecialize(model::EoSModel))
    MODEL = typeof(model)
    pure = Vector{MODEL}(undef,0)
    for comp ∈ model.components
        push!(pure,MODEL([comp]))
    end
    return pure
end

export split_model