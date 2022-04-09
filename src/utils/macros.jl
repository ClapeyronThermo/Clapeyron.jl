"""
    arbitraryparam(params)

returns the first field in the struct that is a subtype of `ClapeyronParam`. errors if it finds none.
"""
function arbitraryparam(params)
    paramstype = typeof(params)
    idx = findfirst(z->z <: ClapeyronParam,fieldtypes(paramstype))
    if isnothing(idx)
        error("The paramater struct ", paramstype, " must contain at least one ClapeyronParam")
    end
     return fieldnames(paramstype)[idx] |> z->getfield(params,z)
end



"""
    @groups

This macro is an alias to

    model.groups.i_flattenedgroups

`iflattenedgroups` is an iterator that goes through all groups in flattenedgroups.
"""
macro groups()
    return :($(esc(:(model.groups.i_flattenedgroups::UnitRange{Int64}))))
end

"""
    @groups(component)

This macro is an alias to

    model.groups.i_groups[component]

`i_groups[component]` is an iterator that goes through all groups in relevent to a given component.
"""
macro groups(component)
    return :($(esc(:(model.groups.i_groups[$(component)]::Vector{Int}))))
end

"""
    @sites(component)

This macro is an alias to

    model.sites.i_sites[component]

`i_sites[component]` is an iterator that goes through all sites relevant to
each group in a GC model, and to each main component in a non-GC model.
"""
macro sites(component)
    return :($(esc(:(model.sites.i_sites[$(component)]::Vector{Int}))))
end

"""
    @f(func,a,b,c,...)

This macro is an alias to
    
    func(model, V, T, z, a, b, c, ...)

where `func` is the name of the function, `model` is the model struct,
`V` is the volume, `T` is the absolute temperature, `z` is an array of number
of moles of each component, and `a`, `b`, `c`, ... are arbitrary parameters
that get passed to `func`.

It is very common for functions that are involved in the models to contain the
`model`, `V`, `T` and `z` parameters, so this macro helps reduce code repetition
as long as the first four parameters in the function are written exactly as above.

"""
macro f(func, args...)
    quote
        $func(model,V,T,z,$(args...))
    end |> esc
end


"""
    @newmodelgc modelname parent paramstype

This is a data type that contains all the information needed to use an EoS model.
It also functions as an identifier to ensure that the right functions are called.

The user is expected to create an outter constructor that takes this signature

    function modelname(components::Array{String,1})

It should then return name(params::paramtype, groups::GroupParam, sites::SiteParam, idealmodel::IdealModel)

= Fields =
The Struct consists of the following fields:

* components: a string lists of components
* icomponents: an iterator that goes through the indices corresponding to each component
* groups: a [`GroupParam`](@ref)
* sites: a [`SiteParam`](@ref)
* params: the Struct paramstype that contains all parameters in the model
* idealmodel: the IdealModel struct that determines which ideal model to use
* assoc_options: struct containing options for the association solver. see [`AssocOptions`](@ref)
* references: reference for this EoS

See the tutorial or browse the implementations to see how this is used.
"""
macro newmodelgc(name, parent, paramstype)
    quote 
    struct $name{T <: IdealModel} <: $parent
        components::Array{String,1}
        icomponents::UnitRange{Int}
        groups::GroupParam
        sites::SiteParam
        params::$paramstype
        idealmodel::T
        assoc_options::AssocOptions
        references::Array{String,1}
    end

    has_sites(::Type{<:$name}) = true
    has_groups(::Type{<:$name}) = true

    function Base.show(io::IO, mime::MIME"text/plain", model::$name)
        return gc_eosshow(io, mime, model)
    end

    function Base.show(io::IO, model::$name)
        return eosshow(io, model)
    end

    Base.length(model::$name) = Base.length(model.icomponents)

    molecular_weight(model::$name,z=SA[1.0]) = group_molecular_weight(model.groups,mw(model),z)

end |> esc
end

"""

    @newmodel name parent paramstype

This is exactly the same as the above but for non-GC models.
All group parameters are absent in this struct.
The sites are associated to the main component rather than the groups,
and the respective fieldnames are named correspondingly.
"""
macro newmodel(name, parent, paramstype)
    quote 
    struct $name{T <: IdealModel} <: $parent
        components::Array{String,1}
        icomponents::UnitRange{Int}
        sites::SiteParam
        params::$paramstype
        idealmodel::T
        assoc_options::AssocOptions
        references::Array{String,1}
    end
    has_sites(::Type{<:$name}) = true
   
    function Base.show(io::IO, mime::MIME"text/plain", model::$name)
        return eosshow(io, mime, model)
    end

    function Base.show(io::IO, model::$name)
        return eosshow(io, model)
    end
    molecular_weight(model::$name,z=SA[1.0]) = comp_molecular_weight(mw(model),z)
    Base.length(model::$name) = Base.length(model.icomponents)
    end |> esc
end

"""
    @newmodelsimple

Even simpler model, primarily for the ideal models.
Contains neither sites nor ideal models.
"""
macro newmodelsimple(name, parent, paramstype)
    quote 
    struct $name <: $parent
        components::Array{String,1}
        icomponents::UnitRange{Int}
        params::$paramstype
        references::Array{String,1}
    end

    function Base.show(io::IO, mime::MIME"text/plain", model::$name)
        return eosshow(io, mime, model)
    end

    function Base.show(io::IO, model::$name)
        return eosshow(io, model)
    end

    Base.length(model::$name) = Base.length(model.icomponents)

    end |> esc
end

const IDEALTYPE = Union{T,Type{T}} where T<:EoSModel

function (::Type{model})(params::EoSParam,
        groups::GroupParam,
        sites::SiteParam,
        idealmodel::IDEALTYPE = BasicIdeal;
        ideal_userlocations::Vector{String}=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    components = groups.components
    icomponents = 1:length(components)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    return model(components, icomponents,
    groups,
    sites,
    params, init_idealmodel, assoc_options, references)
end

function (::Type{model})(params::EoSParam,
        groups::GroupParam,
        idealmodel::IDEALTYPE = BasicIdeal;
        ideal_userlocations::Vector{String}=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    sites = SiteParam(groups.components)
    return model(params,groups,sites,idealmodel;ideal_userlocations,references,assoc_options,verbose)
end

#non GC
function (::Type{model})(params::EoSParam,
        sites::SiteParam,
        idealmodel::IDEALTYPE = BasicIdeal;
        ideal_userlocations::Vector{String}=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel
    
    components = sites.components
    icomponents = 1:length(components)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    return model(components, icomponents,
    sites, params, init_idealmodel, assoc_options, references)
end


#normal macro model
function (::Type{model})(params::EoSParam,
        idealmodel::IDEALTYPE;
        ideal_userlocations::Vector{String}=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    arbparam = arbitraryparam(params)
    components = arbparam.components
    sites = SiteParam(components)
    return model(params,sites,idealmodel;ideal_userlocations,references,assoc_options,verbose)
end

function (::Type{model})(params::EoSParam;
    references::Vector{String}=String[],
    verbose::Bool = false) where model <:EoSModel
    #if there isnt any params, just put empty values.
    if Base.issingletontype(typeof(params))
        components = String[]
        icomponents =1:0
    else
        arbparam = arbitraryparam(params)
        components = arbparam.components
        icomponents = 1:length(components)
    end
    return model(components,icomponents,params,references)
end

function init_model(idealmodel::EoSModel,components,userlocations,verbose)
    return idealmodel
end

function init_model(::Nothing,components,userlocations,verbose)
    return nothing
end


"""
    @registermodel(model)

given an existing model, composed of Clapeyron EoS models, ClapeyronParams or EoSParams, it will generate 
the necessary traits to make the model compatible with Clapeyron routines.

"""
macro registermodel(model)
    _model = getfield(@__MODULE__(),model)
    ∅ = :()

    _has_components = hasfield(_model,:components)
    _has_icomponents = hasfield(_model,:icomponents)
    _has_sites = hasfield(_model,:sites)
    _has_groups = hasfield(_model,:groups)
    
    _sites = _has_sites ? :(has_sites(::Type{<:$model}) = true) : ∅
    _groups = _has_groups ? :(has_groups(::Type{<:$model}) = true) : ∅

    _eos_show = 
    if _has_components
        if _has_groups
            quote
                function Base.show(io::IO, mime::MIME"text/plain", model::$model)
                    return gc_eosshow(io, mime, model)
                end
            
                function Base.show(io::IO, model::$model)
                    return gc_eosshow(io, model)
                end
            end
        else
            quote
                function Base.show(io::IO, mime::MIME"text/plain", model::$model)
                    return eosshow(io, mime, model)
                end
            
                function Base.show(io::IO, model::$model)
                    return eosshow(io, model)
                end
            end
        end
    else
        ∅
    end
  

    _length =
    if _has_icomponents
    :(Base.length(model::$model) = Base.length(model.icomponents))
    elseif _has_components
        :(Base.length(model::$model) = Base.length(model.components))
    else
        ∅
    end

    _molecular_weight = 
    if _has_components
        if _has_groups 
            :(molecular_weight(model::$model,z=SA[1.0]) =group_molecular_weight(model.groups,mw(model),z))
        else
            :(molecular_weight(model::$model,z=SA[1.0]) =comp_molecular_weight(mw(model),z))
        end
    else
        ∅
    end

return quote 
    $_eos_show
    $_sites
    $_groups
    $_length
    $_molecular_weight
    end |> esc
end

export @newmodel, @f, @newmodelgc, @newmodelsimple
