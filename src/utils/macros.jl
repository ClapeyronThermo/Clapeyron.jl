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
    @comps

This macro is an alias to

    model.icomponents

The caveat is that `model` has to exist in the local namespace.
`model` is expected to be an EoSModel type that contains the `icomponents` field.
`icomponents` is an iterator that goes through all component indices.
"""
macro comps()
    return :($(esc(:(model.icomponents))))
end

"""
    @groups

This macro is an alias to

    model.groups.i_flattenedgroups

`iflattenedgroups` is an iterator that goes through all groups in flattenedgroups.
"""
macro groups()
    return :($(esc(:(model.groups.i_flattenedgroups))))
end

"""
    @groups(component)

This macro is an alias to

    model.groups.i_groups[component]

`i_groups[component]` is an iterator that goes through all groups in relevent to a given component.
"""
macro groups(component)
    return :($(esc(:(model.groups.i_groups[$(component)]))))
end

"""
    @sites(component)

This macro is an alias to

    model.sites.i_sites[component]

`i_sites[component]` is an iterator that goes through all sites relevant to
each group in a GC model, and to each main component in a non-GC model.
"""
macro sites(component)
    return :($(esc(:(model.sites.i_sites[$(component)]))))
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
    args = [esc(arg) for arg ∈ args]
    return :($(func)($(esc(:model)),$(esc(:V)),$(esc(:T)),$(esc(:z)),$(args...)))
end

"""
    @newmodelgc name parent paramstype

This is a data type that contains all the information needed to use an EoS model.
It also functions as an identifier to ensure that the right functions are called.

The user is expected to create an outter constructor that takes this signature

    function name(components::Array{String,1})

It should then return name(params::paramtype, groups::GroupParam, sites::SiteParam, idealmodel::IdealModel)

= Fields =
The Struct consists of the following fields:

* components: a string lists of components
* lengthcomponents: the number of components
* icomponents: an iterator that goes through the indices corresponding to each component

* allcomponentgroups: a list of groups for each component
* lengthallcomponentgroups: a list containing the number of groups for each component
* allcomponentngroups: a list of the group multiplicity of each group corresponding to each group in allcomponentsgroup
* igroups: an iterable that contains a list of group indices corresponding to flattenedgroups for each component

* flattenedgroups: a list of all unique groups--the parameters correspond to this list
* lengthflattenedgroups: the number of unique groups
* allcomponentnflattenedgroups: the group multiplicities corresponding to each group in flattenedgroups
* iflattenedgroups: an iterator that goes through the indices for each flattenedgroup

* allgroupsites: a list containing a list of all sites corresponding to each group in flattenedgroups
* lengthallgroupsites: a list containing the number of unique sites for each group in flattenedgroups
* allgroupnsites: a list of the site multiplicities corresponding to each group in flattenedgroups
* isites: an iterator that goes through the indices corresponding to each group in flattenedgroups

* params: the Struct paramstype that contains all parameters in the model
* idealmodel: the IdealModel struct that determines which ideal model to use
* absolutetolerance: the absolute tolerance for solvers; the default value is 1E-12
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
        absolutetolerance::Float64
        references::Array{String,1}
    end

    function $name(params::$paramstype, groups::GroupParam, sites::SiteParam, idealmodel=BasicIdeal;
                    ideal_userlocations=String[],
                    references=String[],
                    absolutetolerance=1E-12,
                    verbose=false)
        return _newmodelgc($name, params, groups, sites, idealmodel;
            ideal_userlocations=ideal_userlocations,
            references=references,
            absolutetolerance=absolutetolerance,
            verbose=verbose)
    end
    function $name(params::$paramstype, groups::GroupParam, idealmodel=BasicIdeal;
                    ideal_userlocations=String[],
                    references=String[],
                    absolutetolerance=1E-12, verbose=false)

        return _newmodelgc($name, params, groups, idealmodel;
                    ideal_userlocations=ideal_userlocations,
                    references=references,
                    absolutetolerance=absolutetolerance,
                    verbose=verbose)
    end
    

    function Base.show(io::IO, mime::MIME"text/plain", model::$name)
        return eosshow(io, mime, model)
    end

    function Base.show(io::IO, model::$name)
        return eosshow(io, model)
    end

    Base.length(model::$name) = Base.length(model.icomponents)

    molecular_weight(model::$name,z=SA[1.0]) = group_molecular_weight(model.groups,mw(model),z)

end |> esc
end

const IDEALTYPE = Type{T} where T<:IdealModel

function _newmodelgc(eostype,params, groups::GroupParam, sites::SiteParam, idealmodel::IDEALTYPE=BasicIdeal;
                    ideal_userlocations=String[],
                    references::Vector{String}=String[],
                    absolutetolerance::Float64=1E-12,
                    verbose::Bool=false)

    components = groups.components
    lengthcomponents = length(components)
    icomponents = 1:lengthcomponents

    verbose && idealmodel != BasicIdeal && @info("Now creating ideal model ", idealmodel, ".")
    idealmodelstruct = idealmodel(components; userlocations=ideal_userlocations, verbose=verbose)

return eostype{idealmodel}(components, icomponents,
                        groups,
                        sites,
       params, idealmodelstruct, absolutetolerance, references)
end


function _newmodelgc(eostype,params, groups, idealmodel=BasicIdeal;
                    references=String[],
                    absolutetolerance=1E-12,
                    verbose=false)

    components = groups.components
    sites = SiteParam(components)
    return _newmodelgc(eostype,params, groups, sites, idealmodel; references=references, absolutetolerance=absolutetolerance, verbose=verbose)
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
        absolutetolerance::Float64
        references::Array{String,1}
    end
    
    function $name(params::$paramstype, sites::SiteParam, idealmodel=BasicIdeal;
                    ideal_userlocations=String[],
                    references=String[],
                    absolutetolerance=1E-12, verbose=false) 

        return _newmodel($name, params, sites, idealmodel;
            ideal_userlocations=ideal_userlocations,
            references=references,
            absolutetolerance=absolutetolerance,
            verbose=verbose)
    end
    function $name(params::$paramstype, idealmodel=BasicIdeal;
                    ideal_userlocations=String[],
                    references=String[],
                    absolutetolerance=1E-12,
                    verbose=false)

        return _newmodel($name, params, idealmodel;
            ideal_userlocations=ideal_userlocations,
            references=references,
            absolutetolerance=absolutetolerance)
    end

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


function _newmodel(eostype, params, sites::SiteParam, idealmodel::IDEALTYPE=BasicIdeal;
                ideal_userlocations=String[],
                references::Array{String,1}=String[],
                absolutetolerance::Float64=1E-12,
                verbose::Bool=false)

    arbparam = arbitraryparam(params)
    components = arbparam.components
    icomponents = 1:length(components)
    verbose && idealmodel != BasicIdeal && @info("Now creating ideal model ", idealmodel, ".")
    idealmodelstruct = idealmodel(components; userlocations=ideal_userlocations, verbose=verbose)

    return eostype{idealmodel}(components, icomponents,
        sites,
        params, idealmodelstruct, absolutetolerance, references)
end

function _newmodel(eostype, params, idealmodel=BasicIdeal;
                ideal_userlocations=String[],
                references=String[],
                absolutetolerance=1E-12,
                verbose=false)
    
    arbparam = arbitraryparam(params)
    components = arbparam.components
    sites = SiteParam(components)
    return _newmodel(eostype, params, sites, idealmodel; references=references, ideal_userlocations=ideal_userlocations, absolutetolerance=absolutetolerance)
end

"""
Even simpler model, primarily for the ideal models.
Contains neither sites nor ideal models.
"""
macro newmodelsimple(name, parent, paramstype)
    quote 
    struct $name <: $parent
        components::Array{String,1}
        icomponents::UnitRange{Int}

        params::$paramstype
        absolutetolerance::Float64
        references::Array{String,1}
    end

    function $name(params::$paramstype;
                references::Array{String,1}=String[],
                absolutetolerance::Float64=1E-12, verbose::Bool=false)
        
        return _newmodelsimple($name, params; references=references, absolutetolerance=absolutetolerance, verbose=verbose)
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

function _newmodelsimple(eostype,params;
                        references::Array{String,1}=String[],
                        absolutetolerance::Float64=1E-12,
                        verbose::Bool=false)

    arbparam = arbitraryparam(params)
    components = arbparam.components
    icomponents = 1:length(components)

    return eostype(components, lengthcomponents, icomponents,
        params, absolutetolerance, references)
end

export @newmodel, @f, _newmodel
