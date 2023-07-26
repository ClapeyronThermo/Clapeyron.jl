#const IDEALTYPE = Union{T,Type{T}} where T<:EoSModel

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

    1:length(model.groups.flattenedgroups)

"""
macro groups()
    return quote
            (1:length(model.groups.flattenedgroups))::UnitRange{Int64}
    end |> esc
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

default_references(M) = String[]
default_locations(M) = String[]
default_getparams_arguments(model,userlocations,verbose) = ParamOptions(;verbose,userlocations)
init_groups(M,components,userlocations,verbose) = GroupParam(components,String[];userlocations,verbose)
transform_params(M,params) = params
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
    struct $name{T <: Clapeyron.IdealModel} <: $parent
        components::Array{String,1}
        groups::Clapeyron.GroupParam
        sites::Clapeyron.SiteParam
        params::$paramstype
        idealmodel::T
        assoc_options::Clapeyron.AssocOptions
        references::Array{String,1}
    end

    function $name(params::$paramstype,
        groups::Clapeyron.GroupParam,
        idealmodel = Clapeyron.BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::Clapeyron.AssocOptions = Clapeyron.AssocOptions(),
        verbose::Bool = false)

        return Clapeyron.build_model($name,params,groups,idealmodel;ideal_userlocations,references,assoc_options,verbose)
    end

    function $name(params::$paramstype,
        groups::Clapeyron.GroupParam,
        sites::Clapeyron.SiteParam,
        idealmodel = BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::Clapeyron.AssocOptions = Clapeyron.AssocOptions(),
        verbose::Bool = false)

        return Clapeyron.build_model($name,params,groups,sites,idealmodel;ideal_userlocations,references,assoc_options,verbose)
    end

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
        sites::Clapeyron.SiteParam
        params::$paramstype
        idealmodel::T
        assoc_options::Clapeyron.AssocOptions
        references::Array{String,1}
    end

    function $name(params::$paramstype,
        sites::Clapeyron.SiteParam,
        idealmodel = Clapeyron.BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::Clapeyron.AssocOptions = Clapeyron.AssocOptions(),
        verbose::Bool = false)

        return Clapeyron.build_model($name,params,sites,idealmodel;ideal_userlocations,references,assoc_options,verbose)
    end

    function $name(params::$paramstype,
        idealmodel = Clapeyron.BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::Clapeyron.AssocOptions = Clapeyron.AssocOptions(),
        verbose::Bool = false)

        return Clapeyron.build_model($name,params,idealmodel;ideal_userlocations,references,assoc_options,verbose)
    end

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
        params::$paramstype
        references::Array{String,1}
    end

    function $name(params::$paramstype;
            references::Vector{String}=String[],
            verbose::Bool = false)

        return Clapeyron.build_model($name,params;references,verbose)
    end
    end |> esc
end

function build_model(::Type{model},params::EoSParam,
        groups::GroupParam,
        sites::SiteParam,
        idealmodel = BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    components = groups.components
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    return model(components,
    groups,
    sites,
    params, init_idealmodel, assoc_options, references)
end

function build_model(::Type{model},params::EoSParam,
        groups::GroupParam,
        idealmodel = BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    sites = SiteParam(groups.components)
    return model(params,groups,sites,idealmodel;ideal_userlocations,references,assoc_options,verbose)
end

#non GC
function build_model(::Type{model},params::EoSParam,
        sites::SiteParam,
        idealmodel = BasicIdeal;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    components = sites.components

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    return model(components,
    sites, params, init_idealmodel, assoc_options, references)
end


#normal macro model
function build_model(::Type{model},params::EoSParam,
        idealmodel;
        ideal_userlocations=String[],
        references::Vector{String}=String[],
        assoc_options::AssocOptions = AssocOptions(),
        verbose::Bool = false) where model <:EoSModel

    arbparam = arbitraryparam(params)
    components = arbparam.components
    sites = SiteParam(components)
    return model(params,sites,idealmodel;ideal_userlocations,references,assoc_options,verbose)
end

function build_model(::Type{model},params::EoSParam;
    references::Vector{String}=String[],
    verbose::Bool = false) where model <:EoSModel
    #if there isnt any params, just put empty values.
    if Base.issingletontype(typeof(params))
        components = String[]
    else
        arbparam = arbitraryparam(params)
        components = arbparam.components
    end
    return model(components,params,references)
end

"""
    init_model(model::EoSModel,components,userlocations=String[],verbose = false)
    init_model(::Type{𝕄},components,userlocations=String[],verbose = false) where 𝕄 <: EoSModel

Utility for building simple models. if a model instance is passed, it will return that instance.
otherwise, it will build the model from the input components and user locations.

It is normally used for models that don't have additional submodels (like ideal models)
or when such submodels are not used at all (like the pure model part of an Activity model when used in an Advanced mixing rule Cubic model)

## Example

julia> Clapeyron.init_model(MonomerIdeal,["methane","ethane"])
MonomerIdeal with 2 components:
 "methane"
 "ethane"
Contains parameters: Mw

```julia-repl
julia> model = Clapeyron.init_model(MonomerIdeal,["methane","ethane"])
MonomerIdeal with 2 components:
 "methane"
 "ethane"
Contains parameters: Mw

julia> model.params.Mw[1] = 1000
1000

julia> model2 = Clapeyron.init_model(model,["methane","ethane"])
MonomerIdeal with 2 components:
 "methane"
 "ethane"
Contains parameters: Mw

julia> model2.params.Mw
SingleParam{Float64}("Mw") with 2 components:
 "methane" => 1000.0
 "ethane" => 30.07
```

"""
function init_model(model::EoSModel,components,userlocations=String[],verbose = false)
    return model
end

function init_model(::Nothing,components,userlocations=String[],verbose = false)
    return nothing
end

function init_model(::Type{𝕄},components,userlocations=String[],verbose = false) where  𝕄 <: EoSModel
    if verbose
        @info "Building an instance of $(info_color(string(𝕄))) with components $components"
    end
    return 𝕄(components;userlocations,verbose)
end

function init_model(f::Function,components,userlocations=String[],verbose = false)
    if verbose
        @info "building an EoS model, using function $(info_color(string(f))) with components $components"
    end
    return f(components;userlocations,verbose)
end
"""
    @registermodel(model)

given an existing model, composed of Clapeyron EoS models, ClapeyronParams or EoSParams, it will generate
the necessary traits to make the model compatible with Clapeyron routines.

"""
macro registermodel(model)
    esc(model)
end

function build_eosmodel(::Type{M},components,idealmodel,userlocations,group_userlocations,ideal_userlocations,verbose,assoc_options = nothing) where T <: EoSModel

    result = Dict{Symbol,Any}()

    options = default_getparams_arguments(T,userlocations,verbose)
    if has_groups(M)
        groups =  init_groups(M,components,group_userlocations,verbose)
        params_in = getparams(groups, default_locations(M),options)
        result[:groups] = groups
        result[:components] = groups.components
    else
        groups = nothing
        params_in = getparams(components, default_locations(M),options)
        result[:components] = components
    end

    params_out = transform_params(M,params_in)
    if has_sites(M)
        assoc_mix!(params_out,assoc_options)
    end

    pkgparam = build_eosparam(M,params_out)
    result[:params] = pkgparam

    if has_sites(M)
        _sites = get(params_out,"sites",nothing)
        if isnothing(_sites)
            _sites = SiteParam(components)
        end
        result[:sites] = _sites
        result[:assoc_options] = assoc_options
    end

    if hasfield(M,:references)
        references = default_references(M)
        result[:references] = references
    end

    if hasfield(M,:idealmodel)
        init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
        result[:idealmodel] = init_idealmodel
    end

    return M((result[k] for k in fieldnames(M))...)
end

export @newmodel, @f, @newmodelgc, @newmodelsimple
#=
function __newmodel(name, parent, paramstype,sites,idealmodel)

    if sites
    struct $name{T <: IdealModel} <: $parent
        components::Array{String,1}
        groups::GroupParam
        sites::SiteParam
        params::$paramstype
        idealmodel::T
        assoc_options::AssocOptions
        references::Array{String,1}
    end

end =#
