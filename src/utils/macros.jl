#const IDEALTYPE = Union{T,Type{T}} where T<:EoSModel

"""
    arbitraryparam(params)

Returns the first field in the struct that is a subtype of `ClapeyronParam`. Errors if it finds none.
"""
function arbitraryparam(params)
    paramstype = typeof(params)
    idx = findfirst(z->z <: ClapeyronParam,fieldtypes(paramstype))
    if isnothing(idx)
        error("The parameter struct ", paramstype, " must contain at least one ClapeyronParam")
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

"""
    @sum(expr)

A macro that can be used to sum over all the variables in an expression. A faster alternative to `sum(@. expr)`.

## Example
```julia
x = [1,2,3]
y = [0.1,0.2,0.3]
z = 2

x1 = @sum(x[i]+y[i]*z)
x2 = sum(@. x+y*z)
x1 ≈ x2 #true
```
"""
macro sum(expr)
    variable_names = Expr(:tuple)
    iterator = Symbol[]
    length_indicator = Symbol[]
    cache = (variable_names.args,iterator,length_indicator)
    if expr.head == :call
        args = expr.args
        for i in 2:length(args)
            __sum_add_variables(cache,args[i])
        end
    elseif expr.head == :ref
        __sum_add_variables(cache,expr)
    else

    end
    iterator = unique!(iterator)
    length(iterator) != 1 && error("@sum: only one iterator index is allowed")
    length(length_indicator) == 0 && error("@sum: no length indicator found")
    idx = iterator[1]
    len = length_indicator[1]
    res_expr = Expr(:call,:(Base.promote_eltype))
    append!(res_expr.args,variable_names.args)
    return quote
        local __sum_result__ = zero($res_expr)
            @inbounds for $idx in 1:first(size($len))
                __sum_result__ += $expr
            end
            __sum_result__
    end  |> esc
end

function __sum_add_variables(cache,expr::Symbol)
    vars,_,_ = cache
    push!(vars,expr)
end

__sum_add_variables(cache,expr::Number) = nothing

function __sum_add_variables(cache,expr::Expr)
    vars,idx,len = cache
    if expr.head == :ref #vector or array
        sym_name = expr.args[1]
        push!(len,sym_name)
        push!(vars,sym_name)
        for j in 2:length(expr.args)
            push!(idx,expr.args[j])
        end
    elseif expr.head == :call
        for i in 2:length(expr.args)
            __sum_add_variables(cache,expr.args[i])
        end     
    end
end

"""
    default_locations(::Type{T}) where T <: EoSModel

Used for models defined via the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros.

Defines the default locations used for parsing the parameters for the input `EoSModel` type, relative to the database location.
"""
default_locations(m::EoSModel) = default_locations(parameterless_type(m))
default_locations(M) = String[]

"""
    default_gclocations(::Type{T}) where T <: EoSModel

Used for models defined via the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros.

Defines the default locations used for parsing groups for the input `EoSModel` type, relative to the database location.
"""
default_gclocations(m::EoSModel) = default_gclocations(parameterless_type(m))
default_gclocations(M) = String[]

"""
    default_ignore_missing_singleparams(::Type{T}) where T <: EoSModel

Used for models defined via the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros.

Defines the default parameters to ignore when constructing a model.
"""
default_ignore_missing_singleparams(M) = String[]

"""
    default_asymmetricparams(::Type{T}) where T <: EoSModel

Used for models defined via the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros.

Defines the default asymmetric parameters when constructing a model.
"""
default_asymmetricparams(M) = String[]

"""
    default_getparams_arguments(::Type{T},userlocations,verbose) where T <: EoSModel

Used for models defined via the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros.

Defines the `ParamsOptions` object that is passed as arguments to `getparams`, when building the input `EoSModel`.
"""
default_getparams_arguments(M,userlocations,verbose) = ParamOptions(;verbose,userlocations, ignore_missing_singleparams=default_ignore_missing_singleparams(M), asymmetricparams=default_asymmetricparams(M))

"""
    transform_params(::Type{T},params) where T <: EoSModel
    transform_params(::Type{T},params,components_or_groups) where T <: EoSModel
    transform_params(::Type{T},params,components_or_groups,verbose) where T <: EoSModel

Used for models defined via the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros.

Given a collection of params, with `(keytype(params)) isa String`, returns a modified collection with all the parameters necessary to build the `params` field contained in the `EoSModel`.

You can overload the 2, 3 or 4-argument version, depending on the need of a components vector (or `GroupParam` in a GC model), or if you want to customize the `verbose` message.

## Example

For the PC-SAFT equation of state, we perform Lorentz-Berthelot mixing of `epsilon` and `sigma`, and we scale the `sigma` parameters:
```julia
function transform_params(::Type{PCSAFT},params)
    segment = params["segment"]
        k = get(params,"k",nothing)
        params["sigma"].values .*= 1E-10
        sigma = sigma_LorentzBerthelot(params["sigma"])
        epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
        params["sigma"] = sigma
        params["epsilon"] = epsilon
        return params
    end
```
"""
function transform_params end

transform_params(M,params) = params
transform_params(M,params,components_or_groups) = transform_params(M,params)
function transform_params(M,params,components_or_groups,verbose)
    verbose && @info "generating parameters for $M"
    transform_params(M,params,components_or_groups)
end


function __struct_expr!(name,paramstype,idealmodel = true)
    if idealmodel
        struct_expr = :($name{I <: IdealModel})
    elseif !idealmodel && paramstype isa Symbol
        struct_expr = name
    else
        struct_expr = :($name{XX})
        popat!(struct_expr.args,length(struct_expr.args))
    end
    if paramstype isa Expr && paramstype.head == :curly
        append!(struct_expr.args,paramstype.args[2:end])
        curly_args = paramstype.args
        for i in 1:length(curly_args)
            if curly_args[i] isa Expr && curly_args[i].head == :<:
                curly_args[i] = curly_args[i].args[1]
            end
        end
    end
    return struct_expr
end

"""
    @newmodelgc modelname parent paramstype [sitemodel = true, use_struct_param = false]

This is a data type that contains all the information needed to use an EoS model.
It also functions as an identifier to ensure that the right functions are called.

You can pass an optional 4th `Bool` argument  indicating if you want to use sites with this model or not. Defaults to `true`.

You can also pass another optional 5th `Bool` argument indicating if a second order GroupParam (`StructGroupParam`) is used or not. Defaults to `false`

= Fields =
The Struct consists of the following fields:

* components: a string lists of components
* groups: a [`GroupParam`](@ref)
* sites: a [`SiteParam`](@ref) (optional)
* params: the Struct paramstype that contains all parameters in the model
* idealmodel: the IdealModel struct that determines which ideal model to use
* assoc_options: struct containing options for the association solver. see [`AssocOptions`](@ref)
* references: reference for this EoS

See the tutorial or browse the implementations to see how this is used.
"""
macro newmodelgc(name, parent, paramstype,sitemodel = true,use_struct_param = false)
    
    if use_struct_param
        grouptype = :StructGroupParam
    else
        grouptype = :GroupParam
    end

    struct_expr = __struct_expr!(name,paramstype)
    res = if sitemodel
        quote
            struct $struct_expr <: $parent
                components::Array{String,1}
                groups::Clapeyron.$grouptype
                sites::Clapeyron.SiteParam
                params::$paramstype
                idealmodel::I
                assoc_options::Clapeyron.AssocOptions
                references::Array{String,1}
            end

            function $name(components;
                idealmodel = Clapeyron.BasicIdeal,
                userlocations = String[],
                group_userlocations = String[],
                ideal_userlocations = String[],
                assoc_options = Clapeyron.default_assoc_options($name),
                reference_state = nothing,
                verbose = false)

                Clapeyron.build_eosmodel($name,components,idealmodel,userlocations,group_userlocations,ideal_userlocations,verbose,assoc_options,reference_state)
            end
        end
    else
        quote
            struct $struct_expr <: $parent
                components::Array{String,1}
                groups::Clapeyron.$grouptype
                params::$paramstype
                idealmodel::I
                references::Array{String,1}
            end

            function $name(components;
                idealmodel = Clapeyron.BasicIdeal,
                userlocations = String[],
                group_userlocations = String[],
                ideal_userlocations = String[],
                reference_state = nothing,
                verbose = false)

                Clapeyron.build_eosmodel($name,components,idealmodel,userlocations,group_userlocations,ideal_userlocations,verbose,nothing,reference_state)
            end
        end
    end

    return res |> esc
end

"""

    @newmodel name parent paramstype [sitemodel = true]

This is exactly the same as the above but for non-GC models.
All group parameters are absent in this struct.
The sites are associated to the main component rather than the groups,
and the respective fieldnames are named correspondingly.

You can pass an optional bool indicating if you want to use sites with this model or not. Defaults to `true`.

## Example
```julia
struct MySAFTParam
    a::SingleParam{Float64}
    b::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

@newmodel MySAFT SAFTModel MySAFTParam #defines a model, with association sites

struct MyModelParam
    a::SingleParam{Float64}
    b::SingleParam{Float64}
end

@newmodel MyModel EoSModel MyModelParam false #defines a model without sites
```
"""
macro newmodel(name, parent, paramstype,sitemodel = true)
    
    struct_expr = __struct_expr!(name,paramstype)

    res = if sitemodel
        quote
            struct $struct_expr <: $parent
                components::Array{String,1}
                sites::Clapeyron.SiteParam
                params::$paramstype
                idealmodel::I
                assoc_options::Clapeyron.AssocOptions
                references::Array{String,1}
            end

            function $name(components;
                idealmodel = Clapeyron.BasicIdeal,
                userlocations = String[],
                ideal_userlocations = String[],
                assoc_options = Clapeyron.default_assoc_options($name),
                reference_state = nothing,
                verbose = false)

                Clapeyron.build_eosmodel($name,components,idealmodel,userlocations,nothing,ideal_userlocations,verbose,assoc_options,reference_state)
            end
        end
    else
        quote
            struct $struct_expr <: $parent
                components::Array{String,1}
                params::$paramstype
                idealmodel::I
                references::Array{String,1}
            end

            function $name(components;
                idealmodel = Clapeyron.BasicIdeal,
                userlocations = String[],
                ideal_userlocations = String[],
                reference_state = nothing,
                verbose = false)

                Clapeyron.build_eosmodel($name,components,idealmodel,userlocations,nothing,ideal_userlocations,verbose,nothing,reference_state)
            end
        end
    end

    return res |> esc
end

"""
    @newmodelsimple name parent paramstype

Even simpler model, primarily for the ideal models.
Contains neither sites nor ideal models.
"""
macro newmodelsimple(name, parent, paramstype)
    struct_expr = __struct_expr!(name,paramstype,false)
    return quote
        struct $struct_expr <: $parent
            components::Array{String,1}
            params::$paramstype
            references::Array{String,1}
        end

        function $name(components;userlocations = String[],reference_state = nothing,verbose = false)
            Clapeyron.build_eosmodel($name,components,nothing,userlocations,nothing,nothing,verbose,nothing,reference_state)
        end
    end |> esc
end

"""
    @newmodelsingleton name parent

A macro that defines an EoSModel without any fields (\"singleton\" struct.). Useful for defining EoS that don't use any parameters, while being composable with other `EoSModels`.
"""
macro newmodelsingleton(name,parent)
    quote
    struct $name <: $parent end
    Clapeyron.is_splittable(::$name) = false
    function $name(components;userlocations = String[],verbose = false,reference_state = nothing)
        reference_state_checkempty($name,reference_state)
        return $name()
    end
    end |> esc
end

"""
    init_model(model::EoSModel,components,userlocations = String[],verbose = false)
    init_model(::Type{𝕄},components,userlocations = String[],verbose = false) where 𝕄 <: EoSModel

Utility for building simple models. If a model instance is passed, it will return that instance.
Otherwise, it will build the model from the input components and user locations.

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
function init_model(model::EoSModel,components,userlocations = String[],verbose = false,reference_state = nothing)
    return model
end

function init_model(::Nothing,components,userlocations = String[],verbose = false,reference_state = nothing)
    return nothing
end

function init_model(::Type{𝕄},components,userlocations = String[],verbose = false,reference_state = nothing) where  𝕄 <: EoSModel
    if verbose
        @info "Building an instance of $(info_color(string(𝕄))) with components $components"
    end
    if has_reference_state(𝕄)
        return 𝕄(components;userlocations,verbose,reference_state)
    else
        return 𝕄(components;userlocations,verbose)
    end
end

function init_model(f::Function,components,userlocations = String[],verbose = false,reference_state = nothing)
    if verbose
        @info "building an EoS model, using function $(info_color(string(f))) with components $components"
    end
    return f(components;userlocations,verbose,reference_state)
end

macro initmodel(modelexpr)
    model = initmodel_macro(modelexpr)
    return quote
        if $model isa EoSModel || $model === nothing
            $model
        else
            $modelexpr
        end
    end |> esc
end

function initmodel_macro(expr::Expr)
    if expr.head == :call
        return expr.args[1]
    else
        throw(error("invalid argument to @initmodel macro"))
    end
end

"""
    @registermodel(model)

Given an existing model, composed of Clapeyron EoS models, ClapeyronParams or EoSParams, it will generate
the necessary traits to make the model compatible with Clapeyron routines.

!!! info
    This macro is a no-op from Clapeyron 0.5 onwards.
"""
macro registermodel(model)
    esc(model)
end

function build_eosmodel(::Type{M},components,idealmodel,userlocations,group_userlocations,ideal_userlocations,verbose,assoc_options = nothing,reference_state = nothing) where M <: EoSModel

    paramtype = fieldtype(M,:params)
    _components = format_components(components)

    #non-splittable
    if paramtype === Nothing
        return M()
    #we don't need to parse params.
    elseif Base.issingletontype(paramtype)
        return M(_components,paramtype(),default_references(M))
    end
    #all fields of the model.
    result = Dict{Symbol,Any}()
    result[:components] = _components
    #parse params from database.
    options = default_getparams_arguments(M,userlocations,verbose)
    if has_groups(M)
        groups = GroupParam(format_gccomponents(components),default_gclocations(M);group_userlocations,verbose)
        params_in = getparams(groups, default_locations(M),options)
        result[:groups] = groups
    else
        groups = nothing
        params_in = getparams(_components, default_locations(M),options)
    end

    #inject reference state if not built
    if has_reference_state(M)
            params_in["reference_state"] = __init_reference_state_kw(reference_state)
    else
        #this could fail when the type is not entirely known.
        #reference_state_checkempty(M,reference_state)
    end

    #put AssocOptions inside params, so it can be used in transform_params
    if has_sites(M)
        if !haskey(params_in,"assoc_options")
            params_in["assoc_options"] = assoc_options
        else
            #throw(error("cannot overwrite \"assoc_options\" key, already exists!"))
        end

        #legacy case: the model has a SiteParam, but it does not have association parameters.
        #we just build an empty one
        if !haskey(params_in,"sites")
            #todo: check how this interact with GC, but i suspect that with our new Approach
            #we always want component-based sites
            params_in["sites"] = SiteParam(_components)
        end
    end

    #perform any transformations, pass components or groups
    if has_groups(M)
        params_out = transform_params(M,params_in,groups,verbose)
    else
        params_out = transform_params(M,params_in,_components,verbose)
    end

    #mix sites
    if has_sites(M)
        assoc_mix!(params_out,_components)
    end
    #build EoSParam
    pkgparam = build_eosparam(paramtype,params_out)
    result[:params] = pkgparam

    #build SiteParam, if needed
    if has_sites(M)
        _sites = get(params_out,"sites",nothing)
        if isnothing(_sites)
            _sites = SiteParam(_components)
        end
        result[:sites] = _sites
        result[:assoc_options] = assoc_options
    end

    #add references, if needed
    if hasfield(M,:references)
        result[:references] = default_references(M)
    end

    #build idealmodel, if needed
    
    if hasfield(M,:idealmodel)
        if has_reference_state(idealmodel)
            #=
            we want to execute set_reference_state!(model) only once (ideal models don't have)
            saturation information so some standard states cannot be initialized.
            
            To avoid this, we set the input reference state to :no_set, and then we reset it to the 
            input value. with this strategy, we can differenciate between standalone ideal models and
            ideal models stored inside a residual model.
            =#
            input_reference_state = __init_reference_state_kw(reference_state)
            std_type = input_reference_state.std_type
            input_reference_state.std_type = :no_set
            init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose,input_reference_state)
            input_reference_state.std_type = std_type
            idmodel_reference_state = Clapeyron.reference_state(init_idealmodel)
            idmodel_reference_state.std_type = std_type
        else
            init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
        end
        result[:idealmodel] = init_idealmodel
    end

    #build model
    model = M((result[k] for k in fieldnames(M))...)
    #fit reference state
    set_reference_state!(model,verbose = verbose)
    return model
end

export @newmodel, @f, @newmodelgc, @newmodelsimple, @newmodelsingleton
