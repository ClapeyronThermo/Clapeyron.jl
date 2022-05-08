"""
    ModelMapping(source, target, transformation)

Part of `ModelOptions`. Used to define how source from csv is mapped to the target struct. The `source` and `target` should be defined in `inputparams` and `params` respecfully in `ModelOptions`.

## Fields
- `source::Union{String,Vector{String}}`: The source headers from csv.
- `target::Union{Symbol,Vector{Symbol}}`: The target parameters.
- `transformation::function = identity`: The transformation from source to target.
- `self_in_args::Bool = false`: If true, the first argument in the transformation function is the model object itself. Note that it is inadvisable to access other params directly via getproperties on the model object as there is no guarantees that the values are up to date at the time of update. The intension is primarily for dispatch. If it is dependent on parameters, explicitly have them as arguments to the transformation function, or ensure that the parameters accessed are never changed.
"""
struct ModelMapping
    source::Vector{Symbol}
    target::Vector{Symbol}
    transformation::Function
    self_in_args::Bool
end
function ModelMapping(
        source::Vector{Symbol},
        target::Vector{Symbol},
        transformation::Function = identity;
        self_in_args::Bool = false
    )
    return ModelMapping(
        source,
        target,
        transformation,
        self_in_args
    )
end


"""
    ModelMember(name, default_type; separate_namespace = false)

Part of `ModelOptions`. This is used to specify the options for member models. Note that only swappable members need to be specified here. If member models contain member models with the same name as one specified here, the same instance will be used. If this is undesired behaviour, explicitly assign the relevant parents in the `restrictparents` field.

## Fields
- `name::Symbol`: The name of this member. It will be used as the fieldname for the created model object.
- `default_type::Type`: The default type for this member. Has to have a constructor with the same function signature.
- `split::Bool = false`: If `true`, create a vector of pure models of this type.
- `separate_namespace::Bool = true`: This is for namespace resolution of input parameter names from user-provided csv in `userlocations`. If `true`, all headers in the input csvs specified in `userlocations` should be prefixed with `{name}__`. Note that these prefixes do not nest (the `name` must be distinct per model, and so there need only be one level).
- `overwritelocations::Union{Vector{String},Nothing} = nothing`: If the model wants to overwrite the array of locations from a location in the Clapeyron database specific for this model, they may do so here. This location is not passed down to subsequent member models.
- `restrictparents::Union{Symbol,Nothing} = nothing`: This is for namespoce resolution of member models. If `nothing`, any model with a member with this name will be assigned the same instantiated object. When ambiguity arises, say if two distinct `activity` models are used, provide a vector of types to specify where this member is applicable, making use of the `nameinparent` field as well for correct asignment.
- `nameinparrent::Symbol = nothing`: If `nothing`, just take the given name. This is for when the name in member model is different from the name in current model, which allows multiple member models of the same type to be present in the main model.

## Separate namespace
To account for the possibility that parameters may share the same names across models, we can enable this option to resolve conflicts. For example, if we had a model

    A{B}

and both `A` and `B` contain distinct parameters `m`. Then we would enable `separate_namespace` for member `B`, and in the parameters file, if there is a column called `B__m`, it will take precedence over the column `m`. If column `B__m` is not present, it will find that column `m` is already initialised when getting parameters for `A`, so it will just take the same reference to that parameter.

## Restrict parents
The order of initialisation is as provided in the `members` field vector. For example, if we have model `A` has `Y` and `X` as member models (note the order) and model `X` is a member of model `Y`. If we specified it as

    A{Y{X},X}

then without any other instructions, it will first initialise `Y` with the default `X`, then intialise a second `X` for the entry in `A`.

If we instead had it written as

    A{X,Y{X}}

Unless otherwise specified, the instance `X` that resides in `A` is the same instance as the one that resides in `Y`.

If we specify `[:A]` in `restrictparents` for `ModelMember` of `X`, then a separate instance of a `X` (the default) will be constructed for `Y`.

The base model `A` will always contain an instance of any member specified in `members`. So in this case, specifying `[:Y]` for the `restrictparents` field here will have the same effect as `nothing` (both referencing the same `X`). In effect, this means that for parameter estimation, we will never need to look into modifying any params deeper than one level of members.

If we wish the `X` in `A` be distinct from the `X` in `Y`, yet have them both modular (swappable from the main constructor interface), we can give them separate names in the base model, but change the `nameinparent` to reflect the correct name. In this case, we will have

    A{X, X2, Y{X}}

where for `X`, we have `restrictparents` as `[:A]`, whereas for `X2`, we have `restrictparents` as `[:Y]` (superfluous here, as per discussion above, but provides clarification on usage), and `nameinparent` as `X`. Then `X2` here is the same instance as the `X` in `Y`. Note that the order of the `restrictparents` can be swapped, so thet he `X` in `Y` is the same instance as the `X` in `A`, instead of `X2` as before. However, in this case, both `restrictparams` have to be specified.
"""
struct ModelMember
    name::Symbol
    default_type::Type
    split::Bool
    separate_namespace::Bool
    overwritelocations::Union{Vector{String},Nothing}
    restrictparents::Union{Vector{Symbol},Nothing}
    nameinparent::Symbol
end
function ModelMember(
        name::Symbol,
        default_type::Type;
        split::Bool = false,
        separate_namespace::Bool = false,
        overwritelocations::Union{Vector{String},Nothing} = nothing,
        restrictparents::Union{Vector{Symbol},Nothing} = nothing,
        nameinparent::Union{Symbol,Nothing} = nothing
    )
    if isnothing(nameinparent)
        nameinparent = name
    end
    return ModelMember(
        name,
        default_type,
        split,
        separate_namespace,
        overwritelocations,
        restrictparents,
        nameinparent
    )
end

"""
    ParamField(name, type)

Part of `ModelOptions`. This is used to specify the relevant parameters and their concrete types in the `inputparams` and `params`. The corresponding structs will be generated to this specification.

## Fields
- `name::Symbol`: The name of the parameter field.
- `type::DataType`: Should by either SingleParam{T}, PairParam{T} or AssocParam{T}.
"""
struct ParamField
    name::Symbol
    type::DataType
end

"""
    ModelOptions(args...)

A complete definition for how the model object will be created in Clapeyron. If the `parent` field is given, it will take the field values of the parent unless explicitly overwritten.

## Fields
- `name::Symbol`: The name of the model. A struct with this name will be generated in to the global/module namespace.
- `supertype::DataType = EoSModel`: An abstract base type to be a subtype of.
- `parent::Union{ModelOptions,Nothing} = nothing`: The `ModelOptions` of the parent model.
- `members::Vector{ModelMember} = ModelMember[]`: A list of ModelMembers. Note that the order reflects the order of initialisation. Read the docs for `ModelMember` for more details.
- `locations::Vector{String} = String[]`: Default locations in Clapeyron database to look for parameters. Note that for user-specified parameters, it is easier to use the `user_locations` parameter.
- `inputparams::Vector{ParamField} = ParamField[]`: A list of relevant source parameters. The model constructor will extract the String headers according to the Symbol name given.
- `params::Vector{ParamField} = ParamField[]`: A list of relevant target parameters.
- `mappings::Vector{ModelMapping} = ModelMapping[]`: Mappings from source to target.
- `has_params::Bool = true`: Whether this model has params. A simple example of this not being the case is `BasicIdeal`.
- `has_components::Bool = true`: Whether this model is dependent on components. A simple example of this not being the case is `BasicIdeal`.
- `has_sites::Bool = false`: Whether this model has association.
- `has_groups::Bool = false`: Whether this model contains groups.
- `param_options::ParamOptions = nothing`: The default ParamOptions. If `has_params` is `true`, but param_options is `nothing`, use create one using the default constructor.
- `assoc_options::AssocOptions = nothing`: The default AssocOptions. If `has_sites` is `true`, but assoc_options is `nothing`, use create one using the default constructor.
- `references::Vector{String} = String[]`: References for this model. Usually DOIs.
- `inputparamstype::Symbol = nothing`: A struct with this name will be generated in the global/module namespace for the input params. If given `nothing`, the constructor will fill it in with `{name}InputParam`. Note that if a struct with the same name is already defined, it will not redefine it.
- `paramstype::Symbol = nothing`: A struct with this name will be generated in the global/module namespace for the target params. If given `nothing`, the constructor will fill it in with `{name}Param`. Note that if a struct with the same name is already defined, it will not redefine it.
"""

struct ModelOptions
    name::Symbol
    supertype::DataType
    parent::Union{ModelOptions,Nothing}
    members::Vector{ModelMember}
    locations::Vector{String}
    inputparams::Vector{ParamField}
    params::Vector{ParamField}
    mappings::Vector{ModelMapping}
    has_params::Bool
    has_components::Bool
    has_sites::Bool
    has_groups::Bool
    param_options::Union{ParamOptions,Nothing}
    assoc_options::Union{AssocOptions,Nothing}
    references::Vector{String}
    inputparamstype::Symbol
    paramstype::Symbol
end

function ModelOptions(
        name::Symbol;
        supertype::Union{DataType,Nothing} = nothing,
        parent::Union{ModelOptions,Nothing} = nothing,
        members::Union{Vector{ModelMember},Nothing} = nothing,
        locations::Union{Vector{String},Nothing} = nothing,
        inputparams::Union{Vector{ParamField},Nothing} = nothing,
        params::Union{Vector{ParamField},Nothing} = nothing,
        mappings::Union{Vector{ModelMapping},Nothing} = nothing,
        has_params::Union{Bool,Nothing} = nothing,
        has_components::Union{Bool,Nothing} = nothing,
        has_sites::Union{Bool,Nothing} = nothing,
        has_groups::Union{Bool,Nothing} = nothing,
        param_options::Union{ParamOptions,Nothing} = nothing,
        assoc_options::Union{AssocOptions,Nothing} = nothing,
        references::Union{Vector{String},Nothing} = nothing,
        inputparamstype::Union{Symbol,Nothing} = nothing,
        paramstype::Union{Symbol,Nothing} = nothing
    )
    if !isnothing(parent)
        if isnothing(supertype)
            supertype = parent.supertype
        end
        if isnothing(members)
            members = parent.members
        end
        if isnothing(locations)
            locations = parent.locations
        end
        if isnothing(inputparams)
            if isnothing(inputparamstype)
                inputparamstype = parent.inputparamstype
            end
            inputparams = parent.inputparams
        end
        if isnothing(params)
            if isnothing(paramstype)
                paramstype = parent.paramstype
            end
            params = parent.params
        end
        if isnothing(mappings)
            mappings = parent.mappings
        end
        if isnothing(has_params)
            has_params = parent.has_params
        end
        if isnothing(has_components)
            has_components = parent.has_components
        end
        if isnothing(has_sites)
            has_sites = parent.has_sites
        end
        if isnothing(has_groups)
            has_groups = parent.has_groups
        end
        if isnothing(param_options)
            param_options = parent.param_options
        end
        if isnothing(assoc_options)
            assoc_options = parent.assoc_options
        end
        if isnothing(references)
            references = parent.references
        end
    else
        # Default values if no parent is given.
        if isnothing(supertype)
            supertype = EoSModel
        end
        if isnothing(members)
            members = ModelMember[]
        end
        if isnothing(locations)
            locations = String[]
        end
        if isnothing(inputparams)
            inputparams = ParamField[]
        end
        if isnothing(params)
            params = ParamField[]
        end
        if isnothing(mappings)
            mappings = ModelMapping[]
        end
        if isnothing(has_params)
            has_params = true
        end
        if isnothing(has_components)
            has_components = true
        end
        if isnothing(has_sites)
            has_sites = false
        end
        if isnothing(has_groups)
            has_groups = false
        end
        if isnothing(references)
            references = String[]
        end
    end

    return ModelOptions(
        name,
        supertype,
        parent,
        members,
        locations,
        inputparams,
        params,
        mappings,
        has_params,
        has_components,
        has_sites,
        has_groups,
        has_params && isnothing(param_options) ? DefaultParamOptions : param_options,
        has_sites && isnothing(assoc_options) ? DefaultAssocOptions : assoc_options,
        references,
        isnothing(inputparamstype) ? Symbol(String(name) * "InputParam") : inputparamstype,
        isnothing(paramstype) ? Symbol(String(name) * "Param") : paramstype,
    )
end

##### Expr generators #####

function _generatecode_param_struct(
        name::Symbol,
        paramfields::Vector{ParamField}
    )::Expr
    block = Expr(:block)
    for paramfield ∈ paramfields
        push!(block.args, :($(paramfield.name)::$(paramfield.type)))
    end
    return Expr(:struct, false, :($name <: EoSParam), block)
end

function _generatecode_model_struct(modeloptions::ModelOptions)::Expr
    if isempty(modeloptions.members)
        defheader = modeloptions.name
    else
        defheader = Expr(:curly, modeloptions.name)
        for i ∈ 1:length(modeloptions.members)
            push!(defheader.args, Symbol(:M, Symbol(i)))
        end
    end

    block = Expr(:block)
    push!(block.args, :(components::Vector{String}))
    if modeloptions.has_groups
        push!(block.args, :(groups::GroupParam))
    end
    if modeloptions.has_sites
        push!(block.args, :(sites::SiteParam))
    end
    if modeloptions.has_params
        push!(block.args, :(inputparams::$(modeloptions.inputparamstype)))
        push!(block.args, :(params::$(modeloptions.paramstype)))
        push!(block.args, :(mappings::Vector{ModelMapping}))
    end
    for (i, member) ∈ enumerate(modeloptions.members)
        if member.split
            push!(block.args, Expr(:(::), member.name, Expr(:curly, :EoSVectorParam, Symbol("M" * string(i)))))
        else
            push!(block.args, Expr(:(::), member.name, Symbol("M" * string(i))))
        end
    end
    if modeloptions.has_sites
        push!(block.args, :(assoc_options::AssocOptions))
    end
    push!(block.args, :(references::Vector{String}))

    return Expr(:struct, false, Expr(:<:, defheader, modeloptions.supertype), block)
end

"""
Generates a model constructor with the folowing function signature:

## Arguments
- `components::Union{String, Vector{String}}`: The components for this model.
_ `userlocations::Vector{String} = String[]`: An array of filepaths to csvs specified by the user.
- `verbose::Bool = false`: Print more information.
- `param_options::ParamOptions`: This is present for models with parameters to change the behaviour of `getparams`.
- `assoc_options::AssocOptions`: This is present for models with sites to change the behaviour of association.
- various member models: For models with modular components, each member model can be swapped out here. For example, the `idealmodal` or `activity`.

## System arguments
These are used by Clayeyron when initialising member models. They should not be modified unless you know what you are doing.
- `_overwritelocations::Union{Vector{String}} = nothing`: If not `nothing`, overwrite the `locations` for this model.
- `_initialisedmodels::Dict{Symbol, Dict{Symbol, Any}} = Dict{Symbol, Dict{Symbol, Any}}()`: If model is already initialised, just take a reference to it unless specified otherwise in the `MemberModel`.
- `_namespace::String = ""`: Give higher priority to columns prepended with `{membermodel_name}__{inputparam_name}` if present.
- `_accumulatedparams::Dict{String, ClapeyronParam} = Dict{String, ClapeyronParam}())`: If there are parameters with the same name, just point to existing reference.
"""
function _generatecode_model_constructor(
        modeloptions::ModelOptions
    )::Expr
    # Creating function args
    func_head = Expr(:call)
    push!(func_head.args, modeloptions.name)
    # Keyword args first
    parameters = Expr(:parameters)
    push!(parameters.args, Expr(:kw, :(userlocations::Vector{String}), :(String[])))
    push!(parameters.args, Expr(:kw, :(verbose::Bool), :false))
    if modeloptions.has_params
        push!(parameters.args, Expr(:kw, :(param_options::ParamOptions), :($(modeloptions.param_options))))
    end
    if modeloptions.has_sites
        push!(parameters.args, Expr(:kw, :(assoc_options::AssocOptions), :($(modeloptions.assoc_options))))
    end
    for member ∈ modeloptions.members
        push!(parameters.args, Expr(:kw, :($(member.name)), :($(member.default_type))))
    end
    push!(parameters.args, Expr(:kw, :(_overwritelocations::Union{Vector{String},Nothing}), :nothing))
    push!(parameters.args, Expr(:kw, :(_initialisedmodels::Dict{Symbol,Dict{Symbol,Any}}), :(Dict{Symbol,Dict{Symbol,Any}}())))
    push!(parameters.args, Expr(:kw, :(_namespace::String), :""))
    push!(parameters.args, Expr(:kw, :(_accumulatedparams::Dict{String,ClapeyronParam}), :(Dict{String,ClapeyronParam}())))
    push!(func_head.args, parameters)
    # Now positional args
    if modeloptions.has_components
        # `components` is mandatory
        push!(func_head.args, :(components::Union{String,Vector{String}}))
    else
        # `components` is optional
        push!(func_head.args, Expr(:kw, :(components::Union{String,Vector{String}}), :(String[])))
    end

    # Creating function body
    block = Expr(:block)
    if modeloptions.has_params
        push!(block.args, :(mappings = $(modeloptions.mappings)))
        push!(block.args, :(locations = $(modeloptions.locations)))
        push!(block.args, Expr(:if, :(!isnothing(_overwritelocations)), :(locations = _overwritelocations)))
        if modeloptions.has_sites
            push!(block.args, :((rawparams, sites) = getparams(components, locations, param_options; userlocations, verbose)))
        else
            push!(block.args, :(rawparams = getparams(components, locations, param_options; userlocations, verbose)))
        end
        if modeloptions.has_groups  # Have to figure out what to do with GC later.
            push!(block.args, :((inputparams, params) = _initparams($(modeloptions.inputparamstype), $(modeloptions.paramstype), rawparams, mappings, _namespace)))
        else
            push!(block.args, :(merge!(_accumulatedparams, merge(rawparams, _accumulatedparams))))
            push!(block.args, :((inputparams, params) = _initparams($(modeloptions.inputparamstype), $(modeloptions.paramstype), _accumulatedparams, mappings, _namespace)))
        end
    end
    for member ∈ modeloptions.members
        if member.split
            push!(block.args, :($(member.name) = _initpuremodel($(member.name), components, Symbol($(modeloptions.name)), Symbol($(member.nameinparent)), userlocations, $(member.overwritelocations), _initialisedmodels, _namespace, _accumulatedparams, verbose)))
        else
            push!(block.args, :($(member.name) = _initmodel($(member.name), components, Symbol($(modeloptions.name)), Symbol($(member.nameinparent)), userlocations, $(member.overwritelocations), _initialisedmodels, _namespace, _accumulatedparams, verbose)))
        end
        if isnothing(member.restrictparents)
            push!(block.args, Expr(:if, :(!haskey(_initialisedmodels, :_)), :(_initialisedmodels[:_] = Dict{Symbol,Any}())))
            push!(block.args, :(_initialisedmodels[:_][Symbol($(member.nameinparent))] = $(member.name)))
        else
            for parent ∈ member.restrictparents
                push!(block.args, Expr(:if, :(!haskey(_initialisedmodels, Symbol($parent))), :(_initialisedmodels[Symbol($parent)] = Dict{Symbol,Any}())))
                push!(block.args, :(_initialisedmodels[Symbol($parent)][Symbol($(member.nameinparent))] = $(member.name)))
            end
        end
    end
    push!(block.args, :(references = $(modeloptions.references)))

    # Create object
    call = Expr(:call)
    push!(call.args, modeloptions.name)
    push!(call.args, :(components))
    if modeloptions.has_groups
        push!(call.args, :(groups))
    end
    if modeloptions.has_sites
        push!(call.args, :(sites))
    end
    if modeloptions.has_params
        push!(call.args, :(inputparams))
        push!(call.args, :(params))
        push!(call.args, :(mappings))
    end
    for (i, member) ∈ enumerate(modeloptions.members)
        push!(call.args, :($(member.name)))
    end
    if modeloptions.has_sites
        push!(call.args, :(assoc_options))
    end
    push!(call.args, :(references))
    push!(block.args, Expr(:(=), :model, call))
    push!(block.args, :(updateparams!(model)))
    push!(block.args, Expr(:return, :model))
    return Expr(:function, func_head, block)
end

#####

function _initmodel(
        model::Union{Type,Function},
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        _overwritelocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false
    )
    if model === Nothing
        return nothing
    end
    if caller ∈ keys(_initialisedmodels)
        if nameinparent ∈ _initialisedmodels[caller]
            return _initialisedmodels[caller][model]
        end
    end
    # Default key if restrictparents is `nothing`.
    if :_ ∈ keys(_initialisedmodels)
        if nameinparent ∈ _initialisedmodels[caller]
            return _initialisedmodels[caller][model]
        end
    end
    verbose && @info("Creating member model: $model")
    return model(
        components;
        userlocations,
        verbose,
        _overwritelocations,
        _initialisedmodels,
        _namespace,
        _accumulatedparams
    )
end

function _initmodel(
        model,
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        _overwritelocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false
    )
    return model
end

function _initpuremodel(
        model::Union{Type,Function},
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        _overwritelocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false
    )
    if model === Nothing
        return nothing
    end
    if caller ∈ keys(_initialisedmodels)
        if nameinparent ∈ keys(_initialisedmodels[caller])
            return _initialisedmodels[caller][model]
        end
    end
    # Default key if restrictparents is `nothing`.
    if :_ ∈ keys(_initialisedmodels)
        if nameinparent ∈ keys(_initialisedmodels[caller])
            return _initialisedmodels[caller][model]
        end
    end
    verbose && @info("Creating member pure models: $puremodels")
    return EoSVectorParam(model(
        components;
        userlocations,
        verbose,
        _overwritelocations,
        _initialisedmodels,
        _namespace,
        _accumulatedparams
    ))
end

function _initpuremodel(
        model,
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        _overwritelocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false
    )
    return model
end


"""
    _initparams(::Type{I}, ::Type{P}, rawparams, mappings, namespace)

Returns the constructed InputParamType and ParamType based on the rawparams and mappings.

Strategy:
- If a mapping is present, and it is not an identity transformation, construct a new param struct with the output param name.
- If a mapping is present, and it is an identity transformation, then it is simply a name change. Create a new param struct, but the `value` is a reference to the inputparam array.
- If no mapping is present, then the parameter is taken as-is from the raw/input params. Just point to the same param object.

## Arguments
- `Type{InputParamType}`: The type for the inputparm struct.
- `Type{ParamType`: The type for the param struct.
- `rawparams::Dict`: The dict returned from `getparams`.
- `mappings::Vector{ModelMapping}`: The list of mappings.
- `namespace::String = ""`: If there exists a header that is prefixed with `namespace * "__"`, this header will take priority over the non-prefixed version.

## Returns
Tuple{InputParamType,ParamType}: The constructed params.
"""
function _initparams(
        ::Type{I}, 
        ::Type{P},
        rawparams::Dict,
        mappings::Vector{ModelMapping},
        namespace::String = ""
    )::Tuple{I,P} where {I,P}
    inputparams_names = collect(fieldnames(I))
    inputparams_types = collect(fieldtypes(I))
    if isempty(namespace)
        inputparams_params = [rawparams[String(param)] for param ∈ inputparams_names]
    else
        inputparams_params = []
        for name ∈ input_params_names
            priority_header = namespace * "__" * String(name)
            if priority_header ∈ rawparams
                push!(inputparams_params, rawparams[priority_header])
            else
                push!(inputparams_params, rawparams[String(name)])
            end
        end
    end
    params_names = collect(fieldnames(P))
    params_types = collect(fieldtypes(P))

    # For targets with transformation
    # Map param with index of mapping
    params_mappings_index = Dict{Symbol,Integer}()
    # For targets without transformations
    # Map param with index of inputparams_params
    params_inputparams_index = Dict{Symbol,Integer}()
    for (i, mapping) ∈ enumerate(mappings)
        # For identity mappings, add to params_input_params_index
        # for making a reference to the input later
        if mapping.transformation === identity
            # Arrays will only have one element
            params_inputparams_index[first(mapping.target)] = findfirst(isequal(first(mapping.source)), inputparams_names)
            continue
        end
        for target ∈ mapping.target
            params_mappings_index[target] = i
        end
    end
    # For params not defined in the mapping, search in inputparams
    for param ∈ params_names
        haskey(params_mappings_index, param) && continue
        haskey(params_inputparams_index, param) && continue
        params_inputparams_index[param] = findfirst(isequal(param), inputparams_names)
    end
    
    # Construct the list of params
    params_params = []
    for param ∈ params_names
        if haskey(params_inputparams_index, param)
            index_inputparam = params_inputparams_index[param]
            inputparam = inputparams_params[index_inputparam]
            if param ∈ inputparams_names
                # Push shared reference with inputparam
                push!(params_params, inputparam)
            else
                # For explicit identity mappings, assume that names change
                # But still reference the same `values` array.
                inputparam_type = inputparams_types[index_inputparam] 
                push!(params_params, Base.typename(inputparam_type).wrapper(inputparam, String(param); isdeepcopy=false))
            end
        else
            # Push a newly initialised ClapeyronParam
            index_mapping = params_mappings_index[param]
            index_param = findfirst(isequal(param), params_names)
            param_type = params_types[index_param]
            indices_inputparams = indexin(mappings[index_mapping].source, inputparams_names)
            inputparams = [inputparams_params[index] for index ∈ indices_inputparams]
            sources = _get_sources(inputparams)
            if param_type <: SingleParam || param_type <: PairParam
                # Initialise empty SingleParam or PairParam
                push!(params_params, param_type(String(param), inputparams[1].components; sources))
            elseif param_type <: AssocParam
                # If AssocParam, just copy from one of the inputparams
                # We assume inputparam is also an AssocParam
                push!(params_params, Base.typename(inputparam_type).wrapper(first(inputparams), String(param); isdeepcopy=true, sources))
            end
        end
    end
    return (I(inputparams_params...), P(params_params...))
end

"""
    updateparams!(model; updatemembers)
    updateparams!(self, inputparams, params, mappings)

Given an `EoSModel`, update the values of `params` with current values of `inputparams` given the `mappings`. Note that for identity transformations, the values share the same references, so they need not be updated explictly.  

## Arguments
- `model::EoSModel`: The model to update params on.
- `updatemembers::Bool = true`: If `true` it will update all direct members of `model`. Note that it does not recursively update members of members. If it is desired for parameters of model members of model members to be updated, say for parameter estimation, ensure that it is also a direct member of the `model`.
"""
function updateparams!(model::EoSModel; updatemembers::Bool = true)::Nothing
    if hasproperty(model, :inputparams)
        updateparams!(model, model.inputparams, model.params, model.mappings, updatemembers)
    end
end

function updateparams!(
        self::EoSModel,
        inputparams::EoSParam,
        params::EoSParam,
        mappings::Vector{ModelMapping},
        updatemembers::Bool
       )::Nothing
    for mapping ∈ filter(m -> !(m.transformation === identity),  mappings)
        # At the moment, we are stil constructing new ClapeyronParams,
        # then broadcast assigning the params to them. There might be
        # some optimisations possible here, but this will maintain
        # maximum compatibility with existing functions for now.
        if mapping.self_in_args
            outputs = mapping.transformation(self, [getfield(inputparams, f) for f ∈ mapping.source]...)
        else
            outputs = mapping.transformation([getfield(inputparams, f) for f ∈ mapping.source]...)
        end
        if typeof(outputs) <: ClapeyronParam
            getfield(params, first(mapping.target)).values .= outputs.values
        else
            toupdates = [getfield(params, f) for f ∈ mapping.target]
            for (output, toupdate) ∈ zip(outputs, toupdates)
                toupdate.values .= output.values
            end
        end
    end
    if updatemembers
        for member in modelmembers(self)
            updateparams!(getfield(self, member); updatemembers=false)
        end
    end
end

"""
    createmodel(modeloptions; verbose)

Create the models in the global (or module) scope using `eval`.

See the tutorial or browse the implementations to see how this is used.
"""
function createmodel(modeloptions::ModelOptions; verbose::Bool = false)
    verbose && @info("Generating model code for " * String(modeloptions.name))
    if !isdefined(@__MODULE__, modeloptions.inputparamstype)
        inputparams = _generatecode_param_struct(modeloptions.inputparamstype, modeloptions.inputparams)
        verbose && @info(inputparams)
        eval(inputparams)
    else
        newfields = [param.name for param in modeloptions.inputparams]
        if !isempty(newfields)  # If not specified, assume users are OK with existing definition.
            oldfields = fieldnames(eval(modeloptions.inputparamstype))
            if !issetequal(newfields, oldfields)
                error("$(modeloptions.inputparamstype) is already defined with fields $oldfileds, so cannot redefined it with fields $newfields.")
            end
        end
    end

    if !isdefined(@__MODULE__, modeloptions.paramstype)
        params = _generatecode_param_struct(modeloptions.paramstype, modeloptions.params)
        verbose && @info(params)
        eval(params)
    else
        newfields = [param.name for param in modeloptions.params]
        if !isempty(newfields)  # If not specified, assume users are OK with existing definition.
            oldfields = fieldnames(eval(modeloptions.paramstype))
            if !issetequal(newfields, oldfields)
                error("$(modeloptions.paramstype) is already defined with fields $oldfileds, so cannot redefined it with fields $newfields.")
            end
        end
    end
    model = _generatecode_model_struct(modeloptions)
    verbose && @info(model)
    eval(model)
    constructor = _generatecode_model_constructor(modeloptions)
    verbose && @info(constructor)
    eval(constructor)

    eval(quote
        has_sites(::Type{<:$(modeloptions.name)}) = $(modeloptions.has_sites)
        function Base.show(io::IO, model::$(modeloptions.name))
            return eosshow(io, model)
        end
        Base.length(model::$(modeloptions.name)) = Base.length(model.components)
    end)
    if modeloptions.has_groups
        eval(quote
            has_groups(::Type{<:$(modeloptions.name)}) = true
            function Base.show(io::IO, mime::MIME"text/plain", model::$(modeloptions.name))
                return gc_eosshow(io, mime, model)
            end
            molecular_weight(model::$(modeloptions.name), z = SA[1.0]) = group_molecular_weight(model.groups, mw(model), z)
        end)
    else
        eval(quote
            has_groups(::Type{<:$(modeloptions.name)}) = false
            function Base.show(io::IO, mime::MIME"text/plain", model::$(modeloptions.name))
                return eosshow(io, mime, model)
            end
            molecular_weight(model::$(modeloptions.name), z = SA[1.0]) = comp_molecular_weight(mw(model), z)
        end)
    end
    eval(:(modelmembers(model::$(modeloptions.name)) = $([member.name for member in modeloptions.members])))
end

export ModelMapping, ModelMember, ParamField, ModelOptions, createmodel, updateparams!
