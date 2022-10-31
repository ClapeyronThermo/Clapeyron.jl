const ∅_expr = quote end |> Base.remove_linenums!

"""
    ModelMapping(source, target, transformation)

Part of `ModelOptions`. Used to define how source from csv is mapped to the target struct. The `source` and `target` should be defined in `inputparams` and `params` respecfully in `ModelOptions`.

## Fields
- `source::Union{String,Vector{String}}`: The source headers from csv. There are special symbols that extract non-parameter inputs, but note that these arguments have to be propagated by the respective functions if they are composed.
    - `:_model`: The model object itself. Note that it is inadvisable to access other params directly via getproperties on the model object as there is no guarantees that the values are up to date at the time of update. The intension is primarily for dispatch. If it is dependent on parameters, explicitly have them as arguments to the transformation function, or ensure that the parameters accessed are never changed.
    - `:_groups`: Takes `model.groups`.
    - `:_sites`: Takes `model.sites`.
    - `:_nothing`: Place `nothing` into argument.
- `target::Union{Symbol,Vector{Symbol}}`: The target parameters.
- `transformation::function = identity`: The transformation from source to target.
- `source_cache::Vector{Any}: Storing references to objects so that `getfield` does not have to be called each time mapping is updated.
- `target_cache::Vector{Any}: Storing references to objects so that `getfield` does not have to be called each time mapping is updated.
"""
struct ModelMapping
    source::Vector{Symbol}
    target::Vector{Symbol}
    transformation::Function
    source_cache::Vector{Any}
    target_cache::Vector{Any}
end
function ModelMapping(
        source::Union{Vector{Symbol},Symbol},
        target::Union{Vector{Symbol},Symbol},
        transformation::Function = identity,
    )
    if source isa Symbol
        source = [source]
    end
    if target isa Symbol
        target = [target]
    end
    return ModelMapping(
        source,
        target,
        transformation,
        Any[],
        Any[],
    )
end

"""
    ModelMember(
            name,
            default_type;
            typeconstraint,
            swappable,
            split,
            groupcontribution_allowed,
            separate_namespace,
            overwritelocations,
            overwritegrouplocations,
            restrictparents,
            nameinparent,
        )

Part of `ModelOptions`. This is used to specify the options for member models. Note that only swappable members need to be specified here. If member models contain member models with the same name as one specified here, the same instance will be used. If this is undesired behaviour, explicitly assign the relevant parents in the `restrictparents` field.

## Fields
- `name::Symbol`: The name of this member. It will be used as the fieldname for the created model object.
- `default_type::Symbol`: The default type for this member. Has to have a constructor with the same function signature. Symbol is used instead of the type to avoid circular reference.
- `typeconstraint::Union{Symbol,Nothing} = nothing`: If not `nothing`, make sure that the type to input is `::Union{typeconstraint,Type{typeconstraint}}`.
- `swappable::Bool = false`: If `true`, put argument `name` into the constructor signature to allow users to swap out member models from the default type. This is useful for, say, idealmodels, but should be set to `false` if the model very specifically depends on a particular combination of member models.
- `split::Bool = false`: If `true`, create a vector of pure models of this type.
- `nothing_allowed::Bool = false`: If `true`, the type constraint can be `Nothing`.
- `groupcontribution_allowed::Bool = false`: If `true`, create fields for injecting group contribution arguments like `usergrouplocations` and `groupdefinitions` for this member. Note that the default does not have to be a group contribution method; just that if there is a possibility of switching to one.
- `separate_namespace::Bool = true`: This is for namespace resolution of input parameter names from user-provided csv in `userlocations`. If `true`, all headers in the input csvs specified in `userlocations` should be prefixed with `{name}__`. Note that these prefixes do not nest (the `name` must be distinct per model, and so there need only be one level).
- `overwritelocations::Union{Vector{String},Nothing} = nothing`: If the model wants to overwrite the array of locations from a location in the Clapeyron database specific for this model, they may do so here. This location is not passed down to subsequent member models. This only applies to the default model.
- `overwritegrouplocations::Union{Vector{String},Nothing} = nothing`: If the model wants to overwrite the array of grouplocations from a location in the Clapeyron database specific for this model, they may do so here. This location is not passed down to subsequent member models. This only applies to the default model.
- `restrictparents::Union{Symbol,Nothing} = nothing`: This is for namespoce resolution of member models. If `nothing`, any model with a member with this name will be assigned the same instantiated object. When ambiguity arises, say if two distinct `activity` models are used, provide a vector of types to specify where this member is applicable, making use of the `nameinparent` field as well for correct asignment.
- `nameinparrent::Symbol = nothing`: If `nothing`, just take the given name. This is for when the name in member model is different from the name in current model, which allows multiple member models of the same type to be present in the main model.

## Separate namespace
To account for the possibility that parameters may share the same names across models, we can enable this option to resolve conflicts. For example, if we had a model

    A{B}

and both `A` and `B` contain distinct parameters `m`. Then we would enable `separate_namespace` for member `B`, and in the parameters file, if there is a column called `B__m`, it will take precedence over the column `m`. If column `B__m` is not present, it will find that column `m` is already initialised when getting parameters for `A`, so it will just take the same reference to that parameter.

Note that it can get very difficult to keep track of which models have overwritten namespaces, especially when they are used as member models themselves, so use this sparingly, and if possible, try to avoid column name conflicts in the input parameters whenever possible.

## Restrict parents
The constructor caches all initialised models, and if a member model later initialises another member model with the same `nameinparent` (set to `name` unless specified otherwise), it will just take the same reference to the previously initialised model. Note that ony models initialised the main non-member model will add to the cache.

The order of initialisation is as provided in the `members` field vector. For example, if we have model `A` has `Y` and `X` as member models (note the order) and model `X` is a member of model `Y`. If we specified it as

    A{Y{X},X}

then without any other instructions, it will first initialise `Y` with the default `X`, then intialise a second `X` for the entry in `A`.

If we instead had it written as

    A{X,Y{X}}

unless otherwise specified, it first initialises `X`, then it initialises `Y`, where it finds that `X` already initialised, so it just takes a reference to it.

If we specify `[:A]` in `restrictparents` for `ModelMember` of `X`, then a separate instance of a `X` (the default) will be constructed for `Y`.

The base model `A` will always contain an instance of any member specified in `members`. So in this case, specifying `[:Y]` for the `restrictparents` field here will have the same effect as `nothing` (both referencing the same `X`). In effect, this means that for parameter estimation, we will never need to look into modifying any params deeper than one level of members.

If we wish the `X` in `A` be distinct from the `X` in `Y`, yet have them both modular (swappable from the main constructor interface), we can give them separate names in the base model, but change the `nameinparent` to reflect the correct name. In this case, we will have

    A{X, X2, Y{X}}

where for `X`, we have `restrictparents` as `[:A]`, whereas for `X2`, we have `restrictparents` as `[:Y]` (superfluous here, as per discussion above, but provides clarification on usage), and `nameinparent` as `X`. Then `X2` here is the same instance as the `X` in `Y`. Note that the order of the `restrictparents` can be swapped, so thet he `X` in `Y` is the same instance as the `X` in `A`, instead of `X2` as before. However, in this case, both `restrictparams` have to be specified.
"""
struct ModelMember
    name::Symbol
    default_type::Symbol
    typeconstraint::Union{Symbol,Nothing}
    swappable::Bool
    split::Bool
    nothing_allowed::Bool
    groupcontribution_allowed::Bool
    separate_namespace::Bool
    overwritelocations::Union{Vector{String},Nothing}
    overwritegrouplocations::Union{Vector{String},Nothing}
    restrictparents::Union{Vector{Symbol},Nothing}
    nameinparent::Symbol
    namespace::String
end
function ModelMember(
        name::Symbol,
        default_type::Symbol;
        typeconstraint::Union{Symbol,Nothing} = nothing,
        swappable::Bool = true,
        split::Bool = false,
        nothing_allowed::Bool = false,
        groupcontribution_allowed::Bool = false,
        separate_namespace::Bool = false,
        overwritelocations::Union{Vector{String},Nothing} = nothing,
        overwritegrouplocations::Union{Vector{String},Nothing} = nothing,
        restrictparents::Union{Vector{Symbol},Nothing} = nothing,
        nameinparent::Union{Symbol,Nothing} = nothing,
    )
    if isnothing(nameinparent)
        nameinparent = name
    end
    return ModelMember(
        name,
        default_type,
        typeconstraint,
        swappable,
        split,
        nothing_allowed,
        groupcontribution_allowed,
        separate_namespace,
        overwritelocations,
        overwritegrouplocations,
        restrictparents,
        nameinparent,
        separate_namespace ? string(name) : "",
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
- `locations::Vector{String} = String[]`: Default locations in Clapeyron database to look for parameters. Note that for user-created models, it is easier to use the `userlocations` parameter.
- `grouplocations::Vector{String} = String[]`: Default grouplocations in Clapeyron database to look for parameters. Note that for user-created models, it is easier to use the `usergrouplocations` parameter.
- `inputparams::Vector{ParamField} = ParamField[]`: A list of relevant source parameters. The model constructor will extract the String headers according to the Symbol name given. If the `inputparams` is not specified and there is no parent, make it equal to `params`.
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
    grouplocations::Vector{String}
    inputparams::Vector{ParamField}
    params::Vector{ParamField}
    mappings::Vector{ModelMapping}
    has_params::Bool
    has_components::Bool
    has_sites::Bool
    has_groups::Bool
    param_options::Union{ParamOptions,Nothing}
    assoc_options::Union{AssocOptions,Nothing}
    references::Union{Vector{String},Nothing}
    inputparamstype::Symbol
    paramstype::Symbol
end

function ModelOptions(
        name::Symbol;
        supertype::Union{DataType,Nothing} = nothing,
        parent::Union{ModelOptions,Nothing} = nothing,
        members::Union{Vector{ModelMember},Nothing} = nothing,
        locations::Union{Vector{String},Nothing} = nothing,
        grouplocations::Union{Vector{String},Nothing} = nothing,
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
        paramstype::Union{Symbol,Nothing} = nothing,
    )
    if !isnothing(parent)
        isnothing(supertype) && (supertype = parent.supertype)
        isnothing(members) && (members = parent.members)
        isnothing(locations) && (locations = parent.locations)
        isnothing(grouplocations) && (grouplocations = parent.grouplocations)
        if isnothing(inputparams)
            isnothing(inputparamstype) && (inputparamstype = parent.inputparamstype)
            inputparams = parent.inputparams
        end
        if isnothing(params)
            isnothing(paramstype) && (paramstype = parent.paramstype)
            params = parent.params
        end
        isnothing(mappings) && (mappings = parent.mappings)
        isnothing(has_params) && (has_params = parent.has_params)
        isnothing(has_components) && (has_components = parent.has_components)
        isnothing(has_sites) && (has_sites = parent.has_sites)
        isnothing(has_groups) && (has_groups = parent.has_groups)
        isnothing(param_options) && (param_options = parent.param_options)
        isnothing(assoc_options) && (assoc_options = parent.assoc_options)
        isnothing(references) && (references = parent.references)
    else
        # Default values if no parent is given.
        isnothing(supertype) && (supertype = EoSModel)
        isnothing(members) && (members = ModelMember[])
        isnothing(locations) && (locations = String[])
        isnothing(grouplocations) && (grouplocations = String[])
        isnothing(params) && (params = ParamField[])
        isnothing(inputparams) && (inputparams = params)
        isnothing(mappings) && (mappings = ModelMapping[])
        isnothing(has_params) && (has_params = true)
        isnothing(has_components) && (has_components = true)
        isnothing(has_sites) && (has_sites = false)
        isnothing(has_groups) && (has_groups = false)
    end
    isnothing(paramstype) && (paramstype = Symbol(name, :Param))
    if isnothing(inputparamstype)
        if inputparams === params
            inputparamstype = paramstype
        else
            inputparamstype = Symbol(name, :InputParam)
        end
    end
    return ModelOptions(
        name,
        supertype,
        parent,
        members,
        locations,
        grouplocations,
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
        inputparamstype,
        paramstype,
    )
end

##### Expr generators #####

"""
    _generatecode_param_struct(modeloptions)

Generates the parameter struct. This can be used for both `inputparams` and `params`.
"""
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

"""
    _generatecode_model_struct(modeloptions)

Generates the model struct.
"""
function _generatecode_model_struct(modeloptions::ModelOptions)::Expr
    if isempty(modeloptions.members)
        defheader = modeloptions.name
        typestatement = modeloptions.supertype
    else
        defheader = Expr(:curly, modeloptions.name)
        for i ∈ 1:length(modeloptions.members)
            push!(defheader.args, Symbol(:M, Symbol(i)))
        end
        typestatement = Expr(:where, modeloptions.supertype)
        for (i, member) ∈ enumerate(modeloptions.members)
            if isnothing(member.typeconstraint)
                push!(typestatement.args, Symbol(:M, Symbol(i)))
            else
                if member.nothing_allowed
                    push!(typestatement.args, Expr(:<:, Symbol(:M, Symbol(i)), :(Union{$(member.typeconstraint),Nothing})))
                else
                    push!(typestatement.args, Expr(:<:, Symbol(:M, Symbol(i)), member.typeconstraint))
                end
            end
        end
    end

    block = Expr(:block)
    modeloptions.has_components && push!(block.args, :(components::Vector{String}))
    modeloptions.has_groups && push!(block.args, :(groups::GroupParam))
    modeloptions.has_sites && push!(block.args, :(sites::SiteParam))
    if modeloptions.has_params
        push!(block.args, :(inputparams::$(modeloptions.inputparamstype)))
        push!(block.args, :(params::$(modeloptions.paramstype)))
        push!(block.args, :(mappings::Vector{ModelMapping}))
    end
    for (i, member) ∈ enumerate(modeloptions.members)
        if member.split
            push!(block.args, Expr(:(::), member.name, Expr(:curly, :EoSVectorParam, Symbol(:M, Symbol(i)))))
        else
            push!(block.args, Expr(:(::), member.name, Symbol(:M, Symbol(i))))
        end
    end
    if modeloptions.has_sites
        push!(block.args, :(assoc_options::AssocOptions))
    end
    if modeloptions.references !== nothing
        push!(block.args, :(references::Vector{String}))
    end
    return Expr(:struct, false, Expr(:<:, defheader, typestatement), block)
end

function is_empty_model(model::ModelOptions)
    return !(model.has_sites | model.has_groups | model.has_components | model.has_params | !isnothing(model.references))
end

"""
    _generate_model_constructor(modeloptions)

Generates a model constructor with the folowing function signature:

## Arguments
- `components::Union{String, Vector{String}}`: The components for this model.
_ `userlocations::Vector{String} = String[]`: An array of filepaths to csvs specified by the user for parameters.
_ `usergrouplocations::Vector{String} = String[]`: An array of filepaths to csvs specified by the user for groups. If this is not a GC model, it is not used.
_ `groupdefinitions::Vector{GroupDefinition} = GroupDefinition[]`: An array of group definitions. If this is not a GC model, it is not used.
- `verbose::Bool = false`: Print more information.
- `param_options::ParamOptions`: This is present for models with parameters to change the behaviour of `getparams`.
- `assoc_options::AssocOptions`: This is present for models with sites to change the behaviour of association.
- various member models: For models with modular components, each member model can be swapped out here. For example, the `idealmodal` or `activity`.
- various {member}_usergrouplocations: If member model is a GC model, this will be present to pass location to csvs for group definitions that can be specified by the user.
- various {member}_groupdefinitions: If member model is a GC model, this will be present to pass group definiitons.

## System arguments
These are used by Clayeyron when initialising member models. They should not be modified unless you know what you are doing.
- `_overwritelocations::Union{Vector{String}} = nothing`: If not `nothing`, overwrite the `locations` for this model.
- `_overwritegrouplocations::Union{Vector{String}} = nothing`: If not `nothing`, overwrite the `grouplocations` for this model.
- `_initialisedmodels::Dict{Symbol, Dict{Symbol, Any}} = Dict{Symbol, Dict{Symbol, Any}}()`: If model is already initialised, just take a reference to it unless specified otherwise in the `MemberModel`.
- `_namespace::String = ""`: Give higher priority to columns prepended with `{membermodel_name}__{inputparam_name}` if present.
- `_accumulatedparams::Dict{String, ClapeyronParam} = Dict{String, ClapeyronParam}())`: If there are parameters with the same name, just point to existing reference.
- `_ismembermodel::Bool = false`: Automatically changed to `true` if called from `_initmodel` or `initpuremodel`. One reason for this is to ensure that only non-member models push to `_initialisedmodels`.
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
    push!(parameters.args, Expr(:kw, :(usergrouplocations::Vector{String}), :(String[])))
    push!(parameters.args, Expr(:kw, :(groupdefinitions::Vector{GroupDefinition}), :(GroupDefinition[])))
    push!(parameters.args, Expr(:kw, :(verbose::Bool), :false))
    if modeloptions.has_params
        push!(parameters.args, Expr(:kw, :(param_options::ParamOptions), :($(modeloptions.param_options))))
    end
    if modeloptions.has_sites
        push!(parameters.args, Expr(:kw, :(assoc_options::AssocOptions), :($(modeloptions.assoc_options))))
    end
    for member ∈ modeloptions.members
        if member.swappable
            if isnothing(member.typeconstraint)
                push!(parameters.args, Expr(:kw, :($(member.name)), :($(member.default_type))))
            else
                if member.nothing_allowed
                    push!(parameters.args, Expr(:kw, :($(member.name)::Union{$(member.typeconstraint),Type{<:Union{$(member.typeconstraint),Nothing}},Nothing}), :($(member.default_type))))
                else
                    push!(parameters.args, Expr(:kw, :($(member.name)::Union{$(member.typeconstraint),Type{<:$(member.typeconstraint)}}), :($(member.default_type))))
                end
            end
        end
        if member.groupcontribution_allowed
            push!(parameters.args, Expr(:kw, Symbol(:($(member.name)), :_usergrouplocations), :(String[])))
            push!(parameters.args, Expr(:kw, Symbol(:($(member.name)), :_groupdefinitions), :(GroupDefinition[])))
        end
    end
    push!(parameters.args, Expr(:kw, :(_overwritelocations::Union{Vector{String},Nothing}), :nothing))
    push!(parameters.args, Expr(:kw, :(_overwritegrouplocations::Union{Vector{String},Nothing}), :nothing))
    push!(parameters.args, Expr(:kw, :(_initialisedmodels::Dict{Symbol,Dict{Symbol,Any}}), :(Dict{Symbol,Dict{Symbol,Any}}(:_ => Dict{Symbol,Any}()))))
    push!(parameters.args, Expr(:kw, :(_namespace::String), :""))
    push!(parameters.args, Expr(:kw, :(_accumulatedparams::Dict{String,ClapeyronParam}), :(Dict{String,ClapeyronParam}())))
    push!(parameters.args, Expr(:kw, :(_ismembermodel::Bool), :false))
    push!(func_head.args, parameters)
    # Now positional args
    
    # `components` is mandatory
    push!(func_head.args, :(components::Union{String,Vector{String}}))
    
    # Creating function body
    block = Expr(:block)
    if modeloptions.has_groups
        push!(block.args, :(grouplocations = $(modeloptions.grouplocations)))
        push!(block.args, Expr(:if, :(!isnothing(_overwritegrouplocations)), :(locations = _overwritegrouplocations)))
        push!(block.args, :(groups = GroupParam(components, grouplocations; groupdefinitions, usergrouplocations, verbose, param_options)))
    end
    if modeloptions.has_params
        push!(block.args, :(mappings = $(modeloptions.mappings)))
        push!(block.args, :(locations = $(modeloptions.locations)))
        push!(block.args, Expr(:if, :(!isnothing(_overwritelocations)), :(locations = _overwritelocations)))
        if modeloptions.has_sites
            if modeloptions.has_groups
                push!(block.args, :((rawparams, sites) = getparams(groups, locations, param_options; userlocations, verbose)))
            else
                push!(block.args, :((rawparams, sites) = getparams(components, locations, param_options; userlocations, verbose)))
            end
        else
            if modeloptions.has_groups
                push!(block.args, :(rawparams = getparams(groups, locations, param_options; userlocations, verbose)))
            else
                push!(block.args, :(rawparams = getparams(components, locations, param_options; userlocations, verbose)))
            end
        end
        push!(block.args, :(_accumulatedparams = merge!(rawparams, _accumulatedparams)))
        if modeloptions.has_groups
            push!(block.args, :((inputparams, params) = _initparams(groups.flattenedgroups, $(modeloptions.inputparamstype), $(modeloptions.paramstype), _accumulatedparams, mappings, _namespace)))
        else
            push!(block.args, :((inputparams, params) = _initparams(components, $(modeloptions.inputparamstype), $(modeloptions.paramstype), _accumulatedparams, mappings, _namespace)))
        end
    end

    # Member model initialisation.
    for member ∈ modeloptions.members
        if !member.swappable
            push!(block.args, :($(member.name) = $(member.default_type)))
        end
        if !isnothing(member.overwritelocations)
            push!(block.args, Expr(:if, :($(member.default_type) === $(member.name)), :(overwritelocations = $(member.overwritelocations)), :(overwritelocations = nothing)))
        else
            push!(block.args, :(overwritelocations = nothing))
        end
        if !isnothing(member.overwritegrouplocations)
            push!(block.args, Expr(:if, :($(member.default_type) === $(member.name)), :(overwritegrouplocations = $(member.overwritegrouplocations)), :(overwritegrouplocations = nothing)))
        else
            push!(block.args, :(overwritegrouplocations = nothing))
        end
        if member.groupcontribution_allowed
            if member.split
                push!(block.args, :($(member.name) = _initpuremodel(
                    $(member.name),
                    components,
                    $(Meta.quot(:($(modeloptions.name)))),
                    $(Meta.quot(:($(member.nameinparent)))),
                    userlocations,
                    $(Symbol(:($(member.name)), :_usergrouplocations)),
                    $(Symbol(:($(member.name)), :_groupdefinitions)),
                    overwritelocations,
                    overwritegrouplocations,
                    _initialisedmodels,
                    $(member.namespace),
                    _accumulatedparams,
                    verbose,
                )))
            else
                push!(block.args, :($(member.name) = _initmodel(
                    $(member.name),
                    components,
                    $(Meta.quot(:($(modeloptions.name)))),
                    $(Meta.quot(:($(member.nameinparent)))),
                    userlocations,
                    $(Symbol(:($(member.name)), :_usergrouplocations)),
                    $(Symbol(:($(member.name)), :_groupdefinitions)),
                    overwritelocations,
                    overwritegrouplocations,
                    _initialisedmodels,
                    $(member.namespace),
                    _accumulatedparams,
                    verbose,
                )))
            end
        else
            if member.split
                push!(block.args, :($(member.name) = _initpuremodel(
                    $(member.name),
                    components,
                    $(Meta.quot(:($(modeloptions.name)))),
                    $(Meta.quot(:($(member.nameinparent)))),
                    userlocations,
                    String[],
                    GroupDefinition[],
                    overwritelocations,
                    overwritegrouplocations,
                    _initialisedmodels,
                    $(member.namespace),
                    _accumulatedparams,
                    verbose,
                )))
            else
                push!(block.args, :($(member.name) = _initmodel(
                    $(member.name),
                    components,
                    $(Meta.quot(:($(modeloptions.name)))),
                    $(Meta.quot(:($(member.nameinparent)))),
                    userlocations,
                    String[],
                    GroupDefinition[],
                    overwritelocations,
                    overwritegrouplocations,
                    _initialisedmodels,
                    $(member.namespace),
                    _accumulatedparams,
                    verbose,
                )))
            end
        end
        # Adding initialised model to `_initilasiedmodels` cache.
        if isnothing(member.restrictparents)
            push!(block.args, Expr(
                :if,
                :(!_ismembermodel),
                :(_initialisedmodels[:_][$(Meta.quot(:($(member.nameinparent))))] = $(member.name))
            ))
        else
            for parent ∈ member.restrictparents
                push!(block.args, Expr(
                    :if,
                    :(!_ismembermodel),
                    Expr(
                        :if,
                        :(!haskey(_initialisedmodels, $(Meta.quot(:($parent))))),
                        :(_initialisedmodels[Symbol($parent)] = Dict{Symbol,Any}())
                    )))
                push!(block.args, Expr(
                    :if,
                    :(!_ismembermodel),
                    :(_initialisedmodels[$(Meta.quot(:($parent)))][$(Meta.quot(:($(member.nameinparent))))] = $(member.name))))
            end
        end
    end
    push!(block.args, :(references = $(modeloptions.references)))

    # Create object
    call = Expr(:call)
    push!(call.args, modeloptions.name)
    if modeloptions.has_components
        push!(call.args, :(components))
    end
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
    if modeloptions.references !== nothing
        push!(call.args, :(references))
    end
    push!(block.args, Expr(:(=), :model, call))
    push!(block.args, :(updateparams!(model; updatemembers=false)))
    push!(block.args, Expr(:return, :model))
    return Expr(:function, func_head, block)
end

#####

"""
    _initmodel(model, components, caller, nameinparent, userlocations, usergrouplocations,
            groupdefinitions, _overwritelocations, _overwritegrouplocations, _initialisedmodels,
            _namespace, _accumulatedparams, verbose)
    _initpuremodel(model, components, caller, nameinparent, userlocations, usergrouplocations,
            groupdefinitions, _overwritelocations, _overwritegrouplocations, _initialisedmodels,
            _namespace, _accumulatedparams, verbose)

These functions will call the construction of a member model, but only if it is not already initialised. Look at documentation for `MemberModel` for more information about how models are cached in the order that they are created.
"""
function _initmodel(
        model::Union{Type,Function},
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        usergrouplocations::Vector{String},
        groupdefinitions::Vector{GroupDefinition},
        _overwritelocations::Union{Vector{String},Nothing},
        _overwritegrouplocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false,
    )
    if model === Nothing
        return nothing
    end
    if caller ∈ keys(_initialisedmodels)
        if nameinparent ∈ keys(_initialisedmodels[caller])
            return _initialisedmodels[caller][nameinparent]
        end
    end
    # Default key if restrictparents is `nothing`.
    if :_ ∈ keys(_initialisedmodels)
        if nameinparent ∈ keys(_initialisedmodels[:_])
            return _initialisedmodels[:_][nameinparent]
        end
    end
    verbose && @info("Creating member model: $model")
    return model(
        components;
        userlocations,
        usergrouplocations,
        groupdefinitions,
        verbose,
        _overwritelocations,
        _overwritegrouplocations,
        _initialisedmodels,
        _namespace,
        _accumulatedparams,
        _ismembermodel=true,
    )
end

function _initmodel(
        model,
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        usergrouplocations::Vector{String},
        groupdefinitions::Vector{GroupDefinition},
        _overwritelocations::Union{Vector{String},Nothing},
        _overwritegrouplocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false,
    )
    return model
end

function _initpuremodel(
        model::Union{Type,Function},
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        usergrouplocations::Vector{String},
        groupdefinitions::Vector{GroupDefinition},
        _overwritelocations::Union{Vector{String},Nothing},
        _overwritegrouplocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false,
    )
    if model === Nothing
        return nothing
    end
    if caller ∈ keys(_initialisedmodels)
        if nameinparent ∈ keys(_initialisedmodels[caller])
            return _initialisedmodels[caller][nameinparent]
        end
    end
    # Default key if restrictparents is `nothing`.
    if :_ ∈ keys(_initialisedmodels)
        if nameinparent ∈ keys(_initialisedmodels[:_])
            return _initialisedmodels[:_][nameinparent]
        end
    end
    verbose && @info("Creating member pure models: $model")
    return EoSVectorParam(model(
        components;
        userlocations,
        verbose,
        _overwritelocations,
        _initialisedmodels,
        _namespace,
        _accumulatedparams,
        _ismembermodel=true,
    ))
end

function _initpuremodel(
        model,
        components::Vector{String},
        caller::Symbol,
        nameinparent::Symbol,
        userlocations::Vector{String},
        usergrouplocations::Vector{String},
        groupdefinitions::Vector{GroupDefinition},
        _overwritelocations::Union{Vector{String},Nothing},
        _overwritegrouplocations::Union{Vector{String},Nothing},
        _initialisedmodels::Dict{Symbol,Dict{Symbol,Any}},
        _namespace::String,
        _accumulatedparams::Dict{String,ClapeyronParam},
        verbose::Bool = false,
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
- `components::Vector{String}`: The list of components (or groups).
- `Type{InputParamType}`: The type for the inputparm struct.
- `Type{ParamType`: The type for the param struct.
- `rawparams::Dict`: The dict returned from `getparams`.
- `mappings::Vector{ModelMapping}`: The list of mappings.
- `namespace::String = ""`: If there exists a header that is prefixed with `namespace * "__"`, this header will take priority over the non-prefixed version.

## Returns
Tuple{InputParamType,ParamType}: The constructed params.
"""
function _initparams(
        components::Vector{String},
        ::Type{I}, 
        ::Type{P},
        rawparams::Dict,
        mappings::Vector{ModelMapping},
        namespace::String = "",
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
            inputparams = [inputparams_params[index] for index ∈ filter(!isnothing, indices_inputparams)]
            sources = _get_sources(inputparams)
            if param_type <: SingleParam || param_type <: PairParam
                # Initialise empty SingleParam or PairParam
                push!(params_params, param_type(String(param), components; sources))
            elseif param_type <: AssocParam
                # If AssocParam, just copy from one of the inputparams
                # We assume inputparam is also an AssocParam
                push!(params_params, Base.typename(inputparam_type).wrapper(components, String(param); isdeepcopy=true, sources))
            end
        end
    end
    return (I(inputparams_params...), P(params_params...))
end

"""
    updateparams!(model; updatemembers)

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

function updateparams!(::Nothing; updatemembers::Bool=false)::Nothing end

function updateparams!(
        model::EoSModel,
        inputparams::EoSParam,
        params::EoSParam,
        mappings::Vector{ModelMapping},
        updatemembers::Bool,
       )::Nothing
    if updatemembers
        for member in modelmembers(model)
            updateparams!(getfield(model, member); updatemembers=false)
        end
    end
    for mapping ∈ filter(m -> !(m.transformation === identity),  mappings)
        # At the moment, we are stil constructing new ClapeyronParams,
        # then broadcast assigning the params to them. There might be
        # some optimisations possible here, but this will maintain
        # maximum compatibility with existing functions for now.
        if isempty(mapping.source_cache)
            for arg in mapping.source
                if arg === :_model
                    push!(mapping.source_cache, model)
                elseif arg === :_groups
                    push!(mapping.source_cache, model.groups)
                elseif arg === :_sites
                    push!(mapping.source_cache, model.sites)
                else
                    push!(mapping.source_cache, getfield(inputparams, arg))
                end
            end
        end
        if isempty(mapping.target_cache)
            append!(mapping.target_cache, [getfield(params, f) for f ∈ mapping.target])
        end
        outputs = mapping.transformation(mapping.source_cache...)
        if outputs isa ClapeyronParam
            first(mapping.target_cache).values .= outputs.values
        elseif outputs isa Union{Tuple, Vector}
            for (output, target) ∈ zip(outputs, mapping.target_cache)
                if output isa ClapeyronParam
                    target.values .= output.values
                else
                    target.values .= output
                end
            end
        else
            first(mapping.target_cache).values .= outputs
        end
    end
end

function _inputparams_expr(modeloptions::ModelOptions,verbose = false)::Expr
    if !isempty(modeloptions.inputparams)
        if !isdefined(@__MODULE__, modeloptions.inputparamstype)
            inputparams = _generatecode_param_struct(modeloptions.inputparamstype, modeloptions.inputparams)
            verbose && @info(inputparams)
            return inputparams
        else
            newfields = [param.name for param in modeloptions.inputparams]
            newtypes = [param.type for param in modeloptions.inputparams]
            oldfields = fieldnames(eval(modeloptions.inputparamstype))
            oldtypes = fieldtypes(eval(modeloptions.inputparamstype))
            if !issetequal(newfields, oldfields)
                error("$(modeloptions.inputparamstype) is already defined with fields $oldfields, so it cannot be redefined with fields $newfields.")
            end
            if !issetequal(newtypes, oldtypes)
                error("$(modeloptions.inputparamstype) has the same fields $oldfields as an existing definition, but it has types $oldtypes, insted of the $newtypes.")
            end
            return ∅_expr
        end
    end
    return ∅_expr
end

function _params_expr(modeloptions::ModelOptions, verbose = false)::Expr
    if !isempty(modeloptions.params)
        if !isdefined(@__MODULE__, modeloptions.paramstype)
            params =  _generatecode_param_struct(modeloptions.paramstype, modeloptions.params)
            verbose && @info(params)
            return params
        else
            newfields = [param.name for param in modeloptions.params]
            newtypes = [param.type for param in modeloptions.params]
            oldfields = fieldnames(eval(modeloptions.paramstype))
            oldtypes = fieldtypes(eval(modeloptions.paramstype))
            if !issetequal(newfields, oldfields)
                error("$(modeloptions.paramstype) is already defined with fields $oldfields, so it cannot be redefined with fields $newfields.")
            end
            if !issetequal(newtypes, oldtypes)
                error("$(modeloptions.paramstype) has the same fields $oldfields as an existing definition, but it has types $oldtypes, insted of the $newtypes.")
            end
            return ∅_expr
        end
    end
    return ∅_expr
end

function _has_groups_expr(modeloptions::ModelOptions)::Expr
    if modeloptions.has_groups
        return quote
            has_groups(::Type{<:$(modeloptions.name)}) = true
            function Base.show(io::IO, mime::MIME"text/plain", model::$(modeloptions.name))
                return gc_eosshow(io, mime, model)
            end
            molecular_weight(model::$(modeloptions.name), z = SA[1.0]) = group_molecular_weight(model.groups, mw(model), z)
        end
    else
        return quote
            function Base.show(io::IO, mime::MIME"text/plain", model::$(modeloptions.name))
                return eosshow(io, mime, model)
            end
            molecular_weight(model::$(modeloptions.name), z = SA[1.0]) = comp_molecular_weight(mw(model), z)
        end
    end
end

function _has_sites_expr(modeloptions::ModelOptions)::Expr
    if modeloptions.has_sites
        return :(has_sites(::Type{<:$(modeloptions.name)}) = $(modeloptions.has_sites))
    else
        return ∅_expr
    end
end

function _length_expr(modeloptions::ModelOptions)::Expr
    if modeloptions.has_components
        return :(Base.length(model::$(modeloptions.name)) = Base.length(model.components))
    else
        return :(Clapeyron.components(model::$(modeloptions.name)) = nothing)
    end
end

function _short_show_expr(modeloptions::ModelOptions)::Expr
    res = quote
        function Base.show(io::IO, model::$(modeloptions.name))
            return eosshow(io, model)
        end
    end
    return res
end

function _modelmember_expr(modeloptions::ModelOptions)::Expr
    tp = ((member.name for member in modeloptions.members)...,)
    return :(modelmembers(model::$(modeloptions.name)) = $tp)
end
"""
    createmodel(modeloptions; verbose)

Create the models in the global (or module) scope using `eval`.

The structs constructed to namespace are
- modeloptions.paramstype (if necessary and not already exist)
- modeloptions.inputparamstype (if necessary and not already exist)
- modeloptions.name

Also defines for this modeloptions.name
- constructor for modeloptions.name
- has_sites
- Base.length
- Base.show
- has_groups
- molecular_weight
- modelmembers

See the tutorial or browse the implementations to see how this is used.
"""
function createmodel(modeloptions::ModelOptions; verbose::Bool = false)
    verbose && @info("Generating model code for " * String(modeloptions.name))
    eval(_inputparams_expr(modeloptions,verbose))
    eval(_params_expr(modeloptions,verbose))
    model = _generatecode_model_struct(modeloptions)
    verbose && @info(model)
    eval(model)
    constructor = _generatecode_model_constructor(modeloptions)
    verbose && @info(constructor)
    eval(constructor)
    eval(_has_sites_expr(modeloptions))
    eval(_short_show_expr(modeloptions))
    eval(_length_expr(modeloptions))
    eval(_has_groups_expr(modeloptions))
    eval(_modelmember_expr(modeloptions))
end

is_splittable(::Vector{ModelMapping}) = false

export ModelMapping, ModelMember, ParamField, ModelOptions, createmodel, updateparams!

macro createmodel2(modeloptions_expr,verbose_expr = false)
    verbose = verbose_val(verbose_expr)
    modeloptions = @eval $modeloptions_expr
    modeloptions isa ModelOptions || throw(error("input model options not a ModelOptions struct."))

    verbose && @info("Generating model code for " * String(modeloptions.name))
    inputparams = _inputparams_expr(modeloptions,verbose) |> Base.remove_linenums!
    params = _params_expr(modeloptions,verbose) |> Base.remove_linenums!
    model = _generatecode_model_struct(modeloptions) |> Base.remove_linenums!
    constructor = _generatecode_model_constructor(modeloptions)  |> Base.remove_linenums!
    sites = _has_sites_expr(modeloptions) |> Base.remove_linenums!
    _length =  _length_expr(modeloptions) |> Base.remove_linenums!
    short_show = _short_show_expr(modeloptions) |> Base.remove_linenums!
    groups = _has_groups_expr(modeloptions) |> Base.remove_linenums!
    modelmember = _modelmember_expr(modeloptions) |> Base.remove_linenums!
    res =  quote
        $inputparams
        $params
        $model
        $constructor
        $sites
        $short_show
        $_length
        $groups
        $modelmember
    end  |> Base.remove_linenums!


    verbose && println(res)
    return ∅_expr
end

function verbose_val(v::Expr)
    sym = v.args[1]
    val = v.args[2]
    if sym !== :verbose
        throw(error("incorrect keyword: expected verbose, got", string(sym)))
    end
    if !isa(val,Bool)
        throw(error("incorrect value: expected verbose::Bool, got", string(typeof(val))))
    end
    return val
end

verbose_val(v::Bool) = v

#= code generated by the createmodel function (via @createmodel2 macro)
begin
    struct sPCSAFT{M1} <: (Clapeyron.sPCSAFTModel where M1 <: IdealModel)
        components::Vector{String}
        sites::SiteParam
        inputparams::PCSAFTInputParam
        params::PCSAFTParam
        mappings::Vector{ModelMapping}
        idealmodel::M1
        assoc_options::AssocOptions
        references::Vector{String}
    end
    function sPCSAFT(components::Union{String, Vector{String}}; 
        userlocations::Vector{String} = String[], 
        usergrouplocations::Vector{String} = String[], 
        groupdefinitions::Vector{GroupDefinition} = GroupDefinition[],
         verbose::Bool = false,
        param_options::ParamOptions = ParamOptions(String[], String[], ["dipprnumber", "smiles"], "species", "source", "site", "groups", true, Dict("e" => "n_e", "e2" => "n_e2", "e1" => "n_e1", "H" => "n_H"),true, "~|~"),
        assoc_options::AssocOptions = AssocOptions(1.0e-12, 1.0e-12, 1000, 0.5, :sparse_nocombining), 
        idealmodel::Union{IdealModel, Type{<:IdealModel}} = BasicIdeal, 
        idealmodel_usergrouplocations = String[], 
        idealmodel_groupdefinitions = GroupDefinition[], 
        _overwritelocations::Union{Vector{String}, Nothing} = nothing, 
        _overwritegrouplocations::Union{Vector{String}, Nothing} = nothing, 
        _initialisedmodels::Dict{Symbol, Dict{Symbol,Any}} = Dict{Symbol, Dict{Symbol, Any}}(:_ => Dict{Symbol, Any}()), _namespace::String = "", 
        _accumulatedparams::Dict{String, ClapeyronParam} = Dict{String, ClapeyronParam}(), 
        _ismembermodel::Bool = false)

        mappings = ModelMapping[ModelMapping([:m], [:segment], identity, Any[], Any[]), ModelMapping([:sigma], [:sigma], Clapeyron.sigma_LorentzBerthelot ∘ Clapeyron.var"#416#417"(), Any[], Any[]), ModelMapping([:epsilon, :k], [:epsilon], Clapeyron.epsilon_LorentzBerthelot, Any[], Any[])]
        
        locations = ["SAFT/PCSAFT/sPCSAFT/", "properties/molarmass.csv"]
        
        if !(isnothing(_overwritelocations))
            locations = _overwritelocations
        end
        (rawparams, sites) = getparams(components, locations, param_options; userlocations, verbose)
        merge!(_accumulatedparams, merge(rawparams, _accumulatedparams))
        (inputparams, params) = _initparams(components, PCSAFTInputParam, PCSAFTParam, _accumulatedparams, mappings, _namespace)
        overwritelocations = nothing
        overwritegrouplocations = nothing
        idealmodel = _initmodel(idealmodel, components, :sPCSAFT, :idealmodel, userlocations, idealmodel_usergrouplocations, idealmodel_groupdefinitions, overwritelocations, overwritegrouplocations, _initialisedmodels, "", _accumulatedparams, verbose)
        if !_ismembermodel
            (_initialisedmodels[:_])[:idealmodel] = idealmodel
        end
        references = ["10.1021/ie020753p"]
        model = sPCSAFT(components, sites, inputparams, params, mappings, idealmodel, assoc_options, references)
        updateparams!(model; updatemembers = false)
        return model
    end
    has_sites(::Type{<:sPCSAFT}) = begin
            true
        end
    begin
        function Base.show(io::IO, model::sPCSAFT)
            return eosshow(io, model)
        end
    end
    Base.length(model::sPCSAFT) = begin
            Base.length(model.components)
        end
        
    begin
        has_groups(::Type{<:sPCSAFT}) = begin
                false
            end
        function Base.show(io::IO, mime::MIME"text/plain", model::sPCSAFT)
            return eosshow(io, mime, model)
        end
        molecular_weight(model::sPCSAFT, z = SA[1.0]) = begin
            comp_molecular_weight(mw(model), z)
        end
    end
    modelmembers(model::sPCSAFT) = begin
            (:idealmodel,)
        end
end

=#