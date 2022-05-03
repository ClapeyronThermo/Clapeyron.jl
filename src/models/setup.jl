"""
    ModelMapping(source, target, transformation)

Part of `ModelOptions`. Used to define how source from csv is mapped to the target struct. The `source` and `target` should be defined in `sourceparams` and `params` respecfully in `ModelOptions`.

## Fields
- `source::Union{String,Vector{String}}`: The source headers from csv.
- `target::Union{Symbol,Vector{Symbol}}`: The target parameters.
- `transformation::function`: The transformation from source to target.
"""
struct ModelMapping
    source::Union{Symbol,Vector{Symbol}}
    target::Union{Symbol,Vector{Symbol}}
    transformation::Function
end

"""
    ModelMember{T}(name, modeloptions; separate_namespace = false)

Part of `ModelOptions`. This is used to specify the options for member models. The parameter `T` should be a `ModelOptions`, but it's defined parametrically because of current limitations in Julia when dealing with circular definitions. Constructor definition is below `ModelOptions`.

## Fields
- `name::Symbol`: The name of this member. It will be used as the fieldname for the created model object.
- `modeloptions::ModelOptions`: The `ModelOptions` for the member object.
- `separate_namespace::Bool`: If `true`, all headers in the input csvs specified in `userlocations` should be prefixed with `{name}__`.
"""
struct ModelMember{T}
    name::Symbol
    modeloptions::T  # ModelOptions; Workaround for circular reference.
    separate_namespace::Bool
    # We probably also want a ParamOptions overwrite here.
end

"""
    ParamField(name, type)

Part of `ModelOptions`. This is used to specify the relevant parameters and their concrete types in the `sourceparams` and `params`. The corresponding structs will be generated to this specification.

## Fields
- `name::Symbol`: The name of the parameter field.
- `type::DataType`: Should by either SingleParam{T}, PairParam{T} or AssocParam{T}
"""
struct ParamField
    name::Symbol
    type::DataType
end

"""
    ModelOptions(args...)

A complete definition for how the model object will be created in Clapeyron. 

## Fields
- `name::Symbol`: The name of the model. A struct with this name will be generated in to the current namespace.
- `supertype::DataType`: An abstract base type to be a subtype of.
- `parent::Union{ModelOptions,Nothing} = nothing`: The `ModelOptions` of the parent if this is a variant model.
- `members::Vector{ModelMember{ModelOptions}} = ModelMember{ModelOptions}[]`: A list of modular components like `IdealModel`s.
- `locations::Vector{String} = String[]`: Default locations in Clapeyron database to look for parameters.
- `sourceparams::Vector{ParamField} = ParamField[]`: A list of relevant source parameters. The model constructor will extract the String headers according to the Symbol name given.
- `params::Vector{ParamField} = ParamField[]`: A list of relevant target parameters.
- `mappings::Vector{ModelMapping} = ModelMapping[]`: Mappings from source to target.
- `has_params::Bool = true`: Whether this model has params. A simple example of this not being the case is `BasicIdeal`.
- `has_components::Bool = true`: Whether this model is dependent on components. A simple example of this not being the case is `BasicIdeal`.
- `has_sites::Bool = false`: Whether this model has association.
- `has_groups::Bool = false`: Whether this model contains groups.
- `references::Vector{String}`: References for this model. Usually DOIs.
- `sourceparamstype::Symbol = nothing`: A struct with this name will be generated in the namespace for the source params. If given `nothing`, the constructor will fill it in with `{name}SourceParam`.
- `paramstype::Symbol = nothing`: A struct with this name will be generated in the namespace for the target params. If given `nothing`, the constructor will fill it in with `{name}Param`.
"""

struct ModelOptions
    name::Symbol
    supertype::DataType
    parent::Union{ModelOptions,Nothing}
    members::Vector{ModelMember{ModelOptions}}
    locations::Vector{String}
    sourceparams::Vector{ParamField}
    params::Vector{ParamField}
    mappings::Vector{ModelMapping}
    has_params::Bool
    has_components::Bool
    has_sites::Bool
    has_groups::Bool
    references::Vector{String}
    sourceparamstype::Symbol
    paramstype::Symbol
end

function ModelOptions(
        name::Symbol;
        supertype::DataType = Any,
        parent::Union{ModelOptions,Nothing} = nothing,
        members::Vector{ModelMember{ModelOptions}} = ModelMember{ModelOptions}[],
        locations::Vector{String} = String[],
        sourceparams::Vector{ParamField} = ParamField[],
        params::Vector{ParamField} = ParamField[],
        mappings::Vector{ModelMapping} = ModelMapping[],
        has_params::Bool = true,
        has_components::Bool = true,
        has_sites::Bool = false,
        has_groups::Bool = false,
        references::Vector{String} = String[],
        sourceparamstype::Union{Symbol,Nothing} = nothing,
        paramstype::Union{Symbol,Nothing} = nothing
    )
    return ModelOptions(
        name,
        supertype,
        parent,
        members,
        locations,
        sourceparams,
        params,
        mappings,
        has_params,
        has_components,
        has_sites,
        has_groups,
        references,
        isnothing(sourceparamstype) ? Symbol(String(name) * "SourceParam") : sourceparamstype,
        isnothing(paramstype) ? Symbol(String(name) * "Param") : paramstype,
    )
end

function ModelMember(
    name::Symbol,
    modeloptions::ModelOptions;
    separate_namespace::Bool = false
   )
    return ModelMember{ModelOptions}(name, modeloptions, separate_namespace)
end


##### Expr generator #####

function _generatecode_param_struct(
        name::Symbol,
        paramfields::Vector{ParamField}
    )::Expr
    block = Expr(:block)
    for paramfield in paramfields
        push!(block.args, :($(paramfield.name)::$(paramfield.type)))
    end
    return Expr(:struct, false, :($name <: EoSParam), block)
end

function _generatecode_model_struct(modeloptions::ModelOptions)::Expr
    # arg[2] of :block. Probably better ways to do this, but
    # I'm just running Meta.parse on a built String for now.
    defheader = String(modeloptions.name)
    if isempty(modeloptions.members)
        defheader *= ""
    else
        defheader *= "{"
        for i in 1:length(modeloptions.members)
            defheader *= "M" * string(i) * ","
        end
        defheader *= "}"
    end
    defheader *= "<:" * string(modeloptions.supertype)

    block = Expr(:block)
    push!(block.args, :(components::Vector{String}))
    if modeloptions.has_groups
        push!(block.args, :(groups::GroupParam))
    end
    if modeloptions.has_sites
        push!(block.args, :(sites::SiteParam))
    end
    if modeloptions.has_params
        push!(block.args, :(params::$(modeloptions.paramstype)))
    end
    for (i, member) in enumerate(modeloptions.members)
        push!(block.args, Expr(:(::), member.name, Symbol("M" * string(i))))
    end
    if modeloptions.has_sites
        push!(block.args, :(assoc_options::AssocOptions))
    end
    push!(block.args, :(references::Vector{String}))

    return Expr(:struct, false, Meta.parse(defheader), block)
end


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
    push!(parameters.args, Expr(:kw, :(namespace::String), :""))
    if modeloptions.has_sites
        push!(parameters.args, Expr(:kw, :(assoc_options::AssocOptions), :(AssocOptions())))
    end
    for member in modeloptions.members
        push!(parameters.args, Expr(:kw, :($(member.name)), :($(member.modeloptions.name))))
    end
    push!(func_head.args, parameters)
    # Now positional args
    if modeloptions.has_components
        push!(func_head.args, :(components::Union{String,Vector{String}}))
    else
        push!(func_head.args, Expr(:kw, :(components::Union{String,Vector{String}}), :(String[])))
    end

    # Creating function body
    block = Expr(:block)
    if modeloptions.has_params
        push!(block.args, :(params, sites = getparams(components, $(modeloptions.locations); userlocations=userlocations, verbose=verbose)))
        push!(block.args, :(packagedparams = assignparams($(modeloptions.paramstype), params, $(modeloptions.mappings))))
    end
    for member in modeloptions.members
        if member.separate_namespace
            push!(block.args, :($(member.name) = _initmodel($(member.modeloptions.name), components; userlocations, namespace=String($(member.name), verbose))))
        else
            push!(block.args, :($(member.name) = _initmodel($(member.modeloptions.name), components; userlocations, verbose)))
        end
    end
    push!(block.args, :(references = $(modeloptions.references)))

    # Generate return command
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
        push!(call.args, :(packagedparams))
    end
    for (i, member) in enumerate(modeloptions.members)
        push!(call.args, :($(member.name)))
    end
    if modeloptions.has_sites
        push!(call.args, :(assoc_options))
    end
    push!(call.args, :(references))
    push!(block.args, Expr(:return, call))
    return Expr(:function, func_head, block)
end

function _initmodel(
        model::DataType,
        components;
        userlocations,
        namespace::String = "",
        verbose::Bool = false
    )
    verbose && @info("Creating member model: $model")
    return model(components; activity, userlocations, verbose)
end

function _initmodel(
        model,
        components;
        userlocations,
        namespace::String = "",
        verbose::Bool = false
    )
    return model
end


"""
    createmodel(modeloptions; verbose)

Create the models in the global (or module) scope using `eval`.

See the tutorial or browse the implementations to see how this is used.
"""
function createmodel(modeloptions::ModelOptions; verbose::Bool = false)
    verbose && println("Generating model code for $(modeloptions.name)")
    sourceparams = _generatecode_param_struct(modeloptions.sourceparamstype, modeloptions.sourceparams)
    verbose && println(sourceparams)
    eval(sourceparams)
    params = _generatecode_param_struct(modeloptions.paramstype, modeloptions.params)
    verbose && println(params)
    eval(params)
    model = _generatecode_model_struct(modeloptions)
    verbose && println(model)
    eval(model)
    constructor = _generatecode_model_constructor(modeloptions)
    verbose && println(constructor)
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
end

export ModelMapping, ModelMember, ParamField, ModelOptions, createmodel
