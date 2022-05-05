"""
    ModelMapping(source, target, transformation)

Part of `ModelOptions`. Used to define how source from csv is mapped to the target struct. The `source` and `target` should be defined in `inputparams` and `params` respecfully in `ModelOptions`.

## Fields
- `source::Union{String,Vector{String}}`: The source headers from csv.
- `target::Union{Symbol,Vector{Symbol}}`: The target parameters.
- `transformation::function`: The transformation from source to target.
"""
struct ModelMapping
    source::Vector{Symbol}
    target::Vector{Symbol}
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

Part of `ModelOptions`. This is used to specify the relevant parameters and their concrete types in the `inputparams` and `params`. The corresponding structs will be generated to this specification.

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
- `name::Symbol`: The name of the model. A struct with this name will be generated in to the global/module namespace.
- `supertype::DataType`: An abstract base type to be a subtype of.
- `parent::Union{ModelOptions,Nothing} = nothing`: The `ModelOptions` of the parent if this is a variant model.
- `members::Vector{ModelMember{ModelOptions}} = ModelMember{ModelOptions}[]`: A list of modular components like `IdealModel`s.
- `locations::Vector{String} = String[]`: Default locations in Clapeyron database to look for parameters. Note that for user-specified parameters, it is easier to use the `user_locations` parameter.
- `inputparams::Vector{ParamField} = ParamField[]`: A list of relevant source parameters. The model constructor will extract the String headers according to the Symbol name given.
- `params::Vector{ParamField} = ParamField[]`: A list of relevant target parameters.
- `mappings::Vector{ModelMapping} = ModelMapping[]`: Mappings from source to target.
- `has_params::Bool = true`: Whether this model has params. A simple example of this not being the case is `BasicIdeal`.
- `has_components::Bool = true`: Whether this model is dependent on components. A simple example of this not being the case is `BasicIdeal`.
- `has_sites::Bool = false`: Whether this model has association.
- `has_groups::Bool = false`: Whether this model contains groups.
- `references::Vector{String}`: References for this model. Usually DOIs.
- `inputparamstype::Symbol = nothing`: A struct with this name will be generated in the global/module namespace for the input params. If given `nothing`, the constructor will fill it in with `{name}InputParam`.
- `paramstype::Symbol = nothing`: A struct with this name will be generated in the global/module namespace for the target params. If given `nothing`, the constructor will fill it in with `{name}Param`.
"""

struct ModelOptions
    name::Symbol
    supertype::DataType
    parent::Union{ModelOptions,Nothing}
    members::Vector{ModelMember{ModelOptions}}
    locations::Vector{String}
    inputparams::Vector{ParamField}
    params::Vector{ParamField}
    mappings::Vector{ModelMapping}
    has_params::Bool
    has_components::Bool
    has_sites::Bool
    has_groups::Bool
    references::Vector{String}
    inputparamstype::Symbol
    paramstype::Symbol
end

function ModelOptions(
        name::Symbol;
        supertype::DataType = Any,
        parent::Union{ModelOptions,Nothing} = nothing,
        members::Vector{ModelMember{ModelOptions}} = ModelMember{ModelOptions}[],
        locations::Vector{String} = String[],
        inputparams::Vector{ParamField} = ParamField[],
        params::Vector{ParamField} = ParamField[],
        mappings::Vector{ModelMapping} = ModelMapping[],
        has_params::Bool = true,
        has_components::Bool = true,
        has_sites::Bool = false,
        has_groups::Bool = false,
        references::Vector{String} = String[],
        inputparamstype::Union{Symbol,Nothing} = nothing,
        paramstype::Union{Symbol,Nothing} = nothing
    )
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
        references,
        isnothing(inputparamstype) ? Symbol(String(name) * "InputParam") : inputparamstype,
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
    # arg[3] of :struct. Probably better ways to do this, but
    # I'm just running Meta.parse on a built String for now.
    defheader = String(modeloptions.name)
    if isempty(modeloptions.members)
        defheader *= ""
    else
        defheader *= "{"
        for i ∈ 1:length(modeloptions.members)
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
        push!(block.args, :(inputparams::$(modeloptions.inputparamstype)))
        push!(block.args, :(params::$(modeloptions.paramstype)))
        push!(block.args, :(mappings::Vector{ModelMapping}))
    end
    for (i, member) ∈ enumerate(modeloptions.members)
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
    for member ∈ modeloptions.members
        push!(parameters.args, Expr(:kw, :($(member.name)), :($(member.modeloptions.name))))
    end
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
        push!(block.args, :((rawparams, sites) = getparams(components, $(modeloptions.locations); userlocations=userlocations, verbose=verbose)))
        push!(block.args, :((inputparams, params) = _initparams($(modeloptions.inputparamstype), $(modeloptions.paramstype), rawparams, mappings)))
    end
    for member ∈ modeloptions.members
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
    push!(block.args, Expr(:return, call))
    return Expr(:function, func_head, block)
end

#####

function _initmodel(
        model::DataType,
        components;
        userlocations,
        namespace::String = "",
        verbose::Bool = false
    )
    verbose && @info("Creating member model: $model")
    return model(components; userlocations, verbose)
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
    _initparams(::Type{I}, ::Type{P}, rawparams, mappings)

Returns the constructed InputParamType and ParamType based on the rawparams and mappings.

Strategy:

- If a mapping is present, and it is not an identity transformation, construct a new param struct with the output param name.
- If a mapping is present, and it is an identity transformation, then it is simply a name change. Create a new param struct, but the `value` is a reference to the inputparam array.
- If no mapping is present, then the parameter is taken as-is from the raw/input params. Just point to the same param object.

# Arguments
- `Type{InputParamType}`: The type for the inputparm struct.
- `Type{ParamType`: The type for the param struct.
- `rawparams::Dict`: The dict returned from `getparams`.
- `mappings::Vector{ModelMapping}`: The list of mappings.

# Returns
Tuple{InputParamType,ParamType}
"""
function _initparams(
        ::Type{I}, 
        ::Type{P},
        rawparams::Dict,
        mappings::Vector{ModelMapping}
    )::Tuple{I,P} where {I,P}
    inputparams_names = collect(fieldnames(I))
    inputparams_types = collect(fieldtypes(I))
    inputparams_params = [rawparams[String(param)] for param ∈ inputparams_names]
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
            # Arrays will only have one elements
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
                # But still reference the same vallues array.
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
    inputparams = I(inputparams_params...)
    params = P(params_params...)
    updateparams!(inputparams, params, mappings)
    return (inputparams, params)
end

"""
    updateparams!(model)
    updateparams!(inputparams, params, mappings)

Given an `EoSModel`, update the values of `params` with current
values of `inputparams` given the `mappings`. Note that for
identity transformations, the values share the same references,
so they need not be updated explictly.
"""
function updateparams!(model::EoSModel)::Nothing
    updateparams!(model.inputparams, model.params, model.mappings)
end

function updateparams!(
        inputparams,
        params,
        mappings::Vector{ModelMapping}
       )::Nothing
    # At the moment, we are stil constructing new ClapeyronParams,
    # then broadcast assigning the params to them. There might be
    # some optimisations possible here, but this will maintain
    # maximum compatibility with existing functions for now.
    for mapping in filter(m -> !(m.transformation === identity),  mappings)
        outputs = mapping.transformation([getfield(inputparams, f) for f in mapping.source]...)
        if typeof(outputs) <: ClapeyronParam
            outputs = [outputs]
        end
        toupdates = [getfield(params, f) for f in mapping.target]
        for (output, toupdate) in zip(outputs, toupdates)
            toupdate.values .= output.values
        end
    end
end

"""
    createmodel(modeloptions; verbose)

Create the models in the global (or module) scope using `eval`.

See the tutorial or browse the implementations to see how this is used.
"""
function createmodel(modeloptions::ModelOptions; verbose::Bool = false)
    verbose && @info("Generating model code for $(modeloptions.name)")
    inputparams = _generatecode_param_struct(modeloptions.inputparamstype, modeloptions.inputparams)
    verbose && @info(inputparams)
    eval(inputparams)
    params = _generatecode_param_struct(modeloptions.paramstype, modeloptions.params)
    verbose && @info(params)
    eval(params)
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
end

export ModelMapping, ModelMember, ParamField, ModelOptions, createmodel, updateparams!
