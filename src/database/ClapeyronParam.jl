"""
    ClapeyronParam

Abstract type corresponding to a Clapeyron parameter.

it has to be splittable (via [`split_model`](@ref)) and have a `components` field
"""
abstract type ClapeyronParam end
abstract type EoSParam end
export EoSParam

const PARSED_GROUP_VECTOR_TYPE =  Vector{Tuple{String, Vector{Pair{String, Int64}}}}

function pack_vectors(x::AbstractVector{<:AbstractVector})
    return PackedVectorsOfVectors.pack(x)
end

function pack_vectors(x::SparseMatrixCSC{<:AbstractVector})
    return SparsePackedMofV(x)
end

include("params/paramvectors.jl")
include("params/SingleParam.jl")
include("params/PairParam.jl")
include("params/AssocParam.jl")
include("params/GroupParam.jl")
include("params/SiteParam.jl")
include("params/AssocOptions.jl")

export SingleParam, SiteParam, PairParam, AssocParam, GroupParam
export AssocOptions


"""
    assignparams(paramstruct, parms, mappings)

    Returns the constructed EoSParam struct using defined transformations in the `mappings` vector. If a field is not present in a target of one of the mappings, look for the equivalent name in the keys of the raw `params` dict. If it is also not found there, throw an error because the EoSParam cannot be constructed.

    # Arguments
    - `paramstruct::DataType`: The `EoSParam` type.
    - `params::Dict`: The Dict of parameters returned from `getparams`.
    - `mappings::Vector{Tuple{Union{String,Vector{String}},Symbol,Function}}`: A vector of ordered triplets: the of the key(s) in params, the target Symbol, and a function to apply mapping.
"""
function assignparams(
        paramstruct::DataType,
        params::Dict,
        mappings::Vector{Tuple{Any,Symbol,Function}}
    )
    rawnames = collect(keys(params))
    mappednames = [mapping[2] for mapping ∈ mappings]
    # List to be splatted into the paramstruct constructor
    transformedparams = []
    for fieldname ∈ fieldnames(paramstruct)
        mappednameindex = findfirst(isequal(fieldname), mappednames)
        if isnothing(mappednameindex)
            # If it's not in target mapping, take from raw params.
            String(fieldname) ∉ rawnames && error(fieldname, " not found.")
            push!(transformedparams, params[String(fieldname)])
        else
            # argnames can be a String or Vector{String}
            argnames = mappings[mappednameindex][1]
            args = []
            if typeof(argnames) <: AbstractString
                push!(args, params[argnames])
            else
                for argname in argnames
                    push!(args, params[argname])
                end
            end
            push!(transformedparams, mappings[mappednameindex][3](args...))
        end
    end
    return paramstruct(transformedparams...)
end
export assignparams
