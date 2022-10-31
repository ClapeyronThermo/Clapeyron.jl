"""
    ClapeyronParam
Abstract type corresponding to a Clapeyron parameter.
it has to be splittable (via [`split_model`](@ref)) and have a `components` field
"""
abstract type ClapeyronParam end

"""
    ClapeyronDataParam
Abstract type corresponding to a Clapeyron data parameter.
it has to be splittable (via [`split_model`](@ref)), have a `components` field and have stored values
"""
abstract type ClapeyronDataParam <: ClapeyronParam end

abstract type EoSParam end
export EoSParam

function Base.show(io::IO, mime::MIME"text/plain", params::EoSParam)
    names = fieldnames(typeof(params))
    if length(names) == 1
        print(io, typeof(params), " for ", getfield(params, first(names)).components, " with ", length(names), " param:")
    else
        print(io, typeof(params), " for ", getfield(params, first(names)).components, " with ", length(names), " params:")
    end
    for name in names
        param = getfield(params, name)
        print(io, "\n ", name, "::", typeof(param))
    end
end
function param(::Type{T},dict::Dict{String,ClapeyronParam}) where T <: EoSParam
    return T((dict[string(name)] for name in fieldnames(T))...)
end

function Base.show(io::IO, params::EoSParam)
    print(io, typeof(params))
end

function Base.copyto!(dest::T,src::T) where T<:EoSParam
    for name in fieldnames(T)
        copyto!(getfield(dest,name),gerfield(dest,name))
    end
    return dest
end

#general fallback that only depends on the the definition of Param(x::InputParam)
#it can be overloaded to use implace values, but it is a performance consideration
function updateparams2!(model,input,output)
    model === nothing && return model #if the model is nothing, return unchanged
    components(model) === nothing && return model #if a model doesn't have components, shouldn't be updated!
    T = typeof(model)
    types = fieldtypes(T)
    #default behaviour, return the same model if we are unable to know if we can update
    if !(typeof(input) ∈ types & typeof(output) ∈ types)
        @warn("cannot update $T with input parameters $input and output parameters $output")
        return model
    end
    output_f = parameterless_type(output) #in utils/core_utils.jl
    calculated_output = output_f(input) #calculate output from input (user defined).
    copyto!(output,calculated_output) #copy the results to the output,we suppose output = model.params, but it could not.
    return model
end

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

const SingleOrPair = Union{<:SingleParameter,<:PairParameter}
function Base.show(io::IO,param::SingleOrPair)
    print(io, typeof(param), "(\"", param.name, "\")")
    show(io,param.components)  
end

include("params/AssocParam.jl")
include("params/GroupParam.jl")
include("params/SiteParam.jl")
include("params/AssocOptions.jl")
include("params/combiningrules_base.jl") #general combining Rules for Single and Pair Params.
include("params/combiningrules_implace.jl") #implace specific rules.
include("params/combiningrules.jl") #out of place versions

const SingleOrPair = Union{<:SingleParameter,<:PairParameter}
function Base.show(io::IO,param::SingleOrPair)
    print(io, typeof(param), "(\"", param.name, "\")")
    show(io,param.components)
end

export SingleParam, SiteParam, PairParam, AssocParam, GroupParam
export AssocOptions

"""
    diagvalues(x)

A common function to retrieve the main diagonal values that work on both SingleParam and PairParam.
"""
function diagvalues end

function diagvalues(x::AbstractMatrix)
    return view(x, diagind(x))
end

function diagvalues(x::AbstractVector)
    return x
end

function diagvalues(x::SingleOrPair)
    return diagvalues(x.values)
end

function _get_sources(x::Vector)::Vector{String}
    return collect(Set(r for y ∈ x for r ∈ y.sources))
end

function _get_sources(x)::Vector{String}
    return copy(x.sources)
end
