"""
    ClapeyronParam
Abstract type corresponding to a Clapeyron parameter.
it has to be splittable (via [`split_model`](@ref)) and have a `components` field
"""
abstract type ClapeyronParam end
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

function Base.show(io::IO, params::EoSParam)
    print(io, typeof(params))
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
include("params/AssocParam.jl")
include("params/GroupParam.jl")
include("params/SiteParam.jl")
include("params/AssocOptions.jl")

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