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


