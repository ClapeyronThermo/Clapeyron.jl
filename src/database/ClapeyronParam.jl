"""
    ClapeyronParam
Abstract type corresponding to a Clapeyron parameter.
it has to be splittable (via [`split_model`](@ref)) and have a `components` field
"""
abstract type ClapeyronParam end
abstract type EoSParam end
export EoSParam

custom_show(param::EoSParam) = _custom_show_param(typeof(param))

function _custom_show_param(::Type{T}) where T <: EoSParam
    types = fieldtypes(T)
    return all(x -> x <: ClapeyronParam,types)
end

function Base.show(io::IO, mime::MIME"text/plain", params::EoSParam)
    !custom_show(params) && return show_default(io,mime,params)
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

function build_eosparam(::Type{T},data) where T <: EoSParam
    names = fieldnames(T)
    return T((data[string(name)] for name in names)...)
end

Base.eltype(p::EoSParam) = Float64

const PARSED_GROUP_VECTOR_TYPE =  Vector{Tuple{String, Vector{Pair{String, Int64}}}}

function pack_vectors(x::AbstractVector{<:AbstractVector})
    return PackedVectorsOfVectors.pack(x)
end

function pack_vectors(x::SparseMatrixCSC{<:AbstractVector})
    return SparsePackedMofV(x)
end

function param_length_check(paramtype,name,comp_length,val_length)
    if comp_length != val_length
        throw(DimensionMismatch(string(paramtype) * "(\"$(name)\"): expected length of components ($comp_length) equal to component length in values ($val_length)"))
    end
end

function _str_to_idx(param,i,j)
    idx_i = findfirst(isequal(i),param.components)
    idx_j = findfirst(isequal(j),param.components)
    if idx_i === nothing && idx_j === nothing
        throw(BoundsError(param,(-1,-1)))
    elseif idx_i === nothing
        throw(BoundsError(param,(-1,idx_j)))
    elseif idx_j === nothing
        throw(BoundsError(param,(idx_i,-1)))
    else
       return (idx_i::Int,idx_j::Int)
    end
end

function _str_to_idx(param,i)
    idx = findfirst(isequal(i),param.components)
    isnothing(idx) && throw(BoundsError(param,-1))
    return idx::Int
end

include("params/paramvectors.jl")
include("params/SingleParam.jl")
include("params/PairParam.jl")
include("params/AssocParam.jl")
include("params/GroupParam.jl")
include("params/SiteParam.jl")
include("params/AssocOptions.jl")
include("params/SpecialComp.jl")
include("params/ReferenceState.jl")

const SingleOrPair = Union{<:SingleParameter,<:PairParameter}
function Base.show(io::IO,param::SingleOrPair)
    print(io, typeof(param), "(\"", param.name, "\")")
    show(io,param.components)
end

#internal utility function
#shortcut for model.params.val, but returns nothing if the val is not found.
@pure function getparam(model::EoSModel,val::Symbol)
    M = typeof(model)
    if hasfield(M,:params)
        if hasfield(typeof(model.params),val)
            return getfield(model.params,val)
        end
    end
    return nothing
end

Base.iterate(param::SingleOrPair) = iterate(param.values) 
Base.iterate(param::SingleOrPair,state) = iterate(param.values,state)

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

function diagvalues(x::Number)
    return x
end

function _get_sources(x::Vector)::Vector{String}
    return collect(Set(r for y ∈ x for r ∈ y.sources))
end

function _get_sources(x)::Vector{String}
    return copy(x.sources)
end