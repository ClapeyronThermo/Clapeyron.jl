"""
    ClapeyronParam
Abstract type corresponding to a Clapeyron parameter.
It requires to define `Base.eltype`, and specify how it is splitted (via defining the adecuate `each_split_model`  method or by defining `Clapeyron.is_splittable(param) = false` to mark as non-splittable)
"""
abstract type ClapeyronParam end

"""
    EoSParam
Abstract type corresponding to a container of `ClapeyronParam`s.
it supposes that all fields are `ClapeyronParam`s.
"""
abstract type EoSParam end
abstract type ParametricEoSParam{T} <: EoSParam end


"""
    OptionsParam <: ClapeyronParam
Abstract type corresponding to a Clapeyron parameter that only contains options.
It is assumed that this parameter is equal to all components.
"""
abstract type OptionsParam <: ClapeyronParam end
Base.eltype(param::OptionsParam) = Bool
is_splittable(::OptionsParam) = false

export EoSParam, ParametricEoSParam

paramtype(m::ClapeyronParam) = eltype(m)
paramtype(::Type{M}) where M <: ClapeyronParam = eltype(m)

custom_show(param::EoSParam) = _custom_show_param(typeof(param))


function build_parametric_param(param::Type{T}, args...) where T <: ParametricEoSParam
    return __build_parametric_param(param,args)
end

#dynamic version
function dyn_build_parametric_param(param::Type{T}, args) where T <: ParametricEoSParam
    TT = mapreduce(paramtype,promote_type,args)
    Param = parameterless_type(T)
    vals = map(Base.Fix1(_convert_param,TT),args)
    return Param{TT}(vals...)
end

#static version
@generated function __build_parametric_param(param::Type{P},args::A) where {P <: ParametricEoSParam,A}
    argtypes = fieldtypes(A)
    TT = mapreduce(paramtype,promote_type,argtypes)
    paramname = Base.typename(P).name
    paramlength = length(fieldnames(P))
    expr = :($paramname{$TT}())
    args = expr.args
    for i in 1:paramlength
        push!(args,:(_convert_param($TT,args[$i])))
    end
    return expr
end

function _custom_show_param(::Type{T}) where T <: EoSParam
    types = fieldtypes(T)
    return all(x -> x <: ClapeyronParam,types)
end

Solvers.primalval(x::EoSParam) = Solvers.primalval_struct(x)

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
    params = map(name -> data[string(name)],names)
    return T(params...)
end

function build_eosparam(::Type{T},data) where T <: ParametricEoSParam
    names = fieldnames(T)
    params = map(name -> data[string(name)],names)
    dyn_build_parametric_param(T,params)
end

Base.eltype(p::EoSParam) = Float64
Base.eltype(p::ParametricEoSParam{T}) where T = T

const PARSED_GROUP_VECTOR_TYPE = Vector{Tuple{String, Vector{Pair{String, Int64}}}}

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


function _convert_param(T::V,val) where V
    return _convert_param(T,parameterless_type(val),val)
end

function _convert_param(T::Type{V},val::SpecialComp) where V
    return val
end

function _convert_param(T::Type{V},val::ReferenceState) where V
    return val
end

function _convert_param(T::Type{V},::Type{SingleParameter},val::SingleParameter{V}) where V
    return val
end

function _convert_param(T::Type{V},::Type{PairParameter},val::PairParameter{V}) where V
    return val
end

function _convert_param(T::Type{V},::Type{AssocParam},val::AssocParam{V}) where V
    return val
end

function _convert_param(T::Type{V},::Type{MixedGCSegmentParam},val::MixedGCSegmentParam{V}) where V
    return val
end

function _convert_param(T::Type{V},::Type{SingleParameter},val) where {V}
    return convert(SingleParam{T},val)
end

function _convert_param(T::Type{V},::Type{PairParameter},val) where {V}
    return convert(PairParam{T},val)
end

function _convert_param(T::Type{V},::Type{AssocParam},val) where {V}
    return convert(AssocParam{T},val)
end

function _convert_param(T::Type{V},::Type{MixedGCSegmentParam},val) where {V}
    return convert(MixedGCSegmentParam{T},val)
end

const SingleOrPair = Union{<:SingleParameter,<:PairParameter}
function Base.show(io::IO,param::SingleOrPair)
    print(io, typeof(param), "(\"", param.name, "\")")
    show(io,param.components)
end

#internal utility function
#shortcut for model.params.val, but returns nothing if the val is not found.
Base.@assume_effects :foldable function getparam(model::EoSModel,val::Symbol)
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
    return view(x, linearidx(x))
end

function diagvalues(x::AbstractVector)
    return view(x, linearidx(x))
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