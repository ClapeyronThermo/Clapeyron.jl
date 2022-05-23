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
        #while the decorative params.name are nice, for accessing, the field name is more useful
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
function inputparams!(model,input,output)
    calculated_output = output(input)
    copyto!(output,calculated_output)
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

function Base.show(io::IO,param::Union{SingleParameter,PairParameter})
    print(io, typeof(param), "(\"", param.name, "\")")
    show(io,param.components)  
end

include("params/AssocParam.jl")
include("params/GroupParam.jl")
include("params/SiteParam.jl")
include("params/AssocOptions.jl")
include("params/combiningrules_base.jl") #general combining Rules for Single and Pair Params.
include("params/combiningrules.jl") #specific rules.

export SingleParam, SiteParam, PairParam, AssocParam, GroupParam
export AssocOptions

"""
    diagvalues(x)

A common function to retrieve the main diagonal values that work on both SingleParam and PairParam.
"""
function diagvalues(x::SingleParam)
    return x.values
end

function diagvalues(x::PairParam)
    return view(x.values, diagind(x.values))
end

function _get_sources(x::Vector)::Vector{String}
    return collect(Set(r for y ∈ x for r ∈ y.sources))
end

function _get_sources(x)::Vector{String}
    return copy(x.sources)
end
