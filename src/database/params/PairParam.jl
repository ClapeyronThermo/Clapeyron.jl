struct PairParameter{T,V<:AbstractMatrix{T}} <: ClapeyronDataParam
    name::String
    components::Array{String,1}
    values::V
    symmetric::Bool
    ismissingvalues::Array{Bool,2}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

"""
    PairParam{T}
Struct designed to contain pair data. used a matrix as underlying data storage.
## Creation:
```julia-repl
julia> kij = PairParam("interaction params",["water","ammonia"],[0.1 0.0;0.1 0.0])
PairParam{Float64}["water", "ammonia"]) with values:
2×2 Matrix{Float64}:
 0.1  0.0
 0.1  0.0
julia> kij.values
2×2 Matrix{Float64}:
 0.1  0.0
 0.1  0.0
julia> kij.diagvalues
2-element view(::Vector{Float64}, 
1:3:4) with eltype Float64:
 0.1
 0.0
```
## Example usage in models:
```julia
#lets compute ∑xᵢxⱼkᵢⱼ
function alpha(model,x)
    kij = model.params.kij.values
    ki = model.params.kij.diagvalues
    res = zero(eltype(molarfrac))
    for i in @comps 
        @show ki[i] #diagonal values
        for j in @comps 
            res += x[i]*x[j]*kij[i,j]
        end
    end
    return res
end
```
"""
const PairParam{T} = PairParameter{T,Matrix{T}} where T

#indexing

Base.@propagate_inbounds Base.getindex(param::PairParameter{T,<:AbstractMatrix{T}},i::Int) where T = param.values[i,i]
Base.@propagate_inbounds Base.getindex(param::PairParameter{T,<:AbstractMatrix{T}},i::Int,j::Int) where T = param.values[i,j]
Base.setindex!(param::PairParameter,val,i) = setindex!(param.values,val,i,i)
function Base.setindex!(param::PairParameter,val,i,j) 
    setindex!(param.values,val,i,j)
    param.symmetric && setindex!(param.values,val,j,i)
end

#Broadcasting

Base.broadcastable(param::PairParameter) = param.values
Base.BroadcastStyle(::Type{<:PairParameter}) = Broadcast.Style{PairParameter}()
function Base.copyto!(param::PairParameter,x)
    Base.copyto!(param.values,x)
    return param
end
Base.size(param::PairParameter) = size(param.values)

components(x::PairParameter) = x.components

PairParam(name,components,values,symmetric, missingvals,src,sourcecsv) = PairParameter(name,components,values,symmetric,missingvals,src,sourcecsv)

#unsafe constructor
function PairParam(name,components,values,symmetric = true)
    missingvals = fill(false,size(values))
    src = String[]
    sourcecsv = String[]
    return PairParameter(name,components,values,symmetric,missingvals,src,sourcecsv)
end

function PairParam(name::String,
                    components::Array{String,1},
                    values::Array{T,2},
                    symmetric = true,
                    ismissingvalues = fill(false,length(components),length(components)),
                    sourcecsvs::Array{String,1} = String[], 
                    sources::Array{String,1} = String[]) where T
    
    _values,_ismissingvalues = defaultmissing(values)
    if !all(ismissingvalues)
        _ismissingvalues = ismissingvalues
    end
    return PairParameter(name, components,_values, symmetric, _ismissingvalues, sourcecsvs, sources)
end

# If no value is provided, just initialise empty param.
function PairParam{T}(
        name::String,
        components::Vector{String};
        sources::Vector{String} = String[]
    ) where T <: AbstractString
    values = fill("", length(components), length(components))
    missingvals = fill(false, size(values))
    return PairParam(name, components, values, missingvals,true, String[], sources)
end

function PairParam{T}(
        name::String,
        components::Vector{String};
        sources::Vector{String} = String[]
    ) where T <: Number
    values = zeros(T, length(components), length(components))
    missingvals = fill(false, size(values))
    return PairParam(name, components, values, missingvals, String[], sources)
end

function PairParam(x::PairParam, name::String = x.name; isdeepcopy = true, sources = x.sources)
    if isdeepcopy
        values = deepcopy(x.values)
        symmetric = x.symmetric
        return PairParam(
            name,
            x.components,
            values,
            symmetric,
            deepcopy(x.ismissingvalues),
            x.sourcecsvs,
            sources
        )
    end
    return PairParam(
        name,
        x.components,
        x.values,
        x.symmetric,
        x.ismissingvalues,
        x.sourcecsvs,
        sources
    )
end

PairParameter(x::PairParam, name::String = x.name; isdeepcopy = true, sources = x.sources) = PairParam(x, name; isdeepcopy, sources)

function PairParam(x::SingleParam,name::String=x.name,symmetric = true)
    pairvalues = singletopair(x.values,missing)
    for i in 1:length(x.values)
        if x.ismissingvalues[i]
            pairvalues[i,i] = missing
        end
    end
    _values,_ismissingvalues = defaultmissing(pairvalues)
    return PairParam(name, x.components, _values,symmetric,_ismissingvalues,x.sourcecsvs, x.sources)
end

function Base.show(io::IO,mime::MIME"text/plain",param::PairParameter) 
    sym = param.symmetric ? "Symmetric" : ""
    print(io,sym," PairParam{",eltype(param.values),"}")
    show(io,param.components)
    println(io,") with values:")
    show(io,mime,param.values)
end

function Base.show(io::IO,param::PairParameter)
    print(io, "PairParam{",eltype(param.values),"}", "(\"", param.name,"\",)[")
    print(io,Base.summary(param.values))
    print(io,"]")
end

#convert utilities
function Base.convert(::Type{PairParam{Float64}},param::PairParam{Int})
    values = Float64.(param.values)
    return PairParam(param.name,param.components,values,param.symmetric,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{PairParam{Bool}},param::PairParam{Int})
    @assert all(z->(isone(z) | iszero(z)),param.values)
    values = Array(Bool.(param.values))
    return PairParam(param.name,param.components,values,param.symmetric,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{PairParam{Int}},param::PairParam{Float64})
    @assert all(z->isinteger(z),param.values)
    values = Int.(param.values)
    return PairParam(param.name,param.components,values,param.symmetric,param.ismissingvalues,param.sourcecsvs,param.sources)
end

#broadcasting utilities

function pack_vectors(param::PairParameter{<:AbstractVector})
    name,components,vals,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    vals = pack_vectors(vals)
    return PairParam(name,components,vals,param.symmetric,missingvals,srccsv,src)
end

const PackedSparsePairParam{T} = Clapeyron.PairParameter{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}, SparsePackedMofV{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, 
true}, PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{T}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}}} where T

# Operations

function Base.:(+)(param::PairParameter, x::Number)
    values = param.values .+ x
    return PairParam(param.name, param.components, values, param.symmetric, param.ismissingvalues, param.sourcecsvs, param.sources)
end

function Base.:(*)(param::PairParameter, x::Number)
    values = param.values .* x
    return PairParam(param.name, param.components, values, param.symmetric, param.ismissingvalues, param.sourcecsvs, param.sources)
end

function Base.:(^)(param::PairParameter, x::Number)
    values = param.values .^ x
    return PairParam(param.name, param.components, values, param.symmetric, param.ismissingvalues, param.sourcecsvs, param.sources)
end
