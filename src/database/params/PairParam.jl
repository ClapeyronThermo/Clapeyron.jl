struct PairParameter{T,V<:AbstractMatrix{T},D} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::V
    diagvalues::D
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
const PairParam{T} = PairParameter{T,Matrix{T},SubArray{T, 1, Vector{T}, Tuple{StepRange{Int64, Int64}}, true}} where T

PairParam(name,components,values,diagvals, missingvals,src,sourcecsv) = PairParameter(name,components,values,diagvals,missingvals,src,sourcecsv)

#unsafe constructor
function PairParam(name,components,values)
    missingvals = fill(false,size(values))
    diagvals = view(values, diagind(values))
    src = String[]
    sourcecsv = String[]
    return PairParam(name,components,values,diagvals, missingvals,src,sourcecsv)
end
function PairParam(name::String,
                    components::Array{String,1},
                    values::Array{T,2},
                    ismissingvalues = fill(false,length(components),length(components)),
                    sourcecsvs::Array{String,1} = String[], 
                    sources::Array{String,1} = String[]) where T
    
    _values,_ismissingvalues = defaultmissing(values)
    diagvalues = view(_values, diagind(_values))
    if !all(ismissingvalues)
        _ismissingvalues = ismissingvalues
    end
    return PairParam(name, components,_values, diagvalues, _ismissingvalues, sourcecsvs, sources)
end

function PairParam(x::PairParameter,name::String=x.name)
    values = deepcopy(x.values)
    diagvalues = view(values,diagind(values))
    return PairParam(name, x.components,values ,diagvalues, deepcopy(x.ismissingvalues), x.sourcecsvs, x.sources)
end

function PairParam(x::SingleParameter,name::String=x.name)
    pairvalues = singletopair(x.values,missing)
    for i in 1:length(x.values)
        if x.ismissingvalues[i]
            pairvalues[i,i] = missing
        end
    end
    _values,_ismissingvalues = defaultmissing(pairvalues)
    diagvalues = view(_values, diagind(_values))
    return PairParam(name, x.components, _values,diagvalues,_ismissingvalues,x.sourcecsvs, x.sources)
end
#=
function PairParam(x::PairParam, v::Matrix,name::String=x.name)
    return PairParam(name, x.components,deepcopy(v), x.sourcecsvs, x.sources)
end
function PairParam(x::SingleParam, v::Vector,name::String=x.name)
    pairvalues = singletopair(v,missing)
    return PairParam(x.name, x.components, pairvalues,x.sourcecsvs, x.sources)
end
function PairParam(x::SingleParam, v::Matrix,name::String=x.name)
    return PairParam(x.name, x.components, deepcopy(v),x.sourcecsvs, x.sources)
end
=#
#barebones constructor by list of pairs.

#convert utilities
function Base.show(io::IO,mime::MIME"text/plain",param::PairParameter) 
    print(io,"PairParam{",eltype(param.values),"}")
    show(io,param.components)
    println(io,") with values:")
    show(io,mime,param.values)
end

function Base.show(io::IO,param::PairParameter)
    print(io, "PairParam{",eltype(param.values),"}", "(\"", param.name, "\")[")
    print(io,Base.summary(param.values))
    print(io,"]")
end

function Base.convert(::Type{PairParam{Float64}},param::PairParam{Int})
    values = Float64.(param.values)
    diagvalues = view(values, diagind(values))
    return PairParam(param.name,param.components,values,diagvalues,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{PairParam{Bool}},param::PairParam{Int})
    @assert all(z->(isone(z) | iszero(z)),param.values)
    values = Array(Bool.(param.values))
    diagvalues = view(values, diagind(values))
    return PairParam(param.name,param.components,values,diagvalues,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{PairParam{Int}},param::PairParam{Float64})
    @assert all(z->isinteger(z),param.values)
    values = Int.(param.values)
    diagvalues = view(values, diagind(values))
    return PairParam(param.name,param.components,values,diagvalues,param.ismissingvalues,param.sourcecsvs,param.sources)
end

#broadcasting utilities
Base.broadcastable(param::PairParameter) = param.values

function pack_vectors(param::PairParameter{<:AbstractVector})
    name,components,vals,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    vals = pack_vectors(vals)
    return PairParam(name,components,vals,nothing,missingvals,srccsv,src)
end

const PackedSparsePairParam{T} = Clapeyron.PairParameter{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}, SparsePackedMofV{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, 
true}, PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{T}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}}, Nothing} where T

