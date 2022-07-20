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
julia> diagvalues(kij)
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
    ki = diagvalues(model.params.kij)
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

#indexing

Base.@propagate_inbounds Base.getindex(param::PairParameter{T,<:AbstractMatrix{T}},i::Int) where T = param.values[i,i]
Base.@propagate_inbounds Base.getindex(param::PairParameter{T,<:AbstractMatrix{T}},i::Int,j::Int) where T = param.values[i,j]
Base.setindex!(param::PairParameter,val,i) = setindex!(param.values,val,i,i,false)
function Base.setindex!(param::PairParameter,val,i,j,symmetric = true) 
    setindex!(param.values,val,i,j)
    symmetric && setindex!(param.values,val,j,i)
end

#Broadcasting

Base.broadcastable(param::PairParameter) = param.values
Base.BroadcastStyle(::Type{<:PairParameter}) = Broadcast.Style{PairParameter}()

#copyto!

function Base.copyto!(param::PairParameter,x)
    Base.copyto!(param.values,x)
    return param
end

function Base.copyto!(dest::PairParameter,src::PairParameter) #used to set params
    #key check
    dest.components == src.components || throw(DimensionMismatch("components of source and destination pair parameters are not the same for $dest"))
    
    #=
    TODO: it does not check that dest.symmetric = src.symmetric, the only solution i see at the moment 
    is to make the Single, Pair and Assoc Params, mutable structs, but i don't really know the performance
    implications of that.
    supposedly, this copyto! is only used internally, and both src and dest params are already the same. but it would
    be good to enforce that.
    =#
    copyto!(dest.values,src.values)
    dest.ismissingvalues .= src.ismissingvalues
    return dest
end

Base.size(param::PairParameter) = size(param.values)

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

# If no value is provided, just initialise empty param.
function PairParam(
        name::String,
        components::Vector{String};
        sources::Vector{String} = String[]
    ) where T <: AbstractString
    values = fill(0.0, length(components), length(components))
    missingvals = fill(false, size(values))
    return PairParam(name, components, values, missingvals, String[], sources)
end

function PairParam(x::PairParam, name::String = x.name; isdeepcopy = true, sources = x.sources)
    if isdeepcopy
        values = deepcopy(x.values)
        diagvalues = view(values, diagind(values))
        return PairParam(
            name,
            x.components,
            values,
            diagvalues,
            deepcopy(x.ismissingvalues),
            x.sourcecsvs,
            sources
        )
    end
    return PairParam(
        name,
        x.components,
        x.values,
        x.diagvalues,
        x.ismissingvalues,
        x.sourcecsvs,
        sources
    )
end

PairParameter(x::PairParam, name::String = x.name; isdeepcopy = true, sources = x.sources) = PairParam(x, name; isdeepcopy, sources)

function PairParam(x::SingleParam,name::String=x.name)
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

function Base.show(io::IO,mime::MIME"text/plain",param::PairParameter) 
    #sym = param.symmetric ? "Symmetric " : ""
    sym = ""
    _size = size(param)
    _size_str = string(_size[1]) * "×" * string(_size[2]) * " "
    print(io,sym,_size_str,"PairParam{",eltype(param.values),"}(")
    show(io,param.components)
    println(io,") with values:")
    Base.print_matrix(IOContext(io, :compact => true),param.values)
end

#convert utilities
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

function pack_vectors(param::PairParameter{<:AbstractVector})
    name,components,vals,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    vals = pack_vectors(vals)
    return PairParam(name,components,vals,nothing,missingvals,srccsv,src)
end

const PackedSparsePairParam{T} = Clapeyron.PairParameter{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}, SparsePackedMofV{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, 
true}, PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{T}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}}, Nothing} where T

# Operations (seems that they are not needed)
#=
function Base.:(+)(param::PairParameter, x::Number)
    values = param.values .+ x
    return PairParam(param.name, param.components, values, param.ismissingvalues, param.sourcecsvs, param.sources)
end

function Base.:(*)(param::PairParameter, x::Number)
    values = param.values .* x
    return PairParam(param.name, param.components, values, param.diagvalues, param.ismissingvalues, param.sourcecsvs, param.sources)
end

function Base.:(^)(param::PairParameter, x::Number)
    values = param.values .^ x
    return PairParam(param.name, param.components, values, param.diagvalues, param.ismissingvalues, param.sourcecsvs, param.sources)
end=#