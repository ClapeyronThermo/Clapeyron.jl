struct PairParameter{T,V<:AbstractMatrix{T}} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::V
    ismissingvalues::Array{Bool,2}
    sourcecsvs::Union{Vector{String},Nothing}
    sources::Union{Vector{String},Nothing}
end

#methods to fill missing sources/sourcescsvs
PairParameter(name,components,values,ismissingvalues) = PairParameter(name,components,values,ismissingvalues,nothing,nothing)
PairParameter(name,components,values,ismissingvalues,src) = PairParameter(name,components,values,ismissingvalues,src,nothing)

"""
    PairParam{T}
Struct designed to contain pair data. used a matrix as underlying data storage.
## Creation:
```julia-repl
julia> kij = PairParam("interaction params",["water","ammonia"],[0.1 0.0;0.1 0.0])
PairParam{Float64}(["water", "ammonia"]) with values:
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
const PairParam{T} = PairParameter{T,Matrix{T}} where T
PairParam(name,components,vals,missingvals,srccsv,src) = PairParameter(name,components,vals,missingvals,srccsv,src)

#indexing
Base.@propagate_inbounds Base.getindex(param::PairParameter{T},i::Int) where T = param.values[i,i]
Base.@propagate_inbounds Base.getindex(param::PairParameter{T},i::Int,j::Int) where T = param.values[i,j]

function Base.getindex(param::PairParameter{T},i::AbstractString) where T
    idx = _str_to_idx(param,i)
    return param[idx]
end

function Base.getindex(param::PairParameter{T},i::AbstractString,j::AbstractString) where T
    idx_i,idx_j = _str_to_idx(param,i,j)
    return param[idx_i,idx_j]
end

Base.setindex!(param::PairParameter,val,i) = setindex!(param,val,i,i,false)

function Base.setindex!(param::PairParameter,val,i,j,symmetric = true) 
    setindex!(param.values,val,i,j)
    symmetric && setindex!(param.values,val,j,i)
end

function Base.setindex!(param::PairParameter,val,i::AbstractString,j::AbstractString,symmetric = true) 
    idx_i,idx_j = _str_to_idx(param,i,j)
    setindex!(param.values,val,idx_i,idx_j)
    symmetric && setindex!(param,val,idx_j,idx_i)
end

function Base.setindex!(param::PairParameter,val,i::AbstractString)
    idx = _str_to_idx(param,i)
    setindex!(param,val,idx::Int)
end

#Broadcasting

Base.broadcastable(param::PairParameter) = param.values
Base.BroadcastStyle(::Type{<:PairParameter}) = Broadcast.Style{PairParameter}()

#copyto!

function Base.copyto!(dest::PairParameter,src::Base.Broadcast.Broadcasted)
    Base.copyto!(dest.values,src)
    return dest
end

function Base.copyto!(dest::PairParameter,src::AbstractArray)
    Base.copyto!(dest.values,src)
    return dest
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

Base.eltype(param::PairParameter{T}) where T = T
Base.eltype(param::Type{<:PairParameter{T}}) where T = T

#primalval
function Solvers.primalval(x::PairParameter)
    return PairParameter(x.name,x.components,Solvers.primalval_eager(x.values),x.ismissingvalues,x.sourcecsvs,x.sources)
end

#barebones constructor, we provide vals and missing vals
function PairParam(name,components,values::Matrix{T},missingvals,src,sourcecsv) where T 
    val_length = LinearAlgebra.checksquare(values)
    param_length_check(PairParam,name,length(components),val_length)
    PairParameter{T,Matrix{T}}(name,components,values,missingvals,src,sourcecsv)
end

function PairParam(name,components,values::Vector{T},missingvals::Vector,src,sourcecsv) where T 
    n = length(values)
    mat_values = zeros(T,(n,n))
    mat_missing = ones(Bool,(n,n))
    for i in 1:n
        mat_values[i,i] = values[i]
        mat_missing[i,i] = false
    end
    PairParam(name,components,mat_values,mat_missing,src,sourcecsv)
end

function PairParam(name,components,values::AbstractMatrix{T},missingvals,src,sourcecsv) where T 
    return PairParam(name,components,convert(Matrix{T},values),missingvals,src,sourcecsv)
end

PairParam(name,components,values,missingvals,src) = PairParam(name,components,values,missingvals,src,nothing)
PairParam(name,components,values,missingvals) = PairParam(name,components,values,missingvals,nothing,nothing)

#constructor in case we provide a normal vector
function PairParam(name, components, values_or_missing::AbstractMatrix{T}) where T
    if nonmissingtype(T) != T
        values,ismissingvalues = defaultmissing(values_or_missing)
    else
        values,ismissingvalues = values_or_missing,fill(false, size(values_or_missing))
    end
    return PairParam(name, components, values, ismissingvalues)
end

function PairParam(name, components, values_or_missing::AbstractVector{T}) where T 
    if nonmissingtype(T) != T
        vec_values,ismissingvalues = defaultmissing(values_or_missing)
        pairvalues = singletopair(vec_values,missing)
        for i in 1:length(vec_values)
            if vec_ismissingvalues[i]
                pairvalues[i,i] = missing
            end
        end
        values,ismissingvalues = defaultmissing(pairvalues)
    else
        values = singletopair(values_or_missing)
        ismissingvalues = fill(true,size(values))
        for i in 1:length(values_or_missing)
            ismissingvalues[i,i] = false
        end
    end
   
    return PairParam(name, components, values, ismissingvalues)
end

# If no value is provided, just initialise empty param.
function PairParam{T}(name,components) where T <: Number
    nc = length(components)
    values = fill(zero(T), (nc,nc))
    ismissingvalues = fill(true,(nc,nc))
    return PairParam(name, components, values, ismissingvalues)
end

PairParam(name,components) = PairParam{Float64}(name,components)

function PairParam(x::SingleParam,name::String=x.name)
    pairvalues = singletopair(x.values,missing)
    for i in 1:length(x.values)
        if x.ismissingvalues[i]
            pairvalues[i,i] = missing
        end
    end
    _values,_ismissingvalues = defaultmissing(pairvalues)
    return PairParam(name, x.components, _values, _ismissingvalues, x.sourcecsvs, x.sources)
end

function PairParam(x::PairParam,name::String = x.name)
    return PairParam(name, x.components, deepcopy(x.values), deepcopy(x.ismissingvalues), x.sourcecsvs, x.sources)
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
function Base.convert(::Type{PairParam{T1}},param::PairParam{T2}) where {T1<:Number,T2<:Number}
    values = convert(Matrix{T1},param.values)
    return PairParam(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end
