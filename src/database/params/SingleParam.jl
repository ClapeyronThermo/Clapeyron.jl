struct SingleParameter{T,V<:AbstractVector{T}} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::V
    ismissingvalues::Array{Bool,1}
    sourcecsvs::Union{Vector{String},Nothing}
    sources::Union{Vector{String},Nothing}
end

#methods to fill missing sources/sourcescsvs
SingleParameter(name,components,values,ismissingvalues) = SingleParameter(name,components,values,ismissingvalues,nothing,nothing)
SingleParameter(name,components,values,ismissingvalues,src) = SingleParameter(name,components,values,ismissingvalues,src,nothing)


"""
    SingleParam{T}
Struct designed to contain single parameters. Basically a vector with some extra info.
## Creation:
```julia-repl
julia> mw = SingleParam("molecular weight",["water","ammonia"],[18.01,17.03])
SingleParam{Float64}("molecular weight") with 2 components:
 "water" => 18.01
 "ammonia" => 17.03
julia> mw.values
2-element Vector{Float64}:
 18.01
 17.03
julia> mw.components
2-element Vector{String}:
 "water"
 "ammonia"
julia> mw2 = SingleParam(mw,"new name")
SingleParam{Float64}("new name") with 2 components:
 "water" => 18.01
 "ammonia" => 17.03
julia> has_oxigen = [true,false]; has_o = SingleParam(mw2,has_oxigen)
SingleParam{Bool}("new name") with 2 components:
 "water" => true
 "ammonia" => false
```
## Example usage in models:
```
function molecular_weight(model,molar_frac)
    mw = model.params.mw.values
    res = zero(eltype(molarfrac))
    for i in @comps #iterating through all components
        res += molar_frac[i]*mw[i]
    end
    return res
end
```
"""
const SingleParam{T} = SingleParameter{T,Vector{T}} where T

#indexing

Base.@propagate_inbounds Base.getindex(param::SingleParameter{T,<:AbstractVector{T}},i::Int) where T = param.values[i]

function Base.getindex(param::SingleParameter{T,<:AbstractVector{T}},i::AbstractString) where T
    idx = _str_to_idx(param,i)
    return param[idx]
end

Base.setindex!(param::SingleParameter,val,i::Integer) = setindex!(param.values,val,i)
function Base.setindex!(param::SingleParameter,val,i::AbstractString)
    idx = _str_to_idx(param,i)
    setindex!(param,val,idx)
end
#broadcasting
Base.size(param::SingleParameter) = size(param.values)
Base.length(param::SingleParameter) = length(param.values)
Base.broadcastable(param::SingleParameter) = param.values
Base.BroadcastStyle(::Type{<:SingleParameter}) = Broadcast.Style{SingleParameter}()

#copyto!
function Base.copyto!(dest::SingleParameter,src::Base.Broadcast.Broadcasted) #general, just copies the values, used in a .= f.(a)
    Base.copyto!(dest.values,src)
    return dest
end

function Base.copyto!(dest::SingleParameter,src::AbstractArray) #general, just copies the values, used in a .= f.(a)
    Base.copyto!(dest.values,src)
    return dest
end

function Base.copyto!(dest::SingleParameter,src::SingleParameter) #used to set params
    #key check
    dest.components == src.components || throw(DimensionMismatch("components of source and destination single parameters are not the same for $dest"))
    copyto!(dest.values,src.values)
    dest.ismissingvalues .= src.ismissingvalues
    return dest
end

#primalval
function Solvers.primalval(x::SingleParameter)
    return SingleParameter(x.name,x.components,Solvers.primalval_eager(x.values),x.ismissingvalues,x.sourcecsvs,x.sources)
end

Base.eltype(param::SingleParameter{T}) where T = T
Base.eltype(param::Type{<:SingleParameter{T}}) where T = T

#linear algebra

LinearAlgebra.dot(param::SingleParameter,x::Union{<:AbstractVector,<:Number}) = dot(param.values,x)
LinearAlgebra.dot(x::Union{<:AbstractVector,<:Number},param::SingleParameter) = dot(x,param.values)

#barebones constructor, we provide vals and missing vals
function SingleParam(name,components,values::Vector{T},missingvals,src,sourcecsv) where T 
    param_length_check(SingleParam,name,length(components),length(values))
    SingleParameter{T,Vector{T}}(name,components,values,missingvals,src,sourcecsv)
end

function SingleParam(name,components,values::AbstractVector{T},missingvals,src,sourcecsv) where T 
    return SingleParam(name,components,convert(Vector{T},values),missingvals,src,sourcecsv)
end

function SingleParam(name,components,values::Matrix{T},missingvals::Matrix{Bool},src,sourcecsv) where T 
    n = length(components)
    vec_values = zeros(T,n)
    vec_missing = zeros(Bool,n)
    for i in 1:n
        vec_values[i] = values[i,i]
        vec_missing[i] = missingvals[i,i]
    end

    return SingleParam(name,components,vec_values,vec_missing,src,sourcecsv)
end

SingleParam(name,components,values,missingvals,src) = SingleParam(name,components,values,missingvals,src,nothing)
SingleParam(name,components,values,missingvals) = SingleParam(name,components,values,missingvals,nothing,nothing)

#constructor in case we provide a normal vector
function SingleParam(name, components, values_or_missing::AbstractVector{T}) where T
    if nonmissingtype(T) != T
        values,ismissingvalues = defaultmissing(values_or_missing)
    else
        values,ismissingvalues = values_or_missing,fill(false, length(values_or_missing))
    end
    return SingleParam(name, components, values, ismissingvalues)
end

# If no value is provided, just initialise empty param.
function SingleParam{T}(name,components) where T<: Number
    nc = length(components)
    values = fill(zero(T),nc)
    ismissingvalues = fill(true,nc)
    return SingleParam(name, components, values, ismissingvalues)
end

SingleParam(name,components) = SingleParam{Float64}(name,components) 

function Base.show(io::IO, ::MIME"text/plain", param::SingleParameter)
    len = length(param.values)
    print(io, "SingleParam{",eltype(param.values), "}(\"", param.name)
    println(io, "\") with ", len, " component", ifelse(len==1, ":", "s:"))
    separator = " => "
    vals = [ifelse(m,missing,v) for (m,v) in zip(param.ismissingvalues, param.values)]
    show_pairs(io,param.components,vals,separator)
end

function SingleParam(oldparam::SingleParameter, v::Vector)
    _values,_ismissingvalues = defaultmissing(v)
    param_length_check(SingleParam,name,length(oldparam.components),length(_values))
    return SingleParam(oldparam.name, oldparam.components, _values, _ismissingvalues , oldparam.sourcecsvs, oldparam.sources)
end

#convert utilities
function Base.convert(::Type{SingleParam{T1}},param::SingleParam{T2}) where {T1<:Number,T2<:Number}
    values = convert(Vector{T1},param.values)
    return SingleParam(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{SingleParam{String}},param::SingleParam{<:AbstractString})
    values = convert(Vector{String},param.values)
    return SingleParameter(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

#pack vectors
const PackedVectorSingleParam{T} = Clapeyron.SingleParameter{PackedSubVector{T}, PackedVector{T}}

function pack_vectors(param::SingleParameter{<:AbstractVector})
    name,components,vals,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    vals = pack_vectors(vals)
    return SingleParam(name,components,vals,missingvals,srccsv,src)
end

function pack_vectors(params::Vararg{SingleParameter{T},N}) where {T<:Number,N}
    param = first(params)
    name,components,_,missingvals,srccsv,src = param.name,param.components,param.values,param.ismissingvalues,param.sourcecsvs,param.sources
    len = length(params)
    vals = [zeros(len) for _ in params]

    for i in 1:length(vals)
        vali = vals[i]
        for (k,par) in pairs(params)
            vali[k] = par.values[i]
        end
    end
    vals = PackedVectorsOfVectors.pack(vals)
    SingleParam(name,components,vals,missingvals,srccsv,src)
end

#= Operations
function Base.:(+)(param::SingleParameter, x::Number)
    values = param.values .+ x
    return SingleParam(param.name, param.components, values, param.ismissingvalues, param.sourcecsvs, param.sources)
end
function Base.:(*)(param::SingleParameter, x::Number)
    values = param.values .* x
    return SingleParam(param.name, param.components, values, param.ismissingvalues, param.sourcecsvs, param.sources)
end
function Base.:(^)(param::SingleParameter, x::Number)
    values = param.values .^ x
    return SingleParam(param.name, param.components, values, param.ismissingvalues, param.sourcecsvs, param.sources)
end =#