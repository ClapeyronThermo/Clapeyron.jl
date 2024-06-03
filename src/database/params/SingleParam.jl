struct SingleParameter{T,V<:AbstractVector{T}} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::V
    ismissingvalues::Array{Bool,1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

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

Base.eltype(param::SingleParameter{T}) where T = T

#linear algebra

LinearAlgebra.dot(param::SingleParameter,x::Union{<:AbstractVector,<:Number}) = dot(param.values,x)
LinearAlgebra.dot(x::Union{<:AbstractVector,<:Number},param::SingleParameter) = dot(x,param.values)

function SingleParam(name,components,values,missingvals,src,sourcecsv) 
    param_length_check(SingleParam,name,length(components),length(values))
    SingleParameter(name,components,values,missingvals,src,sourcecsv)
end

function Base.show(io::IO, ::MIME"text/plain", param::SingleParameter)
    len = length(param.values)
    print(io, "SingleParam{",eltype(param.values), "}(\"", param.name)
    println(io, "\") with ", len, " component", ifelse(len==1, ":", "s:"))
    separator = " => "
    vals = [ifelse(m,missing,v) for (m,v) in zip(param.ismissingvalues, param.values)]
    show_pairs(io,param.components,vals,separator)
end

function SingleParam(x::SingleParam, name::String = x.name; isdeepcopy::Bool = true, sources::Vector{String} = x.sources)
    if isdeepcopy
        return SingleParam(
            name,
            x.components,
            deepcopy(x.values),
            deepcopy(x.ismissingvalues),
            x.sourcecsvs,
            sources
        )
    end
    return SingleParam(
        name,
        x.components,
        x.values,
        x.ismissingvalues,
        x.sourcecsvs,
        sources
    )
end

SingleParameter(x::SingleParam, name::String = x.name; isdeepcopy::Bool = true, sources::Vector{String} = x.sources) = SingleParam(x, name; isdeepcopy, sources)

#a barebones constructor, in case we dont build from csv
function SingleParam(
        name::String,
        components::Vector{String},
        values::Vector{T},
        sourcecsvs = String[],
        sources = String[]
    ) where T
    param_length_check(SingleParam,name,length(components),length(values))
    if any(ismissing, values)
        _values,_ismissingvalues = defaultmissing(values)
        TT = eltype(_values)
    else
        _values = values
        _ismissingvalues = fill(false, length(values))
        TT = T
    end
    return  SingleParam{TT}(name, components, _values, _ismissingvalues, sourcecsvs, sources)
end

# If no value is provided, just initialise empty param.
function SingleParam(
        name::String,
        components::Vector{String};
        sources = String[]
    )
    values = fill(0.0, length(components))
    missingvalues = fill(true,length(components))
    return SingleParam(name, components, values,missingvalues, String[], sources)
end

function SingleParam(oldparam::SingleParameter, v::Vector)
    _values,_ismissingvalues = defaultmissing(v)
    param_length_check(SingleParam,name,length(oldparam.components),length(_values))
    return SingleParam(oldparam.name, oldparam.components,_values, _ismissingvalues , oldparam.sourcecsvs, oldparam.sources)
end

#convert utilities
function Base.convert(::Type{SingleParam{T1}},param::SingleParam{T2}) where {T1<:Number,T2<:Number}
    values = T1.(param.values)
    return SingleParam(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{SingleParam{Bool}},param::SingleParam{<:Union{Int,Float64,Bool}})
    #@assert all(z->(isone(z) | iszero(z)),param.values)
    values = Array(Bool.(param.values))
    return SingleParam(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{SingleParam{String}},param::SingleParam{<:AbstractString})
    values = String.(param.values)
    return SingleParameter(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

#trying to break stack overflow on julia 1.6
function Base.convert(::Type{SingleParam{String}},param::SingleParam{String})
    return param
end


#pack vectors

const PackedVectorSingleParam{T} = Clapeyron.SingleParameter{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}, PackedVectorsOfVectors.PackedVectorOfVectors{Vector{Int64}, Vector{T}, SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true}}}

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