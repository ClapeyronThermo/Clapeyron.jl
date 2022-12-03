#=
This file contains functionalities present in newer versions of julia, but not on LTS,
that are used by this package.
=#


#Copied from ArrayInterfaceCore.jl
Base.@pure __parameterless_type(T) = Base.typename(T).wrapper

"""
    parameterless_type(x)
given a type or instance of type, it will return the type without any parameters attached.
## Examples:
```julia-repl
julia> Clapeyron.parameterless_type(Vector{Float64})
Array
julia> Clapeyron.parameterless_type(rand(1:5,10))
Array
```
"""
parameterless_type(x) = parameterless_type(typeof(x))
parameterless_type(x::Type) = __parameterless_type(x)

#this is never used in a critical path, so we just use a default copying method
if !isdefined(Base,:keepat!)
    function keepat!(a,inds)
        b = a[inds]
        resize!(a,length(b))
        a .= b
        return a
    end
end

if !isdefined(Base,:eachsplit)
    eachsplit(str::AbstractString, dlm; limit::Integer=0, keepempty::Bool=true) = split(str,dlm;limit,keepempty)
    eachsplit(str::AbstractString; limit::Integer=0, keepempty::Bool=false)  = split(str;limit,keepempty)
end

split_2(str) = NTuple{2}(eachsplit(str, limit=2))
split_2(str,dlm) = NTuple{2}(eachsplit(str,dlm, limit=2))

#=
"""
    concrete(x)

Given an array of heterogeneous values, return an array of concrete values.
"""
concrete(x::Vector{Float64}) = x
concrete(x::Vector{Int64}) = x
concrete(x::Vector{String}) = x
concrete(x::Vector{Bool}) = x
concrete(x) = convert(AbstractArray{mapreduce(typeof, promote_type, x)}, x)=#