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

"""
    concrete(x)

Given an array of heterogeneous values, return an array of concrete values.
"""
concrete(x::Vector{Float64}) = x 
concrete(x::Vector{Int64}) = x
concrete(x::Vector{String}) = x
concrete(x::Vector{Bool}) = x 
concrete(x) = convert(AbstractArray{mapreduce(typeof, promote_type, x)}, x);