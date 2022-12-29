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


function show_pairs(io,keys,vals=nothing,separator="",f_print = print;quote_string = true,pair_separator = '\n',prekey = ifelse(pair_separator === '\n'," ",""))
    if length(keys) == 0
        return nothing
    end
    if vals === nothing #useful for printing only keys
        vals = Iterators.repeated("")
    end
    i = 0
    for (k,v) in zip(keys,vals)
        i += 1
        if i > 1
            print(io,pair_separator)
        end
        if quote_string
            quot = '\"'
            print(io,prekey,quot,k,quot,separator)
        else
            print(io,prekey,k,separator)
        end
        f_print(io,v)
    end
end

function _vecparser_eltype(vals)
    for val in eachsplit(vals,' ')
        if isnothing(tryparse(Int,val))
            return Float64 
        end
    end
    return Int
end

function _vecparser(T::Type{X},vals::String,dlm = ' ') where X <: Union{Int,Float64}
    strip_vals = strip(vals,('[',']'))
    res = Vector{T}(undef,0)
    for strval in eachsplit(strip_vals,dlm,keepempty = false)
        val = tryparse(T,strval)
        if !isnothing(val)
            push!(res,val)
        else
            colors = Base.text_colors
            red = colors[:bold] * colors[:red]
            reset = colors[:normal]
            errval =  red * strval * reset
            error("cannot parse $errval as a number in $vals")
        end
    end
    return res
end

function _vecparser(vals::String,dlm = ' ')
    strip_vals = strip(vals,('[',']'))
    T = _vecparser_eltype(strip_vals)
    return _vecparser(T,vals,dlm)
end
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
